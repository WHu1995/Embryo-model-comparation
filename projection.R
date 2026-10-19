library(Seurat)
library(ggplot2)
library(dplyr)

  
query <- readRDS("C:/MZGlab/GAN paper/revise/3rd_edition/GANfigures/Cebpa_OE_real_72h_GSM8678833_2000_matched_genes.rds")
# 重新建立Seurat v4对象
counts <- GetAssayData(query, assay = "RNA", slot = "counts")
meta <- query@meta.data[colnames(counts), , drop = FALSE]
query <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = ""
)

table(query$condition)


folder <- "C:/MZGlab/GAN paper/revise/3rd_edition/cebpa_wei/12h/"
data <- Read10X(folder)

query <- CreateSeuratObject(
  counts = data[["Gene Expression"]],
  min.cells = 3
)
query<- subset(query, orig.ident == "real")


query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^mt-")
VlnPlot(query, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,raster=FALSE,pt.size = 0,group.by = "orig.ident")

query <- subset(query, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 5& nCount_RNA < 15000)

query <- NormalizeData(query)
query <- ScaleData(query)
query<- FindVariableFeatures(query, nfeatures = 2000)
query <- RunPCA(query, verbose = FALSE)

options(future.globals.maxSize = 10 * 1024^3)

###calculate anchors
anchors <- FindTransferAnchors(
  reference = Mouseref,
  query = query,
  dims = 1:40,
  reference.reduction = "pca")#optional set:,k.anchor =10

###label transfer
A_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = Mouseref,
  refdata = list(
    celltype = "cell_type",
    cellstage = "orig.ident",
    lineage = "lineage_group"
  ),
  reference.reduction = "pca",
  reduction.model = "umap"
)
# river plot

library(ggplot2)
library(ggalluvial)

# 1. 根据 lineage score 判断 LowConf
hist(
  A_mapped$predicted.celltype.score,
  breaks = 30,
  main = "Predicted Lineage Group Score",
  xlab = "Prediction Score",
  ylab = "Number of Cells"
)

low <- is.na(A_mapped$predicted.celltype.score) |
  A_mapped$predicted.celltype.score < 0.4

celltype <- ifelse(low, "LowConf",
                   as.character(A_mapped$predicted.lineage))
stage <- ifelse(low, "LowConf",
                as.character(A_mapped$predicted.cellstage))

# 2. 设置从上到下的顺序
lineage_order <- c(
  "Cardiovascular_system", "Ectoderm", "Endoderm",
  "Epiblast", "Mesoderm", "Trophectoderm",
  "PGC", "LowConf"
)

stage_order <- c(
  "E3.5", "E4.5", "E4.75", "E5.0", "E5.25", "E5.5",
  "E6.25", "E6.5", "E6.75", "E7.0", "E7.25",
  "E7.5", "E7.75", "E8.0", "E8.25", "E8.5",
  "E8.75", "E9.5", "LowConf"
)

# ggalluvial 默认反向堆叠
celltype <- factor(celltype, levels = rev(lineage_order))
stage <- factor(stage, levels = rev(stage_order))

# 3. 统计
df <- as.data.frame(table(celltype, stage))
colnames(df) <- c("Lineage", "Stage", "Count")
df <- df[df$Count > 0, ]

# 4. 颜色
lineage_cols <- c(
  Cardiovascular_system = "#FAB6BE",
  Ectoderm = "#6FB6D9",
  Endoderm = "#F7971D",
  Epiblast = "#88A7B1",
  Mesoderm = "#80C090",
  Trophectoderm = "#B19CCA",
  PGC = "#E76BF3",
  LowConf = "#999999"
)

stage_cols <- c(
  E3.5="#F8766D", E4.5="#F3C389", E4.75="#FCBBB6",
  E5.0="#D39200", E5.25="#E9C980", E5.5="#7CAE00",
  E6.25="#80DD9C", E6.5="#AFD980", E6.75="#00C19F",
  E7.0="#80DD9C", E7.25="#CE9AF6", E7.5="#B983FF",
  E7.75="#619CFF", E8.0="#B983FF", E8.25="#ED68ED",
  E8.5="#FF61CC", E8.75="#D63384",
  E9.5="black"
)

# 5. 画河流图
pdf("cebpa_wei_12h.pdf",width = 4, height = 6)
ggplot(df, aes(axis1=Lineage, axis2=Stage, y=Count)) +
  geom_alluvium(
    aes(fill=Lineage),
    width=0.2, alpha=0.75,
    reverse=FALSE
  ) +
  geom_stratum(
    aes(fill=after_stat(stratum)),
    width=0.2, color="black",
    reverse=FALSE
  ) +
  geom_text(
    stat="stratum",
    aes(label=after_stat(stratum)),
    size=3, reverse=FALSE
  ) +
  scale_fill_manual(values=c(lineage_cols, stage_cols)) +
  scale_x_discrete(
    limits=c("Lineage", "Stage"),
    expand=c(0.2,0.2)
  ) +
  theme_void() +
  theme(legend.position="none")
dev.off()


