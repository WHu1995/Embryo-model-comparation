# projection.R
# Project query scRNA-seq data onto the integrated mouse reference atlas

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(plotly)
library(tidyr)

###############################################################################
# Step 1. Set input objects
###############################################################################

# Set your query Seurat object
query <- your_query_seurat_object

# Set your reference Seurat object
Mouseref <- your_reference_map
##Or download pre-built Mouseref in https://zenodo.org/records/19476835

###############################################################################
# Step 2. Optional QC filtering for query data
###############################################################################
VlnPlot(
  query,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  raster = FALSE,
  pt.size = 0,
  group.by = "orig.ident"
)

# Modify thresholds as needed for your dataset.
# Higher-quality cells (e.g., lower mitochondrial content, appropriate feature counts)
# generally lead to more accurate projection results.

query <- subset(
  query,
  subset = nFeature_RNA > 700 & nFeature_RNA < 10000 & percent.mt < 10
)

###############################################################################
# Step 3. Calculate transfer anchors
###############################################################################
options(future.globals.maxSize = 10 * 1024^3)
anchors <- FindTransferAnchors(
  reference = Mouseref,
  query = query,
  dims = 1:40,
  reference.reduction = "pca"
  # optional: k.anchor = 10
)

###############################################################################
# Step 4. label transfer
###############################################################################

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

# Keep original query annotation for downstream comparison
A_mapped$original_annotation <- query$orig.ident



###############################################################################
# Visualization 
###############################################################################

###############################################################################
# Step 5. Define LowConf cells
###############################################################################
#according to celltype score 
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


###############################################################################
# Step 6. Define lineage groups in the reference
###############################################################################
Mouseref$lineage_group <- dplyr::case_when(
  Mouseref$cell_type %in% c("EPI", "Pre_EPI", "ICM", "Caudal.epi") ~ "Epiblast",
  Mouseref$cell_type %in% c("PGC") ~ "PGC",
  Mouseref$cell_type %in% c(
    "TE", "polar TE", "mural TE", "ExE", "EPC", "EPC.mig",
    "Pri.pTGC", "Sec.pTGC", "SpT", "SynTI.p"
  ) ~ "Trophectoderm",
  Mouseref$cell_type %in% c(
    "Brain", "Neuralectoderm", "Spinal cord",
    "Neural Crest", "Sur.ecto", "NMP", "Axial midline"
  ) ~ "Ectoderm",
  Mouseref$cell_type %in% c(
    "Mesoderm", "Caudal.meso", "Nas.meso", "Para.meso",
    "Phary.meso", "Mesenchyme", "ExM", "PS", "Allantois"
  ) ~ "Mesoderm",
  Mouseref$cell_type %in% c(
    "Def.endo", "Gut", "PE", "VE", "DVE", "EmVE",
    "ExVE", "PaE"
  ) ~ "Endoderm",
  Mouseref$cell_type %in% c(
    "Blood.prog", "Erythroid", "Cardiomyocytes",
    "Endothelium", "HEPs"
  ) ~ "Cardiovascular",
  TRUE ~ "Uncategorized"
)

table(Mouseref$lineage_group)


###############################################################################
# Step 7. Define cell type colors
###############################################################################
lineage_colors <- list(
  Cardiovascular = colorRampPalette(c("#FFC0CB", "#D2665A"))(6),
  Epiblast = colorRampPalette(c("#727D73", "#B3C8CF", "#89A8B2", "#9BA8A0"))(3),
  Mesoderm = c(
    "#D0F0C0", "#A8D8A5", "#80C090", "#50A080", "#20806A",
    "#187060", "#71BBB2", "#5FA9A6", "#4D979B"
  ),
  Ectoderm = c(
    "#D8EFFF", "#B0E2FF", "#48B9FF", "#6FB6D9",
    "#5A9EC7", "#326698", "#1E4A7C"
  ),
  Endoderm = c(
    "#FFF3E0", "#FFE0B2", "#FFCC80", "#FFB74D",
    "#FF9800", "#EF6C00", "#D84315", "#BF360C"
  ),
  Trophectoderm = c(
    "#2E1A47", "#673AB7", "#D1C4E9", "#F8BBD0", "#F48FB1",
    "#EC407A", "#AD1457", "#880E4F", "#AB47BC", "#E1BEE7"
  ),
  PGC = colorRampPalette(c("#1FB6C4"))(2),
  Uncategorized = "#000000"
)

color_mapping <- unlist(lapply(names(lineage_colors), function(lineage) {
  types <- Mouseref@meta.data %>%
    filter(lineage_group == lineage) %>%
    pull(cell_type) %>%
    unique()
  
  if (length(types) == 0) return(NULL)
  setNames(lineage_colors[[lineage]][seq_along(types)], types)
}))
###############################################################################
# Step 8. UMAP visualization of projected cells
###############################################################################
# Reference UMAP
ref_umap <- Embeddings(Mouseref, reduction = "umap")
ref_umap_df <- data.frame(
  UMAP_1 = ref_umap[, 1],
  UMAP_2 = ref_umap[, 2],
  source = "reference"
)

# Query projected UMAP
query_umap <- Embeddings(A_mapped, reduction = "ref.umap")
query_umap_df <- data.frame(
  UMAP_1 = query_umap[, 1],
  UMAP_2 = query_umap[, 2],
  source = "query"
)

query_umap_df$original_annotation <- A_mapped$original_annotation
query_umap_df$predicted.celltype <- A_mapped$predicted.celltype
query_umap_df$predicted.cellstage <- A_mapped$predicted.cellstage

ref_umap_df$original_annotation <- NA
ref_umap_df$predicted.celltype <- NA
ref_umap_df$predicted.cellstage <- NA

combined_umap_df <- rbind(ref_umap_df, query_umap_df)

ggplot(combined_umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    data = combined_umap_df[combined_umap_df$source == "reference", ],
    color = "lightgrey",
    alpha = 0.4,
    size = 0.5
  ) +
  geom_point(
    data = combined_umap_df[combined_umap_df$source == "query", ],
    aes(fill = predicted.celltype),
    shape = 21,
    color = "black",
    size = 3,
    stroke = 0.5
  ) +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +
  labs(fill = "Predicted cell type")

###############################################################################
# Step 9. River plot: original annotation -> predicted cell type
###############################################################################
# 2. 设置从上到下的顺序
lineage_order <- c(
  "Cardiovascular", "Ectoderm", "Endoderm",
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
  Cardiovascular = "#FAB6BE",
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
pdf("river_plot.pdf",width = 4, height = 6)
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

###############################################################################
# Step 10. Optional interactive Sankey plot
###############################################################################
sankey_df <- A_mapped@meta.data %>%
  count(original_annotation, predicted.celltype) %>%
  as.data.frame()

labels <- c(unique(sankey_df$original_annotation), unique(sankey_df$predicted.celltype))

sankey_df$source <- match(sankey_df$original_annotation, labels) - 1
sankey_df$target <- match(sankey_df$predicted.celltype, labels) - 1

fig <- plot_ly(
  type = "sankey",
  arrangement = "snap",
  node = list(
    label = labels,
    pad = 15,
    thickness = 20,
    line = list(color = "black", width = 0.5)
  ),
  link = list(
    source = sankey_df$source,
    target = sankey_df$target,
    value = sankey_df$n
  )
)

fig

###############################################################################
# Step 11. Heatmap of predicted cell type vs predicted developmental stage
###############################################################################
plot_df <- A_mapped@meta.data %>%
  dplyr::count(predicted.cellstage, predicted.celltype, name = "cell_number") %>%
  tidyr::complete(predicted.cellstage, predicted.celltype, fill = list(cell_number = 0))

ggplot(plot_df, aes(x = predicted.celltype, y = predicted.cellstage, fill = cell_number)) +
  geom_tile(color = "grey50", linewidth = 0.4) +
  geom_text(aes(label = cell_number), size = 5) +
  scale_fill_gradient(low = "white", high = "#5B8DBB") +
  labs(
    x = "Predicted cell type",
    y = "Predicted cell stage",
    fill = "Cell\nnumber"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key = element_blank()
  )
