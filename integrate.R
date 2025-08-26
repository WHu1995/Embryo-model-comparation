install.packages("matrixStats")
library(batchelor)
library(Seurat)
library(SingleCellExperiment)
library(scuttle)


mat1 <- GetAssayData(test, assay = "RNA", slot = "count")
mat2 <- GetAssayData(post, assay = "RNA", slot = "count")
mat3 <- GetAssayData(placenta_subset, assay = "RNA", slot = "count")
mat4 <- GetAssayData(Cheng, assay = "RNA", slot = "count")
mat5 <- GetAssayData(Mohammed, assay = "RNA", slot = "count")
mat6 <- GetAssayData(Nowotschin, assay = "RNA", slot = "count")

#################
library(Matrix)
#Take the overlapping genes, as these datasets share a large intersection
common_genes <- Reduce(intersect, list(rownames(mat1), rownames(mat2),rownames(mat3), 
                                      rownames(mat4)))
common_genes <- Reduce(intersect, list(rownames(mat1), rownames(mat2),rownames(mat3), 
                                       rownames(mat4),rownames(mat5),rownames(mat6)))

mat1_common <- mat1[common_genes, ]
mat2_common <- mat2[common_genes, ]
mat3_common <- mat3[common_genes, ]
mat4_common <- mat4[common_genes, ]
merged_mat <- cbind(mat1_common, mat2_common, mat3_common, mat4_common)

#Based on the predefined set of common genes, we extracted corresponding expression values 
#from the Nowotschin or Mohammed dataset, assigning zero to genes absent from the dataset.
counts_5 <- GetAssayData(Mohammed, assay = "RNA", slot = "counts")
existing_genes <- intersect(common_genes, rownames(counts_5))
counts_5_filtered <- counts_5[existing_genes, ]
missing_genes <- setdiff(common_genes, existing_genes)
if (length(missing_genes) > 0) {
  zero_block <- Matrix(0, 
                       nrow = length(missing_genes),
                       ncol = ncol(counts_5),
                       dimnames = list(missing_genes, colnames(counts_5)))
  counts_5_aligned <- rbind(counts_5_filtered, zero_block)
} else {
  counts_5_aligned <- counts_5_filtered
}

counts_5_aligned <- counts_5_aligned[common_genes, ]
#reconstruct seurat obj
Mohammed <- CreateSeuratObject(
  counts = counts_5_aligned,
  assay = "RNA",
  meta.data = Mohammed@meta.data  # 保留原metadata
)


mat5 <- GetAssayData(Mohammed, assay = "RNA", slot = "count")
mat6 <- GetAssayData(Nowotschin, assay = "RNA", slot = "count")

#Merge total expression matrices
final_merged <- cbind(merged_mat,mat5,mat6)



metadata1 <- test@meta.data
metadata2 <- post@meta.data
metadata3 <- placenta_subset@meta.data
metadata4 <- Cheng@meta.data
metadata5 <- Mohammed@meta.data
metadata6 <- Nowotschin@meta.data

metadata1<- metadata1[, c("orig.ident", "cell_type","dataset")]
metadata2 <- metadata2[, c("orig.ident", "cell_type","dataset")]
metadata3 <- metadata3[, c("orig.ident", "cell_type","dataset")]
metadata4 <- metadata4[, c("orig.ident", "cell_type","dataset")]
metadata5 <- metadata5[, c("orig.ident", "cell_type","dataset")]
metadata6 <- metadata6[, c("orig.ident", "cell_type","dataset")]
#Merge total metadata
merged_metadata <- rbind(metadata1, metadata2,metadata3,metadata4,metadata5,metadata6)

#re-construct total Seurat obj
save(final_merged,file="merged_raw_count.RData")
save(merged_metadata,file="merged_metadata.RData")

Moss <- CreateSeuratObject(counts = final_merged, meta.data = merged_metadata)
Moss <- NormalizeData(Moss)
Moss <- ScaleData(Moss)
Moss<- FindVariableFeatures(Moss, nfeatures = 2000)

GetAssayData(Moss, assay = "RNA", slot = "counts")[1:5, 1:5]
################integrat dataset using FastMNN##########################################
#select 2000 integration features
sceList <- list(test,post,placenta_subset,Mohammed,Cheng,Nowotschin)
sce_filtered_list <- lapply(sceList, function(obj) {
  subset(obj, features = common_genes)
})
features <- SelectIntegrationFeatures(object.list = sce_filtered_list,nfeatures = 2000)

sce1 <- SingleCellExperiment(list(counts = mat1), colData = DataFrame(metadata1))
sce2 <- SingleCellExperiment(list(counts = mat2), colData = DataFrame(metadata2))
sce3 <- SingleCellExperiment(list(counts = mat3), colData = DataFrame(metadata3))
sce4 <- SingleCellExperiment(list(counts = mat4), colData = DataFrame(metadata4))
sce5 <- SingleCellExperiment(list(counts = mat5), colData = DataFrame(metadata5))
sce6 <- SingleCellExperiment(list(counts = mat6), colData = DataFrame(metadata6))

sce1 <- sce1[features, ]
sce2 <- sce2[features, ]
sce3 <- sce3[features, ]
sce4 <- sce4[features, ]
sce5 <- sce5[features, ]
sce6 <- sce6[features, ]

sce1 <- logNormCounts(sce1)
sce2 <- logNormCounts(sce2)
sce3 <- logNormCounts(sce3)
sce4 <- logNormCounts(sce4)
sce5 <- logNormCounts(sce5)
sce6 <- logNormCounts(sce6)

#integrat using fastMNN
mnn_result <- fastMNN(
  sce1, sce2, sce3,  sce4, sce5, sce6,      
  d = 50,            
  k = 20             
)
assay(mnn_result, "counts") <- 2^assay(mnn_result, "logcounts") - 1
mmn <- as.Seurat(mnn_result, data = "logcounts")


table(mnn_result)

#add metadata info
save(mnn_result,file="merged_mnn_result.RData")

colData(mnn_result) <- cbind(colData(mnn_result), merged_metadata)

corrected_pca <- reducedDim(mnn_result, "corrected")


#add pca embeddings
Moss[["pca"]] <- CreateDimReducObject(embeddings = corrected_pca, key = "PC_", assay = "RNA")
#run UMAP and save
Moss <- RunUMAP(Moss, reduction = "pca", dims = 1:40,return.model = TRUE,seed.use =666)


save(Moss, file = "Moss.RData")

#Visualization
stage_colors <- c(
  "E3.5" = "#8b0000", 
  "E4.5" = "#bd0026",
  "E4.75" = "#e31a1c",
  "E5.0" = "#fc4e2a",
  "E5.25" = "#fd6e33",
  "E5.5" = "#fd8d3c",
  "E6.25" = "#feb24c",
  "E6.5" = "#fed976",  
  "E6.75" = "#FFE38B",
  "E7.0" = "#c7e9b4",  
  "E7.25" = "#7fcdbb", 
  "E7.5" = "#41b6c4",
  "E7.75" = "#1d91c0",
  "E8.0" = "#6baed6",
  "E8.25" = "#4292c6",
  "E8.5" = "#266abd",
  "E8.75" = "#225ea8"
)

DimPlot(Moss, reduction = "umap", group.by = "orig.ident", cols = stage_colors, raster = FALSE)
DimPlot(Moss, reduction = "umap", label = T, group.by = "cell_type",raster=FALSE) + NoLegend()
DimPlot(Moss, reduction = "umap", label = F, group.by = "dataset",raster=FALSE)

FeaturePlot(Moss, features = "Pou5f1",raster=FALSE)

table(Moss$cell_type)


####################Visualization with ggplot2##############
library(ggplot2)
# Extract data from DimPlot
umap_data <- as.data.frame(Embeddings(Moss, "umap"))
umap_data$cell_type <- Moss$cell_type
umap_data$orig.ident <- Moss$orig.ident
umap_data$dataset <- Moss$dataset

set.seed(123)
##optional: downsample to speed up the process
umap_sub <- umap_data[sample(nrow(umap_data), size = 50000), ]
library(RColorBrewer)
set2_colors <- brewer.pal(8, "Set2")
num_cell_types <- length(unique(umap_data$cell_type))
additional_colors <- colorRampPalette(colors = c("#FF6347", "#4682B4", "#32CD32", "#FFD700", "#8A2BE2", "#00CED1", "#DC143C", "#ADFF2F", "#FF4500"))(num_cell_types - length(set2_colors))
color_palette <- c(set2_colors, additional_colors)


ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = cell_type)) +
  geom_point(shape = 21, size = 2, color = "white", stroke = 0.01,alpha=1) +
  scale_fill_manual(values = color_palette) +  # Use a large Moss color palette
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  NoLegend()

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = dataset)) +
  geom_point(shape = 21, size = 2, color = "white", stroke = 0.01,alpha=1) +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2, fill = dataset)) +
  geom_point(shape = 21, size = 2, color = "white", stroke = 0.01, alpha=1) +
  scale_fill_manual(values = c("#4E79A7", "#EDC948", "#CC79A7", "#F28E2B", "#E15759", "#59A14F")) + # 任选上述方案
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right")



ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2, fill = orig.ident)) +
  geom_point(shape = 21, size = 1.8, color = "white", stroke = 0.01,alpha=1) +
  scale_fill_manual(values = stage_colors) +  # Use a large Moss color palette
  theme_minimal() +
  theme(panel.grid = element_blank()) 

##color the cell type by lineage:
#
library(Seurat)
library(dplyr)
library(ggplot2)
table(Moss$cell_type)
# annotated cell lineage
Moss@meta.data <- Moss@meta.data %>%
  mutate(
    lineage_group = case_when(
      grepl("Blood.prog|Erythroid|HEPs|Cardiomyocyte|Endotheli|CMPs", cell_type) ~ "circulatory_system",
      grepl("ICM|Pre_EPI|Epiblast|Caudal.epi", cell_type) ~ "Epiblast",
      grepl("Mesoderm|Caudal.meso|Inter.meso|Para.meso|Mixed.meso|ExM|Nas.meso|Mesenchyme|Phary.meso|PS", cell_type) ~ "Mesoderm",
      grepl("Neural|neurecto|Brain|Spinal|NMP|Notochord|Midline|Pr.ecto|Sur.ecto|Pri.ecto", cell_type) ~ "Ectoderm",
      grepl("endo|Gut|Def.endo|VE|EmVE|ExE.endo|Visceral|PrE|PaE|Allantois", cell_type) & 
        !grepl("Trophectoderm", cell_type) ~ "Endoderm",
      grepl("ExE|TE|pTGC|SynTI|EPC|SpT", cell_type) ~ "Trophectoderm",
      cell_type %in% c("PGC") ~ "PGC",
      TRUE ~ "Uncategorized"
    )
  )
table(Moss@meta.data$lineage_group)
lineage_composition <- Moss@meta.data %>%
  group_by(lineage_group) %>%
  summarise(
    CellTypes = paste(unique(cell_type), collapse = ", "),
    TotalCells = n()
  ) %>%
  arrange(desc(TotalCells))
View(lineage_composition)

##prepare lineage colors
lineage_colors <- list(
  circulatory_system  = colorRampPalette(c("#FFC0CB", "#D2665A"))(10),  
  Epiblast = colorRampPalette(c("#727D73", "#B3C8CF","#89A8B2"))(3),       
  Mesoderm = c(
    "#D0F0C0",  
    "#A8D8A5",  
    "#80C090",  
    "#68B08B",  
    "#50A080",  
    "#389075",  
    "#20806A",  
    "#187060",  
    "#71BBB2",  
    "#5FA9A6",  
    "#4D979B"   
  ),      
  Ectoderm =  c(
    "#D8EFFF",
    "#B0E2FF", 
    "#87CEEB", 
    "#60C2FF",
    "#48B9FF",
    "#6FB6D9", 
    "#5A9EC7",
    "#4682B4",  
    "#326698", 
    "#1E4A7C"
  ),        
  Endoderm = c(
    "#FFF3E0",
    "#FFE0B2",  
    "#FFCC80",  
    "#FFB74D",  
    "#FF9800",  
    "#EF6C00", 
    "#D84315",  
    "#BF360C", 
    "#8D6E63",
    "#780C28"
  ),       
  Ectoderm = colorRampPalette(c("#FFD700", "#FAFAD2"))(2),       
  Trophectoderm = c(
    "#2E1A47", 
    "#4A256E", 
    "#673AB7", 
    "#9575CD", 
    "#B39DDB", 
    "#D1C4E9",
    "#F8BBD0", 
    "#F48FB1", 
    "#EC407A", 
    "#D81B60", 
    "#AD1457", 
    "#880E4F"
  ),  
  PGC = colorRampPalette(c("#1FB6C4"))(1))

celltype_order <- levels(factor(Moss$cell_type))
color_mapping <- unlist(lapply(names(lineage_colors), function(lineage){
  types <- Moss@meta.data %>% 
    filter(lineage_group == lineage) %>% 
    pull(cell_type) %>% 
    unique()
  setNames(lineage_colors[[lineage]][1:length(types)], types)
}))
save(color_mapping,file = "color_mapping.RData")

# get UMAP Embeddings
umap_df <- data.frame(
  UMAP_1 = Embeddings(Moss, "umap")[,1],
  UMAP_2 = Embeddings(Moss, "umap")[,2],
  cell_type = Moss$cell_type,
  lineage_group = Moss$lineage_group
)
table(Moss$cell_type)
# Visualization
set.seed(123)  
umap_sub <- umap_df[sample(nrow(umap_df), size = 50000), ]
ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = cell_type), 
             shape = 21, 
             size = 2, 
             color = "white", 
             stroke = 0.2, 
             alpha = 0.9) +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.key.size = unit(5, "mm"),
    legend.text = element_text(size = 7),
    axis.line = element_line(linewidth = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(
    ncol = 2,
    override.aes = list(size = 3, alpha = 1),
    title = "Cell Type"
  ))



genes <-c("Sox2", "Gata3", "T", "Nanos3",
          "Pou5f1", "Eomes", "Mesp1", "Dppa3",
          "Nanog", "Elf5", "Tbx6", "Dnd1",
          "Otx2", "Gata4", "Klf1", "Ascl2",
          "Pax6", "Sox17", "Etv2", "Prl3d1",
          "Lhx2", "Foxa2", "Gypa", "Plac8"
)

library(Seurat)
library(scRNAtoolVis)
#Visualization scRNAtoolVish
featurePlot(Moss,genes = genes,ncol = 6)




Moss@meta.data$orig.ident <- NULL
Moss@meta.data$orig.ident.x <- NULL
Moss@meta.data$cell_type <- NULL
Moss@meta.data$orig.ident.x <- NULL

