library(Seurat)
library(ggplot2)
library(dplyr)
###set your query seurat obj
query <-  your.seurat.obj

options(future.globals.maxSize = 10 * 1024^3)

###calculate anchors
anchors <- FindTransferAnchors(
  reference = Moss,
  query = query,
  dims = 1:40,
  reference.reduction = "pca")#optional set:,k.anchor =10

##################projection by cell type###################
A_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = Moss,
  refdata = list(celltype = "cell_type"),
  reference.reduction = "pca",
  reduction.model = "umap"
)
#add original annotations
A_mapped$cell_type <- query$orig.ident

##check annotations after projection
table(A_mapped$predicted.celltype)


##Visualization
#get reference UMAP Embeddings
test_umap <- Embeddings(Moss, reduction = "umap")
test_umap_df <- data.frame(UMAP_1 = test_umap[,1], UMAP_2 = test_umap[,2], dataset = "test")
#get projected UMAP Embeddings
a_mapped_umap <- Embeddings(A_mapped, reduction = "ref.umap")
a_mapped_umap_df <- data.frame(UMAP_1 = a_mapped_umap[,1], UMAP_2 = a_mapped_umap[,2], dataset = "A_mapped")

#add cell type annotations, both original and projected
a_mapped_umap_df$celltypeorig <- A_mapped@meta.data$cell_type
a_mapped_umap_df$celltype <- A_mapped@meta.data$predicted.celltype
test_umap_df$celltype <- NA
test_umap_df$celltypeorig <- NA
#combine two df
Moss_reduced_umap_df <- rbind(test_umap_df, a_mapped_umap_df)

#Visualization with ggplot2
# UMAP plot
library(ggplot2)
ggplot(Moss_reduced_umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    data = Moss_reduced_umap_df[Moss_reduced_umap_df$dataset == "test",],
    color = "lightgrey", 
    alpha = 0.4,
    size = 0.5  
  ) +
  geom_point(
    data = Moss_reduced_umap_df[Moss_reduced_umap_df$dataset == "A_mapped",],
    aes(fill = celltype),  
    shape = 21,            
    color = "black",      
    size = 3,              
    stroke = 0.5          
  ) +
  scale_fill_manual(values = color_mapping) +  
  theme_classic() +
  labs(fill = "Cell Type")

# river plot
library(ggplot2)
library(ggalluvial)
celltype <- A_mapped$orig.ident
predicted_celltype <- A_mapped$predicted.celltype
data_summary <- as.data.frame(table(celltype, predicted_celltype))
colnames(data_summary) <- c("CellType", "PredictedCellType", "Count")

ggplot(data_summary, aes(axis1 = CellType, axis2 = PredictedCellType, y = Count)) +
  geom_alluvium(aes(fill = PredictedCellType), width = 0.2) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.2, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CellType", "PredictedCellType"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none")


#########projection by embryonic stage#################

A_mapped <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = Moss,
  refdata = list(cellstage = "orig.ident"),
  reference.reduction = "pca",
  reduction.model = "umap"
)


table(A_mapped$predicted.cellstage)
table(A_mapped$cell_type)

##river plot
library(ggplot2)
library(ggalluvial)

celltype <- A_mapped$cell_type
predicted_cellstage <- A_mapped$predicted.cellstage

data_summary <- as.data.frame(table(celltype, predicted_cellstage))
colnames(data_summary) <- c("CellType", "predicted_cellstage", "Count")

ggplot(data_summary, aes(axis1 = CellType, axis2 = predicted_cellstage, y = Count)) +
  geom_alluvium(aes(fill = predicted_cellstage), width = 0.2) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.2, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CellType", "predicted_cellstage"), expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none")

