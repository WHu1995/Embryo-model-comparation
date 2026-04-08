# visualization.R
# Visualization examples for the integrated reference atlas

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

###############################################################################
# Step 1. Load reference object
###############################################################################
# Replace with the path to your downloaded reference object
Mouseref <- readRDS("Mouseref.rds")

###############################################################################
# Step 2. Define colors
###############################################################################
# Developmental stage colors
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

###############################################################################
# Step 3. Basic Seurat visualization
###############################################################################
# UMAP colored by developmental stage
DimPlot(
  Mouseref,
  reduction = "umap",
  group.by = "orig.ident",
  cols = stage_colors,
  raster = FALSE
)

# UMAP colored by cell type
DimPlot(
  Mouseref,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  raster = FALSE
) + NoLegend()

# UMAP colored by dataset
DimPlot(
  Mouseref,
  reduction = "umap",
  group.by = "dataset",
  raster = FALSE
)

# Marker visualization
FeaturePlot(
  Mouseref,
  features = "Pou5f1",
  raster = FALSE
)

# Check cell type composition
table(Mouseref$cell_type)

###############################################################################
# Step 4. Extract UMAP coordinates for ggplot2-based visualization
###############################################################################
umap_data <- as.data.frame(Embeddings(Mouseref, "umap"))
umap_data$cell_type <- Mouseref$cell_type
umap_data$orig.ident <- Mouseref$orig.ident
umap_data$dataset <- Mouseref$dataset

###############################################################################
# Step 5. Optional downsampling for faster plotting
###############################################################################
set.seed(123)
n_plot <- min(50000, nrow(umap_data))
umap_sub <- umap_data[sample(nrow(umap_data), size = n_plot), ]

###############################################################################
# Step 6. Generate a color palette for cell types
###############################################################################
num_cell_types <- length(unique(umap_data$cell_type))

base_colors <- brewer.pal(8, "Set2")

if (num_cell_types <= length(base_colors)) {
  celltype_colors <- setNames(base_colors[1:num_cell_types], sort(unique(umap_data$cell_type)))
} else {
  extra_colors <- colorRampPalette(
    c("#FF6347", "#4682B4", "#32CD32", "#FFD700", "#8A2BE2",
      "#00CED1", "#DC143C", "#ADFF2F", "#FF4500")
  )(num_cell_types - length(base_colors))
  
  full_palette <- c(base_colors, extra_colors)
  celltype_colors <- setNames(full_palette, sort(unique(umap_data$cell_type)))
}

###############################################################################
# Step 7. ggplot2 visualization
###############################################################################
# UMAP colored by cell type
ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2, fill = cell_type)) +
  geom_point(
    shape = 21,
    size = 1.8,
    color = "white",
    stroke = 0.01,
    alpha = 1
  ) +
  scale_fill_manual(values = celltype_colors) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# UMAP colored by dataset
ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2, fill = dataset)) +
  geom_point(
    shape = 21,
    size = 1.8,
    color = "white",
    stroke = 0.01,
    alpha = 1
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# UMAP colored by developmental stage
ggplot(umap_sub, aes(x = UMAP_1, y = UMAP_2, fill = orig.ident)) +
  geom_point(
    shape = 21,
    size = 1.8,
    color = "white",
    stroke = 0.01,
    alpha = 1
  ) +
  scale_fill_manual(values = stage_colors) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

###############################################################################
# Step 8. Optional: visualize a panel of marker genes
###############################################################################
marker_genes <- c(
  "Sox2", "Gata3", "T", "Nanos3",
  "Pou5f1", "Eomes", "Mesp1", "Dppa3",
  "Nanog", "Elf5", "Tbx6", "Dnd1",
  "Otx2", "Gata4", "Klf1", "Ascl2",
  "Pax6", "Sox17", "Etv2", "Prl3d1",
  "Lhx2", "Foxa2", "Gypa", "Plac8"
)

FeaturePlot(
  Mouseref,
  features = marker_genes,
  raster = FALSE,
  ncol = 6
)
