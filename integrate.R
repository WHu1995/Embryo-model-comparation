# integrate.R
# Build an integrated Seurat reference atlas from multiple scRNA-seq datasets

library(Seurat)
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(batchelor)
library(S4Vectors)

###############################################################################
# Step 1. Prepare input datasets
###############################################################################
# Please replace these example object names with your own Seurat objects.
# Each object should contain:
#   1) RNA counts
#   2) metadata columns: orig.ident(embryostage), cell_type, dataset

dataset_list <- list(
  your_dataset1,
  your_dataset2,
  your_dataset3,
  your_dataset4,
  your_dataset5,
  your_dataset6,...
)

dataset_names <- c(
  "your_dataset1",
  "your_dataset2",
  "your_dataset3",
  "your_dataset4",
  "your_dataset5",
  "your_dataset6",...
)

names(dataset_list) <- dataset_names

###############################################################################
# Step 2. Extract count matrices
###############################################################################
count_list <- lapply(dataset_list, function(obj) {
  GetAssayData(obj, assay = "RNA", slot = "counts")
})

###############################################################################
# Step 3. Identify common genes across datasets
###############################################################################
common_genes <- Reduce(intersect, lapply(count_list, rownames))

###############################################################################
# Step 4. Align all datasets to the same gene set
###############################################################################
# For genes missing in an individual dataset, fill with zeros.
align_counts_to_common_genes <- function(count_mat, common_genes) {
  existing_genes <- intersect(common_genes, rownames(count_mat))
  count_filtered <- count_mat[existing_genes, , drop = FALSE]

  missing_genes <- setdiff(common_genes, existing_genes)

  if (length(missing_genes) > 0) {
    zero_block <- Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(count_mat),
      sparse = TRUE,
      dimnames = list(missing_genes, colnames(count_mat))
    )
    count_aligned <- rbind(count_filtered, zero_block)
  } else {
    count_aligned <- count_filtered
  }

  count_aligned <- count_aligned[common_genes, , drop = FALSE]
  return(count_aligned)
}

aligned_count_list <- lapply(count_list, align_counts_to_common_genes, common_genes = common_genes)

###############################################################################
# Step 5. Reconstruct Seurat objects with aligned counts
###############################################################################
dataset_list_aligned <- mapply(
  FUN = function(obj, aligned_counts) {
    CreateSeuratObject(
      counts = aligned_counts,
      assay = "RNA",
      meta.data = obj@meta.data
    )
  },
  obj = dataset_list,
  aligned_counts = aligned_count_list,
  SIMPLIFY = FALSE
)

names(dataset_list_aligned) <- dataset_names

###############################################################################
# Step 6. Extract and merge metadata
###############################################################################
required_metadata <- c("orig.ident", "cell_type", "dataset")

metadata_list <- lapply(dataset_list_aligned, function(obj) {
  meta <- obj@meta.data

  # check required columns
  missing_cols <- setdiff(required_metadata, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing metadata columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  meta[, required_metadata, drop = FALSE]
})

merged_metadata <- do.call(rbind, metadata_list)

###############################################################################
# Step 7. Merge raw count matrices and build merged Seurat object
###############################################################################
final_merged_counts <- do.call(cbind, aligned_count_list)

save(final_merged_counts, file = "merged_raw_count.RData")
save(merged_metadata, file = "merged_metadata.RData")

Ref <- CreateSeuratObject(
  counts = final_merged_counts,
  meta.data = merged_metadata
)

Ref <- NormalizeData(Ref)
Ref <- ScaleData(Ref)
Ref <- FindVariableFeatures(Ref, nfeatures = 2000)

###############################################################################
# Step 8. Select integration features
###############################################################################
dataset_list_subset <- lapply(dataset_list_aligned, function(obj) {
  subset(obj, features = common_genes)
})

features <- SelectIntegrationFeatures(
  object.list = dataset_list_subset,
  nfeatures = 2000
)

###############################################################################
# Step 9. Convert Seurat objects to SingleCellExperiment
###############################################################################
sce_list <- mapply(
  FUN = function(count_mat, meta) {
    sce <- SingleCellExperiment(
      list(counts = count_mat),
      colData = DataFrame(meta)
    )
    sce <- sce[features, ]
    sce <- logNormCounts(sce)
    return(sce)
  },
  count_mat = aligned_count_list,
  meta = metadata_list,
  SIMPLIFY = FALSE
)

###############################################################################
# Step 10. Integrate datasets using fastMNN
###############################################################################
mnn_result <- do.call(
  fastMNN,
  c(
    sce_list,
    list(
      d = 50,
      k = 20
    )
  )
)

save(mnn_result, file = "merged_mnn_result.RData")

###############################################################################
# Step 11. Add metadata and corrected embeddings back to Seurat object
###############################################################################
colData(mnn_result) <- cbind(colData(mnn_result), merged_metadata)

corrected_pca <- reducedDim(mnn_result, "corrected")

Ref[["pca"]] <- CreateDimReducObject(
  embeddings = corrected_pca,
  key = "PC_",
  assay = "RNA"
)

###############################################################################
# Step 12. Run UMAP on corrected space
###############################################################################
Ref <- RunUMAP(
  Ref,
  reduction = "pca",
  dims = 1:40,
  return.model = TRUE,
  seed.use = 666
)

###############################################################################
# Step 13. Save final reference object
###############################################################################
saveRDS(Ref, file = "Moss.rds")
