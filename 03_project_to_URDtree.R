####################################################################################
####Download pre-constructed tree file from https://zenodo.org/records/23484968#####
####################################################################################
tree_file <- "path/to/tree_file"
tree <- readRDS(tree_file)
query_name <- "query"
dims <- 1:40



# 1. Load the final URD tree.
urd_tree <- readRDS(tree_file)
tree_cells <- as.data.table(urd_tree@tree$walks.fan.rotated.layout)
tree_cells <- tree_cells[
  is.finite(plot_x) & is.finite(plot_y) &
    is.finite(pseudotime) & is.finite(segment)
]
tree_cells[, `:=`(
  cell_type=as.character(urd_tree@meta[cell_id,"cell_type"]),
  lineage=as.character(urd_tree@meta[cell_id,"lineage_group"]),
  stage=as.numeric(urd_tree@meta[cell_id,"stage_num"])
)]
skeleton <- as.data.table(urd_tree@tree$fan.rotated.skeleton)
tips <- as.data.table(urd_tree@tree$fan.tip.labels)
roots <- data.table(label=c("E3.5 TE","E3.5 ICM"), plot_x=0, plot_y=c(5,-5))
pt_range <- range(tree_cells$pseudotime,na.rm=TRUE)
tree_cells[,progress_score:=(pseudotime-pt_range[1])/diff(pt_range)]

# Shared tree frame and labels.
tree_base <- ggplot()+
  geom_segment(data=skeleton,aes(x,y,xend=xend,yend=yend),
               color="#65757c",linewidth=.58,alpha=.65)+
  coord_cartesian(xlim=c(-3,43),ylim=c(-85,38),clip="off")+
  theme_void(base_size=13)+
  theme(plot.title=element_text(face="bold",hjust=.5),legend.position="bottom",
        plot.margin=margin(15,75,12,70))
tree_labels <- list(
  geom_point(data=tips,aes(plot_x,plot_y),color="#34454d",size=1.8),
  geom_text(data=tips,aes(plot_x,plot_y,label=label),nudge_x=.45,
            hjust=0,size=3.6,color="#26383f"),
  geom_label(data=roots,aes(plot_x,plot_y,label=label),nudge_x=-.55,
             hjust=1,size=4.8,fontface="bold",label.size=0,fill="white")
)

# 2. Visualize the reference tree.
tree_plot <- tree_base+
  geom_point(data=tree_cells[order(progress_score)],
             aes(plot_x,plot_y,color=progress_score),size=.72,alpha=.78)+
  scale_color_viridis_c(option="plasma",limits=c(0,1),name="Progress score")+
  tree_labels+
  labs(title="Moss reference URD tree")

stage_levels <- c(3.5,4.5,4.75,5,5.25,5.5,6.25,6.5,6.75,
                  7,7.25,7.5,7.75,8,8.25,8.5,8.75)
stage_colors <- c(
  "3.5"="#ef4b3f","4.5"="#f46d43","4.75"="#f88d52","5"="#fdae61",
  "5.25"="#fed481","5.5"="#fee8a4","6.25"="#e6f598","6.5"="#c4e687",
  "6.75"="#a6d96a","7"="#7fca64","7.25"="#50b76b","7.5"="#2ca25f",
  "7.75"="#1b9e77","8"="#238b8d","8.25"="#2b8cbe","8.5"="#2878b5",
  "8.75"="#345995"
)
tree_cells[,stage_f:=factor(stage,levels=stage_levels)]
tree_stage_plot <- tree_base+
  geom_point(data=tree_cells,aes(plot_x,plot_y,color=stage_f),size=.72,alpha=.78)+
  scale_color_manual(values=stage_colors,drop=FALSE,name="Embryonic stage")+
  tree_labels+
  labs(title="Moss reference URD tree")

tree_plot
tree_stage_plot


ggsave(file.path(out_dir,"URD_tree.pdf"),tree_plot,
       width=11,height=16,bg="white")

ggsave(file.path(out_dir,"URD_tree_stage.pdf"),tree_stage_plot,
       width=11,height=16,bg="white")


# 3. Project A_mapped using its ref.pca coordinates.
ref_pca <- Embeddings(Mouseref,"pca")[tree_cells$cell_id,dims,drop=FALSE]
query_pca <- Embeddings(A_mapped,"ref.pca")[,dims,drop=FALSE]
pca_sd <- Stdev(Mouseref,"pca")[dims]
pca_sd[!is.finite(pca_sd)|pca_sd<1e-8] <- 1
ref_pca <- sweep(ref_pca,2,pca_sd,"/")
query_pca <- sweep(query_pca,2,pca_sd,"/")

meta <- A_mapped[[]]
result <- data.table(
  cell_id=rownames(query_pca),
  predicted_celltype=as.character(meta[rownames(query_pca),"predicted.celltype"]),
  prediction_score=as.numeric(meta[rownames(query_pca),"predicted.celltype.score"]),
  predicted_stage=as.character(meta[rownames(query_pca),"predicted.cellstage"]),
  predicted_lineage=as.character(meta[rownames(query_pca),"predicted.lineage"]),
  segment=NA_integer_, segment_confidence=NA_real_,
  tree_x=NA_real_, tree_y=NA_real_, progress_score=NA_real_
)

# First vote for one branch segment, then interpolate among five real reference
# cells on that segment. This keeps query points inside the tree branches.
for (celltype in unique(result$predicted_celltype)) {
  query_index <- which(result$predicted_celltype==celltype)
  reference_index <- which(tree_cells$cell_type==celltype)
  if (!length(reference_index)) {
    query_lineage <- unique(na.omit(result$predicted_lineage[query_index]))
    reference_index <- which(tree_cells$lineage %in% query_lineage)
  }
  if (!length(reference_index)) reference_index <- seq_len(nrow(tree_cells))
  
  k <- min(30L,length(reference_index))
  nearest <- nn2(ref_pca[reference_index,,drop=FALSE],
                 query_pca[query_index,,drop=FALSE],k=k)
  
  for (j in seq_along(query_index)) {
    neighbors <- reference_index[nearest$nn.idx[j,]]
    weights <- 1/(nearest$nn.dists[j,]+1e-8)
    segment_weight <- tapply(weights,tree_cells$segment[neighbors],sum)
    selected_segment <- as.integer(names(segment_weight)[which.max(segment_weight)])
    take <- head(which(tree_cells$segment[neighbors]==selected_segment),5L)
    selected <- neighbors[take]
    selected_weight <- weights[take]/sum(weights[take])
    q <- query_index[j]
    
    result$segment[q] <- selected_segment
    result$segment_confidence[q] <- max(segment_weight)/sum(segment_weight)
    result$tree_x[q] <- weighted.mean(tree_cells$plot_x[selected],selected_weight)
    result$tree_y[q] <- weighted.mean(tree_cells$plot_y[selected],selected_weight)
    raw_pseudotime <- weighted.mean(tree_cells$pseudotime[selected],selected_weight)
    result$progress_score[q] <- (raw_pseudotime-pt_range[1])/diff(pt_range)
  }
}

A_mapped$URD_tree_x <- result$tree_x[match(colnames(A_mapped),result$cell_id)]
A_mapped$URD_tree_y <- result$tree_y[match(colnames(A_mapped),result$cell_id)]
A_mapped$URD_progress_score <- result$progress_score[match(colnames(A_mapped),result$cell_id)]
A_mapped$URD_segment <- result$segment[match(colnames(A_mapped),result$cell_id)]
A_mapped$URD_segment_confidence <- result$segment_confidence[match(colnames(A_mapped),result$cell_id)]

projection_plot <- tree_base+
  geom_point(data=tree_cells,aes(plot_x,plot_y),color="grey82",size=.4,alpha=.4)+
  geom_point(data=result[order(progress_score)],
             aes(tree_x,tree_y,color=progress_score),size=.85,alpha=.75)+
  scale_color_viridis_c(option="plasma",limits=c(0,1),name="Progress score")+
  tree_labels+
  labs(title=paste0(query_name," projected onto the Moss URD tree"))

projection_plot

ggsave("Projected_on_URD_tree.pdf",projection_plot,
       width=11,height=16,bg="white")


