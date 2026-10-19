# This is the example 

# Rebuild the final lineage-constrained Moss reference URD tree from a
# Seurat reference object, then project an independent Seurat query onto it.
#
# Edit only the configuration section before running:
#   Rscript Moss_URD_tree_rebuild_and_query_projection.R
#
# Required reference metadata: dataset, orig.ident, cell_type, lineage_group.
# Data-derived steps: lineage-aware ComBat, diffusion map, dual-root flood
# pseudotime, directed random walks, and terminal-cell similarity.
# The high-level display topology is explicitly lineage constrained and must
# not be described as a completely unsupervised topology discovery.
# Query projection uses Seurat MapQuery, inverse-distance segment voting, and
# interpolation among up to five real reference cells on the winning segment.
#
# Tested with R 4.4.1, Seurat 4.4.0, sva 3.54.0, and URD 1.1.1 at


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
reference_rds <- "data/mouse_ref.rds"
query_rds <- "data/query.rds"
output_dir <- "results"
query_name <- "query"

seed <- 20260902
pca_dims <- 1:40
sampling_cap_first <- 50
sampling_cap_final <- 30
minimum_tip_cells <- 8
datasets_to_stage_calibrate <- c("Nowotschin,2019", "Jiang,2023")
minimum_prediction_score_for_figures <- 0.4

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
options(future.globals.maxSize=32*1024^3)

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(URD)
  library(Matrix)
  library(sva)
  library(data.table)
  library(RANN)
  library(igraph)
  library(ggplot2)
  library(ggnewscale)
})

require_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) stop(object_name, " is missing: ", paste(missing, collapse=", "))
}

sample_cap <- function(meta, cap) {
  unique(unlist(lapply(split(rownames(meta), meta$stratum), function(x) {
    if (length(x)>cap) sample(x, cap) else x
  })))
}

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

# -----------------------------------------------------------------------------
# Prepare Reference
# -----------------------------------------------------------------------------
mouse_ref <- readRDS(reference_rds)
stopifnot(inherits(mouse_ref, "Seurat"))
meta <- mouse_ref[[]]
require_columns(meta, c("dataset", "orig.ident", "cell_type", "lineage_group"), "mouse_ref metadata")

meta$stage_num <- as.numeric(gsub("[^0-9.]", "", as.character(meta$orig.ident)))
meta$stage <- as.character(meta$orig.ident)
if (anyNA(meta$stage_num)) stop("Could not parse all embryonic stages from orig.ident.")
meta$stratum <- paste(meta$dataset, meta$stage, meta$cell_type, sep=" | ")

te_types <- c("TE","mural TE","polar TE","ExE","EPC","EPC.mig","SpT",
              "Pri.pTGC","Sec.pTGC","SynTI.p")
meta$tree_compartment <- ifelse(meta$cell_type %in% te_types, "TE", "ICM")
meta$bio_lineage <- "Other"
meta$bio_lineage[meta$cell_type %in% te_types] <- "Trophoblast"
meta$bio_lineage[meta$cell_type %in% c("PE","VE","ExVE","EmVE","PaE","DVE")] <- "Primitive_endoderm"
meta$bio_lineage[meta$cell_type %in% c("Def.endo","Gut")] <- "Definitive_endoderm"
meta$bio_lineage[meta$cell_type %in% c(
  "Mesoderm","PS","ExM","Mesenchyme","Para.meso","Nas.meso","Caudal.meso",
  "Phary.meso","Allantois","Axial midline","HEPs","Blood.prog","Erythroid",
  "Endothelium","Cardiomyocytes"
)] <- "Mesoderm"
meta$bio_lineage[meta$cell_type %in% c("Sur.ecto","Neuralectoderm","Spinal cord","Brain","Neural Crest")] <- "Ectoderm"
meta$bio_lineage[meta$cell_type=="NMP"] <- "Posterior_progenitor"
meta$bio_lineage[meta$cell_type %in% c("EPI","ICM","Pre_EPI","Caudal.epi")] <- "Epiblast"
meta$bio_lineage[meta$cell_type=="PGC"] <- "PGC"
meta$bio_lineage_model <- factor(meta$bio_lineage)

root_stage <- min(meta$stage_num)
set.seed(seed)
cells_first <- sample_cap(meta, sampling_cap_first)
roots <- rownames(meta)[meta$stage_num==root_stage & meta$cell_type %in% c("TE","ICM")]
rare <- rownames(meta)[meta$cell_type %in% names(which(table(meta$cell_type)<=200))]
cells_first <- sort(unique(c(cells_first, roots, rare)))

set.seed(seed)
meta_first <- meta[cells_first, , drop=FALSE]
cells <- sample_cap(meta_first, sampling_cap_final)
roots <- rownames(meta_first)[meta_first$stage_num==root_stage & meta_first$cell_type %in% c("TE","ICM")]
rare <- rownames(meta_first)[meta_first$cell_type %in% names(which(table(meta_first$cell_type)<=200))]
cells <- sort(unique(c(cells, roots, rare)))
meta <- meta_first[cells, , drop=FALSE]
counts <- GetAssayData(mouse_ref, assay="RNA", slot="counts")[, cells, drop=FALSE]

message("Reference cells retained: ", ncol(counts))

# -----------------------------------------------------------------------------
# Diffusion Map
# -----------------------------------------------------------------------------
seu <- CreateSeuratObject(counts, meta.data=meta, min.cells=3, min.features=100)
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(seu)
hvg <- hvg[!grepl("^(mt-|Mt-|Rpl|Rps|Gm[0-9]+|Hist[0-9]|Hba-|Hbb-)", hvg)]
markers <- c(
  "Cdx2","Eomes","Elf5","Tfap2c","Gata3","Id2","Bmp4","Ascl2","Tpbpa",
  "Pcdh12","Hand1","Plac8","Gcm1","Hbegf","Prl2c2","Prl4a1","Krt19",
  "Lgals3","Pou5f1","Nanog","Sox2","Gata6","Sox17","Foxa2","Cer1",
  "Lefty1","T","Mesp1","Mixl1","Tbx6","Msgn1","Pdgfra","Kdr","Tal1",
  "Etv2","Gata1","Klf1","Nkx2-5","Tbx5","Myl7","Sox1","Sox3","Pax6",
  "Sox10","Tfap2a","Tubb3","Dppa3","Prdm1","Nanos3"
)
hvg <- unique(c(hvg, intersect(markers, rownames(counts))))
counts <- counts[hvg, , drop=FALSE]
meta <- meta[colnames(counts), , drop=FALSE]
rm(seu); gc()

urd <- createURD(count.data=counts, meta=meta, min.cells=3, min.genes=100,
                 min.counts=1, gene.max.cut=1e9, max.genes.in.ram=5000, verbose=TRUE)
urd@group.ids$stage <- urd@meta[rownames(urd@group.ids), "stage"]
urd@group.ids$cell_type <- urd@meta[rownames(urd@group.ids), "cell_type"]
urd@group.ids$dataset <- urd@meta[rownames(urd@group.ids), "dataset"]
urd@group.ids$bio_lineage <- urd@meta[rownames(urd@group.ids), "bio_lineage"]
urd@group.ids$compartment <- urd@meta[rownames(urd@group.ids), "tree_compartment"]
urd@var.genes <- intersect(hvg, rownames(urd@logupx.data))

urd <- batchComBat(urd, batch.column="dataset",
                   adjustment.columns=c("stage_num","bio_lineage_model"),
                   genes.correct=urd@var.genes, verbose=TRUE)
urd@var.genes <- rownames(urd@logupx.data)
profiles <- as.matrix(t(urd@logupx.data))
urd <- urdSubset(urd, cells.keep=rownames(profiles)[!duplicated(profiles)])
rm(profiles, counts); gc()

urd <- calcDM(urd, genes.use=urd@var.genes, knn=100, sigma.use="local",
              dcs.store=40, verbose=TRUE)

# -----------------------------------------------------------------------------
# Pseudotime
# -----------------------------------------------------------------------------
root_te <- rownames(urd@meta)[urd@meta$stage_num==root_stage & urd@meta$cell_type=="TE"]
root_icm <- rownames(urd@meta)[urd@meta$stage_num==root_stage & urd@meta$cell_type=="ICM"]
root_cells <- c(root_te, root_icm)
urd@meta$root_identity <- ifelse(rownames(urd@meta) %in% root_te, "E3.5 TE",
                                 ifelse(rownames(urd@meta) %in% root_icm, "E3.5 ICM", NA))
urd@group.ids$root_identity <- urd@meta[rownames(urd@group.ids), "root_identity"]

set.seed(seed)
floods <- floodPseudotime(urd, root.cells=root_cells, n=50,
                          minimum.cells.flooded=2, verbose=FALSE)
urd <- floodPseudotimeProcess(urd, floods, floods.name="pseudotime")
raw_pt <- setNames(urd@pseudotime$pseudotime, rownames(urd@pseudotime))
pt <- as.data.table(urd@meta[names(raw_pt), ], keep.rownames="cell_id")
pt[, `:=`(stage=as.numeric(stage_num), raw_pseudotime=raw_pt[cell_id])]

calibrate <- intersect(datasets_to_stage_calibrate, unique(pt$dataset))
corrected_pt <- raw_pt
if (length(calibrate)) {
  target <- pt[!dataset %in% calibrate, .(target=median(raw_pseudotime)), by=stage][order(stage)]
  target[, target:=isoreg(stage, target)$yf]
  late_slope <- max(diff(tail(target$target,2))/diff(tail(target$stage,2)), 0.02)
  target_at <- function(x) {
    y <- approx(target$stage, target$target, xout=x, rule=2)$y
    late <- x>max(target$stage)
    y[late] <- tail(target$target,1)+late_slope*(x[late]-max(target$stage))
    y
  }
  observed <- pt[dataset %in% calibrate,
                 .(observed=median(raw_pseudotime)), by=.(dataset,stage)][order(dataset,stage)]
  observed[, correction:=target_at(stage)-observed]
  for (ds in calibrate) {
    z <- observed[dataset==ds]
    correction_function <- if (nrow(z)>1) approxfun(z$stage,z$correction,rule=2) else function(x) rep(z$correction[1],length(x))
    use <- pt$dataset==ds
    corrected_pt[pt$cell_id[use]] <- pt$raw_pseudotime[use]+correction_function(pt$stage[use])
  }
}
corrected_pt[pt$stage==root_stage] <- 0
urd@pseudotime$pseudotime_flood_original <- raw_pt[rownames(urd@pseudotime)]
urd@pseudotime$pseudotime <- pmax(corrected_pt[rownames(urd@pseudotime)],0)

# -----------------------------------------------------------------------------
# Urd Topology
# -----------------------------------------------------------------------------
tip_definitions <- list(
  "Primary pTGC"="Pri.pTGC", "Secondary pTGC"="Sec.pTGC",
  "SpT / migrating EPC"=c("SpT","EPC.mig"), "SynTI progenitor"="SynTI.p",
  "Extraembryonic VE"="ExVE", "Embryonic VE"="EmVE", "Parietal endoderm"="PaE",
  "Definitive endoderm"="Def.endo", "Gut endoderm"="Gut",
  "Paraxial mesoderm"="Para.meso", "Nascent mesoderm"="Nas.meso",
  "Pharyngeal mesoderm"="Phary.meso", "Caudal mesoderm"="Caudal.meso",
  "Extraembryonic mesoderm"="ExM", "Mesenchyme"="Mesenchyme", "Allantois"="Allantois",
  "Axial midline"="Axial midline", "Cardiomyocytes"="Cardiomyocytes",
  "Endothelium"="Endothelium", "Erythroid"="Erythroid", "Brain"="Brain",
  "Neural ectoderm"="Neuralectoderm", "Spinal cord"="Spinal cord",
  "Neural crest"="Neural Crest", "Surface ectoderm"="Sur.ecto",
  "Primordial germ cells"="PGC"
)

tip_group <- setNames(rep(NA_character_,nrow(urd@meta)),rownames(urd@meta))
for (i in seq_along(tip_definitions)) {
  z <- rownames(urd@meta)[urd@meta$cell_type %in% tip_definitions[[i]]]
  tab <- table(urd@meta[z,"stage_num"])
  eligible <- as.numeric(names(tab)[tab>=minimum_tip_cells])
  tip_stage <- if (length(eligible)) max(eligible) else max(as.numeric(names(tab)))
  z <- z[urd@meta[z,"stage_num"]==tip_stage]
  tip_group[z] <- as.character(i)
}
urd@group.ids$tip.clusters <- tip_group[rownames(urd@group.ids)]

logistic <- pseudotimeDetermineLogistic(urd,"pseudotime",optimal.cells.forward=30,
                                        max.cells.back=60,do.plot=FALSE)
transition_matrix <- as.matrix(pseudotimeWeightTransitionMatrix(
  urd,"pseudotime",logistic.params=logistic,verbose=TRUE))
set.seed(seed)
walks <- simulateRandomWalksFromTips(
  urd,tip.group.id="tip.clusters",root.cells=root_cells,
  transition.matrix=transition_matrix,n.per.tip=5000,
  root.visits=1,max.steps=6000,verbose=TRUE)
rm(transition_matrix); gc()
urd <- processRandomWalksFromTips(urd,walks,n.subsample=10,verbose=TRUE)

tree <- loadTipCells(urd,"tip.clusters")
tree <- buildTree(
  tree,pseudotime="pseudotime",tips.use=1:26,
  divergence.method="preference",weighted.fusion=TRUE,
  cells.per.pseudotime.bin=18,bins.per.pseudotime.window=5,
  minimum.visits=8,visit.threshold=.7,save.all.breakpoint.info=TRUE,
  p.thresh=1e-4,min.cells.per.segment=8,min.pseudotime.per.segment=.005,
  dendro.node.size=90,dendro.cell.jitter=.16,dendro.cell.dist.to.tree=.05,
  verbose=TRUE)
tree <- nameSegments(tree,segments=as.character(1:26),
                     segment.names=names(tip_definitions),short.names=names(tip_definitions))
tree@tree$unconstrained.segment.joins <- tree@tree$segment.joins

# -----------------------------------------------------------------------------
# Constrained Layout
# -----------------------------------------------------------------------------
tip_names <- names(tip_definitions)
tip_types <- unname(tip_definitions)
tip_lineage <- c(rep("Trophoblast",4),rep("Endoderm",5),rep("Mesoderm",11),rep("Ectoderm",5),"PGC")
tip_y <- c(34,30,24,20, 12,8,3,-3,-7, -14,-22,-26,-19,-30,-34,-38,-10,-43,-47,-51,
           -58,-62,-70,-66,-74,-55)

pt_tree <- setNames(tree@pseudotime$pseudotime,rownames(tree@pseudotime))
tip_cells <- vector("list",26); tip_end <- numeric(26)
for (i in 1:26) {
  z <- rownames(tree@meta)[tree@meta$cell_type %in% tip_types[[i]]]
  tab <- table(tree@meta[z,"stage_num"])
  eligible <- as.numeric(names(tab)[tab>=minimum_tip_cells])
  tip_stage <- if (length(eligible)) max(eligible) else max(as.numeric(names(tab)))
  z <- z[tree@meta[z,"stage_num"]==tip_stage]
  tip_cells[[i]] <- z
  tip_end[i] <- as.numeric(quantile(pt_tree[z],.97,na.rm=TRUE))
}

nodes <- data.table(
  node=c("TE_root","TE_line","pTGC_pair","SpT_pair","ICM_root","ICM_split","Primitive","VE_pair",
         "Epiblast","Germ_layer","Defendo_pair","Mesoderm","Axial_mesoderm","Para_axial_pair",
         "Nascent_caudal_pair","Extra_mesoderm","ExM_mesen_pair","Cardio_hemato","Endo_ery_pair",
         "Ectoderm","Neural","Brain_neural_pair","Surface_spinal_pair",paste0("tip",1:26)),
  pseudotime=c(0,.28,.52,.48,0,.16,.26,.52,.36,.43,.50,.44,.48,.56,.53,.49,.57,.52,.58,
               .47,.53,.61,.50,tip_end)
)
edges <- data.table(
  from=c("TE_root","TE_line","TE_line","pTGC_pair","pTGC_pair","SpT_pair","SpT_pair",
         "ICM_root","ICM_split","ICM_split","Primitive","Primitive","VE_pair","VE_pair",
         "Epiblast","Germ_layer","Germ_layer","Germ_layer","Germ_layer","Defendo_pair","Defendo_pair",
         "Mesoderm","Mesoderm","Mesoderm","Axial_mesoderm","Axial_mesoderm","Axial_mesoderm",
         "Para_axial_pair","Para_axial_pair","Nascent_caudal_pair","Nascent_caudal_pair",
         "Extra_mesoderm","Extra_mesoderm","ExM_mesen_pair","ExM_mesen_pair","Cardio_hemato",
         "Cardio_hemato","Endo_ery_pair","Endo_ery_pair","Ectoderm","Ectoderm","Neural","Neural",
         "Brain_neural_pair","Brain_neural_pair","Surface_spinal_pair","Surface_spinal_pair"),
  to=c("TE_line","pTGC_pair","SpT_pair","tip1","tip2","tip3","tip4",
       "ICM_split","Primitive","Epiblast","tip7","VE_pair","tip5","tip6",
       "Germ_layer","Defendo_pair","Mesoderm","Ectoderm","tip26","tip8","tip9",
       "Axial_mesoderm","Extra_mesoderm","Cardio_hemato","Para_axial_pair","Nascent_caudal_pair","tip12",
       "tip10","tip17","tip11","tip13","ExM_mesen_pair","tip16","tip14","tip15",
       "tip18","Endo_ery_pair","tip19","tip20","Neural","Surface_spinal_pair","Brain_neural_pair","tip24",
       "tip21","tip22","tip23","tip25")
)
edges[,segment:=c(27,28,29,1,2,3,4,30,31,33,7,32,5,6,48,34,35,43,26,8,9,36,39,41,
                  37,38,12,10,17,11,13,40,16,14,15,18,42,19,20,44,46,45,24,21,22,23,25)]

graph <- graph_from_data_frame(edges[,.(from,to)])
desc_tips <- function(node) as.integer(sub("tip","",grep("^tip",as_ids(subcomponent(graph,node,mode="out")),value=TRUE)))
nodes[,y:=vapply(node,function(x){z<-desc_tips(x);if(length(z))mean(tip_y[z])else NA_real_},0.0)]
nodes[node=="TE_root",y:=5]; nodes[node=="ICM_root",y:=-5]
nodes[node=="Ectoderm",y:=-59]; nodes[grep("^tip",node),y:=tip_y]
nodes[,x:=pseudotime/.70*36]
edges[nodes,`:=`(x=i.x,y=i.y,start_pt=i.pseudotime),on=.(from=node)]
edges[nodes,`:=`(xend=i.x,yend=i.y,end_pt=i.pseudotime),on=.(to=node)]

dm <- scale(as.matrix(tree@dm@eigenvectors[,1:20,drop=FALSE]))
tree_cells <- rownames(dm)
centroids <- do.call(rbind,lapply(tip_cells,function(z)colMeans(dm[intersect(z,tree_cells),,drop=FALSE])))
allowed <- list(Trophoblast=1:4,Primitive_endoderm=5:7,Definitive_endoderm=8:9,
                Mesoderm=10:20,Ectoderm=21:25,PGC=26,Posterior_progenitor=10:25,
                Epiblast=8:26,Other=8:26)
assigned_tip <- integer(length(tree_cells)); names(assigned_tip) <- tree_cells
for (i in 1:26) assigned_tip[intersect(tree_cells,rownames(tree@meta)[tree@meta$cell_type %in% tip_types[[i]]])] <- i
for (k in which(assigned_tip==0)) {
  allowed_tips <- allowed[[as.character(tree@meta[tree_cells[k],"bio_lineage"])]]
  if (is.null(allowed_tips)) allowed_tips <- 8:26
  assigned_tip[k] <- allowed_tips[which.min(rowSums(
    (centroids[allowed_tips,,drop=FALSE]-matrix(dm[k,],length(allowed_tips),ncol(dm),byrow=TRUE))^2))]
}

paths <- lapply(1:26,function(i){
  root <- if(i<=4)"TE_root" else "ICM_root"
  vertices <- as_ids(shortest_paths(graph,root,paste0("tip",i),mode="out")$vpath[[1]])
  edges[data.table(from=head(vertices,-1),to=tail(vertices,-1)),on=.(from,to),nomatch=0][order(start_pt)]
})
d <- data.table(cell_id=tree_cells,pseudotime=pt_tree[tree_cells],stage=tree@meta[tree_cells,"stage_num"],
                dataset=tree@meta[tree_cells,"dataset"],cell_type=tree@meta[tree_cells,"cell_type"],
                bio_lineage=tree@meta[tree_cells,"bio_lineage"],assigned_tip=assigned_tip[tree_cells])
d[,`:=`(segment=NA_integer_,plot_x=NA_real_,plot_y=NA_real_,center_x=NA_real_,center_y=NA_real_,progress=NA_real_)]
for (i in 1:26) {
  use <- which(d$assigned_tip==i); path <- paths[[i]]
  for (k in use) {
    edge_index <- which(d$pseudotime[k]<=path$end_pt+1e-10)[1]
    if (is.na(edge_index)) edge_index <- nrow(path)
    q <- pmin(1,pmax(0,(d$pseudotime[k]-path$start_pt[edge_index])/
                           max(path$end_pt[edge_index]-path$start_pt[edge_index],1e-6)))
    d$segment[k] <- path$segment[edge_index]; d$progress[k] <- q
    d$center_x[k] <- path$x[edge_index]+q*(path$xend[edge_index]-path$x[edge_index])
    d$center_y[k] <- path$y[edge_index]+q*(path$yend[edge_index]-path$y[edge_index])
  }
}

d[,raw_width:=dm[cell_id,2]]
d[,z_width:={z<-raw_width-median(raw_width);sc<-c(IQR(z),mad(z),sd(z),1);sc<-sc[is.finite(sc)&sc>1e-6][1];pmin(2.5,pmax(-2.5,z/sc))},by=segment]
d[edges,`:=`(dx=i.xend-i.x,dy=i.yend-i.y),on=.(segment)]
d[,edge_len:=sqrt(dx^2+dy^2)]
d[,width:=.62*pmin(1,.20+progress/.10)]
d[segment<=26,width:=width*(1-.82*(pmax(0,(progress-.25)/.75)^2*(3-2*pmax(0,(progress-.25)/.75))))]
d[segment %in% c(27,30,31,33),width:=pmin(width,.13)]
d[,offset:=z_width*width]
d[,`:=`(plot_x=center_x-dy/edge_len*offset,plot_y=center_y+dx/edge_len*offset)]

skeleton <- edges[,.(segment,x,y,xend,yend,start_pt,end_pt)]
tips <- nodes[match(paste0("tip",1:26),node),.(segment=1:26,plot_x=x,plot_y=y)]
tips[,`:=`(label=tip_names,lineage=tip_lineage)]
roots <- data.table(label=c("E3.5 TE","E3.5 ICM"),plot_x=0,plot_y=c(5,-5))

incoming_segment <- setNames(edges$segment,edges$to)
nonroot_edges <- edges[!segment %in% c(27,30)]
joins <- rbind(data.table(parent=47L,child=c(27L,30L),pseudotime=0),
               nonroot_edges[,.(parent=as.integer(incoming_segment[from]),child=segment,pseudotime=start_pt)])

tree@tree$lineage.constrained.segment.joins <- as.data.frame(joins)
tree@tree$lineage.constrained.nodes <- as.data.frame(nodes)
tree@tree$lineage.constrained.edges <- as.data.frame(edges)
tree@tree$walks.fan.rotated.layout <- as.data.frame(d)
tree@tree$fan.rotated.skeleton <- as.data.frame(skeleton)
tree@tree$fan.tip.labels <- as.data.frame(tips)
tree@tree$active_display_layout <- "walks.fan.rotated.layout (lineage-constrained)"
tree@tree$lineage_constraint_method <- list(
  high_level_topology="Known embryonic lineage constraints; no cross-lineage terminal fusion",
  cell_assignment="Nearest terminal centroid in lineage-balanced diffusion coordinates, restricted to biological lineage",
  pseudotime="Dual-root URD flood pseudotime with optional late-dataset stage calibration"
)
tree@group.ids$segment_lineage_constrained <- as.character(setNames(d$segment,d$cell_id)[rownames(tree@group.ids)])

# -----------------------------------------------------------------------------
# Save Tree
# -----------------------------------------------------------------------------
stage_levels <- c(3.5,4.5,4.75,5,5.25,5.5,6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75)
stage_colors <- c("3.5"="#ef4b3f","4.5"="#f46d43","4.75"="#f88d52","5"="#fdae61",
                  "5.25"="#fed481","5.5"="#fee8a4","6.25"="#e6f598","6.5"="#c4e687",
                  "6.75"="#a6d96a","7"="#7fca64","7.25"="#50b76b","7.5"="#2ca25f",
                  "7.75"="#1b9e77","8"="#238b8d","8.25"="#2b8cbe","8.5"="#2878b5","8.75"="#345995")
lineage_colors <- c(Trophoblast="#ef4b3f",Endoderm="#169c68",Mesoderm="#237dbd",
                    Ectoderm="#7650b3",PGC="#d43f85")
d[,stage_f:=factor(stage,levels=stage_levels)]

tree_plot <- ggplot()+
  geom_segment(data=skeleton,aes(x,y,xend=xend,yend=yend),color="#65757c",linewidth=.68,alpha=.72)+
  geom_point(data=d,aes(plot_x,plot_y,color=stage_f),size=.78,alpha=.78)+
  geom_segment(data=skeleton[segment %in% c(27,30,31,33)],aes(x,y,xend=xend,yend=yend),
               color="#65757c",linewidth=.34,alpha=.82)+
  scale_color_manual(values=stage_colors,drop=FALSE,name="Embryonic stage")+
  ggnewscale::new_scale_color()+
  geom_point(data=tips,aes(plot_x,plot_y,color=lineage),size=2.5)+
  geom_text(data=tips,aes(plot_x,plot_y,label=label,color=lineage),nudge_x=.45,hjust=0,size=4)+
  scale_color_manual(values=lineage_colors,guide="none")+
  geom_label(data=roots,aes(plot_x,plot_y,label=label),nudge_x=-.55,hjust=1,size=5.2,
             fontface="bold",label.size=0,fill="white")+
  coord_cartesian(xlim=c(-3,43),ylim=c(-85,38),clip="off")+
  labs(title="Moss reference: lineage-constrained URD tree")+
  theme_void(base_size=13)+
  theme(plot.title=element_text(face="bold",hjust=.5),legend.position="bottom",
        plot.margin=margin(15,80,12,75))

tree_prefix <- file.path(output_dir,"Moss_reference_lineage_constrained_URD")
saveRDS(tree,paste0(tree_prefix,".rds"),compress="gzip")
fwrite(d,file.path(output_dir,"lineage_constrained_cell_layout.csv.gz"))
fwrite(skeleton,file.path(output_dir,"lineage_constrained_skeleton.csv"))
fwrite(joins,file.path(output_dir,"lineage_constrained_segment_joins.csv"))
ggsave(paste0(tree_prefix,".pdf"),tree_plot,width=11,height=16,device=pdf_device,bg="white")
ggsave(paste0(tree_prefix,".png"),tree_plot,width=11,height=16,dpi=320,bg="white")
tree_plot
