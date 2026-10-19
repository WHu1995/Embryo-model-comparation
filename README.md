🧬 Embryo Reference Mapping & In Silico Perturbation Framework**

🚀 Workflow

····································································

Step 1 — Build Reference Atlas

You can either build the reference from scratch or download the precomputed version:
👉 (https://zenodo.org/records/22019550)

If building from scratch, please refer to:

01_build_refmap.R

·Contains step-by-step annotated code

·Covers preprocessing, integration, and dimensional reduction

····································································

Step 1.5 — Visualize the reference atlas

Please see 01.5_Mouseref_visualization.R for example code to load the prebuilt reference object and generate basic UMAP and marker-expression plots. you can skip if you use your own reference map.

····································································

Step 2 — Projection (Label Transfer)

To project your query dataset onto the reference, see:

projection.R

·Implements FindTransferAnchors and MapQuery

·Outputs predicted cell types and developmental stages


····································································

Step 3 — Train GAN Model

To train the generative model, see:

sccGAN_code.ipynb

·Trains GAN on selected embryonic stages

·Learns distribution of cell states



····································································

Step 4 — In Silico Overexpression (OE)

To perform in silico perturbation (e.g., OE simulation), see:

GCN_OE_code.ipynb

·Simulates gene perturbations (e.g., transcription factor OE)

·Can be combined with projection to evaluate lineage shifts


····································································

Session info:
R and R packages:
R version 4.4.1
batchelor: 1.22.0
Seurat: 4.4.0
SeuratObject: 5.2.0
SingleCellExperiment: 1.28.1
scuttle: 1.16.0
Matrix: 1.7-0
matrixStats: 1.4.1
ggplot2: 3.5.1
dplyr: 1.1.4
RColorBrewer: 1.1-3
ggalluvial: 0.12.5
scRNAtoolVis: 0.1.0

Python:
python: 3.10.20
pandas: 2.3.3
numpy: 2.2.6
keras: 3.12.1
torch: 2.10.0+cu128
sklearn: 1.7.2
matplotlib: 3.10.8

