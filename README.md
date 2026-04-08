🧬 Embryo Reference Mapping & In Silico Perturbation Framework**

🚀 Workflow

····································································

Step 1 — Build Reference Atlas

You can either build the reference from scratch or download the precomputed version:
👉 https://zenodo.org/records/19476835

If building from scratch, please refer to:

integrate.R

·Contains step-by-step annotated code

·Covers preprocessing, integration, and dimensional reduction

····································································

Step 1.5 — Visualize the reference atlas

Please see visualization.R for example code to load the prebuilt reference object and generate basic UMAP and marker-expression plots. you can skip if you use your own reference map.

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
