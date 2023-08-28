# Single cell RNA sequencing of PBMCs treated with ADP-heptose
Code I have used to analyze PBMC stimulation data for my Ph.D. The single-cell RNA sequencing data generated using the 10X platform was analyzed using Seurat packages and I have
utilized the tutorial provided by HBCTraining modules.

- Batched sequencing data that was deconvoluted and loaded onto R to be analyzed using Seurat.
- Data was processed and quality controlled to remove doublets, dying cells, low-quality droplets.
- Cell cycle was assessed to determine its contribution to variability.
- Normalized data using SCtransform and merged treated and untreated conditions.
- Determined clusters in the dataset and determined cluster identities manually.
- Performed differential gene expression analysis in each cluster, treated vs. untreated.
- Using DEG dataset, plotted data using Volcano plots and analyzed DEG dataset using GO-term enrichment to drive hypothesis generation.
