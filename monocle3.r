#trajectory inference of different cell populations

#load packages

library(Seurat)
library(tidyverse)
library(monocle3)
library(SeuratWrappers)


#load in the data

seurat_integrated <- readRDS(file = "integrated_seurat.rds")

#determine cell identities
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()


#monocle on all of the PBMC population


seurat_integrated@active.assay = 'RNA'

pbmc_monocle <- as.cell_data_set(seurat_integrated)
pbmc_monocle <- estimate_size_factors(pbmc_monocle)
pbmc_monocle@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(seurat_integrated[["RNA"]])
pbmc_monocle <- cluster_cells(cds = pbmc_monocle, reduction_method = "UMAP")
pbmc_monocle <- learn_graph(pbmc_monocle, use_partition = TRUE)
pbmc_monocle <- order_cells(pbmc_monocle,reduction_method = "UMAP")
plot_cells(pbmc_monocle,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           cell_size = 1,
           show_trajectory_graph = TRUE)

#genes that contribute to the trajectories

pbmc_pr_test_res <- graph_test(pbmc_monocle, neighbor_graph="principal_graph", cores=4)

#export the Moran's test as csv to better wrangle the data

write.csv(pbmc_pr_test_res, "pbmc.morans.csv")

pr_deg_ids <- row.names(subset(pbmc_pr_test_res, q_value < 0.05))



plot_cells(pbmc_monocle, genes=c(
  "GP9", "GP1BB", "EPB41L3", "KYNU", "IFI30", "CYP1B1", "SLC7A11", "FCER1G", 
  "MAFB", "HLA-DRA", "TYROBP", "CAVIN2", "NKG7", "C15orf48", "GNLY", 
  "CST3", "PF4V1", "GZMB", "C5AR1", "IDO1"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
#---





#---

#CD14 monocyte trajectory inference

#subset CD14 monocytes

cd14_mono_seurat <- subset(seurat_integrated, idents = "CD14 monocyte")

DimPlot(cd14_mono_seurat, reduction = "umap", pt.size = 2)

#perform reclustering

#integrated analysis since working with the clusters

DefaultAssay(cd14_mono_seurat) <- "integrated"

# Run the standard workflow for visualization and clustering
cd14_mono_seurat <- ScaleData(cd14_mono_seurat, verbose = FALSE)
cd14_mono_seurat <- RunPCA(cd14_mono_seurat, npcs = 30, verbose = FALSE)
# UMAP and Clustering
cd14_mono_seurat <- RunUMAP(cd14_mono_seurat, reduction = "pca", dims = 1:10)
cd14_mono_seurat <- FindNeighbors(cd14_mono_seurat, reduction = "pca", dims = 1:20)
cd14_mono_seurat <- FindClusters(cd14_mono_seurat, resolution = 0.5)

DimPlot(cd14_mono_seurat, reduction ="umap",label = F, pt.size = 2)

DimPlot(cd14_mono_seurat, reduction ="umap",label = F, pt.size = 2, group.by = "stim")

#convert re-clustered CD14 monocytes to monocle obj

cd14_mono_seurat@active.assay = 'RNA'

cd14_monocle <- as.cell_data_set(cd14_mono_seurat)
cd14_monocle <- estimate_size_factors(cd14_monocle)
cd14_monocle@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cd14_mono_seurat[["RNA"]])
cd14_monocle <- cluster_cells(cds = cd14_monocle, reduction_method = "UMAP")
cd14_monocle <- learn_graph(cd14_monocle, use_partition = TRUE)
cd14_monocle <- order_cells(cd14_monocle,reduction_method = "UMAP")
plot_cells(cd14_monocle,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           genes = c("GBP1", "IDO1", "STAT1"),
           cell_size = 2,
           show_trajectory_graph = TRUE)

#expression of top 5 genes within the monocyte cluster

plot_cells(cd14_monocle, genes=c("GBP1", "STAT1", "APOL6", "GBP5", "IRF1", "GBP2", "IDO1", 
                                 "LAP3", "TAP1", "PSMB9", "SOD2", "PARP14", 
                                 "SOCS3", "B2M", "WARS", "RNF213", "GBP4", "PSME2", "PDE4B"),
  show_trajectory_graph=FALSE,
  label_cell_groups=FALSE,
  label_leaves=FALSE)


#genes that contribute to the trajectories

cd14_pr_test_res <- graph_test(cd14_monocle, neighbor_graph="principal_graph", cores=4)

#export the Moran's test as csv to better wrangle the data

write.csv(cd14_pr_test_res, "cd14.morans.csv")

cd14_pr_deg_ids <- row.names(subset(cd14_pr_test_res, q_value < 0.05))



#plot the data based on the Moran's test

plot_cells(cd14_monocle, genes=c(
  "CXCR6", "PPP1R16B", "MYC", "SH2D1B", "FABP5", "SLC25A37",
  "AUTS2", "NFATC2", "ATP6V0B", "PFDN5", "RPL22", "AIF1",
  "PASK", "TMA7", "JAML", "PRKACB", "RPL27", "SH2D1A", "MIF",
  "PRF1", "YWHAB"
),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)


#---

#NK cells

NK_seurat <- subset(seurat_integrated, idents = "NK cell")

#perform reclustering

#integrated analysis since working with the clusters

DefaultAssay(NK_seurat) <- "integrated"

# Run the standard workflow for visualization and clustering
NK_seurat <- ScaleData(NK_seurat, verbose = FALSE)
NK_seurat <- RunPCA(NK_seurat, npcs = 30, verbose = FALSE)
# UMAP and Clustering
NK_seurat <- RunUMAP(NK_seurat, reduction = "pca", dims = 1:10)
NK_seurat <- FindNeighbors(NK_seurat, reduction = "pca", dims = 1:20)
NK_seurat <- FindClusters(NK_seurat, resolution = 0.5)

DimPlot(NK_seurat, reduction ="umap",label = T)

#convert re-clustered CD14 monocytes to monocle obj

NK_seurat@active.assay = "RNA"

NK_monocle <- as.cell_data_set(NK_seurat)
NK_monocle <- estimate_size_factors(NK_monocle)
NK_monocle@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(NK_seurat[["RNA"]])
NK_monocle <- cluster_cells(cds = NK_monocle, reduction_method = "UMAP")
NK_monocle <- learn_graph(NK_monocle, use_partition = TRUE)
NK_monocle <- order_cells(NK_monocle,reduction_method = "UMAP")
plot_cells(NK_monocle,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           genes = c("IFNG", "CD69"),
           cell_size = 2,
           show_trajectory_graph = TRUE)


#CD4 memory T-cells

CD4_mem_seurat <- subset(seurat_integrated, idents = "Memory CD4 T-cell")

#perform reclustering

#integrated analysis since working with the clusters

DefaultAssay(CD4_mem_seurat) <- "integrated"

# Run the standard workflow for visualization and clustering
CD4_mem_seurat <- ScaleData(CD4_mem_seurat, verbose = FALSE)
CD4_mem_seurat <- RunPCA(CD4_mem_seurat, npcs = 30, verbose = FALSE)
# UMAP and Clustering
CD4_mem_seurat <- RunUMAP(CD4_mem_seurat, reduction = "pca", dims = 1:10)
CD4_mem_seurat <- FindNeighbors(CD4_mem_seurat, reduction = "pca", dims = 1:20)
CD4_mem_seurat <- FindClusters(CD4_mem_seurat, resolution = 0.5)

DimPlot(CD4_mem_seurat, reduction ="umap",label = T, pt.size = 2)

DimPlot(CD4_mem_seurat, reduction ="umap",label = F, pt.size = 2, group.by = "stim")

#convert re-clustered CD14 monocytes to monocle obj

CD4_mem_seurat@active.assay = "RNA"

CD4_mem_monocle <- as.cell_data_set(CD4_mem_seurat)
CD4_mem_monocle <- estimate_size_factors(CD4_mem_monocle)
CD4_mem_monocle@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(CD4_mem_seurat[["RNA"]])
CD4_mem_monocle <- cluster_cells(cds = CD4_mem_monocle, reduction_method = "UMAP")
CD4_mem_monocle <- learn_graph(CD4_mem_monocle, use_partition = TRUE)
CD4_mem_monocle <- order_cells(CD4_mem_monocle,reduction_method = "UMAP")
plot_cells(CD4_mem_monocle,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           cell_size = 2,
           show_trajectory_graph = TRUE)

#CD4 memory top 20 DEG plotted with monocle3

plot_cells(CD4_mem_monocle,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           genes = c(
             "TNFRSF4", "DDIT4", "KLF2", "HSP90AB1", "CD52", "CREM",
             "NFKBIA", "TNFRSF18", "SNHG16", "IRF1", "ANXA1", "PSME2",
             "CD3G", "JUNB", "SRGN", "TRAC", "TRAF4", "ENO1", "C17orf49",
             "BHLHE40"
           ),
           cell_size = 1,
           show_trajectory_graph = TRUE)

#genes that contribute to the trajectories

cd4_mem_pr_test_res <- graph_test(CD4_mem_monocle, neighbor_graph="principal_graph", cores=4)

#export the Moran's test as csv to better wrangle the data

write.csv(cd4_mem_pr_test_res, "CD4.mem.morans.csv")

cd4.mem_pr_deg_ids <- row.names(subset(cd4_mem_pr_test_res, q_value < 0.05))



#plot the data based on the Moran's test

plot_cells(CD4_mem_monocle, genes=c(
  "DLGAP5", "UBE2C", "RRM2", "GTSE1", "MT-ATP6", "BIRC5",
  "CCNB2", "MT-ND3", "MKI67", "MT-CYB", "CEP55", "KIF14",
  "HIST1H3G", "SHCBP1", "MT-CO2", "HIST1H3C", "POLQ", "MT-CO3",
  "SPC25", "HJURP"
),
show_trajectory_graph=FALSE,
label_cell_groups=FALSE,
label_leaves=FALSE)
 