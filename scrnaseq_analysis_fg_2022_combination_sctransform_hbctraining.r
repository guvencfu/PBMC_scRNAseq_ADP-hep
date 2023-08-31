#Single cell RNA sequencing data analysis
#Investigator: Furkan Guvenc
#Date: 2022

Dataset integration and analysis was done based on the vignette provided in https://satijalab.org/seurat/archive/v3.1/immune_alignment.html
https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
https://github.com/hbctraining/scRNA-seq/tree/master/lessons




```{r}
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggdark)
library(tidyverse)

```


#if the data has been analyzed already, start from here

```{r}
seurat_integrated <- readRDS(file = "integrated_seurat.rds")
```


#If the data is not analyzed, start from scratch here;
#load the 10X Cellranger output data
```{r}
#untreated data
unt <- Read10X(data.dir = "pbmc_untreated_6h/count/sample_feature_bc_matrix/")

#ADP-heptose treated data
ah <- Read10X(data.dir = "pbmc_adphep_6h/count/sample_feature_bc_matrix/")

#Process the datasets prior to integration

#untreated samples
pbmc.unt <- CreateSeuratObject(counts = unt$`Gene Expression`, project = "unt", min.cells = 3, min.features = 200)

pbmc.unt$stim <- "untreated"
# 
# pbmc.unt[["percent.mt"]] <- PercentageFeatureSet(pbmc.unt, pattern = "^MT-")
# 
# pbmc.unt <- subset(pbmc.unt, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
# 
# pbmc.unt <- NormalizeData(pbmc.unt, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# pbmc.unt <- FindVariableFeatures(pbmc.unt, selection.method = "vst", nfeatures = 2000)

#---

pbmc.ah <- CreateSeuratObject(counts = ah$`Gene Expression`, project = "ah", min.cells = 3, min.features = 200)

pbmc.ah$stim <- "adp-hep"
# 
# pbmc.ah[["percent.mt"]] <- PercentageFeatureSet(pbmc.ah, pattern = "^MT-")
# 
# pbmc.ah <- subset(pbmc.ah, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 10)
# 
# pbmc.ah <- NormalizeData(pbmc.ah, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# pbmc.ah <- FindVariableFeatures(pbmc.ah, selection.method = "vst", nfeatures = 2000)

```


#Merge data for quality control. Merging makes QC of the datasets easier.
```{r}
merged_seurat <- merge(x = pbmc.unt, 
                       y = pbmc.ah, 
                       add.cell.id = c("unt", "ah"))


View(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)


# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^unt_"))] <- "unt"
metadata$sample[which(str_detect(metadata$cells, "^ah_"))] <- "ah"

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
                           
# Create .RData object to load at any time
save(merged_seurat, file="merged_filtered_seurat.rds")

# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")


# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)


# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)


#based on the above graphs, remove if nUMI >= 500, ngene >= 250, log10genes >0.8 and mitoratio <0.2;

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

save(filtered_seurat, file="seurat_filtered.rds")

```


Check cell cycles to see if it is a source of variation
```{r}

seurat_phase <- NormalizeData(filtered_seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
#The cell cycle expression matrix and how to use it is from https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html

exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
    as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)    


# Identify the most variable genes in the cell cycle phase dataset
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

#cell cycle is not a major source of variation in the dataset, so it doesnt need to be regressed out.

```


#perform cell cycle scoring and apply scTransform on all samples

```{r}

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("unt", "ah")]

for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
    split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.genes, s.features=s.genes)
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
    }


```



#Perform integration
```{r}

integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)


# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Save integrated seurat object
saveRDS(seurat_integrated, "integrated_seurat.rds")


# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  


# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
			     reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)    


```

#identidfy significant PCs

```{r}

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)


# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)

```


```{r}

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:40)
                                
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))



# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


#Trying a lower resolution
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3)


#continue with a resolution of 0.4


saveRDS(seurat_integrated, file = "integrated_seurat.rds")

```

Segragation of clusters by sample
```{r}

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
        dplyr::count(ident, orig.ident) %>%
        tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()

```


#segragation of samples by cell cycle
```{r}

# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

```

Segragation of samples by uninteresting variation
```{r}

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

```

#exploration of the PCs that drive different clusters

```{r}

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
            "ident",
            "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)


# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]

# Plotting a UMAP plot for each of the PCs
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
  
# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
        ggplot(pc_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=pc), 
                           alpha = 0.7) +
                scale_color_gradient(guide = FALSE, 
                                     low = "grey90", 
                                     high = "blue")  +
                geom_text(data=umap_label, 
                          aes(label=ident, x, y)) +
                ggtitle(pc)
}) %>% 
        plot_grid(plotlist = .)


# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

```



```{r}
seurat_integrated <- readRDS(file = "integrated_seurat.rds")
```



#Identify conserved cell type markers

This function (FindConservedMarkers) will allow analysis of genes that are conserved across the clusters irrespective of the treatment condition

```{r}

#to use the FindConservedMarkers function, install these packages'
#uncomment to install

#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
#BiocManager::install("limma")

library(multtest)
library(metap)
library(limma)



DefaultAssay(seurat_integrated) <- "RNA"

# Create an empty list to store the clusters
clusters_list <- list()

# Loop through the different ident.1 values and find the conserved markers
for (i in 0:14) {
  clusters_list[[paste0("cluster", i)]] <- FindConservedMarkers(seurat_integrated, 
                                                                 ident.1 = i, 
                                                                 grouping.var = "stim", 
                                                                 verbose = FALSE)
}



#Top 20 genes in each cluster


library(dplyr)


# Create an empty list to store the top 20 genes for each cluster
top20_genes_list <- list()

# Loop through the clusters_list and extract the top 20 genes for each cluster
for (i in 0:14) {
  # Get the cluster object from clusters_list
  cluster <- clusters_list[[paste0("cluster", i)]]
  
  # Sort the cluster by untreated_avg_log2FC in descending order and extract top 20 genes
  top20_genes <- cluster %>%
    arrange(desc(untreated_avg_log2FC)) %>%
    rownames() %>%
    as.vector() %>%
    head(20)
  
  # Store the top 20 genes in the top20_genes_list
  top20_genes_list[[paste0("cluster", i, ".top20")]] <- top20_genes
}




#export the files

write.csv(cluster0, "output/combined.data/conserved.cluster.genes/cluster0.csv")
write.csv(cluster1, "output/combined.data/conserved.cluster.genes/cluster1.csv")
write.csv(cluster2, "output/combined.data/conserved.cluster.genes/cluster2.csv")
write.csv(cluster3, "output/combined.data/conserved.cluster.genes/cluster3.csv")
write.csv(cluster4, "output/combined.data/conserved.cluster.genes/cluster4.csv")
write.csv(cluster5, "output/combined.data/conserved.cluster.genes/cluster5.csv")
write.csv(cluster6, "output/combined.data/conserved.cluster.genes/cluster6.csv")
write.csv(cluster7, "output/combined.data/conserved.cluster.genes/cluster7.csv")
write.csv(cluster8, "output/combined.data/conserved.cluster.genes/cluster8.csv")
write.csv(cluster9, "output/combined.data/conserved.cluster.genes/cluster9.csv")
write.csv(cluster10, "output/combined.data/conserved.cluster.genes/cluster10.csv")
write.csv(cluster11, "output/combined.data/conserved.cluster.genes/cluster11.csv")
write.csv(cluster12, "output/combined.data/conserved.cluster.genes/cluster12.csv")
write.csv(cluster13, "output/combined.data/conserved.cluster.genes/cluster13.csv")

```



#finding differentially expressed features (cluster biomarkers)

This is done to define the markers that define the clusters via differential expression.

```{r}

markers.to.plot <- c(cluster0.top20, cluster1.top20, cluster2.top20, cluster3.top20, cluster4.top20,
                     cluster5.top20, cluster6.top20, cluster7.top20, cluster8.top20, cluster9.top20,
                     cluster10.top20, cluster11.top20, cluster12.top20, cluster13.top20)

DotPlot(seurat_integrated, features = rev(top20_genes_list$cluster0.top20), cols = c("blue", "red"), dot.scale = 1, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
DotPlot(immune.combined, features = rev(cluster0.top20), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()



```


#save the project

```{r}
saveRDS(immune.combined, file = "../per_sample_outs/pbmc.combined.rds")

```




#Determine cluster identities

```{r}

DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

#CD14 monocyte markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = T)

#FCGR3A+ markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#Macrophages
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#conventional DCs

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#plasmacytoid dendritic cells

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
#B-cells

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD79A", "MS4A1", "SIGLEC2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#T-cells


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D", "CD8A", "IL7R", "SELL", "CCR7"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#T-cell activation markers

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD3D", "CD8A", "CD69", "IL2", "EGR1", "FOS"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

#platellets

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PPBP"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


#marker genes from the previous analysis
FeaturePlot(seurat_integrated, features = c("CD3D", "SELL", "CD8A", "GNLY", "CD79A", "FCGR3A", 
    "CCL2", "PPBP", "CCR7", "CD4", "CD14", "LYZ"), cols = c("grey", "blue"), ncol = 3)


#rename the clusters based on the marker genes

# Define the list with the new cluster names
new_cluster_names <- c("CD14 monocyte", "Memory CD4 T-cell", "Naïve CD4 T-cell",
                       "Naïve CD8 T-cell", "CD14 monocyte", "NK cell", "B-cells",
                       "Memory CD8 T-cell", "Memory CD8 T-cell", "FCGR3A+ monocytes",
                       "Activated CD4 T-cell", "NKT cell", "Conventional DCs", "pDC",
                       "platelets")

# Rename the cluster identities in the seurat_integrated object
seurat_integrated <- RenameIdents(seurat_integrated, 
                                  "0" = new_cluster_names[1],
                                  "1" = new_cluster_names[2],
                                  "2" = new_cluster_names[3],
                                  "3" = new_cluster_names[4],
                                  "4" = new_cluster_names[1],
                                  "5" = new_cluster_names[6],
                                  "6" = new_cluster_names[7],
                                  "7" = new_cluster_names[8],
                                  "8" = new_cluster_names[9],
                                  "9" = new_cluster_names[10],
                                  "10" = new_cluster_names[11],
                                  "11" = new_cluster_names[12],
                                  "12" = new_cluster_names[13],
                                  "13" = new_cluster_names[14],
                                  "14" = new_cluster_names[15])

# Optional: Check the updated identities
table(Idents(seurat_integrated))

saveRDS(seurat_integrated, file = "integrated_seurat.rds")



```

#Perform DEG analysis on the clusters

```{r}


DefaultAssay(seurat_integrated) <- "RNA" #integrated is used for the visualization following integration. Use "RNA" for DE analysis
seurat_integrated$celltype.stim <- paste(Idents(seurat_integrated), seurat_integrated$sample, sep = "_")
seurat_integrated$celltype <- Idents(seurat_integrated)
Idents(seurat_integrated) <- "celltype.stim"


# List of clusters
clusters <- c("CD14 monocyte", "Memory CD4 T-cell", "Naïve CD4 T-cell", 
              "Naïve CD8 T-cell", "NK cell", "B-cells", "Memory CD8 T-cell",
              "FCGR3A+ monocytes", "Activated CD4 T-cell", "NKT cell",
              "Conventional DCs", "pDC", "platelets")
# Create an empty list to store the results
marker_results_list <- list()

# Loop through each cluster to perform FindMarkers and save results to CSV files
for (cluster_name in clusters) {
  # Form the treatment condition names based on the cluster name
  ident_ah <- paste(cluster_name, "_ah", sep = "")
  ident_unt <- paste(cluster_name, "_unt", sep = "")
  
  # Perform FindMarkers
  marker_results <- FindMarkers(seurat_integrated, ident.1 = ident_ah, ident.2 = ident_unt, verbose = TRUE)
  
  # Extract gene names
  gene_names <- rownames(marker_results)
  
  # Add gene names to the marker results data.frame
  marker_results <- cbind(Gene = gene_names, marker_results)
  
  # Save results to CSV file
  file_name <- paste(cluster_name, "_response.csv", sep = "")
  write.csv(marker_results, file = file_name, row.names = FALSE)
}


```




Plotting the data
```{r}

library(ggdark)
library(viridis)
library(ggplot2)


DefaultAssay(seurat_integrated) <- "RNA"

#dimplots for the original clusters


Idents(seurat_integrated) <- "celltype"
p1 <- DimPlot(seurat_integrated, label = T, label.color = "white", label.box = T, repel = T) & dark_theme_gray()

#plot based on stimulation

p2 <- DimPlot(immune.combined, group.by = "stim") & dark_theme_classic()


plot_grid(p1, p2)

#IDO1 comparison in treated and untreated samples

p3 <- FeaturePlot(seurat_integrated, features = "IDO1") & scale_color_viridis(option = "D") & dark_theme_gray()

plot_grid(p1, p3)
#IFNg expression in the cells

p4 <- FeaturePlot(immune.combined, features = "IFNG") & scale_color_viridis(option = "D") & dark_theme_gray()

plot_grid(p1, p4)


#vlnplot of TIFA and ALPK1

VlnPlot(seurat_integrated, features = "TIFA")

VlnPlot(seurat_integrated, features = "ALPK1")

FeaturePlot(seurat_integrated, features = c("TIFA", "ALPK1"), blend = T, pt.size = 2)

```

#various plots
```{r}

library(RColorBrewer)

DefaultAssay(seurat_integrated) <- "RNA"

#dimplot with labels

dimplot.all <- DimPlot(seurat_integrated, reduction = "umap", label = T)

ggsave(filename = "dimplot.all.png", plot = dimplot.all, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 12, height = 7, units = "in")

dimplot.nolegend <- DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 3) + NoLegend()

ggsave(filename = "dimplot.nolegend.png", plot = dimplot.nolegend, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")



#dimplot of the data grouped by stim

umap.stim <- DimPlot(seurat_integrated, reduction = "umap", group.by = "stim")

ggsave(filename = "umap_stim.png", plot = umap.stim, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

#IFNg expression in the dataset

infg.featureplot <- FeaturePlot(seurat_integrated, reduction = "umap", features = "IFNG")

ggsave(filename = "ifng_featureplot.png", plot = infg.featureplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

#CD69 in the population

cd69.featureplot <- FeaturePlot(seurat_integrated, reduction = "umap", features = "CD69")

ggsave(filename = "CD69_featureplot.png", plot = cd69.featureplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

CD69.vlnplot <- VlnPlot(seurat_integrated, features = "CD69", split.by = "stim", split.plot = T)

ggsave(filename = "CD69_vlnplot.png", plot = CD69.vlnplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

cd69.tcells.vlnplot <- CD69.vlnplot <- VlnPlot(seurat_integrated, features = "CD69", split.by = "stim", split.plot = T,
                                               idents = c("Naïve CD4 T-cell", "Memory CD4 T-cell", "Naïve CD8 T-cell", 
                                                          "Memory CD8 T-cell"))



#expression of NK cell inducing cytokines
FeaturePlot(seurat_integrated, features = c("IL2", "IL12A", "IL15", "IL18"))

#IL15 expression in monocytes

VlnPlot(seurat_integrated, features = "IL15", split.by = "stim", idents = c("CD14 monocyte", "FCGR3A+ monocytes", "Conventional DCs"))

library(ggsignif)

VlnPlot(immune.combined, features = "IL15", split.by = "stim", idents = "Monocytes")

#IDO1 plots

ido1.feature <- FeaturePlot(seurat_integrated, features = "IDO1", reduction = "umap", label = T, repel = T)

ggsave(filename = "ido1_featureplot.png", plot = ido1.feature, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

ido1.vlnplot <- VlnPlot(seurat_integrated, features = "IDO1", split.by = "stim", split.plot = T)

ggsave(filename = "ido1_vlnplot.png", plot = ido1.vlnplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")


#dotplot for ifng stimulated genes

dotplot_features <- c("IFNGR1", "GBP1", "IDO1", "NFKBIA","JUNB", "TNFAIP3", "CCL4", "SOCS3")
RidgePlot(immune.combined, features = c("IFNGR1", "GBP1", "IDO1", "NFKBIA","JUNB", "TNFAIP3", "CCL4", "SOCS3"), ncol = 2)

DotPlot(immune.combined, features = c("IFNGR1", "GBP1", "IDO1", "NFKBIA","JUNB", "TNFAIP3", "CCL4", "SOCS3"), split.by = "stim")

#IFNg and CD69 blend plots
FeaturePlot(seurat_integrated, features = c("CD69", "IFNG"), blend = T, pt.size = 1)

#CD69 and CD24 blend plots

FeaturePlot(seurat_integrated, features = c("CD69", "IL2RA"), blend = T, pt.size = 1, blend.threshold = 1)
  

#heatmap of top 20 markers in each cluster

pbmc.markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DefaultAssay(seurat_integrated) <- "RNA"
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

scaled.seurat <- ScaleData(seurat_integrated, verbose = F)

DoHeatmap(scaled.seurat, features = top10$gene) + NoLegend()

#variable features across all clusters heatmap

DefaultAssay(seurat_integrated) <- "integrated"
variable.heatmap <- DoHeatmap(seurat_integrated, features = VariableFeatures(seurat_integrated)[1:75], cells = 1:1000, size = 3, angle = 90) + NoLegend()


#heatmaps of DEG in all clusters

BiocManager::install("dittoSeq")
library(dittoSeq)

genes <- c("")

dittoHeatmap(sce, genes,
    annot.by = c("label", "donor"))

```

