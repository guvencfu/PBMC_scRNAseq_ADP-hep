#gene set enrichment analysis

#based on the tutorials provided in https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# BiocManager::install("clusterProfiler")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(scales)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


#prepare inputs

#load in the files

# Load required library for file manipulation
library(readr)

# Set the path to the folder containing the CSV files
folder_path <- "output/SCTransform.manual.annotation.output/DEG.genes.list/"

# Get a list of all the CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Loop through each CSV file and read it as a data frame
for (csv_file in csv_files) {
  # Extract the filename without the extension to use as the object name
  file_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Replace spaces with underscores in the filename
  file_name <- gsub(" ", "_", file_name)
  
  # Read the CSV file and assign it to a data frame with the modified filename as the object name
  assign(file_name, read_csv(csv_file))
}


#GSE for B-cells
{
  #prepare input
  
  b.cell.genelist <- B_cells_response$avg_log2FC
  
  #name the vector
  
  names(b.cell.genelist) <- B_cells_response$Gene
  
  #omit NA values
  
  b.cell.genelist <- na.omit(b.cell.genelist)
  
  #sort the list in decreasing order
  
  b.cell.genelist <- sort(b.cell.genelist, decreasing = T)
  
  #gse analysis
  gse.b.cell <- gseGO(geneList=b.cell.genelist, 
               ont ="BP", 
               keyType = "SYMBOL",
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = F, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
  
  #plot the data
  b.cell.dotplot <- dotplot(gse.b.cell, showCategory=10, split = ".sign", title = "B-cells, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "b_cell_dotplot.png", plot = b.cell.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cnetplot(gse.b.cell, categorySize="pvalue", foldChange=b.cell.genelist, showCategory = 10)
  
  bcell.p1 <- cnetplot(gse.b.cell, categorySize="pvalue", foldChange=b.cell.genelist, showCategory = 5)
  
  #install.packages(scales)
  #library(scales)
  
  min.value <- floor(min(bcell.p1$data$color, na.rm = T))
  max.value <- ceiling(max(bcell.p1$data$color, na.rm = T))
  
  bcell.p1 + scale_color_gradientn(name = "fold change",
                             colours = c("blue","white","red"),
                             values = rescale(c(min.value, 0, max.value)),
                             limits= c(min.value, max.value), 
                             breaks=c(min.value , 0, max.value) )
}

#GSE for CD4 naive T-cells
{
  #prepare input
  
  cd4.cell.genelist <- `Naive_CD4_T-cell_response`$avg_log2FC
  
  #name the vector
  
  names(cd4.cell.genelist) <- `Naive_CD4_T-cell_response`$Gene  
  
  #omit NA values
  
  cd4.cell.genelist <- na.omit(cd4.cell.genelist)
  
  #sort the list in decreasing order
  
  cd4.cell.genelist <- sort(cd4.cell.genelist, decreasing = T)
  
  #gse analysis
  gse.cd4 <- gseGO(geneList=cd4.cell.genelist, 
                   ont ="BP", 
                   keyType = "SYMBOL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = F, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")
  
  #plot the data
  cd4.naive.dotplot <- dotplot(gse.cd4, showCategory=8, split = ".sign", title = "CD4 Naive T-cells, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd4_naive_dotplot.png", plot = cd4.naive.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cnetplot(gse.cd4, categorySize="pvalue", foldChange=cd4.cell.genelist, showCategory = 10)

}

#GSE for CD4 Memory T-cells
{
  #prepare input
  
  cd4.memory.genelist <- `Memory_CD4_T-cell_response`$avg_log2FC
  
  #name the vector
  
  names(cd4.memory.genelist) <- `Memory_CD4_T-cell_response`$Gene
  
  #omit NA values
  
  cd4.memory.genelist <- na.omit(cd4.memory.genelist)
  
  #sort the list in decreasing order
  
  cd4.memory.genelist <- sort(cd4.memory.genelist, decreasing = T)
  
  #gse analysis
  gse.cd4.memory <- gseGO(geneList=cd4.memory.genelist, 
                   ont ="BP", 
                   keyType = "SYMBOL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = F, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")
  
  #plot the data
  cd4.memory.dotplot <- dotplot(gse.cd4.memory, showCategory=10, split = ".sign", title = "CD4 Memory T-cells, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd4_memory_dotplot.png", plot = cd4.memory.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cd4.cnetplot <- cnetplot(gse.cd4.memory, categorySize="pvalue", foldChange=cd4.memory.genelist, showCategory = 10)
  
  ggsave(filename = "cd4_cnetplot.png", plot = cd4.cnetplot, path = "figures/", dpi = 320, 
         width = 7, height = 7, units = "in")
}


#GSE for Naive CD8 T-cells
{
  #prepare input
  
  cd8.cell.genelist <- `Naive_CD8_T-cell_response`$avg_log2FC
  
  #name the vector
  
  names(cd8.cell.genelist) <- `Naive_CD8_T-cell_response`$Gene
  
  #omit NA values
  
  cd8.cell.genelist <- na.omit(cd8.cell.genelist)
  
  #sort the list in decreasing order
  
  cd8.cell.genelist <- sort(cd8.cell.genelist, decreasing = T)
  
  #gse analysis
  gse.cd8 <- gseGO(geneList=cd8.cell.genelist, 
                   ont ="BP", 
                   keyType = "SYMBOL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = F, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")
  
  #plot the data
  cd8.naive.dotplot<- dotplot(gse.cd8, showCategory=10, split = ".sign", title = "CD8 naive T-cells, untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd8_naive_dotplot.png", plot = cd8.naive.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cnetplot(gse.cd8, categorySize="pvalue", foldChange=cd8.cell.genelist, showCategory = 10)
  
}

#GSE for memory CD8 T-cells

{
  #prepare input
  
  cd8.memory.cell.genelist <- `Memory_CD8_T-cell_response`$avg_log2FC
  
  #name the vector
  
  names(cd8.memory.cell.genelist) <- `Memory_CD8_T-cell_response`$Gene
  
  #omit NA values
  
  cd8.memory.cell.genelist <- na.omit(cd8.memory.cell.genelist)
  
  #sort the list in decreasing order
  
  cd8.memory.cell.genelist <- sort(cd8.memory.cell.genelist, decreasing = T)
  
  #gse analysis
  gse.cd8.memory <- gseGO(geneList=cd8.memory.cell.genelist, 
                   ont ="BP", 
                   keyType = "SYMBOL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = F, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")
  
  #plot the data
  memory_cd8_dotplot <- dotplot(gse.cd8.memory, showCategory=10, split = ".sign", title = "Memory CD8 T-cells, untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd8_memory_dotplot.png", plot = memory_cd8_dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cnetplot(gse.cd8.memory, categorySize="pvalue", foldChange=cd8.cell.genelist, showCategory = 10)
}

#GSE for Activated CD4 T-cells
{
  #prepare input
  
  t.cell.genelist <- `Activated_CD4_T-cell_response`$avg_log2FC
  
  #name the vector
  
  names(t.cell.genelist) <- `Activated_CD4_T-cell_response`$Gene
  
  #omit NA values
  
  t.cell.genelist <- na.omit(t.cell.genelist)
  
  #sort the list in decreasing order
  
  t.cell.genelist <- sort(t.cell.genelist, decreasing = T)
  
  #gse analysis
  gse.t.cell <- gseGO(geneList=t.cell.genelist, 
                   ont ="BP", 
                   keyType = "SYMBOL",
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.001, 
                   verbose = F, 
                   OrgDb = org.Hs.eg.db, 
                   pAdjustMethod = "none")
  
  #plot the data
  active.cd4.dotplot <- dotplot(gse.t.cell, showCategory=10, split = ".sign", title = "T-cells, unt vs ah")  + facet_grid(.~.sign)
  
  ggsave(filename = "cd4_active_dotplot.png", plot = active.cd4.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  cnetplot(gse.t.cell, categorySize="pvalue", foldChange=t.cell.genelist, showCategory = 7)
  
}

#GSE for CD16 (FCGR3A) monocytes
{
  #prepare input
  
  cd16.genelist <- FCGR3A_monocytes_response$avg_log2FC
  
  #name the vector
  
  names(cd16.genelist) <- FCGR3A_monocytes_response$Gene
  
  #omit NA values
  
  cd16.genelist <- na.omit(cd16.genelist)
  
  #sort the list in decreasing order
  
  cd16.genelist <- sort(cd16.genelist, decreasing = T)
  
  #gse analysis
  gse.cd16 <- gseGO(geneList=cd16.genelist, 
                      ont ="BP", 
                      keyType = "SYMBOL",
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = F, 
                      OrgDb = org.Hs.eg.db, 
                      pAdjustMethod = "none")
  
  #dotplot
  cd16.dotplot <- dotplot(gse.cd16, showCategory=10, split = ".sign" ,title = "CD16 Monocytes, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd16_dotplot.png", plot = cd16.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  #category netplot
  cnetplot(gse.cd16, categorySize="pvalue", foldChange=cd16.genelist, showCategory = 10)
  
  #hierchical clustering of the enriched genes
  
  treeplot(gse.cd16, showCategory = 10)
  
  #enrichment map
  
  pairwise.cd16 <- pairwise_termsim(gse.cd16)
  
  emapplot(pairwise.cd16)
  
  upsetplot(pairwise.cd16)
  
}


#GSE for CD14 monocytes
{
  #prepare input
  
  cd14.genelist <- CD14_monocyte_response$avg_log2FC
  
  #name the vector
  
  names(cd14.genelist) <- CD14_monocyte_response$Gene
  
  #omit NA values
  
  cd14.genelist <- na.omit(cd14.genelist)
  
  #sort the list in decreasing order
  
  cd14.genelist <- sort(cd14.genelist, decreasing = T)
  
  #gse analysis
  gse.cd14 <- gseGO(geneList=cd14.genelist, 
                    ont ="BP", 
                    keyType = "SYMBOL",
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = F, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "none")
  
  #dotplot
  cd14.dotplot <- dotplot(gse.cd14, showCategory=8, split = ".sign", title = "CD14 Monocytes, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  ggsave(filename = "cd14_GO_dp.png", plot = cd14.dp, path = "figures/", dpi = 320, 
         width = 7, height = 9, units = "in")
  
  ggsave(filename = "cd14_dotplot.png", plot = cd14.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  #category netplot
  cd14.netplot <- cnetplot(gse.cd14, categorySize="pvalue", foldChange=cd14.genelist, showCategory = 5)
  
  ggsave(filename = "cd14_cnetplot.png", plot = cd14.netplot, path = "figures/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
  #hierchical clustering of the enriched genes
  
  cd14.tree <- treeplot(gse.cd14, showCategory = 10)
  

}

#GSE for Activated CD14 monocytes
{
  #prepare input
  
  cd14.active.genelist <- cd14.active$avg_log2FC
  
  #name the vector
  
  names(cd14.active.genelist) <- cd14.active$X
  
  #omit NA values
  
  cd14.active.genelist <- na.omit(cd14.active.genelist)
  
  #sort the list in decreasing order
  
  cd14.active.genelist <- sort(cd14.active.genelist, decreasing = T)
  
  #gse analysis
  gse.cd14.active <- gseGO(geneList=cd14.active.genelist, 
                    ont ="BP", 
                    keyType = "SYMBOL",
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = F, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "none")
  
  #dotplot
  dotplot(gse.cd14.active, showCategory=10, split = ".sign",title = "Activated CD14 Monocytes, Untreated vs ADP-heptose") + facet_grid(.~.sign)
  
  #category netplot
  cnetplot(gse.cd14.active, categorySize="pvalue", foldChange=cd14.active.genelist, showCategory = 10)
  
  #hierchical clustering of the enriched genes
  
  treeplot(gse.cd16, showCategory = 10)
  
  
}


#GSE for NK cells
{
  #prepare input
  
  nk.genelist <- nk$avg_log2FC
  
  #name the vector
  
  names(nk.genelist) <- nk$X
  
  #omit NA values
  
  nk.genelist <- na.omit(nk.genelist)
  
  #sort the list in decreasing order
  
  nk.genelist <- sort(nk.genelist, decreasing = T)
  
  #gse analysis
  gse.nk <- gseGO(geneList=nk.genelist, 
                    ont ="BP", 
                    keyType = "SYMBOL",
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = F, 
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "none")
  
  #plot the data
  dotplot(gse.nk, showCategory=30, title = "NK cells, Untreated vs ADP-heptose")
  
  cnetplot(gse.nk, categorySize="pvalue", foldChange=nk.genelist, showCategory = 7)
  
}

#GSE for NKT cells
{
  #prepare input
  
  nkt.genelist <- nkt$avg_log2FC
  
  #name the vector
  
  names(nkt.genelist) <- nkt$X
  
  #omit NA values
  
  nkt.genelist <- na.omit(nkt.genelist)
  
  #sort the list in decreasing order
  
  nkt.genelist <- sort(nkt.genelist, decreasing = T)
  
  #gse analysis
  gse.nkt <- gseGO(geneList=nkt.genelist, 
                  ont ="BP", 
                  keyType = "SYMBOL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = F, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "none")
  
  #plot the data
  dotplot(gse.nkt, showCategory=30, title = "NKT cells, Untreated vs ADP-heptose")
  
  cnetplot(gse.nkt, categorySize="pvalue", foldChange=nkt.genelist, showCategory = 7)
}
#GSE for DC
{
  #prepare input
  
  dc.genelist <- Conventional_DCs_response$avg_log2FC
  
  #name the vector
  
  names(dc.genelist) <- Conventional_DCs_response$Gene
  
  #omit NA values
  
  dc.genelist <- na.omit(dc.genelist)
  
  #sort the list in decreasing order
  
  dc.genelist <- sort(dc.genelist, decreasing = T)
  
  #gse analysis
  gse.dc <- gseGO(geneList=dc.genelist, 
                  ont ="ALL", 
                  keyType = "SYMBOL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.001, 
                  verbose = F, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "none")
  
  #plot the data
  dc.dotplot <- dotplot(gse.dc, showCategory=20, title = "Dendritic cells, Untreated vs ADP-heptose", split = ".sign") + facet_grid(.~.sign)
  
  ggsave(filename = "dc_dotplot.png", plot = dc.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 9, units = "in")
}

#GSE for pDC
{
  #prepare input
  
  pdc.genelist <- pDC_response$avg_log2FC
  
  #name the vector
  
  names(pdc.genelist) <- pDC_response$Gene
  
  #omit NA values
  
  pdc.genelist <- na.omit(pdc.genelist)
  
  #sort the list in decreasing order
  
  pdc.genelist <- sort(pdc.genelist, decreasing = T)
  
  #gse analysis
  gse.pdc <- gseGO(geneList=pdc.genelist, 
                  ont ="ALL", 
                  keyType = "SYMBOL",
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.001, 
                  verbose = F, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "none")
  
  #plot the data
  pdc.dotplot <- dotplot(gse.pdc, showCategory=10, title = "pDC cells, Untreated vs ADP-heptose", split = ".sign") + facet_grid(.~.sign)
  
  ggsave(filename = "pdc_dotplot.png", plot = pdc.dotplot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}
