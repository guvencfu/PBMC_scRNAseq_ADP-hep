#Volcano plots for DEG in each cell cluster

##############################

#Blighe, K, S Rana, and M Lewis. 2018. "EnhancedVolcano: 
#Publication-ready volcano plots with enhanced colouring and labeling." 
#https://github.com/kevinblighe/EnhancedVolcano.

##############################

#BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
library(tidyverse)

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



#plot the data

#CD14 Monocytes
{
  CD14 <- CD14 <- data.frame(CD14_monocyte_response[, -1], row.names = CD14_monocyte_response$Gene)
   #subset the dataframe to include avg_log2fc and adj.pval
  CD14 <- CD14[,c(-1, -3, -4)]
  
  #plot the data
  cd14.plot <- EnhancedVolcano(CD14,
                  lab = rownames(CD14),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  pCutoff = 0.05,
                  FCcutoff = 0.25,
                  legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                 'p-value & Log (base 2) FC'),
                  title = "CD14 Monocyte response",
                  subtitle = "Untreated vs ADP-heptose, 6h")
  cd14.plot
  
  ggsave(filename = "cd14_volcano.png", plot = cd14.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")

}


#FCGR3A+ Monocytes
{
  FCGR3A <- data.frame(FCGR3A_monocytes_response[, -1], row.names = FCGR3A_monocytes_response$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  FCGR3A <- FCGR3A[,c(-1, -3, -4)]
  
  #plot the data
  FCGR3A.plot <- EnhancedVolcano(FCGR3A,
                               lab = rownames(FCGR3A),
                               x = 'avg_log2FC',
                               y = 'p_val_adj',
                               pCutoff = 0.05,
                               FCcutoff = 0.25,
                               legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                              'p-value & Log (base 2) FC'),
                               title = "FCGR3A+ Monocyte response",
                               subtitle = "Untreated vs ADP-heptose, 6h")
  FCGR3A.plot
  
  
  ggsave(filename = "FCGR3A_volcano.png", plot = FCGR3A.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
}


#Naive CD4 T-cells
{
  naive.cd4 <- data.frame(`Naïve_CD4_T-cell_response`[, -1], row.names = `Naïve_CD4_T-cell_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  naive.cd4 <- naive.cd4[,c(-1, -3, -4)]
  
  #plot the data
  naive.cd4.plot <- EnhancedVolcano(naive.cd4,
                               lab = rownames(naive.cd4),
                               x = 'avg_log2FC',
                               y = 'p_val_adj',
                              xlim = c(-1, 1.5),
                               pCutoff = 0.05,
                               FCcutoff = 0.25,
                               legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                              'p-value & Log (base 2) FC'),
                               title = "CD4 Naive T-cell response",
                               subtitle = "Untreated vs ADP-heptose, 6h")
  
  naive.cd4.plot
  
  ggsave(filename = "naivecd4_volcano.png", plot = naive.cd4.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}

#Activated CD4 T-cells
{
  active.cd4 <- data.frame(`Activated_CD4_T-cell_response`[, -1], row.names =`Activated_CD4_T-cell_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  active.cd4 <- active.cd4[,c(-1, -3, -4)]
  
  #plot the data
  active.cd4.plot <- EnhancedVolcano(active.cd4,
                                    lab = rownames(active.cd4),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj',
                                    xlim = c(-1, 1.5),
                                    pCutoff = 0.05,
                                    FCcutoff = 0.25,
                                    legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                                   'p-value & Log (base 2) FC'),
                                    title = "Active CD4 T-cell response",
                                    subtitle = "Untreated vs ADP-heptose, 6h")
  
  active.cd4.plot
  
  ggsave(filename = "activecd4_volcano.png", plot = active.cd4.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}



#CD4 memory T-cells
{
  memory.cd4 <- data.frame(`Memory_CD4_T-cell_response`[,-1], row.names = `Memory_CD4_T-cell_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  memory.cd4 <- memory.cd4[,c(-1, -3, -4)]
  
  #plot the data
  cd4.memory.plot <- EnhancedVolcano(memory.cd4,
                              lab = rownames(memory.cd4),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              xlim = c(-1, 1.5),
                              pCutoff = 0.05,
                              FCcutoff = 0.25,
                              legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                             'p-value & Log (base 2) FC'),
                              title = "CD4 Memory T-cell response",
                              subtitle = "Untreated vs ADP-heptose, 6h")
  
  cd4.memory.plot
  
  ggsave(filename = "cd4_memory_volcano.png", plot = cd4.memory.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}


#CD8 T-cells
{
  cd8.naive <- data.frame(`Naïve_CD8_T-cell_response`[,-1], row.names = `Naïve_CD8_T-cell_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  cd8.naive <- cd8.naive[,c(-1, -3, -4)]
  
  #plot the data
  cd8.plot <- EnhancedVolcano(cd8.naive,
                              lab = rownames(cd8.naive),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              pCutoff = 0.05,
                              FCcutoff = 0.25,
                              legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                             'p-value & Log (base 2) FC'),
                              title = "CD8 T-cell response",
                              subtitle = "Untreated vs ADP-heptose, 6h")
  cd8.plot
  
  ggsave(filename = "cd8_naive_volcano.png", plot = cd8.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}


#Memory CD8 T-cells
{
  cd8.memory <- data.frame(`Memory_CD8_T-cell_response`[,-1], row.names = `Memory_CD8_T-cell_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  cd8.memory <- cd8.memory[,c(-1, -3, -4)]
  
  #plot the data
  cd8.memory.plot <- EnhancedVolcano(cd8.memory,
                              lab = rownames(cd8.memory),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              pCutoff = 0.05,
                              FCcutoff = 0.25,
                              legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                             'p-value & Log (base 2) FC'),
                              title = "Memory CD8 T-cell response",
                              subtitle = "Untreated vs ADP-heptose, 6h")
  cd8.memory.plot
  
  ggsave(filename = "cd8_memory_volcano.png", plot = cd8.memory.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
}


#B-cells

{
  b.cell <- data.frame(`B-cells_response`[,-1], row.names = `B-cells_response`$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  b.cell <- b.cell[,c(-1, -3, -4)]
  
  #plot the data
  b.cell.plot <- EnhancedVolcano(b.cell,
                              lab = rownames(b.cell),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              pCutoff = 0.05,
                              FCcutoff = 0.25,
                              legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                             'p-value & Log (base 2) FC'),
                              title = "B-cell response",
                              subtitle = "Untreated vs ADP-heptose, 6h")
  
  b.cell.plot
  
  ggsave(filename = "b-cell_volcano.png", plot = b.cell.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
}

#NK-cells

{
  nk <- data.frame(NK_cell_response[,-1], row.names = NK_cell_response$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  nk <- nk[,c(-1, -3, -4)]
  
  #plot the data
  nk.plot <- EnhancedVolcano(nk,
                                 lab = rownames(nk),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                pCutoff = 0.05,
                                FCcutoff = 0.25,
                                 legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                                'p-value & Log (base 2) FC'),
                                 title = "NK cell response",
                                 subtitle = "Untreated vs ADP-heptose, 6h")
  
  nk.plot
}


#NKT cells

{
  nkt <- data.frame(NKT_cell_response[,-1], row.names = NKT_cell_response$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  nkt <- nkt[,c(-1, -3, -4)]
  
  #plot the data
  nkt.plot <- EnhancedVolcano(nkt,
                             lab = rownames(nkt),
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             pCutoff = 0.05,
                             FCcutoff = 0.25,
                             legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                            'p-value & Log (base 2) FC'),
                             title = "NKT cell response",
                             subtitle = "Untreated vs ADP-heptose, 6h")
  
  nkt.plot
}

#Dendritic cells
{
  dc <- data.frame(Conventional_DCs_response[, -1], row.names = Conventional_DCs_response$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  dc <- dc[,c(-1, -3, -4)]
  
  #plot the data
  dc.plot <- EnhancedVolcano(dc,
                              lab = rownames(dc),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              pCutoff = 0.05,
                              FCcutoff = 0.25,
                              legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                             'p-value & Log (base 2) FC'),
                              title = "Conventional DC response",
                              subtitle = "Untreated vs ADP-heptose, 6h")
  dc.plot
  
  ggsave(filename = "dc_volcano.png", plot = dc.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
}

#pDC
{
  pdc <- pdc <- data.frame(pDC_response[, -1], row.names = pDC_response$Gene)
  #subset the dataframe to include avg_log2fc and adj.pval
  pdc <- pdc[,c(-1, -3, -4)]
  
  #plot the data
  pdc.plot <- EnhancedVolcano(pdc,
                               lab = rownames(pdc),
                               x = 'avg_log2FC',
                               y = 'p_val_adj',
                               pCutoff = 0.05,
                               FCcutoff = 0.25,
                               legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                                              'p-value & Log (base 2) FC'),
                               title = "pDC response",
                               subtitle = "Untreated vs ADP-heptose, 6h")
  pdc.plot
  
  ggsave(filename = "pdc_volcano.png", plot = pdc.plot, path = "output/SCTransform.manual.annotation.output/", dpi = 320, 
         width = 7, height = 7, units = "in")
  
}



