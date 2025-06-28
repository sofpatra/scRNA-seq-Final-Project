
#load needed libraries
library(Seurat)
library(tidyverse)
library(DoubletFinder)


samples <- c("HB17_background", "HB17_PDX", "HB17_tumor", "HB30_PDX", "HB30_tumor", "HB53_background",
             "HB53_tumor")

seurat_objects <- readRDS("seurat_objects.rds")

# Create list to hold filtered objects (singlets)
filtered_seurat_objects <- list()

# Loop through each sample
for (sample in samples) {
  cat("Processing sample:", sample, "\n")
  
  obj <- seurat_objects[[sample]]
  
  # Basic preprocessing
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  
  # Clustering + UMAP
  obj <- FindNeighbors(obj, dims = 1:10, reduction = "pca")
  obj <- FindClusters(obj, resolution = 1, cluster.name = "unintegrated_clusters")
  obj <- RunUMAP(obj, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")
  
  # DoubletFinder parameter sweep
  sweep.res.list <- paramSweep(obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nCells <- ncol(obj)
  doublet_rate <- 0.01
  nExp <- round(nCells * doublet_rate)
  
  # Run DoubletFinder
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp, sct = FALSE)
  
  # Automatically find the classification column (it’s usually the last one)
  df_col <- tail(colnames(obj@meta.data), 1)
  
  # Plot and table
  print(DimPlot(obj, reduction = "umap.unintegrated", group.by = df_col) + ggtitle(paste(sample, "- Doublets")))
  print(table(obj@meta.data[[df_col]]))
  
  # Subset singlets
  singlet_cells <- rownames(obj@meta.data)[obj@meta.data[[df_col]] == "Singlet"]
  obj <- obj[, singlet_cells]
  
  # Store filtered object
  filtered_seurat_objects[[sample]] <- obj
  
  # Clean up
  rm(obj)
  gc()
}


#save filtered objects 
saveRDS(filtered_seurat_objects, file = "doublet_filtered_seurat_objects.rds")
