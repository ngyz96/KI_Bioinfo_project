library(Seurat)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_01_data.rds')
#processing
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, npcs = 50, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30)

#SNN graph and clustering
data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
data <- FindClusters(data, resolution = 0.4)

saveRDS(data, file = 'R/Data/combined_02f_data.rds')