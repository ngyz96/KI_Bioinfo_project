library(Seurat)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_01_data.rds')

#normalisation
data <- SCTransform(object = data)

#run dimensionality reduction
data <- RunPCA(data, npcs = 50, verbose = FALSE)
data <- RunUMAP(data, dims = 1:30)

#SNN graph and clustering
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = 0.4)

saveRDS(data, file = 'R/Data/combined_02a_data.rds')