library(Seurat)
library(harmony)
library(SeuratWrappers)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_01_data.rds')

data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data, npcs = 50)
data <- RunHarmony(data, group.by.vars = "cell_ident")
data <- RunUMAP(data, reduction = "harmony", dims = 1:30)
data <- FindNeighbors(data, reduction = "harmony", dims = 1:30)
data <- FindClusters(data, resolution = 0.4)

saveRDS(data, file = 'R/Data/combined_02e_harmony_data.rds')