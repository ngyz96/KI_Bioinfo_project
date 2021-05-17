library(Seurat)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_01_data.rds')

#normalisation for each individual sample
data.list <- SplitObject(data, split.by = 'cell_ident')
for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
}

data.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:50)
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:50)

DefaultAssay(data.integrated) <- "integrated"

data.integrated <- ScaleData(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, npcs = 50, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, dims = 1:30)

#SNN graph and clustering
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 0.4)

saveRDS(data.integrated, file = 'R/Data/combined_02b_data.rds')