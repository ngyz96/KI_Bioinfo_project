library(Seurat)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_01_data.rds')

#normalisation for each individual sample
data.list <- SplitObject(data, split.by = 'cell_ident')
rm(data)
for (i in 1:length(data.list)) {
    data.list[[i]] <- SCTransform(data.list[[i]], verbose = FALSE)
}

#prepare and integrate sample data
options(future.globals.maxSize = 3800 * 1024^2)
data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, 
                                    verbose = T)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                           anchor.features = data.features, verbose = FALSE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

#run dimensionality reduction
data.integrated <- RunPCA(data.integrated, npcs = 50, verbose = FALSE)
data.integrated <- RunUMAP(data.integrated, dims = 1:30)

#SNN graph and clustering
data.integrated <- FindNeighbors(data.integrated,, dims = 1:30)
data.integrated <- FindClusters(data.integrated, resolution = 0.4)

saveRDS(data.integrated, file = 'R/Data/combined_02c_data.rds')