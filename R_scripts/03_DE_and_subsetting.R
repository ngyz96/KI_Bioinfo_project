library(Seurat)
library(dplyr)

#set directory
working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load previous raw data
data <- readRDS('R/Data/combined_02d_mnn_data.rds')

#find cluster markers for all clusters
markers_all <- FindAllMarkers(data, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers_all, 'R/markers_all.csv')
#subsetting astrocyte like population
astrocytes <- subset(data, idents = c(0, 5, 6, 13))
saveRDS(astrocytes, file = 'R/Data/astrocytes_mnn_data.rds')

#performing DE on cluster 13 from 0,5,6
markers_13 <- FindMarkers(astrocytes, ident.1 = 13, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers_13, 'R/markers_13.csv')