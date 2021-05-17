library(Seurat)
library(dplyr)

working_directory <- '/proj/uppstore2019102/ESCG_data/10X_19_ML_13_R38'
setwd(working_directory)

#load in data from old rats, create seurat object
old1.data <- Read10X(data.dir = '10X_19_057/outs/filtered_feature_bc_matrix')
old2.data <- Read10X(data.dir = '10X_19_059/outs/filtered_feature_bc_matrix')
old3.data <- Read10X(data.dir = '10X_19_061/outs/filtered_feature_bc_matrix')

old1 <- CreateSeuratObject(counts=old1.data, project="SO1", min.cells = 5, min.features = 200)
old2 <- CreateSeuratObject(counts=old2.data, project="SO2", min.cells = 5, min.features = 200)
old3 <- CreateSeuratObject(counts=old3.data, project="SO3", min.cells = 5, min.features = 200)

#merge and remove used objects from memory
old <- merge(old1, y=c(old2,old3), project = 'old')
rm(old1,old1.data,old2,old2.data,old3,old3.data)

#repeat the same with the young rats data
young1.data <- Read10X(data.dir = '10X_19_058/outs/filtered_feature_bc_matrix')
young2.data <- Read10X(data.dir = '10X_19_060/outs/filtered_feature_bc_matrix')
young3.data <- Read10X(data.dir = '10X_19_062/outs/filtered_feature_bc_matrix')

young1 <- CreateSeuratObject(counts=young1.data, project="SY1", min.cells = 5, min.features = 200)
young2 <- CreateSeuratObject(counts=young2.data, project="SY2", min.cells = 5, min.features = 200)
young3 <- CreateSeuratObject(counts=young3.data, project="SY3", min.cells = 5, min.features = 200)

young <- merge(young1, y=c(young2,young3), project = 'young')
rm(young1,young1.data,young2,young2.data,young3,young3.data)

#merge all data
data <- merge(old, young, project='youngvsold')
rm(old, young)
#save cell identity for downstream analysis
data <- StashIdent(object = data, save.name = 'cell_ident')
#save age group of rats
data <- RenameIdents(data, 'SO1' = 'old', 'SO2' = 'old', 'SO3' = 'old', 'SY1' = 'young', 'SY2' = 'young', 'SY3' = 'young')
data[['age']] <- Idents(data)
#calculate mitochondria gene
data <- PercentageFeatureSet(data, pattern = '^Mt-', col.name = 'percent.mt')
#QC filtering
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 40)

#save data file
saveRDS(data, file = 'R/Data/combined_01_data.rds')