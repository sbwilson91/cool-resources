# Script used for the umap dimensional reduction for single cell datasets.

library(reticulate)
library(Seurat)
data <- readRDS("S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/Spence_HFK_D96-D108_combined.rds")
data <- RunUMAP(data, dims = 1:30, seed.use = 25)
saveRDS(data, "S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/Spence_HFK_D96-D108_combined.rds")
DefaultAssay(data) <- "RNA"

library(reticulate)
library(Seurat)
data <- readRDS("S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/RETD7.rds")
data <- RunUMAP(data, dims = 1:30, seed.use = 8)
saveRDS(data, file = "S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/RETD7.rds")


library(reticulate)
library(Seurat)
data <- readRDS("S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/UretericEp_Cult_Combined.rds")
data <- RunUMAP(data, dims = 1:30, seed.use = 8, dim.embed = 3)
saveRDS(data, file = "S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/UretericEp_Cult_Combined.rds")


library(reticulate)
library(Seurat)
data <- readRDS("S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/Lindstrom_wk16_HFK_cortex.rds")
data <- RunUMAP(data, dims = 1:30, seed.use = 15, dim.embed = 3)
saveRDS(data, file = "S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/Lindstrom_wk16_HFK_cortex.rds")

library(reticulate)
library(Seurat)
data <- readRDS("S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/UE_merge_anchors.rds")
data <- RunUMAP(data, dims = 1:30, seed.use = 2, dim.embed = 3)
saveRDS(data, file = "S:/PhD/R-projects/UretericEpithelialCulture/output/seurat/UE_merge_anchors.rds")


