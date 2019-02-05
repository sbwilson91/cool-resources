# GPL: Platform, describes list of elements in an array 
# GSM: Samples, describes conditions each sample was handled, manipulated etc.
# GSE: Series, defines sets of related samples in a group, how related, how ordered

# BiocManager::install("GEOquery")
#install.packages("Seurat")
library(GEOquery)
library(Seurat)
#library(Biobase)


# downloaded geo data
getGEOSuppFiles(GEO = "GSE102596")

untar("GSE102596/GSE102596_RAW.tar", exdir = "GSE102596/")

geo.data <- read.table("GSE102596/GSM2741551_count-table-human16w.tsv.gz", sep = "\t")

geo <- CreateSeuratObject(raw.data = geo.data, project = "GEO_DOWNLOAD")

geo <- NormalizeData(geo)

geo <- FindVariableGenes(object = geo, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
geo <- ScaleData(geo, display.progress = F)
geo <- RunPCA(object = geo, pc.genes = geo@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
PCAPlot(object = geo, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = geo)
geo <- FindClusters(object = geo, reduction.type = "pca", dims.use = 1:10, 
                    resolution = seq(0, 1, 0.1), print.output = 0, save.SNN = TRUE)
geo <- RunTSNE(geo, reduction.use = "pca", dims.use = 1:10, 
               do.fast = T)
TSNEPlot(object = geo)

FeaturePlot(geo, features.plot = "CDH6", cols.use = c("grey", "red"))

# P1 kidney
getGEOSuppFiles(GEO = "GSE94333")

untar("GSE94333/GSE94333_RAW.tar", exdir = "GSE94333/")

p1.data <- read.table("GSE94333/GSM2473320_run1871_out_gene_exon_tagged.dge.txt.gz", sep = "\t",header = T, row.names = 1)

p1 <- CreateSeuratObject(raw.data = p1.data, project = "P1_KIDNEY")
cols
p1 <- NormalizeData(p1)

p1 <- FindVariableGenes(object = p1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
p1 <- ScaleData(p1, display.progress = F)
p1 <- RunPCA(object = p1, pc.genes = p1@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
PCAPlot(object = p1, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = p1)
p1 <- FindClusters(object = p1, reduction.type = "pca", dims.use = 1:10, 
                    resolution = seq(0, 1, 0.1), print.output = 0, save.SNN = TRUE)
p1 <- RunTSNE(p1, reduction.use = "pca", dims.use = 1:10, 
               do.fast = T)
TSNEPlot(object = p1)

FeaturePlot(p1, features.plot = "CDH6", cols.use = c("grey", "red"))
