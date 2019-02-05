# GPL: Platform, describes list of elements in an array 
# GSM: Samples, describes conditions each sample was handled, manipulated etc.
# GSE: Series, defines sets of related samples in a group, how related, how ordered

# BiocManager::install("GEOquery")
#install.packages("Seurat")
library(GEOquery)
library(Seurat)
#library(Biobase)


# downloaded geo data: GSE118184 single cell data from humphries comparison paper
gse <- getGEO(GEO = "GSE118184", destdir = "./GSE118184/")

getGEOSuppFiles(GEO = "GSE118184")
gse.files <- list.files(path = "./GSE118184/")

#untar(tarfile = "./GSE118184/GSE118184_Human_kidney_snRNA.dge.txt.gz", exdir = "./GSE118184/extracted/")
# single nuclear rna human kidney dataset
snrna.data <- read.table("GSE118184/GSE118184_Human_kidney_snRNA.dge.txt.gz", sep = "\t")

snrna <- CreateSeuratObject(raw.data = snrna.data, project = "Human_Kidney_snRNA")

snrna <- NormalizeData(snrna)

snrna <- FindVariableGenes(object = snrna, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
snrna <- ScaleData(snrna, display.progress = F)
snrna <- RunPCA(object = snrna, pc.genes = snrna@var.genes, do.print = F)
#PCAPlot(object = geo, dim.1 = 1, dim.2 = 2)
#PCElbowPlot(object = geo)
snrna <- FindClusters(object = snrna, reduction.type = "pca", dims.use = 1:10, 
                    resolution = seq(0, 1, 0.2), print.output = 0, save.SNN = TRUE)
snrna <- RunTSNE(snrna, reduction.use = "pca", dims.use = 1:10, 
               do.fast = T)
clustree::clustree(snrna)
TSNEPlot(object = snrna, group.by = "res.0.6")
snrna <- SetIdent(snrna, ident.use = snrna@meta.data$res.0.6)

FeaturePlot(snrna, features.plot = c("TMEM213", "SLC12A1", "SLC12A3", "LRP2"), cols.use = c("grey", "red"))

markers <- FindAllMarkers(snrna, min.pct = 0.1, logfc.threshold = 1)

markers.list <- lapply(0:(length(unique(snrna@ident))-1), function(x) {
  markers %>%
    dplyr::filter(cluster == x, p_val_adj < 0.05, avg_logFC > 0) %>%
    dplyr::arrange(-avg_logFC) %>%
    select(Gene = gene, LogFC = avg_logFC, pVal = p_val_adj)
})
write_csv(markers, "./GSE118184/markers_snrna.csv")
top <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

write_rds(x = snrna, path = "./GSE118184/seurat/scrna_seurat.rds")

# takasato ipsc day 34
mt.ipsc.d34.data <- read.table("GSE118184/GSE118184_Takasato_iPS_day34.dge.txt.gz", sep = "\t")

mt.ipsc.d34 <- CreateSeuratObject(raw.data = mt.ipsc.d34.data, project = "Takasato_ipsc_d34")

mt.ipsc.d34 <- NormalizeData(mt.ipsc.d34)

mt.ipsc.d34 <- FindVariableGenes(object = mt.ipsc.d34, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
mt.ipsc.d34 <- ScaleData(mt.ipsc.d34, display.progress = F)
mt.ipsc.d34 <- RunPCA(object = mt.ipsc.d34, pc.genes = mt.ipsc.d34@var.genes, do.print = F)
#PCAPlot(object = geo, dim.1 = 1, dim.2 = 2)
#PCElbowPlot(object = geo)
mt.ipsc.d34 <- FindClusters(object = mt.ipsc.d34, reduction.type = "pca", dims.use = 1:10, 
                      resolution = seq(0, 1, 0.2), print.output = 0, save.SNN = TRUE)
mt.ipsc.d34 <- RunTSNE(mt.ipsc.d34, reduction.use = "pca", dims.use = 1:10, 
                 do.fast = T)
#clustree::clustree(mt.ipsc.d34)
TSNEPlot(object = mt.ipsc.d34, group.by = "res.0.6", do.label = T)
mt.ipsc.d34 <- SetIdent(mt.ipsc.d34, ident.use = mt.ipsc.d34@meta.data$res.0.6)

FeaturePlot(mt.ipsc.d34, features.plot = c("COL1A1", "SLC12A1", "LYPD1", "EPCAM"), cols.use = c("grey", "red"))
FeaturePlot(mt.ipsc.d34, features.plot = c("PAX2", "PAX8", "COL3A1", "CRABP2"), cols.use = c("grey", "red"))

markers <- FindAllMarkers(mt.ipsc.d34, min.pct = 0.1, logfc.threshold = 1)

markers.list <- lapply(0:(length(unique(mt.ipsc.d34@ident))-1), function(x) {
  markers %>%
    dplyr::filter(cluster == x, p_val_adj < 0.05, avg_logFC > 0) %>%
    dplyr::arrange(-avg_logFC) %>%
    select(Gene = gene, LogFC = avg_logFC, pVal = p_val_adj)
})
write_csv(markers, "./GSE118184/markers_mt_ipsc_d34.csv")
top <- markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

write_rds(x = mt.ipsc.d34, path = "./GSE118184/seurat/mt_ipsc_d34_seurat.rds")
