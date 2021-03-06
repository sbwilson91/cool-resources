# install scmap
#BiocManager::install("scmap", version = "3.8")
library(scmap)
library(SingleCellExperiment)
library(GEOquery)
library(Seurat)
library(clustree)
library(tidyverse)


# download human fetal kidney GEO dataset
getGEOSuppFiles(GEO = "GSE102596")

untar("GSE102596/GSE102596_RAW.tar", exdir = "GSE102596/")

hfk.data <- read.table("GSE102596/GSM2741551_count-table-human16w.tsv.gz", sep = "\t")

hfk.data[1:4, 1:4]

hfk.sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(hfk.data)))
rowD

hfk.sce <- selectFeatures(hfk.sce, suppress_plot = F)

hfk <- CreateSeuratObject(raw.data = hfk.data, project = "Human_Kidney_hfk")

hfk <- NormalizeData(hfk)

hfk <- FindVariableGenes(object = hfk, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
hfk <- ScaleData(hfk, display.progress = F)
hfk <- RunPCA(object = hfk, pc.genes = hfk@var.genes, do.print = F)
#PCAPlot(object = geo, dim.1 = 1, dim.2 = 2)
PCElbowPlot(object = hfk)
hfk <- FindClusters(object = hfk, reduction.type = "pca", dims.use = 1:10, 
                      resolution = seq(0, 2, 0.2), print.output = 0, save.SNN = TRUE)
hfk <- RunTSNE(hfk, reduction.use = "pca", dims.use = 1:10, 
                 do.fast = T)
clustree(hfk)
TSNEPlot(object = hfk, group.by = "res.0.6")
hfk <- SetIdent(hfk, ident.use = hfk@meta.data$res.0.6)

FeaturePlot(hfk, features.plot = c("TMEM213", "GATA3", "CD14", "TCF21"), cols.use = c("grey", "red"))

markers <- FindAllMarkers(hfk, min.pct = 0.1, logfc.threshold = 0.5)

markers.list <- lapply(0:(length(unique(hfk@ident))-1), function(x) {
  markers %>%
    dplyr::filter(cluster == x, p_val_adj < 0.05, avg_logFC > 0) %>%
    dplyr::arrange(-avg_logFC) %>%
    select(Gene = gene, LogFC = avg_logFC, pVal = p_val_adj)
})

head(markers.list[[10]],20)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10)
new.cluster.ids <- c("Stroma 1", "Stroma 2", "Neph Progenitors", "Stroma DCN", "Endothelium", "Stroma 3", 
                     "Macrophage", "Coll Duct/Dist Tubule", "Prox Tubule/Gloms", "Cell Cycle", "Blood")
hfk@ident <- plyr::mapvalues(hfk@ident,from = current.cluster.ids, to=new.cluster.ids)


TSNEPlot(object = hfk, do.label = T)

hfk.sce <- Convert(from = hfk, to = "sce")

rowData(hfk.sce)$feature_symbol <- rownames(hfk.sce)
hfk.sce <- hfk.sce[!duplicated(rownames(hfk.sce)),]
hfk.sce <- selectFeatures(hfk.sce, suppress_plot = F)
r_data <- as.data.frame(rowData(hfk.sce))

log_count <- as.matrix(logcounts(hfk.sce))
cols <- ncol(log_count)

count <- as.matrix(counts(hfk.sce))
dropouts <- rowSums(count == 0)/cols * 100
tmp <- linearModel(hfk.sce, 500)
r_data$scmap_features <- tmp$scmap_features
r_data$scmap_scores <- tmp$scmap_scores
rowData(hfk.sce) <- r_data
p <- ggplot_features(tmp$for_plotting, tmp$fit)
plot(p)
table(rowData(hfk.sce)$scmap_features)
hfk.sce <- indexCluster(hfk.sce, cluster_col = "ident")

tmp <- hfk.sce[rowData(hfk.sce)$scmap_features,]
gene <- cell_class <- exprs <- NULL
exprs_mat <- as.matrix(logcounts(tmp))
rownames(exprs_mat) <- as.data.frame(rowData(tmp))$feature_symbol
colnames(exprs_mat) <- as.data.frame(colData(tmp))[["ident"]]

hfk.sce <- indexCluster.SingleCellExperiment(hfk.sce, cluster_col = "ident")

head(metadata(hfk.sce)$scmap_cluster_index)

heatmap(as.matrix(metadata(hfk.sce)$scmap_cluster_index))


# Projection of dataset back onto itself

scmapCluster_results <- scmapCluster.SingleCellExperiment(projection = hfk.sce,
                                     index_list = list(hfk = metadata(hfk.sce)$scmap_cluster_index),
                                     threshold = 0.2)
head(scmapCluster_results$scmap_cluster_labs)

plot(
  getSankey(
    colData(hfk.sce)$ident, 
    scmapCluster_results$scmap_cluster_labs[,'hfk'],
    plot_height = 800,
    plot_width = 1200
  )
)
head(scmapCluster_results$scmap_cluster_labs)
hfk@meta.data$scmap <- scmapCluster_results$scmap_cluster_labs[,'hfk']
# Plot scmap assigned cluster identity to tsne plot of data
plot1 <- TSNEPlot(hfk, do.return=T, no.legend = T, do.label = T)
plot2 <- TSNEPlot(hfk, do.label =F, group.by = "scmap")
plot_grid(plot1, plot2)


# Compare our organoid datasets to the hfk dataset

reporters <- readRDS("../../../../../GROUPS/Stem cells and Regeneration/Single Cell/Reporters_18vs25Days/output/seurat/D18vsD25_combined_20clusters_seurat.rds")
TSNEPlot(reporters, do.label = T, group.by = "res.1.6")
reporters <- SetIdent(reporters, ident.use = reporters@meta.data$res.1.6)
TSNEPlot(reporters)

rep.sce <- Convert(from = reporters, to = "sce")
rowData(rep.sce)$feature_symbol <- rownames(rep.sce)
rep.sce <- selectFeatures.SingleCellExperiment(rep.sce, n_features = 500, suppress_plot = F)
rep.sce <- indexCluster.SingleCellExperiment(rep.sce, cluster_col = "res.1.6")
heatmap(as.matrix(metadata(rep.sce)$scmap_cluster_index))

scmapCluster_results_rep <- scmapCluster.SingleCellExperiment(
  projection = rep.sce, 
  index_list = list(
    reporters = metadata(hfk.sce)$scmap_cluster_index
  ), threshold = 0.7
)

head(scmapCluster_results_rep$scmap_cluster_labs)
head(scmapCluster_results_rep$scmap_cluster_siml)


plot(
  getSankey(
    colData(rep.sce)$"res.1.6", 
    scmapCluster_results_rep$scmap_cluster_labs[,'reporters'],
    plot_height = 800,
    plot_width = 800
  )
)

qplot(x = scmapCluster_results_rep$scmap_cluster_siml[,'reporters'])
reporters@meta.data$scmap <- scmapCluster_results_rep$scmap_cluster_labs[,'reporters']
# Plot scmap assigned cluster identity to tsne plot of data
plot1 <- TSNEPlot(reporters, do.return=T, no.legend = F, do.label = T)
plot2 <- TSNEPlot(reporters, do.label =T, group.by = "scmap")
plot_grid(plot1, plot2)

reporters <- SetIdent(reporters, ident.use = reporters@meta.data$scmap)
neph.prog <- SubsetData(reporters, ident.use = "Neph Progenitors")
tubules <- SubsetData(reporters, ident.use = c("Prox Tubule/Gloms", "Coll Duct/Dist Tubule"))
unassigned <- SubsetData(reporters, ident.use = c("unassigned"))
plot1 <- TSNEPlot(neph.prog, do.return=T, no.legend = T, do.label = T)
plot2 <- TSNEPlot(tubules, do.return = T, no.legend=T, do.label = T)
TSNEPlot(unassigned, do.return = T, no.legend = F)
plot_grid(plot1, plot2)


