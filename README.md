# cool-resources
Will use this to keep track of packages/tools/resources that are useful throughout my PhD and beyond
## AccessingData-Bioconductor.R
Contains R script allowing for the download of GEO data
## GEO-Seurat.R
Also contains R script to download GEO data then convert to seurat
## scmap
A tool for unsupervised projection of single cell RNA-seq data:
https://github.com/hemberg-lab/scmap

## UMAP.R
This script was made due to the lack of compatibility with the MCRI Rstudio server and using the umap-learn python package.
It contains the lines I used to import an .rds seurat file, perform the UMAP dimensional reduction, and resave the .rds file with the new UMAP computations. Seurat allows for mapping functions and including which reduction to use, so as along as the file contains the coordinates the umap package does not need to be accessible.

## reproducible_env.R
Rmarkdown file containing instructions on using renv for reproducible environments. I couldn't get it to properly work, so will stick with packrat for now.

# R graph Gallery

http://r-graph-gallery.com/

Location for formatted graphing scripts for all various graphs