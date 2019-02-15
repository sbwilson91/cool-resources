linearModel <- function(object, n_features) {
  log_count <- as.matrix(logcounts(object))
  cols <- ncol(log_count)
  if (!"counts" %in% assayNames(object)) {
    warning("Your object does not contain counts() slot. Dropouts were calculated using logcounts() slot...")
    dropouts <- rowSums(log_count == 0)/cols * 100
  } else {
    count <- as.matrix(counts(object))
    dropouts <- rowSums(count == 0)/cols * 100
  }
  # do not consider spikes and genes with 0 and 100 dropout rate
  dropouts_filter <- dropouts != 0 & dropouts != 100
  if (!is.null(isSpike(object))) {
    for (spikes in spikeNames(object)) {
      dropouts_filter <- as.logical(dropouts_filter * (!isSpike(object, spikes)))
    }
  }
  dropouts_filter <- which(dropouts_filter)
  dropouts <- log2(dropouts[dropouts_filter])
  expression <- rowSums(log_count[dropouts_filter, ])/cols
  
  fit <- lm(dropouts ~ expression)
  gene_inds <- fit$residuals
  names(gene_inds) <- 1:length(gene_inds)
  gene_inds <- as.numeric(names(head(sort(gene_inds, decreasing = TRUE), n_features)))
  
  scmap_features <- rep(FALSE, nrow(object))
  scmap_features[dropouts_filter[gene_inds]] <- TRUE
  
  scmap_scores <- rep(NA, nrow(object))
  scmap_scores[dropouts_filter] <- fit$residuals
  
  d <- as.data.frame(cbind(expression, dropouts))
  d$Gene <- rownames(object)[dropouts_filter]
  d$Features <- "Other"
  d$Features[gene_inds] <- "Selected"
  d$Features <- factor(d$Features, levels = c("Selected", "Other"))
  
  return(list(scmap_features = scmap_features, scmap_scores = scmap_scores, for_plotting = d, fit = fit))
}


ggplot_features <- function(d, fit) {
  dropouts <- Features <- NULL
  cols <- c("#d73027", "#4575b4")
  p <- ggplot(d, aes(x = expression, y = dropouts, colour = Features)) + geom_point(size = 0.7) + 
    scale_colour_manual(values = cols) + labs(x = "log2(Expression)", y = "log2(% of dropouts)") + 
    geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1]) + theme_classic(base_size = 12)
  return(p)
}


indexCluster.SingleCellExperiment <- function(object, cluster_col) {
  if(!checks_for_index(object)) {
    return(object)
  }
  if (!cluster_col %in% colnames(colData(object))) {
    stop("Please define an existing cluster column of the `colData` slot of the input object using the `cluster_col` parameter!")
    return(object)
  }
  tmp <- object[rowData(object)$scmap_features, ]
  gene <- cell_class <- exprs <- NULL
  exprs_mat <- as.matrix(logcounts(tmp))
  rownames(exprs_mat) <- as.data.frame(rowData(tmp))$feature_symbol
  colnames(exprs_mat) <- as.data.frame(colData(tmp))[[cluster_col]]
  
  # calculate median feature expression in every cell class of object
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c("gene", "cell_class", "exprs")
  exprs_mat <- exprs_mat %>% group_by(gene, cell_class) %>% summarise(med_exprs = median(exprs))
  exprs_mat <- reshape2::dcast(exprs_mat, gene ~ cell_class, value.var = "med_exprs")
  rownames(exprs_mat) <- exprs_mat$gene
  index <- exprs_mat[, 2:ncol(exprs_mat), drop = FALSE]
  index <- index[order(rownames(index)), , drop = FALSE]
  index <- index[, colSums(index) > 0, drop = FALSE]
  if (ncol(index) == 0) {
    stop("scmap index is empty because the median expression in the selected features is 0 in every cell cluster! Try to increase the number of selected features!")
    return(object)
  }
  metadata(object)$scmap_cluster_index <- index
  return(object)
}


checks_for_index <- function(object) {
  if (is.null(object)) {
    stop("Please provide a `SingleCellExperiment` object using the `object` parameter!")
    return(FALSE)
  }
  if (!"SingleCellExperiment" %in% is(object)) {
    stop("Input object is not of `SingleCellExperiment` class! Please provide an object of the correct class!")
    return(FALSE)
  }
  if (is.null(rowData(object)$scmap_features)) {
    stop("Features are not selected! Please run `selectFeatures()` or `setFeatures()` first!")
    return(FALSE)
  }
  if (is.null(rowData(object)$feature_symbol)) {
    stop("There is no `feature_symbol` column in the `rowData` slot of the `reference` dataset! Please write your gene/transcript names to this column!")
    return(FALSE)
  }
  return(TRUE)
}

checks_for_projection <- function(projection, index_list) {
  if (is.null(projection)) {
    stop("Please provide a `SingleCellExperiment` object for the `projection` parameter!")
    return(FALSE)
  }
  if (!"SingleCellExperiment" %in% is(projection)) {
    stop("`projection` dataset has to be of the `SingleCellExperiment` class!")
    return(FALSE)
  }
  if (is.null(rowData(projection)$feature_symbol)) {
    stop("There is no `feature_symbol` column in the `rowData` slot of the `projection` dataset! Please write your gene/transcript names to this column!")
    return(FALSE)
  }
  if (is.null(index_list)) {
    stop("Please provide a list of precomputed scmap indexes as the `reference` parameter!")
    return(FALSE)
  }
  if (!"list" %in% is(index_list)) {
    stop("Please provide a list of precomputed scmap indexes as the `reference` parameter!")
    return(FALSE)
  }
  return(TRUE)
}


scmapCluster.SingleCellExperiment <- function(projection, index_list, threshold) {
  if(!checks_for_projection(projection, index_list)) {
    return(projection)
  }
  labels <- list()
  simls <- list()
  for (n in seq_len(length(index_list))) {
    index <- index_list[[n]]
    # find and select only common features, then subset both datasets
    tmp <- setFeatures(projection, rownames(index))
    index <- index[rownames(index) %in% rowData(tmp)$feature_symbol[rowData(tmp)$scmap_features], , drop = FALSE]
    tmp <- tmp[rowData(tmp)$scmap_features, ]
    
    if (nrow(index) < 10) {
      warning("There are less than ten features in common between the `reference` and `projection` datasets. Most probably they come from different organisms! Please redefine your query!")
      return(projection)
    }  
    
    # get expression values of the projection dataset
    proj_exprs <- as.matrix(logcounts(tmp))
    rownames(proj_exprs) <- rowData(tmp)$feature_symbol
    
    # prepare projection dataset
    proj_exprs <- proj_exprs[order(rownames(proj_exprs)), ]
    
    # calculate similarities and correlations
    tmp <- t(index)
    res <- proxy::simil(tmp, t(proj_exprs), method = "cosine")
    res <- matrix(res, ncol = nrow(tmp), byrow = TRUE)
    max_inds1 <- max.col(res)
    maxs1 <- rowMaxs(res)
    
    res <- cor(index, proj_exprs, method = "pearson")
    max_inds2 <- max.col(t(res))
    maxs2 <- colMaxs(res)
    
    res <- cor(index, proj_exprs, method = "spearman")
    max_inds3 <- max.col(t(res))
    maxs3 <- colMaxs(res)
    
    cons <- cbind(colnames(index)[max_inds1], colnames(index)[max_inds2], 
                  colnames(index)[max_inds3])
    
    maximums <- cbind(maxs1, maxs2, maxs3)
    
    # cells with at least one NA correlation value become unassigned
    non_na_inds <- !is.na(max_inds1) & !is.na(max_inds2) & !is.na(max_inds3)
    
    # create labels
    maxs <- rep(NA, nrow(cons))
    labs <- rep("unassigned", nrow(cons))
    unique_labs <- unlist(apply(cons, 1, function(x) {
      length(unique(x))
    }))
    
    ## all similarities agree
    if (length(which(unique_labs == 1 & non_na_inds)) > 0) {
      labs[unique_labs == 1 & non_na_inds] <- cons[unique_labs == 1 & non_na_inds, 1]
      maxs_tmp <- rowMaxs(maximums[unique_labs == 1 & non_na_inds, , drop = FALSE])
      maxs[unique_labs == 1 & non_na_inds] <- maxs_tmp
    }
    
    ## only two similarities agree
    if (length(which(unique_labs == 2 & non_na_inds)) > 0) {
      tmp <- cons[unique_labs == 2 & non_na_inds, , drop = FALSE]
      inds <- unlist(apply(tmp, 1, function(x) {
        which(duplicated(x))
      }))
      labs[unique_labs == 2 & non_na_inds] <- tmp[cbind(seq_along(inds), inds)]
      
      ## calculate maximum similarity in case of two agreeing similarities
      tmp1 <- matrix(apply(tmp, 2, `==`, labs[unique_labs == 2 & non_na_inds]), ncol = 3)
      inds <- t(apply(tmp1, 1, which))
      maxs_tmp <- cbind(maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 
                                                                                                      1]), inds[, 1])], maximums[unique_labs == 2 & non_na_inds, , drop = FALSE][cbind(seq_along(inds[, 
                                                                                                                                                                                                      1]), inds[, 2])])
      maxs_tmp <- rowMaxs(maxs_tmp)
      maxs[unique_labs == 2 & non_na_inds] <- maxs_tmp
    }
    
    ## check the similarity threshold
    labs[!is.na(maxs) & maxs < threshold] <- "unassigned"
    labels[[n]] <- labs
    simls[[n]] <- maxs
  }
  names(labels) <- names(index_list)
  names(simls) <- names(index_list)
  return(order_and_combine_labels(labels, simls))}
order_and_combine_labels <- function(labels, simls) {
  unassigned_rate_order <- order(
    unlist(
      lapply(labels, function(x) {
        length(x[x == "unassigned"])/length(x)
      })
    )
  )
  labels <- labels[unassigned_rate_order]
  simls <- simls[unassigned_rate_order]
  labels <- do.call(cbind, labels)
  simls <- do.call(cbind, simls)
  max_simls_inds <- apply(
    simls, 1, function(x) {
      if(!all(is.na(x))) {
        return(which.max(x))
      } else {
        return(NA)
      }
    }
  )
  inds <- which(!is.na(max_simls_inds))
  cons_labels <- rep("unassigned", nrow(labels))
  cons_labels[inds] <- labels[cbind(inds, max_simls_inds[inds])]
  return(list(scmap_cluster_labs = labels, scmap_cluster_siml = simls, combined_labs = cons_labels))
}
selectFeatures.SingleCellExperiment <- function(object, n_features, suppress_plot) {
  if (is.null(as.data.frame(rowData(object))$feature_symbol)) {
    stop("There is no feature_symbol column in the rowData slot! 
                Please create one and then run this function again. Please note
                that feature symbols in the reference dataset must correpond 
                to the feature symbols in the mapping dataset, otherwise the 
                mapping will not work!.")
    return(object)
  }
  
  r_data <- as.data.frame(rowData(object))
  tmp <- linearModel(object, n_features)
  r_data$scmap_features <- tmp$scmap_features
  r_data$scmap_scores <- tmp$scmap_scores
  rowData(object) <- r_data
  
  if (!suppress_plot) {
    p <- ggplot_features(tmp$for_plotting, tmp$fit)
    print(p)
  }
  
  return(object)
}
