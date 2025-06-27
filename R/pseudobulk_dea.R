pseudobulk_profiles <- function(seurat_obj, sample_names, cluster_col = "integrated_snn_res.0.5") {
  # Determine the union of genes present in any sample's counts layer
  LayerNames <- Layers(so.split)
  LayerNames <- LayerNames[grepl("counts", LayerNames)]
  
  all_genes <- character()
  for (i in seq_along(LayerNames)) {
    layer_name <- LayerNames[i]
    counts_layer <- LayerData(seurat_obj, assay = "RNA", layer = layer_name)
    if (!is.null(counts_layer)) {
      all_genes <- union(all_genes, rownames(counts_layer))
    }
  }
  
  # Get unique clusters from metadata
  clusters <- unique(seurat_obj@meta.data[[cluster_col]])
  
  # Initialize a list to store pseudobulk matrices per cluster
  pseudobulk_list <- list()
  
  # Loop over each cluster
  for (clust in clusters) {
    
    # Identify cells belonging to the current cluster
    cells_in_cluster <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[cluster_col]] == clust]
    
    # Initialize the pseudobulk matrix: union of genes x samples
    pseudobulk_mat <- matrix(0, nrow = length(all_genes), ncol = length(LayerNames))
    rownames(pseudobulk_mat) <- all_genes
    colnames(pseudobulk_mat) <- gsub("counts\\.", "", LayerNames) # Remove prefix from sample names
    
    # Loop over each sample/layer
    for (i in seq_along(LayerNames)) {
      layer_name <- LayerNames[i]
      counts_layer <- LayerData(seurat_obj, assay = "RNA", layer = layer_name)
      
      # Ensure that the layer exists
      if (is.null(counts_layer)) next
      
      # Identify cells in this sample that belong to the current cluster
      cells_in_sample_cluster <- intersect(colnames(counts_layer), cells_in_cluster)
      
      if (length(cells_in_sample_cluster) > 0) {
        # Sum the counts per gene for the cells in this sample and cluster
        sample_counts <- rowSums(counts_layer[, cells_in_sample_cluster, drop = FALSE])
        # Create a full vector for all genes; missing genes remain zero
        full_counts <- setNames(rep(0, length(all_genes)), all_genes)
        full_counts[names(sample_counts)] <- sample_counts
        pseudobulk_mat[, i] <- full_counts
      }
    }
    
    pseudobulk_list[[as.character(clust)]] <- pseudobulk_mat
  }
  
  return(pseudobulk_list)
}



run_deseq2 <- function(pbp,
                       cluster_col,
                       colData, 
                       design = ~ condition,
                       pseudocount = 0) {
  
  # Ensure provided colData is ordered by sample names
  colData <- colData[order(rownames(colData)), ]
  
  # Initialize list to store DESeq2 datasets for each cluster
  dds_list <- list()
  cell_types <- names(pbp)
  
  for (ct in cell_types) {
    # Retrieve the pseudobulk counts matrix as a data frame
    curr_counts <- as.data.frame(pbp[[ct]])
    # Reorder columns to match the rownames of colData
    curr_counts <- curr_counts[, rownames(colData), drop = FALSE]
    
    # Remove genes with all zeros across all samples
    nonzero_genes <- rowSums(curr_counts) > 0
    curr_counts <- curr_counts[nonzero_genes, , drop = FALSE]
    
    # Identify samples (columns) with nonzero counts; drop those that are entirely zero
    nonzero_samples <- colSums(curr_counts) > 0
    if (any(!nonzero_samples)) {
      # Optionally, you could issue a warning or message:
      message("Cluster ", ct, ": Dropping samples with all-zero counts: ",
              paste(colnames(curr_counts)[!nonzero_samples], collapse = ", "))
      curr_counts <- curr_counts[, nonzero_samples, drop = FALSE]
      local_colData <- colData[nonzero_samples, , drop = FALSE]
    } else {
      local_colData <- colData
    }
    
    # Optionally add a pseudocount (if pseudocount != 0, this ensures no zeros)
    if (pseudocount != 0) {
      curr_counts <- curr_counts + pseudocount
    }
    
    # Create DESeq2 dataset using the specified design formula
    dds <- DESeqDataSetFromMatrix(countData = curr_counts,
                                  colData = local_colData,
                                  design = design)
    dds <- DESeq(dds)
    dds_list[[ct]] <- dds
  }
  
  return(dds_list)
}

extract_contrasts <- function(dds_list, ct) {
  # Retrieve the DESeq2 object for the given cell type
  dds <- dds_list[[ct]]
  
  # Contrast: CS vs WT
  contrast_CS_WT <- contraster(dds,
                               group1 = list(c("condition", "CS")),
                               group2 = list(c("condition", "WT")))
  res_CS_WT <- results(dds, contrast = contrast_CS_WT)
  
  # Contrast: CS vs KO
  contrast_CS_KO <- contraster(dds,
                               group1 = list(c("condition", "CS")),
                               group2 = list(c("condition", "KO")))
  res_CS_KO <- results(dds, contrast = contrast_CS_KO)
  
  # Contrast: KO vs WT
  contrast_KO_WT <- contraster(dds,
                               group1 = list(c("condition", "KO")),
                               group2 = list(c("condition", "WT")))
  res_KO_WT <- results(dds, contrast = contrast_KO_WT)
  
  # Return a named list of results
  return(list("CS_vs_WT" = res_CS_WT,
              "CS_vs_KO" = res_CS_KO,
              "KO_vs_WT" = res_KO_WT))
  
  
  
}

add_fdrtool <- function(contrast_results_list) {
  result_list <- list()
  
  for (contrast_name in names(contrast_results_list)) {
    # Convert DESeq2 results object to a data frame
    res <- contrast_results_list[[contrast_name]]
    df <- as.data.frame(res)
    
    # Extract the test statistic column instead of the p-value column
    stats <- df$stat
    fdr_p <- rep(NA, length(stats))
    fdr_q <- rep(NA, length(stats))
    
    # Identify indices with non-NA statistic values
    valid_idx <- which(!is.na(stats))
    
    if (length(valid_idx) > 0) {
      # Use fdrtool with statistic = "normal" on the 'stat' column
      ft <- fdrtool(stats[valid_idx], statistic = "normal", plot = FALSE)
      
      fdr_p[valid_idx] <- ft$pval
      fdr_q[valid_idx] <- ft$qval
    }
    
    # Add the new fdrtool columns to the data frame
    df$fdrtool_pval <- fdr_p
    df$fdrtool_qval <- fdr_q
    
    result_list[[contrast_name]] <- df
  }
  
  return(result_list)
}

#’ Run DESeq2 on pseudobulk profiles, either per-cluster or in one combined model
#’
#’ @param pbp Named list of pseudobulk matrices (genes × samples), one entry per cluster
#’ @param colData Data.frame of sample metadata; rownames must match colnames of each matrix
#’ @param condition_col String; name of column in colData giving the condition factor (e.g. "condition" or other)
#’ @param cluster_mode String, one of "perCluster" or "combined"
#’ @param design Optional formula for DESeq2 design; if NULL, defaults to:
#’   - perCluster: ~ condition_col
#’   - combined: ~ cluster + condition + cluster:condition
#’ @param reference_level String; reference level for the primary condition factor (default "WT")
#’ @param min_reps Integer; minimum samples per condition required (default 2)
#’ @param filter_min_counts Integer; drop genes with fewer than this total count (default 10)
#’ @param lfc_shrink Logical; whether to apply DESeq2’s lfcShrink() to results (default TRUE)
#’ @param shrink_type String, one of "apeglm", "ashr", or "normal" (default "apeglm"); type of shrinkage to use
#’ @param BPPARAM BiocParallelParam (optional) for parallel per-cluster runs
#’ @return
#’ - If `cluster_mode = "perCluster"`, a named list of tibbles, one per cluster × contrast  
#’ - If `cluster_mode = "combined"`, a named list of tibbles of within-cluster contrasts  
#’ @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink resultsNames
#’ @importFrom tibble rownames_to_column as_tibble
#’ @importFrom dplyr filter pull
#’ @importFrom BiocParallel bplapply
run_pseudobulk_deseq2 <- function(pbp,
                                  colData,
                                  condition_col,
                                  cluster_mode      = c("perCluster", "combined"),
                                  design            = NULL,
                                  reference_level   = "WT",
                                  min_reps          = 2,
                                  filter_min_counts = 10,
                                  lfc_shrink        = TRUE,
                                  shrink_type       = c("apeglm","ashr","normal"),
                                  BPPARAM           = NULL) {
  library(DESeq2); library(dplyr); library(tibble)
  cluster_mode <- match.arg(cluster_mode)
  shrink_type  <- match.arg(shrink_type)
  
  # 1. Validate inputs
  stopifnot(is.list(pbp), is.data.frame(colData))
  stopifnot(condition_col %in% colnames(colData))
  
  # 2. Set default designs
  if (is.null(design)) {
    per_design  <- as.formula(paste("~", condition_col))
    comb_design <- ~ cluster + condition + cluster:condition
  } else {
    per_design  <- design
    comb_design <- design
  }
  
  # 3. Prepare combined data if needed
  if (cluster_mode == "combined") {
    # 3.1 Ensure all pseudobulk matrices share the same genes
    common_genes <- Reduce(intersect, lapply(pbp, rownames))
    pbp <- lapply(pbp, function(m) m[common_genes, , drop=FALSE])
    
    # 3.2 Prefix columns with cluster name to keep unique sample IDs
    pbp_named <- mapply(function(mat, cl) {
      cn <- colnames(mat)
      new_cn <- paste0(cl, "__", cn)
      colnames(mat) <- new_cn
      mat
    }, pbp, names(pbp), SIMPLIFY = FALSE)
    
    # 3.3 Combine count matrices
    combined_counts <- do.call(cbind, pbp_named)
    
    # 3.4 Build combined metadata by splitting the new colnames
    cn <- colnames(combined_counts)
    parts <- do.call(rbind, strsplit(cn, "__", fixed = TRUE))
    clusters <- parts[,1]
    samples  <- parts[,2]
    
    # 3.5 Subset and repeat colData for each cluster-sample
    if (!all(samples %in% rownames(colData))) {
      stop("Sample names from pseudobulk do not match colData rownames")
    }
    combined_meta <- colData[samples, , drop=FALSE]
    rownames(combined_meta) <- cn
    
    # 3.6 Add 'cluster' and relevel 'condition'
    combined_meta$cluster   <- factor(clusters, levels = names(pbp))
    combined_meta$condition <- relevel(
      factor(combined_meta[[condition_col]]),
      ref = reference_level
    )
    
    # 3.7 Check replicates per cluster-condition and filter genes
    if (any(table(combined_meta$cluster, combined_meta$condition) < min_reps)) {
      warning("Some cluster–condition combos have fewer than ", min_reps, " samples")
    }
    keep_genes       <- rowSums(combined_counts) >= filter_min_counts
    combined_counts <- combined_counts[keep_genes, , drop=FALSE]
  }
  
  # 4. Helper: get results with shrinkage and robust conversion to data.frame
  get_tidy_res <- function(dds, coef) {
    # Single coefficient shrink if applicable
    if (length(coef) == 1 && coef %in% resultsNames(dds)) {
      res_obj <- if (lfc_shrink) lfcShrink(dds, coef = coef, type = shrink_type) else results(dds, name = coef)
    } else {
      res_obj <- results(dds, contrast = coef)
      if (lfc_shrink && shrink_type == "ashr") {
        # apply ashr shrink for contrasts
        library(ashr)
        ash_out <- ash(res_obj$log2FoldChange, res_obj$lfcSE)
        res_obj$log2FoldChange <- ash_out$PosteriorMean
        res_obj$lfcSE          <- ash_out$PosteriorSD
      }
    }
    res_df <- as.data.frame(res_obj, row.names = rownames(res_obj))
    as_tibble(rownames_to_column(res_df, var = "gene"))
  }
  
  # 5a. Per-cluster analysis
  if (cluster_mode == "perCluster") {
    run_one <- function(clust) {
      mat <- pbp[[clust]][, rownames(colData), drop=FALSE]
      mat <- mat[rowSums(mat) >= filter_min_counts, , drop=FALSE]
      mat <- mat[, colSums(mat) > 0, drop=FALSE]
      meta <- colData[colnames(mat), , drop=FALSE]
      meta$condition <- relevel(factor(meta[[condition_col]]), ref = reference_level)
      dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = per_design)
      dds <- DESeq(dds)
      conds <- levels(dds$condition)
      pairs <- combn(conds, 2, simplify = FALSE)
      res_list <- setNames(lapply(pairs, function(p) {
        coef <- c(condition_col, p[1], p[2])
        out  <- get_tidy_res(dds, coef)
        attr(out, "contrast") <- paste(p[1], "vs", p[2])
        out
      }), vapply(pairs, function(p) paste(p[1], "vs", p[2]), ""))
      res_list
    }
    if (!is.null(BPPARAM)) {
      library(BiocParallel)
      result <- bplapply(names(pbp), run_one, BPPARAM = BPPARAM)
    } else {
      result <- lapply(names(pbp), run_one)
    }
    names(result) <- names(pbp)
    return(result)
  }
  
  # 5b. Combined mode analysis
  dds <- DESeqDataSetFromMatrix(countData = combined_counts, colData = combined_meta, design = comb_design)
  dds <- DESeq(dds)
  # boost maxit for convergence
  dds <- nbinomWaldTest(dds, maxit = 50)
  resnames     <- resultsNames(dds)
  
  # explicitly build the primary condition contrast name
  numerator   <- levels(combined_meta$condition)[2]
  denominator <- reference_level
  cond_contrast <- paste0(condition_col, "_", numerator, "_vs_", denominator)
  if (!(cond_contrast %in% resnames)) {
    stop("Primary condition contrast not found: ", cond_contrast)
  }
  
  # Build results per cluster
  results_list <- setNames(lapply(levels(combined_meta$cluster), function(cl) {
    if (cl == levels(combined_meta$cluster)[1]) {
      # baseline cluster uses main effect only
      get_tidy_res(dds, coef = cond_contrast)
    } else {
      # build interaction term name
      int_name <- paste0("cluster", cl, ".condition", numerator, "_vs_", denominator)
      if (!(int_name %in% resnames)) {
        stop("Interaction term not found for cluster ", cl, ": ", int_name)
      }
      # sum main effect + interaction
      get_tidy_res(dds, coef = list(cond_contrast, int_name))
    }
  }), levels(combined_meta$cluster))
  
  return(results_list)
}




#’ Prepare a DESeq2‐style colData table from a Seurat object
#’
#’ @param seurat_obj     A Seurat v5 object
#’ @param sample_col     Name of the metadata column holding your sample IDs (e.g. “orig.ident”)
#’ @param condition_col  Name of the metadata column holding your condition (e.g. “condition”)
#’ @param extra_cols     Optional character vector of additional metadata columns to carry along
#’ @param reference      Optional character: level of `condition_col` to use as the reference
#’ @return A data.frame whose rownames are sample IDs and columns are condition (and any `extra_cols`), with `condition` as a factor
prepare_colData <- function(seurat_obj,
                            sample_col    = "orig.ident",
                            condition_col = "condition",
                            extra_cols    = NULL,
                            reference     = NULL) {
  
  # 1) Grab the raw metadata
  meta <- as.data.frame(seurat_obj@meta.data)
  
  # 2) Check requested columns exist
  req_cols <- unique(c(sample_col, condition_col, extra_cols))
  missing  <- setdiff(req_cols, colnames(meta))
  if (length(missing)) {
    stop("These columns are missing from @meta.data: ", paste(missing, collapse = ", "))
  }
  
  # 3) Subset and dedupe
  colData <- meta[, req_cols, drop = FALSE]
  colData <- unique(colData)  # one row per sample
  
  # 4) Use sample IDs as rownames, then drop that column
  rownames(colData) <- colData[[sample_col]]
  colData[[sample_col]] <- NULL
  
  # 5) Turn condition into a factor (and re‐level if requested)
  colData[[condition_col]] <- factor(colData[[condition_col]])
  if (!is.null(reference)) {
    colData[[condition_col]] <- relevel(colData[[condition_col]], ref = reference)
  }
  
  return(colData)
}

#’ Build DESeqDataSet objects for pseudobulk profiles
#’
#’ @param pbp Named list of pseudobulk count matrices (genes × samples) per cluster
#’ @param colData Sample metadata data.frame; rownames match sample IDs
#’ @param condition_col Column name in colData for condition factor
#’ @param cluster_mode "perCluster" or "combined"
#’ @param design Optional DESeq2 design formula; default uses condition or cluster+condition+interaction
#’ @param reference_level Reference level for condition (default "WT")
#’ @param min_reps Minimum replicates per condition (default 2)
#’ @param filter_min_counts Genes with total counts < this are filtered (default 10)
#’ @return A list:
#’   - If perCluster: named list of DESeqDataSet, one per cluster
#’   - If combined: one DESeqDataSet for all clusters
#’ @importFrom DESeq2 DESeqDataSetFromMatrix DESeq nbinomWaldTest
build_pseudobulk_dds <- function(pbp,
                                 colData,
                                 condition_col,
                                 cluster_mode      = c("perCluster","combined"),
                                 design            = NULL,
                                 reference_level   = "WT",
                                 min_reps          = 2,
                                 filter_min_counts = 10) {
  library(DESeq2)
  cluster_mode <- match.arg(cluster_mode)
  
  # set design formulas
  if (is.null(design)) {
    per_design  <- as.formula(paste("~", condition_col))
    comb_design <- ~ cluster + condition + cluster:condition
  } else {
    per_design  <- design
    comb_design <- design
  }
  
  # per-cluster: build one DDS per cluster
  if (cluster_mode == "perCluster") {
    dds_list <- lapply(names(pbp), function(cl) {
      mat  <- pbp[[cl]]
      # subset samples and genes
      mat  <- mat[rowSums(mat) >= filter_min_counts, , drop=FALSE]
      ok   <- colSums(mat) > 0
      mat  <- mat[, ok, drop=FALSE]
      meta <- colData[colnames(mat), , drop=FALSE]
      meta$condition <- relevel(factor(meta[[condition_col]]), ref=reference_level)
      dds <- DESeqDataSetFromMatrix(mat, meta, design = per_design)
      DESeq(dds)
    })
    names(dds_list) <- names(pbp)
    return(dds_list)
  }
  
  # combined: merge matrices, build one DDS
  common_genes <- Reduce(intersect, lapply(pbp, rownames))
  pbp         <- lapply(pbp, function(m) m[common_genes, , drop=FALSE])
  # prefix cluster and bind
  mats <- mapply(function(m, cl) {
    colnames(m) <- paste(cl, colnames(m), sep="__")
    m
  }, pbp, names(pbp), SIMPLIFY=FALSE)
  combined_counts <- do.call(cbind, mats)
  # build combined metadata
  parts <- do.call(rbind, strsplit(colnames(combined_counts), "__", fixed=TRUE))
  samples <- parts[,2]; clusters <- parts[,1]
  meta <- colData[samples, , drop=FALSE]
  rownames(meta) <- colnames(combined_counts)
  meta$cluster   <- factor(clusters, levels = names(pbp))
  meta$condition <- relevel(factor(meta[[condition_col]]), ref = reference_level)
  # filter genes
  keep <- rowSums(combined_counts) >= filter_min_counts
  combined_counts <- combined_counts[keep, , drop=FALSE]
  
  dds <- DESeqDataSetFromMatrix(combined_counts, meta, design = comb_design)
  dds <- DESeq(dds)
  # increase iterations
  nbinomWaldTest(dds, maxit = 50)
}
library(DESeq2)
library(tibble)
library(ashr)
library(DESeq2)
library(tibble)
library(ashr)

extract_pseudobulk_results <- function(dds_list,
                                       condition_col,
                                       cluster_mode      = c("perCluster","combined"),
                                       lfc_shrink        = TRUE,
                                       shrink_type       = c("apeglm","ashr","normal")) {
  cluster_mode <- match.arg(cluster_mode)
  shrink_type  <- match.arg(shrink_type)
  
  tidy_dds <- function(dds) {
    rn       <- resultsNames(dds)
    lv       <- levels(dds[[condition_col]])
    ref      <- lv[1]
    non_ref  <- lv[-1]
    
    # 1) main-effect contrasts
    main_contrasts <- setNames(
      lapply(non_ref, function(l) c(condition_col, l, ref)),
      paste0("cond_", non_ref, "_vs_", ref)
    )
    
    results <- list()
    for (nm in names(main_contrasts)) {
      ctr  <- main_contrasts[[nm]]
      res  <- results(dds, contrast = ctr)
      if (lfc_shrink && shrink_type != "normal") {
        coef_name <- paste0(condition_col, "_", ctr[2], "_vs_", ctr[3])
        res        <- lfcShrink(dds, coef = coef_name, type = shrink_type)
      }
      results[[nm]] <- as_tibble(rownames_to_column(as.data.frame(res), "gene"))
    }
    
    # 2) combined‐mode: add cluster‐specific sums
    if (cluster_mode=="combined" && "cluster" %in% names(colData(dds))) {
      cls <- levels(dds$cluster)
      for (cl in cls[-1]) {
        for (lvl in non_ref) {
          main_name <- paste0(condition_col, "_", lvl, "_vs_", ref)
          int_name  <- paste0("cluster", cl, ".condition", lvl)
          if (!all(c(main_name, int_name) %in% rn)) {
            stop("Cannot find coefficients: ", main_name, " or ", int_name)
          }
          # ── build as TWO character vectors of names ──
          contrast_list <- list(
            c(main_name),  # group-1: main effect
            c(int_name)    # group-2: interaction effect
          )
          # inside tidy_dds(), combined-mode branch:
          rn <- resultsNames(dds)
          
          # Loop over clusters & levels...
          weights <- setNames(numeric(length(rn)), rn)
          weights[main_name] <- 1
          weights[int_name]  <- 1
          
          res_comb <- results(dds, contrast = weights)
          
          if (lfc_shrink && shrink_type == "ashr") {
            ash_out                 <- ash(res_comb$log2FoldChange, res_comb$lfcSE)
            res_comb$log2FoldChange <- ash_out$PosteriorMean
            res_comb$lfcSE          <- ash_out$PosteriorSD
          }
          
          nm <- paste0("cl_", cl, "_", lvl, "_vs_", ref)
          results[[nm]] <- as_tibble(rownames_to_column(as.data.frame(res_comb), "gene"))
          
        }
      }
    }
    
    results
  }
  
  if (cluster_mode == "perCluster") {
    lapply(dds_list, tidy_dds)
  } else {
    tidy_dds(dds_list)
  }
}

extract_pseudobulk_results_wc <- function(dds_list,
                                          condition_col,
                                          cluster_mode      = c("perCluster","combined"),
                                          lfc_shrink        = TRUE,
                                          shrink_type       = c("apeglm","ashr","normal")) {
  cluster_mode <- match.arg(cluster_mode)
  shrink_type  <- match.arg(shrink_type)
  
  tidy_dds <- function(dds) {
    rn      <- resultsNames(dds)
    lv      <- levels(dds[[condition_col]])
    ref     <- lv[1]
    non_ref <- lv[-1]
    
    results <- list()
    
    # 1) main-effect contrasts
    for (lvl in non_ref) {
      ctr_name <- paste0("cond_", lvl, "_vs_", ref)
      res      <- results(dds, contrast = c(condition_col, lvl, ref))
      if (lfc_shrink && shrink_type != "normal") {
        coef_nm <- paste0(condition_col, "_", lvl, "_vs_", ref)
        res     <- lfcShrink(dds, coef = coef_nm, type = shrink_type)
      }
      results[[ctr_name]] <- as_tibble(rownames_to_column(as.data.frame(res), "gene"))
    }
    
    # 2) combined-mode: cluster-specific sums via contraster()
    if (cluster_mode == "combined" && "cluster" %in% names(colData(dds))) {
      cls <- levels(dds$cluster)
      for (cl in cls) {
        for (lvl in non_ref) {
          # Build the numeric contrast for cluster=cl, condition=lvl vs ref
          grp1 <- list(c("cluster", cl), c(condition_col, lvl))
          grp2 <- list(c("cluster", cl), c(condition_col, ref))
          weights <- contraster(dds, group1 = grp1, group2 = grp2, weighted = FALSE)
          
          # Extract results for this cluster–condition pair
          res_comb <- results(dds, contrast = weights)
          if (lfc_shrink && shrink_type == "ashr") {
            ash_out                 <- ash(res_comb$log2FoldChange, res_comb$lfcSE)
            res_comb$log2FoldChange <- ash_out$PosteriorMean
            res_comb$lfcSE          <- ash_out$PosteriorSD
          }
          
          nm <- paste0("cl_", cl, "_", lvl, "_vs_", ref)
          results[[nm]] <- as_tibble(rownames_to_column(as.data.frame(res_comb), "gene"))
        }
      }
    }
    results
  }
  
  if (cluster_mode == "perCluster") {
    lapply(dds_list, tidy_dds)
  } else {
    tidy_dds(dds_list)
  }
}


extract_pseudobulk_results_all_pairs <- function(dds_list,
                                                 condition_col,
                                                 cluster_mode      = c("perCluster","combined"),
                                                 shrink_type       = c("ashr","normal")) {
  cluster_mode <- match.arg(cluster_mode)
  shrink_type  <- match.arg(shrink_type)
  
  tidy_dds <- function(dds) {
    lv         <- levels(colData(dds)[[condition_col]])
    cond_pairs <- combn(lv, 2, simplify = FALSE)
    results    <- list()
    
    # 1) Global pairwise contrasts
    for (pair in cond_pairs) {
      a <- pair[1]; b <- pair[2]
      nm <- paste0("cond_", a, "_vs_", b)
      
      # raw results
      res_raw <- results(dds, contrast = c(condition_col, a, b))
      
      # shrink with ASHR (supports any contrast)
      if (shrink_type == "ashr") {
        res <- lfcShrink(dds,
                         contrast = c(condition_col, a, b),
                         res      = res_raw,
                         type     = "ashr")
      } else {
        res <- lfcShrink(dds,
                         contrast = c(condition_col, a, b),
                         res      = res_raw,
                         type     = "normal")
      }
      
      results[[nm]] <- as_tibble(rownames_to_column(as.data.frame(res), "gene"))
    }
    
    # 2) Cluster-specific pairwise contrasts
    if (cluster_mode == "combined" && "cluster" %in% names(colData(dds))) {
      cls <- levels(dds$cluster)
      for (cl in cls) {
        for (pair in cond_pairs) {
          a <- pair[1]; b <- pair[2]
          nm <- paste0("cl_", cl, "_", a, "_vs_", b)
          
          # numeric contrast weights via contraster()
          weights <- contraster(dds,
                                group1   = list(c("cluster", cl), c(condition_col, a)),
                                group2   = list(c("cluster", cl), c(condition_col, b)),
                                weighted = FALSE)
          
          # raw results
          res_raw <- results(dds, contrast = weights)
          
          # ASHR shrink
          if (shrink_type == "ashr") {
            res <- lfcShrink(dds,
                             contrast = weights,
                             res      = res_raw,
                             type     = "ashr")
          } else {
            res <- lfcShrink(dds,
                             contrast = weights,
                             res      = res_raw,
                             type     = "normal")
          }
          
          results[[nm]] <- as_tibble(rownames_to_column(as.data.frame(res), "gene"))
        }
      }
    }
    
    results
  }
  
  if (cluster_mode == "perCluster") {
    lapply(dds_list, tidy_dds)
  } else {
    tidy_dds(dds_list)
  }
}
#’ Run DESeq2 on pseudobulk profiles, per‐cluster only
#’
#’ @param pbp Named list of pseudobulk count matrices (genes × samples), one entry per cluster
#’ @param colData Data.frame of sample metadata; rownames must match colnames of each matrix
#’ @param condition_col String; name of column in colData giving the condition factor
#’ @param design Optional formula for DESeq2 design; defaults to ~ condition_col
#’ @param reference_level String; reference level for the condition factor (default "WT")
#’ @param filter_min_counts Integer; drop genes with fewer than this total count (default 10)
#’ @param lfc_shrink Logical; whether to apply DESeq2’s lfcShrink() to results (default TRUE)
#’ @param shrink_type String, one of "apeglm", "ashr", or "normal" (default "apeglm")
#’ @param BPPARAM BiocParallelParam (optional) for parallel execution
#’ @return Named list (per cluster) of lists of tibbles, one tibble per pairwise contrast
#’ @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink resultsNames
#’ @importFrom tibble rownames_to_column as_tibble
#’ @importFrom dplyr relevel
#’ @importFrom BiocParallel bplapply
run_pseudobulk_deseq2_perCluster <- function(pbp,
                                             colData,
                                             condition_col,
                                             design            = NULL,
                                             reference_level   = "WT",
                                             filter_min_counts = 10,
                                             lfc_shrink        = TRUE,
                                             shrink_type       = c("apeglm","ashr","normal"),
                                             BPPARAM           = NULL) {
  library(DESeq2); library(dplyr); library(tibble)
  shrink_type <- match.arg(shrink_type)
  
  # Validate inputs
  stopifnot(is.list(pbp), is.data.frame(colData))
  stopifnot(condition_col %in% colnames(colData))
  
  # Set design
  per_design <- if (is.null(design)) 
    as.formula(paste("~", condition_col)) 
  else 
    design
  
  # Helper: extract and shrink results
  get_tidy_res <- function(dds, coef) {
    # Single coefficient vs contrast
    if (length(coef) == 1 && coef %in% resultsNames(dds)) {
      res_obj <- if (lfc_shrink) 
        lfcShrink(dds, coef = coef, type = shrink_type) 
      else 
        results(dds, name = coef)
    } else {
      res_obj <- results(dds, contrast = coef)
      if (lfc_shrink && shrink_type == "ashr") {
        library(ashr)
        ash_out <- ash(res_obj$log2FoldChange, res_obj$lfcSE)
        res_obj$log2FoldChange <- ash_out$PosteriorMean
        res_obj$lfcSE          <- ash_out$PosteriorSD
      }
    }
    df <- as.data.frame(res_obj, row.names = rownames(res_obj))
    as_tibble(rownames_to_column(df, var = "gene"))
  }
  
  # Per‐cluster analysis
  run_one <- function(cl) {
    mat <- pbp[[cl]][, intersect(colnames(pbp[[cl]]), rownames(colData)), drop = FALSE]
    mat <- mat[rowSums(mat) >= filter_min_counts, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    
    meta <- colData[colnames(mat), , drop = FALSE]
    meta$condition <- relevel(factor(meta[[condition_col]]), ref = reference_level)
    
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData   = meta,
                                  design    = per_design)
    dds <- DESeq(dds)
    
    conds <- levels(dds$condition)
    pairs <- combn(conds, 2, simplify = FALSE)
    
    # Extract each pairwise contrast
    res_list <- setNames(lapply(pairs, function(p) {
      contrast_vec <- c(condition_col, p[1], p[2])
      out <- get_tidy_res(dds, contrast_vec)
      attr(out, "contrast") <- paste(p[1], "vs", p[2])
      out
    }), vapply(pairs, function(p) paste(p[1], "vs", p[2]), ""))
    
    res_list
  }
  
  # Execute per‐cluster (optionally parallel)
  if (!is.null(BPPARAM)) {
    library(BiocParallel)
    results <- bplapply(names(pbp), run_one, BPPARAM = BPPARAM)
  } else {
    results <- lapply(names(pbp), run_one)
  }
  names(results) <- names(pbp)
  results
}
#’ Run DESeq2 on pseudobulk profiles, per-cluster only (returning DDS and results)
#’
#’ @param pbp Named list of pseudobulk count matrices (genes × samples), one entry per cluster
#’ @param colData Data.frame of sample metadata; rownames must match colnames of each matrix
#’ @param condition_col String; name of column in colData giving the condition factor
#’ @param design Optional formula for DESeq2 design; defaults to ~ condition_col
#’ @param reference_level String; reference level for the condition factor (default "WT")
#’ @param filter_min_counts Integer; drop genes with fewer than this total count (default 10)
#’ @param lfc_shrink Logical; whether to apply DESeq2’s lfcShrink() to results (default TRUE)
#’ @param shrink_type String, one of "apeglm", "ashr", or "normal" (default "apeglm")
#’ @param BPPARAM BiocParallelParam (optional) for parallel execution
#’ @return Named list (per cluster) of lists, each with components:
#’   - dds: the DESeqDataSet after DESeq()
#’   - results: a named list of tibbles, one per pairwise contrast
#’ @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink resultsNames
#’ @importFrom tibble rownames_to_column as_tibble
#’ @importFrom dplyr relevel
#’ @importFrom BiocParallel bplapply
run_pseudobulk_deseq2_perCluster_withDDS <- function(pbp,
                                                     colData,
                                                     condition_col,
                                                     design            = NULL,
                                                     reference_level   = "WT",
                                                     filter_min_counts = 10,
                                                     lfc_shrink        = TRUE,
                                                     shrink_type       = c("apeglm","ashr","normal"),
                                                     BPPARAM           = NULL) {
  library(DESeq2)
  library(dplyr)
  library(tibble)
  shrink_type <- match.arg(shrink_type)
  
  # Validate inputs
  stopifnot(is.list(pbp), is.data.frame(colData))
  stopifnot(condition_col %in% colnames(colData))
  
  # Set design
  per_design <- if (is.null(design)) {
    as.formula(paste("~", condition_col))
  } else {
    design
  }
  
  # Helper: extract and shrink results
  get_tidy_res <- function(dds, coef) {
    if (length(coef) == 1 && coef %in% resultsNames(dds)) {
      res_obj <- if (lfc_shrink) {
        lfcShrink(dds, coef = coef, type = shrink_type)
      } else {
        results(dds, name = coef)
      }
    } else {
      res_obj <- results(dds, contrast = coef)
      if (lfc_shrink && shrink_type == "ashr") {
        library(ashr)
        ash_out <- ash(res_obj$log2FoldChange, res_obj$lfcSE)
        res_obj$log2FoldChange <- ash_out$PosteriorMean
        res_obj$lfcSE          <- ash_out$PosteriorSD
      }
    }
    as_tibble(rownames_to_column(as.data.frame(res_obj), var = "gene"))
  }
  
  # Per-cluster analysis: returns a list with components 'dds' and 'results'
  run_one <- function(cl) {
    mat <- pbp[[cl]][, intersect(colnames(pbp[[cl]]), rownames(colData)), drop = FALSE]
    mat <- mat[rowSums(mat) >= filter_min_counts, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    
    meta <- colData[colnames(mat), , drop = FALSE]
    meta$condition <- relevel(factor(meta[[condition_col]]), ref = reference_level)
    
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData   = meta,
                                  design    = per_design)
    dds <- DESeq(dds)
    
    conds <- levels(dds$condition)
    pairs <- combn(conds, 2, simplify = FALSE)
    
    res_list <- setNames(lapply(pairs, function(p) {
      contrast_vec <- c(condition_col, p[1], p[2])
      out <- get_tidy_res(dds, contrast_vec)
      attr(out, "contrast") <- paste(p[1], "vs", p[2])
      out
    }), vapply(pairs, function(p) paste(p[1], "vs", p[2]), ""))
    
    list(dds = dds, results = res_list)
  }
  
  # Execute per-cluster (optionally parallel)
  if (!is.null(BPPARAM)) {
    library(BiocParallel)
    out <- bplapply(names(pbp), run_one, BPPARAM = BPPARAM)
  } else {
    out <- lapply(names(pbp), run_one)
  }
  names(out) <- names(pbp)
  
  out
}


contraster <- function(dds,    # A DESeqDataSet object containing colData (sample metadata) and a design formula
                       group1, # A list of character vectors; each vector specifies a column in colData and the values defining group1
                       group2, # A list of character vectors; each vector specifies a column in colData and the values defining group2
                       weighted = FALSE # Boolean flag: If FALSE, remove duplicate rows from model matrix; if TRUE, keep them
){
  
  # Step 1: Generate the design model matrix based on the experimental design
  # This matrix encodes the statistical design as numerical values, with columns representing different factor levels.
  design_matrix <- model.matrix(design(dds), colData(dds))
  
  # Step 2: Initialize lists to store boolean masks (TRUE/FALSE) indicating which samples belong to each group
  group1_filter_mask <- list()
  group2_filter_mask <- list()
  
  # Step 3: Generate filtering masks for group1
  # Iterate over each filtering condition in group1 and create a logical vector.
  for(i in 1:length(group1)){
    
    # Extract the column name specified in group1[[i]][1] and filter based on the provided values.
    # The result is a logical vector where TRUE indicates the sample belongs to group1.
    group1_filter_mask[[i]] <- colData(dds)[[group1[[i]][1]]] %in% group1[[i]][2:length(group1[[i]])]
    
  }
  
  # Step 4: Generate filtering masks for group2
  for(i in 1:length(group2)){
    
    # Similar to the logic for group1, extract column name and check for matching values.
    group2_filter_mask[[i]] <- colData(dds)[[group2[[i]][1]]] %in% group2[[i]][2:length(group2[[i]])]
    
  }
  
  # Step 5: Combine multiple conditions within each group using logical AND (`&`)
  # This ensures that samples must satisfy **all** conditions in the list.
  group1_filter_mask <- Reduce(function(x, y) x & y, group1_filter_mask)
  group2_filter_mask <- Reduce(function(x, y) x & y, group2_filter_mask)
  
  # Step 6: Subset the design matrix to extract rows corresponding to each group
  # This filters the design matrix to only include samples that belong to group1 or group2.
  group1_design_matrix <- design_matrix[group1_filter_mask, , drop=FALSE]
  group2_design_matrix <- design_matrix[group2_filter_mask, , drop=FALSE]
  
  # Step 7: Remove duplicate rows if `weighted = FALSE`
  # Duplicate rows can arise due to technical replicates or redundant samples.
  if(!weighted){
    
    group1_design_matrix <- group1_design_matrix[!duplicated(group1_design_matrix), , drop=FALSE]
    group2_design_matrix <- group2_design_matrix[!duplicated(group2_design_matrix), , drop=FALSE]
    
  }
  
  # Step 8: Compute the contrast (difference between mean values in the two groups)
  # This provides a numerical difference between the model matrix means of group1 and group2.
  return(colMeans(group1_design_matrix) - colMeans(group2_design_matrix))
  
}
