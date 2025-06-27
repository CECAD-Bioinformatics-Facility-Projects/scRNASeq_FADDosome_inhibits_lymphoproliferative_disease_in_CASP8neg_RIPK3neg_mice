plot_pheatmap <- function(expr_mat,
                          log_transform     = TRUE,
                          scale_rows        = TRUE,
                          cluster_rows      = TRUE,
                          cluster_cols      = TRUE,
                          color_palette     = colorRampPalette(
                            c("green","green","black",
                              "red",  "red"))(100),
                          border_color      = "grey60",
                          show_rownames     = TRUE,
                          show_colnames     = TRUE,
                          mark_nan          = TRUE,
                          include_nan_genes = TRUE,
                          gaps_col          = NULL,    # e.g. c(4, 6),
                          title
) {
  
  ## 1. Input checks
  if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) {
    stop("expr_mat must have row names (genes) and column names (samples).")
  }
  
  ## 2. Turn into matrix + optional log transform
  mat <- as.matrix(expr_mat)
  if (log_transform) mat <- log(mat + 1)
  
  ## 3. Row-scaling (Z-score), preserving NA
  if (scale_rows) {
    mat <- t(apply(mat, 1, function(x) {
      mu  <- mean(x, na.rm = TRUE)
      sdx <- sd(x,   na.rm = TRUE)
      (x - mu) / sdx
    }))
  }
  
  ## 4. Mark & zero-out non-finite cells
  number_matrix <- NULL
  if (mark_nan) {
    number_matrix <- matrix("", nrow = nrow(mat), ncol = ncol(mat),
                            dimnames = dimnames(mat))
    bad_idx       <- which(!is.finite(mat), arr.ind = TRUE)
    number_matrix[bad_idx] <- "X"
    mat[bad_idx]          <- 0
  }
  
  ## 5. Optionally drop genes that have any ✗
  if (!include_nan_genes && mark_nan) {
    # Keep only rows where number_matrix is all ""
    good_rows <- apply(number_matrix, 1, function(r) all(r == ""))
    mat        <- mat[good_rows, , drop = FALSE]
    number_matrix <- number_matrix[good_rows, , drop = FALSE]
  }
  
  ## 6. Plot
  library(pheatmap)
  ph <- pheatmap(
    mat,
    color            = color_palette,
    cluster_rows     = cluster_rows,
    cluster_cols     = cluster_cols,
    border_color     = border_color,
    show_rownames    = show_rownames,
    show_colnames    = show_colnames,
    display_numbers  = number_matrix,
    number_color     = "black",
    fontsize_number  = 12,
    fontsize_row     = 12,
    fontsize_col     = 12,
    treeheight_row   = 20,
    treeheight_col   = 20,
    gaps_col         = gaps_col,
    main             = title
  ) 
  
  ph |> ggplotify::as.ggplot()
  
  # invisible(rownames(mat))
}


plot_pp_heatmap <- function(pp,
                            cell_type             = "T",
                            selected_samples = NULL,
                            genes,
                            sample_order,
                            include_nan_genes = TRUE,
                            ...  # additional args to pass to plot_pheatmap()
) {
  #— 1. Extract the assay data.frame
  if(!is.null(selected_samples)){
    #— 1.1. Subset to the selected samples
    pp[[cell_type]] <- pp[[cell_type]][, selected_samples, drop = FALSE]
  } else {
    #— 1.2. Use all samples
    pp[[cell_type]] <- pp[[cell_type]][, colnames(pp), drop = FALSE]
  }
  df <- as.data.frame(pp[[cell_type]])
  
  #— 2. Filter to only the requested genes
  #    Use base subsetting to preserve rownames exactly :contentReference[oaicite:3]{index=3}
  df <- df[rownames(df) %in% genes, , drop = FALSE]
  
  #— 3. Safe reordering of genes
  #    Keep only genes present, in the user’s specified order :contentReference[oaicite:4]{index=4}
  good_genes <- genes[genes %in% rownames(df)]
  df <- df[good_genes, , drop = FALSE]
  
  #— 4. Safe reordering of samples
  #    Warn if any requested sample is missing :contentReference[oaicite:5]{index=5}
  missing_samps <- setdiff(sample_order, colnames(df))
  if (length(missing_samps)) {
    warning("These samples were not found and will be dropped: ",
            paste(missing_samps, collapse = ", "))
  }
  good_samps <- sample_order[sample_order %in% colnames(df)]
  df <- df[, good_samps, drop = FALSE]
  
  #— 5. Call your existing plot_pheatmap()
  #    Pass through additional parameters (e.g. include_nan_genes) :contentReference[oaicite:6]{index=6}
  plot_pheatmap(df,
                include_nan_genes = include_nan_genes,
                ...)
  
  #— 6. Return the final gene order invisibly
  #  invisible(rownames(df))
}

#' Plot pseudobulk heatmap from DESeq2-per-cluster results
#'
#' This wrapper takes the output of `run_pseudobulk_deseq2_perCluster_withDDS()`
#' and generates a heatmap for a given cluster, gene set, and sample order,
#' mirroring `plot_pp_heatmap()` behavior.
#'
#' @param pbpdds Named list-of-lists returned by `run_pseudobulk_deseq2_perCluster_withDDS()`.
#'   Each element must have a `$dds` component (DESeqDataSet).
#' @param cell_type String; name of the cluster to plot (must exist in `names(pbpdds)`).
#' @param genes Character vector of gene IDs to include.
#' @param sample_order Character vector specifying ordering of samples (columns).
#' @param include_nan_genes Logical; passed through to `plot_pheatmap()` to control
#'   inclusion of genes with all-NA values.
#' @param norm Logical; whether to use normalized counts (`counts(dds, normalized=TRUE)`)
#'   or raw counts (default `TRUE`).
#' @param ... Additional arguments passed to `plot_pheatmap()`.
#'
#' @return Invisibly, the ordered gene vector (rownames of the plotted matrix).
#' @importFrom DESeq2 counts
#' @importFrom SummarizedExperiment assay
#' @export
plot_pp_heatmap_from_pbpdds <- function(
    pbpdds,
    cell_type,
    genes,
    sample_order,
    include_nan_genes = TRUE,
    norm = TRUE,
    ...
) {
  if (!cell_type %in% names(pbpdds)) {
    stop("Cluster '", cell_type, "' not found in pseudobulk results")
  }
  
  # 1. Extract DESeqDataSet and counts matrix
  dds <- pbpdds[[cell_type]]$dds
  mat <- if (norm) {
    DESeq2::vst(dds, blind = FALSE) |> assay()
  } else {
    SummarizedExperiment::assay(dds)
  }
  df  <- as.data.frame(mat)
  
  # 2. Filter to requested genes
  df <- df[rownames(df) %in% genes, , drop = FALSE]
  
  # 3. Reorder genes safely
  good_genes <- genes[genes %in% rownames(df)]
  df <- df[good_genes, , drop = FALSE]
  
  # 4. Reorder samples safely, warn if missing
  missing_samps <- setdiff(sample_order, colnames(df))
  if (length(missing_samps)) {
    warning("These samples were not found and will be dropped: ",
            paste(missing_samps, collapse = ", "))
  }
  good_samps <- sample_order[sample_order %in% colnames(df)]
  df <- df[, good_samps, drop = FALSE]
  
  # 5. Plot via existing heatmap function
  plot_pheatmap(
    df,
    include_nan_genes = include_nan_genes, title="Robust normalization",
    ...
  )
}