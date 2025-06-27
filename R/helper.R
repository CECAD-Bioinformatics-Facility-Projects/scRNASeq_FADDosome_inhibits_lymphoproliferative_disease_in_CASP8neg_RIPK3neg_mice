dt_table <- function(data.frame, title="", filename="", rownames=F){
  
  table <- DT::datatable(as.data.frame(data.frame), rownames = rownames, filter = 'top', caption = title, extensions = 'Buttons',
                         options = list(dom = 'Bfrtip', scrollX=T,
                                        buttons =
                                          list(list(extend='copy',title=filename),
                                               list(extend='csv',      title=filename),
                                               list(extend='excel', filename=filename),
                                               list(extend='pdf',   filename=filename),
                                               list(extend='print', filename=filename)),
                                        autoWidth = F, pageLength = 10))
  return(table)
}

# make a function that removes axis lines, axis ticks and axis.text from ggplot
theme_void.2 <- function() {
  theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) + NoLegend()
}

# 
downsample_seurat_object <- function(seurat_object, sample_size=NULL){
  # Assuming `seurat_object` is your integrated Seurat object
  # and `sample` is a metadata column in `seurat_object@meta.data`
  # that identifies the sample each cell belongs to.
  
  # Step 1: Determine the fixed number of cells per sample
  
  # If `sample_size` is not provided, use the minimum number of cells per sample
  if(is.null(sample_size)) {
    sample_size <- min(table(seurat_object$orig.ident))
  }  
  
  fixed_cell_number <- sample_size
  
  # Initialize a vector to store cells to keep
  cells_to_keep <- c()
  
  # Step 2: Loop through each sample
  for(sample_name in unique(seurat_object$orig.ident)) {
    # Identify cells belonging to the current sample
    cells_in_sample <- rownames(seurat_object@meta.data[
      seurat_object@meta.data$orig.ident == sample_name,])
    
    # Randomly select cells to keep, if the sample has more cells than the fixed number
    if(length(cells_in_sample) > fixed_cell_number) {
      set.seed(42) # For reproducibility
      cells_to_keep <- c(cells_to_keep, sample(cells_in_sample, fixed_cell_number))
    } else {
      # If the sample has fewer than the fixed number, keep all cells
      cells_to_keep <- c(cells_to_keep, cells_in_sample)
    }
  }
  
  # Step 3: Subset the Seurat object to only keep the selected cells
  seurat_object_subset <- subset(seurat_object, cells = cells_to_keep)
}

#’ Run DESeq2 on pseudobulk profiles, per‐cluster only
#’
#’ @param pbp Named list of pseudobulk matrices (genes × samples), one entry per cluster
#’ @param colData Data.frame of sample metadata; rownames must match colnames of each matrix
#’ @param condition_col String; name of column in colData giving the condition factor
#’ @param design Optional formula for DESeq2 design; defaults to ~ condition_col
#’ @param reference_level String; reference level for the condition factor (default "WT")
#’ @param filter_min_counts Integer; drop genes with fewer than this total count (default 10)
#’ @param lfc_shrink Logical; whether to apply lfcShrink() (default TRUE)
#’ @param shrink_type String: "apeglm", "ashr", or "normal" (default "apeglm")
#’ @param BPPARAM BiocParallelParam for optional parallel execution
#’ @return Named list (per cluster) of lists of tibbles, one tibble per pairwise contrast
#’ @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink resultsNames
#’ @importFrom tibble rownames_to_column as_tibble
#’ @importFrom BiocParallel bplapply
run_pseudobulk_deseq_perCluster <- function(pbp,
                                            colData,
                                            condition_col,
                                            design            = NULL,
                                            reference_level   = "WT",
                                            filter_min_counts = 10,
                                            lfc_shrink        = TRUE,
                                            shrink_type       = c("apeglm","ashr","normal"),
                                            BPPARAM           = NULL) {
  library(DESeq2); library(tibble)
  shrink_type <- match.arg(shrink_type)
  
  # Validate inputs
  stopifnot(is.list(pbp), is.data.frame(colData))
  stopifnot(condition_col %in% colnames(colData))
  
  # Default design
  if (is.null(design)) {
    per_design <- as.formula(paste("~", condition_col))
  } else {
    per_design <- design
  }
  
  # Helper to extract and shrink results
  get_tidy_res <- function(dds, contrast_vec) {
    # raw results
    res_raw <- if (length(contrast_vec)==1 && contrast_vec %in% resultsNames(dds)) {
      results(dds, name = contrast_vec)
    } else {
      results(dds, contrast = contrast_vec)
    }
    
    # apply shrinkage
    if (lfc_shrink) {
      if (length(contrast_vec)==1 && shrink_type=="apeglm") {
        res_obj <- lfcShrink(dds, coef = contrast_vec, type = "apeglm")
      } else if (shrink_type=="ashr") {
        res_obj <- lfcShrink(dds,
                             contrast = contrast_vec,
                             res      = res_raw,
                             type     = "ashr")
      } else {
        res_obj <- lfcShrink(dds,
                             contrast = contrast_vec,
                             res      = res_raw,
                             type     = "normal")
      }
    } else {
      res_obj <- res_raw
    }
    
    # to tibble
    as_tibble(rownames_to_column(as.data.frame(res_obj), "gene"))
  }
  
  # Per‐cluster DESeq2
  run_one_cluster <- function(clust) {
    mat  <- pbp[[clust]]
    # keep only samples present in colData and with enough counts
    mat  <- mat[ , intersect(colnames(mat), rownames(colData)), drop=FALSE]
    mat  <- mat[rowSums(mat)>=filter_min_counts, , drop=FALSE]
    mat  <- mat[, colSums(mat)>0, drop=FALSE]
    
    meta <- colData[colnames(mat), , drop=FALSE]
    meta$condition <- relevel(factor(meta[[condition_col]]), ref = reference_level)
    
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData   = meta,
                                  design    = per_design)
    dds <- DESeq(dds)
    
    # all pairwise condition contrasts
    conds <- levels(dds$condition)
    pairs <- combn(conds, 2, simplify = FALSE)
    
    # extract results
    res_list <- lapply(pairs, function(p) {
      contrast_vec <- c(condition_col, p[1], p[2])
      out <- get_tidy_res(dds, contrast_vec)
      attr(out, "contrast") <- paste(p[1], "vs", p[2])
      out
    })
    names(res_list) <- vapply(pairs, function(p) paste(p[1], "vs", p[2]), "")
    
    res_list
  }
  
  # Run per‐cluster, optionally in parallel
  if (!is.null(BPPARAM)) {
    library(BiocParallel)
    results <- bplapply(names(pbp), run_one_cluster, BPPARAM = BPPARAM)
  } else {
    results <- lapply(names(pbp), run_one_cluster)
  }
  names(results) <- names(pbp)
  
  results
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