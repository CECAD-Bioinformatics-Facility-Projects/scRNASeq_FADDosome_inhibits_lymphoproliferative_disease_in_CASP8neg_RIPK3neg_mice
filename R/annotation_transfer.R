# Transfer cell_type / cell_type.fine annotation from an original (published)
# Seurat object onto a freshly clustered object, without assuming cluster IDs
# are stable across reruns.
#
# Strategy (cell-first, cluster-fallback):
#   1. For every barcode that exists in both objects, copy cell_type and
#      cell_type.fine verbatim from the original (hand-curated ground truth).
#   2. For cells present only in the new object, use majority-vote from
#      shared cells in the same (sub-)cluster. Sub-clustering is applied only
#      to clusters where shared-cell purity is below `purity_threshold`, so
#      the fallback vote is drawn from a more homogeneous neighbourhood.
#
# Used by both relabel_from_original.R (standalone relabeling) and
# initial_analysis_updated.R (in-pipeline annotation of a fresh rerun).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

transfer_annotation_from_original <- function(new, orig, cluster_col,
                                              purity_threshold = 0.90,
                                              min_shared        = 30,
                                              sub_resolutions   = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.2),
                                              graph_name        = "integrated_snn") {
  if (!"cell_type" %in% colnames(orig@meta.data) ||
      !"cell_type.fine" %in% colnames(orig@meta.data)) {
    stop("Original object missing cell_type / cell_type.fine.", call. = FALSE)
  }
  if (!cluster_col %in% colnames(new@meta.data)) {
    stop("cluster_col '", cluster_col, "' not found in new object.", call. = FALSE)
  }

  common      <- intersect(colnames(orig), colnames(new))
  only_in_new <- setdiff(colnames(new), common)
  message(sprintf("Shared barcodes: %d / %d new (%.2f%%)",
                  length(common), ncol(new), 100 * length(common) / ncol(new)))
  message(sprintf("Cells present only in the new object: %d (fallback needed)",
                  length(only_in_new)))

  orig_bt   <- setNames(as.character(orig$cell_type),      colnames(orig))
  orig_fine <- setNames(as.character(orig$cell_type.fine), colnames(orig))

  DefaultAssay(new) <- "RNA"
  parent_ids <- as.character(new@meta.data[[cluster_col]])
  names(parent_ids) <- colnames(new)

  cluster_purity <- function(cells) {
    shared <- intersect(cells, common)
    n <- length(shared)
    if (!n) return(list(purity = NA_real_, n_shared = 0L,
                        label_broad = NA_character_, label_fine = NA_character_))
    tb <- sort(table(orig_bt[shared]),   decreasing = TRUE)
    tf <- sort(table(orig_fine[shared]), decreasing = TRUE)
    list(purity      = unname(tb[1]) / n,
         n_shared    = n,
         label_broad = names(tb)[1],
         label_fine  = names(tf)[1])
  }

  subcluster_ids <- function(so, target, resolution) {
    tmp_col <- paste0("__sub__", target, "__", resolution)
    Idents(so) <- cluster_col
    so <- FindSubCluster(so, cluster = target, resolution = resolution,
                         graph.name = graph_name, subcluster.name = tmp_col)
    setNames(as.character(so@meta.data[[tmp_col]]), colnames(so))
  }

  baseline <- tibble(cluster = sort(unique(parent_ids))) |>
    rowwise() |>
    mutate(s = list(cluster_purity(names(parent_ids)[parent_ids == cluster]))) |>
    mutate(purity = s$purity, n_shared = s$n_shared,
           label_broad = s$label_broad, label_fine = s$label_fine) |>
    select(-s) |>
    ungroup() |>
    arrange(purity)

  message("\n=== Baseline shared-cell purity per new cluster ===")
  print(as.data.frame(baseline), row.names = FALSE)

  impure <- baseline$cluster[!is.na(baseline$purity) &
                             baseline$purity < purity_threshold &
                             baseline$n_shared >= min_shared]
  message(sprintf("\n%d cluster(s) below purity %.2f -> sub-clustering: %s",
                  length(impure), purity_threshold,
                  paste(impure, collapse = ", ")))

  final_ids <- parent_ids
  split_log <- list()
  for (cl in impure) {
    cells_cl <- names(parent_ids)[parent_ids == cl]
    chosen_res <- NA_real_; chosen_ids <- NULL
    for (res in sub_resolutions) {
      sub_ids <- subcluster_ids(new, target = cl, resolution = res)[cells_cl]
      sub_summary <- tibble(sub = unique(sub_ids)) |>
        rowwise() |>
        mutate(s = list(cluster_purity(names(sub_ids)[sub_ids == sub]))) |>
        mutate(purity = s$purity, n_shared = s$n_shared) |>
        select(-s) |>
        ungroup()
      ok <- all(is.na(sub_summary$purity) |
                sub_summary$n_shared < min_shared |
                sub_summary$purity >= purity_threshold)
      message(sprintf("  cluster %s @ res=%.2f -> %d pieces, min purity=%.2f, %s",
                      cl, res, nrow(sub_summary),
                      min(sub_summary$purity, na.rm = TRUE),
                      if (ok) "ACCEPT" else "try higher"))
      chosen_res <- res; chosen_ids <- sub_ids
      if (ok) break
    }
    final_ids[cells_cl] <- chosen_ids
    split_log[[cl]] <- list(resolution = chosen_res,
                            n_pieces   = length(unique(chosen_ids)))
  }

  # Per-cell label assignment: verbatim copy for shared cells; per-(sub-)cluster
  # majority-vote for the rest.
  cell_type_final      <- rep(NA_character_, ncol(new))
  cell_type_fine_final <- rep(NA_character_, ncol(new))
  label_source         <- rep(NA_character_, ncol(new))
  names(cell_type_final) <- names(cell_type_fine_final) <- names(label_source) <- colnames(new)

  cell_type_final[common]      <- orig_bt[common]
  cell_type_fine_final[common] <- orig_fine[common]
  label_source[common]         <- "copied_from_original"

  sub_summary_map <- tibble(cluster = unique(final_ids)) |>
    rowwise() |>
    mutate(s = list(cluster_purity(names(final_ids)[final_ids == cluster]))) |>
    mutate(fallback_broad = s$label_broad,
           fallback_fine  = s$label_fine,
           fallback_purity = s$purity,
           fallback_n_shared = s$n_shared) |>
    select(-s) |>
    ungroup()

  fb_broad <- setNames(sub_summary_map$fallback_broad, sub_summary_map$cluster)
  fb_fine  <- setNames(sub_summary_map$fallback_fine,  sub_summary_map$cluster)

  for (bc in only_in_new) {
    cl <- final_ids[[bc]]
    cell_type_final[bc]      <- fb_broad[[cl]]
    cell_type_fine_final[bc] <- fb_fine[[cl]]
    label_source[bc]         <- "cluster_majority_vote"
  }

  message(sprintf("\nAssignment source: %d copied, %d fallback",
                  sum(label_source == "copied_from_original"),
                  sum(label_source == "cluster_majority_vote")))

  new_col <- paste0(cluster_col, ".resplit")
  new@meta.data[[new_col]]        <- unname(final_ids[colnames(new)])
  new@meta.data$cell_type         <- unname(cell_type_final[colnames(new)])
  new@meta.data$cell_type.fine    <- unname(cell_type_fine_final[colnames(new)])
  new@meta.data$annotation_source <- unname(label_source[colnames(new)])

  comp <- tibble(cluster = new@meta.data[[new_col]],
                cell_type = new@meta.data$cell_type) |>
    dplyr::count(cluster, cell_type) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(pct = 100 * n / sum(n), total = sum(n)) |>
    dplyr::ungroup()

  list(object = new, cluster_col_resplit = new_col,
      composition = comp, split_log = split_log)
}
