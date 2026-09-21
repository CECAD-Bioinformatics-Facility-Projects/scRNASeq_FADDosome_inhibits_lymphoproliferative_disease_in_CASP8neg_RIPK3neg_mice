# Canonical marker-based cell-type annotation.
#
# Assigns each cluster the identity whose canonical marker panel it expresses
# most strongly, then refines the T / NK compartment for the fine annotation.
# Because the label of a cluster is decided by *marker expression* rather than
# by its Louvain ID, the result does not depend on run-specific cluster
# numbering -- readers who re-run the pipeline get the same labels even when
# Seurat happens to number the clusters differently.
#
# This is the single source of truth for the marker panels: it is used both to
# ASSIGN labels (reproduce_seurat_object.R) and to VERIFY them
# (validate_annotation.R), so annotation and validation can never drift apart.
#
# The broad panels and the map from panel -> cell type mirror the identities
# used in supplementary/initial_analysis.Rmd
# (broad_annotation_map / T_refined_annotation_map).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ---- marker panels ----------------------------------------------------------

# Broad cell types: each cluster is labelled by whichever panel it expresses
# most (mean fraction of cells with non-zero counts, averaged over the panel).
broad_marker_panels <- list(
  T         = c("Cd3d", "Cd3e", "Cd3g", "Cd2"),
  B         = c("Cd19", "Cd79a", "Cd79b", "Ms4a1", "Pax5"),
  NK        = c("Ncr1", "Klrb1c", "Nkg7", "Klre1"),
  DC        = c("Itgax", "H2-Ab1", "H2-Aa", "Flt3"),
  Mono      = c("Ly6c1", "Ly6c2", "Ccr2", "Cx3cr1", "Csf1r"),
  MF        = c("Adgre1", "Cd68", "Marco"),
  Neutro    = c("S100a8", "S100a9", "Ly6g", "Mpo", "Elane"),
  Erythroid = c("Hba-a1", "Hbb-bs", "Hbb-bt", "Gypa", "Klf1")
)

# Extra markers used only to split the T / NK compartment for the fine labels.
fine_marker_sets <- list(
  Cd4    = c("Cd4"),
  Cd8    = c("Cd8a", "Cd8b1"),
  prolif = c("Mki67", "Top2a")
)

# Thresholds for the fine T / NK refinement (all on the 0-1 detection-rate
# scale, i.e. fraction of cells in the cluster with non-zero counts). They are
# deliberately exposed here so a reader can adjust them for their own data.
fine_thresholds <- list(
  prolif   = 0.40,  # a T cluster with >= this Mki67/Top2a detection -> "prolif. T"
  cd_min   = 0.10,  # if both Cd4 and Cd8 detection are below this -> "DN-T"
  nkt_tcr  = 0.25   # an NK-panel cluster with T-panel detection >= this -> "NKT"
)

# ---- internal helpers -------------------------------------------------------

# Return the RNA count matrix, tolerant of both Assay and Assay5 storage.
.annotation_counts <- function(so) {
  DefaultAssay(so) <- "RNA"
  if (inherits(so[["RNA"]], "Assay5")) {
    layer <- if ("counts" %in% Layers(so[["RNA"]])) "counts" else Layers(so[["RNA"]])[1]
  } else {
    layer <- "counts"
  }
  GetAssayData(so, assay = "RNA", layer = layer)
}

# Per-cluster mean detection rate of a set of genes.
.detect_by_cluster <- function(cnts, cluster, genes) {
  genes <- intersect(genes, rownames(cnts))
  if (!length(genes)) {
    return(setNames(rep(NA_real_, length(unique(cluster))), unique(cluster)))
  }
  per_cell <- Matrix::colMeans(cnts[genes, , drop = FALSE] > 0)
  tapply(per_cell, cluster, mean)
}

# ---- public API -------------------------------------------------------------

# Per-cluster mean detection for every broad panel (also used by validation).
panel_detection_by_cluster <- function(so, cluster_col,
                                        panels = broad_marker_panels) {
  cnts    <- .annotation_counts(so)
  cluster <- as.character(so@meta.data[[cluster_col]])
  out <- vapply(panels, function(m) .detect_by_cluster(cnts, cluster, m),
                FUN.VALUE = numeric(length(unique(cluster))))
  out <- as.data.frame(out)
  out$cluster <- rownames(out)
  rownames(out) <- NULL
  out[, c("cluster", names(panels))]
}

# cluster -> broad cell type (named character vector).
assign_broad_by_markers <- function(so, cluster_col,
                                     panels = broad_marker_panels) {
  byc <- panel_detection_by_cluster(so, cluster_col, panels)
  mat <- as.matrix(byc[, names(panels)])
  mat[is.na(mat)] <- -Inf
  labs <- names(panels)[max.col(mat, ties.method = "first")]
  setNames(labs, byc$cluster)
}

# cluster -> fine cell type, refining the T / NK compartment of a broad map.
assign_fine_by_markers <- function(so, cluster_col, broad_map,
                                    panels = broad_marker_panels,
                                    thresholds = fine_thresholds) {
  cnts    <- .annotation_counts(so)
  cluster <- as.character(so@meta.data[[cluster_col]])

  cd4    <- .detect_by_cluster(cnts, cluster, fine_marker_sets$Cd4)
  cd8    <- .detect_by_cluster(cnts, cluster, fine_marker_sets$Cd8)
  prolif <- .detect_by_cluster(cnts, cluster, fine_marker_sets$prolif)
  tcr    <- .detect_by_cluster(cnts, cluster, panels$T)

  fine <- broad_map
  for (cl in names(broad_map)) {
    lab <- broad_map[[cl]]
    if (lab == "T") {
      if (!is.na(prolif[[cl]]) && prolif[[cl]] >= thresholds$prolif) {
        fine[[cl]] <- "prolif. T"
      } else if (max(cd4[[cl]], cd8[[cl]], na.rm = TRUE) < thresholds$cd_min) {
        fine[[cl]] <- "DN-T"
      } else if (isTRUE(cd8[[cl]] > cd4[[cl]])) {
        fine[[cl]] <- "CD8+ T"
      } else {
        fine[[cl]] <- "CD4+ T"
      }
    } else if (lab == "NK") {
      fine[[cl]] <- if (!is.na(tcr[[cl]]) && tcr[[cl]] >= thresholds$nkt_tcr) "NKT" else "NK"
    }
  }
  fine
}

# Convenience wrapper: annotate an object in place and return it with
# `cell_type` and `cell_type.fine` columns filled from the cluster labels.
annotate_by_markers <- function(so, cluster_col,
                                panels = broad_marker_panels,
                                thresholds = fine_thresholds) {
  broad_map <- assign_broad_by_markers(so, cluster_col, panels)
  fine_map  <- assign_fine_by_markers(so, cluster_col, broad_map, panels, thresholds)
  cl <- as.character(so@meta.data[[cluster_col]])
  so@meta.data$cell_type      <- unname(broad_map[cl])
  so@meta.data$cell_type.fine <- unname(fine_map[cl])
  list(object = so, broad_map = broad_map, fine_map = fine_map)
}
