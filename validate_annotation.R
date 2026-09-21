#!/usr/bin/env Rscript
# Sanity-check the annotation transferred onto the re-run object:
#
#   1. Per-cluster label composition + purity (does each cluster carry a
#      single dominant cell_type / cell_type.fine?).
#   2. Canonical marker-panel enrichment: each cluster's mean detection
#      frequency of each panel should be highest for the panel that matches
#      its assigned label (broad type).
#   3. Summary plot combining per-cluster purity and marker-panel agreement,
#      for a single at-a-glance read of the validation result.
#   4. UMAP plots: one coloured by cluster ID, one by cell_type, one by
#      cell_type.fine.
#
# Run from the project root:
#   Rscript validate_annotation.R
#
# Env vars:
#   NEW_QS      annotated Seurat object to validate (default:
#               results/rerun_20260916_origdoublets/seurat_objects.combined.cleansed.annotated.250428.qs)
#   OUT_DIR     output directory (default: results/annotation_validation)
#   CLUSTER_COL cluster column to use (default:
#               integrated_snn_res.0.1_and_cluster_4_and_cluster_9.resplit,
#               falling back to integrated_snn_res.0.1_and_cluster_4_and_cluster_9)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(qs)
  library(Matrix)
})

# replace the reproduced_20260921_083232 with the folder generated in you run.

NEW_QS      <- Sys.getenv("NEW_QS",
                          unset = "results/reproduced_20260921_083232/seurat_objects.combined.cleansed.annotated.qs")
OUT_DIR     <- Sys.getenv("OUT_DIR", unset = "results/reproduced_20260921_083232/annotation_validation")
CLUSTER_COL <- Sys.getenv("CLUSTER_COL", unset = "")

if (!file.exists(NEW_QS)) stop("NEW_QS does not exist: ", NEW_QS, call. = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading: ", NEW_QS)
so <- qs::qread(NEW_QS)

# Pick the cluster column: prefer the .resplit produced by the transfer step.
candidates <- c(
  CLUSTER_COL,
  "integrated_snn_res.0.1_and_cluster_4_and_cluster_9.resplit",
  "integrated_snn_res.0.1_and_cluster_4_and_cluster_9",
  "seurat_clusters"
)
candidates <- candidates[nzchar(candidates)]
cluster_col <- candidates[candidates %in% colnames(so@meta.data)][1]
if (is.na(cluster_col)) stop("No usable cluster column found.", call. = FALSE)
message("Using cluster column: ", cluster_col)

meta <- so@meta.data
meta$cluster <- as.character(meta[[cluster_col]])

# --- 1. Per-cluster label composition + purity ------------------------------



comp <- meta %>%
  dplyr::count(cluster, cell_type, cell_type.fine) %>%
  group_by(cluster) %>%
  mutate(share = n / sum(n)) %>%
  arrange(cluster, desc(n)) %>%
  ungroup()

top_per_cluster <- comp %>%
  group_by(cluster) %>%
  summarise(
    n_cells         = sum(n),
    top_broad       = cell_type[which.max(n)],
    top_broad_share = max(tapply(n, cell_type, sum)) / sum(n),
    top_fine        = cell_type.fine[which.max(n)],
    top_fine_share  = max(tapply(n, cell_type.fine, sum)) / sum(n),
    .groups         = "drop"
  ) %>%
  arrange(top_broad, cluster)

write.csv(comp,           file.path(OUT_DIR, "cluster_label_composition.csv"), row.names = FALSE)
write.csv(top_per_cluster, file.path(OUT_DIR, "cluster_top_labels.csv"),        row.names = FALSE)

cat("\n=== Cluster label summary (top broad / fine) ===\n")
print(as.data.frame(top_per_cluster), row.names = FALSE)

impure_b <- top_per_cluster[top_per_cluster$top_broad_share < 0.90, ]
impure_f <- top_per_cluster[top_per_cluster$top_fine_share  < 0.90, ]
cat(sprintf("\nClusters with < 90%% broad purity: %d\n", nrow(impure_b)))
if (nrow(impure_b)) print(as.data.frame(impure_b), row.names = FALSE)
cat(sprintf("Clusters with < 90%% fine  purity: %d\n", nrow(impure_f)))
if (nrow(impure_f)) print(as.data.frame(impure_f), row.names = FALSE)

# --- 2. Marker-panel sanity check -------------------------------------------

# Use the same canonical panels that reproduce_seurat_object.R annotates with,
# so annotation and validation can never drift apart.
source("R/marker_annotation.R")
panels <- broad_marker_panels

DefaultAssay(so) <- "RNA"
if (inherits(so[["RNA"]], "Assay5")) {
  layer <- if ("counts" %in% Layers(so[["RNA"]])) "counts" else Layers(so[["RNA"]])[1]
} else {
  layer <- "counts"
}
cnts <- GetAssayData(so, assay = "RNA", layer = layer)

panel_detect <- lapply(panels, function(m) {
  m <- intersect(m, rownames(cnts))
  if (!length(m)) return(rep(NA_real_, ncol(cnts)))
  Matrix::colMeans(cnts[m, , drop = FALSE] > 0)
})
panel_detect <- as.data.frame(panel_detect, check.names = FALSE)
rownames(panel_detect) <- colnames(cnts)

by_cluster <- panel_detect %>%
  tibble::rownames_to_column("bc") %>%
  mutate(cluster = meta[bc, "cluster"]) %>%
  group_by(cluster) %>%
  summarise(across(all_of(names(panels)), mean), .groups = "drop")

top_panel_per_cluster <- by_cluster %>%
  tidyr::pivot_longer(-cluster, names_to = "panel", values_to = "detect") %>%
  group_by(cluster) %>%
  slice_max(detect, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::rename(top_panel = panel, top_panel_detect = detect)

check <- top_per_cluster %>%
  left_join(top_panel_per_cluster, by = "cluster") %>%
  mutate(panel_matches_broad = top_panel == top_broad)

write.csv(check,     file.path(OUT_DIR, "cluster_marker_check.csv"),   row.names = FALSE)
write.csv(by_cluster, file.path(OUT_DIR, "cluster_panel_detection.csv"), row.names = FALSE)

cat("\n=== Marker-panel sanity: top-scoring panel per cluster vs. assigned broad label ===\n")
print(as.data.frame(check[, c("cluster", "n_cells", "top_broad", "top_panel",
                              "top_panel_detect", "panel_matches_broad")]),
      row.names = FALSE)
disagree <- check[!check$panel_matches_broad, ]
cat(sprintf("\nClusters where the top marker panel disagrees with the transferred broad label: %d\n",
            nrow(disagree)))
if (nrow(disagree)) print(as.data.frame(disagree), row.names = FALSE)

# --- 3. Summary plot for interpretation -------------------------------------
# One figure combining both checks: per-cluster broad/fine purity (bars,
# against the 90% threshold used above) and whether the marker panel agrees
# with the transferred broad label (point shape), ordered by broad purity so
# problem clusters sort to one end.

summary_df <- check %>%
  dplyr::select(cluster, n_cells, top_broad, top_broad_share, top_fine_share) %>%
  tidyr::pivot_longer(c(top_broad_share, top_fine_share),
                       names_to = "purity_type", values_to = "purity") %>%
  dplyr::mutate(purity_type = ifelse(purity_type == "top_broad_share", "broad", "fine"))

marker_flag <- dplyr::distinct(check, cluster, panel_matches_broad)

cluster_order <- top_per_cluster$cluster[order(top_per_cluster$top_broad_share)]
summary_df$cluster  <- factor(summary_df$cluster,  levels = cluster_order)
marker_flag$cluster  <- factor(marker_flag$cluster, levels = cluster_order)

p_summary <- ggplot(summary_df, aes(x = cluster, y = purity, fill = purity_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "red") +
  geom_point(data = marker_flag,
             aes(x = cluster, y = 1.05, shape = panel_matches_broad),
             inherit.aes = FALSE, size = 2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                      labels = c(`TRUE` = "marker matches", `FALSE` = "marker disagrees"),
                      name = "Marker-panel check") +
  coord_flip() +
  labs(x = "Cluster", y = "Label purity (share of top label)",
       fill = "Purity type",
       title = "Validation summary: per-cluster label purity & marker-panel agreement",
       subtitle = "Dashed line = 90% purity threshold") +
  theme_minimal()

ggsave(file.path(OUT_DIR, "validation_summary.pdf"), p_summary,
       width = 8, height = max(4, 0.3 * nrow(top_per_cluster) + 2))
message("Summary plot written to ", file.path(OUT_DIR, "validation_summary.pdf"))

# --- 4. UMAP plots ----------------------------------------------------------

if (!"umap" %in% names(so@reductions)) {
  stop("No UMAP embedding present on the object; expected 'umap'.", call. = FALSE)
}

so$cluster <- meta$cluster
p_cluster <- DimPlot(so, group.by = "cluster",         label = TRUE, repel = TRUE) +
  ggtitle("New object: clusters") + NoLegend()
p_broad   <- DimPlot(so, group.by = "cell_type",       label = TRUE, repel = TRUE) +
  ggtitle("New object: cell_type (transferred from original)")
p_fine    <- DimPlot(so, group.by = "cell_type.fine",  label = TRUE, repel = TRUE) +
  ggtitle("New object: cell_type.fine (transferred from original)")

ggsave(file.path(OUT_DIR, "umap_clusters.pdf"),        p_cluster, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "umap_cell_type.pdf"),       p_broad,   width = 8, height = 6)
ggsave(file.path(OUT_DIR, "umap_cell_type_fine.pdf"),  p_fine,    width = 8, height = 6)

message("\nOutputs written to ", OUT_DIR)
