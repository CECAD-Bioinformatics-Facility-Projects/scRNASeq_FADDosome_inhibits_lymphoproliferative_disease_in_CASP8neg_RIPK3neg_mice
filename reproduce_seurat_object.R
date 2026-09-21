#!/usr/bin/env Rscript
# =============================================================================
# Reproduce the integrated, annotated Seurat object from the raw data.
# =============================================================================
#
# This script is written for READERS who want to reproduce the single-cell
# object from scratch. It takes the nine per-sample Seurat objects created by the
# targets pipeline and produces the final integrated, clustered and annotated 
# object.
#
# Pipeline (nothing is discarded beyond standard doublet / QC filtering):
#
#   1. RPCA integration of the nine samples
#   2. Multi-resolution graph clustering
#   3. Doublet calling (scDblFinder)      <- optional: run it or 
#											[default] use the called doublets from figshare.
#   4. QC filtering (singlets, %mt, nFeature)      
#   5. Re-clustering (resolution 0.1) + UMAP
#   6. Sub-clustering of the two mixed clusters
#   7. Marker-based cell-type annotation         
#
#   -> results/<out>/seurat_objects.combined.cleansed.annotated.qs
#
# -----------------------------------------------------------------------------
# STEP 0 (prerequisite): build the per-sample objects from the raw matrices.
#
#   The nine per-sample objects are produced by the targets pipeline in this
#   repository from the raw Cell Ranger / Singleron count matrices:
#
#       Rscript -e 'targets::tar_make()'
#
#   That writes results/individual_seurat_objs/<sample>.qs. Point this script
#   at that directory with PER_SAMPLE_DIR if you keep them elsewhere.
# -----------------------------------------------------------------------------
#
# Run from the project root:
#
#   Rscript reproduce_seurat_object.R
#
# Configuration (all optional, via environment variables):
#
#   PER_SAMPLE_DIR  directory with the nine per-sample .qs objects
#                   (default: results/individual_seurat_objs)
#   OUT_DIR         output directory
#                   (default: results/reproduced_<timestamp>)
#
#   DOUBLET_MODE    "cached" (default) or "compute".
#                     cached  = reuse the exact doublet calls from the
#                               published run (CACHED_DOUBLET_QS). This is the
#                               default so the QC step is bit-for-bit
#                               reproducible.
#                     compute = run scDblFinder from scratch (DOUBLET_SEED).
#   CACHED_DOUBLET_QS  barcode-keyed doublet calls used when DOUBLET_MODE=cached
#                   (default: data/doublets_original.qs)
#   DOUBLET_SEED    seed for scDblFinder when DOUBLET_MODE=compute (default 42)
#
#   ANNOTATION_MODE "markers" (default) or "reference".
#                     markers   = assign cell types from canonical marker
#                                 panels (see R/marker_annotation.R). Fully
#                                 self-contained; robust to cluster renumbering.
#                     reference = copy cell_type / cell_type.fine by barcode
#                                 from a provided published object
#                                 (REFERENCE_QS), for an exact label match.
#   REFERENCE_QS    published object used when ANNOTATION_MODE=reference
#                   (default: data/seurat_objects.combined.cleansed.annotated.250428.qs)
#
# The script never overwrites an existing output; every stage is cached so a
# re-run resumes where it stopped. Delete a stage file to recompute it.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
})

options(future.globals.maxSize = 10 * 1024^3)

# ---- configuration ----------------------------------------------------------

# PS. "unset" defines the fallback value, so if this script is run using rstudio
# rather than rscript with specification of environment variables, then the 
# values of the unset parameter are used as fallback. 

PER_SAMPLE_DIR    <- Sys.getenv("PER_SAMPLE_DIR", unset = "results/individual_seurat_objs")
OUT_DIR           <- Sys.getenv("OUT_DIR",
                                unset = file.path("results",
                                                  paste0("reproduced_", format(Sys.time(), "%Y%m%d_%H%M%S"))))
DOUBLET_MODE      <- match.arg(Sys.getenv("DOUBLET_MODE", unset = "cached"),
                               choices = c("cached", "compute"))
CACHED_DOUBLET_QS <- Sys.getenv("CACHED_DOUBLET_QS", unset = "data/doublets_original.qs")
DOUBLET_SEED      <- as.integer(Sys.getenv("DOUBLET_SEED", unset = "42"))
ANNOTATION_MODE   <- match.arg(Sys.getenv("ANNOTATION_MODE", unset = "markers"),
                               choices = c("markers", "reference"))
REFERENCE_QS      <- Sys.getenv("REFERENCE_QS",
                                unset = "data/seurat_objects.combined.cleansed.annotated.250428.qs")

sample_names <- c("CSR1", "CSR2", "WTR1", "WTR2", "WTR3", "WTR4", "CSR3", "KOR1", "KOR2")
per_sample_qs <- file.path(PER_SAMPLE_DIR, paste0(sample_names, ".qs"))

log_step <- function(msg) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))

log_step("Configuration:")
message("  PER_SAMPLE_DIR  : ", PER_SAMPLE_DIR)
message("  OUT_DIR         : ", OUT_DIR)
message("  DOUBLET_MODE    : ", DOUBLET_MODE)
if (DOUBLET_MODE == "cached")  message("  CACHED_DOUBLET_QS: ", CACHED_DOUBLET_QS)
if (DOUBLET_MODE == "compute") message("  DOUBLET_SEED    : ", DOUBLET_SEED)
message("  ANNOTATION_MODE : ", ANNOTATION_MODE)
if (ANNOTATION_MODE == "reference") message("  REFERENCE_QS    : ", REFERENCE_QS)

# ---- input checks -----------------------------------------------------------

missing <- per_sample_qs[!file.exists(per_sample_qs)]
if (length(missing)) {
  stop("Missing per-sample input(s):\n  ", paste(missing, collapse = "\n  "),
       "\n\nBuild them first with:  Rscript -e 'targets::tar_make()'", call. = FALSE)
}
if (DOUBLET_MODE == "cached" && !file.exists(CACHED_DOUBLET_QS)) {
  stop("DOUBLET_MODE=cached but CACHED_DOUBLET_QS not found: ", CACHED_DOUBLET_QS,
       "\nEither provide the CACHED_DOUBLET_QS or set DOUBLET_MODE=compute to call doublets from scratch instead.", call. = FALSE)
}
if (ANNOTATION_MODE == "reference" && !file.exists(REFERENCE_QS)) {
  stop("ANNOTATION_MODE=reference but REFERENCE_QS not found: ", REFERENCE_QS,
       "\nEither provie the REFERENCE_QS or use the default ANNOTATION_MODE=markers to annotate without it.", call. = FALSE)
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

f_integrated <- file.path(OUT_DIR, "01_integrated.qs")
f_clustered  <- file.path(OUT_DIR, "02_clustered.qs")
f_doublet    <- file.path(OUT_DIR, "03_doublets.qs")
f_filtered   <- file.path(OUT_DIR, "04_filtered.qs")
f_reclust    <- file.path(OUT_DIR, "05_reclustered.qs")
f_final      <- file.path(OUT_DIR, "seurat_objects.combined.cleansed.annotated.qs")

qsave_safe <- function(obj, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing file: ", path, call. = FALSE)
  qs::qsave(obj, path)
  message("  saved: ", path)
}

# =============================================================================
# 1. RPCA integration
# =============================================================================
# All seeds below are Seurat's defaults, matching the original analysis:
#   RunPCA / RunUMAP  seed.use    = 42 (Seurat default)
#   FindClusters      random.seed = 0  (Seurat default; not set explicitly)

if (!file.exists(f_integrated)) {
  log_step("[1/7] Integrating the nine samples (RPCA) ...")
  seurat_objects <- lapply(per_sample_qs, qs::qread)
  names(seurat_objects) <- sample_names

  # R/sc_workflow.R creates every per-sample object with a hardcoded
  # project = "SeuratAnalysis", so orig.ident is not sample-specific at this
  # point. So we set it here.
  for (nm in names(seurat_objects)) seurat_objects[[nm]]$orig.ident <- nm

  seurat_objects <- lapply(seurat_objects, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x
  })

  features <- SelectIntegrationFeatures(object.list = seurat_objects)
  seurat_objects <- lapply(seurat_objects, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
  })

  anchors <- FindIntegrationAnchors(object.list    = seurat_objects,
                                    anchor.features = features,
                                    reduction       = "rpca")
  combined <- IntegrateData(anchorset = anchors)
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)

  qsave_safe(combined, f_integrated)
} else {
  log_step("[1/7] Loading cached integrated object")
  combined <- qs::qread(f_integrated)
}

# =============================================================================
# 2. Multi-resolution clustering
# =============================================================================

if (!file.exists(f_clustered)) {
  log_step("[2/7] Multi-resolution clustering ...")
  for (res in c(0.0125, 0.025, 0.05, 0.1, 0.5)) {
    message("  resolution: ", res)
    combined <- FindClusters(combined, resolution = as.numeric(res))
  }
  qsave_safe(combined, f_clustered)
} else {
  log_step("[2/7] Loading cached clustered object")
  combined <- qs::qread(f_clustered)
}

# =============================================================================
# 3. Doublet detection (scDblFinder)
# =============================================================================

if (!file.exists(f_doublet)) {
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
  DefaultAssay(combined) <- "RNA"

  if (DOUBLET_MODE == "cached") {
    log_step(paste0("[3/7] Loading published doublet calls: ", CACHED_DOUBLET_QS))
    cached      <- qs::qread(CACHED_DOUBLET_QS)
    cached_meta <- if (is(cached, "Seurat")) cached@meta.data else as.data.frame(cached)
    dbl_cols    <- grep("^scDblFinder", colnames(cached_meta), value = TRUE)
    if (!length(dbl_cols)) stop("CACHED_DOUBLET_QS has no scDblFinder.* columns.", call. = FALSE)

    our_bc  <- colnames(combined)
    absent  <- setdiff(our_bc, rownames(cached_meta))
    if (length(absent)) {
      stop(sprintf("Cached doublet calls miss %d of our cells (e.g. %s).",
                   length(absent), paste(head(absent, 3), collapse = ", ")), call. = FALSE)
    }
    combined <- AddMetaData(combined, metadata = cached_meta[our_bc, dbl_cols, drop = FALSE])
    rm(cached, cached_meta); gc(verbose = FALSE)
  } else {
    log_step(paste0("[3/7] Running scDblFinder (seed = ", DOUBLET_SEED, ") ..."))
    sce <- as.SingleCellExperiment(combined)
    set.seed(DOUBLET_SEED)
    sce <- scDblFinder(sce, samples = "orig.ident", clusters = "seurat_clusters")
    doublet_info <- as.data.frame(colData(sce))
    doublet_info <- doublet_info[, grepl("^scDblFinder", colnames(doublet_info)), drop = FALSE]
    rownames(doublet_info) <- colnames(combined)
    combined <- AddMetaData(combined, metadata = doublet_info)

    # Persist just the calls so a later run can reuse them via DOUBLET_MODE=cached.
    qs::qsave(doublet_info, file.path(OUT_DIR, "doublet_calls.qs"))
  }

  qsave_safe(combined, f_doublet)
} else {
  log_step("[3/7] Loading cached doublet-annotated object")
  combined <- qs::qread(f_doublet)
}

# =============================================================================
# 4. QC filtering + PCA  (no cluster is removed)
# =============================================================================
# Keep singlets with reasonable library complexity and low mitochondrial
# content. Unlike the original run, NO low-quality cluster is dropped by ID:
# Louvain cluster numbers are run-specific, so an ID-based drop is not
# reproducible. The result is therefore a strict superset of the published
# object (~17 extra low-count cells, 0.017%); downstream results are unaffected.

if (!file.exists(f_filtered)) {
  log_step("[4/7] QC filtering (singlets, %mt <= 5, 200 < nFeature < 5000) ...")
  DefaultAssay(combined) <- "RNA"
  combined@meta.data$percent.mt <- PercentageFeatureSet(combined, pattern = "^mt-")

  combined <- subset(combined, subset = scDblFinder.class == "singlet")
  combined <- subset(combined,
                     subset = percent.mt <= 5 &
                              nFeature_RNA > 200 &
                              nFeature_RNA < 5000)

  DefaultAssay(combined) <- "integrated"
  combined <- RunPCA(combined, npcs = 100)

  qsave_safe(combined, f_filtered)
} else {
  log_step("[4/7] Loading cached filtered object")
  combined <- qs::qread(f_filtered)
}

# =============================================================================
# 5. Re-clustering (resolution 0.1) + UMAP
# =============================================================================

if (!file.exists(f_reclust)) {
  log_step("[5/7] Re-clustering at resolution 0.1 and computing UMAP ...")
  combined <- FindNeighbors(combined, dims = 1:60)
  combined <- FindClusters(combined, resolution = 0.1)
  combined <- RunUMAP(combined, dims = 1:60)
  qsave_safe(combined, f_reclust)
} else {
  log_step("[5/7] Loading cached reclustered object")
  combined <- qs::qread(f_reclust)
}

# =============================================================================
# 6. Sub-clustering of the two mixed clusters
# =============================================================================
# Clusters 4 and 7 each hold more than one lineage at res 0.1; splitting them
# yields the granularity used for annotation. The resulting column name keeps
# the original report's (historical) label.

log_step("[6/7] Sub-clustering clusters 4 and 7 ...")
Idents(combined) <- "integrated_snn_res.0.1"
DefaultAssay(combined) <- "RNA"

combined <- FindSubCluster(
  combined, cluster = 4, resolution = 0.05,
  graph.name = "integrated_snn",
  subcluster.name = "integrated_snn_res.0.1_and_cluster_4"
)
Idents(combined) <- "integrated_snn_res.0.1_and_cluster_4"
combined <- FindSubCluster(
  combined, cluster = 7, resolution = 0.05,
  graph.name = "integrated_snn",
  subcluster.name = "integrated_snn_res.0.1_and_cluster_4_and_cluster_9"
)

cluster_col <- "integrated_snn_res.0.1_and_cluster_4_and_cluster_9"

# =============================================================================
# 7. Cell-type annotation
# =============================================================================

if (ANNOTATION_MODE == "markers") {
  log_step("[7/7] Annotating cell types from canonical markers ...")
  source("R/marker_annotation.R")
  ann <- annotate_by_markers(combined, cluster_col = cluster_col)
  combined <- ann$object

  message("\n  cluster -> broad / fine label (marker-derived):")
  map_tbl <- data.frame(cluster = names(ann$broad_map),
                        cell_type = unname(ann$broad_map),
                        cell_type.fine = unname(ann$fine_map[names(ann$broad_map)]))
  print(map_tbl, row.names = FALSE)
} else {
  log_step("[7/7] Annotating by barcode transfer from the published object ...")
  source("R/annotation_transfer.R")
  orig <- qs::qread(REFERENCE_QS)
  xfer <- transfer_annotation_from_original(combined, orig, cluster_col = cluster_col)
  combined <- xfer$object
  rm(orig); gc(verbose = FALSE)
}

# Batch / technology metadata (used by the downstream reports).
combined@meta.data <- combined@meta.data %>%
  mutate(
    batch = case_when(
      orig.ident %in% c("WTR1", "WTR2", "CSR1")         ~ "batch1",
      orig.ident %in% c("WTR3", "CSR2")                 ~ "batch2",
      orig.ident %in% c("WTR4", "CSR3", "KOR1", "KOR2") ~ "batch3"
    ),
    technology = case_when(
      orig.ident %in% c("WTR1", "WTR2", "CSR1") ~ "singleron",
      TRUE                                      ~ "tenX"
    )
  )

qsave_safe(combined, f_final)

log_step(paste0("Done. Final annotated object: ", f_final))
message("\nNext: verify the annotation with\n",
        "  NEW_QS=", f_final, " Rscript validate_annotation.R")
