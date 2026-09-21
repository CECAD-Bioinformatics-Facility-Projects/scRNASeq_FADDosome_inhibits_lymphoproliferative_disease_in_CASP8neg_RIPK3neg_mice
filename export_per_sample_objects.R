#!/usr/bin/env Rscript
# Export per-sample Seurat objects from the targets pipeline (_targets/) to
# standalone .qs files that reproduce_seurat_object.R can read.
#
# The targets pipeline (_targets.R) builds each sample as a target named
# SO_UMAP_<id>, storing it inside the _targets/ cache. This script reads
# each target and exports it as a plain .qs file in results/individual_seurat_objs/.
#
# Run from the project root:
#   Rscript export_per_sample_objects.R
#
# Requires:
#   - _targets/ cache populated (run: Rscript -e 'targets::tar_make()' first)
#   - targets, Seurat, qs packages installed

suppressPackageStartupMessages({
  library(targets)
  library(qs)
})

sample_names <- c("CSR1", "CSR2", "WTR1", "WTR2", "WTR3", "WTR4", "CSR3", "KOR1", "KOR2")
out_dir <- "results/individual_seurat_objs"

if (!dir.exists("_targets")) {
  stop("_targets/ directory not found. Run  Rscript -e 'targets::tar_make()'  first.", call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
message("Exporting per-sample objects to: ", out_dir)

for (name in sample_names) {
  target_name <- paste0("SO_UMAP_", name)
  out_file <- file.path(out_dir, paste0(name, ".qs"))
  
  if (file.exists(out_file)) {
    message("  ", name, " ... already exists, skipping")
    next
  }
  
  message("  ", name, " ... ", appendLF = FALSE)
  tryCatch({
    obj <- tar_read_raw(target_name)
    qsave(obj, out_file)
    message("saved")
  }, error = function(e) {
    message("ERROR: ", e$message)
  })
}

message("Done. Per-sample objects are ready for reproduce_seurat_object.R")
message("Next: Rscript reproduce_seurat_object.R")
