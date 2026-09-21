# Reproducing the integrated Seurat object from raw data

In the following we describe for readers who would like to rebuild the
integrated, annotated Seurat object **from the raw per-sample count
matrices**, rather than starting from the published annotated object (see
the main [README](README.md)).

If you only need to reproduce the paper's figures, you don't need this
document. Instead, just follow the main README, starting from the 
published annotated Seurat object on Figshare.

---


## 0. System Requirements
See [Data Availability](README.md#system-requirements) section for most of system
requirements. However as this analysis is too resource consuming you need to
increased hardware capabilities. 

### Hardware
- **RAM**: 128 GB recommended (64 GB minimum)
- **Storage**: ~50 GB free
  - ~5 GB for the input Seurat object
  - ~10 GB for the loaded Docker image
  - remainder for intermediate results and outputs
- **CPU**: 4+ cores recommended
- **Network**: required only for the one-time download of the Docker image
  tarball and the input data


## 1. Prerequisites

### 1.1 Raw count matrices

Nine count matrices are required, one per sample in the known mtx format. 
After downloading them, make sure, they are placed under `data/` in the cloned repository, 
using exactly this layout (each folder holds the usual `matrix.mtx.gz`,
`barcodes.tsv.gz`, `features.tsv.gz` triplet or the uncompressed versions):

```
data/
├── Casp8_CS_REP1/
├── Casp8_CS_REP2/
├── Casp8_CS_REP3/
├── Casp8_KO_REP1/
├── Casp8_KO_REP2/
├── Casp8_WT_REP1/
├── Casp8_WT_REP2/
├── Casp8_WT_REP3/
└── Casp8_WT_REP4/
```
This means also that you need to rename the downloaded files once you put them
in the according folders as they have prefixes on GEO.

> **Raw data availability:** The raw fastq files and the raw counts are available
> via Gene Expression Omnibus repository (GEO) via accession number GSE347887

P.S. To reproduce the raw counts from the fastq files refer to the versions of
cellranger and celescope mentioned in the method description. In the following
we assume you are starting from the raw counts.

### 1.2 Cached doublet calls (optional but recommended)

You have two options regarding doublet calling and removal. You can use
the published run's doublet calls, which are provided on Figshare 
(see the main README's
[Data Availability](README.md#data-availability) section) so the QC step
can be exactly reproduced. This is the default `DOUBLET_MODE=cached`. 
If you prefer to call doublets from scratch instead (e.g. you don't
have this file, or want to test the `scDblFinder` step itself), set
`DOUBLET_MODE=compute` when running `reproduce_seurat_object.R` (see
below). In this latter case, no cached file is needed. Please expect 
increased stochastic effect on called doublets in that case. If you go for the
cached version download the file and place it under data. 

### 1.3 Environment

Follow the main README's [Installation Guide](README.md#2-installation-guide)
(Steps 1–5) first, so you have a running container with `renv::restore()`
already applied.

---

## 2. Steps

Run the following from the project root, either
1. in the RStudio Server console or 
2. from the terminal via `docker compose exec -T rstudio Rscript ...`.

### Step 1 — Build the nine per-sample objects (`targets`)

```r
targets::tar_make()
```

This runs QC, normalization, variable-feature selection, PCA, neighbor
graph construction and clustering independently for each of the nine
samples (81 targets total). This depends on your system's resources. 
On average expect roughly 15–20 minutes.

### Step 2 — Export the per-sample objects

Either by the running the script using RStudio's GUI or by running

```bash
Rscript export_per_sample_objects.R
```
in the terminal.

This steps reads targets objects of all nine samples 
(`CSR1, CSR2, CSR3, KOR1, KOR2, WTR1, WTR2, WTR3, WTR4`) and export them to 
`results/individual_seurat_objs/<sample>.qs` 

### Step 3 — Integrate, cluster, filter and annotate

```bash
Rscript reproduce_seurat_object.R
```

This step take most of the running time (in our case it took roughly 1-1.5 
hours on 4 cores, mainly due RPCA anchor-finding across all sample pairs). 

It performs, in order:
1. RPCA integration of the nine samples
2. Multi-resolution graph clustering
3. Doublet calling (cached or computed - see section 1.2)
4. QC filtering (singlets, `%mt <= 5`, `200 < nFeature < 5000`) — no
   cluster is removed at this stage
5. Re-clustering (resolution 0.1) + UMAP
6. Sub-clustering of the two mixed clusters
7. Marker-based cell-type annotation

Every stage is cached under `results/reproduced_<timestamp>/`, so a
re-run resumes where it left off. If you one to recompute a step, you need to
delete the corresponding stage file from the cache. Configuration is via 
environment variables-see the header comment of
`reproduce_seurat_object.R` for the full list (`PER_SAMPLE_DIR`,
`OUT_DIR`, `DOUBLET_MODE`, `CACHED_DOUBLET_QS`, `DOUBLET_SEED`,
`ANNOTATION_MODE`, `REFERENCE_QS`). However in each case there is a default 
value. So in case your settings match the defaults it will run without 
setting the environment variables.

The final object is written to:

```
results/reproduced_<timestamp>/seurat_objects.combined.cleansed.annotated.qs
```

### Step 4 — Validate the annotation

If you like you can run a annotation validation script using

```bash
NEW_QS=results/reproduced_<timestamp>/seurat_objects.combined.cleansed.annotated.qs \
  Rscript validate_annotation.R
```

Or using RStudio's GUI.

This script checks cluster label purity and cross-checks each cluster's assigned
broad cell type against canonical marker-panel scores. A successful run
reports 0 clusters below 90% purity and ideally 0 marker-panel disagreements, 
and writes detailed tables and plots to `results/annotation_validation/`.

---

## 3. Expected result

The pipeline includes stochastic steps, particularly doublet detection with 
scDblFinder, exact agreement between runs is not expected. In our tests, 
rebuilding the Seurat object from raw data and rerunning scDblFinder with 
different seeds yielded the consistent population annotations, and highly similar 
downstream results. Using the cached doublet calls yielded almost no differences
across runs. Using the published object yielded exactly the same results.

## 4. Troubleshooting

- **`Error: DOUBLET_MODE=cached but CACHED_DOUBLET_QS not found`** — you
  are using the default cached mode but `data/doublets_original.qs` is
  missing. Either add the file (§1.2) or re-run with
  `DOUBLET_MODE=compute Rscript reproduce_seurat_object.R`.
- **`renv::status()` reports packages "not installed"`** — run
  `renv::restore(prompt = FALSE)` to bring the library back in sync with
  `renv.lock` (this happens if the container was recreated after a
  package was installed manually instead of via `renv.lock` + rebuild).
