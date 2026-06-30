# Nature Research Software Policy - Compliance Checklist

## Required content

- [x] **Compiled standalone software and/or source code**
  Docker image (`faddosome-casp8-ripk3:4.4.1`) ships all compiled
  packages; source code in `R/`, `generate_heatmaps.Rmd`,
  `propeller_composition_report.Rmd`, `demo/demo_of_generate_heatmaps.Rmd`.

- [x] **A small (simulated or real) dataset to demo the software/code**
  `demo_seurat_objects.combined.cleansed.annotated.250428.qs` on Figshare
  (DOI: 10.6084/m9.figshare.29425877); used by
  `demo/demo_of_generate_heatmaps.Rmd` and can be use by propeller t-test
  script as well. See README for details.

---

## README - required sections

### 1. System requirements

- [x] All software dependencies listed (Docker Engine 28.1.1+, R 4.4.1,
  RStudio Server rocker:4.4.1, Git) -> README §1
- [x] Operating systems listed with version numbers (Ubuntu 22.04.5 LTS
  tested; macOS/Windows WSL2 noted as compatible) -> README §1
- [x] Versions the software has been tested on -> README §1
- [x] Required non-standard hardware listed (RAM 16 GB, ~50 GB storage,
  4+ CPU cores) -> README §1

### 2. Installation guide

- [x] Step-by-step installation instructions (clone -> download data ->
  SETUP.bash -> docker load / build -> docker compose up) -> README §2
- [x] Typical install time on a "normal" desktop computer
  (Option A: a few minutes; Option B: 40–60 min) -> README §2

### 3. Demo

- [x] Instructions to run on data (`demo/demo_of_generate_heatmaps.Rmd`,
  click Knit) -> README §3
- [x] Expected output (`results/demo_hms_pp_selected.pdf`) -> README §3
- [x] Expected run time for demo on a "normal" desktop (a few seconds)
  -> README §3

### 4. Instructions for use

- [x] How to run the software on your data (full heatmap and propeller
  workflows with Knit / command-line instructions) -> README §4
- [x] *(OPTIONAL)* Reproduction instructions — full quantitative
  reproduction of all manuscript figures covered in README §4 and §5,
  including exact package versions
  via `renv.lock`, and prebuilt Docker image -> README §4–5
