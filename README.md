# scRNASeq Analysis: FADDosome inhibits lymphoproliferative disease in CASP8neg RIPK3neg mice

This repository contains the code and environment needed to reproduce
figures from the paper on FADDosome inhibition effects on lymphoproliferative
disease in CASP8-negative, RIPK3-negative mice.

Two independent, self-contained analyses are provided:

| Analysis | Script | Output |
|----------|--------|--------|
| **Pseudobulk heatmaps** | `generate_heatmaps.Rmd` (demo: `demo/demo_of_generate_heatmaps.Rmd`) | Heatmap figures |
| **Erythroid compositional analysis (propeller)** | `propeller_composition_report.Rmd` | Propeller bar plot (PDF/PNG/TIFF/SVG) + source-data Excel |

Both use the same Docker environment, renv lockfile, and input data.

## 1. System Requirements

### Operating Systems 
- **Linux**: Ubuntu 22.04.5 LTS (tested)
- You should be able to run the code on any operating system that supports the Docker version we used to generate the Docker container and any other compatible versions. However, this was tested on the aforementioned Linux System.
- Ideally, it should be possible to run the code outside the container if the libraries are installed on the host system as well. But this was not tested.

### Software Dependencies
It was tested with
- **R**: 4.4.0+ 
- **RStudio-Server**: 2025.05.0+496 (Mariposa Orchid) for Ubuntu Focal
- **Docker**: 28.1.1+

### R Package Dependencies (managed by renv)
See renv.lock file.

### Hardware Requirements
To be able to run it you need the following hardware:
- **RAM**: 8 GB RAM
- **Storage**: around 50 GB free space at least 
[For the input, you need around 5GB, and the remaining space is needed for the installation of system libraries and R packages for Docker]
- **CPU**: 4+ cores recommended
- **Network**: Internet connection for installation

## 2. Installation Guide

### Prerequisites
Install Docker and Git on your system. 

### Installation Steps

First, clone the repository (needed for both options below):

```bash
# Clone repository
git clone https://github.com/CECAD-Bioinformatics-Facility-Projects/scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice.git
cd scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice
```

There are two ways to obtain the analysis environment. **Option A (recommended)**
uses the prebuilt Docker image, so you reproduce the *exact* environment used for
the paper without compiling any packages and without depending on package
repositories (CRAN/Bioconductor) staying online. **Option B** rebuilds the image
from source using the `Dockerfile` and `renv.lock`.

Both options use the same `SETUP.bash` helper to configure mounts and a port.
Run it now and follow the prompts. When asked
`How many additional directories would you like to mount?` type `1` and provide
as the source directory the directory where you downloaded the data, and as
destination `/home/rstudio/project/data/`. When finishing the setup a random
port number will be generated and stored in `compose.yml`.

```bash
./SETUP.bash
```

---

#### Option A — Use the prebuilt image (recommended)

The prebuilt image `faddosome-casp8-ripk3:4.4.1` contains R and all required
packages already installed. Download the image tarball from Figshare, load it,
then start the container — no build step and no package installation are
required.

> Replace the download URL below with the Figshare link for the image tarball
> once it is deposited.

**A1. Download and load the image from Figshare:**

```bash
# Download faddosome-casp8-ripk3-4.4.1.tar.gz from Figshare, then load it:
docker load < faddosome-casp8-ripk3-4.4.1.tar.gz
```

Once the image is loaded, start the container **without building**:

```bash
# Start Docker environment (uses the loaded image, skips the build)
docker compose up -d
```

Because `compose.yml` defines both an `image:` and a `build:` section, make sure
to use `docker compose up` (not `docker compose build`) so the prebuilt image is
used as-is.

---

#### Option B — Build the image from source (alternative)

This path rebuilds the environment from the `Dockerfile` and `renv.lock`. Use it
if you cannot obtain the prebuilt image, need a different CPU architecture, or
want to regenerate the environment from the recipe.

```bash
# Build Docker image
docker-compose build
```

By default, the Dockerfile runs the following commands:
`RUN renv::restore(prompt=FALSE)`
This will install all libraries used in this project. This can take up to 30-45 minutes.

```bash
# Start Docker environment
docker-compose up -d
```

---

#### Accessing RStudio (both options)

Then you can access RStudio of the container at http://localhost:50362 on your host 
[REPLACE 50362 with your randomly generated port number, which you can find
in compose.yml]

Username: rstudio
Password: 1rstudio


### Installation Time

- **Option A (prebuilt image)**: a few minutes (download + `docker load`/`pull`); no package installation.
- **Option B (build from source)**: 40-50 minutes the first time (compiles all R packages).
- **Generating Heatmaps**: about 1 minute.


## 3. Demo

### Quick Demo
A subset of the Seurat object is provided to test the generated heatmap workflow

#### Run Demo
Open the demo_of_generate_heatmaps.Rmd file within rstudio-server and click
the Run button, and select Run all. It should usually take a few seconds until it
runs through. This could also be run outside the container if the corresponding
packages are installed.


#### Expected Output
- `results/demo_seurat_objects.combined.cleansed.annotated.250428.qs`: PDF
containing all the heatmaps you should get.

## 4. Expected Runtime
- **Demo**: A Few seconds on a standard desktop
- **Fully reproducible workflow of full dataset**: 51 minutes

#### Verify Success
```r
# Check outputs exist
file.exists("results/demo_seurat_objects.qs")
# Load results
demo_obj <- qs::qread("results/demo_seurat_objects.combined.cleansed.annotated.250428.qs")
print(demo_obj)
```

## 4. Instructions for Use and Reproducibility

### Input Data Format

Organize the Seurat objects in a folder of your preference, and when running
SETUP.bash Do as described in the installation steps.

### Reproducibility

For exact reproducibility, the analyses use a fixed random seed
(`set.seed(42)`), so jitter/label placement and any stochastic steps are
deterministic across runs.

### Run Analysis
After running the container open http://localhost:50362 [according port in compose.yml]
Click in the Files pane in RStudio on the demo folder and open the demo_of_generate_heatmaps.Rmd
file, then run the whole Rmd file by clicking on Run and then Run AL,L or as a shortcut
Ctrl+Alt+R

To reproduce the paper's heatmaps, open the file generate_heatmaps.Rmd and follow
same steps.

### Outputs [Figures]
`hm_pp_selected.pdf` for the full version and `demo_hms_pp_selected.pdf` for the demo version.


## 5. Erythroid Compositional Analysis (Propeller Plot)

This analysis tests whether the proportion of erythroid cells at each
maturation stage (progenitor, erythroblast, late erythroblast) differs
across the three CASP8 genotypes (WT, KO, CS), using the propeller
compositional test (Phipson et al., *Bioinformatics* 2022).

### Run

After the container is running (steps 1–2 above), either:

**Option A — RStudio (interactive):**
Open `propeller_composition_report.Rmd` in RStudio Server and click
**Knit**.

**Option B — command line:**
```bash
docker compose exec -u rstudio rstudio \
  Rscript -e 'rmarkdown::render("/home/rstudio/project/propeller_composition_report.Rmd")'
```

### Expected Output

| File | Description |
|------|-------------|
| `propeller_composition_report.html` | Full report with interactive tables (Excel download buttons) |
| `figures/erythroid_propeller_composition.pdf` | Vector figure (publication) |
| `figures/erythroid_propeller_composition.tiff` | 300 dpi raster (journal submission) |
| `figures/erythroid_propeller_composition.png` | 300 dpi raster (presentations) |
| `figures/erythroid_propeller_tables.xlsx` | Source data + propeller statistics (two sheets) |

### Expected Runtime

About 2–3 minutes on a standard desktop (dominated by loading the Seurat
object).

---

## Data Availability

The annotated Seurat object required for both the heatmap and the propeller
compositional analysis is available on Figshare:
DOI: 10.6084/m9.figshare.29425877

- **Demo data**: Navigate the link and download: demo_seurat_objects.combined.cleansed.annotated.250428.qs 
or click directly on https://figshare.com/ndownloader/files/55758923

- **Full data object**:  Navigate the link and download: seurat_objects.combined.cleansed.annotated.250428.qs 

## Contact

For issues or questions:
- **Technical problems**: Create GitHub issue
- **Scientific questions**: Contact corresponding author
- Include session info (`sessionInfo()`) when reporting bugs

## License

This project is released under the **MIT License** (see `LICENSE`). The custom
analysis code is available without restriction. If you use it, please cite via
the `CITATION.cff` file.
