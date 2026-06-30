# scRNASeq Analysis: FADDosome inhibits lymphoproliferative disease in CASP8neg RIPK3neg mice

This repository contains the code and environment needed to reproduce
figures from the paper on FADDosome inhibition effects on lymphoproliferative
disease in CASP8-negative, RIPK3-negative mice.

Two independent, self-contained analyses are provided:

| Analysis | Script | Output |
|----------|--------|--------|
| **Pseudobulk heatmaps** | `generate_heatmaps.Rmd` (demo: `demo/demo_of_generate_heatmaps.Rmd`) | Heatmap figures |
| **Erythroid compositional analysis (propeller)** | `propeller_composition_report.Rmd` | Propeller bar plot (PDF/PNG/TIFF/SVG) + source-data Excel |

Both use the same Docker environment, the same `renv.lock`, and the same
input data (the annotated Seurat object on Figshare).

---

## 1. System Requirements

### Operating systems
- **Linux**: Ubuntu 22.04.5 LTS (tested)
- Should work on any OS that supports a recent Docker Engine (macOS, Windows
  with WSL2), but only the Linux configuration above has been tested.
- It is also possible to run the analyses outside the container if every
  R package in `renv.lock` is installed on the host — not tested.

### Software dependencies
Tested with:
- **Docker Engine**: 28.1.1+ (with `docker compose` v2 plugin)
- **R**: 4.4.1 (inside the container)
- **RStudio Server**: rocker/rstudio:4.4.1 base image
- **Git**: any recent version (for cloning)

### R package dependencies
All R package versions are pinned in `renv.lock` and shipped inside the
prebuilt Docker image (Option A) or installed automatically when the image
is built (Option B).

### Hardware
- **RAM**: 16 GB recommended (8 GB minimum for the demo only; the full
  Seurat object is large)
- **Storage**: ~50 GB free
  - ~5 GB for the input Seurat object
  - ~10 GB for the loaded Docker image
  - remainder for intermediate results and outputs
- **CPU**: 4+ cores recommended
- **Network**: required only for the one-time download of the Docker image
  tarball and the input data

---

## 2. Installation Guide

### Prerequisites
Install Docker (with the `docker compose` v2 plugin) and Git on your system.

### Step 1 — Clone the repository

```bash
git clone https://github.com/CECAD-Bioinformatics-Facility-Projects/scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice.git
cd scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice
```

### Step 2 — Download the input data

Download the annotated Seurat object(s) from Figshare (see
[Data Availability](#data-availability) for direct links) into a directory
on your host machine — for example `~/faddosome_data/`. You will mount this
directory into the container in Step 3.

### Step 3 — Configure mounts and port with `SETUP.bash`

`SETUP.bash` interactively generates `compose.yml` with your project name,
data mount, host port, and password.

```bash
./SETUP.bash
```

Answer the prompts as follows:

| Prompt | Answer |
|--------|--------|
| `Project Name` | any name, e.g. `faddosome-repro` |
| `Mount the working directory…?` | `y` |
| `Mount the Default R .cache…?` | **`n`** (mounting the host `.cache` shadows the prebuilt R library and will break Option A) |
| `How many additional directories would you like to mount?` | `1` |
| `Source` | absolute path to your data directory, e.g. `/home/youruser/faddosome_data` |
| `Destination` | `/home/rstudio/project/data` |
| `password` | `y` to keep the default `1rstudio`, or `c` to set your own |
| `Do you accept…?` | `y` |
| `Would you like me to start the container?` | `n` (you will start it explicitly in Step 4) |

`SETUP.bash` writes `compose.yml` with `image: faddosome-casp8-ripk3:4.4.1`
and a randomly chosen free host port. Note the port — you will need it to
open RStudio in your browser.

### Step 4 — Obtain the analysis environment

There are two ways to obtain the environment.
**Option A (recommended)** uses the prebuilt Docker image so you reproduce
the *exact* environment used for the paper without compiling any packages
and without depending on CRAN/Bioconductor staying online.
**Option B** rebuilds the image from `Dockerfile` + `renv.lock`.

---

#### Option A — Use the prebuilt image (recommended)

The prebuilt image `faddosome-casp8-ripk3:4.4.1` contains R 4.4.1 and all
required packages already installed. Download the image tarball from
Figshare, load it, then start the container — no build step and no package
installation are required.

> Download the image tarball from
> <https://doi.org/10.6084/m9.figshare.32835902>
> (private review link available on request — see [Data Availability](#data-availability)).

```bash
# Download faddosome-casp8-ripk3-4.4.1.tar.gz from Figshare into the
# project root, then load it into your local Docker:
docker load < faddosome-casp8-ripk3-4.4.1.tar.gz

# Verify the image is present:
docker images | grep faddosome-casp8-ripk3
# Expected:  faddosome-casp8-ripk3   4.4.1   <id>   <date>   ~9 GB

# Start the container using the loaded image (DO NOT run `docker compose build`):
docker compose up -d
```

Because `compose.yml` defines both an `image:` and a `build:` section, use
`docker compose up` (not `docker compose build`) so the prebuilt image is
used as-is.

---

#### Option B — Build the image from source (alternative)

Use this path if you cannot obtain the prebuilt image, need a different
CPU architecture, or want to regenerate the environment from the recipe.

```bash
DOCKER_BUILDKIT=1 docker compose build
docker compose up -d
```

The build installs all R packages listed in `renv.lock` into an external
library at `/home/rstudio/renv-library` inside the image. Expect 40–60
minutes the first time (most of it compiling Bioconductor packages).

---

### Step 5 — Open RStudio Server in your browser

```
http://localhost:<PORT>
```

Replace `<PORT>` with the port `SETUP.bash` wrote into `compose.yml` (look
for the `ports:` entry, e.g. `"50362:8787"` → use `50362`).

| | |
|---|---|
| **Username** | `rstudio` |
| **Password** | `1rstudio` (or whatever you set in `SETUP.bash`) |

You should land in the `/home/rstudio/project` directory, with the project
files in the Files pane.

### Step 6 — Verify the environment is working

In RStudio's R console:

```r
source("renv/activate.R")
.libPaths()[1]
# Expected: "/home/rstudio/renv-library/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu"

library(Seurat); library(rmarkdown); library(openxlsx); library(qs)
# All four should load without errors.
```

### Stopping the container

```bash
docker compose down
```

This stops and removes the container but leaves the image and your data
intact. To restart later, just run `docker compose up -d` again.

### Installation time

- **Option A (prebuilt image)**: a few minutes (download + `docker load`).
- **Option B (build from source)**: 40–60 minutes the first time.

---

## 3. Demo

The demo satisfies the Nature Research software policy requirement for a
quick end-to-end test that verifies the environment and code are working
correctly. It runs the **exact same analysis pipeline** as the full
workflow but on a small bundled subset of the data, so it completes in
**seconds rather than ~1.5 minutes**.

> If you are a reviewer or editor checking software availability: running
> the demo is all that is needed to verify reproducibility. The full
> dataset is available on Figshare (see [Data Availability](#data-availability))
> for complete reproduction of the paper figures.

### Run the demo

1. In RStudio Server, open `demo/demo_of_generate_heatmaps.Rmd`.
2. Click **Knit** (or **Run → Run All**, <kbd>Ctrl</kbd>+<kbd>Alt</kbd>+<kbd>R</kbd>).

### Expected output

| File | Description |
|------|-------------|
| `results/demo_hms_pp_selected.pdf` | PDF containing the demo heatmaps |

### Verify success

```r
file.exists("results/demo_hms_pp_selected.pdf")
# TRUE
```

### Expected runtime
**A few seconds** on a standard desktop.

---

## 4. Reproducing the full analyses

### 4.1 Pseudobulk heatmaps (full dataset)

1. Make sure the full Seurat object (`seurat_objects.combined.cleansed.annotated.250428.qs`,
   see [Data Availability](#data-availability)) is present in the mounted
   data directory.
2. In RStudio Server, open `generate_heatmaps.Rmd`.
3. Click **Run → Run All** or **Knit**.

**Expected output**: `results/hms_pp_selected.pdf` — publication
heatmaps matching the manuscript figures.

**Expected runtime**: ~1.5 minutes on a standard desktop (previously
~51 minutes before the prebuilt image was available, as most of that
time was `renv::restore` compiling packages during the build).

### 4.2 Erythroid compositional analysis (propeller plot)

This analysis tests whether the proportion of erythroid cells at each
maturation stage (progenitor, erythroblast, late erythroblast) differs
across the three CASP8 genotypes (WT, KO, CS), using the propeller
compositional test (Phipson et al., *Bioinformatics* 2022).

> **Note on demo:** A separate demo version of this analysis is not
> provided. The propeller test requires sufficient cells per
> sample/condition to produce statistically meaningful results — a
> reduced subset would give degenerate output. The full analysis
> completes in 2–3 minutes, which is short enough to serve as its own
> demonstration.

**Option A — RStudio (interactive):**
Open `propeller_composition_report.Rmd` and click **Knit**.

**Option B — command line:**
```bash
docker compose exec -u rstudio rstudio \
  Rscript -e 'rmarkdown::render("/home/rstudio/project/propeller_composition_report.Rmd")'
```

**Expected output:**

| File | Description |
|------|-------------|
| `propeller_composition_report.html` | Full report with interactive tables (Excel download buttons) |
| `figures/erythroid_propeller_composition.pdf` | Vector figure (publication) |
| `figures/erythroid_propeller_composition.tiff` | 300 dpi raster (journal submission) |
| `figures/erythroid_propeller_composition.png` | 300 dpi raster (presentations) |
| `figures/erythroid_propeller_composition.svg` | Editable vector |
| `figures/erythroid_propeller_qq.{pdf,png,svg,tiff}` | Diagnostic QQ plots |
| `figures/erythroid_propeller_tables.xlsx` | Source data + propeller statistics (two sheets) |

**Expected runtime**: 2–3 minutes on a standard desktop (dominated by
loading the Seurat object).

---

## 5. Reproducibility

Both analyses set a fixed random seed (`set.seed(42)`) so jitter, label
placement, and any stochastic steps are deterministic across runs. The
exact R version (4.4.1) and exact package versions are pinned by the
Docker image and `renv.lock`.

### Input data format

The mounted data directory must contain the Seurat object(s) listed below
under [Data Availability](#data-availability), as plain `.qs` files at the
top level (no subdirectory).

---

## Data Availability

The annotated Seurat objects are deposited on Figshare:
**DOI: <https://doi.org/10.6084/m9.figshare.29425877>**

> **Note for editors and reviewers:** The Figshare deposit is currently
> under embargo pending publication. A private review link will be shared
> directly with the editor and reviewers upon request. Once the paper is
> accepted and published, all files will become publicly accessible via
> the DOI above and the direct download links below.

- **Demo data** — `demo_seurat_objects.combined.cleansed.annotated.250428.qs`
  - Public link (live now): <https://figshare.com/ndownloader/files/55758923>
- **Full data** — `seurat_objects.combined.cleansed.annotated.250428.qs`
  - Public link: *(to be activated on publication — navigate the DOI above)*
  - Private review link: *(to be shared directly with editors/reviewers on request)*
- **Prebuilt Docker image tarball** — `faddosome-casp8-ripk3-4.4.1.tar.gz`
  (~3.7 GB) — DOI: <https://doi.org/10.6084/m9.figshare.32835902>
  - Public link: *(to be activated on publication)*
  - Private review link: *(to be shared directly with editors/reviewers on request)*

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `there is no package called 'rmarkdown'` (or any other package) | Host `.cache` mounted into the container shadowed the prebuilt library, or you ran `docker compose build` against an unpatched fork. | Edit `compose.yml` and remove any bind mount whose target is `/home/rstudio/.cache`, then `docker compose down && docker compose up -d`. With Option A the image already ships the library. |
| RStudio login fails | Wrong password | Default is `1rstudio` unless changed in `SETUP.bash`. Check `PASSWORD=` in `compose.yml`. |
| `Cannot find Seurat object at /home/rstudio/project/data/…` | Data directory not mounted, or filename mismatch | Re-run `SETUP.bash` and double-check the data source path and the exact filename. |
| Browser shows port not reachable | Different port than expected | Look at the `ports:` line in `compose.yml`; the host port may be different from 50362. |

---

## Contact

For issues or questions:
- **Technical problems**: open a GitHub issue (include `sessionInfo()` and
  the `docker compose logs` output)
- **Scientific questions**: contact the corresponding author

## License

This project is released under the **MIT License** (see `LICENSE`). The
custom analysis code is available without restriction. If you use it,
please cite via the `CITATION.cff` file.
