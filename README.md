# scRNASeq Analysis: FADDosome inhibits lymphoproliferative disease in CASP8neg RIPK3neg mice

This repository contains the complete analysis pipeline for single-cell RNA sequencing data of the paper FADDosome inhibition effects on lymphoproliferative disease in CASP8-negative, RIPK3-negative mice, which is needed to reproduce the heatmaps in the corresponding paper.

## 1. System Requirements

### Operating Systems 
- **Linux**: Ubuntu 22.04.5 LTS (tested)
Ideally you should be able to run the pipeline on any operating system which supports the docker version we used to generate the docker container and any other compatible versions. However this was tested on the aforementioned Linux System.

### Software Dependencies
- **R**: 4.4.0+ (tested o4.4.1)
- **Rstudio-Server**: Tested with 2025.05.0+496 (Mariposa Orchid) for Ubuntu Focal
- **Docker**: 28.1.1+

### R Package Dependencies (managed by renv)
See renv.lock file.

### Hardware Requirements
To be able to run it through from scratch without interruption you need the following hardware:
- **RAM**: 128 GB RAM
- **Storage**: 50 GB free space (100 GB recommended)
- **CPU**: 4+ cores recommended
- **Network**: Internet connection for installation

To be able to run it starting from the Seurat object containing the expression data and the annotation of cell types you need
- **RAM**: 32 GB RAM
- **Storage**: About 10 GB free space (50 Recommended)
- **CPU**: 

## 2. Installation Guide

### Prerequisites
Install Docker and Git on your system. 

### Installation Steps

```bash
# Clone repository
git clone https://github.com/CECAD-Bioinformatics-Facility-Projects/scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice.git
cd scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice
```
Then you need to run ./SETUP.bash and follow the instructions. This will adjust
the docker environment to be compatible with yours. When asked
`How many additional directories would you like to mount?` type 1 and provide
as source directory the directory where you downloaded the data. As destination
provide /home/rstudio/project/data/

When finishing the setup a random port number will be generated. This port will
be stored in compose.yml

Now you can run:
```bash
# Build Docker image
docker-compose build
```

This will mostly take quite a few time (see below)

```bash
# Start Docker environment
docker-compose up -d

# Access RStudio at http://localhost:50362 [REPLACE 50362 with your randomly generated port number]
# Username: rstudio, Password: 1rstudio
```

By default the dockerfile runs the following commands:
RUN renv::restore(prompt=FALSE)

### Installation Time

When running from scratch [counts matrices]:
- **Docker and R packages**: Installing the container [40-50 minutes first time]

When running the seurat object 

## 3. Demo

### Quick Demo
A simulated dataset is provided in `data/demo/` for testing.

#### Run Demo
```r
# Load libraries
library(targets)
library(Seurat)

# Execute demo pipeline
tar_make()
```

#### Expected Output
- `results/demo_seurat_objects.qs`: Processed data
- `results/plots/qc_metrics.pdf`: Quality control plots
- `results/plots/umap_clusters.pdf`: UMAP visualization
- `results/de_results.csv`: Differential expression results

#### Expected Runtime
- **Demo**: 5-10 minutes on standard desktop
- **Full analysis**: 2-4 hours

#### Verify Success
```r
# Check outputs exist
file.exists("results/demo_seurat_objects.qs")

# Load results
demo_obj <- qs::qread("results/demo_seurat_objects.qs")
print(demo_obj)
```

## 4. Instructions for Use

### Input Data Format
Organize your 10X Genomics data as:
```
data/your_experiment/
├── sample1/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── sample2/
    └── ...
```

### Configuration
Edit `_targets.R` to specify your samples:
```r
# Update sample paths
sample_info <- list(
  "sample1" = "data/your_experiment/sample1/",
  "sample2" = "data/your_experiment/sample2/"
)

# Adjust parameters
min_genes_per_cell <- 200
max_genes_per_cell <- 5000
max_mitochondrial_percent <- 20
clustering_resolution <- 0.5
```

### Run Analysis
Open rstudio-server
```r
# Complete pipeline
tar_make()
```

### Outputs
- **Seurat objects and other R objects**: `results/*.qs`
- **Plots**: `results/plots/*.pdf`

## Reproduction Instructions

### Reproduce Manuscript Figures

## Data Availability

- **Preprocessed data**: https://figshare.com/s/b0aa00e0cb8e70a68429
- **Demo data**: `data/demo/`

## Contact

For issues or questions:
- **Technical problems**: Create GitHub issue
- **Scientific questions**: Contact corresponding author
- Include session info (`sessionInfo()`) when reporting bugs
