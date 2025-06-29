# scRNASeq Analysis: FADDosome inhibits lymphoproliferative disease in CASP8neg RIPK3neg mice

This repository contains the code needed to reproduce the heatmap from single-cell RNA sequencing data of the paper FADDosome inhibition effects on lymphoproliferative disease in CASP8-negative, RIPK3-negative mice.

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

By default the dockerfile runs the following commands:
`RUN renv::restore(prompt=FALSE)`
Which will install all libraries used in this project. This can take up to 30-45 minutes. 

Once the container is build you can run:
```bash
# Start Docker environment
docker-compose up -d

Then you can access RStudio of the container at http://localhost:50362 on you host 
[REPLACE 50362 with your randomly generated port number which you can find
in compose.yml]

Username: rstudio
Password: 1rstudio


### Installation Time

- **Docker and R packages**: Installing the container [40-50 minutes first time]
- **Generating Heatmaps**: about 1 minute.


## 3. Demo

### Quick Demo
A subset of the seurat object is provided to test generate heatmap workflow

#### Run Demo
Open the demo_of_generate_heatmaps.Rmd file within rstudio-server and click
the Run button and select Run all. It should usually take few seconds until it
run through. This could also be run outside the container if the corresponding
packages are installed.


#### Expected Output
- `results/demo_seurat_objects.combined.cleansed.annotated.250428.qs`: PDF
containing all the heatmaps you should get.

## 4. Expected Runtime
- **Demo**: Few seconds on standard desktop
- **Full workflow**: 51 minutes

#### Verify Success
```r
# Check outputs exist
file.exists("results/demo_seurat_objects.qs")

# Load results
demo_obj <- qs::qread("results/demo_seurat_objects.combined.cleansed.annotated.250428.qs")
print(demo_obj)
```

## 4. Instructions for Use

### Input Data Format

Organize the seurat object in a folder of oyur preference and when running
SETUP.bash do as described in the installation steps.
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

- **Preprocessed data**: 10.6084/m9.figshare.29425877 
seurat_objects.combined.cleansed.annotated.250428
- **Demo data**: `data/demo/`

## Contact

For issues or questions:
- **Technical problems**: Create GitHub issue
- **Scientific questions**: Contact corresponding author
- Include session info (`sessionInfo()`) when reporting bugs
