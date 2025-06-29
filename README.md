# scRNASeq Analysis: FADDosome inhibits lymphoproliferative disease in CASP8neg RIPK3neg mice

This repository contains the code and environment needed to reproduce the scRNA-Seq pseudobulk heatmaps from the paper FADDosome inhibition effects on lymphoproliferative disease in CASP8-negative, RIPK3-negative mice.

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

```bash
# Clone repository
git clone https://github.com/CECAD-Bioinformatics-Facility-Projects/scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice.git
cd scRNASeq_FADDosome_inhibits_lymphoproliferative_disease_in_CASP8neg_RIPK3neg_mice
```
Then you need to run ./SETUP.bash and follow the instructions. This will adjust
the Docker environment to be compatible with yours. When asked
`How many additional directories would you like to mount?` Type 1 and provide
as the source directory the directory where you downloaded the data. As destination,
provide /home/rstudio/project/data/

When finishing the setup a random port number will be generated. This port will
be stored in compose.yml

Now you can run:
```bash
# Build Docker image
docker-compose build
```

By default, the Dockerfile runs the following commands:
`RUN renv::restore(prompt=FALSE)`
This will install all libraries used in this project. This can take up to 30-45 minutes. 

Once the container is built, you can run:
```bash
# Start Docker environment
docker-compose up -d
```

Then you can access RStudio of the container at http://localhost:50362 on your host 
[REPLACE 50362 with your randomly generated port number, which you can find
in compose.yml]

Username: rstudio
Password: 1rstudio


### Installation Time

- **Docker and R packages**: Installing the container [40-50 minutes first time]
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

### Run Analysis
After running the container open http://localhost:50362 [according port in compose.yml]
Click in the Files pane in RStudio on the demo folder and open the demo_of_generate_heatmaps.Rmd
file, then run the whole Rmd file by clicking on Run and then Run AL,L or as a shortcut
Ctrl+Alt+R

To reproduce the paper's heatmaps, open the file generate_heatmaps.Rmd and follow
same steps.

### Outputs [Figures]
`hm_pp_selected.pdf` for the full version and `demo_hms_pp_selected.pdf` for the demo version.


## Data Availability

- Seurat objects needed to generate the heatmap for the demo version, as well as
for the full dataset, are available under
- DOI: 10.6084/m9.figshare.29425877 [unpublished yet]
- Private Link For Reviewer: [you got from the journal]

- **Demo data**: Navigate the link and download: demo_seurat_objects.combined.cleansed.annotated.250428.qs 
or click directly on https://figshare.com/ndownloader/files/55758923 [this option doesn't work because item not published yet]

- ** Full data object**:  Navigate the link and download: seurat_objects.combined.cleansed.annotated.250428.qs 

## Contact

For issues or questions:
- **Technical problems**: Create GitHub issue
- **Scientific questions**: Contact corresponding author
- Include session info (`sessionInfo()`) when reporting bugs
