# Demo Dataset

This directory contains a small simulated single-cell RNA-seq dataset for demonstration and testing purposes.

## Contents

- `demo_counts.rds`: Simulated count matrix (1000 cells x 2000 genes)
- `demo_metadata.csv`: Cell metadata including sample IDs and experimental conditions
- `demo_features.csv`: Gene information and annotations

## Usage

This demo dataset is automatically used when running the demo pipeline:

```r
# Load the demo data
demo_data <- readRDS("data/demo/demo_counts.rds")
demo_metadata <- read.csv("data/demo/demo_metadata.csv")

# Or run the complete demo pipeline
tar_make(demo_analysis)
```

## Dataset Characteristics

- **Cells**: 1,000 simulated cells
- **Genes**: 2,000 genes (including mitochondrial and ribosomal genes)
- **Conditions**: 3 experimental groups (Control, Treatment1, Treatment2)
- **Replicates**: 2 biological replicates per condition
- **Expected clusters**: 5-8 cell clusters
- **Runtime**: ~5-10 minutes for complete analysis

This dataset is designed to test all major pipeline components including quality control, normalization, clustering, and differential expression analysis.
