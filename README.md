# Red Swamp Crayfish Gene Flow and Invasion Dynamics

This repository contains code and data processing pipelines for:

[**Adams et al. (2025)** "Genomic analyses reveal multiple introduction events and fine-scale geographic barriers influencing gene flow in a widespread aquatic invader." *Ecology and Evolution* 15(12): e72550.](https://onlinelibrary.wiley.com/doi/10.1002/ece3.72550)

## Overview

This project uses population genomic approaches to investigate the introduction history and genetic structure of invasive red swamp crayfish (*Procambarus clarkii*) populations in southeastern Michigan, USA. Using whole-genome sequencing data from 763 individuals across 20  waterbodies, we use examine genetic diversity and population structure, employ Approximate Bayesian Computation (ABC) to infer demographic histories, and landscape genetics to identify barriers to gene flow.

## Repository Structure

```
.
├── README.md                    # This file
├── code/                        # Data processing and simulation scripts
│   ├── ABC/                     # Approximate Bayesian Computation simulations
│   ├── dataProcessing*.Rmd      # Raw sequence processing pipelines
│   ├── genotyping_and_filtering.Rmd  # SNP calling and filtering
│   └── cleaning_metadata.Rmd    # Metadata standardization
├── analysis/                    # Statistical analyses and visualizations
│   ├── ABC/                     # ABC model selection and parameter estimation
│   ├── genetic_diversity.Rmd    # Population genetic diversity analyses
│   ├── PCA_and_DAPC.Rmd        # Multivariate analyses
│   ├── directionality_index.Rmd # Range expansion direction tests
│   └── RSC_landGen.Rmd         # Landscape genetics analyses
├── data/                        # Input data files
│   ├── ABC/                     # Empirical summary statistics
│   ├── directionality_index/    # Coordinate and genotype files
│   ├── landscape_genetics/      # Environmental raster layers
│   ├── pca_dapc/               # PCA/DAPC results
│   └── metadata_N1278.csv      # Sample metadata
├── output/                      # Analysis results
│   ├── ABC_output/             # Model selection and parameter estimates
│   ├── genetic_diversity_output/  # Diversity metrics
│   ├── directionality_index_output/  # ψ statistics
│   └── landgen_output/         # MLPE results
├── figures/                     # Publication figures
├── tables/                      # Publication tables
└── sensitive/                   # Non-public location data (see Data Availability)
```

## Workflow

### 1. Data Processing

**Sequence Processing → Genotyping → Filtering**

Raw Illumina reads were processed through adapter trimming, quality filtering, mapping to the *P. clarkii* reference genome, and SNP calling. See:

- `code/dataProcessing_fastq2bams_Adams.Rmd` - Adams 2022 sequences
- `code/dataProcessing_fastq2bams_Sard-Homola.Rmd` - Sard & Homola 2021 sequences
- `code/genotyping_and_filtering.Rmd` - SNP calling and quality filters

### 2. Population Structure

**PCA, DAPC, Clustering**

Multivariate analyses see:

- `analysis/PCA_and_DAPC.Rmd` - Principal components and discriminant analysis

### 3. Genetic Diversity

**Heterozygosity, Allelic Richness, F-statistics, etc.**

Population-level diversity metrics, pairwise differentiation, and neighbor-joining tree. See:

- `analysis/genetic_diversity.Rmd` - Comprehensive diversity analyses

### 4. Approximate Bayesian Computation

**Simulations → Model Selection → Parameter Estimation → Cross-Validation**

ABC workflow to infer demographic history using 100,000 simulations per model across 10 competing demographic scenarios. See:

- `code/ABC/RSC_ABCsims_20250722_anon.R` - Forward-time coalescent simulations
- `analysis/ABC/` - ABC pipeline including:
  - Random Forest, multinomial logistic regression, and neural network model selection
  - Cross-validation for model selection
  - Parameter estimation for best-fit models
  - Cross-validation for parameter estimation


### 5. Directionality Analysis

**ψ-statistics → Z-scores**

Tests for directional expansion patterns. See:

- `analysis/directionality_index.Rmd` - Range expansion directionality tests

### 6. Landscape Genetics

**Resistance Surface Optimization → Model Selection**

Optimization of landscape resistance surfaces using genetic algorithms to identify which environmental features best explain patterns of genetic differentiation. See:

- `analysis/RSC_landGen.Rmd` - Landscape genetics analyses using ResistanceGA

## Data Availability

### Public Data

Raw FASTQ files and barcodes for demultiplexing RAD data are available on Dryad [https://doi.org/10.5061/dryad.djh9w0wcq](https://doi.org/10.5061/dryad.djh9w0wcq)

Individual sequences after demultiplexing libraries are available on the NCBI sequence read archive (SRA):

- **BioProject:** [PRJNA1148680](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1148680)

Processed data in this repository:

- Filtered VCF files (`data/`)
- Summary statistics (`output/`)
- Analysis results (`output/`, `tables/`)

### Sensitive Data

Due to privacy concerns, precise geographic coordinates and detailed location information are not publicly available. Researchers requiring access to precise coordinates for legitimate scientific purposes should contact the corresponding author.

## Requirements

### Software Dependencies

**Sequence Processing:**
- Stacks v2.x (genotyping)
- BCFtools (VCF manipulation)
- PLINK v1.9 (data formatting)
- BWA-MEM (read mapping)
- SAMtools (BAM processing)

**R Packages:**
```r
# Population genetics
vcfR, adegenet, hierfstat, strataG, poppr, mmod, ape, pegas

# Multivariate statistics  
ade4, MASS

# Landscape genetics
ResistanceGA, corMLPE

# ABC
abc, nnet, randomForest, caret

# Spatial analysis
raster, sp, sf, gdistance

# Visualization and data handling
tidyverse, ggplot2, patchwork, RColorBrewer
```

**Python/Other:**
- Custom holoSimCell simulation framework (included in `code/ABC/`)

### Computing Requirements

ABC simulations require high-performance computing resources

## License

Code is provided as-is for reproducibility and educational purposes.