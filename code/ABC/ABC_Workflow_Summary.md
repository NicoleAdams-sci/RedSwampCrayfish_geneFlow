# Red Swamp Crayfish ABC Analyses - Workflow Summary

## Overview

This project uses Approximate Bayesian Computation (ABC) to infer population genetic parameters and demographic history for Red Swamp Crayfish populations in SE Michigan, USA. The workflow involves simulation, model selection, parameter estimation, and cross-validation.

## Key Data Files

- **modList.txt**: List of demographic models to test
- **rsc_100k_subset_new.csv**: Combined simulation results (100k per model)
- **rsc.empirical.stats.csv**: Observed summary statistics from empirical data

---

## Workflow Stage 1: ABC Simulations

### 1. RSC_ABCsims_20250722_anon.R

**Location**: `code/ABC/RSC_ABCsims_20250722_anon.R`

**Purpose**: ABC simulation engine with demographic models

**Dependencies**:
- Event files: `*_sourceSink.txt` (demographic scenarios)
- Helper functions: `holoSimCell_helper-functions.R`, `holoSimCell_additional_stats_functions_NEA.R`
- Missing data pattern: `jy22yLoA_adults_baseFilter_imiss99_dp5miss50_imiss75_hz.ab_1pRt_100kb_rmdups_orderedMissMat.txt`

**Inputs**: Command line args: `node`, `nreps`, `runBy`

**Outputs**: `paramYstats_[runBy].[node].csv`

**Function**:
- Simulates genetic data under 10 different demographic models
- Calculates summary statistics (FST, heterozygosity, private alleles, etc.)
- Runs on SLURM cluster as array jobs

**Run Command**:
```bash
sbatch rsc_reUp.sbatch
# Once reach desired number of replicates across all models, you must manually kill the rsc_reUp.sbatch script 
```

**Note**: The script is called via `ABCsims.sbatch` or `rsc_reUp.sbatch`, which runs it as an array job, each performing multiple replications per model.

---

## Workflow Stage 2: Data Processing & Combination

### 2. subset_combine_old_new_models.sh

**Location**: `sensitive/ABC/subset_combine_old_new_models.sh`

**Purpose**: Combines old and new simulation results

**Output**: `rsc_100k_subset_new.csv`

### 3. extract_params4plots.sh

**Location**: `analysis/ABC/extract_params4plots.sh`

**Purpose**: Extracts parameter columns for plotting

**Output**: `rsc_100k_params4PEplots_new.csv`

---

## Workflow Stage 3: Model Selection

### 4. Random Forest Model Selection

**Script**: `rsc_RF_pred_subset.R`

**Location**: `analysis/ABC/rsc_RF_pred_subset.R`

**Dependencies**:
- `rsc_100k_subset_new.csv`
- `rsc.empirical.stats.csv`

**Outputs**:
- `rsc_RF_forest_100ksubset.rda`
- `RF_rscpred_100ksubset.Rws`

**Run Command**:
```bash
sbatch rsc_RF_pred_subset.sh
```

### 5. Multinomial Logistic Model Selection

**Script**: `rsc_modsel_subset_mnlog.R`

**Location**: `analysis/ABC/rsc_modsel_subset_mnlog.R`

**Dependencies**: Same as RF script

**Outputs**: `rsc_modsel_100ksubset_MNLOG-[tol].rda`

**Run Command**:
```bash
sbatch rsc_mnlog.sh [tolerance]
```

### 6. Neural Network Model Selection

**Script**: `rsc_modsel_subset_nnetTest.R`

**Location**: `analysis/ABC/rsc_modsel_subset_nnetTest.R`

**Dependencies**: Same as RF script

**Outputs**: `rsc_modsel_subset_NNET-[tol].rda`

**Run Command**:
```bash
sbatch rsc_nnet.sbatch [tolerance]
```

### 7. Model Selection Summary Table

**Script**: `mnlog_nnet_modsel_table.R`

**Location**: `analysis/ABC/mnlog_nnet_modsel_table.R`

**Purpose**: Combine model selection results into summary table

**Dependencies**: All `.rda` files from model selection

**Outputs**: `model_selection_results_table.csv`

**Run Command**:
```bash
Rscript mnlog_nnet_modsel_table.R
```

---

## Workflow Stage 4: Cross-Validation

### 8. Multinomial Logistic Cross-Validation

**Script**: `rsc_CV_mnlog.sh`

**Location**: `analysis/ABC/rsc_CV_mnlog.sh`

**Purpose**: Submit multinomial logistic cross-validation jobs

**Dependencies**:
- `modList.txt`
- `rsc_CV_subset.R`

**Outputs**: `rsc_100k_CV_mnlogistic_[MODEL]-[tol]-[node].csv` files

**Run Command**:
```bash
bash rsc_CV_mnlog.sh
# This creates and submits individual sbatch scripts for each model
# Example: sbatch -a 1-10 rsc_100k_CV_mnlog_[MODEL].sh
```

### 9. Neural Network Cross-Validation

**Script**: `rsc_CV_nnet.sh`

**Location**: `analysis/ABC/rsc_CV_nnet.sh`

**Purpose**: Submit neural network cross-validation jobs

**Dependencies**:
- `modList.txt`
- `rsc_CV_subset.R`

**Outputs**: `rsc_100k_CV_neuralnet_[MODEL]-[tol]-[node].csv` files

**Run Command**:
```bash
bash rsc_CV_nnet.sh
# This creates and submits individual sbatch scripts for each model
# Example: sbatch -a 1-20 rsc_100k_CV_nnet_[MODEL].sh
```

### 10. Core Cross-Validation Engine

**Script**: `rsc_CV_subset.R`

**Location**: `analysis/ABC/rsc_CV_subset.R`

**Purpose**: Core cross-validation engine (called by both shell scripts above)

**Dependencies**: `rsc_100k_subset_new.csv`

### 11. Combine Cross-Validation Results

**Script**: `combine_cv_files.sh`

**Location**: `analysis/ABC/combine_cv_files.sh`

**Purpose**: Combines all CV results into single file

**Output**: `rsc_100k_CV_combo_new.csv`

**Run Command**:
```bash
bash combine_cv_files.sh
```

---

## Workflow Stage 5: Parameter Estimation

### 12. Parameter Estimation

**Script**: `rsc_PE_subset.R`

**Location**: `analysis/ABC/rsc_PE_subset.R`

**Purpose**: Parameter estimation for best models

**Dependencies**:
- `rsc_100k_subset_new.csv`
- `rsc.empirical.stats.csv`

**Outputs**: `rsc_PE_subset_new_[model]_[method]-[tol].rda`

**Run Command**:
```bash
sbatch rsc_PE.sbatch [model] [method] [tolerance]
```

### 13. Parameter Estimation Cross-Validation

**Script**: `rsc_CV4PE.R`

**Location**: `analysis/ABC/rsc_CV4PE.R`

**Purpose**: Cross-validation for parameter estimation accuracy

**Dependencies**: `rsc_100k_subset_new.csv`

**Outputs**: `rsc_CV4PE_new_[model]_[method]_[tol].csv`

**Run Command**:
```bash
sbatch rsc_CV4PE.sbatch [model] [method] [tolerance]
```

---

## Workflow Stage 6: Visualization & Analysis

### 14. Cross-Validation Plots

**Script**: `abc_CV_plots.R`

**Location**: `analysis/ABC/abc_CV_plots.R`

**Purpose**: Plot cross-validation results

**Dependencies**: `rsc_100k_CV_combo_new.csv`

**Outputs**: Model selection accuracy plots

**Run Command**:
```bash
Rscript abc_CV_plots.R
```

### 15. Parameter Estimation CV Plots

**Script**: `abc_CV4PE_plots.R`

**Location**: `analysis/ABC/abc_CV4PE_plots.R`

**Purpose**: Plot parameter estimation cross-validation

**Dependencies**: `rsc_CV4PE_*.csv` files

**Outputs**: Parameter estimation accuracy plots

**Run Command**:
```bash
Rscript abc_CV4PE_plots.R
```

### 16. Posterior Distribution Plots

**Script**: `abc_posteriors_plots.R`

**Location**: `analysis/ABC/abc_posteriors_plots.R`

**Purpose**: Plot posterior parameter distributions

**Dependencies**: `rsc_PE_subset_*.rda` files

**Outputs**: Posterior distribution plots

**Run Command**:
```bash
Rscript abc_posteriors_plots.R
```

### 17. Confusion Matrix

**Script**: `abc_confusion_matrix.Rmd`

**Location**: `analysis/ABC/abc_confusion_matrix.Rmd`

**Purpose**: Access confusion matrix for model selection

**Dependencies**: `rsc_RF_forest_100ksubset.rda`

---

## Dependency Graph

```
RSC_ABCsims_20250722_anon.R (via rsc_reUp.sbatch)
  ↓ (produces individual CSV files)
subset_combine_old_new_models.sh
  ↓ (produces rsc_100k_subset_new.csv)
  │
  ├── Model Selection:
  │   ├── rsc_RF_pred_subset.R (via rsc_RF_pred_subset.sh)
  │   ├── rsc_modsel_subset_mnlog.R (via rsc_mnlog.sh)
  │   ├── rsc_modsel_subset_nnetTest.R (via rsc_nnet.sbatch)
  │   └── mnlog_nnet_modsel_table.R
  │
  ├── Cross-Validation:
  │   ├── rsc_CV_mnlog.sh → rsc_CV_subset.R
  │   ├── rsc_CV_nnet.sh → rsc_CV_subset.R
  │   └── combine_cv_files.sh → rsc_100k_CV_combo_new.csv
  │
  ├── Parameter Estimation:
  │   ├── rsc_PE_subset.R (via rsc_PE.sbatch)
  │   └── rsc_CV4PE.R (via rsc_CV4PE.sbatch)
  │
  └── Visualization:
      ├── abc_CV_plots.R
      ├── abc_CV4PE_plots.R
      ├── abc_posteriors_plots.R
      └── abc_confusion_matrix.Rmd
```

---

## Typical Execution Order

1. **Run simulations**: 
   ```bash
   sbatch rsc_reUp.sbatch
   ```

2. **Combine simulation outputs**: 
   ```bash
   bash subset_combine_old_new_models.sh
   bash extract_params4plots.sh rsc_100k_subset_new.csv rsc_100k_params4PEplots_new.csv
   ```

3. **Model selection**: 
   ```bash
   sbatch rsc_RF_pred_subset.sh
   sbatch rsc_mnlog.sh [tolerance]
   sbatch rsc_nnet.sbatch [tolerance]
   Rscript mnlog_nnet_modsel_table.R
   ```

4. **Cross-validation**: 
   ```bash
   bash rsc_CV_mnlog.sh  # Creates and submits array jobs
   bash rsc_CV_nnet.sh   # Creates and submits array jobs
   ```

5. **Combine CV results**: 
   ```bash
   bash combine_cv_files.sh
   ```

6. **Parameter estimation**: 
   ```bash
   sbatch rsc_PE.sbatch [model] [method] [tolerance]
   sbatch rsc_CV4PE.sbatch [model] [method] [tolerance]
   ```

7. **Generate plots and summaries**:
   ```bash
   Rscript abc_CV_plots.R
   Rscript abc_CV4PE_plots.R
   Rscript abc_posteriors_plots.R
   Rscript -e "rmarkdown::render('abc_confusion_matrix.Rmd')"
   ```

---

## Repository Structure

```
RedSwampCrayfish_geneFlow/
├── code/
│   └── ABC/
│       ├── RSC_ABCsims_20250722_anon.R
│       ├── rsc_reUp.sbatch
│       ├── ABCsims.sbatch
│       └── [helper functions and event files]
├── analysis/
│   └── ABC/
│       ├── [all analysis and plotting scripts]
│       └── Red_Swamp_Crayfish_ABC_Workflow_Summary.pdf
├── data/
│   └── ABC/
│       ├── rsc.empirical.stats.csv
│       └── modList.txt
├── output/
│   └── ABC_output/
│       └── [all generated output files]
└── sensitive/
    └── ABC/
        └── subset_combine_old_new_models.sh
```

---


## Notes

- All SLURM scripts use the MSU HPCC environment with R/4.2.2
- Cross-validation shell scripts (`rsc_CV_mnlog.sh`, `rsc_CV_nnet.sh`) are wrappers that create and submit individual sbatch scripts for each model
- The `rsc_reUp.sbatch` script is the primary way to run ABC simulations
- File paths in scripts reference the MSU HPCC file system structure
