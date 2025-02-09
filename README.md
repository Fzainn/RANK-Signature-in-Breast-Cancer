# RANK Overexpression Signature Analysis

## OverView
This project contains R scripts for analyzing RANK overexpression signatures and their impact on survival outcomes in breast cancer patients. The analysis is conducted using two datasets:
  1. *METABRIC Dataset* - Examines the prognostic value of the RANK signature in breast cancer survival.
  2. *DCIS Dataset* - Investigates the role of the RANK signature in predicting disease recurrence in ductal carcinoma in situ (DCIS).

# Requirements
## R Libraries
Ensure the following R libraries are installed before running the scripts:
  ```r
install.packages(c("dplyr", "survival", "survminer", "readxl", "tidyr", "pheatmap", "ggplot2"))`

# Load the libraries:
    library(dplyr)
    library(survival)
    library(survminer)
    library(readxl)
    library(tidyr)
    library(pheatmap)
    library(ggplot2)
```


Data Source

* METABRIC Data
    * Gene expression data: data_mrna_illumina_microarray.txt 
    * Clinical data: data/data_clinical_patient.txt
    * [Download](https://www.synapse.org/Synapse:syn1688369/wiki/27311)
 
* DCIS Data
    * Gene expression data: RAHBT_LCM_reads_matrix_Epi265_Str196_VSTnorm.csv
    * Clinical data: DCISClinical.tsv
    * [Download](https://humantumoratlas.org/)
      
* RANK Signature Gene List : mmc2.xlsx
* [Download](https://www.ncbi.nlm.nih.gov/gds?linkname=pubmed_gds&from_uid=34004159)


# Workflow
## METABRIC Analysis
1. Load METABRIC gene expression and clinical data.
2. Remove duplicate genes and missing values.
3. Filter data to include only RANK signature genes.
4. Compute RANK signature scores per patient.
5. Merge gene expression data with clinical data.
6. Define high and low RANK signature score groups.
7. Perform survival analysis:
   * Kaplan-Meier survival curves
   * Cox proportional hazards regression
8. Save results:
   * Merged dataset: merged_data.csv
   * Kaplan-Meier plot: KM_plot_RANK_signature.png
   * Cox model results: cox_results_signature_score.csv

## DCIS Analysis
1. Load DCIS gene expression and clinical data.
2. Process and filter gene expression data for RANK signature genes.
3. Compute RANK signature scores per patient.
4. Merge gene expression data with clinical data.
5. Define high and low RANK signature score groups.
6. Perform recurrence analysis:
   * Kaplan-Meier survival curves.
7. Save results:
   * Merged dataset: merged_data_DCIS.csv
   * Kaplan-Meier plot: Combined_KM_Plot_DCIS_RANK_Signature.png

