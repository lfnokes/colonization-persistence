# Colonization-Persistence Trade-offs in the Human Microbiome (Currently being QC'd)

**Author:** Liam F. Nokes  
**Institution:** Dartmouth College, Environmental Studies Senior Honors Thesis  
**Advisor:** Dr. Bala Chaudhary  
**Published:** June 2025  
**Citation:** Nokes, Liam F., "Colonization-Persistence Trade-offs in the Human Microbiome" (2025). *Environmental Studies Senior Theses*. 13. https://digitalcommons.dartmouth.edu/environmental_studies_senior_theses/13

---

## Overview

This repository contains all code, data, and figures supporting the thesis. We test for a colonization-persistence trade-off across human-associated microbial communities using a participant-level meta-analysis of longitudinal microbiome datasets from the [MGnify database](https://www.ebi.ac.uk/metagenomics/). Island biogeography-based models are applied to calculate effective colonization and persistence rates for microbial families and genera across 606 individuals in 26 studies spanning gut, skin, oral, respiratory, and vaginal microbiomes.

---

## Repository Structure

```
colonization-persistence/
├── scripts/                         # All analysis scripts (run in order listed below)
│   ├── API_query.R                  # Step 1: Study identification via MGnify API
│   ├── Data_metadata_integration.R  # Step 2: Data download, cleaning, and metadata matching
│   ├── Combining_all_data.R         # Step 3: Combining cleaned data across all studies
│   ├── Cleanish_results.Rmd         # Step 4: C/E rate calculation and effect size computation
│   └── Cleanish_metaanalysis.Rmd    # Step 5: Meta-analysis, meta-regression, and all figures
│
├── data/
│   ├── raw_inputs/                  # Raw input files for Data_metadata_integration.R
│   │   ├── ENA_metadata_tables/     # Sample metadata fetched from ENA (35 files)
│   │   ├── MGYS00005378/            # Manually downloaded OTU table and metadata for MGYS5378
│   │   ├── MGYS00006120/            # Manually downloaded taxonomic profiles for MGYS6120
│   │   ├── MGYS00004729_meta.xlsx   # Publication metadata for MGYS4729
│   │   ├── meta0308.csv             # Publication metadata for MGYS0308
│   │   ├── tax_incomplete1070.xlsx  # Publication taxonomy for MGYS1070
│   │   └── tax_1818.txt             # Publication taxonomy for MGYS1818
│   ├── S1_study_list/               # Supplementary File 1: Study list and filtering info
│   ├── S2_raw_taxonomy/             # Supplementary File 2: Cleaned taxonomy tables per study
│   ├── S3_artifact_taxa/            # Supplementary File 3: Removed low-variability taxa
│   ├── S4_subject_metadata/         # Supplementary File 4: Subject metadata and disease status
│   ├── S5_leave_one_out/            # Supplementary File 5: Leave-one-out analysis data
│   ├── S6_variance_components/      # Supplementary File 6: Multilevel model variance components
│   ├── S7_family_rates/             # Supplementary File 7: Family-level C/E rates and ranks
│   └── results/                     # Key intermediate and final result CSVs
│
└── figures/
    ├── S5_leave_one_out/            # Supplementary File 5: Leave-one-out sensitivity figures
    └── main_analysis/               # Final published figures
```

---

## Supplementary File Guide

| File | Location | Description |
|------|----------|-------------|
| **S1** | `data/S1_study_list/` + `scripts/API_query.R` | Full list of 26 included studies with sequencing platform, bioinformatic pipeline, and MGnify accession IDs. `API_query.R` contains the keyword search and filtering logic used to identify candidate studies. |
| **S2** | `data/S2_raw_taxonomy/` + `scripts/Data_metadata_integration.R` + `scripts/Combining_all_data.R` | Cleaned taxonomy tables (`usable_tax_[MGYS_ID].csv`) for each study, produced after downloading MAE objects from MGnify, applying the 0.01% abundance filter, and matching to subject/timepoint metadata. Raw MAE objects are not included due to file size but are freely available via the MGnify API (see Data Availability below). |
| **S3** | `data/S3_artifact_taxa/all_artifact_taxa.csv` | Complete list of taxa removed from each subject/body-habitat combination due to limited presence/absence variability (fewer than 2 transitions). |
| **S4** | `data/S4_subject_metadata/final_biome_list_true.0001.csv` | Full subject list with associated body habitat, disease status, antibiotic use, probiotic use, infant status, and pregnancy status for all 606 participants. |
| **S5** | `data/S5_leave_one_out/` + `figures/S5_leave_one_out/` + `scripts/Cleanish_metaanalysis.Rmd` | Leave-one-out sensitivity analysis conducted at both the subject/body-habitat scale and the study scale. CSV files contain effect sizes with each unit removed; figures (`leaveout_plot_1.pdf` through `leaveout_plot_16.pdf`) show the resulting distributions. See the figure legend in `Cleanish_metaanalysis.Rmd` for the mapping of plot numbers to model types. |
| **S6** | `data/S6_variance_components/` | Variance components (sigma²) from the three-level `rma.mv()` multilevel meta-analytic models, showing the relative contribution of study-level, subject-level, and body-site-level variance to overall heterogeneity. **To be generated:** run the `extract_sigma2` code chunk in `scripts/Cleanish_metaanalysis.Rmd`. |
| **S7** | `data/S7_family_rates/` | Colonization and extinction rates and standardized ranks for every microbial family (`fam_ce_values_no_art.csv`) and genus (`gen_ce_values_no_art.csv`) in every subject/body-habitat combination included in the analysis. |

---

## How to Reproduce Results

### Prerequisites

Install the following R packages:

```r
# From CRAN
install.packages(c("tidyverse", "ggplot2", "island", "FAVA"))

# From Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("MultiAssayExperiment")

# From GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("wviechtb/metafor")
devtools::install_github("daniel1noble/orchaRd")
devtools::install_github("EBI-Metagenomics/MGnifyR")
```

> **Note on package versions:** All analyses were run using the most current versions of each package available in April–May 2025. Exact versions were not recorded at the time. If you encounter unexpected errors, this is the most likely source — consider checking GitHub release histories for breaking changes around that period. A future update to this repository will include a `sessionInfo()` output or `renv.lock` file for full reproducibility.

### Running the Pipeline

Run scripts in the following order. It may make sense to skip over certain steps that could be affected by changes in MGnify API structure (Steps 1-3).

1. **`API_query.R`** — Queries the MGnify API to identify longitudinal human microbiome studies. Outputs a filtered list of candidate studies. Requires an internet connection and the `MGnifyR` package. **Note:** most users should skip this step and start at Step 2 (see below).

2. **`Data_metadata_integration.R`** — Downloads MAE objects for each of the 26 included studies from MGnify and integrates with subject/timepoint metadata. All supporting input files (ENA metadata tables, manually downloaded OTU tables and taxonomic profiles for studies that required non-API retrieval) are included in `data/raw_inputs/`. Note: MAE downloads are large and were originally run on Dartmouth's Discovery HPC cluster; this script may take several hours on a standard machine. **Note:** if MAE downloads prevent the execution of this script, skip Steps 2 and 3 and start with Step 4 using the `usable_tax**.csv` files included in `data/results/`.

3. **`Combining_all_data.R`** — Applies the 0.01% abundance filter, converts to presence/absence, splits by subject and body site, and compiles the cleaned taxonomy tables found in `data/S2_raw_taxonomy/`.

4. **`Cleanish_results.Rmd`** — Fits the island biogeography model (`irregular_single_dataset` from the `island` package) to each subject/body-habitat combination at the family and genus level to estimate colonization and extinction rates. Also calculates Spearman correlation and regression slope between log-transformed colonization and persistence rates. Outputs the C/E rate tables and effect size CSVs in `data/results/`.

5. **`Cleanish_metaanalysis.Rmd`** — Runs all meta-analyses, meta-regressions, bias tests, leave-one-out analyses, and generates all figures. Input files are in `data/results/`. Outputs figures to `figures/`.

> **Recommended starting point:** For most purposes, we recommend starting at **Step 2** (`Data_metadata_integration.R`) rather than Step 1. The API query (Step 1) was used to identify candidate studies at the time of this analysis, but MGnify's API structure and study contents may change over time. All raw input data needed to run Step 2 onward — including ENA sample metadata and manually retrieved files for studies that required non-API access — are archived in `data/raw_inputs/`, making Steps 2–5 fully reproducible without relying on the current state of the API.
>
> **Shortcut:** If you want to reproduce only the meta-analysis and figures without re-running the full data pipeline, start at Step 5 (`Cleanish_metaanalysis.Rmd`) using the pre-computed files already provided in `data/results/`.

---

## Data Availability

Raw microbiome sequencing data and taxonomic profiles for all 26 studies are publicly available through the [MGnify database](https://www.ebi.ac.uk/metagenomics/) (European Nucleotide Archive). Study accession IDs are listed in `data/S1_study_list/FINAL_Study_list.csv`. Data can be downloaded programmatically using the [MGnifyR](https://github.com/EBI-Metagenomics/MGnifyR) R package, as demonstrated in `scripts/Data_metadata_integration.R`.

Raw MAE (MultiAssayExperiment) R objects are not included in this repository due to file size, but can be regenerated by running `scripts/Data_metadata_integration.R`. All other inputs required by that script — including ENA sample metadata (fetched via the ENA API) and manually downloaded files for two studies (MGYS00005378 and MGYS00006120) where the standard MGnifyR download was not available — are archived in `data/raw_inputs/`.

---

## Key R Packages

| Package | Source | Purpose |
|---------|--------|---------|
| `island` | CRAN | Colonization/extinction rate estimation |
| `metafor` | GitHub (`wviechtb/metafor`) | Multilevel meta-analysis and meta-regression |
| `orchaRd` | GitHub (`daniel1noble/orchaRd`) | Orchard plot visualization |
| `FAVA` | CRAN | Temporal compositional variability (F_ST-based) |
| `MGnifyR` | GitHub (`EBI-Metagenomics/MGnifyR`) | MGnify API data access |
| `tidyverse` | CRAN | Data manipulation |
| `MultiAssayExperiment` | Bioconductor | MAE object handling |

---

## Contact

Liam F. Nokes — liam.f.nokes.25@dartmouth.edu  
For questions about the data or code, please open a GitHub issue.
