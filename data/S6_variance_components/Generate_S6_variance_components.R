# =============================================================================
# Generate_S6_variance_components.R
#
# Supplementary File 6: Variance components from multilevel meta-analytic models
#
# This script fits the same intercept-only rma.mv() models used in
# Cleanish_metaanalysis.Rmd and extracts the three sigma² variance components
# (study-level, subject-level, body-habitat-level) for each model. Results are
# written to data/S6_variance_components/variance_components.csv.
#
# INPUT:  All files in data/results/ (pre-computed; use the Step 5 shortcut)
# OUTPUT: data/S6_variance_components/variance_components.csv
#
# Run from the repository root, or open the .Rproj if present.
# Requires: metafor (GitHub: wviechtb/metafor), tidyverse, here
# =============================================================================

library(metafor)
library(tidyverse)
library(here)

# ---------------------------------------------------------------------------
# 1. Load data  (mirrors the "Data" chunk in Cleanish_metaanalysis.Rmd)
# ---------------------------------------------------------------------------

parent_dir <- 'C:/Users/liamn/Desktop/colonization-persistence_repo_management/'
results_dir  <- 'data/results'
metadata_dir <- 'data/S4_subject_metadata'

my_data.0001 <- read.csv(file.path(metadata_dir, "final_biome_list_true.0001.csv"))
my_data.0001$true_biome[my_data.0001$true_biome == "respiratory_oral"] <- "oral"

stability_data   <- read.csv(file.path(results_dir, "stability_df.csv"))
merged_with_slope <- read.csv(file.path(results_dir, "merged_metrics.csv"))
merged_data.0001  <- inner_join(my_data.0001, merged_with_slope, by = join_by(ID == ID))

amount_of_families <- read.csv(file.path(results_dir, "distinct_families.csv"))
merged_data.0001   <- inner_join(merged_data.0001, amount_of_families, by = join_by(ID == ID))

pc_slope_data    <- read.csv(file.path(results_dir, "pc_slopes.csv"))
merged_data.0001 <- inner_join(merged_data.0001, pc_slope_data, by = join_by(ID == ID))

fava_data        <- read.csv(file.path(results_dir, "fava_values.csv"))
merged_data.0001 <- inner_join(merged_data.0001, fava_data, by = join_by(ID == ID))

# Slim to the columns used in model fitting
no_rho_no_slope <- merged_data.0001[, c(2, 11, 15:26, 31:39, 47)]
no_rho_no_slope$pub_year_simple <- no_rho_no_slope$pub_year - mean(no_rho_no_slope$pub_year)

# Effect size data frames
family.0001   <- read.csv(file.path(results_dir, "merged_data_nobeta.0001.csv"))
family_no_art <- read.csv(file.path(results_dir, "merged_data_nobeta_no_art.csv"))
genus.0001    <- read.csv(file.path(results_dir, "merged_gen.0001_final.csv"))
genus_no_art  <- read.csv(file.path(results_dir, "merged_gen_no_art_final.csv"))

# Merge with subject metadata and square SE (matches Cleanish_metaanalysis.Rmd)
prep_dataset <- function(es_df, meta_df) {
  df <- inner_join(es_df, meta_df, by = join_by(ID == ID))
  df <- df %>% mutate(study_ID = ifelse(study_ID == "MT5147", "M5147", study_ID))
  df$ste <- df$ste^2
  df
}

family.0001_complete   <- prep_dataset(family.0001,   no_rho_no_slope)
family_no_art_complete <- prep_dataset(family_no_art, no_rho_no_slope)
genus.0001_complete    <- prep_dataset(genus.0001,    no_rho_no_slope)
genus_no_art_complete  <- prep_dataset(genus_no_art,  no_rho_no_slope)

datasets <- list(
  family.0001   = family.0001_complete,
  family.no_art = family_no_art_complete,
  genus.0001    = genus.0001_complete,
  genus.no_art  = genus_no_art_complete
)

# Apply minimum-taxa filter (>9 families or genera)
datasets_slim <- datasets
datasets_slim$family.0001   <- datasets_slim$family.0001[datasets_slim$family.0001$distinct_families > 9, ]
datasets_slim$family.no_art <- datasets_slim$family.no_art[datasets_slim$family.no_art$distinct_families > 9, ]
datasets_slim$genus.0001    <- datasets_slim$genus.0001[datasets_slim$genus.0001$distinct_genera > 9, ]
datasets_slim$genus.no_art  <- datasets_slim$genus.no_art[datasets_slim$genus.no_art$distinct_genera > 9, ]

# ---------------------------------------------------------------------------
# 2. Fit intercept-only multilevel models  (mirrors run_meta_model())
# ---------------------------------------------------------------------------

run_meta_model <- function(data, level = "family", effect_size_col = "z", variance_col = "v") {
  min_taxa_col  <- ifelse(level == "family", "distinct_families", "distinct_genera")
  data_filtered <- data[data[[min_taxa_col]] > 9, ]

  model <- rma.mv(
    yi     = data_filtered[[effect_size_col]],
    V      = data_filtered[[variance_col]],
    random = list(~1 | study_ID,
                  ~1 | subject_id,
                  ~1 | biome_ID),
    method = "REML",
    data   = data_filtered,
    sparse = TRUE
  )
  return(model)
}

meta_models <- list()

for (name in names(datasets_slim)) {
  dat   <- datasets_slim[[name]]
  level <- ifelse(grepl("genus", name), "genus", "family")

  meta_models[[paste0("mod_rho_",   name)]] <- run_meta_model(dat, level = level, effect_size_col = "z",     variance_col = "v")
  meta_models[[paste0("mod_slope_", name)]] <- run_meta_model(dat, level = level, effect_size_col = "slope", variance_col = "ste")

  cat("Fitted:", name, "\n")
}

# ---------------------------------------------------------------------------
# 3. Extract sigma² components
#
# rma.mv() returns sigma2 in the order the random effects were specified:
#   sigma2[1] = study_ID   (between-study)
#   sigma2[2] = subject_id (between-subject, within-study)
#   sigma2[3] = biome_ID   (between-body-habitat, within-subject)
# ---------------------------------------------------------------------------

extract_sigma2 <- function(model_list) {
  do.call(rbind, lapply(names(model_list), function(name) {
    m   <- model_list[[name]]
    s2  <- m$sigma2          # named numeric vector of length 3
    tot <- sum(s2)

    data.frame(
      model              = name,
      taxon_level        = ifelse(grepl("genus", name), "genus", "family"),
      artifact_filter    = ifelse(grepl("no_art", name), "artifacts_removed", "0.01pct_threshold"),
      effect_type        = ifelse(grepl("rho", name), "rho_z", "slope_beta"),
      k                  = m$k,
      sigma2_study       = round(s2[1], 6),
      sigma2_subject     = round(s2[2], 6),
      sigma2_biome       = round(s2[3], 6),
      sigma2_total       = round(tot,   6),
      pct_study          = round(s2[1] / tot * 100, 2),
      pct_subject        = round(s2[2] / tot * 100, 2),
      pct_biome          = round(s2[3] / tot * 100, 2),
      stringsAsFactors   = FALSE
    )
  }))
}

sigma2_table <- extract_sigma2(meta_models)

## sigma2_table is the desired output and contains the variance amonng studies, subjects, and subject/biome combinations.
