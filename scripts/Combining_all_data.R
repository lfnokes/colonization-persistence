#####This script is for creating one big data structure with a list of all
#####studies, subjects, biomes, with OTU tables ordered by timepoint. 
#####This data structure will be used to analyze all of the data all at once
##Loading some packages I might need
library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)
library(strucchange)
library(readr)
library(rootSolve)
library(readxl)
library(MGnifyR)
library(biomformat)
library(island)
library(readxl)
library(stringr)
setwd("Studies_scripts")

##Bringing in all of the data
usable0278 <- read.csv("usable_tax0278.csv")
usable0308 <- read.csv("usable_tax0308.csv")
usable1070 <- read.csv("usable_tax1070.csv")
usable1278 <- read.csv("usable_tax1278.csv")
usable1818 <- read.csv("usable_tax1818.csv")
usable2184 <- read.csv("usable_tax2184.csv")
usable3619 <- read.csv("usable_tax3619.csv")
usable3659 <- read.csv("usable_tax3659.csv")
usable4708 <- read.csv("usable_tax4708.csv")
usable4729 <- read.csv("usable_tax4729.csv")
usable5142 <- read.csv("usable_tax5142.csv")
usable5147 <- read.csv("usable_tax5147.csv")
usable5176 <- read.csv("usable_tax5176.csv")
usable5196 <- read.csv("usable_tax5196.csv")
usable5261 <- read.csv("usable_tax5261.csv")
usable5378 <- read.csv("usable_tax5378.csv")
usable5791 <- read.csv("usable_tax5791.csv")
usable5816 <- read.csv("usable_tax5816.csv")
usable5860 <- read.csv("usable_tax5860.csv")
usable6120 <- read.csv("usable_tax6120.csv")
usable6263 <- read.csv("usable_tax6263.csv")
usable6273 <- read.csv("usable_tax6273.csv")
usable6469 <- read.csv("usable_tax6469.csv")
usable6530 <- read.csv("usable_tax6530.csv")
usable6638 <- read.csv("usable_tax6638.csv")

#the one for 1068 is already split up
split_tax1068 <- read_rds("split_tax1068.rds")

###Now, I am going to split these all up into lists by subject and biome.
###To do this, I have to make sure the taxonomies are somewhat standardized and in fact are the rownames of the OTU tables
###I like the format of sk__Bacteria;p__...;c__... etc. so I am going to try making everything like that

##0278
rownames(usable0278) <- usable0278$X
standardize_taxonomy <- function(taxa_list) {
  tax_levels <- c("k", "p", "c", "o", "f", "g", "s")  # up to species if needed
  
  sapply(taxa_list, function(taxa) {
    if (tolower(taxa) == "unusigned" || taxa == "") {
      return("Unassigned")
    }
    
    parts <- unlist(strsplit(taxa, ":"))
    # Make sure it doesn't exceed 7 levels
    parts <- parts[1:min(length(parts), length(tax_levels))]
    
    # Fill missing levels with empty strings
    parts <- c(parts, rep("", length(tax_levels) - length(parts)))
    
    paste(paste0(tax_levels, "__", parts), collapse = ";")
  })
}

new_rownames0278 <- standardize_taxonomy(rownames(usable0278))
rownames(usable0278) <- new_rownames0278
usable0278 <- usable0278[-44,]
##Done

##0308
rownames(usable0308) <- usable0308$X
new_rownames0308 <- standardize_taxonomy(rownames(usable0308))
rownames(usable0308) <- new_rownames0308
usable0308 <- usable0308[-107,]
##Done

##1070
rownames(usable1070) <- usable1070$X
new_rownames1070 <- paste0("k__;", rownames(usable1070))
rownames(usable1070) <- new_rownames1070
##Done

##1278
rownames(usable1278) <- usable1278$X
##Done

##1818
rownames(usable1818) <- usable1818$X
##Done

##2184
rownames(usable2184) <- usable2184$X
##Done

##3619
standardize_from_phylum_down_full <- function(taxa_list) {
  tax_levels <- c("p", "c", "o", "f", "g", "s")  # "st" gets handled separately
  
  sapply(taxa_list, function(t) {
    parts <- unlist(strsplit(t, "_"))
    
    # Handle "Candidatus_Something" as one unit
    if (length(parts) >= 2 && parts[1] == "Candidatus") {
      parts[1:2] <- paste(parts[1], parts[2], sep = "_")
      parts <- parts[-2]
    }
    
    # First 6 parts go to p..s
    tax_core <- head(parts, length(tax_levels))
    tax_core <- c(tax_core, rep("", length(tax_levels) - length(tax_core)))  # pad if needed
    
    # Remaining parts (after 6) go into "st"
    if (length(parts) > length(tax_levels)) {
      st <- paste(parts[(length(tax_levels) + 1):length(parts)], collapse = "_")
    } else {
      st <- ""
    }
    
    # Combine with prefixes
    full <- c(paste0(tax_levels, "__", tax_core), paste0("st__", st))
    paste(full, collapse = ";")
  })
}
new_rownames3619 <- standardize_from_phylum_down_full(usable3619$X)
rownames(usable3619) <- new_rownames3619
##Done

##3659
rownames(usable3659) <- usable3659$X
#column names are also a bit messed up since some have additional underscores in biome name
colnames(usable3659)
fix_colnames <- function(names_vec) {
  sapply(names_vec, function(name) {
    # Split on underscores
    parts <- strsplit(name, "_")[[1]]
    
    if (length(parts) == 4) {
      # Replace third underscore with period
      paste(parts[1:3], sep = "_", collapse = "_")
    } else {
      # Leave unchanged if there isn't a third underscore
      name
    }
  })
}
new_colnames3659 <- fix_colnames(colnames(usable3659))
colnames(usable3659) <- new_colnames3659
##Done

##4708
rownames(usable4708) <- usable4708$X
##Done

##4729
rownames(usable4729) <- usable4729$X
##Done

##5142
new_rownames5142 <- standardize_from_phylum_down_full(usable5142$X)
rownames(usable5142) <- new_rownames5142
##Done

##5147
new_rownames5147 <- standardize_from_phylum_down_full(usable5147$X)
rownames(usable5147) <- new_rownames5147
##Done

##5147 Malawi Twins
usable_malawi5147 <- read.csv("usable_malawi5147.csv")
new_rownames_malawi5147 <- standardize_from_phylum_down_full(usable_malawi5147$X)
rownames(usable_malawi5147) <- new_rownames_malawi5147

##5176
rownames(usable5176) <- usable5176$X
##Done

##5196
new_rownames5196 <- standardize_from_phylum_down_full(usable5196$X)
rownames(usable5196) <- new_rownames5196
##Done

##5261
new_rownames5261 <- standardize_from_phylum_down_full(usable5261$X)
rownames(usable5261) <- new_rownames5261
##Done

##5378
standardize_semicolon_taxa_with_otu <- function(taxa_list) {
  tax_levels <- c("k", "p", "c", "o", "f", "g", "s")  # standard ranks
  max_taxa_len <- length(tax_levels)
  
  sapply(taxa_list, function(t) {
    parts <- unlist(strsplit(t, ";"))
    
    otu_id <- parts[1]
    taxonomy <- parts[-1]
    
    # Replace "unclassified" with "" and trim
    taxonomy <- gsub("^unclassified$", "", taxonomy)
    
    # Pad taxonomy to 7 levels if shorter
    taxonomy <- head(c(taxonomy, rep("", max_taxa_len)), max_taxa_len)
    
    # Add prefixes
    tax_string <- paste(paste0(tax_levels, "__", taxonomy), collapse = ";")
    
    # Append OTU as ID
    paste0(tax_string, ";st__", otu_id)
  })
}
new_rownames5378 <- standardize_semicolon_taxa_with_otu(usable5378$X)
rownames(usable5378) <- new_rownames5378
##Done

##5791
rownames(usable5791) <- usable5791$X
##Done

##5816
rownames(usable5816) <- usable5816$X
##Done

##5860
standardize_weird_taxonomy <- function(taxa_list) {
  sapply(taxa_list, function(t) {
    # Split into underscore and semicolon parts (if present)
    split_parts <- unlist(strsplit(t, "_(?=Bacteria;|Archaea;)", perl = TRUE))
    
    # First part: underscore-separated taxonomy from k
    main_part <- split_parts[1]
    levels <- unlist(strsplit(main_part, "_"))
    
    # Taxonomic ranks starting from kingdom
    tax_levels <- c("k", "p", "c", "o", "f", "g")
    
    # Pad or trim the levels to match length
    levels <- head(c(levels, rep("", length(tax_levels))), length(tax_levels))
    
    # If second part exists, extract last field after semicolon
    if (length(split_parts) > 1) {
      semi_parts <- unlist(strsplit(split_parts[2], ";"))
      species <- tail(semi_parts, n = 1)
    } else {
      species <- ""
    }
    
    # Build the full taxonomy string
    full <- c(paste0(tax_levels, "__", levels), paste0("s__", species))
    paste(full, collapse = ";")
  })
}
new_rownames5860 <- standardize_weird_taxonomy(usable5860$X)
view(new_rownames5860)
##Eukaryotes are messing everything up and appear to be mostly artifacts
usable5860 <- usable5860[-c(355:362),]
new_rownames5860 <- new_rownames5860[-c(355:362)]
##weird pipeline created some duplicates that don't make sense, will combine them
usable5860$new_rownames <- new_rownames5860 
usable5860 <- usable5860[,-1]
unique(usable5860$new_rownames)
usable5860_agg <- aggregate(. ~ new_rownames, data = usable5860, FUN = sum)
rownames(usable5860_agg) <- usable5860_agg$new_rownames
usable5860 <- usable5860_agg
##Done

##6120
new_rownames6120 <- str_replace_all(usable6120$X, "\\|", ";")
rownames(usable6120) <- new_rownames6120
##Done

##6263
rownames(usable6263) <- usable6263$X
##Done

##6273
rownames(usable6273) <- usable6273$X
##Done

##6469
rownames(usable6469) <- usable6469$X
##Done

##6530
rownames(usable6530) <- usable6530$X
##Done

##6638
rownames(usable6638) <- usable6638$X
##Done

##Making all of it into a big list
standardized_list <- list(usable0278,
usable0308,
usable1070,
usable1278,
usable1818,
usable2184,
usable3619,
usable3659,
usable4708,
usable4729,
usable5142,
usable5147,
usable5176,
usable5196,
usable5261,
usable5378,
usable5791,
usable5816,
usable5860,
usable6120,
usable6263,
usable6273,
usable6469,
usable6530,
usable6638)
write_rds(standardized_list, "standardized_list.rds")

write_rds(usable_malawi5147, "standardized_malawi_twins.rds")

####Now, to split all of them up by biome and subject
big_table_list <- read_rds("standardized_list.rds")
#Need to add malawi twins
MT5147 <- read_rds("standardized_malawi_twins.rds")

view(big_table_list[[1]])
study_accessions <- c("M0278","M0308","M1070","M1278","M1818","M2184","M3619",
                      "M3659","M4708","M4729","M5142","M5147","M5176","M5196",
                      "M5261","M5378","M5791","M5816","M5860","M6120","M6263",
                      "M6273","M6469","M6530","M6638")
names(big_table_list) <- study_accessions
split_by_subject_and_biome <- function(df) {
  # Split all column names
  col_parts <- strsplit(colnames(df), "_")
  
  # Filter to only columns with exactly 3 parts
  valid_idx <- sapply(col_parts, length) == 3
  
  if (!any(valid_idx)) return(list())  # Return empty list if no valid columns
  
  col_parts <- col_parts[valid_idx]
  valid_cols <- colnames(df)[valid_idx]
  
  # Create a matrix of subject, timepoint, biome
  col_info <- do.call(rbind, col_parts)
  colnames(col_info) <- c("subject", "timepoint", "biome")
  
  # Group by subject-biome
  subject_biome <- paste0(col_info[, "subject"], "_", col_info[, "biome"])
  split_indices <- split(seq_along(valid_cols), subject_biome)
  
  # Create sub-dataframes
  lapply(names(split_indices), function(name) {
    df[, valid_cols[split_indices[[name]]], drop = FALSE]
  }) |> `names<-`(names(split_indices))
}

nested_table_list <- lapply(big_table_list, split_by_subject_and_biome)
#Malawi twins
MT5147_split <- split_by_subject_and_biome(MT5147)

flat_OG_list <- flatten_nested_list(nested_table_list)

filtered_nested_list <- lapply(nested_table_list, function(study) {
  Filter(function(df) ncol(df) > 9, study)
})
#Malawi twins
MT5147_filtered <- Filter(function(df) ncol(df) > 9, MT5147_split)

###I am noticing duplicates
#it appears that the duplicates happen when there is an exact overlap in the sample ID, date, and biome
#I have looked over the datasets where it happens and there are three causes:
#1. I misidentified the biome or number of biomes present, so it conflated two or more biomes (2184)
#2. The duplicates occur do to a sampling scheme that seems to resample on the same day sometimes, e.g. if there is an illness (or for no documented reason)
#3. there are actual systematic duplicates
# in the first case (2184) I need to go back to edit the metadata merge process
# in the second case, I am ignoring the duplicates for now because except in the case of 5147 there is no rhyme or reason to them (6273, 6530)
# in the third case, it makes sense to combine them by adding them together since every single sample has a duplicate (4729, 5816, 5860, 5147(in some cases))
# I am not worried about adding them together when there are duplicates of every sample because it is presence/absence
# when there are not duplicates of every sample, then some samples might get artificially inflated in diversity relative to others
# The solution to the duplicates problem is to create a function that adds dataframes componentwise when they contain the same timepoints for the same individual/biome combo


resolve_duplicates_in_study <- function(study_df_list) {
  # Extract base names (without trailing .1, .2, etc.)
  base_names <- sub("\\.\\d+$", "", names(study_df_list))
  
  # Group by base name
  grouped <- split(seq_along(study_df_list), base_names)
  
  # Output list
  cleaned_list <- list()
  
  for (base in names(grouped)) {
    indices <- grouped[[base]]
    dfs <- study_df_list[indices]
    
    if (length(dfs) == 1) {
      # Only one entry — keep as-is
      cleaned_list[[base]] <- dfs[[1]]
    } else {
      # Multiple versions — check for match
      ref_df <- dfs[[1]]
      combined <- ref_df
      for (i in 2:length(dfs)) {
        curr_df <- dfs[[i]]
        if (ncol(curr_df) == ncol(ref_df)) {
          # Same shape — add values
          combined <- combined + curr_df
        } else {
          # Different shape — save separately as extra
          cleaned_list[[paste0(base, ".extra")]] <- curr_df
        }
      }
      # Save the combined dataframe as the base
      cleaned_list[[base]] <- combined
    }
  }
  
  cleaned_list
}

cleaned_nested_list <- lapply(filtered_nested_list, resolve_duplicates_in_study)
write_rds(cleaned_nested_list, "cleaned_nested_list.rds")

#Malawi twins
MT5147_dup <- resolve_duplicates_in_study(MT5147_filtered)
write_rds(MT5147_dup, "cleaned_MT5147_twins.rds")

####OK! we have a list of dataframes with all of the data that we are going to analyze (plus some labeled "extra" data that needs to me matched with metadata later)
####Now, we want to convert this list into presence/absence
list_for_conversion <- read_rds("cleaned_nested_list.rds")

#Malawi twins
MT5147_for_conversion <- read_rds("cleaned_MT5147_twins.rds")

make_numeric <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(replace(x, is.na(x), 0)))
  return(df)
}
numeric_list <- lapply(list_for_conversion, function(study) {
  lapply(study, make_numeric)
})
#Malawi twins
MT5147_numeric <- lapply(MT5147_for_conversion, make_numeric)
##
############################################################
#before making a presence-absence list I should do a prevalence filtered pa_list
#going to filter out OTUs with less than 0.01% relative abundance
filter_by_abundance <- function(otu_table, min_rel_abundance = 0.0001) {
  # Compute row sums and total sum of abundance values
  row_totals <- rowSums(otu_table, na.rm = FALSE)
  table_total <- sum(row_totals)
  
  # Compute relative abundance for each OTU
  rel_abundance <- row_totals / table_total
  
  # Filter rows by relative abundance threshold
  otu_table[rel_abundance >= min_rel_abundance, ]
}

numeric_list_.0001 <- lapply(numeric_list, function(study) {
  lapply(study, filter_by_abundance)
})

##Malawi twins
MT5147_.0001 <- lapply(MT5147_numeric, filter_by_abundance)

##1068
for_filtering_1068 <- read_rds("M1068_Tax.rds")
M1068_tax <- Filter(function(df) ncol(df) > 9, for_filtering_1068)
M1068_.0001 <- lapply(M1068_tax, filter_by_abundance)
M1068_.0001 <- lapply(M1068_.0001, as.data.frame)

newnames1068 <- str_replace(names(M1068_.0001), "SA", "saliva")
newnames1068 <- str_replace(newnames1068, "SP", "supragingival.placque")
newnames1068 <- str_replace(newnames1068, "T", "tongue.dorsum")

newnames1068 <- sapply(newnames1068, function(old_name) {
  number <- sub(".*?(\\d+)$", "\\1", old_name)
  biome <- sub("(.*?)(\\d+)$", "\\1", old_name)
  paste0("S", number, "_", biome)
})

names(M1068_.0001) <- newnames1068

##combining all of them
tax_list_.0001 <- numeric_list_.0001
tax_list_.0001$MT5147 <- MT5147_.0001
tax_list_.0001$M1068 <- M1068_.0001
write_rds(tax_list_.0001, "prevalence_filtered_nested.rds")
##

############################################################
to_presence_absence <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(x > 0))
  return(df)
}
pa_nested_list <- lapply(numeric_list, function(study) {
  lapply(study, to_presence_absence)
})

#Malawi twins
MT5147_pa <- lapply(MT5147_numeric, to_presence_absence)
##

##prevalence filtering
tax_list_.0001_pa <- lapply(tax_list_.0001, function(study) {
  lapply(study, to_presence_absence)
})
##

##Confirming it worked
view(pa_nested_list[["M0278"]][[1]])
####Now we have a list that is in the form of presence/absence
####So we need to remove all of the OTUs that have all 0s in each column of a data frame
remove_all_zero_rows <- function(df) {
  df[rowSums(df) > 0, , drop = FALSE]
}

cleaned_nested_pa_list <- lapply(pa_nested_list, function(study) {
  lapply(study, remove_all_zero_rows)
})

#Malawi twins
MT5147_pa_clean <- lapply(MT5147_pa, remove_all_zero_rows)
##

####Now, I want to convert the column names to the timepoint number and reorder them according to the timepoint
find_bad_timepoints <- function(df) {
  raw_parts <- strsplit(colnames(df), "_")
  timepoints <- sapply(raw_parts, function(x) if (length(x) >= 2) x[2] else NA)
  bad <- timepoints[is.na(suppressWarnings(as.numeric(timepoints)))]
  return(unique(bad))
}
bad_timepoints_by_study <- lapply(cleaned_nested_pa_list, function(study) {
  lapply(study, find_bad_timepoints)
})

#Malawi twins (no bad timepoints here)
lapply(MT5147_pa_clean, find_bad_timepoints)

###I've encountered an issue!!!! there are some "null" columns it seems (I am not super worried)
###The main issue is that for 5147, there were some weirdly duplicated strings in the relevant columns
###So I need to go back and fix it (just a step I appear to have forgotten when making the "usable_tax.csv")
###This issue has been fixed
##I have now gone through all 7 persisting "bad timepoints" and have determined that they can be removed
drop_bad_timepoint_columns <- function(df) {
  # Split column names and extract the timepoint portion
  timepoints <- sapply(strsplit(colnames(df), "_"), function(x) if (length(x) >= 2) x[2] else NA)
  
  # Check which are valid numeric values
  is_good <- !is.na(suppressWarnings(as.numeric(timepoints)))
  
  # Keep only the good columns
  df[, is_good, drop = FALSE]
}
good_nested_pa_list <- lapply(cleaned_nested_pa_list, function(study) {
  lapply(study, drop_bad_timepoint_columns)
})

##prevalence-filtered (no bad timepoints in the last two lists)
good_tax_.0001_pa <- lapply(tax_list_.0001_pa[-c(26:27)], function(study) {
  lapply(study, drop_bad_timepoint_columns)
})
good_tax_.0001_pa$MT5147 <- tax_list_.0001_pa$MT5147

#Adding back 1068 later because it is already ordered
##

convert_colnames_to_timepoints <- function(df) {
  colnames(df) <- as.numeric(sapply(strsplit(colnames(df), "_"), function(x) x[2]))
  df
}
numeric_nested_pa_list <- lapply(good_nested_pa_list, function(study) {
  lapply(study, convert_colnames_to_timepoints)
})

#Malawi twins
numeric_MT5147_pa <- lapply(MT5147_pa_clean, convert_colnames_to_timepoints)

#prevalence filtering
numeric_tax_.0001 <- lapply(good_tax_.0001_pa, function(study) {
  lapply(study, convert_colnames_to_timepoints)
})
numeric_tax_.0001$M1068 <- tax_list_.0001_pa$M1068
##

##removing the dataframes that fell below 10 timepoints with the removal of the bad names
clean_numeric_nested_pa_list <- lapply(numeric_nested_pa_list, function(study) {
  Filter(function(df) ncol(df) > 9, study)
})

##prevalence filtered
clean_tax_.0001 <- lapply(numeric_tax_.0001, function(study) {
  Filter(function(df) ncol(df) > 9, study)
})

#Lastly for this section, ordering the columns by the numbers
sort_columns_by_timepoint <- function(df) {
  df[, order(as.numeric(colnames(df))), drop = FALSE]
}

ordered_nested_pa_list <- lapply(clean_numeric_nested_pa_list, function(study) {
  lapply(study, sort_columns_by_timepoint)
})

#Malawi twins
MT5147_ordered_pa <- lapply(numeric_MT5147_pa ,sort_columns_by_timepoint)

#prevalence filtered
ordered_.0001 <- lapply(clean_tax_.0001, function(study) {
  lapply(study, sort_columns_by_timepoint)
})
##

####Now, let's convert the taxonomy rownames into taxonomy columns for each level of classification
tax_levels <- c("sk__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "st__")
col_names  <- c("superkingdom", "kingdom", "phylum", "class", "order",
                "family", "genus", "species", "strain")
parse_taxonomy <- function(tax_string) {
  # Split and extract level keys and values
  parts <- unlist(strsplit(tax_string, ";"))
  
  # Only keep parts that match expected taxonomic labels
  valid_parts <- parts[grepl("^[a-z]{1,2}__", parts)]
  
  # If no valid taxonomic parts are found, return empty row
  if (length(valid_parts) == 0) {
    return(as.data.frame(as.list(setNames(rep("", length(tax_levels)), col_names)), stringsAsFactors = FALSE))
  }
  
  level_keys <- sub("^([^_]+__).*", "\\1", valid_parts)
  level_values <- sub("^[^_]+__", "", valid_parts)
  
  named_tax <- setNames(level_values, level_keys)
  
  # Fill in any missing levels with blanks
  full_row <- setNames(rep("", length(tax_levels)), tax_levels)
  full_row[names(named_tax)] <- named_tax
  
  as.data.frame(as.list(setNames(full_row, col_names)), stringsAsFactors = FALSE)
}

add_taxonomy_columns <- function(df) {
  rownames_vec <- rownames(df)
  
  # Catch NULL rownames early
  if (is.null(rownames_vec)) {
    stop("Dataframe has NULL rownames — cannot parse taxonomy.")
  }
  
  # Parse taxonomy
  taxonomy_parsed <- lapply(rownames_vec, function(x) {
    result <- parse_taxonomy(x)
    if (ncol(result) != 9) {
      print(paste("Problematic rowname:", x))
      print(result)
    }
    return(result)
  })
  
  taxonomy_df <- do.call(rbind, taxonomy_parsed)
  df <- cbind(taxonomy_df, df)
  rownames(df) <- NULL
  return(df)
}

###Removing some problematic row name additions I made earlier
removek__ <- function(x){
  rownames(x) <- str_remove(rownames(x), "k__;")
  return(x)
}
ordered_nested_pa_list$M1070 <- lapply(ordered_nested_pa_list$M1070, removek__)
ordered_.0001$M1070 <- lapply(ordered_.0001$M1070, removek__)

###Some dataset has "t__" instead of "st__" for strain. Let's find it
find_t_taxonomies <- function(nested_list) {
  lapply(nested_list, function(study) {
    lapply(study, function(df) {
      any(grepl(";t__", rownames(df), fixed = TRUE))
    })
  })
}

t_taxonomy_flags <- find_t_taxonomies(ordered_nested_pa_list)
##The culprit is 6120
replace_t_with_st <- function(x){
  rownames(x) <- str_replace(rownames(x),"t__","st__")
  return(x)
}
ordered_nested_pa_list$M6120 <- lapply(ordered_nested_pa_list$M6120, replace_t_with_st)
ordered_.0001$M6120 <- lapply(ordered_.0001$M6120, replace_t_with_st)
##Fixed

###Now converting taxonomy
tax_nested_list <- lapply(ordered_nested_pa_list, function(study) {
  lapply(study, add_taxonomy_columns)
})

#Malawi twins
MT5147_tax <- lapply(MT5147_ordered_pa, add_taxonomy_columns)

#prevalence filtered
taxonomy_.0001 <- lapply(ordered_.0001, function(study) {
  lapply(study, add_taxonomy_columns)
})
##

##1068 REMOVE ALL OF THE DATAFRAMES IN THIS SPLIT_TAX THAT HAVE FEWER THAN 18 columns
standardize_and_filter_study <- function(study_list, min_cols = 19) {
  # Process each dataframe
  cleaned <- lapply(study_list, function(df) {
    # Rename 'Domain' to 'superkingdom'
    if ("Domain" %in% colnames(df)) {
      colnames(df)[colnames(df) == "Domain"] <- "superkingdom"
    }
    
    # Add 'strain' column if missing
    if (!"strain" %in% colnames(df)) {
      df$strain <- ""
    }
    
    return(df)
  })
  
  # Filter out dataframes with too few columns
  cleaned <- Filter(function(df) ncol(df) >= min_cols, cleaned)
  
  return(cleaned)
}


clean1068 <- standardize_and_filter_study(split_tax1068)
reorder_tax_columns <- function(df_list) {
  lapply(df_list, function(df) {
    if (ncol(df) >= 19) {
      df <- df[, c(1:8, 19, 9:18)]
    }
    return(df)
  })
}
M1068 <- reorder_tax_columns(clean1068)

#next, I have to modify the names
newnames1068 <- str_replace(names(M1068), "SA", "saliva")
newnames1068 <- str_replace(newnames1068, "SP", "supragingival.placque")
newnames1068 <- str_replace(newnames1068, "T", "tongue.dorsum")

newnames1068 <- sapply(newnames1068, function(old_name) {
  number <- sub(".*?(\\d+)$", "\\1", old_name)
  biome <- sub("(.*?)(\\d+)$", "\\1", old_name)
  paste0("S", number, "_", biome)
})

names(M1068) <- newnames1068

##Now I need to add it to my list of studies
tax_nested_list$M1068 <- M1068

####DONE DONE DONE DONE DONE DONE DONE DONE
flatten_nested_list <- function(nested_list) {
  flat_list <- list()
  
  for (study_name in names(nested_list)) {
    study <- nested_list[[study_name]]
    
    for (sample_name in names(study)) {
      new_name <- paste0(study_name, "_", sample_name)
      flat_list[[new_name]] <- study[[sample_name]]
    }
  }
  
  return(flat_list)
}

final_tax_list <- flatten_nested_list(tax_nested_list)

final_.0001 <- flatten_nested_list(taxonomy_.0001)
####SAVE THIS FILE, IT CAN BE USED FOR ALL CALCULATIONS
##Saving both just in case anyway
write_rds(final_tax_list, "final_tax_list.rds")
write_rds(tax_nested_list, "nested_tax_list.rds")
write_rds(final_.0001, "final_flat_.0001.rds")

##Malawi twins
write_rds(MT5147_tax, "final_MT5147.rds")

###DONE





