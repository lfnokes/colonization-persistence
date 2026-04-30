#####This script is about putting together the OTU tables with the metadata necessary
#####to split them up by subject/biome and analyze them by timepoint

#Getting all possibly necessary packages
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

#Getting all of the multi-assay-experiment data (every thing from MGnify)
mae0278 <- read_rds("MGYS00000278/mae_MGYS00000278.rds")#Done
mae0308 <- read_rds("MGYS00000308/mae_MGYS00000308.rds")#Done
#mae1068 <- read_rds("MGYS00001068/mae_MGYS00001068.rds")#Done
mae1278 <- read_rds("mae_MGYS00001278.rds")#Done
mae1809 <- read_rds("mae_MGYS00001809.rds")
mae2184 <- read_rds("mae_MGYS00002184.rds")#Done
mae2425 <- read_rds("MGYS00002425/mae_MGYS00002425.rds")#Done
mae3619 <- read_rds("mae_MGYS00003619.rds")#Done
mae3659 <- read_rds("mae_MGYS00003659.rds")#Done
mae4708 <- read_rds("mae_MGYS00004708.rds")#Done
mae4729 <- read_rds("mae_MGYS00004729.rds")
mae5142_1 <- read_rds("mae_MGYS00005142.rds")
mae5142_2 <- read_rds("mae_MGYS00005142_2.rds")
mae5147 <- read_rds("mae_MGYS00005147.rds")
mae5176 <- read_rds("mae_MGYS00005176.rds")
mae5196 <- read_rds("mae_MGYS00005196.rds")#Done
mae5261 <- read_rds("mae_MGYS00005261.rds")#Done
mae5791 <- read_rds("mae_MGYS00005791.rds")
mae5816 <- read_rds("mae_MGYS00005816.rds")
mae6263 <- read_rds("mae_MGYS00006263.rds")
mae6273 <- read_rds("mae_MGYS00006273.rds")
mae6469 <- read_rds("mae_MGYS00006469.rds")#Done
mae6530 <- read_rds("mae_MGYS00006530.rds")#Done
mae6638 <- read_rds("mae_MGYS00006638.rds")#Done

#Some studies have metadata in ENA:
meta1278 <- read_rds("ENA_metadata_tables/PRJEB6518_sample_metadata.rds")
meta2184 <- read_rds("ENA_metadata_tables/PRJEB19825_sample_metadata.rds")
meta3619 <- read_rds("ENA_metadata_tables/PRJNA290380_sample_metadata.rds")#Different
meta3659 <- read_rds("ENA_metadata_tables/PRJEB6292_sample_metadata.rds")
meta4708 <- read_rds("ENA_metadata_tables/PRJEB32062_sample_metadata.rds")
meta5196 <- read_rds("ENA_metadata_tables/PRJNA393472_sample_metadata.rds")
meta5261 <- read_rds("ENA_metadata_tables/PRJNA423040_sample_metadata.rds")
meta5860 <- read_rds("ENA_metadata_tables/PRJNA301903_sample_metadata.rds")
meta6469 <- read_rds("ENA_metadata_tables/PRJNA360332_sample_metadata.rds")
meta6530 <- read_rds("ENA_metadata_tables/PRJNA520924_sample_metadata.rds")
meta6638 <- read_rds("ENA_metadata_tables/PRJNA790843_sample_attributes.rds")

#Some studies have metadata in their publications:
meta4729 <- read_xlsx("MGYS00004729_meta.xlsx")

####Now, to begin merging the data with the metadata (where metadata is separate)
##1278
mae_meta1278 <- as.data.frame(mae1278@colData@listData)#there are duplicates here, use the ones from pipeline 2.0
mae_tax1278 <- mae1278@ExperimentList$taxonomy@assays@data@listData$counts
mae_meta1278_5.0 <- mae_meta1278 %>% filter(analysis_pipeline.version=="5.0")
merged_metadata1278 <- meta1278 %>%
  inner_join(mae_meta1278_5.0, by = c("sample_accession" = "sample_biosample"))
df1278 <- merged_metadata1278
df1278 <- df1278 %>%
  mutate(
    donor_id = str_extract(`host subject id`, "Donor[A-Z]"),  # Extract "DonorA" or "DonorB"
    donor_id = str_remove(donor_id, "Donor"),  # Remove "Donor", leaving just "A" or "B"
    biome = ifelse(`common_name` == "human gut metagenome", "gut", "oral"),  # Convert biome names
    new_colname = paste(donor_id, collection_day, biome, sep = "_")  # Create the new name
  )
otu_col_mapping1278 <- df %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax1278) <- otu_col_mapping1278[colnames(mae_tax1278), "new_colname"]
write.csv(mae_tax1278, "usable_tax1278.csv")
##now, mae_tax1278 has all of the necessary metadata in the column names (ID_day_biome)

##2184
mae_meta2184 <- as.data.frame(mae2184@colData@listData)
mae_tax2184 <- mae2184@ExperimentList$taxonomy@assays@data@listData$counts
merged_metadata2184 <- meta2184 %>%
  inner_join(mae_meta2184, by = c("sample_accession" = "sample_biosample"))
df2184 <- merged_metadata2184
unique(df2184$organism)
df2184 <- df2184 %>%
  mutate(
    donor_id = host_subject_id,
    biome = ifelse(organism == "human gut metagenome", "gut", ifelse(organism == "human skin metagenome", "skin", "oral")),  # Convert biome names
    new_colname = paste(donor_id, days_since_experiment_start, biome, sep = "_")  # Create the new name
  )
otu_col_mapping2184 <- df2184 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax2184) <- otu_col_mapping2184[colnames(mae_tax2184), "new_colname"]
write.csv(mae_tax2184, "usable_tax2184.csv")
##Done

##3619
mae_meta3619 <- as.data.frame(mae3619@colData@listData)
mae_meta3619_5.0 <- mae_meta3619 %>% filter(analysis_pipeline.version=="5.0")
mae_tax3619 <- mae3619@assays@data@listData$counts #its got some weird organisms (plants?)
mae_tax3619 <- mae_tax3619[-c(1:51),] #removing some obvious contaminants (viridiplantae, aves, etc.)
merged_metadata3619 <- meta3619 %>%
  inner_join(mae_meta3619_5.0, by = c("sample_accession" = "sample_biosample"))
df3619 <- merged_metadata3619
df3619 <- df3619 %>%
  mutate(
    donor_id = host_subject_id,
    biome = ifelse(organism == "human gut metagenome", "gut", "gut"),  # Convert biome names
    new_colname = paste(donor_id, host_age, biome, sep = "_")  # Create the new name
  )
otu_col_mapping3619 <- df3619 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax3619) <- otu_col_mapping3619[colnames(mae_tax3619), "new_colname"]
mae_tax3619 <- mae_tax3619[,-c(1:759)] #removing all of the ones from the previous pipeline
write.csv(mae_tax3619, "usable_tax3619.csv")
##Done

##3659
mae_meta3659 <- as.data.frame(mae3659@colData@listData)
mae_tax3659 <- mae3659@ExperimentList$taxonomy@assays@data@listData$counts
merged_metadata3659 <- meta3659 %>%
  inner_join(mae_meta3659, by = c("sample_accession" = "sample_biosample"))
df3659 <- merged_metadata3659
df3659 <- df3659 %>%
  mutate(
    donor_id = str_replace(subject,"_","-"),
    biome = str_replace(surface,"_","-"),  # Convert biome names
    daynumber = str_remove(day, "D"),
    new_colname = paste(donor_id, daynumber, biome, sep = "_")  # Create the new name
  )
otu_col_mapping3659 <- df3659 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax3659) <- otu_col_mapping3659[colnames(mae_tax3659), "new_colname"]
write.csv(mae_tax3659, "usable_tax3659.csv")
##Done

##4708
mae_meta4708 <- as.data.frame(mae4708@colData@listData)
mae_tax4708 <- mae4708@ExperimentList$taxonomy@assays@data@listData$counts
merged_metadata4708 <- meta4708 %>%
  inner_join(mae_meta4708, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata4708)
merged_metadata4708$numeric_date <- as.numeric(as.Date(unlist(merged_metadata4708$`collection date`)))
df4708 <- merged_metadata4708
df4708 <- df4708 %>%
  mutate(
    donor_id = `host subject id`,
    biome = `environment (biome)`,  # Convert biome names
    new_colname = paste(donor_id, numeric_date, biome, sep = "_")  # Create the new name
  )
otu_col_mapping4708 <- df4708 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax4708) <- otu_col_mapping4708[colnames(mae_tax4708), "new_colname"]
write.csv(mae_tax4708, "usable_tax4708.csv")
##Done

##5196
mae_meta5196 <- as.data.frame(mae5196@colData@listData)
mae_tax5196 <- mae5196@assays@data@listData$counts
merged_metadata5196 <- meta5196 %>%
  inner_join(mae_meta5196, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5196)
df5196 <- merged_metadata5196
df5196 <- df5196 %>%
  mutate(
    donor_id = host_subject_id,
    biome = ifelse(host_tissue_sampled == "vaginal epithelial surface", "vaginal", "reagent"),  # Convert biome names
    new_colname = paste(donor_id, gest_day_collection, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5196 <- df5196 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5196) <- otu_col_mapping5196[colnames(mae_tax5196), "new_colname"]
write.csv(mae_tax5196, "usable_tax5196.csv")
##Done

##5261
mae_meta5261 <- as.data.frame(mae5261@colData@listData)
mae_tax5261 <- mae5261@assays@data@listData$counts
mae_meta5261_5.0 <- mae_meta5261 %>% filter(analysis_pipeline.version=="5.0")
merged_metadata5261 <- meta5261 %>%
  inner_join(mae_meta5261_5.0, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5261)
df5261 <- merged_metadata5261
df5261 <- df5261 %>%
  mutate(
    donor_id = isolate,
    biome = ifelse(organism == "human lung metagenome", "lung", "lung"),  # Convert biome names
    new_colname = paste(donor_id, host_age, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5261 <- df5261 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5261) <- otu_col_mapping5261[colnames(mae_tax5261), "new_colname"]
mae_tax5261 <- mae_tax5261[,-c(1:478)]
write.csv(mae_tax5261, "usable_tax5261.csv")
##Done

##5860/2425
mae_meta5860 <- as.data.frame(mae2425@colData@listData)
mae_tax5860 <- mae2425@assays@data@listData$counts
merged_metadata5860 <- meta5860 %>%
  inner_join(mae_meta5860, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5860)
df5860 <- merged_metadata5860
df5860$host_subject_id_clean <- sapply(df5860$host_subject_id, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
df5860$host_day_of_life_clean <- sapply(df5860$host_day_of_life, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
df5860 <- df5860 %>%
  mutate(
    donor_id = df5860$host_subject_id_clean,
    biome = "gut",  # Convert biome names
    day = df5860$host_day_of_life_clean,
    new_colname = paste(donor_id, day, biome, sep = "_")  # Create the new name
  )
df5860$new_colname
otu_col_mapping5860 <- df5860 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5860) <- otu_col_mapping5860[colnames(mae_tax5860), "new_colname"]
write.csv(mae_tax5860, "usable_tax5860.csv")
##Done

##6469
mae_meta6469 <- as.data.frame(mae6469@colData@listData)
mae_tax6469 <- mae6469@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata6469 <- meta6469 %>%
  inner_join(mae_meta6469, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata6469)
merged_metadata6469$numeric_date <- as.numeric(as.Date(unlist(merged_metadata6469$collection_date)))
df6469 <- merged_metadata6469
df6469 <- df6469 %>%
  mutate(
    donor_id = str_remove(isolation_source, "participant_"),
    biome = "lung",  # Convert biome names
    new_colname = paste(donor_id, numeric_date, biome, sep = "_")  # Create the new name
  )
otu_col_mapping6469 <- df6469 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax6469) <- otu_col_mapping6469[colnames(mae_tax6469), "new_colname"]
write.csv(mae_tax6469, "usable_tax6469.csv")
##Done

##6530
mae_meta6530 <- as.data.frame(mae6530@colData@listData)
mae_tax6530 <- mae6530@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata6530 <- meta6530 %>%
  inner_join(mae_meta6530, by = c("sample_accession" = "sample_biosample"))
unique(merged_metadata6530$organism)
view(merged_metadata6530)
merged_metadata6530$Subject_ID_Complete <- coalesce(merged_metadata6530$`Subject ID`, 
                                                    merged_metadata6530$Subject.ID)
merged_metadata6530$Day_Complete <- coalesce(merged_metadata6530$Day,
                                             merged_metadata6530$Days)
view(merged_metadata6530$sample_accession)
df6530 <- merged_metadata6530
df6530 <- df6530 %>%
  mutate(
    donor_id = Subject_ID_Complete,
    biome = ifelse(organism == "human lung metagenome", "lung", "reagent"),  # Convert biome names
    new_colname = paste(donor_id, Day_Complete, biome, sep = "_")  # Create the new name
  )
otu_col_mapping6530 <- df6530 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax6530) <- otu_col_mapping6530[colnames(mae_tax6530), "new_colname"]
write.csv(mae_tax6530, "usable_tax6530.csv")
##Done

##6638
mae_meta6638 <- as.data.frame(mae6638@colData@listData)
mae_tax6638 <- mae6638@ExperimentList$taxonomy@assays@data@listData$counts
merged_metadata6638 <- meta6638 %>%
  inner_join(mae_meta6638, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata6638)
df6638 <- merged_metadata6638
df6638 <- df6638 %>%
  mutate(
    donor_id = host_subject_id,
    biome = "nasopharyngeal",  # Convert biome names
    new_colname = paste(donor_id, host_age_days, biome, sep = "_")  # Create the new name
  )
otu_col_mapping6638 <- df6638 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax6638) <- otu_col_mapping6638[colnames(mae_tax6638), "new_colname"]
write.csv(mae_tax6638, "usable_tax6638.csv")
##Done

####Now, to begin merging data with metadata (where metadata is included)
##0278
#mae_meta0278 <- as.data.frame(mae0278@colData@listData)
mae_tax0278 <- mae0278@ExperimentList$taxonomy@assays@data@listData$counts
view(mae_meta0278)
analysis_accession0278 <- c("MGYA00000223", "MGYA00000224", "MGYA00000225", "MGYA00000226", 
                            "MGYA00000227", "MGYA00000228", "MGYA00000229", "MGYA00000230", 
                            "MGYA00000231", "MGYA00000232", "MGYA00000233", "MGYA00000234",
                            "MGYA00000235", "MGYA00000236", "MGYA00000237", "MGYA00000238",
                            "MGYA00000239")
day0278 <- as.character(c(371, 118, 413, 432, 441, 454, 6, 85, 92, 98, 3, 100, 623, 745, 831, 835, 838))
new_colname0278 <- paste("Infant1", day, "gut", sep = "_")
df0278 <- data.frame(analysis_accession = analysis_accession0278, new_colname = new_colname0278)
otu_col_mapping0278 <- df0278 %>% 
  select("analysis_accession", "new_colname") %>% 
  column_to_rownames("analysis_accession")
colnames(mae_tax0278) <- otu_col_mapping0278[colnames(mae_tax0278), "new_colname"]
view(mae_tax0278)
write.csv(mae_tax0278, "usable_tax0278.csv")
##Done

##0308
mae_tax0308 <- mae0308@ExperimentList$taxonomy@assays@data@listData$counts
meta0308 <- read.csv("meta0308.csv")
df0308 <- meta0308
df0308 <- df0308 %>%
  mutate(
    donor_id = host_subject_id,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, sampling_day, biome, sep = "_")  # Create the new name
  )
otu_col_mapping0308 <- df0308 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax0308) <- otu_col_mapping0308[colnames(mae_tax0308), "new_colname"]
write.csv(mae_tax0308, "usable_tax0308.csv")
##Done (remember for this one and 0278 the taxonomy tables are non-numeric when there is a 0, it is NA instead)

##1068
mae_tax1068 <- read_rds("MGYS00001068/filtered_PresTax_MGYS00001068.rds")
write_rds(mae_tax1068, "split_tax1068.rds")
##Done. This has already been split up by subject and biome (biome is the two letter code, subject is the number)

###I am adding all of the ENA metadata because it makes this process a lot easier and more standardized
meta1809 <- read_rds("ENA_metadata_tables/PRJEB21196_sample_metadata.rds")
meta5147 <- read_rds("ENA_metadata_tables/PRJEB5482_sample_metadata.rds")
meta5176 <- read_rds("ENA_metadata_tables/PRJNA390646_sample_metadata.rds")
meta5791 <- read_rds("ENA_metadata_tables/PRJNA385949_sample_metadata.rds")
meta5816 <- read_rds("ENA_metadata_tables/PRJEB8347_sample_metadata.rds")
meta6273 <- read_rds("ENA_metadata_tables/PRJEB36034_sample_metadata.rds")

##1809 EXCLUDED DUE TO NOT ENOUGH LONGITUDINAL SAMPLES

##5147/PRJEB5482
mae_meta5147 <- as.data.frame(mae5147@colData@listData)
mae_tax5147 <- mae5147@assays@data@listData$counts
merged_metadata5147 <- meta5147 %>%
  inner_join(mae_meta5147, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5147)
unique(merged_metadata5147$`host subject id`)
unique(merged_metadata5147$organism)
df5147 <- merged_metadata5147
df5147$host_subject_id <- sapply(df5147$`host subject id`, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
df5147$age_clean <- sapply(df5147$age, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
which(df5147$sample_sample.alias)
df5147 <- df5147 %>%
  mutate(
    donor_id = host_subject_id,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, age_clean, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5147 <- df5147 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5147) <- otu_col_mapping5147[colnames(mae_tax5147), "new_colname"]
view(mae_tax5147)
write.csv(mae_tax5147, "usable_tax5147.csv")
##Done

##5147.2
#There is an issue with the names in "host_subject_id" that merges all of the 
#Malawian twin pairs. I am making a separate dataset for them so I can add it for the analysis
malawi5147 <- df5147 %>%
  filter(str_detect(sample_sample.alias, "Bgtw"))
unique(malawi5147$sample_sample.alias)
malawi5147$twin_id <- str_replace(malawi5147$sample_sample.alias, "^([^.]+\\.[^.]+)\\..*$", "\\1")
malawi5147 <- malawi5147 %>%
  mutate(
    donor_id = twin_id,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, age_clean, biome, sep = "_")  # Create the new name
  )
malawi_col_mapping5147 <- malawi5147 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
malawi_tax5147 <- mae5147@assays@data@listData$counts
colnames(malawi_tax5147) <- malawi_col_mapping5147[colnames(malawi_tax5147), "new_colname"]
view(malawi_tax5147)
write.csv(malawi_tax5147, "usable_malawi5147.csv")

##5176/PRJNA390646
mae_meta5176 <- as.data.frame(mae5176@colData@listData)
mae_tax5176 <- mae5176@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata5176 <- meta5176 %>%
  inner_join(mae_meta5176, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5176)
df5176 <- merged_metadata5176
df5176$host_id <- str_extract(df5176$sample_sample.alias, ".*ID")
df5176$host_id <- str_replace(df5176$host_id, "ID", "")
df5176$sample_day_clean <- sapply(df5176$Target_Collection_Day, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
df5176 <- df5176 %>%
  mutate(
    donor_id = host_id,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, sample_day_clean, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5176 <- df5176 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5176) <- otu_col_mapping5176[colnames(mae_tax5176), "new_colname"]
view(mae_tax5176)
write.csv(mae_tax5176, "usable_tax5176.csv")
##Done

##5791/PRJNA385949 #Different ENA number
mae_meta5791 <- as.data.frame(mae5791@colData@listData)
mae_tax5791 <- mae5791@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata5791 <- meta5791 %>%
  inner_join(mae_meta5791, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5791)
df5791 <- merged_metadata5791
df5791 <- df5791 %>%
  mutate(
    donor_id = subject,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, month, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5791 <- df5791 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5791) <- otu_col_mapping5791[colnames(mae_tax5791), "new_colname"]
view(mae_tax5791)
write.csv(mae_tax5791, "usable_tax5791.csv")
##Done

##5816/PRJEB8347 #Different ENA number
mae_meta5816 <- as.data.frame(mae5816@colData@listData)
mae_tax5816 <- mae5816@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata5816 <- meta5816 %>%
  inner_join(mae_meta5816, by = c("sample_accession" = "sample_biosample"))
view(merged_metadata5816)
unique(merged_metadata5816$organism)
df5816 <- merged_metadata5816
df5816$subject_id_clean <- sapply(df5816$`host subject id`, function(x) {
  paste(unique(strsplit(x, ",")[[1]]), collapse = ",")
})
df5816$day <- sub(".*?-.*?-(\\d+)-.*", "\\1", df5816$sample_sample.alias)
df5816 <- df5816 %>%
  mutate(
    donor_id = subject_id_clean,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, day, biome, sep = "_")  # Create the new name
  )
otu_col_mapping5816 <- df5816 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5816) <- otu_col_mapping5816[colnames(mae_tax5816), "new_colname"]
view(mae_tax5816)
write.csv(mae_tax5816, "usable_tax5816.csv")
##Done

##6273/PRJEB36034
mae_meta6273 <- as.data.frame(mae6273@colData@listData)
mae_tax6273 <- mae6273@ExperimentList@listData$taxonomy@assays@data@listData$counts
merged_metadata6273 <- meta6273 %>%
  inner_join(mae_meta6273, by = c("sample_accession" = "sample_biosample"))
reindeer_meta6273 <- merged_metadata6273[which(merged_metadata6273$host_subject_id == "reindeer"),]
view(merged_metadata6273)
df6273 <- merged_metadata6273
df6273 <- df6273 %>%
  mutate(
    donor_id = host_subject_id,
    biome = ifelse(organism == "human lung metagenome", "lung", "reagent"),  # Convert biome names
    new_colname = paste(donor_id, days_since_origin, biome, sep = "_")  # Create the new name
  )
otu_col_mapping6273 <- df6273 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax6273) <- otu_col_mapping6273[colnames(mae_tax6273), "new_colname"]
view(mae_tax6273)
write.csv(mae_tax6273, "usable_tax6273.csv")
##Done

####Now, I am doing the ones without ANY good metadata in ENA, only OK metadata in MGnify

##6263
mae_meta6263 <- as.data.frame(mae6263@colData@listData)
mae_tax6263 <- mae6263@ExperimentList@listData$taxonomy@assays@data@listData$counts
view(mae_meta6263)
df6263 <- mae_meta6263
split_parts6263 <- do.call(rbind, strsplit(df6263$sample_sample.desc, "\\.", fixed = FALSE))
df6263$sample_id <- split_parts6263[,1]
df6263$sample_day <- split_parts6263[,2]
df6263 <- df6263 %>%
  mutate(
    donor_id = sample_id,
    biome = "oropharyngeal",  # Convert biome names
    new_colname = paste(donor_id, sample_day, biome, sep = "_")  # Create the new name
  )
otu_col_mapping6263 <- df6263 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax6263) <- otu_col_mapping6263[colnames(mae_tax6263), "new_colname"]
view(mae_tax6263)
write.csv(mae_tax6263, "usable_tax6263.csv")
##Done

##5142 (This one might require the publication metadata table - the Sample Alias should work though)
mae_meta5142_1 <- as.data.frame(mae5142_1@colData@listData)
mae_meta5142_2 <- as.data.frame(mae5142_2@colData@listData)
mae_meta5142 <- rbind(mae_meta5142_1, mae_meta5142_2)
view(mae_meta5142)
mae_tax5142_1 <- mae5142_1@assays@data@listData$counts
mae_tax5142_2 <- mae5142_2@assays@data@listData$counts
intersection <- intersect(rownames(mae_tax5142_1),rownames(mae_tax5142_2))
union <- union(rownames(mae_tax5142_1), rownames(mae_tax5142_2))
#Merging the taxonomy tables
otu1 <- mae_tax5142_1
otu2 <- mae_tax5142_2
otu1_full <- matrix(0,
                    nrow = length(union),
                    ncol = ncol(otu1),
                    dimnames = list(union, colnames(otu1)))

otu2_full <- matrix(0,
                    nrow = length(union),
                    ncol = ncol(otu2),
                    dimnames = list(union, colnames(otu2)))
otu1_full[rownames(otu1), ] <- as.matrix(otu1)
otu2_full[rownames(otu2), ] <- as.matrix(otu2)
combined_otu <- cbind(otu1_full, otu2_full)
mae_tax5142 <- combined_otu
view(mae_meta5142)
#the samples with JP02 are mice samples, so I am removing them
str
remove <- which(str_detect(mae_meta5142$sample_sample.alias, "JP02")==TRUE)
mae_meta5142_filtered <- mae_meta5142[-remove,]
df5142 <- mae_meta5142_filtered
mae_meta5142_filtered$sample_sample.alias
#extracting the host ids from the 1st-3rd part of the sample alias
df5142$sample_id <- str_extract(mae_meta5142_filtered$sample_sample.alias, "^([^\\.]+\\.){2}[^\\.]+")
#extracting the sample timepoint from the 4th part of the sample alias
df5142$sample_day <- str_match(mae_meta5142_filtered$sample_sample.alias, "^([^\\.]+\\.){3}([^\\.]+)")[, 3]
df5142$sample_day <- str_remove(df5142$sample_day, "d")
#Extracting fraction info from the last part of the sample alias
df5142$biome <- str_extract(mae_meta5142_filtered$sample_sample.alias, "[^\\.]+$")
#these biomes are all variations on "fecal" but refer to slightly different fractions/extraction methods
df5142 <- df5142 %>%
  mutate(
    donor_id = sample_id,
    biome = biome,  # Convert biome names
    new_colname = paste(donor_id, sample_day, biome, sep = "_")  # Create the new name
  )
rownames(df5142) <- NULL
otu_col_mapping5142 <- df5142 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax5142) <- otu_col_mapping5142[colnames(mae_tax5142), "new_colname"]
view(mae_tax5142)
write.csv(mae_tax5142, "usable_tax5142.csv")
##Done

####Now, for the ones where the only good metadata is in the publication

##4729
mae_meta4729 <- as.data.frame(mae4729@colData@listData)
mae_tax4729 <- mae4729@ExperimentList@listData$taxonomy@assays@data@listData$counts
meta4729_filtered <- meta4729[-c(1:10),]
colnames(meta4729_filtered) <- meta4729[9,]
meta4729_filtered <- as.data.frame(meta4729_filtered[,-c(9:22)])
meta4729_filtered$`Accession no. ¥`
###need to remove the second accession number in the strings separated by a comma in the above column
meta4729_filtered$accession_no <- sub(",.*", "", meta4729_filtered$`Accession no. ¥`)
meta4729_filtered <- meta4729_filtered[,-11]
merged_metadata4729 <- meta4729_filtered %>%
  inner_join(mae_meta4729, by = c("accession_no" = "sample_accession"))
view(merged_metadata4729)
df4729 <- merged_metadata4729
df4729 <- df4729 %>%
  mutate(
    donor_id = Infant,
    biome = "nasopharyngeal",  # Convert biome names
    new_colname = paste(donor_id, `Date collected`, biome, sep = "_")  # Create the new name
  )
otu_col_mapping4729 <- df4729 %>%
  select(analysis_accession, new_colname) %>%
  column_to_rownames("analysis_accession")
colnames(mae_tax4729) <- otu_col_mapping4729[colnames(mae_tax4729), "new_colname"]
#This study has duplicates - unclear what the differences are between the duplicates but its every other analysis is a duplicate
write.csv(mae_tax4729, "usable_tax4729.csv")
##Done

####Now, to add studies where all data in their publications

biom6120 <- read_biom("MGYS00006120/taxonomic_profiles.biom")

##1070 (these taxa lack the full lineage so I will have to convert)
#install.packages("taxize")
library(taxize)
tax_incomplete1070 <- read_xlsx("tax_incomplete1070.xlsx")
transposed_incomplete1070 <- t(tax_incomplete1070)
tax_list <- as.list(transposed_incomplete1070[,3])
tax_list <- tax_list[-c(1:2)]
tax_formatted <- tax_incomplete1070[-c(1:3),-c(1:2)]
colnames(tax_formatted) <- tax_list
rownames(tax_formatted) <- tax_incomplete1070$`Additional File 3. Relative proportions of taxa observed in all samples.`[-c(1:3)]
taxa_names1070 <- colnames(tax_formatted)
taxa_chunks <- split(taxa_names1070, ceiling(seq_along(taxa_names1070) / 20))  # 20 at a time

all_results <- list()

for (i in seq_along(taxa_chunks)) {
  cat("Querying chunk", i, "\n")
  chunk <- taxa_chunks[[i]]
  
  res <- classification(chunk, db = "ncbi", messages = FALSE)
  all_results <- c(all_results, res)
  
  Sys.sleep(2)  # Pause for 2 seconds to avoid NCBI rate limits
}

format_lineage <- function(tax_df) {
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  prefix <- c("sk__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  lineage <- sapply(ranks, function(r) {
    entry <- tax_df$name[tax_df$rank == r]
    if (length(entry) > 0) paste0(prefix[which(ranks == r)], entry) else NA
  })
  
  paste(na.omit(lineage), collapse = ";")
}

lineage_strings <- sapply(all_results, function(x) {
  if (inherits(x, "data.frame")) {
    format_lineage(x)
  } else {
    NA  # if the entry is not a dataframe (e.g. failed lookup)
  }
})

make_unique_with_suffix <- function(strings, suffixes = LETTERS) {
  # Create a copy
  new_strings <- strings
  counts <- table(strings)
  
  # Find duplicated names
  dup_names <- names(counts[counts > 1])
  
  for (name in dup_names) {
    idx <- which(strings == name)
    suffix <- suffixes[seq_along(idx)]
    new_strings[idx] <- paste0(name, "_", suffix)
  }
  
  return(new_strings)
}
tax_strings1070 <- make_unique_with_suffix(lineage_strings)
old_colnames1070 <- make_unique_with_suffix(colnames(tax_formatted))

names(tax_strings1070) <- old_colnames1070
tax_withnames1070 <- tax_formatted
colnames(tax_withnames1070) <- old_colnames1070
colnames(tax_withnames1070) <- tax_strings1070[colnames(tax_withnames1070)]
final_tax1070 <- t(tax_withnames1070)
view(final_tax1070)
#NOW I need to go through and edit the column names so they can be used in the same way
# Extract current column names
col_names <- colnames(final_tax1070)

# Use regular expression to extract components
parsed <- stringr::str_match(col_names, "^s(\\d+)\\.w(\\d+)d(\\d+)$")
# Columns:      [,1]       [,2]   [,3]   [,4]
# Example: "s1.w3d2" → "s1.w3d2", "1",   "3",  "2"

# Convert week and day to numeric
subject <- paste0("s", parsed[, 2])
week <- as.numeric(parsed[, 3])
day <- as.numeric(parsed[, 4])

# Calculate study day
study_day <- 7 * (week - 1) + day

# Format new names
new_col_names <- paste0(subject, "_", study_day, "_vaginal")

# Assign new names back to the OTU table
colnames(final_tax1070) <- new_col_names
write.csv(final_tax1070, "usable_tax1070.csv")
##Done

##1818
tax1818 <- read_tsv("tax_1818.txt", skip = 1)
colnames(tax1818) <- str_replace(colnames(tax1818), "A", "Armpit")
colnames(tax1818) <- str_replace(colnames(tax1818), "B", "Back")
colnames(tax1818) <- str_replace(colnames(tax1818), "C", "Chest")
colnames(tax1818) <- str_replace(colnames(tax1818), "G", "Groin")
colnames(tax1818) <- str_replace(colnames(tax1818), "L", "Leg")
row_names <- tax1818$`#OTU ID` 
tax1818 <- tax1818[,-1]
rownames(tax1818) <- row_names
#Renaming the columns in correct order
col_names1818 <- colnames(tax1818)

split_parts1818 <- strsplit(col_names1818, "_")

new_names1818 <- sapply(split_parts1818, function(x) {
  paste0("Astronaut_", x[2], "_", x[1])
})

colnames(tax1818) <- new_names1818
view(tax1818)
write.csv(tax1818, "usable_tax1818.csv")
##Done

##5378
otu_table5378 <- read.csv("MGYS00005378/otu_table5378.csv", row.names = 1, check.names = FALSE)
taxonomy5378 <- read.csv("MGYS00005378/SILVA_taxonomy5378.csv", row.names = 1, check.names = FALSE)
taxonomy5378
taxonomy5378$full_tax <- apply(taxonomy5378, 1, function(x) paste(na.omit(x), collapse = ";"))
taxonomy5378$full_tax
#mapping taxonomy to OTUs
otu_to_tax5378 <- taxonomy5378$full_tax
names(otu_to_tax5378) <- taxonomy5378$otu
view(otu_to_tax5378)
otu_colnames5378 <- colnames(otu_table5378)
mapped_names5378 <- otu_to_tax5378[otu_colnames5378]
#converting OTUs to taxonomy in the table
colnames(otu_table5378) <- mapped_names5378
#correcting rownames
rownames(otu_table5378) <- otu_table5378[,1]
ttax5378 <- otu_table5378[,-c(1:3)]
view(tax5378)
tax5378 <- t(ttax5378)
columns5378 <- paste0(colnames(tax5378), "_gut")
colnames(tax5378) <- columns5378
view(tax5378)
write.csv(tax5378, "usable_tax5378.csv")
##Done

##6120
tax6120 <- as.data.frame(as.matrix(biom_data(biom6120)))
meta6120 <- read.csv("MGYS00006120/metadata_6120.csv")
df6120 <- meta6120
df6120 <- df6120[which(df6120$data_type=="metagenomics"),]
df6120 <- df6120 %>%
  mutate(
    donor_id = Participant.ID,
    biome = "gut",  # Convert biome names
    new_colname = paste(donor_id, interval_days, biome, sep = "_")  # Create the new name
  )
rownames(df6120) <- NULL
otu_col_mapping6120 <- df6120 %>%
  select(External.ID, new_colname) %>%
  column_to_rownames("External.ID")
colnames(tax6120) <- otu_col_mapping6120[colnames(tax6120), "new_colname"]
view(tax6120)
write.csv(tax6120, "usable_tax6120.csv")
##Done

####The only metadata associated with all of these taxonomy tables is the subject, timepoint, and biome. 
####Further metadata is too much to add into the columns now and must be added later subject by subject