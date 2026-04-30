setwd("Studies_scripts//ENA_metadata_tables")
library(httr)
library(readr)
library(dplyr)
library(xml2)
library(rvest)
library(tidyr)

#PRJNA790843

url_samples <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA790843&result=read_run&fields=sample_accession&format=tsv&download=true"
response <- GET(url_samples)
samples_df <- read_tsv(httr::content(response, "text"), col_types = cols())

fetch_sample_metadata_xml <- function(sample_accession) {
  url_metadata <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", sample_accession)
  
  response <- GET(url_metadata)
  
  if (status_code(response) == 200) {
    xml_content <- httr::content(response, "text")  # Retrieve XML as text
    return(xml_content)
  } else {
    print(paste("Failed to fetch metadata for", sample_accession))
    return(NULL)
  }
}

parse_sample_attributes <- function(xml_data) {
  xml_doc <- read_xml(xml_data)
  
  # Extract all <SAMPLE_ATTRIBUTE> nodes
  attributes <- xml_doc %>% xml_find_all(".//SAMPLE_ATTRIBUTE")
  
  # Convert XML attributes into a data frame
  df <- tibble(
    TAG = attributes %>% xml_find_first("TAG") %>% xml_text(),
    VALUE = attributes %>% xml_find_first("VALUE") %>% xml_text()
  )
  
  return(df)
}


# Fetch and parse metadata for all samples
all_sample_attributes <- lapply(samples_df$sample_accession, function(sample) {
  xml_data <- fetch_sample_metadata_xml(sample)
  if (!is.null(xml_data)) {
    attributes <- parse_sample_attributes(xml_data)
    attributes$sample_accession <- sample  # Add sample ID for reference
    return(attributes)
  } else {
    return(NULL)
  }
})

all_sample_attributes_clean <- all_sample_attributes[!sapply(all_sample_attributes, is.null)]
# Combine all samples into one dataframe

all_sample_attributes_clean <- lapply(all_sample_attributes_clean, as.data.frame)

#library(data.table)
#metadata_df <- rbindlist(all_sample_attributes_clean, fill = TRUE)

metadata_df <- do.call(rbind, all_sample_attributes_clean)

metadata_df_wide <- metadata_df %>%
  pivot_wider(names_from = TAG, values_from = VALUE)



# Save to file
write_rds(metadata_df_wide, "ENA_metadata_tables/PRJNA790843_sample_attributes.rds")

library(httr)
library(readr)
library(dplyr)
library(xml2)
library(rvest)
library(tidyr)
library(data.table)

# Function to fetch and process metadata for a given project
fetch_project_metadata <- function(project_accession) {
  # Step 1: Retrieve Sample Accessions for the Project
  url_samples <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                        project_accession, "&result=read_run&fields=sample_accession&format=tsv&download=true")
  
  response <- GET(url_samples)
  
  if (status_code(response) != 200) {
    print(paste("Failed to fetch sample accessions for project:", project_accession))
    return(NULL)
  }
  
  samples_df <- read_tsv(httr::content(response,"text", encoding = "UTF-8"), col_types = cols())
  
  if (nrow(samples_df) == 0) {
    print(paste("No samples found for project:", project_accession))
    return(NULL)
  }
  
  # Step 2: Fetch XML Metadata for Each Sample
  fetch_sample_metadata_xml <- function(sample_accession) {
    url_metadata <- paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", sample_accession)
    response <- GET(url_metadata)
    
    if (status_code(response) == 200) {
      return(httr::content(response, "text", encoding = "UTF-8"))  # Retrieve XML as text
    } else {
      print(paste("Failed to fetch metadata for sample:", sample_accession))
      return(NULL)
    }
  }
  
  # Step 3: Parse Sample Attributes from XML
  parse_sample_attributes <- function(xml_data) {
    xml_doc <- read_xml(xml_data)
    
    # Extract <SAMPLE_ATTRIBUTE> nodes
    attributes <- xml_doc %>% xml_find_all(".//SAMPLE_ATTRIBUTE")
    
    if (length(attributes) == 0) {
      return(NULL)  # No attributes found
    }
    
    # Convert XML attributes into a data frame
    df <- tibble(
      TAG = xml_find_first(attributes, "TAG") %>% xml_text(trim = TRUE),
      VALUE = xml_find_first(attributes, "VALUE") %>% xml_text(trim = TRUE)
    )
    
    return(df)
  }
  
  # Step 4: Process Each Sample
  all_sample_attributes <- lapply(samples_df$sample_accession, function(sample) {
    xml_data <- fetch_sample_metadata_xml(sample)
    if (!is.null(xml_data)) {
      attributes <- parse_sample_attributes(xml_data)
      if (!is.null(attributes)) {
        attributes$sample_accession <- sample  # Add sample ID for reference
        return(attributes)
      }
    }
    return(NULL)
  })
  
  # Step 5: Clean and Merge Data
  all_sample_attributes_clean <- all_sample_attributes[!sapply(all_sample_attributes, is.null)]
  
  if (length(all_sample_attributes_clean) == 0) {
    print(paste("No metadata extracted for project:", project_accession))
    return(NULL)
  }
  
  metadata_df <- do.call(rbind, all_sample_attributes_clean)
  
  # Step 6: Reshape Data (Samples as Rows, Tags as Columns)
  metadata_df_wide <- metadata_df %>%
    pivot_wider(names_from = TAG, values_from = VALUE)
  
  # Step 7: Save to File
  output_file <- paste0(project_accession, "_sample_metadata.rds")
  write_rds(metadata_df_wide, output_file)
  
  print(paste("Metadata saved to:", output_file))
  
  return(metadata_df_wide)
}

fetch_project_metadata("PRJNA290380")

#https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB26925&result=read_run&fields=sample_accession&format=tsv&download=true

accession_numbers <- c("PRJNA290380","PRJEB6518","PRJNA520924","PRJNA360332",
                       "PRJNA301903","PRJNA526551","PRJNA196801","PRJNA393472",
                       "PRJEB32062","PRJEB19825","PRJEB13117","PRJEB6292", 
                       "PRJNA423040")

lapply(accession_numbers, fetch_project_metadata)

accession_numbers_2 <- c("PRJEB21196", "PRJEB5482", "PRJNA390646", "PRJEB43973", 
                         "PRJEB41655", "PRJEB36034")

lapply(accession_numbers_2, fetch_project_metadata)

##Didn't work for PRJEB43973 or PRJEB41655 because they correspond to the assemblies but not the real samples
fetch_project_metadata("PRJNA385949")

fetch_project_metadata("PRJEB8347")
