#### In this script I will be running the API request to find relevant projects
## Step 1: Getting packages to get the data from the MGnify API
library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)
#remotes::install_github("EBI-Metagenomics/MGnifyR")
library(MGnifyR)


##Step 2: compiling a master dataset with all of the relevant studies and metadata
# Base URL for the API query
base_url <- "https://www.ebi.ac.uk/metagenomics/api/latest/studies"

# Initialize an empty list to store all results
all_studies <- list()

# Pagination parameters
page <- 1
page_size <- 25  # Default page size
more_pages <- TRUE

while (more_pages) {
  # Make the API request with pagination
  res <- GET(
    url = base_url,
    query = list(
      accession = "", 
      biome_name = "", 
      lineage = "root:Host-associated:Human",
      page = page
    )
  )
  
  # Parse the response
  data <- fromJSON(rawToChar(res$content))
  
  # Append current page data to the results list
  all_studies <- append(all_studies, data$data)
  
  # Check if there's a next page
  if (!is.null(data$links$`next`)) {
    page <- page + 1  # Increment to the next page
  } else {
    more_pages <- FALSE  # No more pages
  }
}

#combining humans studies into master df
master_df <- all_studies[[3]]
for (i in 1:35){
  master_df <- rbind(master_df, all_studies[[3+5*i]])
}

#getting data for built environment studies
built_res <- GET(
  url = base_url,
  query = list(
    accession = "", 
    biome_name = "", 
    lineage = "root:Engineered:Built environment"
  )
)

built_data <- fromJSON(rawToChar(built_res$content))
built_atts <- built_data$data$attributes

master_df <- rbind(master_df, built_atts)

##Discovered the mixed biome option on MGnify, adding those datasets here (2/10/25)
page <- 1
page_size <- 25
more_pages <- TRUE
mixed_studies <- list()

while(more_pages){
  mixed_res <- GET(
    url = base_url,
    query = list(
      accession = "", 
      biome_name = "", 
      lineage = "root:Mixed",
      page = page
    )
  )
  
  # Parse the response
  mixed_data <- fromJSON(rawToChar(mixed_res$content))
  
  # Append current page data to the results list
  mixed_studies <- append(mixed_studies, mixed_data$data)
  
  # Check if there's a next page
  if (!is.null(mixed_data$links$`next`)) {
    page <- page + 1  # Increment to the next page
  } else {
    more_pages <- FALSE  # No more pages
  }
}

#only two pages so just combining those
mixed_list <- rbind(mixed_studies[[3]], mixed_studies[[8]])

master_df <- rbind(master_df, mixed_list)

##Now we have a list of all of the attributes we need to search through these studies
# Defining keywords
keywords <- c("time", "time series", "time-series", "dynamic", "dynamics", "evolution", 
               "longitudinal", "time point", "time-point", "time points", "time-points", 
               "year", "years", "yearly", "month", "months", "monthly", "week", "weeks", 
               "weekly", "day", "days", "daily", "long term", "long-term", "continued", 
               "continuous", "monitoring", "monitored", "temporal", "fluctuation", 
               "fluctuations", "interval", "intervals", "differential abundance")  # Add your desired keywords
#Try adding "sequential" and "succession" to see if it increases the number of studies at all

# Create a regex pattern for matching any of the keywords
pattern <- paste(keywords, collapse = "|")  # "time|time-series"

# Filter the dataframe for matches in the Title or Abstract columns
filtered_df <- master_df[grepl(pattern, master_df$`study-name`, ignore.case = TRUE) |
                           grepl(pattern, master_df$`study-abstract`, ignore.case = TRUE), ]

##Now we have narrowed down from 6000 studies to 954 to 235. 
## The next step is parsing through these ones to decide whether or not they fit the criteria for inclusion in my study
## Specifically, longitudinal sampling of human microbiome data for 4 or more timepoints.
#write_csv(filtered_df, "filtered_df.csv") #1/10/2025
##COMMENTING THIS SO I DO NOT OVERWRITE MY CURRENT FILE

###Trying with 3 more keywords
keywords2 <- c("time", "time series", "time-series", "dynamic", "dynamics", "evolution", 
              "longitudinal", "time point", "time-point", "time points", "time-points", 
              "year", "years", "yearly", "month", "months", "monthly", "week", "weeks", 
              "weekly", "day", "days", "daily", "long term", "long-term", "continued", 
              "continuous", "monitoring", "monitored", "temporal", "fluctuation", 
              "fluctuations", "interval", "intervals", "differential abundance",
              "sequential", "succession", "successional")  # Add your desired keywords
#Try adding "sequential" and "succession/al" to see if it increases the number of studies at all (NOT DONE YET)

pattern2 <- paste(keywords2, collapse = "|")  # "time|time-series"

# Filter the dataframe for matches in the Title or Abstract columns
filtered_df2 <- master_df[grepl(pattern2, master_df$`study-name`, ignore.case = TRUE) |
                           grepl(pattern2, master_df$`study-abstract`, ignore.case = TRUE), ]

##Adding these three key words added 2 studies to the filtered_df list (we went from 235 to 237).
#Which studies were those?
unique_to_df2 <- anti_join(filtered_df2, filtered_df)

write.csv(unique_to_df2, "new_keywords_df.csv")
view(filtered_df2)
###This was revised 2/10/25. Needed to include the Mixed biome, which included 15 additional candidate studies
write.csv(filtered_df2, "Incl_mixed_df.csv")
##(I am going to add all of the mixed ones to my filtering document from this .csv file)

#The final_filter dataframe combines the filtered_df with the new_keywords_df manually
##I HAVE NOW MANUALLY GONE THROUGH ALL 252 studies presented here and determined 
##if they have longitudinal data for 4 or more time points (based on abstract or associated
##publication) and if these longitudinal samples include human microbiome data. 
##These data are stored in the final_filter.csv. The column "longitudinal-and-human"
##contains a 1 if the study fits these requirements and a 0 if not, along with an explanation
##why it does or does not fulfil the requirements for potential inclusion.

#Here we load in the final_filter dataset and filter it to get the 72 viable studies and 
#extract their additional metadata from MGnify
final_filter <- read.csv("final_filter.csv")
accession_list <- final_filter$accession[final_filter$longitudinal.and.human==1]
justification_list <- final_filter$notes[final_filter$longitudinal.and.human==1]
bioproject_list <- final_filter$bioproject[final_filter$longitudinal.and.human==1]

#pulling metadata for these studies
# Function to query MGnify metadata by accession
get_mgnify_metadata <- function(accession, resource = "studies") {
  # Base URL for MGnify API
  base_url <- "https://www.ebi.ac.uk/metagenomics/api/v1"
  
  # Construct the endpoint
  url <- paste0(base_url, "/", resource, "/", accession)
  
  # Make the API request
  response <- GET(url)
  
  data <- fromJSON(rawToChar(response$content))
}

# Retrieve metadata for each accession
metadata_list <- lapply(accession_list, function(acc) {
  get_mgnify_metadata(acc, resource = "studies")  # Replace "studies" with "projects", "samples", etc., as needed
})

###NOW I NEED TO EXTRACT THE RELEVANT METADATA!!!
metadata_df <- list()
metadata_df$accession <- metadata_list[[1]]$data$attributes$accession
metadata_df$`samples-count` <- metadata_list[[1]]$data$attributes$`samples-count`
metadata_df$`study-abstract` <- metadata_list[[1]]$data$attributes$`study-abstract`
metadata_df$`study-name` <- metadata_list[[1]]$data$attributes$`study-name`
metadata_df$biome <- metadata_list[[1]]$data$relationships$biomes$data$id
metadata_df <- as.data.frame(metadata_df)
for (i in 2:72) {
  metadata_temp <- list()
  metadata_temp$accession <- metadata_list[[i]]$data$attributes$accession
  metadata_temp$`samples-count` <- metadata_list[[i]]$data$attributes$`samples-count`
  metadata_temp$`study-abstract` <- metadata_list[[i]]$data$attributes$`study-abstract`
  metadata_temp$`study-name` <- metadata_list[[i]]$data$attributes$`study-name`
  metadata_temp$biome <- metadata_list[[i]]$data$relationships$biomes$data$id
  metadata_temp <- as.data.frame(metadata_temp)
  metadata_df <- rbind(metadata_df, metadata_temp)
}

#Some entries are repeated because there are two nested biomes reported in the biome
# column, so I am getting rid of one of each of those (the second iteration has the more
# specific biome)
final_meta_df <- metadata_df %>%
  group_by(accession) %>%
  filter(n() == 1 | row_number() == 2) %>%
  ungroup()

#confirming that the order of accession numbers is the same:
identical(final_meta_df$accession, accession_list)

##Adding back in my notes from before:
final_meta_df$notes <- justification_list
final_meta_df$bioproject <- bioproject_list

#Exporting this .csv for future reference
write.csv(final_meta_df, "final_metadata_df.csv")

#Converting the final metadata filter to include only the studies with 10+ timepoints, good metadata (manually collected data)
further <- read.csv("final_metadata_df.csv")
meta_analysis_list <- further[further$Inclusion_num_timepoints != 0, ]
meta_analysis_list <- meta_analysis_list[, -c(11, 16, 17)]
write.csv(meta_analysis_list, "meta_analysis_list.csv")
