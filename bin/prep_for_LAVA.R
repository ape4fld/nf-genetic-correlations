# Description: 
# 1. Obtain the sample overlaps between pairs of traits for LAVA.
# 2. Obtain the info file for LAVA from the metadata file.

# Load packages -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readr)

# Set arguments -----------------------------------------------------------

args <- commandArgs(TRUE)
metadata_file <- as.character(args[1]) # metadata path
formatted_dir <- as.character(args[2]) # directory for formatted summary statistics
ldsc_rg <- as.character(args[3]) # results dir from LDSC rg
run_id <- as.character(args[4]) # run identifier for output file prefixes

# Set output file prefix based on run_id (adds underscore separator automatically)
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# -----------------------------------------------------------
# -----------------------------------------------------------
# Format metadata file into the info file for LAVA:

metadata <- read.table(metadata_file, sep = "\t", header = TRUE) 

info <- metadata %>%
  mutate(suffix = stringr::str_remove(filename, ".gz$") %>% stringr::str_remove(., ".[:alpha:]+$"),
        formatted_path = stringr::str_c(formatted_dir, "/formatted_", suffix, ".tsv")) %>%
  dplyr::select(phenotype = suffix, cases, controls, filename = formatted_path)
write.table(info, paste0(file_prefix, "info_file.txt"), sep = "\t", row.names = F, quote = F)

# -----------------------------------------------------------
# -----------------------------------------------------------

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    ldsc_rg,
    pattern = "*.rg_results",
    full.names = TRUE
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){
  
  files[[i]] <-
    read_table(
      file = file_paths[i]
    )
  
}

# extract rg values in matrix
all_rg <-
  files %>%
  rbindlist(., fill = TRUE) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz"),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz")
  ) %>%
  select(
    p1,
    p2,
    rg,
    se,
    z,
    p,
    h2_obs,
    h2_obs_se,
    h2_int,
    h2_int_se,
    gcov_int,
    gcov_int_se
  ) %>% distinct(., p1, p2, .keep_all = TRUE) %>%
  dplyr::filter(p1 %in% info$phenotype) %>%
  dplyr::filter(p2 %in% info$phenotype)

phenotypes <- info$phenotype

###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

n <- length(phenotypes)
covar_matrix <- matrix(NA,n,n)

rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

for (k in 1:length(phenotypes)){
  
  for(i in phenotypes) {
    for(j in phenotypes) {

      cat("Getting genetic covariance intercept for", i, "and", j, ".\n")

      if (i == j) {
        gcov_int <- 1
      } else {
        subset <- dplyr::filter(all_rg, p1 == i, p2 == j)
        if (nrow(subset) > 0) {
          gcov_int <- dplyr::filter(all_rg, p1 == i, p2 == j) %>% .[["gcov_int"]]
        } else {
          subset <- dplyr::filter(all_rg, p1 == j, p2 == i)
          if (nrow(subset) > 0) {
          gcov_int <- dplyr::filter(all_rg, p1 == j, p2 == i) %>% .[["gcov_int"]]
        } else {
          cat("Error when filling in covariance matrix with LDSC outputs for traits", i, " and ", j, ".\n")
          next
          }
        }
      }

      covar_matrix[i,j] <-
        gcov_int   
         
      cat("Done.\n")
    }
  }
  
  # Standardise the matrix
  covar_matrix <-
    covar_matrix %>%
    cov2cor() %>%
    round(digits = 5)
  
}

# Save data ---------------------------------------------------------------

write.table(
  covar_matrix,
  file = paste0(file_prefix, "sample_overlap.txt"),
  quote = F,
  row.names = phenotypes,
  sep = "\t"
)

write.table(
  all_rg,
  file = paste0(file_prefix, "all_rg_results.tsv"),
  quote = F,
  row.names = F,
  sep = "\t"
)


