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

metadata_link <- metadata %>%
  mutate(filename_base = stringr::str_remove(filename, ".gz$") %>% stringr::str_remove(., ".[:alpha:]+$"))

info <- metadata_link %>%
  mutate(formatted_path = stringr::str_c(formatted_dir, "/formatted_", filename_base, ".tsv")) %>%
  dplyr::select(phenotype = dataset, cases, controls, filename = formatted_path)

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
  ) %>% distinct(., p1, p2, .keep_all = TRUE) 
  
# keep only traits in metadata file:
all_rg <- all_rg %>%
  dplyr::filter(p1 %in% metadata_link$filename_base) %>%
  dplyr::filter(p2 %in% metadata_link$filename_base) %>%
  dplyr::rename(p1_old = p1, p2_old = p2)

# link p1 and p2 to replace by metadata$dataset column instead of the basename of filename
for (i in 1:nrow(all_rg)) {
  sub1 <- metadata_link %>% dplyr::filter(., filename_base == all_rg$p1_old[i])
  all_rg$p1[i] <- sub1$dataset[1]
  sub2 <- metadata_link %>% dplyr::filter(., filename_base == all_rg$p2_old[i])
  all_rg$p2[i] <- sub2$dataset[1]
}

all_rg <- all_rg %>%
  dplyr::select(p1,p2,rg,se,z,p,h2_obs,h2_obs_se,h2_int, h2_int_se, gcov_int, gcov_int_se)

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


