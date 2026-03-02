# Load packages -----------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

# Arguments -----------------------------------------------------------
args <- commandArgs(TRUE)
metadata_file <- as.character(args[1])
munged_dir <- as.character(args[2])
run_id <- as.character(args[3]) # run identifier for output file prefix

# Set output file prefix based on run_id
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# Read metadata and derive suffixes -------------------------------------------

metadata <- read.table(metadata_file, sep = "\t", header = TRUE)

suffixes <- metadata$dataset

# Get all pairwise combinations of suffixes
pairwise_combos <- combn(suffixes, 2)

# Convert to data frame with file paths and suffixes
df_out <- data.frame(
    suffix1 = pairwise_combos[1, ],
    suffix2 = pairwise_combos[2, ]
)

for (i in 1:nrow(df_out)) {
    sub1 <- metadata %>% dplyr::filter(., dataset == df_out$suffix1[i])
    sub_filename1 <- sub1$filename[1] %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".[:alpha:]+$")
    df_out$file1[i] <- stringr::str_c(munged_dir, "/", sub_filename1, ".sumstats.gz")
    sub2 <- metadata %>% dplyr::filter(., dataset == df_out$suffix2[i])
    sub_filename2 <- sub2$filename[1] %>% stringr::str_remove(".gz$") %>% stringr::str_remove(".[:alpha:]+$")
    df_out$file2[i] <- stringr::str_c(munged_dir, "/", sub_filename2, ".sumstats.gz")
}

df_out <- df_out %>% dplyr::select(file1, suffix1, file2, suffix2)


write.table(df_out, paste0(file_prefix, "pairs_to_test.tsv"),
            sep = "\t", row.names = F, quote = F, col.names = T)
