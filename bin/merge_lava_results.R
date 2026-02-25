# Description: Merge LAVA results from partitioned runs

# Load packages -----------------------------------------------------------

library(stringr)
library(dplyr)

# Set arguments -----------------------------------------------------------

args <- commandArgs(TRUE)
lava_output_dir <- as.character(args[1])
run_id <- as.character(args[2])

# Set file prefix based on run_id
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# Find all partition files ------------------------------------------------

univ_files <- list.files(
  path = lava_output_dir,
  pattern = str_c(file_prefix, ".*_part[0-9]+\\.univ\\.lava\\.rds$"),
  full.names = TRUE
)

bivar_files <- list.files(
  path = lava_output_dir,
  pattern = str_c(file_prefix, ".*_part[0-9]+\\.bivar\\.lava\\.rds$"),
  full.names = TRUE
)

print(str_c("Found ", length(univ_files), " univariate partition files"))
print(str_c("Found ", length(bivar_files), " bivariate partition files"))

# Merge univariate results ------------------------------------------------

univar_merged <- list()

for (f in univ_files) {
  print(str_c("Reading: ", basename(f)))
  partition_data <- readRDS(f)

  # Append non-null elements to merged list
  for (i in seq_along(partition_data)) {
    if (!is.null(partition_data[[i]])) {
      univar_merged[[length(univar_merged) + 1]] <- partition_data[[i]]
    }
  }
}

# Convert list to data frame
univar_df <- bind_rows(univar_merged)
print(str_c("Merged ", nrow(univar_df), " univariate results"))

# Merge bivariate results -------------------------------------------------

bivar_merged <- list()

for (f in bivar_files) {
  print(str_c("Reading: ", basename(f)))
  partition_data <- readRDS(f)

  # Append non-null elements to merged list
  for (i in seq_along(partition_data)) {
    if (!is.null(partition_data[[i]])) {
      bivar_merged[[length(bivar_merged) + 1]] <- partition_data[[i]]
    }
  }
}

# Convert list to data frame
bivar_df <- bind_rows(bivar_merged)
print(str_c("Merged ", nrow(bivar_df), " bivariate results"))

# Save merged results -----------------------------------------------------

# Extract base output name from first file (remove partition suffix)
base_name <- str_replace(basename(univ_files[1]), "_part[0-9]+\\.univ\\.lava\\.rds$", "")

# Save as RDS and TSV:
if (nrow(univar_df) > 0) {
  saveRDS(
  univar_df,
  file = str_c(base_name, ".univ.lava.rds"))

  write.table(
  univar_df,
  file = str_c(base_name, ".univ.lava.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE)
}

if (nrow(bivar_df) > 0) {
  saveRDS(
    bivar_df,
    file = str_c(base_name, ".bivar.lava.rds")
  )

  write.table(
    bivar_df,
    file = str_c(base_name, ".bivar.lava.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

cat("LAVA results merged successfully!")
