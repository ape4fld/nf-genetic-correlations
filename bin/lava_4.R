# Description: run univariate and bivariate tests for GWAS traits using Nextflow
# Partition 4: loci 751 to 1000

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(dplyr)
library(tidyr)
library(stringr)

# Set arguments -----------------------------------------------------------

args <- commandArgs(TRUE)
ref_ld <- as.character(args[1])
loci_file = as.character(args[2])
lava_data_dir = as.character(args[3])
run_id = as.character(args[4]) # run identifier for input/output file prefixes

# Set file prefix based on run_id
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# Partition settings
start_locus <- 751
end_locus <- 1000

# -----------------------------------------------------------
# -----------------------------------------------------------
# Get array of phenotypes:

info <- read.table(here(lava_data_dir, paste0(file_prefix, "info_file.txt")), sep = "\t", header = TRUE)
phenotypes <- info$phenotype

# -----------------------------------------------------------
# -----------------------------------------------------------
# Set LAVA arguments:

args <-
  list(
    ref_prefix = ref_ld,
    loc_file = loci_file,
    info_file = here(lava_data_dir, paste0(file_prefix, "info_file.txt")),
    sample_overlap_file = here(lava_data_dir, paste0(file_prefix, "sample_overlap.txt")),
    phenotypes = phenotypes,
    output_filename = str_c(file_prefix, "local_rg_", str_c(phenotypes, collapse = ":"), "_part4")
  )

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)
input <-
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Adjust end_locus if it exceeds total loci
end_locus <- min(end_locus, n_loci)

# Main --------------------------------------------------------------------

# Print progress
print(str_c("Starting LAVA analysis for loci ", start_locus, " to ", end_locus))
n_loci_partition <- end_locus - start_locus + 1
progress <-
  quantile(
    x = 1:n_loci_partition,
    probs = seq(.05,1,.05)
  ) %>%
  ceiling()

# Set univariate threshold to 0.05/n_loci (use total loci for threshold)
univar_threshold <-
  0.05/n_loci

univar = bivar = list()

for (i in start_locus:end_locus) {

  local_i <- i - start_locus + 1
  if (local_i %in% progress) print(str_c("..", names(progress[which(progress==local_i)])))     # (printing progress)

  # Process locus
  locus <-
    LAVA::process.locus(
      loci[i,],
      input
    )

  # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs),
  # The !is.null(locus) check is necessary before calling the analysis functions.
  if (!is.null(locus)) {

    # extract some general locus info for the output
    loc_info <-
      data.frame(
        locus = locus$id,
        chr = locus$chr,
        start = locus$start,
        stop = locus$stop,
        n_snps = locus$n.snps,
        n_pcs = locus$K
      )

    # Run the univariate and bivariate tests
    loc_out <-
      LAVA::run.univ.bivar(
        locus,
        univ.thresh = univar_threshold
      )

    # Bind
    univar[[i]] <-
      loc_info %>%
      dplyr::bind_cols(loc_out$univ)

    if(!is.null(loc_out$bivar)){

      bivar[[i]] <-
        loc_info %>%
        dplyr::bind_cols(loc_out$bivar)

    }

  }

}

# Save data ---------------------------------------------------------------

saveRDS(
  univar,
  file = file.path(str_c(args$output_filename, ".univ.lava.rds"))
)
saveRDS(
  bivar,
  file = file.path(str_c(args$output_filename, ".bivar.lava.rds"))
)

cat("LAVA analysis partition 4 done!")
