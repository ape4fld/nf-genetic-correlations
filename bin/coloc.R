# Description: run colocalization analysis as a follow-up for LAVA signals

# Packages -------------------------------------------------------

library(coloc)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# Set arguments -----------------------------------------------------------

args <- commandArgs(TRUE)
metadata_file <- as.character(args[1]) # metadata path
sumstats_dir <- as.character(args[2]) # directory for initial summary statistics
lava_bivar_dir <- as.character(args[3]) # bivariate results dir from LAVA
run_id <- as.character(args[4]) # run identifier for output file prefixes
pval_cutoff <- as.numeric(args[5]) # p-value cutoff on LAVA bivariate results for performing coloc

# Set file prefix based on run_id
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# Coloc priors
p1 = 1e-04
p2 = 1e-04
p12 = 5e-06

coloc_results_summ = list()
coloc_results_res = list()
trait_pairs <- array()

# Read datasets -------------------------------------------------------

# read metadata file to get sample sizes
metadata <- fread(metadata_file)

# read LAVA bivariate results
lava_bivar_path <- list.files(lava_bivar_dir, pattern = paste0("^", file_prefix, ".*\\.bivar.lava.tsv$"), full = TRUE)
lava_bivar = list()

if (length(lava_bivar_path) == 0) {
    stop("No LAVA bivariate result files found in directory. This is possibly true if there weren't any pairs of traits to test for local genetic correlations.")
}

for (i in seq_along(lava_bivar_path)) {
    tmp <- tryCatch(
        read.table(lava_bivar_path[i], sep = "\t", header = TRUE),
        error = function(e) NULL
    )
    if (!is.null(tmp) && nrow(tmp) > 0) {
        lava_bivar[[length(lava_bivar) + 1]] <- tmp
    }
}

test_loci <- rbindlist(lava_bivar) %>%
    dplyr::filter(., p < pval_cutoff)

if (nrow(test_loci) == 0) {
    stop("There are no significant LAVA bivariate results. Finishing analysis without creating outputs.")
}

for (i in 1:nrow(test_loci)) {

    locus = test_loci$locus[i]
    chr = test_loci$chr[i]
    start_pos = test_loci$start[i]
    end_pos = test_loci$stop[i]
    phen1 = test_loci$phen1[i]
    phen2 = test_loci$phen2[i]

    metadata_phen1 <- dplyr::filter(metadata, dataset == phen1)
    metadata_phen2 <- dplyr::filter(metadata, dataset == phen2)

    # Validate metadata exists                                                                              
      if (nrow(metadata_phen1) == 0 || is.na(metadata_phen1$filename[1])) {                                   
          warning(paste("Skipping: no metadata found for phenotype", phen1))                                  
          next                                                                                                
      }                                                                                                       
      if (nrow(metadata_phen2) == 0 || is.na(metadata_phen2$filename[1])) {                                   
          warning(paste("Skipping: no metadata found for phenotype", phen2))                                  
          next                                                                                                
      }                                                                                                       
                                                                                                              
      # Build file paths                                                                                      
      phen1_file <- file.path(sumstats_dir, metadata_phen1$filename[1])                                       
      phen2_file <- file.path(sumstats_dir, metadata_phen2$filename[1])                                       
                                                                                                              
      if (!file.exists(phen1_file)) {                                                                         
          warning(paste("Skipping: file not found:", phen1_file))                                             
          next                                                                                                
      }                                                                                                       
      if (!file.exists(phen2_file)) {                                                                         
          warning(paste("Skipping: file not found:", phen2_file))                                             
          next                                                                                                
      }

    # check columns in summary statistics:
    cols_expected = c("chromosome", "base_pair_location", "variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

    # phen1
    phen1_sumstats <- fread(phen1_file) 
    check_colnames_phen1 <- cols_expected %in% colnames(phen1_sumstats)
    
    if ("FALSE" %in% check_colnames_phen1 == TRUE) {
        stop(stringr::str_c("The summary statistics for ", phen1," do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, variant_id, effect_allele, other_allele, beta, standard_error, p_value.\n"))
    }

    phen1_sumstats <- phen1_sumstats %>%
        dplyr::filter(., chromosome == chr & base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
        mutate(., N = metadata_phen1$N[1],
                varbeta = standard_error^2) %>%
        dplyr::filter(., !is.na(variant_id)) %>%
        distinct(., variant_id, .keep_all = TRUE)

    # phen2
    phen2_sumstats <- fread(phen2_file)
    check_colnames_phen2 <- cols_expected %in% colnames(phen2_sumstats)
    
    if ("FALSE" %in% check_colnames_phen2 == TRUE) {
        stop(stringr::str_c("The summary statistics for ", phen2," do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, variant_id, effect_allele, other_allele, beta, standard_error, p_value.\n"))
    }

    phen2_sumstats <- phen2_sumstats %>%
        dplyr::filter(., chromosome == chr & base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
        mutate(., N = metadata_phen2$N[1],
                varbeta = standard_error^2) %>%
        dplyr::filter(., !is.na(variant_id)) %>%
        distinct(., variant_id, .keep_all = TRUE)

    # Run coloc depending on which trait is quantitative or case/control
    # Assumes variance (dY) for quantitative traits = 1

    if (is.na(metadata_phen1$cases[1]) == TRUE & is.na(metadata_phen2$cases[1]) == TRUE) {
        # 1) Both phen1 and phen2 are quantitative:
        coloc_results <- coloc.abf(dataset1 = list(type = "quant",
                                                snp = phen1_sumstats$variant_id,
                                                beta = phen1_sumstats$beta,
                                                varbeta = phen1_sumstats$varbeta,
                                                pvalues = phen1_sumstats$p_value,
                                                N = phen1_sumstats$N,
                                                sdY = 1),
                                dataset2 = list(type = "quant",
                                                snp = phen2_sumstats$variant_id,
                                                beta = phen2_sumstats$beta,
                                                varbeta = phen2_sumstats$varbeta,
                                                pvalues = phen2_sumstats$p_value,
                                                N = phen2_sumstats$N,
                                                sdY = 1),
                                p1 = p1, p2 = p2, p12 = p12
        )
    } 
    
    if (is.na(metadata_phen1$cases[1]) == TRUE & is.na(metadata_phen2$cases[1]) == FALSE) {
        # 2) Phen1 is quantitative and phen2 is case control:
        s_phen2 = metadata_phen2$cases[1] / (metadata_phen2$cases[1] + metadata_phen2$controls[1])

        coloc_results <- coloc.abf(dataset1 = list(type = "quant",
                                                snp = phen1_sumstats$variant_id,
                                                beta = phen1_sumstats$beta,
                                                varbeta = phen1_sumstats$varbeta,
                                                pvalues = phen1_sumstats$p_value,
                                                N = phen1_sumstats$N,
                                                sdY = 1),
                                dataset2 = list(type = "cc",
                                                snp = phen2_sumstats$variant_id,
                                                beta = phen2_sumstats$beta,
                                                varbeta = phen2_sumstats$varbeta,
                                                pvalues = phen2_sumstats$p_value,
                                                N = phen2_sumstats$N,
                                                s = s_phen2),
                                p1 = p1, p2 = p2, p12 = p12
        )
    }
    
    if (is.na(metadata_phen1$cases[1]) == FALSE & is.na(metadata_phen2$cases[1]) == TRUE) {
        # 3) Phen1 is case control and phen2 is quantitative:
        s_phen1 = metadata_phen1$cases[1] / (metadata_phen1$cases[1] + metadata_phen1$controls[1])

        coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                                snp = phen1_sumstats$variant_id,
                                                beta = phen1_sumstats$beta,
                                                varbeta = phen1_sumstats$varbeta,
                                                pvalues = phen1_sumstats$p_value,
                                                N = phen1_sumstats$N,
                                                s = s_phen1),
                                dataset2 = list(type = "quant",
                                                snp = phen2_sumstats$variant_id,
                                                beta = phen2_sumstats$beta,
                                                varbeta = phen2_sumstats$varbeta,
                                                pvalues = phen2_sumstats$p_value,
                                                N = phen2_sumstats$N,
                                                sdY = 1),
                                p1 = p1, p2 = p2, p12 = p12
        )
    }
    
    if (is.na(metadata_phen1$cases[1]) == FALSE & is.na(metadata_phen2$cases[1]) == FALSE) {
        # 4) Phen1 and Phen2 are both case control:
        s_phen1 = metadata_phen1$cases[1] / (metadata_phen1$cases[1] + metadata_phen1$controls[1])
        s_phen2 = metadata_phen2$cases[1] / (metadata_phen2$cases[1] + metadata_phen2$controls[1])

        coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                                snp = phen1_sumstats$variant_id,
                                                beta = phen1_sumstats$beta,
                                                varbeta = phen1_sumstats$varbeta,
                                                pvalues = phen1_sumstats$p_value,
                                                N = phen1_sumstats$N,
                                                s = s_phen1),
                                dataset2 = list(type = "cc",
                                                snp = phen2_sumstats$variant_id,
                                                beta = phen2_sumstats$beta,
                                                varbeta = phen2_sumstats$varbeta,
                                                pvalues = phen2_sumstats$p_value,
                                                N = phen2_sumstats$N,
                                                s = s_phen2),
                                p1 = p1, p2 = p2, p12 = p12
        )
    }

    coloc_results_summ[[i]] <- coloc_results$summary
    coloc_results_res[[i]] <- coloc_results$results
    trait_pairs[i] = stringr::str_c(locus, "::", chr, "::", start_pos, "::", end_pos, "::", phen1, "::", phen2)

}

coloc_summ_df <- setNames(coloc_results_summ, nm = trait_pairs)
coloc_res_df <- setNames(coloc_results_res, nm = trait_pairs)

coloc_summ_df %>% dplyr::bind_rows(., .id = "Traits") %>%
    tidyr::separate(., Traits, c("locus", "chr", "start_locus", "end_locus", "phen1", "phen2"), sep = "::") %>%
    write.table(., stringr::str_c(file_prefix, "coloc_summary.txt"), sep = "\t", row.names = F, quote = F)

coloc_res_df %>% dplyr::bind_rows(., .id = "Traits") %>%
    tidyr::separate(., Traits, c("locus", "chr", "start_locus", "end_locus", "phen1", "phen2"), sep = "::") %>%
    write.table(., stringr::str_c(file_prefix, "coloc_all.txt"), sep = "\t", row.names = F, quote = F)
