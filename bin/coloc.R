# Description: run colocalization analysis as a follow-up for LAVA signals

# Packages -------------------------------------------------------

library(coloc)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(BSgenome)

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
dbSNP144_GRCh37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
dbSNP144_GRCh38 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

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

##### Helper functions:

# 1) rename columns to default
rename_to_default <- function(df, aliases) {
    lower_cols <- tolower(colnames(df))
    for (canonical in names(aliases)) {
        if (canonical %in% colnames(df)) next
        match_idx <- which(lower_cols %in% tolower(aliases[[canonical]]))[1]
        if (!is.na(match_idx)) {
        message(sprintf("Column '%s' matched as '%s'", colnames(df)[match_idx], canonical))
        colnames(df)[match_idx] <- canonical
        lower_cols[match_idx] <- canonical
        }
    }
    df
}

# 2) Function to convert rsIDs to chr:bp GRCh37 (Source: https://github.com/RHReynolds/colochelpR/blob/master/R/convert_rs_to_loc.R)

convert_rs_to_loc <- function(df, SNP_column, dbSNP){
  rs <- BSgenome::snpsById(dbSNP, df[[SNP_column]], ifnotfound = "drop") %>%
    as.data.frame() %>%
    tidyr::unite(col = "loc", seqnames, pos, sep = ":", remove = T) %>%
    dplyr::rename(rs = RefSNP_id) %>%
    dplyr::select(rs, loc)
  filter_vector <- c("rs")
  names(filter_vector) <- SNP_column
  df <- df %>%
    dplyr::inner_join(., rs, by = filter_vector) %>%
    tidyr::separate(., loc, c("chromosome", "base_pair_location"), sep = ":") %>%
    dplyr::mutate(., base_pair_location = as.numeric(base_pair_location),
                     chromosome = case_when(
                        chromosome == "X" ~ as.integer(23),
                        TRUE ~ as.integer(chromosome)
                     ))

  return(df)
}

# 3) Function to get rsIDs from chr:bp GRCh37 (Source: https://github.com/RHReynolds/colochelpR/blob/master/R/convert_rs_to_loc.R)

convert_loc_to_rs <- function(df, dbsnp){
  if(stringr::str_detect(df$CHR[1], "chr")){
    df <-
      df %>%
      dplyr::mutate(CHR = stringr::str_replace(CHR, "chr", ""))
  }
  df <-
    df %>%
    dplyr::mutate(CHR = as.factor(CHR),
                  BP = as.integer(BP),
                  CHR = case_when(
                    CHR == "23" ~ "X",
                    CHR == "x" ~ "X",
                    TRUE ~ CHR
                  ))
  df_gr <-
    GenomicRanges::makeGRangesFromDataFrame(df,
                                            keep.extra.columns = FALSE,
                                            ignore.strand = TRUE,
                                            seqinfo = NULL,
                                            seqnames.field = "CHR",
                                            start.field = "BP",
                                            end.field = "BP",
                                            starts.in.df.are.0based = FALSE)
  df_gr <-
    BSgenome::snpsByOverlaps(dbsnp, df_gr, minoverlap = 1L) %>%
    as.data.frame()
  combined <-
    df_gr %>%
    dplyr::rename(
      variant_id = RefSNP_id,
      CHR = seqnames,
      BP = pos,
      gr_strand = strand) %>%
    dplyr::right_join(df, by = c("CHR", "BP")) %>%
    dplyr::mutate(CHR = case_when(
                        CHR == "X" ~ as.integer(23),
                        TRUE ~ as.integer(CHR)
                     )) %>%
    dplyr::select(-gr_strand, -alleles_as_ambig) %>%
    dplyr::filter(., !(is.na(variant_id))) %>%
    dplyr::distinct(., variant_id, .keep_all = TRUE)
  return(combined)
}

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

    # Normalize column names to default form (case-insensitive alias matching) -----------
    col_aliases <- list(
        variant_id      = c("variant_id", "snp", "rsid", "rs_id", "snpid", "markername",
                            "marker", "id", "name", "variantid"),
        effect_allele   = c("effect_allele", "effectallele", "a1", "ea", "alt",
                            "allele1", "tested_allele", "coded_allele"),
        other_allele    = c("other_allele", "otherallele", "a2", "oa", "ref",
                            "allele2", "non_effect_allele", "nea"),
        beta            = c("beta", "b", "effect", "effect_size", "log_odds", "log_or"),
        standard_error  = c("standard_error", "se", "stderr", "std_err", "std_error"),
        p_value         = c("p_value", "p", "pval", "pvalue", "p_val", "p.value", "p-value")
    )

    df <- rename_to_default(df, col_aliases)

    # phen1
    phen1_sumstats <- fread(phen1_file)

    if (!(metadata_phen1$genome_version[1] %in% c("GRCh37", "GRCh38", "none"))) {
        warning(paste("Skipping: genome version for", phen1, "should be either: 'GRCh37', 'GRCh38' or 'none'."))
        next
    }
    
    # Three scenarios for genome version: none, GRCh37, or GRCh38

    # 1) none - means it doesn't have CHR and BP, so should map the rsIDs to GRCh37 positions.
    # 2) GRCh37 - does it have variant_id? If yes, keep as is, if not then get rsIDs from CHR and BP.
    # 3) GRCh38 - lift down positions to GRCh37, and assess if it has variant_id? If yes, keep as is, if not then get rsIDs from CHR and BP.

    chr_bp_aliases <- list(
        chromosome = c("chromosome", "chr", "chrom", "seqnames"),
        base_pair_location = c("base_pair_location", "bp", "pos", "position",
                               "basepairlocation", "base_pair", "genpos")
    )
    cols_check <- c("chromosome", "base_pair_location", "effect_allele",
                    "other_allele", "beta", "standard_error", "p_value")

    phen1_sumstats <- rename_to_default(phen1_sumstats, col_aliases)

    if (metadata_phen1$genome_version[1] == "none") {
        # No CHR/BP: map rsIDs to GRCh37 positions
        phen1_sumstats <- phen1_sumstats %>%
            convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
            dplyr::filter(., chromosome == chr &
                             base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
            dplyr::mutate(., N = metadata_phen1$N[1], varbeta = standard_error^2) %>%
            dplyr::filter(., !is.na(variant_id)) %>%
            dplyr::distinct(., variant_id, .keep_all = TRUE)
    }

    if (metadata_phen1$genome_version[1] == "GRCh37") {
        phen1_sumstats <- rename_to_default(phen1_sumstats, chr_bp_aliases)
        if (any(!(cols_check %in% colnames(phen1_sumstats)))) {
            stop(stringr::str_c("The summary statistics for ", phen1, " do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value.\n"))
        }
        if ("variant_id" %in% colnames(phen1_sumstats)) {
            phen1_sumstats <- phen1_sumstats %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen1$N[1], varbeta = standard_error^2) %>%
                dplyr::filter(., !is.na(variant_id)) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        } else {
            phen1_sumstats <- phen1_sumstats %>%
                dplyr::rename(CHR = chromosome, BP = base_pair_location) %>%
                convert_loc_to_rs(., dbSNP144_GRCh37) %>%
                dplyr::filter(., CHR == chr &
                                 BP >= start_pos & BP <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen1$N[1], varbeta = standard_error^2) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        }
    }

    if (metadata_phen1$genome_version[1] == "GRCh38") {
        if ("variant_id" %in% colnames(phen1_sumstats)) {
            # Use rsIDs to get GRCh37 positions directly
            phen1_sumstats <- phen1_sumstats %>%
                convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen1$N[1], varbeta = standard_error^2) %>%
                dplyr::filter(., !is.na(variant_id)) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        } else {
            # Get rsIDs from GRCh38 CHR/BP, then get GRCh37 positions
            phen1_sumstats <- rename_to_default(phen1_sumstats, chr_bp_aliases)
            if (any(!(cols_check %in% colnames(phen1_sumstats)))) {
                stop(stringr::str_c("The summary statistics for ", phen1, " do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value.\n"))
            }
            phen1_sumstats <- phen1_sumstats %>%
                dplyr::rename(CHR = chromosome, BP = base_pair_location) %>%
                convert_loc_to_rs(., dbSNP144_GRCh38) %>%
                convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos &
                                 base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen1$N[1], varbeta = standard_error^2) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        }
    }

    # phen2
    phen2_sumstats <- fread(phen2_file)

    if (!(metadata_phen2$genome_version[1] %in% c("GRCh37", "GRCh38", "none"))) {
        warning(paste("Skipping: genome version for", phen2, "should be either: 'GRCh37', 'GRCh38' or 'none'."))
        next
    }

    phen2_sumstats <- rename_to_default(phen2_sumstats, col_aliases)

     if (metadata_phen2$genome_version[1] == "none") {
        # No CHR/BP: map rsIDs to GRCh37 positions
        phen2_sumstats <- phen2_sumstats %>%
            convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
            dplyr::filter(., chromosome == chr &
                             base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
            dplyr::mutate(., N = metadata_phen2$N[1], varbeta = standard_error^2) %>%
            dplyr::filter(., !is.na(variant_id)) %>%
            dplyr::distinct(., variant_id, .keep_all = TRUE)
    }

    if (metadata_phen2$genome_version[1] == "GRCh37") {
        phen2_sumstats <- rename_to_default(phen2_sumstats, chr_bp_aliases)
        if (any(!(cols_check %in% colnames(phen2_sumstats)))) {
            stop(stringr::str_c("The summary statistics for ", phen1, " do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value.\n"))
        }
        if ("variant_id" %in% colnames(phen2_sumstats)) {
            phen2_sumstats <- phen2_sumstats %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen2$N[1], varbeta = standard_error^2) %>%
                dplyr::filter(., !is.na(variant_id)) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        } else {
            phen2_sumstats <- phen2_sumstats %>%
                dplyr::rename(CHR = chromosome, BP = base_pair_location) %>%
                convert_loc_to_rs(., dbSNP144_GRCh37) %>%
                dplyr::filter(., CHR == chr &
                                 BP >= start_pos & BP <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen2$N[1], varbeta = standard_error^2) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        }
    }

    if (metadata_phen2$genome_version[1] == "GRCh38") {
        if ("variant_id" %in% colnames(phen2_sumstats)) {
            # Use rsIDs to get GRCh37 positions directly
            phen2_sumstats <- phen2_sumstats %>%
                convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos & base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen2$N[1], varbeta = standard_error^2) %>%
                dplyr::filter(., !is.na(variant_id)) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        } else {
            # Get rsIDs from GRCh38 CHR/BP, then get GRCh37 positions
            phen2_sumstats <- rename_to_default(phen2_sumstats, chr_bp_aliases)
            if (any(!(cols_check %in% colnames(phen2_sumstats)))) {
                stop(stringr::str_c("The summary statistics for ", phen1, " do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\nchromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value.\n"))
            }
            phen2_sumstats <- phen2_sumstats %>%
                dplyr::rename(CHR = chromosome, BP = base_pair_location) %>%
                convert_loc_to_rs(., dbSNP144_GRCh38) %>%
                convert_rs_to_loc(., "variant_id", dbSNP144_GRCh37) %>%
                dplyr::filter(., chromosome == chr &
                                 base_pair_location >= start_pos &
                                 base_pair_location <= end_pos) %>%
                dplyr::mutate(., N = metadata_phen2$N[1], varbeta = standard_error^2) %>%
                dplyr::distinct(., variant_id, .keep_all = TRUE)
        }
    }

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
    fwrite(., stringr::str_c(file_prefix, "coloc_summary.tsv"), sep = "\t")

coloc_res_df %>% dplyr::bind_rows(., .id = "Traits") %>%
    tidyr::separate(., Traits, c("locus", "chr", "start_locus", "end_locus", "phen1", "phen2"), sep = "::") %>%
    fwrite(., stringr::str_c(file_prefix, "coloc_all.tsv"), sep = "\t")
