# Load packages -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)
library(BSgenome)

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
dbSNP144_GRCh37 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
dbSNP144_GRCh38 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

# Arguments -----------------------------------------------------------
args <- commandArgs(TRUE)
input <- as.character(args[1])
N_samples <- as.character(args[2])
genome_version <- as.character(args[3])
output <- as.character(args[4])

# Function to get rsIDs from chr:bp GRCh37 (Source: https://github.com/RHReynolds/colochelpR/blob/master/R/convert_rs_to_loc.R)

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
    dplyr::select(-gr_strand, -alleles_as_ambig) %>%
    dplyr::filter(., !(is.na(variant_id))) %>%
    dplyr::distinct(., variant_id, .keep_all = TRUE)
  return(combined)
}

# Read harmonized summary stats -----------------------------------------------------------
df <- fread(input)

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

df <- rename_to_default(df, col_aliases)

# If variant_id is still missing, attempt to generate RSIDs from chromosome/position ----                                                              
if (!"variant_id" %in% colnames(df)) {                                                                                                                 
  chr_bp_aliases <- list(                                                                                                                              
  CHR = c("chr", "chromosome", "chrom", "seqnames"),                                                                                                 
  BP  = c("bp", "pos", "position", "base_pair_location", "basepairlocation", "base_pair", "genpos", "base_pair_position")                                                            
  )                                                                                                                                                    
  df <- rename_to_default(df, chr_bp_aliases)                                                                                                          
                                                                                                                                                      
  if (all(c("CHR", "BP") %in% colnames(df))) {
    message("variant_id not found. Generating RSIDs from CHR and BP using SNPlocs dbSNP version 144.")
    if (genome_version == "GRCh37") {
      df <- convert_loc_to_rs(df, dbSNP144_GRCh37)
    } else if (genome_version == "GRCh38") {
      df <- convert_loc_to_rs(df, dbSNP144_GRCh38)
    } else {
      stop(sprintf("Unsupported genome_version '%s'. Expected 'GRCh37' or 'GRCh38'.", genome_version))
    }
  }                                                                                                                                                    
}

# Check for required columns -----------------------------------------------------------

cols_expected = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

check_colnames <- cols_expected %in% colnames(df)

if ("FALSE" %in% check_colnames == TRUE) {
  missing <- cols_expected[!check_colnames]
  stop(sprintf(
    "Could not find column(s): %s\nAccepted aliases per column:\n%s\nColumn names in input: %s\n",
    paste(missing, collapse = ", "),
    paste(sapply(names(col_aliases), function(n) sprintf("  %s: %s", n, paste(col_aliases[[n]], collapse = ", "))), collapse = "\n"),
    paste(colnames(df), collapse = ", ")
  ))
} else {
  df <- df %>%
    dplyr::select(
      variant_id,
      effect_allele,
      other_allele,
      beta,
      standard_error,
      p_value
    )
}

# Check if rsids are indeed rsids -----------------------------------------------------------
df <- df %>%
  filter(grepl("^rs", variant_id))

if (nrow(df) == 0) {
  warning("No rows with rsid starting with 'rs' found. No output will be generated.")
}

# Include sample size in sumstats -----------------------------------------------------------
df <- df %>% mutate(N = N_samples)

# Calculate Z-score and select required columns -----------------------------------------------------------
df_formatted <- df %>%
  mutate(Z = beta / standard_error) %>%
  transmute(
    SNP = variant_id,
    N = N,
    Z = Z,
    A1 = effect_allele,
    A2 = other_allele,
    P = p_value
  )

# Write to output -----------------------------------------------------------
if (nrow(df_formatted) > 0) {
  fwrite(df_formatted, output, sep = "\t")
}