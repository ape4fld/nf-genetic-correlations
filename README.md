# 🧬 nf-genetic-correlations

**Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics**  
Supports **LDSC** for genome-wide correlations and **LAVA** for local (regional) genetic correlations. The users can also follow-up the significant LAVA loci with Bayesian colocalization with **coloc**, to assess if there is a single shared causal variant.

---

## 📖 Overview

Take a look at the [workflow diagram](https://github.com/ape4fld/nf-genetic-correlations/blob/main/workflow.png) for a visual overview.

This pipeline processes **harmonized GWAS summary statistics** (restricted to **European ancestry** for now) and computes:

- **Global genetic correlations** using [LDSC](https://github.com/bulik/ldsc) (also computes SNP-based heritability)
- **Local genetic correlations** using [LAVA](https://github.com/josefin-werme/LAVA) (computes univariate and bivariate tests)
- **Bayesian colocalization** (optional; if enabled by user) using [coloc](https://chr1swallace.github.io/coloc/) (optional; across loci with significant regional genetic correlations)

There are several advantages of using the pipeline:
1) Given that it uses an LDSC .sif image, there is no need to load old python versions (< v3) to run LDSC.
2) The pipeline formats and adapts the GWAS summary statistics for each tool.
3) The user does not need to prepare additional files, other than a metadata file.
4) It partitions LAVA loci so that they all run in parallel, which significantly reduces the running time.
5) It is reproducible and the user can easily re-run the analysis by adding/removing GWAS datasets from the metadata file.

---

## 🚀 Getting Started

### 1. Install

```bash
git clone https://github.com/ape4fld/nf-genetic-correlations.git
cd nf-genetic-correlations
```

### 2. Dependencies

### a) 📦 LDSC Apptainer/Singularity Image

---

Download the LDSC container image from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15920751.svg)](https://doi.org/10.5281/zenodo.15920751)

```bash
# Download the image (1.2GB)
wget https://zenodo.org/records/15920751/files/ldsc_latest.sif
# Place it in the bin/ directory
mv ldsc_latest.sif bin/
```

### b) 📦 R environment Apptainer/Singularity Image

---

Download the R container image from Zenodo (v2):

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20817347.svg)](https://doi.org/10.5281/zenodo.20817347)

```bash
# Download the image (420Mb)
wget https://zenodo.org/records/18683118/files/r_packages.sif
# Place it in the bin/ directory
mv r_packages.sif bin/
```

### 3. Inputs Required

---

#### 📁 a) GWAS Summary Statistics

- Accepted formats: `.tsv`, `.csv`, `.txt`, etc.
- Required columns - only one of each - (order of columns can vary and additional columns will be ignored):

| Column          | Aliases accepted (agnostic to lower/upper case) |
|-----------------|----------------------------------------------------|
| `variant_id`    | `snp`, `rsid`, `rs_id`, `snpid`, `markername`, `marker`, `id`, `name`, `variantid` |
| `chromosome`    | `chr`, `chrom`, `seqnames` |
| `base_pair_location` | `pos`, `position`, `base_pair_location`, `basepairlocation`, `base_pair`, `genpos`, `base_pair_position` |
| `effect_allele` | `effectallele`, `a1`, `ea`, `alt`, `allele1`, `tested_allele`, `coded_allele`      |
| `other_allele`  | `otherallele`, `a2`, `oa`, `ref`, `allele2`, `non_effect_allele`, `nea`  |
| `beta`          | `b`, `effect`, `effect_size`, `log_odds`, `log_or`   |
| `standard_error`| `se`, `stderr`, `std_err`, `std_error` |
| `p_value`       | `p`, `pval`, `pvalue`, `p_val`, `p.value`, `p-value` |
  
⚠️ Notes: 
- **`variant_id`** must be rsIDs 
- If `variant_id` is absent, **`chromosome`** and **`base_pair_location`** must be included (and specify the genome reference version in the `metadata.txt` file). Make sure you select the correct genome version for each summary statistics (either `GRCh37`, `GRCh38`).
- Similarly, if `chromosome` and `base_pair_location` are absent, then `variant_id` must be present. The pipeline will map rsIDs to GRCh37 genomic positions (and `genome_version` in the metadata file should be then set to `none`).

- **Create a directory for the summary statistics:**
```bash
mkdir data/sumstats/
# Store the files here:
/nf-genetic-correlations/data/sumstats/
```

#### 📝 b) Metadata File

Create a single file named `metadata.txt`, tab-separated, with the following columns:

| Column     | Description                                        |
|------------|----------------------------------------------------|
| `dataset`  | Short name for each dataset                        |
| `filename` | File name of the GWAS summary statistics file      |
| `N`        | Total sample size (use max if per-variant varies)  |
| `cases`    | Number of cases (use `NA` for continuous traits)   |
| `controls` | Number of controls (use `NA` for continuous traits)|
| `genome_version` | One of the following: `GRCh37`, `GRCh38` or `none` (if there are no CHR and BP locations) |

- Store the metadata file in:
 ```bash
/nf-genetic-correlations/data/
```
**Note:** an example of [metadata.txt](https://github.com/ape4fld/nf-genetic-correlations/blob/main/data/metadata.txt) is included, which can be edited. Additionally, if the user wants to run the analysis across a subset of the GWAS datasets, it is possible to do so by creating a new metadata file including only those datasets (and specify the file name with the --metadata flag - see  below '🚀 Running the Pipeline').

#### 📦 c) LD Reference Files

1. **LD Scores for LDSC**  
Download the LD scores from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18749273.svg)](https://doi.org/10.5281/zenodo.18749273)

```bash
# Download the compressed directory (65.9Mb)
wget -O eur_w_ld_chr.tar.gz https://zenodo.org/records/18749273/files/eur_w_ld_chr.tar.gz
# Uncompress the directory
tar -xf eur_w_ld_chr.tar.gz
```

Place /eur_w_ld_chr in ld_reference directory:
 ```bash
mv eur_w_ld_chr/ ./data/ld_reference/
```

2. **1000 Genomes Reference or UK Biobank reference (for LAVA)**
Download European PLINK reference files as described in the [LAVA reference guide](https://github.com/josefin-werme/LAVA/blob/main/REFERENCE.md) or download UK Biobank reference files as described in the [LAVA reference guide](https://github.com/josefin-werme/LAVA/blob/main/REFERENCE.md). Note that LAVA developers [highly recommend to use the UK Biobank reference file](https://www.preprints.org/manuscript/202507.0966).

Place 1000 Genomes Reference contents in:
 ```bash
/nf-genetic-correlations/data/ld_reference/g1000_eur/
```
Place UK Biobank reference contents in:
 ```bash
/nf-genetic-correlations/data/ld_reference/ukb_eur/
```
Note: The default LD reference file that is used is the UK Biobank, but the user can specify the LD source with the --lava_ref flag (options: 1KGP_EUR or UKB) - see ```run_nextflow.sh```).

### 4. ⚙️ Nextflow Configuration

---

The pipeline uses **relative paths** by default, making it portable across different systems. The configuration is set up for **[Digital Research Alliance Canada](https://www.alliancecan.ca/en) clusters** but can be adapted for other environments.

#### Minimal Configuration Required:

1. **For Digital Research Alliance Canada users**, update the SLURM account in the [Nextflow config file](https://github.com/ape4fld/nf-genetic-correlations/blob/main/nextflow.config):
```nextflow
process.clusterOptions = '--account=def-xxxxx'  // Replace with your allocation
```

2. **For other HPC/local systems**, you may need to:
- Change the `executor` from 'slurm' to your system (e.g., 'local', 'sge', 'pbs')
- Adjust resource allocations (memory, CPUs, time)

#### ⏱️ Time Considerations for LAVA:

The LAVA process is currently set to 1 hour, which works well for 4-5 phenotypes. However, **running time increases** with more datasets due to pairwise comparisons:
- 3 datasets = 3 pairs
- 5 datasets = 10 pairs  
- 10 datasets = 45 pairs

To adjust the time limit, modify in the [Nextflow config file](https://github.com/ape4fld/nf-genetic-correlations/blob/main/nextflow.config):
```nextflow
withLabel: lava {
    time = "1h"  // Increase for larger analyses
}
```

#### Default Directory Structure:
The pipeline expects this structure relative to where your [main_full.nf](https://github.com/ape4fld/nf-genetic-correlations/blob/main/main_full.nf) file is located.

---

## 🚀 Running the Pipeline

Once you've completed the setup and configuration, you can run the pipeline:

### For Alliance Canada Users:

1. **Edit the SLURM script** [run_nextflow.sh](https://github.com/ape4fld/nf-genetic-correlations/blob/main/run_nextflow.sh):
   - Replace `def-xxxxx` with your compute allocation
   - Options in Nextflow command (see ```run_nextflow.sh``` for an example).
   - All flags are optional.
     
     
   | Flag                | Description                                      | Default    |
   |---------------------|--------------------------------------------------|------------|
   | --run_id            | Give the specific run a prefix | no prefix |
   | --metadata          | Provide a different name to the metadata file | metadata.txt |
   | --lava-ref          | Specify LD reference for LAVA ('1KGP_EUR' or 'UKB') | 'UKB' |
   | --coloc             | Include colocalization analysis (true or false) | false |
   | --pvalue_LAVA_coloc | Provide p-value cutoff for a significant local genetic correlation (for use with --coloc) | 0.05 |
   | --clean_files_only  | Delete intermediate files generated, keep only final results | true |

3. **Submit the job**:
   ```bash
   sbatch run_nextflow.sh
   ```

### For Other HPC/Local Systems:

Run Nextflow directly:
```bash
nextflow run main_full.nf -profile <your_profile> -resume \
     --run_id analysis1 \
     --metadata ./data/metadata.txt \
     --lava_ref 'UKB'
```

The nf-genetic-correlations pipeline will:
- Process your GWAS summary statistics
- Calculate global genetic correlations using LDSC
- Calculate local genetic correlations using LAVA (the bivariate test is only performed for loci that passed a Bonferroni-corrected univariate test (i.e., pvalue < 0.05/2,945))
- Optionally perform colocalization analysis using coloc R package (performed for pairs of traits with significant local genetic correlations; significance defined by the user)
- Output results to the `results/` directory
  
---

## 📊 Expected Outputs

The pipeline generates results in the following directory structure:

```
results/
├── formatted/                 # Formatted summary statistics
│   └── formatted_*.tsv        # One file per GWAS dataset (kept only if --clean_files_only = false)
├── munged/                    # LDSC-ready files
│   └── *.sumstats.gz          # Munged summary statistics (kept only if --clean_files_only = false)
├── ldsc_h2/                   # Heritability estimates
│   └── *.h2_results           # SNP-heritability for each trait
├── ldsc_rg/                   # Global genetic correlations
│   ├── *.rg_results           # Pairwise genetic correlations (kept only if --clean_files_only = false)
│   └── *_all_rg_results.tsv     # Combined results table
└── LAVA/                      # Local genetic correlations
│   ├── *univ.lava.tsv         # Univariate test results (one line per trait)
│   └── *bivar.lava.tsv        # Bivariate test results (one line per trait pair)
└── coloc/                     # Colocalization results, if enabled by user
    ├── *_coloc_all.txt         # Coloc results (one line per variant assessed across traits and loci)
    └── *_coloc_summary.txt     # Coloc summary results (one line per locus)

data/LAVA/                     # LAVA input files
├── info_file.txt              # Trait information
└── sample_overlap.txt         # Sample overlap matrix
```

### Questions or issues:
Please [open an issue](https://github.com/ape4fld/nf-genetic-correlations/issues) if you have questions/suggestions about this Nextflow pipeline, if you encounter problems or if you find a bug! Remember to include your input and output for easier debugging; the more information the better.
