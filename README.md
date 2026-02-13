# 🧬 nf-genetic-correlations

**Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics**  
Supports **LDSC** for genome-wide correlations and **LAVA** for local (regional) genetic correlations. The users can also follow-up the significant LAVA loci with Bayesian colocalization with **coloc**, to assess if there is a single shared causal variant.

---

## 📖 Overview

Take a look at the [workflow diagram](https://github.com/ape4fld/nf-genetic-correlations/blob/main/workflow.jpg) for a visual overview.

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

Some R packages need to be pre-installed in R version 4.3.1:

- Tidyverse: dplyr, tidyr, stringr, readr
- Others: here, data.table
- LAVA (```R via remotes::install_github()```)
- coloc (optional)

In Alliance Canada, you can follow these steps:

```bash
module load StdEnv/2023 r/4.3.1
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
R -e 'install.packages(c("dplyr", "tidyr", "stringr", "readr", "here", "data.table", "coloc"), repos="https://cloud.r-project.org/")'
R -e 'remotes::install_github("josefin-werme/LAVA")'
```

### 3. Inputs Required

---

#### 📁 a) GWAS Summary Statistics

- Accepted formats: `.tsv`, `.csv`, `.txt`, etc.
- Required columns (**names must match exactly**, order can vary and other columns will be ignored):
variant_id, effect_allele, other_allele, beta, standard_error, p_value, N

> ⚠️ `variant_id` must be rsIDs. This pipeline is optimized for harmonized summary stats from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

- Store your files in:
```bash
/genetic_correlations/data/sumstats/
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

- Store the metadata file at:
 ```bash
/genetic_correlations/data/
```
**Note:** an example of [metadata.txt](https://github.com/ape4fld/nf-genetic-correlations/blob/main/data/metadata.txt) is included, which can be edited. Additionally, if the user wants to run the analysis across a subset of the GWAS datasets, it is possible to do so by creating a new metadata file including only those datasets (and specify the file name with the --metadata flag - see  below '🚀 Running the Pipeline').

#### 📦 c) LD Reference Files

1. **LD Scores for LDSC**  
 Download and extract:

 ```bash
 wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
 tar -xvjf eur_w_ld_chr.tar.bz2
```

Place contents in:
 ```bash
/genetic_correlations/data/ld_reference/eur_w_ld_chr/
```

2. **1000 Genomes Reference or UK Biobank reference (for LAVA)**
Download European PLINK reference files as described in the LAVA reference guide or download UK Biobank reference files as described in the LAVA reference guide. Note that LAVA developers [highly recommend to use the UK Biobank reference file](https://www.preprints.org/manuscript/202507.0966).

Place 1000 Genomes Reference contents in:
 ```bash
/genetic_correlations/data/ld_reference/g1000_eur/
```
Place UK Biobank reference contents in:
 ```bash
/genetic_correlations/data/ld_reference/ukb_eur/
```
Note: The default LD reference file that is used is the UK Biobank, but the user can specify the LD source with the --lava_ref flag (options: 1KGP_EUR or UKB) - see ```run_nextflow.sh```).

### 4. 📦 LDSC Apptainer/Singularity Image

---

Download the LDSC container image from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15920751.svg)](https://doi.org/10.5281/zenodo.15920751)

```bash
# Download the image (1.2GB)
wget https://zenodo.org/records/15920751/files/ldsc_latest.sif
# Place it in the bin/ directory
mv ldsc_latest.sif bin/
```

### 5. ⚙️ Nextflow Configuration

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

The LAVA process is currently set to 5 hours, which works well for 4-5 phenotypes. However, **running time increases** with more datasets due to pairwise comparisons:
- 3 datasets = 3 pairs
- 5 datasets = 10 pairs  
- 10 datasets = 45 pairs

To adjust the time limit, modify in the [Nextflow config file](https://github.com/ape4fld/nf-genetic-correlations/blob/main/nextflow.config):
```nextflow
withLabel: lava {
    time = "5h"  // Increase for larger analyses
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
   - Options in Nextflow command (see ```run_nextflow.sh``` for an example):
     
   | Option              | Description                                      |
   |---------------------|--------------------------------------------------|
   | --run_id            | Give the specific run a prefix (optional)        |
   | --metadata          | Provide a different name to the metadata file (optional - default: metadata.txt) |
   | --lava-ref          | Specify LD reference for LAVA (optional; UKB or 1KGP_EUR - default: UKB) |
   | --coloc             | To include colocalization analysis follow-up (optional; true or false - default: false) |
   | --pvalue_LAVA_coloc | Provide p-value cutoff to define a significant local genetic correlation (for use with --coloc) (optional; - default: 0.05) |

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
- Optionally perform colocalization analysis using coloc R package (performed for pairs of traits with significant local genetic correlations (significance defined by the user)
- Output results to the `results/` directory
  
---

## 📊 Expected Outputs

The pipeline generates results in the following directory structure:

```
results/
├── formatted/                 # Formatted summary statistics
│   └── formatted_*.tsv        # One file per GWAS dataset
├── munged/                    # LDSC-ready files
│   └── *.sumstats.gz          # Munged summary statistics
├── ldsc_h2/                   # Heritability estimates
│   └── *.h2_results           # SNP-heritability for each trait
├── ldsc_rg/                   # Global genetic correlations
│   ├── *.rg_results           # Pairwise genetic correlations
│   └── all_rg_results.tsv     # Combined results table
└── LAVA/                      # Local genetic correlations
│   ├── *univ.lava.tsv         # Univariate test results (one line per trait)
│   └── *bivar.lava.tsv        # Bivariate test results (one line per trait pair)
└── coloc/                     # Colocalization results, if enabled by user
    ├── *coloc_all.txt         # Coloc results (one line per variant assessed across traits and loci)
    └── *coloc_summary.txt     # Coloc summary results (one line per locus)

data/LAVA/                     # LAVA input files
├── info_file.txt              # Trait information
└── sample_overlap.txt         # Sample overlap matrix
```

### Key Output Files:

- **`all_rg_results.tsv`**: Summary table with all global genetic correlations (rg), standard errors, p-values, and heritability estimates
- **`local_rg_*_univ.lava.tsv`**: Local heritability and association p-values for each genomic locus per trait
- **`local_rg_*_bivar.lava.tsv`**: Local genetic correlations between trait pairs at specific loci where both traits show significant univariate signals
