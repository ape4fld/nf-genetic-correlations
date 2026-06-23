**Example for running the pipeline nf-genetic-correlations**

The first step to run the pipeline using the example here is to follow the instructions in the [main documentation page](https://github.com/ape4fld/nf-genetic-correlations) for cloning the repository and installing dependencies.

The example uses six GWAS summary statistics that can be downloaded from the GWAS catalog and placed in the `data/` directory:

1. [Alzheimer's disease](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST013001-GCST014000/GCST013196/GCST013196.tsv.gz)
2. [Parkinson's disease](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/GCST009325.tsv)
3. [Amyotrophic Lateral Sclerosis](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz)
4. [Major Depressive Disorder](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726344/GCST90726344.tsv)
5. [Schizophrenia](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435854/GCST90435854.tsv.gz)
6. [Bipolar Disorder](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90042001-GCST90043000/GCST90042719/GCST90042719_buildGRCh37.tsv.gz)

A `metadata.txt` ready-to-use file is provided within `example/` for the GWAS summary statistics above (user must place this file within data/ too).

Next, the user can download the LD reference files for LDSC and LAVA, as indicated in the [main documentation page](https://github.com/ape4fld/nf-genetic-correlations).

Once all data has been downloaded and placed in the correct directories of the cloned repository, the user may adapt the [nextflow.config](https://github.com/ape4fld/nf-genetic-correlations/blob/main/nextflow.config) file to their own HPC system:

1. Modify the process executor scheduler accordingly (line 37) - default: slurm.
2. Modify the custom_profile (lines 93-97), with its own memory, time and CPU limits, according to the user's HPC system.

This example was run using an HPC system with the default memory parameters and the whole process took approximately two hours.

Finally, the user can execute the pipeline using the following command in a bash terminal:

```bash
nextflow run main_full.nf -profile custom_profile \
     --run_id example \
     --metadata ./data/metadata.txt \
     --lava_ref 'UKB' \
     --coloc true \
     --pvalue_LAVA_coloc 0.05
```

**Expected outputs**
After finishing running the pipeline with the parameters set above, the user should be able to see in their `results/` directory the following subdirectories, each with corresponding files, as explained in the [main documentation page](https://github.com/ape4fld/nf-genetic-correlations):
1) ldsc_h2
2) ldsc_rg
3) LAVA
4) coloc

