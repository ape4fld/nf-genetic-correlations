nextflow.enable.dsl=2

// Run identifier - use this to namespace outputs when running different datasets
// Usage: nextflow run main_full.nf --run_id my_analysis --metadata path/to/metadata.txt --lava-ref UKB
params.run_id = params.run_id ?: ''
run_id = params.run_id

// Metadata file - can be overridden per run
metadata_file = params.metadata ? file(params.metadata) : file(params.data_dir).resolve('metadata.txt')

// LAVA reference panel (options: --lava_ref UKB or 1KGP_EUR, default: UKB)
params.lava_ref = params.lava_ref ?: 'UKB'

if (params.lava_ref == 'UKB') {
    ref_ld_chr = "${params.data_dir}/ld_reference/ukb_eur/lava-ukb-v1.1"
} else if (params.lava_ref == '1KGP_EUR') {
    ref_ld_chr = "${params.data_dir}/ld_reference/g1000_eur/g1000_eur"
} else {
    error "Invalid lava_ref value: ${params.lava_ref}. Must be 'UKB' or '1KGP_EUR'"
}

// Coloc: disabled by default, enable with --coloc true
params.coloc = params.coloc ?: false

// Coloc p-value threshold for LAVA bivariate results (default: 0.05)
params.pvalue_LAVA_coloc = params.pvalue_LAVA_coloc ?: 0.05
pvalue_LAVA_coloc = params.pvalue_LAVA_coloc

// Clean intermediate files: enabled by default, disable with --clean_files_only false
params.clean_files_only = params.clean_files_only ?: true

// Build channel from metadata file - only processes files listed in metadata
Channel
    .from(metadata_file.readLines().drop(1))  // skip header before creating channel
    .map { line ->
        def row = line.split('\t')
        def filename = row[1]
        def sampleN = row[2] as Integer
        def filepath = file("${params.data_dir}/sumstats/${filename}")
        def suffixName = filepath.getBaseName(filename.endsWith('.gz') ? 2 : 1)
        tuple(filepath, sampleN, suffixName)
    }
    .set { raw_sumstats_with_N }

// Create channel of LDSC rg results directory
Channel
    .value(file("${params.output_dir}/ldsc_rg"))
    .set { rg_dir_ch }

process FormatSumstats {
    label 'r_env'

    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("formatted_${suffixName}.tsv"), val(sampleN), val(suffixName)

    publishDir "${params.output_dir}/formatted", mode: 'copy'

    script:
    """
    Rscript ${params.bin_dir}/format_sumstats.R \
        ${file} \
        ${sampleN} \
        formatted_${suffixName}.tsv
    """
}

process MungeSumstats {
    label 'ldsc'

    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("${suffixName}.sumstats.gz"), val(suffixName)

    publishDir "${params.output_dir}/munged", mode: 'copy'

    script:
    """
        ldsc munge_sumstats \
            --sumstats ${file} \
            --N ${sampleN} \
            --merge-alleles ${params.hapmap_ref} \
            --chunksize 500000 \
            --out ${suffixName}
    """
}

process RunLDSC_h2 {
    label 'ldsc'

    input:
    tuple path(file), val(suffixName)

    output:
    path("ldsc_h2_${suffixName}.h2_results")

    publishDir "${params.output_dir}/ldsc_h2", mode: 'copy'

    script:
    """
        ldsc ldsc h2 \
            --h2 ${file} \
            --ref-ld-chr ${params.w_ld_chr} \
            --w-ld-chr ${params.w_ld_chr} \
            --out ldsc_h2_${suffixName}
    """
}

process PrepRg {
    label 'r_env'

    input:
    val ready  // signal that munging is complete
    val run_id

    output:
    path("*pairs_to_test.tsv")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
    Rscript ${params.bin_dir}/prep_pairs_for_rg.R \
        ${metadata_file} \
        ${params.output_dir}/munged \
        ${run_id}
    """
}

process RunLDSC_rg {
    label 'ldsc'

    input:
    tuple path(file1), val(suffix1), path(file2), val(suffix2)

    output:
    path("ldsc_rg_${suffix1}_${suffix2}.rg_results")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
        ldsc ldsc rg \
            --rg ${file1} \
            --rg ${file2} \
            --ref-ld-chr ${params.w_ld_chr} \
            --w-ld-chr ${params.w_ld_chr} \
            --out ldsc_rg_${suffix1}_${suffix2}
    """
}

process PrepLAVA {
    label 'r_env'

    input:
    path rg_dir
    val run_id

    output:
    path("*.txt"), emit: data_files
    path("*_all_rg_results.tsv"), emit: rg_all_results

    publishDir "${params.data_dir}/LAVA", mode: 'copy', pattern: "*.txt"
    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy', pattern: "*_all_rg_results.tsv"

script:
    """
    Rscript ${params.bin_dir}/prep_for_LAVA.R \
        ${metadata_file} \
        ${params.output_dir}/formatted \
        ${rg_dir} \
        ${run_id}
    """

}

process RunLAVA_1 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_1.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_2 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_2.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_3 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_3.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_4 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_4.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_5 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_5.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_6 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_6.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_7 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_7.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_8 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_8.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_9 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_9.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process RunLAVA_10 {
    label 'lava'

    input:
    val ready
    val run_id

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/lava_10.R \
            ${ref_ld_chr} \
            ${params.locus_file} \
            ${params.data_dir}/LAVA \
            ${run_id}
    """
}

process MergeLAVA {
    label 'r_env'

    input:
    path lava_rds_files
    val run_id

    output:
    path("*.rds"), emit: rds_files
    path("*.tsv"), emit: tsv_files

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
        Rscript ${params.bin_dir}/merge_lava_results.R \
            ${params.output_dir}/LAVA \
            ${run_id}
    """
}

process Coloc {
    label 'coloc'

    input:
    val ready
    val run_id

    output:
    path("*.txt"), emit: data_files, optional: true

    publishDir "${params.output_dir}/coloc", mode: 'copy', pattern: "*.txt"

    script:
    """
    Rscript ${params.bin_dir}/coloc.R \
        ${metadata_file} \
        ${params.data_dir}/sumstats \
        ${params.output_dir}/LAVA \
        ${run_id} \
        ${pvalue_LAVA_coloc}
    """
}

process Cleanup {

    input:
    val ready

    script:
    """
    rm -f ${params.output_dir}/LAVA/*_part*.rds
    rm -f ${params.output_dir}/formatted/formatted_*.tsv
    rm -f ${params.output_dir}/ldsc_rg/ldsc_rg_*.rg_results
    rm -f ${params.output_dir}/ldsc_rg/*_pairs_to_test.tsv
    rm -f ${params.output_dir}/munged/*.sumstats.gz
    """
}

workflow {
    
    // Step 1: Format
   formatted_sumstats = raw_sumstats_with_N | FormatSumstats

    // Step 2: Munge
    munged_sumstats = formatted_sumstats | MungeSumstats

    // Step 3: Estimate heritability (LDSC)
    munged_sumstats | RunLDSC_h2

    // Step 4. Prepare pairs of traits for RunLDSC_rg
    // Wait for munging to complete, then generate pairs from metadata
    munging_done = munged_sumstats.collect()

    pairs_test = PrepRg(munging_done, run_id)

    // Step 5: Run LDSC rg
    pairs_ch = pairs_test
        .map { file -> file.readLines().drop(1) }  // skip header
        .flatten()
        .map { line -> 
            def row = line.split('\t')
            tuple(row[0], row[1], row[2], row[3])
        }

    output_rg = pairs_ch | RunLDSC_rg

    // Step 6: Get sample overlap and info files for LAVA
    all_rg_results_ch = output_rg
    .collect()
    .map { it[0] }

    lava_input_ch = all_rg_results_ch
        .combine(rg_dir_ch)
        .map { results, rg_dir -> rg_dir }

    prep_lava = PrepLAVA(lava_input_ch, run_id)

    // Step 7: Run LAVA (10 partitions in parallel)
    lava_ready = prep_lava.data_files.collect()

    lava_1_out = RunLAVA_1(lava_ready, run_id)
    lava_2_out = RunLAVA_2(lava_ready, run_id)
    lava_3_out = RunLAVA_3(lava_ready, run_id)
    lava_4_out = RunLAVA_4(lava_ready, run_id)
    lava_5_out = RunLAVA_5(lava_ready, run_id)
    lava_6_out = RunLAVA_6(lava_ready, run_id)
    lava_7_out = RunLAVA_7(lava_ready, run_id)
    lava_8_out = RunLAVA_8(lava_ready, run_id)
    lava_9_out = RunLAVA_9(lava_ready, run_id)
    lava_10_out = RunLAVA_10(lava_ready, run_id)

    // Step 8: Merge LAVA results from all partitions
    all_lava_outputs = lava_1_out
        .mix(lava_2_out, lava_3_out, lava_4_out, lava_5_out,
             lava_6_out, lava_7_out, lava_8_out, lava_9_out, lava_10_out)
        .collect()

    merged_lava = MergeLAVA(all_lava_outputs, run_id)

    // Step 9: Run coloc on significant LAVA bivariate results (if enabled)
    merged_tsv_ready = merged_lava.tsv_files.collect()

    if (params.coloc == true) {
        coloc_ready = Coloc(merged_tsv_ready, run_id)
    }

    // Step 10: Cleanup intermediate files generated, keeping only final results in each results/ directory

    if (params.clean_files_only == true) {

        if (params.coloc == false) {
            merged_rds_ready = merged_lava.rds_files.collect()
            Cleanup(merged_rds_ready)
        }

        if (params.coloc == true) {
            Cleanup(coloc_ready)
        }
    }
}
