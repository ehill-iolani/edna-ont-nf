process BUILD_REPORT {
    label 'process_medium'  // nodes.dmp is held in memory when --taxdump is set
    container 'quay.io/biocontainers/pandas:2.2.1'
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path hit_files
    path consensus_files
    path taxdump    // NCBI taxdump dir/tar.gz, or [] to skip lineage columns

    output:
    path "abundance_table.tsv", emit: report
    path "run_qc_summary.html"

    script:
    def taxdump_arg = taxdump ? "--taxdump ${taxdump}" : ''
    """
    build_report.py \\
        --hits ${hit_files} \\
        --consensus ${consensus_files} \\
        --out-table abundance_table.tsv \\
        --out-html run_qc_summary.html \\
        --min-pident ${params.min_pident} \\
        ${taxdump_arg}
    """
}
