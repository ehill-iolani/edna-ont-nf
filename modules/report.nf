process BUILD_REPORT {
    label 'process_medium'  // nodes.dmp is held in memory when --taxdump is set
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir "${params.outdir}/final_report", mode: 'copy', pattern: '{abundance_table.tsv,run_qc_summary.html}'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'consensus_by_confidence'

    input:
    path hit_files
    path consensus_files
    path taxdump    // NCBI taxdump dir/tar.gz, or [] to skip lineage columns

    output:
    path "abundance_table.tsv", emit: report
    path "run_qc_summary.html"
    path "consensus_by_confidence"

    script:
    def taxdump_arg = taxdump ? "--taxdump ${taxdump}" : ''
    """
    build_report.py \\
        --hits ${hit_files} \\
        --consensus ${consensus_files} \\
        --out-table abundance_table.tsv \\
        --out-html run_qc_summary.html \\
        --out-sorted consensus_by_confidence \\
        --min-pident ${params.min_pident} \\
        --min-abundance ${params.min_abundance} \\
        --min-rel-abundance ${params.min_rel_abundance} \\
        ${taxdump_arg}
    """
}
