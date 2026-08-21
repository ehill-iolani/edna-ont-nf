#!/usr/bin/env nextflow
/*
 * eDNA ONT amplicon classification pipeline
 * Modern replacement for decona: quality-aware clustering (isONclust) +
 * spoa/racon/medaka consensus + reference/BLAST taxonomy assignment.
 */

nextflow.enable.dsl = 2

include { EDNA_AMPLICON } from './workflows/edna_amplicon.nf'
include { MERGE_FASTQ   } from './modules/merge_fastq.nf'

// ---- top-level params (override via -params-file or --flag) ----
// min_len / max_len / min_qual / cluster_id / min_cluster / enable_medaka /
// min_pident defaults all live in nextflow.config, not here -- anything read
// directly inside an included module or workflow script (rather than only
// inside this file's own `workflow` block) must be set there to be reliably
// visible by the time that module/workflow script runs (see nextflow.config)
params.input        = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir       = "results"
params.fwd_primer   = null
params.rev_primer   = null
params.taxdb        = null   // path to reference sequences FASTA; a BLAST db is built from this at runtime
params.help         = false

def helpMessage() {
    log.info """
    eDNA ONT amplicon pipeline
    ---------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --taxdb refs.fasta --outdir results

    Required:
      --input       CSV: sample,fastq_path
      --taxdb       Reference sequences FASTA (BLAST db is built from this each run)

    Key optional:
      --fwd_primer / --rev_primer   primer sequences for cutadapt trimming
      --min_len / --max_len / --min_qual   chopper filtering thresholds
      --cluster_id   isONclust/vsearch similarity threshold (default ${params.cluster_id})
      --min_cluster  minimum reads to polish a cluster (default ${params.min_cluster})
      --enable_medaka   use medaka-polished consensus instead of racon consensus (default ${params.enable_medaka})
      --min_pident   BLAST %identity below which a hit is low-confidence (default ${params.min_pident})
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

workflow {
    if (params.help || !params.input || !params.taxdb) {
        helpMessage()
        exit 0
    }

    taxdb_fasta_ch = Channel.fromPath(params.taxdb, checkIfExists: true)

    EDNA_AMPLICON(samplesheetToChannel(params.input), taxdb_fasta_ch)
}

// dev/debug entry points -- run a single step in isolation, e.g.:
//   nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
workflow MERGE_ONLY {
    if (!params.input) {
        log.error "MERGE_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
}
