include { MAKEBLASTDB      } from '../modules/makeblastdb.nf'
include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from '../modules/read_stats.nf'
include { CUTADAPT         } from '../modules/cutadapt.nf'
include { VSEARCH_CLUSTER  } from '../modules/vsearch_cluster.nf'
include { SPOA_CONSENSUS   } from '../modules/spoa_consensus.nf'
include { MINIMAP2_ALIGN   } from '../modules/minimap2_align.nf'
include { RACON            } from '../modules/racon.nf'
include { MEDAKA           } from '../modules/medaka.nf'
include { MERGE_CONSENSUS  } from '../modules/merge_consensus.nf'
include { BLAST_TAX        } from '../modules/blast_tax.nf'
include { BUILD_REPORT     } from '../modules/report.nf'

workflow EDNA_AMPLICON {

    take:
    reads_ch    // tuple(sample, fastq)
    ref_fasta   // path to reference sequences fasta for BLAST taxonomy
    taxdump     // NCBI taxdump dir/tar.gz, or [] to skip lineage columns

    main:
    // 1. build a BLAST db from the reference fasta (once per run, independent
    //    of per-sample steps below)
    MAKEBLASTDB(ref_fasta)

    // 2. merge multi-part fastq(.gz) files per barcode into one file per sample
    MERGE_FASTQ(reads_ch)

    // 3. length/quality filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 3b. read length / Q-score summary, before vs. after filtering -- a
    // side branch, nothing downstream depends on it
    if (params.enable_read_stats) {
        READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
        READ_STATS_REPORT(
            READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
            READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
        )
    }

    // 4. primer trimming (skipped internally if no primers supplied)
    CUTADAPT(CHOPPER.out.filtered)

    // 5. de novo identity clustering (decona's CD-HIT step, replaced by vsearch)
    VSEARCH_CLUSTER(CUTADAPT.out.trimmed)

    // Everything from here to BLAST_TAX runs once per sample, not once per
    // cluster: VSEARCH_CLUSTER emits all of a sample's cluster fastqs as one
    // list, and each step below loops over them inside its task. A task per
    // cluster meant thousands of Batch VMs (start-up, image pull, GCS staging)
    // each doing seconds of work. Channels are tuple(sample, [files]) throughout.

    // 6. draft consensus per cluster
    SPOA_CONSENSUS(VSEARCH_CLUSTER.out.clusters)

    // 7. alignment-based refinement
    MINIMAP2_ALIGN(SPOA_CONSENSUS.out.draft.join(VSEARCH_CLUSTER.out.clusters))
    RACON(MINIMAP2_ALIGN.out.aligned)

    // 8. ONT-specific polish -- opt-in via --enable_medaka; off by default, in
    // which case the racon output above is used as the consensus directly
    if (params.enable_medaka) {
        MEDAKA(RACON.out.polished)
        consensus_ch = MEDAKA.out.consensus
    } else {
        consensus_ch = RACON.out.polished
            .map { sample, racon_fastas, cluster_fastqs -> tuple(sample, racon_fastas) }
    }

    // 8b. second vsearch pass: fold clusters whose consensus sequences are
    // near-identical back together (one species split across several read
    // clusters by raw ONT error) -- skipped with --merge_id 0
    if (params.merge_id) {
        MERGE_CONSENSUS(consensus_ch)
        // one fasta per merged group, named for its representative cluster
        consensus_ch = MERGE_CONSENSUS.out.merged
    }

    // 9. taxonomy assignment against the freshly built BLAST db. One blastn
    // per sample over all its consensus sequences, so the db is staged once
    // per sample rather than per cluster.
    // .first() turns the (single-emission) db channels into value channels so
    // they're reused for every sample instead of being consumed after one
    BLAST_TAX(consensus_ch, MAKEBLASTDB.out.db_files.first(), MAKEBLASTDB.out.db_name.first())

    // 10. per-run abundance table + QC report. Also gathers every consensus
    // fasta into confident / low_confidence / no_hit dirs, using the same
    // best-hit call as the abundance table.
    BUILD_REPORT(
        BLAST_TAX.out.hits.map { sample, hits -> hits }.collect(),
        consensus_ch
            .flatMap { sample, fastas -> [fastas].flatten() }
            .collect(),
        taxdump
    )

    emit:
    consensus = consensus_ch
    report    = BUILD_REPORT.out.report
}
