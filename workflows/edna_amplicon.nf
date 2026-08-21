include { MAKEBLASTDB      } from '../modules/makeblastdb.nf'
include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER          } from '../modules/chopper.nf'
include { CUTADAPT         } from '../modules/cutadapt.nf'
include { ISONCLUST        } from '../modules/isonclust.nf'
include { SPOA_CONSENSUS   } from '../modules/spoa_consensus.nf'
include { MINIMAP2_ALIGN   } from '../modules/minimap2_align.nf'
include { RACON            } from '../modules/racon.nf'
include { MEDAKA           } from '../modules/medaka.nf'
include { BLAST_TAX        } from '../modules/blast_tax.nf'
include { SORT_CONSENSUS   } from '../modules/sort_consensus.nf'
include { BUILD_REPORT     } from '../modules/report.nf'

workflow EDNA_AMPLICON {

    take:
    reads_ch    // tuple(sample, fastq)
    ref_fasta   // path to reference sequences fasta for BLAST taxonomy

    main:
    // 1. build a BLAST db from the reference fasta (once per run, independent
    //    of per-sample steps below)
    MAKEBLASTDB(ref_fasta)

    // 2. merge multi-part fastq(.gz) files per barcode into one file per sample
    MERGE_FASTQ(reads_ch)

    // 3. length/quality filter
    CHOPPER(MERGE_FASTQ.out.merged)

    // 4. primer trimming (skipped internally if no primers supplied)
    CUTADAPT(CHOPPER.out.filtered)

    // 5. quality-aware de novo clustering (decona's CD-HIT step, replaced)
    ISONCLUST(CUTADAPT.out.trimmed)

    // ISONCLUST emits one fastq per cluster per sample; flatten and tag
    clusters_ch = ISONCLUST.out.clusters
        .flatMap { sample, cluster_fastqs ->
            cluster_fastqs.collect { fq -> tuple(sample, fq.baseName, fq) }
        }

    // 6. draft consensus per cluster
    SPOA_CONSENSUS(clusters_ch)

    // 7. alignment-based refinement
    MINIMAP2_ALIGN(SPOA_CONSENSUS.out.draft.join(clusters_ch, by: [0, 1]))
    RACON(MINIMAP2_ALIGN.out.aligned)

    // 8. ONT-specific polish -- opt-in via --enable_medaka; off by default, in
    // which case the racon output above is used as the consensus directly
    if (params.enable_medaka) {
        MEDAKA(RACON.out.polished)
        consensus_ch = MEDAKA.out.consensus
    } else {
        consensus_ch = RACON.out.polished
            .map { sample, cluster_id, racon_fasta, cluster_fastq -> tuple(sample, cluster_id, racon_fasta) }
    }

    // 9. taxonomy assignment against the freshly built BLAST db
    // .first() turns the (single-emission) db channels into value channels so
    // they're reused for every cluster instead of being consumed after one
    BLAST_TAX(consensus_ch, MAKEBLASTDB.out.db_files.first(), MAKEBLASTDB.out.db_name.first())

    // 10. per-run abundance table + QC report
    BUILD_REPORT(
        BLAST_TAX.out.hits.map { sample, cluster_id, hits -> hits }.collect(),
        consensus_ch.map { sample, cluster_id, fasta -> fasta }.collect()
    )

    // 11. gather every consensus fasta into confident / low_confidence / no_hit
    // dirs, using the same best-hit classification as the abundance table
    SORT_CONSENSUS(consensus_ch.join(BLAST_TAX.out.hits, by: [0, 1]))

    emit:
    consensus = consensus_ch
    report    = BUILD_REPORT.out.report
}
