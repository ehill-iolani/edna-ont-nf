process BLAST_TAX {
    tag "$sample $chunk"
    // process_medium, not process_low: blastn scales with threads, and a
    // several-thousand-query soil sample runs for hours -- the time limit is
    // only there to catch a hung task, since a timeout restarts from zero
    label 'process_medium'
    container "${params.container_registry}/biocontainers/blast:2.15.0--pl5321h6f7f691_1"

    input:
    tuple val(sample), val(chunk), path(consensus_fastas)
    path db_files
    val db_name

    output:
    tuple val(sample), path("${sample}.chunk${chunk}.hits.tsv"), emit: hits

    script:
    // one blastn over a chunk of a sample's consensus sequences (the workflow
    // splits them by --blast_chunk_size), not one task per cluster -- each
    // task stages the whole BLAST db (hundreds of MB for a large reference), so a single
    // sequence per task would spend far longer staging than searching. The
    // file is sample- and chunk-prefixed since these all get collected
    // together in BUILD_REPORT, which joins each sample's chunks back into
    // one published <sample>.hits.tsv.
    //
    // `awk 1` rather than cat so a fasta missing its trailing newline can't
    // run its last line into the next file's header
    """
    awk 1 ${consensus_fastas} > queries.fasta

    blastn -query queries.fasta -db ${db_name} \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs stitle" \\
        -max_target_seqs 50 -evalue 1e-10 \\
        -out blast.tsv

    # 50 targets, not 5: BUILD_REPORT ranks a cluster's hits by identity and
    # looks for near-ties among them, and BLAST's own cut to the first N
    # targets goes by score, which favours longer references. eDNA reference
    # sets mix whole mitogenomes with short partial genes, so with 5 the
    # closest record could be cut before BUILD_REPORT saw it.
    #
    # qcovs (how much of the consensus the hit covers, all HSPs together) is
    # what BUILD_REPORT's --min-qcov checks and ranks on first: a reference
    # shorter than the amplicon can only match part of it, so its identity
    # says nothing about the rest. It goes before stitle, which stays the
    # last column.
    #
    # unmatched clusters (no hit) are flagged, not dropped -- important for
    # undescribed / poorly represented Hawaiian endemic sequence variants
    #
    # qseqid here must match the id blastn used (the consensus fasta's header
    # token, "sample_clusterid") -- otherwise this row can't be joined back to
    # its sample in BUILD_REPORT. The FILENAME test (not NR == FNR) keeps this
    # right when blast.tsv is empty, i.e. when nothing at all hit.
    awk -F'\\t' '
        FILENAME == ARGV[1] { seen[\$1] = 1; next }
        /^>/ {
            id = substr(\$0, 2)
            sub(/[ \t].*/, "", id)
            if (!(id in seen)) printf "%s\\tNO_HIT\\tNA\\tNA\\tNA\\tNA\\tNA\\tflag_for_manual_review\\n", id
        }
    ' blast.tsv queries.fasta > no_hit.tsv

    cat blast.tsv no_hit.tsv > ${sample}.chunk${chunk}.hits.tsv
    """
}
