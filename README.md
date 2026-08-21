# edna-ont-nf

A Nextflow (DSL2) pipeline for classifying Oxford Nanopore amplicon eDNA data:
quality-aware clustering (isONclust) + spoa/racon/medaka consensus +
reference/BLAST taxonomy assignment. Written as a modern replacement for
[decona](https://github.com/Saskia-Oosterbroek/decona).

| decona step | this pipeline | why changed |
|---|---|---|
| qcat demux | (assumes Dorado demux upstream) | qcat unmaintained |
| (implicit) | `MERGE_FASTQ` | Dorado/MinKNOW split each barcode into many part-files; merged into one fastq per sample before filtering |
| NanoFilt | chopper | same author, actively maintained, faster |
| cutadapt | cutadapt | unchanged, still solid |
| CD-HIT clustering | isONclust | quality-aware clustering, not a blind identity cutoff |
| minimap2 + Racon | `MINIMAP2_ALIGN` + `RACON` | split into two processes -- the racon biocontainer doesn't ship minimap2 |
| Medaka | Medaka (optional, off by default) | racon consensus is used directly unless `--enable_medaka` is set |
| BLAST (optional) | BLAST (default) + explicit no-hit/low-identity flagging | undescribed/endemic sequences shouldn't silently drop |
| (pre-built db) | `MAKEBLASTDB` from `--taxdb refs.fasta` | db is (re)built from a plain reference FASTA each run, so refs can be updated without a separate `makeblastdb` step |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.10.0
- One of Docker, Singularity, or Conda (every process is containerized; pick
  a profile below)

## Quickstart

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --taxdb refs.fasta \
  --fwd_primer <seq> --rev_primer <seq> \
  -profile docker
```

`samplesheet.csv` is a two-column CSV:

```csv
sample,fastq
sample_a,data/barcode17/*.fastq.gz
sample_b,data/barcode18/*.fastq.gz
```

`fastq` may be a single file or a glob matching multiple part-files (they're
merged per sample before filtering). Glob paths resolve relative to the
launch directory, not the samplesheet's location.

`refs.fasta` is a plain FASTA of reference sequences; a BLAST db is built
from it at the start of every run.

Dev/debug entry point to run just the merge step in isolation, e.g. while
iterating on a module:

```bash
nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | *(required)* | Samplesheet CSV (`sample,fastq`) |
| `--taxdb` | *(required)* | Reference sequences FASTA; a BLAST db is built from this each run |
| `--outdir` | `results` | Output directory |
| `--fwd_primer` / `--rev_primer` | `null` | Primer sequences for cutadapt trimming; trimming is skipped if unset |
| `--min_len` / `--max_len` / `--min_qual` | `150` / `300` / `10` | chopper length/quality filtering thresholds |
| `--cluster_id` | `0.86` | isONclust similarity threshold -- tune per amplicon/primer set |
| `--min_cluster` | `20` | Minimum reads in a cluster to attempt consensus |
| `--enable_medaka` | `false` | Use medaka-polished consensus instead of the racon consensus downstream |
| `--min_pident` | `90` | BLAST hits below this %identity are flagged `low_identity`, not dropped |

All defaults live in `nextflow.config`, not `main.nf` (see the comments
there if you're adding a new one).

## Pipeline

1. `MAKEBLASTDB` -- build a BLAST db from `--taxdb` (once per run)
2. `MERGE_FASTQ` -- merge multi-part fastq(.gz) files per sample
3. `CHOPPER` -- length/quality filter
4. `CUTADAPT` -- primer trimming (skipped if no primers supplied)
5. `ISONCLUST` -- quality-aware de novo clustering
6. `SPOA_CONSENSUS` -- draft consensus per cluster
7. `MINIMAP2_ALIGN` + `RACON` -- alignment-based polish
8. `MEDAKA` -- ONT-specific polish, opt-in via `--enable_medaka`
9. `BLAST_TAX` -- taxonomy assignment against the run's BLAST db
10. `BUILD_REPORT` -- per-run abundance table + QC summary html
11. `SORT_CONSENSUS` -- gather every consensus fasta by BLAST-hit confidence

Consensus fasta headers are stamped as
`>{sample}_{cluster_id} sample={sample} cluster_size={n_reads}` by whichever
of RACON/MEDAKA produces the final consensus for a cluster.

## Output layout

```
results/
  {sample}/
    00_merged/                merged fastq
    01_filtered/               chopper output
    02_trimmed/                 cutadapt output
    03_clusters/                isONclust clusters
    04_draft/{cluster_id}/      spoa draft consensus
    05_racon/{cluster_id}/      minimap2 alignment + racon consensus
    06_consensus/               medaka consensus (only if --enable_medaka)
    07_taxonomy/{cluster_id}/   BLAST hits per cluster
  blastdb/                     BLAST db built from --taxdb
  consensus_by_confidence/
    confident/                 best BLAST hit >= --min_pident
    low_confidence/            best BLAST hit < --min_pident
    no_hit/                    no BLAST hit at all
  final_report/
    abundance_table.tsv        one row per cluster: sample, cluster_size, best hit, flag_reason
    run_qc_summary.html        cluster counts, per-sample flagged counts
  pipeline_info/                Nextflow timeline/report/trace
```

## Repo structure

```
main.nf                     entry point, samplesheet parsing, --help
workflows/edna_amplicon.nf  subworkflow chaining all steps
modules/*.nf                one process per tool, one container each
bin/build_report.py         abundance table + QC html
nextflow.config              param defaults, profiles, resource labels
conf/test.config             -profile test overrides (small synthetic dataset)
assets/                      your own local samplesheet/taxdb go here (gitignored)
tests/data/                  small synthetic dataset used by -profile test
.github/workflows/ci.yml     runs -profile test on push/PR
```

## Testing

`tests/data/` holds a small synthetic dataset (two samples, a handful of
reads each) checked into the repo -- see `conf/test.config` for the exact
params it's run with. Use it to sanity-check the pipeline end-to-end without
needing your own data:

```bash
nextflow run main.nf -profile test,docker
```

This populates `results/` exactly like a normal run (see Output layout
above), just against the small checked-in fixtures instead of your own data
(eyeball the abundance table, poke around `consensus_by_confidence/`, etc).

Note: the `medaka` container is ONT's own multi-arch image
(`ontresearch/medaka`), not biocontainers. The biocontainers build is
amd64-only and SIGILLs under Docker's emulation on Apple Silicon.
