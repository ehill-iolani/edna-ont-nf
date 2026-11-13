#!/usr/bin/env python3
"""
Aggregate per-cluster BLAST hits + consensus fastas into a single
sample x taxon abundance table, plus a minimal QC/summary HTML.

This is a functional starting point, not the final report -- extend with
per-cluster read counts (from the vsearch clusters.uc), top-hit
filtering by pident/evalue, and a proper MultiQC-style layout once the
sample sheet format and reference taxonomy fields are finalized.
"""
import argparse
import os
import shutil
import tarfile

import pandas as pd

COLS = ["seq_id", "subject_id", "pident", "length", "evalue", "bitscore", "qcovs", "stitle"]
SKIPPED_COLS = ["sample", "trimmed_reads", "largest_cluster", "reason"]


def load_hits(paths):
    frames = []
    for p in paths:
        try:
            df = pd.read_csv(p, sep="\t", names=COLS)
            df["source_file"] = p
            frames.append(df)
        except pd.errors.EmptyDataError:
            continue
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=COLS)


def species_label(stitle):
    # same rule as the frontend's shortTaxonLabel: BLAST's stitle starts with
    # the binomial, so the first two tokens name the species
    toks = stitle.split() if isinstance(stitle, str) else []
    return " ".join(toks[:2])


def rank_hits(hits_df, min_qcov):
    """A cluster's hits, best first: those covering at least min_qcov % of the
    consensus ahead of the rest, ranked by identity, then bitscore.

    Not by bitscore alone: bitscore grows with alignment length, and eDNA
    reference sets mix whole mitogenomes with short partial genes, so a longer
    but more distant reference outscored a closer, shorter one (a Phycis
    blennoides consensus went to a 98.8% P. chesteri mitogenome over 1,930 bp
    instead of the 99.9% P. blennoides record over 1,504 bp). Same rule as
    ITS-nf. Hits under the coverage floor keep their bitscore order, so a
    cluster with none above it still reports its longest hit (and is flagged
    low_coverage). The sort is stable, so BLAST's own order breaks remaining
    ties."""
    hits = hits_df.copy()
    qcov = pd.to_numeric(hits["qcovs"], errors="coerce")
    hits["_covered"] = qcov >= min_qcov
    hits["_identity"] = hits["pident"].where(hits["_covered"])
    return hits.sort_values(["_covered", "_identity", "bitscore"], ascending=False,
                            kind="stable", na_position="last")


def tied_species(ranked, tie_pident):
    """seq_id -> "|"-joined distinct species among the hits within tie_pident
    percentage points of the cluster's best identity (covering hits only), or
    "" when the best hit is unambiguous.

    Short markers often can't tell close species apart (12S across many
    flatfish and ray congeners), and the reference then holds several species
    within a fraction of a percent of the consensus. Naming the first of them
    would be a coin flip, so they're recorded here and the frontend reports
    the cluster at genus level (same genus) or as ambiguous (different
    genera) instead."""
    out = {}
    hits = ranked[(ranked["subject_id"] != "NO_HIT") & ranked["_covered"]]
    for seq_id, grp in hits.groupby("seq_id", sort=False):
        near = grp[grp["pident"] >= grp["pident"].max() - tie_pident]
        names = list(dict.fromkeys(n for n in map(species_label, near["stitle"]) if n))
        out[seq_id] = "|".join(names) if len(names) > 1 else ""
    return out


def load_skipped(paths):
    """VSEARCH_CLUSTER's one-line <sample>.skipped.tsv files (a sample with no
    reads after trimming, or no cluster of --min_cluster reads)."""
    rows = []
    for p in paths or []:
        with open(p) as fh:
            rows += [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return pd.DataFrame(rows, columns=SKIPPED_COLS).sort_values("sample", kind="stable")


def load_consensus_meta(paths):
    # consensus fasta headers are stamped by RACON/MEDAKA as:
    #   >{sample}_{cluster_id} sample={sample} cluster_size={n_reads}
    rows = []
    for p in paths:
        with open(p) as fh:
            header = fh.readline().strip()
        if not header.startswith(">"):
            continue
        tokens = header[1:].split()
        seq_id = tokens[0]
        fields = dict(tok.split("=", 1) for tok in tokens[1:] if "=" in tok)
        rows.append({
            "path": p,
            "seq_id": seq_id,
            "sample": fields.get("sample"),
            "cluster_size": int(fields["cluster_size"]) if "cluster_size" in fields else None,
        })
    return pd.DataFrame(rows, columns=["path", "seq_id", "sample", "cluster_size"])


# flag_reason -> the consensus_by_confidence/ folder a cluster's fasta goes in
CONFIDENCE_DIR = {"": "confident", "low_identity": "low_confidence", "low_coverage": "low_confidence", "no_hit": "no_hit"}


def sort_consensus(best, out_dir):
    """Copy every consensus fasta into out_dir/{confident,low_confidence,no_hit}/
    as {sample}.{cluster_id}.consensus.fasta, by the same call the abundance
    table makes (flag_reason), so the two can't disagree."""
    for d in CONFIDENCE_DIR.values():
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    for row in best.dropna(subset=["path", "sample"]).itertuples():
        cluster_id = row.seq_id[len(row.sample) + 1:]  # seq_id is {sample}_{cluster_id}
        dest = os.path.join(out_dir, CONFIDENCE_DIR[row.flag_reason], f"{row.sample}.{cluster_id}.consensus.fasta")
        shutil.copyfile(row.path, dest)


RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def _open_dmp(taxdump, name):
    # taxdump is either the extracted NCBI taxdump directory or taxdump.tar.gz
    if os.path.isdir(taxdump):
        return open(os.path.join(taxdump, name))
    tf = tarfile.open(taxdump, "r:gz")
    return (line.decode() for line in tf.extractfile(name))


def candidate_names(stitle):
    # reference headers are "ACC Genus species ...", and BLAST's stitle drops
    # the accession, so the binomial is the first two tokens; fall back to the
    # genus alone (e.g. "Genus sp.") so the lineage is at least filled to genus
    if not isinstance(stitle, str):
        return []
    toks = stitle.split()
    return [" ".join(toks[:2]), toks[0]] if len(toks) >= 2 else toks[:1]


def add_lineage(best, taxdump):
    """Append taxid + one column per RANKS entry, resolved from the NCBI taxdump."""
    wanted = {n for s in best["stitle"] for n in candidate_names(s)}

    # a name can map to >1 taxid (homonyms across kingdoms); take the lowest
    # taxid, preferring a scientific name over a synonym
    name_to_taxid = {}
    for line in _open_dmp(taxdump, "names.dmp"):
        taxid, name, _, cls = [f.strip() for f in line.rstrip("\n").split("|")][:4]
        if name in wanted and cls in ("scientific name", "synonym"):
            key = (cls != "scientific name", int(taxid))
            if name not in name_to_taxid or key < name_to_taxid[name][0]:
                name_to_taxid[name] = (key, int(taxid))

    parent, rank = {}, {}
    for line in _open_dmp(taxdump, "nodes.dmp"):
        f = [x.strip() for x in line.split("|")]
        parent[int(f[0])], rank[int(f[0])] = int(f[1]), f[2]

    sci = {}
    need = set()

    def lineage(taxid):
        out = {}
        while True:
            if rank.get(taxid) in RANKS:
                out[rank[taxid]] = taxid
            if parent.get(taxid, taxid) == taxid:
                return out
            taxid = parent[taxid]

    rows = []
    for stitle in best["stitle"]:
        taxid = next((name_to_taxid[n][1] for n in candidate_names(stitle) if n in name_to_taxid), None)
        lin = lineage(taxid) if taxid else {}
        need.update(lin.values())
        rows.append((taxid, lin))

    for line in _open_dmp(taxdump, "names.dmp"):
        f = [x.strip() for x in line.rstrip("\n").split("|")]
        if f[3] == "scientific name" and int(f[0]) in need:
            sci[int(f[0])] = f[1]

    best = best.copy()
    best["taxid"] = [t if t else "" for t, _ in rows]
    for r in RANKS:
        best[r] = [sci.get(lin.get(r), "") for _, lin in rows]
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", nargs="+", required=True)
    ap.add_argument("--consensus", nargs="+", required=True)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-sorted", help="directory to gather consensus fastas into, by confidence")
    ap.add_argument("--taxdump", default=None,
                     help="NCBI taxdump dir or taxdump.tar.gz; adds taxid + full lineage columns")
    ap.add_argument("--min-pident", type=float, default=90,
                     help="hits below this identity are flagged low_identity, not dropped")
    ap.add_argument("--min-qcov", type=float, default=0,
                     help="hits covering less than this %% of the consensus (BLAST qcovs) are flagged "
                          "low_coverage, not dropped; 0 = off")
    ap.add_argument("--tie-pident", type=float, default=1.0,
                     help="hits within this many percentage points of a cluster's best identity "
                          "count as tied with it (see tied_species)")
    ap.add_argument("--min-abundance", type=int, default=0,
                     help="merged clusters with fewer reads than this get low_abundance=true, not dropped")
    ap.add_argument("--min-rel-abundance", type=float, default=0,
                     help="merged clusters holding a smaller fraction of their sample's clustered reads "
                          "than this get low_abundance=true, not dropped")
    ap.add_argument("--skipped", nargs="*", default=[],
                     help="VSEARCH_CLUSTER's <sample>.skipped.tsv files; listed in skipped_samples.tsv")
    ap.add_argument("--out-skipped", default="skipped_samples.tsv",
                     help="written only when --skipped names at least one sample")
    args = ap.parse_args()

    hits_df = load_hits(args.hits)
    consensus_meta = load_consensus_meta(args.consensus)

    # keep the best-ranked hit per cluster as the working call (rank_hits:
    # covering hits first, then identity, then bitscore). Ties are common
    # (several reference records with identical scores), so the sort must be
    # stable: BLAST lists a cluster's hits best-first, and a stable sort keeps
    # that order among ties. pandas' default sort is not stable and picked an
    # arbitrary tied hit, which gave clusters with identical hits different
    # species.
    ranked = rank_hits(hits_df, args.min_qcov)
    # drop_duplicates rather than groupby().first(), which takes each column's
    # first non-null value and could mix fields from different hits
    best = ranked.drop_duplicates("seq_id")[COLS].sort_values("seq_id").reset_index(drop=True)
    best = best.merge(consensus_meta, on="seq_id", how="left")
    best["tied_taxa"] = best["seq_id"].map(tied_species(ranked, args.tie_pident)).fillna("")

    # clusters are flagged, not dropped: for no BLAST hit at all, for a best
    # hit that covers too little of the consensus to be a match for the whole
    # amplicon (qcovs < min-qcov -- typically a reference shorter than the
    # amplicon, or a hit on a conserved stretch only), or for one too
    # divergent to call with confidence (pident < min-pident). Coverage is
    # checked first: a short hit's identity says nothing about the rest of
    # the sequence.
    best["flag_reason"] = ""
    best.loc[best["subject_id"] == "NO_HIT", "flag_reason"] = "no_hit"
    qcovs = pd.to_numeric(best["qcovs"], errors="coerce")
    low_qcov = qcovs.notna() & (qcovs < args.min_qcov)
    best.loc[low_qcov & (best["flag_reason"] == ""), "flag_reason"] = "low_coverage"
    low_pident = best["pident"].notna() & (best["pident"] < args.min_pident)
    best.loc[low_pident & (best["flag_reason"] == ""), "flag_reason"] = "low_identity"

    # Thin support is a separate axis from the BLAST call above, so it gets its
    # own true/false column rather than another flag_reason value: flag_reason
    # stays exactly '' | 'no_hit' | 'low_coverage' | 'low_identity', which downstream readers
    # (e.g. the web frontend) take as the whole story and use as a taxon label.
    # Judged on the *merged* cluster size (reads from every cluster folded into
    # this one), so a species split into several small read clusters isn't
    # penalised for the split. The relative form is a fraction of the sample's
    # clustered reads, which keeps the cutoff meaningful across barcodes of
    # very different depth.
    sample_total = best.groupby("sample")["cluster_size"].transform("sum")
    is_low_abundance = best["cluster_size"].notna() & (
        (best["cluster_size"] < args.min_abundance)
        | (best["cluster_size"] / sample_total < args.min_rel_abundance)
    )
    best["low_abundance"] = is_low_abundance.map({True: "true", False: "false"})

    tax_cols = []
    if args.taxdump:
        best = add_lineage(best, args.taxdump)
        tax_cols = ["taxid"] + RANKS

    if args.out_sorted:
        sort_consensus(best, args.out_sorted)

    best = best[["seq_id", "sample", "cluster_size", "subject_id", "pident",
                 "length", "evalue", "bitscore", "qcovs", "stitle", "tied_taxa"] + tax_cols + ["flag_reason", "low_abundance"]]
    best.to_csv(args.out_table, sep="\t", index=False)

    # "flagged" is the BLAST-confidence flag; low_abundance is counted
    # separately since it says how much support a call has, not whether the
    # call itself is doubtful
    n_clusters = best.shape[0]
    n_flagged = (best["flag_reason"] != "").sum()
    n_low_abundance = (best["low_abundance"] == "true").sum()
    per_sample = (
        best.groupby("sample")
        .agg(clusters=("seq_id", "count"),
             flagged=("flag_reason", lambda s: (s != "").sum()),
             low_abundance=("low_abundance", lambda s: (s == "true").sum()))
        .reset_index()
    )

    # samples VSEARCH_CLUSTER had nothing to cluster for: in no table above,
    # so listed on their own (a file only when there are any, so a run without
    # them has exactly its old outputs)
    skipped = load_skipped(args.skipped)
    if not skipped.empty:
        skipped.to_csv(args.out_skipped, sep="\t", index=False)

    with open(args.out_html, "w") as fh:
        fh.write("<html><body><h2>Run summary</h2>")
        fh.write(f"<p>Clusters processed: {n_clusters}</p>")
        fh.write(f"<p>Clusters flagged for manual review (no hit, query coverage &lt; {args.min_qcov:g}% "
                 f"or pident &lt; {args.min_pident}): {n_flagged}</p>")
        cutoffs = []
        if args.min_abundance:
            cutoffs.append(f"&lt; {args.min_abundance} reads")
        if args.min_rel_abundance:
            cutoffs.append(f"&lt; {args.min_rel_abundance:g} of the sample's clustered reads")
        if cutoffs:
            fh.write(f"<p>Clusters marked low_abundance ({' or '.join(cutoffs)}): {n_low_abundance}</p>")
        fh.write("<h3>Per-sample</h3>")
        fh.write(per_sample.to_html(index=False))
        if not skipped.empty:
            fh.write("<h3>Skipped samples</h3>")
            fh.write("<p>Nothing to identify: no reads left after primer trimming (no_reads), or no "
                     "cluster reached --min_cluster reads (no_cluster). Expected for blanks and "
                     "negative controls.</p>")
            fh.write(skipped.to_html(index=False))
        fh.write("</body></html>")


if __name__ == "__main__":
    main()
