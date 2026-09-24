#!/usr/bin/env python3
"""
Aggregate per-cluster BLAST hits + consensus fastas into a single
sample x taxon abundance table, plus a minimal QC/summary HTML.

This is a functional starting point, not the final report -- extend with
per-cluster read counts (from isONclust final_clusters.tsv), top-hit
filtering by pident/evalue, and a proper MultiQC-style layout once the
sample sheet format and reference taxonomy fields are finalized.
"""
import argparse
import os
import tarfile

import pandas as pd

COLS = ["seq_id", "subject_id", "pident", "length", "evalue", "bitscore", "stitle"]


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
            "seq_id": seq_id,
            "sample": fields.get("sample"),
            "cluster_size": int(fields["cluster_size"]) if "cluster_size" in fields else None,
        })
    return pd.DataFrame(rows, columns=["seq_id", "sample", "cluster_size"])


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
    ap.add_argument("--taxdump", default=None,
                     help="NCBI taxdump dir or taxdump.tar.gz; adds taxid + full lineage columns")
    ap.add_argument("--min-pident", type=float, default=90,
                     help="hits below this identity are flagged low_identity, not dropped")
    args = ap.parse_args()

    hits_df = load_hits(args.hits)
    consensus_meta = load_consensus_meta(args.consensus)

    # keep best hit per cluster (highest bitscore) as the working call.
    # Ties are common (several reference records with identical bitscores),
    # so the sort must be stable: BLAST lists a cluster's hits best-first, and
    # a stable sort keeps that order among tied bitscores, so the hit picked
    # here is the first line of the cluster's hits file -- the same one
    # SORT_CONSENSUS reads to choose the confidence folder. pandas' default
    # sort is not stable and picked an arbitrary tied hit, which could disagree
    # with SORT_CONSENSUS (a cluster flagged low_identity but filed under
    # confident/) and gave clusters with identical hits different species.
    best = (
        hits_df.sort_values("bitscore", ascending=False, kind="stable")
        .groupby("seq_id", as_index=False)
        .first()
    )
    best = best.merge(consensus_meta, on="seq_id", how="left")

    # clusters are flagged, not dropped, either for no BLAST hit at all or for
    # a best hit too divergent to call with confidence (pident < min-pident)
    best["flag_reason"] = ""
    best.loc[best["subject_id"] == "NO_HIT", "flag_reason"] = "no_hit"
    low_pident = best["pident"].notna() & (best["pident"] < args.min_pident)
    best.loc[low_pident & (best["flag_reason"] == ""), "flag_reason"] = "low_identity"

    tax_cols = []
    if args.taxdump:
        best = add_lineage(best, args.taxdump)
        tax_cols = ["taxid"] + RANKS

    best = best[["seq_id", "sample", "cluster_size", "subject_id", "pident",
                 "length", "evalue", "bitscore", "stitle"] + tax_cols + ["flag_reason"]]
    best.to_csv(args.out_table, sep="\t", index=False)

    n_clusters = best.shape[0]
    n_flagged = (best["flag_reason"] != "").sum()
    per_sample = (
        best.groupby("sample")
        .agg(clusters=("seq_id", "count"), flagged=("flag_reason", lambda s: (s != "").sum()))
        .reset_index()
    )

    with open(args.out_html, "w") as fh:
        fh.write("<html><body><h2>Run summary</h2>")
        fh.write(f"<p>Clusters processed: {n_clusters}</p>")
        fh.write(f"<p>Clusters flagged for manual review (no hit or pident &lt; {args.min_pident}): {n_flagged}</p>")
        fh.write("<h3>Per-sample</h3>")
        fh.write(per_sample.to_html(index=False))
        fh.write("</body></html>")


if __name__ == "__main__":
    main()
