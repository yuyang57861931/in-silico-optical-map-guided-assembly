#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reference-guided pseudomolecule builder with contig splitting.

For each contig:
1) Choose the primary alignment (largest query alignment length).
2) Split the contig into Left / Primary / Right by query coordinates.
3) For Left/Right segments:
   - If segment length < min_unaligned_keep, keep with primary ref.
   - Else, if coverage on another ref >= move_ratio, move to that ref;
     otherwise keep with primary ref.

Outputs:
  - pseudomolecules.fasta: concatenated per-ref sequences with N gaps.
  - contig_split.tsv: per-piece coordinates, length, assigned ref, reason, etc.
  - assembly_stats.tsv: summary statistics.

Example:
  python 02.new-build_pseudomolecules.py \
    -x alignments.xmap \
    -c contigs.fa \
    -o pseudomolecules.fasta \
    --gap-size 1000 \
    --move-ratio 0.4 \
    --min-unaligned-keep 100000 \
    --split-tsv contig_split.tsv \
    --stats-tsv assembly_stats.tsv \
    --write-unplaced \
    --write-unaligned --unaligned-fasta unaligned_contigs.fasta
"""

import sys
import gzip
import argparse
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ------------------ utils ------------------

def open_text_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def natural_sort_refs(keys):
    def key_func(k):
        s = str(k)
        for p in ("chr", "CHR", "Chr"):
            if s.startswith(p):
                s = s[len(p):]
                break
        try:
            return (0, int(s))
        except ValueError:
            return (1, str(k))
    return sorted(keys, key=key_func)


# ------------------ IO ------------------

def load_contigs(fa: str) -> dict:
    print("Loading contig sequences ...")
    with open_text_maybe_gzip(fa, "rt") as h:
        contigs = SeqIO.to_dict(SeqIO.parse(h, "fasta"))
    if not contigs:
        raise ValueError("No sequences found in FASTA.")
    print(f"Loaded {len(contigs)} contigs.")
    return contigs


def load_xmap(xmap_path: str) -> pd.DataFrame:
    """
    Load XMAP and keep: qry, ref, qry_start, qry_end, ref_start, ref_end, strand.
    """
    print("Parsing XMAP file ...")
    df = pd.read_csv(
        open_text_maybe_gzip(xmap_path, "rt"),
        sep="\t",
        comment="#",
        header=None,
        usecols=[1, 2, 3, 4, 5, 6, 7],
        names=["qry", "ref", "qry_start", "qry_end", "ref_start", "ref_end", "strand"],
        dtype=str,
        engine="python"
    )
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    # keep contig names as-is
    for c in ["qry_start", "qry_end", "ref_start", "ref_end"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["qry", "ref", "qry_start", "qry_end", "ref_start", "ref_end"])

    df["strand"] = df["strand"].map(lambda s: "-" if str(s).startswith("-") else "+")
    m = df["qry_start"] > df["qry_end"]
    df.loc[m, ["qry_start", "qry_end"]] = df.loc[m, ["qry_end", "qry_start"]].values
    m2 = df["ref_start"] > df["ref_end"]
    df.loc[m2, ["ref_start", "ref_end"]] = df.loc[m2, ["ref_end", "ref_start"]].values

    # alignment length on query (inclusive)
    df["align_len"] = df["qry_end"] - df["qry_start"] + 1
    print(f"Valid alignments in XMAP: {len(df)}")
    return df


# ------------------ splitting logic ------------------

def choose_primary_hit(group: pd.DataFrame) -> pd.Series:
    """Choose primary hit by maximum align_len."""
    return group.loc[group["align_len"].idxmax()]


def coverage_by_ref(piece_start: int, piece_end: int, hits: pd.DataFrame):
    """
    Compute coverage of [piece_start, piece_end] on each ref.
    Return dict {ref: (covered_len, best_hit_row)}.
    """
    cov = {}
    L = piece_end - piece_start + 1
    if hits.empty or L <= 0:
        return cov

    for ref_id, sub in hits.groupby("ref"):
        covered = 0
        best_len = -1
        best_row = None
        for _, r in sub.iterrows():
            s = max(piece_start, int(r["qry_start"]))
            e = min(piece_end, int(r["qry_end"]))
            olap = e - s + 1
            if olap > 0:
                covered += olap
                if olap > best_len:
                    best_len = olap
                    best_row = r
        if covered > 0:
            cov[ref_id] = (covered, best_row)
    return cov


def split_contig_by_rule(contig_id: str,
                         contig_len: int,
                         all_hits_one: pd.DataFrame,
                         move_ratio: float = 0.4,
                         min_unaligned_keep: int = 500_000):
    """
    Split one contig into pieces and assign each piece to a reference.
    """
    out = []

    if all_hits_one.empty:
        out.append(dict(
            contig=contig_id, start=1, end=contig_len, piece_len=contig_len,
            assigned_ref="unplaced", strand="+", reason="unaligned_all",
            support_ref_start=pd.NA, support_ref_end=pd.NA,
            support_cov_len=0, support_cov_frac=0.0,
            is_primary=0, order_key=float("inf")
        ))
        return out

    primary = choose_primary_hit(all_hits_one)
    p_qs, p_qe = int(primary["qry_start"]), int(primary["qry_end"])
    primary_ref = primary["ref"]
    p_ref_start = float(primary["ref_start"])
    p_strand = primary["strand"]
    if p_qs > p_qe:
        p_qs, p_qe = p_qe, p_qs

    out.append(dict(
        contig=contig_id, start=p_qs, end=p_qe, piece_len=p_qe - p_qs + 1,
        assigned_ref=primary_ref, strand=p_strand, reason="primary",
        support_ref_start=primary["ref_start"], support_ref_end=primary["ref_end"],
        support_cov_len=p_qe - p_qs + 1, support_cov_frac=1.0,
        is_primary=1, order_key=p_ref_start
    ))

    others = all_hits_one[all_hits_one.index != primary.name].copy()

    left_piece = (1, p_qs - 1) if p_qs > 1 else None
    right_piece = (p_qe + 1, contig_len) if p_qe < contig_len else None

    def keep_with_primary(s, e, L, side_sign, reason):
        out.append(dict(
            contig=contig_id, start=s, end=e, piece_len=L,
            assigned_ref=primary_ref, strand=p_strand, reason=reason,
            support_ref_start=pd.NA, support_ref_end=pd.NA,
            support_cov_len=0, support_cov_frac=0.0,
            is_primary=0, order_key=p_ref_start + 0.5 * side_sign
        ))

    def handle_side(piece, side_sign):
        if not piece:
            return
        s, e = piece
        if e < s:
            return
        L = e - s + 1

        # short piece: keep with primary
        if L < min_unaligned_keep:
            keep_with_primary(s, e, L, side_sign,
                              f"kept_with_primary(small<{min_unaligned_keep})")
            return

        # compute alternative coverage
        cov = coverage_by_ref(s, e, others)
        if not cov:
            keep_with_primary(s, e, L, side_sign, "kept_with_primary(<ratio, no_alt_hit)")
            return

        best_ref, (covered, best_row) = max(cov.items(), key=lambda kv: kv[1][0])
        frac = covered / L

        if frac >= move_ratio:
            out.append(dict(
                contig=contig_id, start=s, end=e, piece_len=L,
                assigned_ref=best_ref, strand=best_row["strand"],
                reason=f"moved(>={int(move_ratio * 100)}%)",
                support_ref_start=best_row["ref_start"],
                support_ref_end=best_row["ref_end"],
                support_cov_len=int(covered),
                support_cov_frac=round(frac, 4),
                is_primary=0,
                order_key=float(best_row["ref_start"])
            ))
        else:
            keep_with_primary(s, e, L, side_sign,
                              f"kept_with_primary(<{int(move_ratio * 100)}%)")

    handle_side(left_piece, -1)
    handle_side(right_piece, +1)
    return out


def split_all_contigs(xmap_df: pd.DataFrame,
                      contigs: dict,
                      move_ratio: float = 0.4,
                      min_unaligned_keep: int = 500_000):
    """
    Split all contigs: first those in XMAP, then add fully unaligned contigs.
    """
    records = []
    xmap_df = xmap_df.sort_values(["qry", "align_len"], ascending=[True, False]).reset_index(drop=True)

    seen = set()
    for contig_id, grp in xmap_df.groupby("qry"):
        if contig_id not in contigs:
            continue
        seen.add(contig_id)
        clen = len(contigs[contig_id])
        pieces = split_contig_by_rule(contig_id, clen, grp,
                                      move_ratio=move_ratio,
                                      min_unaligned_keep=min_unaligned_keep)
        records.extend(pieces)

    all_ids = set(contigs.keys())
    unaligned = sorted(all_ids - seen)
    for cid in unaligned:
        clen = len(contigs[cid])
        records.append(dict(
            contig=cid, start=1, end=clen, piece_len=clen,
            assigned_ref="unplaced", strand="+", reason="unaligned_all",
            support_ref_start=pd.NA, support_ref_end=pd.NA,
            support_cov_len=0, support_cov_frac=0.0, is_primary=0, order_key=float("inf")
        ))

    split_df = pd.DataFrame.from_records(records, columns=[
        "contig", "start", "end", "piece_len", "assigned_ref", "strand", "reason",
        "support_ref_start", "support_ref_end", "support_cov_len", "support_cov_frac", "is_primary", "order_key"
    ])
    return split_df


# ------------------ build pseudomolecules ------------------

def build_pseudomolecules(split_df: pd.DataFrame,
                          contigs: dict,
                          gap_size: int):
    grouped = defaultdict(list)  # ref -> list of (row, seq)
    used_contigs = set()

    for row in split_df.itertuples(index=False):
        if row.assigned_ref == "unplaced":
            continue
        cid = row.contig
        s, e = int(row.start), int(row.end)
        seq = contigs[cid].seq[s - 1:e]
        if row.strand == "-":
            seq = seq.reverse_complement()
        grouped[row.assigned_ref].append((row, str(seq)))
        used_contigs.add(cid)

    gap = "N" * gap_size
    records = []
    for ref in natural_sort_refs(grouped.keys()):
        items = grouped[ref]
        items.sort(key=lambda x: (
            float(x[0].order_key),
            -int(x[0].is_primary),
            str(x[0].contig),
            int(x[0].start),
        ))
        seqs = [seq for _, seq in items]
        joined = gap.join(seqs)
        records.append(SeqRecord(Seq(joined), id=f"chr{ref}", description=f"pseudomolecule;gap={gap_size}N"))
        print(f"Built chr{ref} with {len(items)} pieces.")
    return records, used_contigs, sum(len(v) for v in grouped.values())


# ------------------ stats ------------------

def write_split_tsv(split_df: pd.DataFrame, path: str):
    cols = [
        "contig", "start", "end", "piece_len", "assigned_ref", "strand",
        "reason", "support_ref_start", "support_ref_end",
        "support_cov_len", "support_cov_frac", "is_primary",
    ]
    split_df[cols].to_csv(path, sep="\t", index=False)
    print(f"Wrote contig split table: {path}")


def write_stats(records,
                contigs,
                used_contigs,
                piece_count,
                split_df,
                out_path: str):
    total_len = sum(len(r.seq) for r in records)
    stats = {
        "total_chromosomes": len(records),
        "total_length": total_len,
        "avg_chrom_length": (total_len / len(records)) if records else 0,
        "total_contigs": len(contigs),
        "used_contigs_unique": len(used_contigs),
        "piece_count": piece_count,
        "primary_piece_count": int((split_df["reason"] == "primary").sum()),
        "moved_piece_count": int(split_df["reason"].str.startswith("moved").sum()),
        "kept_small_piece_count": int(split_df["reason"].str.startswith("kept_with_primary(small").sum()),
        "kept_ratio_piece_count": int(split_df["reason"].str.startswith("kept_with_primary(<").sum()),
        "unaligned_all_contigs": int((split_df["reason"] == "unaligned_all").sum()),
    }
    pd.DataFrame([stats]).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote summary stats: {out_path}")


# ------------------ main ------------------

def main(args):
    contigs = load_contigs(args.contigs_fasta)
    xmap_df = load_xmap(args.xmap_file)

    split_df = split_all_contigs(
        xmap_df, contigs,
        move_ratio=args.move_ratio,
        min_unaligned_keep=args.min_unaligned_keep,
    )
    write_split_tsv(split_df, args.split_tsv)

    records, used_contigs, piece_count = build_pseudomolecules(
        split_df, contigs, gap_size=args.gap_size
    )
    SeqIO.write(records, args.output_fasta, "fasta")
    print(f"Wrote pseudomolecule FASTA: {args.output_fasta}")

    if args.write_unaligned:
        unaligned_ids = sorted(split_df.loc[split_df["reason"] == "unaligned_all", "contig"].unique())
        if unaligned_ids:
            with open(args.unaligned_fasta, "w") as out:
                for cid in unaligned_ids:
                    SeqIO.write(contigs[cid], out, "fasta")
            print(f"Wrote unaligned contigs ({len(unaligned_ids)}): {args.unaligned_fasta}")
        else:
            print("No fully unaligned contigs found.")

    all_ids = set(contigs.keys())
    used_ids_from_split = set(split_df["contig"].unique())
    unused = all_ids - used_ids_from_split
    if args.write_unplaced and unused:
        with open(args.unplaced_fasta, "w") as out:
            for cid in sorted(unused):
                SeqIO.write(contigs[cid], out, "fasta")
        print(f"Wrote unused contigs ({len(unused)}): {args.unplaced_fasta}")

    write_stats(records, contigs, used_contigs, piece_count, split_df, args.stats_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build pseudomolecules with contig splitting.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-x", "--xmap-file", required=True, help=".xmap or .xmap.gz file")
    parser.add_argument("-c", "--contigs-fasta", required=True, help="contigs FASTA (optionally gzipped)")
    parser.add_argument("-o", "--output-fasta", default="pseudomolecules.fasta", help="output pseudomolecule FASTA")
    parser.add_argument("--gap-size", type=int, default=1000, help="number of Ns between adjacent pieces")
    parser.add_argument("--move-ratio", type=float, default=0.4, help="coverage ratio threshold to move a piece")
    parser.add_argument("--min-unaligned-keep", type=int, default=500000, help="unaligned piece length below this is kept with primary ref")
    parser.add_argument("--write-unplaced", action="store_true", help="write contigs that are never used")
    parser.add_argument("--unplaced-fasta", default="unplaced_contigs.fasta", help="output FASTA for never-used contigs")
    parser.add_argument("--write-unaligned", action="store_true", help="write contigs with reason=unaligned_all")
    parser.add_argument("--unaligned-fasta", default="unaligned_contigs.fasta", help="output FASTA for fully unaligned contigs")
    parser.add_argument("--split-tsv", default="contig_split.tsv", help="TSV file of contig pieces and placement")
    parser.add_argument("--stats-tsv", default="assembly_stats.tsv", help="TSV file of assembly statistics")
    args = parser.parse_args()

    try:
        main(args)
    except Exception as e:
        print(f"Error during execution: {e}", file=sys.stderr)
        sys.exit(1)