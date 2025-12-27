#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reference-guided pseudomolecule builder without contig splitting.
Each contig is anchored as a whole using its primary XMAP alignment.
Example:
python 22.new-build_pseudomolecules.no_split.py \
  -x 21.20.lm15_11X_2_CS_IAAS_lo7.xmap.processed \
  -c ../12.split_genome.fa.gz \
  -o 22.lm15_11x_pseudomolecules.fasta \
  --gap-size 1000 \
  --assign-tsv contig_assign.tsv \
  --stats-tsv assembly_stats.tsv \
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


def open_text_maybe_gzip(path, mode="rt"):
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


def load_contigs(fa):
    print("Loading contig sequences...")
    with open_text_maybe_gzip(fa, "rt") as h:
        contigs = SeqIO.to_dict(SeqIO.parse(h, "fasta"))
    if not contigs:
        raise ValueError("No sequences found in FASTA.")
    print(f"Loaded {len(contigs)} contigs.")
    return contigs


def load_xmap(xmap_path):
    print("Parsing XMAP file...")
    df = pd.read_csv(
        open_text_maybe_gzip(xmap_path, "rt"),
        sep="\t",
        comment="#",
        header=None,
        usecols=[1, 2, 3, 4, 5, 6, 7],
        names=["qry", "ref", "qry_start", "qry_end",
               "ref_start", "ref_end", "strand"],
        dtype=str,
        engine="python"
    )

    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    for c in ["qry_start", "qry_end", "ref_start", "ref_end"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["qry", "ref", "qry_start",
                           "qry_end", "ref_start", "ref_end"])

    df["strand"] = df["strand"].map(lambda s: "-" if str(s).startswith("-") else "+")
    m = df["qry_start"] > df["qry_end"]
    df.loc[m, ["qry_start", "qry_end"]] = df.loc[m, ["qry_end", "qry_start"]].values
    m2 = df["ref_start"] > df["ref_end"]
    df.loc[m2, ["ref_start", "ref_end"]] = df.loc[m2, ["ref_end", "ref_start"]].values

    df["align_len"] = df["qry_end"] - df["qry_start"] + 1
    print(f"Valid alignments: {len(df)}")
    return df


def choose_primary_hit(group):
    return group.loc[group["align_len"].idxmax()]


def assign_all_contigs(xmap_df, contigs):
    records = []
    xmap_df = xmap_df.sort_values(
        ["qry", "align_len"], ascending=[True, False]
    ).reset_index(drop=True)

    seen = set()
    for contig_id, grp in xmap_df.groupby("qry"):
        if contig_id not in contigs:
            continue
        seen.add(contig_id)
        primary = choose_primary_hit(grp)
        clen = len(contigs[contig_id])

        records.append(dict(
            contig=contig_id,
            contig_len=clen,
            assigned_ref=str(primary["ref"]),
            strand=str(primary["strand"]),
            support_ref_start=float(primary["ref_start"]),
            support_ref_end=float(primary["ref_end"]),
            support_qry_start=int(primary["qry_start"]),
            support_qry_end=int(primary["qry_end"]),
            support_align_len=int(primary["align_len"]),
            reason="primary_hit_whole_contig",
            order_key=float(primary["ref_start"])
        ))

    for cid in sorted(set(contigs.keys()) - seen):
        clen = len(contigs[cid])
        records.append(dict(
            contig=cid,
            contig_len=clen,
            assigned_ref="unplaced",
            strand="+",
            support_ref_start=pd.NA,
            support_ref_end=pd.NA,
            support_qry_start=pd.NA,
            support_qry_end=pd.NA,
            support_align_len=0,
            reason="unaligned_all",
            order_key=float("inf")
        ))

    return pd.DataFrame.from_records(records)


def build_pseudomolecules(assign_df, contigs, gap_size):
    grouped = defaultdict(list)
    used_contigs = set()

    for row in assign_df.itertuples(index=False):
        if row.assigned_ref == "unplaced":
            continue
        seq = contigs[row.contig].seq
        if row.strand == "-":
            seq = seq.reverse_complement()
        grouped[row.assigned_ref].append((row, str(seq)))
        used_contigs.add(row.contig)

    gap = "N" * gap_size
    records = []
    piece_count = 0

    for ref in natural_sort_refs(grouped.keys()):
        items = grouped[ref]
        items.sort(key=lambda x: (float(x[0].order_key), str(x[0].contig)))
        joined = gap.join(seq for _, seq in items)
        records.append(
            SeqRecord(
                Seq(joined),
                id=f"chr{ref}",
                description=f"pseudomolecule;gap={gap_size}N"
            )
        )
        piece_count += len(items)
        print(f"Built chr{ref} with {len(items)} contigs.")

    return records, used_contigs, piece_count


def write_assign_tsv(assign_df, path):
    cols = [
        "contig", "contig_len", "assigned_ref", "strand",
        "support_ref_start", "support_ref_end",
        "support_qry_start", "support_qry_end",
        "support_align_len", "reason"
    ]
    assign_df[cols].to_csv(path, sep="\t", index=False)
    print(f"Wrote contig assignment table: {path}")


def write_stats(records, contigs, used_contigs, piece_count, assign_df, out_path):
    total_len = sum(len(r.seq) for r in records)
    stats = {
        "total_chromosomes": len(records),
        "total_length": total_len,
        "avg_chrom_length": total_len / len(records) if records else 0,
        "total_contigs": len(contigs),
        "used_contigs_unique": len(used_contigs),
        "contigs_in_pseudomolecules": piece_count,
        "unaligned_all_contigs": int(
            (assign_df["reason"] == "unaligned_all").sum()
        )
    }
    pd.DataFrame([stats]).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote summary statistics: {out_path}")


def main(args):
    contigs = load_contigs(args.contigs_fasta)
    xmap_df = load_xmap(args.xmap_file)

    assign_df = assign_all_contigs(xmap_df, contigs)
    write_assign_tsv(assign_df, args.assign_tsv)

    records, used_contigs, piece_count = build_pseudomolecules(
        assign_df, contigs, gap_size=args.gap_size
    )
    SeqIO.write(records, args.output_fasta, "fasta")
    print(f"Wrote pseudomolecule FASTA: {args.output_fasta}")

    if args.write_unaligned:
        ids = sorted(assign_df.loc[
            assign_df["reason"] == "unaligned_all", "contig"
        ].unique())
        if ids:
            with open(args.unaligned_fasta, "w") as out:
                for cid in ids:
                    SeqIO.write(contigs[cid], out, "fasta")
            print(f"Wrote unaligned contigs: {args.unaligned_fasta}")

    if args.write_unplaced:
        ids = sorted(assign_df.loc[
            assign_df["assigned_ref"] == "unplaced", "contig"
        ].unique())
        if ids:
            with open(args.unplaced_fasta, "w") as out:
                for cid in ids:
                    SeqIO.write(contigs[cid], out, "fasta")
            print(f"Wrote unplaced contigs: {args.unplaced_fasta}")

    write_stats(
        records, contigs, used_contigs,
        piece_count, assign_df, args.stats_tsv
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build pseudomolecules without contig splitting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-x", "--xmap-file", required=True)
    parser.add_argument("-c", "--contigs-fasta", required=True)
    parser.add_argument("-o", "--output-fasta",
                        default="pseudomolecules.fasta")
    parser.add_argument("--gap-size", type=int, default=1000)
    parser.add_argument("--assign-tsv",
                        default="contig_assign.tsv")
    parser.add_argument("--stats-tsv",
                        default="assembly_stats.tsv")
    parser.add_argument("--write-unplaced", action="store_true")
    parser.add_argument("--unplaced-fasta",
                        default="unplaced_contigs.fasta")
    parser.add_argument("--write-unaligned", action="store_true")
    parser.add_argument("--unaligned-fasta",
                        default="unaligned_contigs.fasta")

    args = parser.parse_args()

    try:
        main(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

