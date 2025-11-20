#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine and adjust three XMAP files (Bionano).
Step 1: Use file1 as base. For contigs with main alignment coverage
        < MIN_COV_RATIO of query length, add the longest non-duplicate
        alignment from file2.
Step 2: Modify chr2 and chr20
"""

import pandas as pd

# ===== User parameters =====
file1 = "lm15_0820_2_T2T_08222025.xmap"    # primary alignments
file2 = "LM15_0820_2_T2T_08242025.xmap"    # additional alignments
file3 = "lm15_0820_2_Lo7_08242025.xmap"    # alignments to lo7
outfile = "20.lm15_11X_2_CS_IAAS_lo7.xmap"
statfile = "20.xmap_stats.txt"

# Chromosome IDs
CHR1 = 1
CHR2 = 2
CHR20 = 20

# Thresholds (Mb)
CHR2_CUT_MB     = 243   # chr2 cut point
CHR20_LIMIT_MB  = 357   # chr20 delete/selection threshold
CHR1_INSERT_MB  = 268   # chr1 length used for chr20 replacement (file3)

# Step 1 coverage threshold
MIN_COV_RATIO = 0.6


# =================== Common helpers ===================

def read_xmap(path):
    header_cols = None
    header_line = None
    with open(path) as f:
        for line in f:
            if line.startswith("#h "):
                header_line = line.rstrip("\n")
                header_cols = line.strip().split()[1:]  # drop "#h"
                break
    if header_cols is None:
        raise ValueError(f"Header line starting with '#h' not found in {path}")

    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=header_cols)
    for c in [
        "RefContigID", "RefStartPos", "RefEndPos",
        "QryStartPos", "QryEndPos", "QryLen", "QryContigID"
    ]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, header_cols, header_line


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


# =================== Step 1: supplement low coverage ===================

def step1_supplement(df1, df2):
    """
    For each query contig in df1, if an alignment covers less than
    MIN_COV_RATIO of query length, add the longest non-duplicate
    alignment from df2 for that contig. Collect per-contig stats.
    """
    df1 = df1.copy()
    df2 = df2.copy()
    df1["AlignLen"] = df1["RefEndPos"] - df1["RefStartPos"]
    df2["AlignLen"] = df2["RefEndPos"] - df2["RefStartPos"]

    out_rows = []
    # contig -> dict(before=, after=, added=[(RefContigID, AlignLen, RefStart-RefEnd, overlap_bool)])
    supplemented = {}

    for contig, group in df1.groupby("QryContigID", sort=False):
        before_count = len(group)
        rows_to_add = []
        added_info = []

        for _, row in group.iterrows():
            # always keep the main alignment from file1
            rows_to_add.append(row)

            align_len = row["AlignLen"]
            qlen = row["QryLen"]
            if pd.isna(align_len) or pd.isna(qlen):
                continue

            if align_len < MIN_COV_RATIO * qlen:
                cand = df2[df2["QryContigID"] == contig].copy()
                if len(cand) == 0:
                    continue

                # drop nearly identical alignment to avoid duplication
                cand = cand[~(
                    (cand["RefContigID"] == row["RefContigID"]) &
                    (abs(cand["RefStartPos"] - row["RefStartPos"]) < 1e5) &
                    (abs(cand["RefEndPos"]   - row["RefEndPos"])   < 1e5)
                )]
                if len(cand) == 0:
                    continue

                # pick alignment with maximum AlignLen
                best = cand.loc[cand["AlignLen"].idxmax()]

                # compute overlap flag (for reporting only)
                ovl = overlap_len(row["RefStartPos"], row["RefEndPos"],
                                  best["RefStartPos"], best["RefEndPos"])
                small = min(max(row["AlignLen"], 0), max(best["AlignLen"], 0))
                overlap_flag = (ovl > 0.0 and small > 0 and ovl / small >= 0.3)

                rows_to_add.append(best)
                added_info.append((
                    int(best["RefContigID"]),
                    int(best["AlignLen"]),
                    f"{int(best['RefStartPos'])}-{int(best['RefEndPos'])}",
                    overlap_flag
                ))

        after_count = len(rows_to_add)
        if added_info:
            key = str(int(contig)) if pd.notna(contig) else str(contig)
            supplemented[key] = {
                "before": before_count,
                "after": after_count,
                "added": added_info
            }
        out_rows.extend(rows_to_add)

    out_df = pd.DataFrame(out_rows).reset_index(drop=True)
    out_df = out_df[df1.columns]
    return out_df, supplemented


# =================== Step 2: chr2 and chr20 transformations ===================

def step2_chr2_replace_and_chr20_delete_replace(df1, df3):
    """
    Step 2: apply structural changes on chr2 and chr20.
    - chr2:
      * remove rows with RefEndPos <= CHR2_CUT_MB;
      * build an insertion block from chr20 (rows <= CHR20_LIMIT_MB plus
        whole contigs that cross CHR20_LIMIT_MB), remap to chr2;
      * shift the remaining chr2 tail by delta = max(L_insert - CHR2_CUT_MB, 0).

    - chr20:
      * delete head by whole-contig rule (rows <= CHR20_LIMIT_MB and all
        rows of contigs crossing CHR20_LIMIT_MB);
      * add chr1 rows from file3 with RefEndPos <= CHR1_INSERT_MB,
        remapped to chr20 (coordinates unchanged).
    """
    cut_bp_chr2    = int(CHR2_CUT_MB    * 1_000_000)
    limit_bp_chr20 = int(CHR20_LIMIT_MB * 1_000_000)
    ins_bp_chr1    = int(CHR1_INSERT_MB * 1_000_000)

    df1 = df1.copy()
    df3 = df3.copy()

    # ----- Build chr20 block to insert into chr2 (with whole-contig rule) -----
    chr20_df = df1.loc[df1["RefContigID"] == CHR20].copy()

    # rows fully ending within the limit
    within = chr20_df.loc[chr20_df["RefEndPos"] <= limit_bp_chr20].copy()

    # contigs crossing the threshold, then take all rows of these contigs
    crossing_rows = chr20_df.loc[
        (chr20_df["RefStartPos"] < limit_bp_chr20) &
        (chr20_df["RefEndPos"]   > limit_bp_chr20)
    ]
    crossing_contigs = set(
        crossing_rows["QryContigID"].dropna().astype(str).tolist()
    )

    whole_contig_part = chr20_df.loc[
        chr20_df["QryContigID"].astype(str).isin(crossing_contigs)
    ].copy()

    chr20_selected = (
        pd.concat([within, whole_contig_part], ignore_index=True)
        .drop_duplicates()
    )

    # length of insertion block (in bp, on chr20 coordinates)
    if len(chr20_selected) > 0:
        L_insert_bp = int(chr20_selected["RefEndPos"].max())
    else:
        L_insert_bp = cut_bp_chr2

    # map selected chr20 block to chr2
    chr20_to_chr2 = chr20_selected.copy()
    chr20_to_chr2["RefContigID"] = CHR2

    # ----- chr2: delete head and shift tail -----
    chr2_mask = (df1["RefContigID"] == CHR2)
    removed_chr2_rows = int(
        (df1.loc[chr2_mask, "RefEndPos"] <= cut_bp_chr2).sum()
    )

    tail = df1.loc[chr2_mask & (df1["RefEndPos"] > cut_bp_chr2)].copy()

    truncated_chr2_contigs = (
        tail.loc[tail["RefStartPos"] < cut_bp_chr2, "QryContigID"]
        .dropna().astype(str).unique().tolist()
    )

    # truncate spanning rows at the cut position
    tail.loc[tail["RefStartPos"] < cut_bp_chr2, "RefStartPos"] = cut_bp_chr2

    # shift tail by delta
    delta_bp = max(L_insert_bp - cut_bp_chr2, 0)
    if delta_bp != 0:
        tail["RefStartPos"] = tail["RefStartPos"] + delta_bp
        tail["RefEndPos"]   = tail["RefEndPos"]   + delta_bp

    # ----- chr20: delete head by whole-contig rule -----
    chr20_del_mask_within = (
        (df1["RefContigID"] == CHR20) &
        (df1["RefEndPos"] <= limit_bp_chr20)
    )

    chr20_crossing_rows = df1.loc[
        (df1["RefContigID"] == CHR20) &
        (df1["RefStartPos"] < limit_bp_chr20) &
        (df1["RefEndPos"]   > limit_bp_chr20)
    ]
    chr20_del_contigs = set(
        chr20_crossing_rows["QryContigID"].dropna().astype(str).tolist()
    )

    chr20_del_mask_contig = (
        (df1["RefContigID"] == CHR20) &
        (df1["QryContigID"].astype(str).isin(chr20_del_contigs))
    )

    chr20_del_mask = chr20_del_mask_within | chr20_del_mask_contig

    removed_chr20_rows = int(chr20_del_mask.sum())
    removed_chr20_contigs_count = len(chr20_del_contigs)
    removed_chr20_contigs_list = sorted(list(chr20_del_contigs))

    # chr20 rows kept after deletion
    chr20_residual = df1.loc[
        (df1["RefContigID"] == CHR20) & (~chr20_del_mask)
    ].copy()

    # ----- file3: chr1 head mapped to chr20 -----
    chr1_ins = df3.loc[
        (df3["RefContigID"] == CHR1) &
        (df3["RefEndPos"] <= ins_bp_chr1)
    ].copy()
    chr1_to_chr20 = chr1_ins.copy()
    chr1_to_chr20["RefContigID"] = CHR20

    # ----- Assemble final dataframe -----
    # 1) all chromosomes except chr2 and chr20
    others = df1.loc[
        (df1["RefContigID"] != CHR2) &
        (df1["RefContigID"] != CHR20)
    ].copy()

    # 2) final chr20 = residual chr20 + chr1->chr20
    chr20_final = pd.concat(
        [chr20_residual, chr1_to_chr20], ignore_index=True
    )

    # 3) final chr2 = chr20->chr2 insertion block + shifted tail
    chr2_final = pd.concat([chr20_to_chr2, tail], ignore_index=True)

    # 4) merge everything and sort
    out = pd.concat([others, chr20_final, chr2_final], ignore_index=True)
    out = out.sort_values(
        ["RefContigID", "RefStartPos", "RefEndPos"]
    ).reset_index(drop=True)

    stats = {
        # chr2
        "removed_chr2_rows": removed_chr2_rows,
        "truncated_chr2_contigs": truncated_chr2_contigs,
        "chr20_selected_rows_for_chr2": int(len(chr20_selected)),
        "chr20_crossing_contigs_for_chr2_count": int(len(crossing_contigs)),
        "chr20_crossing_contigs_for_chr2": sorted(list(crossing_contigs)),
        "L_insert_bp": int(L_insert_bp),
        "L_insert_Mb": round(L_insert_bp / 1_000_000, 3),
        "delta_bp": int(delta_bp),
        "delta_Mb": round(delta_bp / 1_000_000, 3),

        # chr20 delete and replace
        "removed_chr20_rows": removed_chr20_rows,
        "removed_chr20_contigs_count": removed_chr20_contigs_count,
        "removed_chr20_crossing_contigs": removed_chr20_contigs_list,
        "chr1_to_chr20_rows": int(len(chr1_to_chr20)),
        "chr20_residual_rows": int(len(chr20_residual)),
        "chr20_gap_Mb_between_268_and_357":
            max(CHR20_LIMIT_MB - CHR1_INSERT_MB, 0)
    }
    return out, stats


# =================== Main ===================

def main():
    df1, header_cols, header_line = read_xmap(file1)
    df2, _, _ = read_xmap(file2)
    df3, _, _ = read_xmap(file3)

    # Step 1: supplement low coverage
    step1_df, supplemented = step1_supplement(df1, df2)

    # Step 2: chr2 replace + shift; chr20 delete + replace
    final_df, step2_stats = step2_chr2_replace_and_chr20_delete_replace(
        step1_df, df3
    )

    # write XMAP result, keep original #h header
    with open(outfile, "w") as fo:
        fo.write(header_line + "\n")
    final_df.to_csv(outfile, sep="\t", index=False, mode="a")

    # write statistics
    with open(statfile, "w") as f:
        f.write("XMAP processing statistics\n\n")
        f.write(f"Total rows in output: {len(final_df)}\n\n")

        # Step 1 stats
        f.write("Step 1: coverage supplement\n")
        f.write(
            f"  Contigs with supplementary alignments: {len(supplemented)}\n"
        )
        if supplemented:
            for contig, info in supplemented.items():
                f.write(
                    f"  - QryContig {contig}: "
                    f"{info['before']} -> {info['after']}\n"
                )
                for refid, alen, span, ovlp in info["added"]:
                    f.write(
                        "      + "
                        f"RefContigID={refid}, "
                        f"AlignLen={alen}, "
                        f"Span={span}, "
                        f"Overlap30%={ovlp}\n"
                    )
        f.write("\n")

        # Step 2 stats
        f.write("Step 2: chr2 replace/shift and chr20 delete/replace\n")
        # chr2
        f.write(
            f"  [chr2] Rows removed with RefEndPos <= "
            f"{CHR2_CUT_MB} Mb: {step2_stats['removed_chr2_rows']}\n"
        )
        f.write(
            "  [chr2] Number of truncated crossing contigs: "
            f"{len(step2_stats['truncated_chr2_contigs'])}\n"
        )
        if step2_stats["truncated_chr2_contigs"]:
            f.write(
                "        Crossing contigs: " +
                ", ".join(step2_stats["truncated_chr2_contigs"]) +
                "\n"
            )
        f.write(
            "  [chr2] Rows in insertion block from chr20: "
            f"{step2_stats['chr20_selected_rows_for_chr2']}\n"
        )
        f.write(
            "  [chr2] chr20 contigs triggering whole-contig inclusion: "
            f"{step2_stats['chr20_crossing_contigs_for_chr2_count']}\n"
        )
        if step2_stats["chr20_crossing_contigs_for_chr2"]:
            f.write(
                "        These contigs: " +
                ", ".join(step2_stats["chr20_crossing_contigs_for_chr2"]) +
                "\n"
            )
        f.write(
            "  [chr2] Insertion block length L_insert: "
            f"{step2_stats['L_insert_bp']} bp "
            f"(~{step2_stats['L_insert_Mb']} Mb)\n"
        )
        f.write(
            "  [chr2] Tail shift delta: "
            f"{step2_stats['delta_bp']} bp "
            f"(~{step2_stats['delta_Mb']} Mb)\n"
        )

        # chr20
        f.write(
            "  [chr20] Rows removed (<= "
            f"{CHR20_LIMIT_MB} Mb, including whole-contig deletions): "
            f"{step2_stats['removed_chr20_rows']}\n"
        )
        f.write(
            "  [chr20] Crossing contigs removed as whole contigs: "
            f"{step2_stats['removed_chr20_contigs_count']}\n"
        )
        if step2_stats["removed_chr20_crossing_contigs"]:
            f.write(
                "        These contigs: " +
                ", ".join(step2_stats["removed_chr20_crossing_contigs"]) +
                "\n"
            )
        f.write(
            "  [chr20] Rows from file3 chr1 <= "
            f"{CHR1_INSERT_MB} Mb mapped to chr20: "
            f"{step2_stats['chr1_to_chr20_rows']}\n"
        )
        f.write(
            "  [chr20] Rows kept on chr20 after deletion "
            f"(> {CHR20_LIMIT_MB} Mb and not crossing): "
            f"{step2_stats['chr20_residual_rows']}\n"
        )
        f.write(
            "  [chr20] Expected gap between 268 Mb and 357 Mb (Mb): "
            f"{step2_stats['chr20_gap_Mb_between_268_and_357']}\n"
        )

    print(f"Output written to: {outfile}")
    print(f"Stats written to: {statfile}")


if __name__ == "__main__":
    main()