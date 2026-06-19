#!/usr/bin/env python3
"""
Filter and subsample E. coli STEC assembly metadata by Stx status, assembly quality, and HC cluster.

Filters:
  - Stx1 == '+' OR Stx2 == '+' (at least one stx gene present)
  - N50 >= 50000
  - Contigs <= 500

Deduplication:
  - One representative per HC{dedup} cluster (removes clones; level configurable)

Subsampling:
  - HC{subsample} clusters with <= max-per-cluster members: keep all
  - HC{subsample} clusters with > max-per-cluster members: random sample

Cluster-frequency report:
  - Always counts/ranks clusters at both HC20 and HC50 in the final output and
    plots each size distribution.  Output filenames, plot titles, and summary
    labels carry the HC level so the two reports are never confused.
"""

import argparse
import logging
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

__version__ = "1.0.0"

# Cluster-frequency reports are always produced at these HC levels.
FREQ_LEVELS = [20, 50]


def build_hc_col(level: int) -> str:
    if level == 0:
        return "HC0 (indistinguishable)(cgMLST V1 + HierCC V1)"
    return f"HC{level}(cgMLST V1 + HierCC V1)"


def setup_logging(log_path: str) -> logging.Logger:
    logger = logging.getLogger("subsample_HC")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("output", help="Output TSV file")
    parser.add_argument("--min-n50", type=int, default=50000, help="Minimum N50 value (default: 50000)")
    parser.add_argument("--max-contigs", type=int, default=500, help="Maximum contig count (default: 500)")
    parser.add_argument("--hc-dedup", type=int, default=0, metavar="LEVEL",
                        help="HC level used for deduplication, e.g. 0, 5, 10 (default: 0)")
    parser.add_argument("--max-per-cluster", type=int, default=10, help="Maximum representatives per subsampling cluster (default: 10)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling (default: 42)")
    parser.add_argument("--hc-subsample", type=int, default=50, metavar="LEVEL",
                        help="HC level used for subsampling, e.g. 20, 50 (default: 50)")
    parser.add_argument("--stx1-col", default="Stx1(Phylotypes)", help="Column name for Stx1 (default: 'Stx1(Phylotypes)')")
    parser.add_argument("--stx2-col", default="Stx2(Phylotypes)", help="Column name for Stx2 (default: 'Stx2(Phylotypes)')")
    args = parser.parse_args()

    base = os.path.splitext(args.output)[0]
    log_path      = base + ".log"
    dropped_path  = base + "_dropped.tsv"
    summary_path  = base + "_summary.txt"
    freq_paths = {lvl: base + f"_HC{lvl}_freq.tsv" for lvl in FREQ_LEVELS}
    plot_paths = {lvl: base + f"_HC{lvl}_freq.png" for lvl in FREQ_LEVELS}

    log = setup_logging(log_path)
    log.info("Started subsample_HC.py v%s", __version__)
    log.info("Input:  %s", args.input)
    log.info("Output: %s", args.output)
    log.info("Parameters: min-n50=%d  max-contigs=%d  hc-dedup=%d  hc-subsample=%d  max-per-cluster=%d  seed=%d",
             args.min_n50, args.max_contigs, args.hc_dedup, args.hc_subsample, args.max_per_cluster, args.seed)
    log.info("Frequency reports (hardcoded): %s", ", ".join(f"HC{lvl}" for lvl in FREQ_LEVELS))

    COL_STX1   = args.stx1_col
    COL_STX2   = args.stx2_col
    COL_N50    = "N50(Assembly stats)"
    COL_CONTIGS = "Contig Number(>=200 bp)(Assembly stats)"
    COL_HC_DEDUP     = build_hc_col(args.hc_dedup)
    COL_HC_SUBSAMPLE = build_hc_col(args.hc_subsample)

    log.debug("Deduplication column: %s", COL_HC_DEDUP)
    log.debug("Subsampling column:   %s", COL_HC_SUBSAMPLE)

    try:
        df = pd.read_csv(args.input, sep='\t', low_memory=False)
    except pd.errors.ParserError as e:
        log.error("Failed to parse input file: %s", e)
        log.error("Hint: if you concatenated files, ensure the second file's header was stripped "
                  "and file1 ends with a newline.  Use: awk 'NR==1 || FNR!=1' file1 file2 > merged.tsv")
        sys.exit(1)
    except Exception as e:
        log.error("Unexpected error reading input file: %s", e)
        sys.exit(1)
    n_input = len(df)
    log.info("Loaded %d rows from %s", n_input, args.input)
    dropped_frames = []

    # Validate required columns
    for col in [COL_STX1, COL_STX2, COL_N50, COL_CONTIGS, COL_HC_DEDUP, COL_HC_SUBSAMPLE]:
        if col not in df.columns:
            log.error("Required column not found: '%s'", col)
            log.error("Available columns: %s", list(df.columns))
            sys.exit(1)

    # Stx filter
    mask = (df[COL_STX1] == "+") | (df[COL_STX2] == "+")
    dropped_frames.append(df[~mask].assign(dropped_at="Stx filter"))
    df = df[mask]
    n_stx = len(df)
    log.info("After Stx filter:      %d  (%d dropped)", n_stx, n_input - n_stx)

    # Assembly quality filter
    mask = (df[COL_N50] >= args.min_n50) & (df[COL_CONTIGS] <= args.max_contigs)
    dropped_frames.append(df[~mask].assign(dropped_at="Quality filter"))
    df = df[mask]
    n_quality = len(df)
    log.info("After quality filter:  %d  (%d dropped)", n_quality, n_stx - n_quality)

    # HC deduplication
    df_dedup = df.drop_duplicates(subset=COL_HC_DEDUP)
    dropped_frames.append(
        df[~df.index.isin(df_dedup.index)].assign(dropped_at=f"HC{args.hc_dedup} deduplication")
    )
    n_dedup = len(df_dedup)
    log.info("After HC%d dedup:       %d  (%d dropped)", args.hc_dedup, n_dedup, n_quality - n_dedup)

    # HC subsampling
    selected_idx = (
        df_dedup.groupby(COL_HC_SUBSAMPLE, group_keys=False)
        .apply(
            lambda x: x if len(x) <= args.max_per_cluster else x.sample(n=args.max_per_cluster, random_state=args.seed),
            include_groups=False,
        )
        .index
    )
    subsampled = df_dedup.loc[selected_idx]
    n_final = len(subsampled)
    dropped_frames.append(
        df_dedup[~df_dedup.index.isin(subsampled.index)].assign(dropped_at=f"HC{args.hc_subsample} subsampling")
    )
    log.info("After HC%d subsample:  %d  (%d dropped)", args.hc_subsample, n_final, n_dedup - n_final)
    log.info("Total dropped: %d  |  Final output: %d", n_input - n_final, n_final)

    # Cluster-frequency analysis on final output, at each hardcoded HC level
    freq_clusters = {}  # level -> n_clusters (or None if column missing)
    for level in FREQ_LEVELS:
        col = build_hc_col(level)
        if col not in subsampled.columns:
            freq_clusters[level] = None
            log.warning("HC%d column not found (%s) — skipping HC%d frequency report", level, col, level)
            continue

        counts = subsampled[col].value_counts().sort_values(ascending=False)
        n_clusters = len(counts)
        freq_clusters[level] = n_clusters
        log.info("HC%d clusters in final output: %d", level, n_clusters)

        freq_df = counts.reset_index()
        freq_df.columns = [f"HC{level}_cluster", "count"]
        freq_df.to_csv(freq_paths[level], index=False, sep='\t')
        log.debug("Written: %s", freq_paths[level])

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.bar(range(len(counts)), counts.values, width=1.0, edgecolor='none')
        ax.set_xlabel(f"HC{level} cluster (ranked by size, largest first)")
        ax.set_ylabel("Number of sequences")
        ax.set_title(
            f"HC{level} cluster size distribution  "
            f"({n_final} sequences  |  {n_clusters} clusters)"
        )
        plt.tight_layout()
        fig.savefig(plot_paths[level], dpi=150)
        plt.close(fig)
        log.debug("Written: %s", plot_paths[level])

    # Write outputs
    subsampled.to_csv(args.output, index=False, sep='\t')
    log.debug("Written: %s", args.output)

    dropped = pd.concat(dropped_frames, ignore_index=True)
    dropped.to_csv(dropped_path, index=False, sep='\t')
    log.debug("Written: %s", dropped_path)

    summary_lines = [
        f"subsample_HC.py v{__version__} summary",
        "=======================",
        f"Input file:              {args.input}",
        f"Output file:             {args.output}",
        f"Dropped file:            {dropped_path}",
        f"Log file:                {log_path}",
    ]
    for level in FREQ_LEVELS:
        summary_lines.append(f"HC{level} freq file:          {freq_paths[level]}")
        summary_lines.append(f"HC{level} plot file:          {plot_paths[level]}")
    summary_lines += [
        "",
        "Parameters:",
        f"  --min-n50              {args.min_n50}",
        f"  --max-contigs          {args.max_contigs}",
        f"  --hc-dedup             {args.hc_dedup}  ({COL_HC_DEDUP})",
        f"  --hc-subsample         {args.hc_subsample}  ({COL_HC_SUBSAMPLE})",
        f"  --max-per-cluster      {args.max_per_cluster}",
        f"  --seed                 {args.seed}",
        f"  --stx1-col             {COL_STX1}",
        f"  --stx2-col             {COL_STX2}",
        "",
        "Filtering steps:",
        f"  Input rows:            {n_input}",
        f"  After Stx filter:      {n_stx}  ({n_input - n_stx} dropped)",
        f"  After quality filter:  {n_quality}  ({n_stx - n_quality} dropped)",
        f"  After HC{args.hc_dedup} dedup:       {n_dedup}  ({n_quality - n_dedup} dropped)",
        f"  After HC{args.hc_subsample} subsample:  {n_final}  ({n_dedup - n_final} dropped)",
        "",
        f"Total dropped:           {n_input - n_final}",
        f"Final output rows:       {n_final}",
    ]
    for level in FREQ_LEVELS:
        note = str(freq_clusters[level]) if freq_clusters[level] is not None else "N/A (column not found)"
        summary_lines.append(f"HC{level} clusters in output: {note}")

    with open(summary_path, "w") as f:
        f.write("\n".join(summary_lines) + "\n")
    log.debug("Written: %s", summary_path)
    log.info("Done.")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        logging.getLogger("subsample_HC").exception("Unhandled exception: %s", e)
        sys.exit(1)
