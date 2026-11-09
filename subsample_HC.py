#!/usr/bin/env python3
"""
Filter and subsample E. coli STEC assembly metadata by Stx status, assembly quality, and HC cluster.

Filters:
  - Stx1 == '+' OR Stx2 == '+' (at least one stx gene present)
  - N50 >= 50000
  - Contigs <= 500

Deduplication:
  - One randomly chosen representative per HC{dedup} cluster (removes clones; level configurable)

Subsampling:
  - HC{subsample} clusters with <= max-per-cluster members: keep all
  - HC{subsample} clusters with > max-per-cluster members: random sample

Within source (--within-source, optional):
  - Deduplication and subsampling run independently within each source type
    (e.g. Human, Bovine), so a clone seen in two hosts keeps one representative
    per host.  Results are combined into a single output.  Rows with a blank
    source are grouped as 'Unknown'.

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


def dedup_and_subsample(df, col_dedup, col_subsample, max_per_cluster, seed):
    """Return (deduplicated, subsampled) frames for one group of rows."""
    # Shuffle first so the representative kept per cluster is random (reproducible via --seed),
    # then restore original row order
    df_dedup = df.sample(frac=1, random_state=seed).drop_duplicates(subset=col_dedup).sort_index()
    if df_dedup.empty:
        return df_dedup, df_dedup
    selected_idx = (
        df_dedup.groupby(col_subsample, group_keys=False)
        .apply(
            lambda x: x if len(x) <= max_per_cluster else x.sample(n=max_per_cluster, random_state=seed),
            include_groups=False,
        )
        .index
    )
    return df_dedup, df_dedup.loc[selected_idx]


def write_freq_reports(subsampled, prefix, log):
    """Write HC cluster-frequency TSV + plot for each FREQ_LEVEL. Returns {level: (n_clusters, tsv, png)}."""
    results = {}
    n_seqs = len(subsampled)
    for level in FREQ_LEVELS:
        col = build_hc_col(level)
        tsv_path = prefix + f"_HC{level}_freq.tsv"
        png_path = prefix + f"_HC{level}_freq.png"
        if col not in subsampled.columns:
            results[level] = (None, tsv_path, png_path)
            log.warning("HC%d column not found (%s) — skipping HC%d frequency report", level, col, level)
            continue

        counts = subsampled[col].value_counts().sort_values(ascending=False)
        n_clusters = len(counts)
        results[level] = (n_clusters, tsv_path, png_path)
        log.info("HC%d clusters in final output: %d", level, n_clusters)

        freq_df = counts.reset_index()
        freq_df.columns = [f"HC{level}_cluster", "count"]
        freq_df.to_csv(tsv_path, index=False, sep='\t')
        log.debug("Written: %s", tsv_path)

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.bar(range(len(counts)), counts.values, width=1.0, edgecolor='none')
        ax.set_xlabel(f"HC{level} cluster (ranked by size, largest first)")
        ax.set_ylabel("Number of sequences")
        ax.set_title(f"HC{level} cluster size distribution  ({n_seqs} sequences  |  {n_clusters} clusters)")
        plt.tight_layout()
        fig.savefig(png_path, dpi=150)
        plt.close(fig)
        log.debug("Written: %s", png_path)
    return results


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
    parser.add_argument("--seed", type=int, default=42, help="Random seed for dedup and subsampling (default: 42)")
    parser.add_argument("--hc-subsample", type=int, default=50, metavar="LEVEL",
                        help="HC level used for subsampling, e.g. 20, 50 (default: 50)")
    parser.add_argument("--stx1-col", default="Stx1(Phylotypes)", help="Column name for Stx1 (default: 'Stx1(Phylotypes)')")
    parser.add_argument("--stx2-col", default="Stx2(Phylotypes)", help="Column name for Stx2 (default: 'Stx2(Phylotypes)')")
    parser.add_argument("--serotype", default=None, metavar="O_ANTIGEN",
                        help="Keep only rows matching this O antigen, e.g. O157 "
                             "(default: no serotype filter)")
    parser.add_argument("--serotype-col", default="O Antigen(Serotype Prediction)",
                        help="Column name for serotype filtering (default: 'O Antigen(Serotype Prediction)')")
    parser.add_argument("--within-source", action="store_true",
                        help="Run deduplication and subsampling separately within each source type "
                             "(host), then combine into one output (default: off)")
    parser.add_argument("--source-col", default="Source Type",
                        help="Column name for source type, used with --within-source (default: 'Source Type')")
    args = parser.parse_args()

    base = os.path.splitext(args.output)[0]
    log_path      = base + ".log"
    dropped_path  = base + "_dropped.tsv"
    summary_path  = base + "_summary.txt"

    log = setup_logging(log_path)
    log.info("Started subsample_HC.py v%s", __version__)
    log.info("Input:  %s", args.input)
    log.info("Output: %s", args.output)
    log.info("Parameters: min-n50=%d  max-contigs=%d  hc-dedup=%d  hc-subsample=%d  max-per-cluster=%d  seed=%d",
             args.min_n50, args.max_contigs, args.hc_dedup, args.hc_subsample, args.max_per_cluster, args.seed)
    log.info("Serotype filter: %s", args.serotype if args.serotype is not None else "(none)")
    log.info("Within source: %s", f"yes ({args.source_col})" if args.within_source else "no")
    log.info("Frequency reports (hardcoded): %s", ", ".join(f"HC{lvl}" for lvl in FREQ_LEVELS))

    COL_STX1   = args.stx1_col
    COL_STX2   = args.stx2_col
    COL_SEROTYPE = args.serotype_col
    COL_SOURCE = args.source_col
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
    required_cols = [COL_STX1, COL_STX2, COL_N50, COL_CONTIGS, COL_HC_DEDUP, COL_HC_SUBSAMPLE]
    if args.serotype is not None:
        required_cols.append(COL_SEROTYPE)
    if args.within_source:
        required_cols.append(COL_SOURCE)
    for col in required_cols:
        if col not in df.columns:
            log.error("Required column not found: '%s'", col)
            log.error("Available columns: %s", list(df.columns))
            sys.exit(1)

    # Serotype prefilter (optional)
    if args.serotype is not None:
        mask = df[COL_SEROTYPE] == args.serotype
        dropped_frames.append(df[~mask].assign(dropped_at="Serotype filter"))
        df = df[mask]
        n_serotype = len(df)
        log.info("After serotype filter: %d  (%d dropped)  [%s == %s]",
                 n_serotype, n_input - n_serotype, COL_SEROTYPE, args.serotype)
    else:
        n_serotype = n_input

    # Stx filter
    mask = (df[COL_STX1] == "+") | (df[COL_STX2] == "+")
    dropped_frames.append(df[~mask].assign(dropped_at="Stx filter"))
    df = df[mask]
    n_stx = len(df)
    log.info("After Stx filter:      %d  (%d dropped)", n_stx, n_serotype - n_stx)

    # Assembly quality filter
    mask = (df[COL_N50] >= args.min_n50) & (df[COL_CONTIGS] <= args.max_contigs)
    dropped_frames.append(df[~mask].assign(dropped_at="Quality filter"))
    df = df[mask]
    n_quality = len(df)
    log.info("After quality filter:  %d  (%d dropped)", n_quality, n_stx - n_quality)

    # Group rows: one group per source type, or a single group of everything
    if args.within_source:
        source_labels = df[COL_SOURCE].fillna("").astype(str).str.strip().replace("", "Unknown")
        groups = [(label, df[source_labels == label]) for label in sorted(source_labels.unique())]
        log.info("Source types (%s): %d  [%s]", COL_SOURCE, len(groups),
                 ", ".join(f"{label}={len(g)}" for label, g in groups))
    else:
        groups = [(None, df)]

    # Dedup + subsample within each group, then combine
    group_results = []
    subsampled_frames = []
    for label, gdf in groups:
        tag = f"[{label}] " if label is not None else ""
        n_in = len(gdf)

        df_dedup, sub = dedup_and_subsample(
            gdf, COL_HC_DEDUP, COL_HC_SUBSAMPLE, args.max_per_cluster, args.seed)
        n_dedup_g = len(df_dedup)
        n_final_g = len(sub)

        dropped_frames.append(
            gdf[~gdf.index.isin(df_dedup.index)].assign(dropped_at=f"{tag}HC{args.hc_dedup} deduplication")
        )
        dropped_frames.append(
            df_dedup[~df_dedup.index.isin(sub.index)].assign(dropped_at=f"{tag}HC{args.hc_subsample} subsampling")
        )
        log.info("%sAfter HC%d dedup:       %d  (%d dropped)", tag, args.hc_dedup, n_dedup_g, n_in - n_dedup_g)
        log.info("%sAfter HC%d subsample:  %d  (%d dropped)", tag, args.hc_subsample, n_final_g, n_dedup_g - n_final_g)

        subsampled_frames.append(sub)
        group_results.append(dict(label=label, n_in=n_in, n_dedup=n_dedup_g, n_final=n_final_g))

    # back to input row order (empty frame if nothing survived the filters)
    subsampled = pd.concat(subsampled_frames).sort_index() if subsampled_frames else df.iloc[0:0]
    n_dedup = sum(r["n_dedup"] for r in group_results)
    n_final = len(subsampled)
    if args.within_source:
        log.info("All sources combined: dedup %d, subsample %d", n_dedup, n_final)
    log.info("Total dropped: %d  |  Final output: %d", n_input - n_final, n_final)

    # Cluster-frequency analysis on final output
    freq = write_freq_reports(subsampled, base, log)

    # Write outputs
    subsampled.to_csv(args.output, index=False, sep='\t')
    log.debug("Written: %s", args.output)

    dropped = pd.concat(dropped_frames, ignore_index=True)
    dropped.to_csv(dropped_path, index=False, sep='\t')
    log.debug("Written: %s", dropped_path)

    # Summary
    summary_lines = [
        f"subsample_HC.py v{__version__} summary",
        "=======================",
        f"Input file:              {args.input}",
        f"Output file:             {args.output}",
        f"Dropped file:            {dropped_path}",
        f"Log file:                {log_path}",
    ]
    for level in FREQ_LEVELS:
        _, tsv_path, png_path = freq[level]
        summary_lines.append(f"HC{level} freq file:          {tsv_path}")
        summary_lines.append(f"HC{level} plot file:          {png_path}")
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
        f"  --serotype             {args.serotype if args.serotype is not None else '(none)'}",
        f"  --serotype-col         {COL_SEROTYPE}",
        f"  --within-source        {'yes' if args.within_source else 'no'}",
        f"  --source-col           {COL_SOURCE}",
        "",
        "Filtering steps:",
        f"  Input rows:            {n_input}",
        f"  After serotype filter: {n_serotype}  ({n_input - n_serotype} dropped)",
        f"  After Stx filter:      {n_stx}  ({n_serotype - n_stx} dropped)",
        f"  After quality filter:  {n_quality}  ({n_stx - n_quality} dropped)",
        f"  After HC{args.hc_dedup} dedup:       {n_dedup}  ({n_quality - n_dedup} dropped)"
        + ("  (within each source)" if args.within_source else ""),
        f"  After HC{args.hc_subsample} subsample:  {n_final}  ({n_dedup - n_final} dropped)"
        + ("  (within each source)" if args.within_source else ""),
    ]
    if args.within_source:
        summary_lines += ["", "Per source (after filters -> dedup -> subsample):"]
        width = max(len(r["label"]) for r in group_results) if group_results else 0
        for r in group_results:
            summary_lines.append(f"  {r['label']:<{width}}  {r['n_in']:>7} -> {r['n_dedup']:>7} -> {r['n_final']:>7}")
    summary_lines += [
        "",
        f"Total dropped:           {n_input - n_final}",
        f"Final output rows:       {n_final}",
    ]
    for level in FREQ_LEVELS:
        n_clusters = freq[level][0]
        note = str(n_clusters) if n_clusters is not None else "N/A (column not found)"
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
