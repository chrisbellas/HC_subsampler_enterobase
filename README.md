# HC_subsampler_enterobase

**Version 1.0.0**

Takes EnteroBase metadata TSV files and filters and subsamples them based on
EnteroBase's HierCC (Hierarchical Clustering). HC clusters organize core-genome
Multi-Locus Sequence Typing (cgMLST) profiles into nested, evolutionary clusters,
where lower HC levels (e.g. HC0, HC5) group near-identical genomes (recent local
outbreak chains) and higher levels group progressively more distantly related
isolates, up to HC1100, which corresponds to legacy MLST ST Complexes in
*E. coli*.

## Uses

Use the script for dereplicating and downsampling data whilst preserving diversity.
Useful for curating genomes for Machine Learning approaches — especially reducing
the O157 sampling bias.

## What it does

`subsample_HC.py` runs the following pipeline on an EnteroBase metadata TSV:

0. **Serotype prefilter** *(optional)* - if `--serotype` is given, keeps only rows
   whose O antigen matches (e.g. `--serotype O157`). Skipped by default.
1. **Stx filter** - keeps rows where Stx1 *or* Stx2 is present (`+`).
2. **Quality filter** - keeps rows with `N50 >= --min-n50` and
   `Contigs <= --max-contigs`.
3. **Deduplication** - keeps one randomly chosen representative per `--hc-dedup`
   cluster (removes clones; reproducible via `--seed`).
4. **Subsampling** - for each `--hc-subsample` cluster, keeps all members if the
   cluster has `<= --max-per-cluster` members, otherwise takes a random sample of
   `--max-per-cluster` (reproducible via `--seed`).
5. **Cluster-frequency report** — counts and ranks clusters at both HC20 and HC50
   in the final output, writing frequency tables and size-distribution plots.

## Input

The input must be the **full metadata TSV exported from EnteroBase** (Data →
download as tab-separated values), not a trimmed subset. Make sure the export
includes:

- **All HierCC HC columns** (HC0, HC2, HC5, HC10, HC20, HC50, HC100, HC200, HC400,
  HC1100, HC1500, HC2000, HC2350) - the script selects HC levels by name at
  runtime, so any level you might pass to `--hc-dedup` / `--hc-subsample`, plus the
  HC20 and HC50 columns used for the frequency reports, must be present.
- **Stx columns** - `Stx1(Phylotypes)` and `Stx2(Phylotypes)`.
- **Serotype** (O and H antigen) columns.
- **Assembly stats** - `N50(Assembly stats)` and
  `Contig Number(>=200 bp)(Assembly stats)`.

Column names must match EnteroBase's export headers exactly. If you concatenate
multiple exports, see the note at the end of this file.

## Usage

```bash
python3 subsample_HC.py input.tsv output.tsv [options]
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--min-n50` | `50000` | Minimum N50 value |
| `--max-contigs` | `500` | Maximum contig count |
| `--hc-dedup LEVEL` | `0` | HC level used for deduplication |
| `--hc-subsample LEVEL` | `50` | HC level used for subsampling |
| `--max-per-cluster` | `10` | Maximum representatives per subsampling cluster |
| `--seed` | `42` | Random seed for deduplication and subsampling |
| `--stx1-col` | `Stx1(Phylotypes)` | Column name for Stx1 |
| `--stx2-col` | `Stx2(Phylotypes)` | Column name for Stx2 |
| `--serotype O_ANTIGEN` | *(none)* | Keep only rows with this O antigen, e.g. `O157` (prefilter, off by default) |
| `--serotype-col` | `O Antigen(Serotype Prediction)` | Column name for serotype filtering |
| `--version` | | Print version and exit |

### Example

```bash
python3 subsample_HC.py metadata.tsv subsampled.tsv \
    --hc-dedup 5 --hc-subsample 50 --max-per-cluster 10
```

## Outputs

Given an output path `out.tsv`, the script writes:

- `out.tsv` - the filtered, subsampled dataset
- `out_dropped.tsv` - all dropped rows, annotated with the step that dropped them
- `out_summary.txt` - run parameters and per-step counts
- `out.log` — full run log
- `out_HC20_freq.tsv`, `out_HC50_freq.tsv` — cluster-frequency tables
- `out_HC20_freq.png`, `out_HC50_freq.png` — cluster size-distribution plots

## Requirements

- Python 3
- pandas
- matplotlib

## Notes

If you concatenated multiple EnteroBase exports, make sure the second file's
header was stripped and the first file ends with a newline, e.g.:

```bash
awk 'NR==1 || FNR!=1' file1.tsv file2.tsv > merged.tsv
```
