# Host-depletion route benchmark — data and code

Data and analysis code for the manuscript:

> **Host fraction changes F1-leading read-level host-depletion routes in simulated
> targeted Oxford Nanopore libraries**

## What this repository contains

| Path | Contents |
|---|---|
| `data/read_level.tsv` | Per-library, per-route read-level confusion matrix (100 libraries × 5 routes = 500 rows): TP/FP/FN/TN against host-positive truth |
| `data/assembly.tsv` | Whole-assembly metrics for 100 libraries × 2 assemblers (contigs, bases, CheckM completeness/contamination) |
| `data/contig_attribution.tsv` | Target-genome recovery after contig attribution (200 rows) |
| `data/binning.tsv` | Best bin per library and its coverage of the strain reference (200 rows) |
| `data/target_library_34_accessions.txt` | The 34-accession target library used for positive selection (the five simulated strains are excluded by name) |
| `scripts/` | Analysis and figure scripts (Python 3.13, matplotlib, scipy, numpy) |

## Experimental design

- **Target genomes:** five *Orientia tsutsugamushi* assemblies (UT76, CH219, Karp, JJOtsu1, UT176).
- **Host:** GRCh38.
- **Host:target ratios:** 9:1, 1:1, 1:9, 1:99 (≈90%, 50%, 10%, 1% host content).
- **Replicates:** five independent simulation seeds per strain × ratio combination.
- **Total:** 5 strains × 4 ratios × 5 seeds = **100 libraries**.
- **Reads:** Oxford Nanopore-style long reads, ~2.57 M reads (~4.2 Gb) per library.

Sample numbering within each seed: samples 0–3 = UT76 (1:9, 9:1, 1:1, 1:99), 4–7 = CH219,
8–11 = Karp, 12–15 = JJOtsu1, 16–19 = UT176.

## Read-level routes

| Route | Strategy | Reference |
|---|---|---|
| noHuman | Kraken2 taxonomic assignment, LCA rule | HPRC human database |
| Deacon negative | minimizer set membership, hit → discard | GRCh38 panhuman index |
| Deacon positive | minimizer set membership, hit → keep | 34-accession *Orientia* library |
| Minimap2 positive | read-to-contig alignment, hit → keep | same 34-accession library |
| Cleanifier | probabilistic k-mer filter | k=29, w=33 human index |

Metrics use **host reads as the positive class**: Accuracy = (TP+TN)/N,
Precision = TP/(TP+FP), Recall = TP/(TP+FN), F1 = 2PR/(P+R).
Target-read retention = TN/(TN+FP).

## Software versions

Deacon 0.19.1 · noHuman 0.5.1 · Cleanifier 1.3.2 · Minimap2 2.31-r1302 ·
Canu master r10522 (commit 86d7d3b) · Flye 2.9.6-b1802 · MetaBAT2 2.18 ·
CheckM 1.2.5 (checkm-genome) · samtools 1.23.1 (htslib 1.23.1) ·
Prodigal 2.6.3 · HMMER 3.4 · Python 3.15.2

## How to reproduce the summary tables

```bash
python3 scripts/summarize.py
```

This recomputes, from `data/read_level.tsv`:

- median (IQR) F1 per route × host fraction
- host recall and target retention per route × host fraction
- Kendall's W and Friedman tests per host fraction
- exact two-sided Wilcoxon signed-rank tests with Holm correction
  (ten pairwise comparisons; SciPy 1.18.1, statsmodels 0.15.0; bootstrap seed 20261115)

and prints the values reported in Table 2 and Figure 1.

## Notes on aggregation

- The **primary analysis uses medians with IQR across all 25 libraries** per cell.
- The **threshold-scan tables report means** (Deacon scan n = 20; Minimap2 scan n = 100).
  These statistics are not interchangeable: for Minimap2 target retention the median is
  64.52% while the mean is 56.28% (IQR 49.09–65.20%), because the distribution is skewed.
- The 25 libraries are the experimental unit; individual reads are not treated as
  independent observations.

## Data availability

Per-library tables are provided here as plain TSV. Larger intermediates
(assemblies, BAM files, CheckM working directories) are available from the
authors on reasonable request.

## License

Code: MIT. Data: CC BY 4.0.
