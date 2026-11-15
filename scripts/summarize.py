#!/usr/bin/env python3
"""Recompute Table 2 and the statistics reported in the manuscript.

Input : data/read_level.tsv
Output: medians (IQR) of F1, host recall, target retention per route x host fraction;
        Friedman tests, Kendall's W, and exact Wilcoxon signed-rank tests (Holm).
"""
import csv, collections, itertools
import numpy as np
from scipy import stats

ROUTES = ["nohuman", "deacon_neg", "deacon_pos34", "mm2_pos34", "cleanifier"]
LEVELS = ["90%", "50%", "10%", "1%"]
RATIO = {"9:1": "90%", "1:1": "50%", "1:9": "10%", "1:99": "1%"}


def level(host_frac):
    if host_frac > 0.7:
        return "90%"
    if host_frac > 0.3:
        return "50%"
    if host_frac > 0.05:
        return "10%"
    return "1%"


def load(path="data/read_level.tsv"):
    f1 = collections.defaultdict(list)
    tret = collections.defaultdict(list)
    for row in csv.DictReader(open(path), delimiter="\t"):
        key = (RATIO[row["ratio"]], row["route"])
        tp, fp, fn, tn = (int(row[k]) for k in ("tp", "fp", "fn", "tn"))
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1[key].append(2 * prec * rec / (prec + rec) if prec + rec else 0.0)
        tret[key].append(tn / (tn + fp) if tn + fp else 0.0)
    return f1, tret


def main():
    f1, tret = load()
    print("route\tcondition\tn\tF1 median (IQR)\ttarget retention median (IQR)")
    for lv in LEVELS:
        for r in ROUTES:
            v = np.array(f1[(lv, r)])
            t = np.array(tret[(lv, r)])
            q1, q3 = np.percentile(v, [25, 75])
            t1, t3 = np.percentile(t, [25, 75])
            print(f"{r}\t{lv}\t{len(v)}\t{np.median(v):.5f} ({q1:.5f}-{q3:.5f})"
                  f"\t{np.median(t)*100:.2f} ({t1*100:.2f}-{t3*100:.2f})")

    print("\ncondition\tFriedman chi2\tKendall W\tglobal P")
    for lv in LEVELS:
        mat = np.array([f1[(lv, r)] for r in ROUTES]).T
        chi2, p = stats.friedmanchisquare(*[mat[:, i] for i in range(len(ROUTES))])
        print(f"{lv}\t{chi2:.2f}\t{chi2/(mat.shape[0]*(len(ROUTES)-1)):.3f}\t{p:.3g}")

    print("\nExact Wilcoxon signed-rank, all ten pairwise comparisons, Holm-corrected")
    for lv in LEVELS:
        pairs = []
        for a, b in itertools.combinations(ROUTES, 2):
            d = np.array(f1[(lv, a)]) - np.array(f1[(lv, b)])
            if np.all(d == 0):
                continue
            pairs.append((stats.wilcoxon(d, method="exact").pvalue, a, b))
        pairs.sort()
        m = len(pairs)
        running = 0.0
        out = []
        for i, (p, a, b) in enumerate(pairs):
            adj = min(1.0, p * (m - i))
            running = max(running, adj)
            out.append((a, b, p, running))
        print(f"[{lv}] m={m}, max adjusted P = {running:.3g}")
        for a, b, p, adj in out:
            print(f"   {a} vs {b}: raw {p:.3g}, adj {adj:.3g}")


if __name__ == "__main__":
    main()
