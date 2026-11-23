#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 6: locus-level dosage-phenotype figure and statistics.

Panel (a): plant height by EasyHap haplotype class (Hap001/Hap002/Hap003 +
           "Other"), with compact letter display from two-sided Mann-Whitney U
           tests and Benjamini-Hochberg correction.
Panel (b): mean relative dosage of the selected k-mers vs plant height, with an
           ordinary least-squares fit, its 95% confidence band, the fitted
           equation, R2 and Spearman's rho.

Statistical conventions (must match the manuscript):
  * the mean dosage of each accession is computed over the markers of the VCF
    (missing values ignored) and STORED TO FOUR DECIMAL PLACES before any
    statistic is calculated;
  * Mann-Whitney U tests are always two-sided and share one BH family within a
    figure panel;
  * pairwise comparisons use the groups displayed in panel (a), i.e. the
    `--top-n` largest classes plus the pooled "Other" class.

Usage:
    python scripts/06_plot_dosage.py \
        --vcf region_bonf.vcf \
        --hapgroup easyhap_bonf/Chr4G_3213203_3231114.HapGroup.tsv \
        --out dosage_vs_height_bonf.png

Environment variables (see config.sh): WORKDIR, SAMPLE_ORDER, TRAITS
"""
import argparse
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

DECIMALS = 4


def bh_adjust(p):
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n, float)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * n / (rank + 1.0))
        q[i] = prev
    return q


def mwu(a, b):
    """Two-sided Mann-Whitney U test (the default alternative differs between
    SciPy versions, so it is always passed explicitly)."""
    return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue


def letter_display(labels, pmat):
    """Compact letter display from a BH-adjusted pairwise P matrix."""
    different = {}
    for i, j in itertools.combinations(range(len(labels)), 2):
        different[(i, j)] = different[(j, i)] = pmat[i][j] < 0.05
    cliques = []
    for k in range(len(labels), 1, -1):
        for sub in itertools.combinations(range(len(labels)), k):
            if all(not different[(a, b)] for a, b in itertools.combinations(sub, 2)):
                if not any(set(sub) < set(c) for c in cliques):
                    cliques.append(sub)
    letters = [""] * len(labels)
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    nxt = 0
    for c in cliques:
        ch = alphabet[nxt]
        nxt += 1
        for g in c:
            letters[g] += ch
    for i in range(len(labels)):
        if not letters[i]:
            letters[i] = alphabet[nxt]
            nxt += 1
    return letters


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="VCF from step 4 (markers are read from the ID column)")
    ap.add_argument("--hapgroup", default=None, help="EasyHap HapGroup.tsv")
    ap.add_argument("--trait-col", default="Plant_height")
    ap.add_argument("--out", default="dosage_vs_height.png")
    ap.add_argument("--top-n", type=int, default=3,
                    help="classes shown in panel (a); the remainder are pooled as Other")
    args = ap.parse_args()

    work = os.environ["WORKDIR"]

    # ---------- markers (VCF ID column) ----------
    markers = []
    for line in open(args.vcf):
        if line.startswith("#"):
            continue
        markers.append(line.rstrip("\n").split("\t")[2])
    markers = list(dict.fromkeys(markers))
    if not markers:
        raise SystemExit("[ERROR] no marker in the VCF")

    # ---------- samples and phenotype ----------
    samples = [l.strip() for l in open(os.environ["SAMPLE_ORDER"]) if l.strip()]
    pheno = {}
    for line in open(os.environ["TRAITS"]):
        p = line.split()
        if len(p) >= 2 and p[0] != "Sample":
            pheno[p[0]] = float(p[1])
    y = np.array([pheno.get(s, np.nan) for s in samples], float)

    # ---------- dosage from BIMBAM ----------
    want = set(markers)
    rows = {}
    for line in open(os.path.join(work, "region_dosage.bimbam")):
        p = line.rstrip("\n").split(", ")
        k = p[0].strip()
        if k in want:
            rows[k] = np.array(
                [np.nan if x in ("NA", "", "nan", "NaN") else float(x) for x in p[3:]],
                float)
    missing = want - set(rows)
    if missing:
        print("[warn] %d k-mer(s) without dosage record, ignored" % len(missing))

    mat = np.vstack([rows[k] for k in markers if k in rows])
    dos = np.round(np.nanmean(mat, axis=0), DECIMALS)   # 0-2 scale, 4 dp
    ok = ~np.isnan(dos) & ~np.isnan(y)

    # ---------- statistics ----------
    rho, p_rho = stats.spearmanr(dos[ok], y[ok])
    fit = stats.linregress(dos[ok], y[ok])
    r2, p_ols = fit.rvalue ** 2, fit.pvalue
    tcrit = stats.t.ppf(0.975, ok.sum() - 2)
    print("[step6] markers: %d   accessions: %d   dosage stored to %d dp"
          % (len(mat), int(ok.sum()), DECIMALS))
    print("[step6] Spearman rho = %+.3f   P = %.2e" % (rho, p_rho))
    print("[step6] OLS  Y = %.1f %+.1f X   R2 = %.3f   P = %.2e   95%% CI [%.1f, %.1f]"
          % (fit.intercept, fit.slope, r2, p_ols,
             fit.slope - tcrit * fit.stderr, fit.slope + tcrit * fit.stderr))

    # ---------- haplotype classes ----------
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6))
    ax = axes[0]
    labels, groups, letters, kw_p = [], [], [], None
    hap_v = None

    if args.hapgroup and os.path.exists(args.hapgroup):
        hap = {}
        with open(args.hapgroup) as fh:
            fh.readline()
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    hap[p[2]] = p[0]
        hap_v = np.array([hap.get(s, "NA") for s in samples])

        counts = {}
        for h in hap_v:
            counts[h] = counts.get(h, 0) + 1
        main_cls = sorted([h for h, _ in sorted(counts.items(), key=lambda kv: -kv[1])
                           if h.startswith("Hap")][:args.top_n])

        for c in main_cls:
            sel = (hap_v == c) & ~np.isnan(y)
            groups.append(y[sel])
            labels.append("%s\n(n = %d)" % (c, int(sel.sum())))
        sel = np.array([h not in main_cls for h in hap_v]) & ~np.isnan(y)
        if sel.sum() > 0:
            groups.append(y[sel])
            labels.append("Other\n(n = %d)" % int(sel.sum()))

        kw = stats.kruskal(*groups)
        kw_p = kw.pvalue
        pairs = list(itertools.combinations(range(len(labels)), 2))
        praw = np.array([mwu(groups[i], groups[j]) for i, j in pairs])
        q = bh_adjust(praw)
        pmat = [[1.0] * len(labels) for _ in labels]
        for (i, j), p, qq in zip(pairs, praw, q):
            pmat[i][j] = pmat[j][i] = qq
        letters = letter_display(labels, pmat)

        overall = float(np.nanmean(y[ok]))
        print("[step6] Kruskal-Wallis H = %.2f, df = %d, P = %.4f"
              % (kw.statistic, len(groups) - 1, kw_p))
        for lab, g, lt in zip(labels, groups, letters):
            print("  %-7s letter=%s  n=%3d  mean=%7.2f  (%+5.2f%%)"
                  % (lab.splitlines()[0], lt, g.size, g.mean(),
                     100.0 * (g.mean() - overall) / overall))
        for (i, j), p, qq in zip(pairs, praw, q):
            print("  %s vs %s  P = %.5f  BH = %.4f  %s"
                  % (labels[i].splitlines()[0], labels[j].splitlines()[0], p, qq,
                     "*" if qq < 0.05 else ""))

        # ---------- panel (a) ----------
        bp = ax.boxplot(groups, patch_artist=True, showfliers=False, widths=0.6,
                        medianprops=dict(color="#1F3864", lw=1.4),
                        whiskerprops=dict(color="#4472C4"),
                        capprops=dict(color="#4472C4"))
        colors = ["#F4B183", "#A9D18E", "#9DC3E6", "#D9D9D9", "#FFE699"]
        for b, c in zip(bp["boxes"], colors):
            b.set_facecolor(c)
            b.set_edgecolor("#7F7F7F")
        rng = np.random.RandomState(7)
        for i, g in enumerate(groups):
            x = (i + 1) + rng.uniform(-0.15, 0.15, g.size)
            ax.scatter(x, g, s=14, color="#404040", alpha=0.55, linewidths=0, zorder=3)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("%s BLUP" % args.trait_col)
        ax.set_title("(a) Haplotype classes   KW P = %.3f" % kw_p, fontsize=10)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        top = max(np.nanmax(g) for g in groups)
        for i, lt in enumerate(letters):
            ax.text(i + 1, top + 8, lt, ha="center", va="bottom",
                    fontsize=11, fontweight="bold")
        ax.set_ylim(min(np.nanmin(g) for g in groups) - 10, top + 26)

    # ---------- panel (b) ----------
    ax2 = axes[1]
    ax2.scatter(dos[ok], y[ok], s=18, color="#4472C4", alpha=0.7, linewidths=0)
    xs = np.linspace(np.nanmin(dos[ok]), np.nanmax(dos[ok]), 100)
    fitted = fit.intercept + fit.slope * xs
    resid = y[ok] - (fit.intercept + fit.slope * dos[ok])
    sxx = ((dos[ok] - dos[ok].mean()) ** 2).sum()
    n = int(ok.sum())
    se_fit = np.sqrt((resid ** 2).sum() / (n - 2)
                     * (1.0 / n + (xs - dos[ok].mean()) ** 2 / sxx))
    ax2.plot(xs, fitted, color="#C00000", lw=1.5)
    ax2.fill_between(xs, fitted - tcrit * se_fit, fitted + tcrit * se_fit,
                     color="#C00000", alpha=0.13, linewidth=0)
    ax2.set_xlabel("Mean k-mer dosage (%d k-mers, 0-2)" % len(mat))
    ax2.set_ylabel("%s BLUP" % args.trait_col)
    ax2.set_title("(b) Dosage vs %s" % args.trait_col.replace("_", " "), fontsize=10)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    txt = ("Y = %.1f %s %.1f X\n"
           "$R^2$ = %.3f,  P = %.1e\n"
           "slope 95%% CI [%.1f, %.1f]\n"
           "Spearman $\\rho$ = %.3f, P = %.1e\n"
           "n = %d" % (fit.intercept, "-" if fit.slope < 0 else "+", abs(fit.slope),
                       r2, p_ols, fit.slope - tcrit * fit.stderr,
                       fit.slope + tcrit * fit.stderr, rho, p_rho, n))
    ax2.text(0.97, 0.96, txt, transform=ax2.transAxes, ha="right", va="top",
             fontsize=8.5, linespacing=1.5,
             bbox=dict(facecolor="white", edgecolor="#D9D9D9", alpha=0.9,
                       boxstyle="round,pad=0.35"))

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close(fig)
    print("[step6] written: " + args.out)

    # ---------- per-accession table (supplementary table source) ----------
    table = os.path.join(work, "dosage_by_accession.tsv")
    with open(table, "w") as fh:
        fh.write("Sample\tHaplotype\tMeanDosage\t%s\n" % args.trait_col)
        for i, s in enumerate(samples):
            fh.write("%s\t%s\t%.4f\t%s\n"
                     % (s, hap_v[i] if hap_v is not None else "NA",
                        dos[i], ("%.6f" % y[i]) if not np.isnan(y[i]) else "NA"))
    print("[step6] written: " + table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
