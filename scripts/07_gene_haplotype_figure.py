#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 7: Fig. 4D -- gene model plus P/A haplotype heatmap.

Cells are labelled with k-mer presence/absence states (A = absent, P = present)
instead of the placeholder VCF alleles (A/A, A/G, G/G). Columns follow the
5'->3' orientation of the minus-strand gene, so genomic coordinates decrease
from left to right.

Inputs (defaults are in data/):
    sample_genotype_tokens_bonf.tsv   EasyHap SampleGenotypeTokens.tsv
    hapgroups_bonf.tsv                EasyHap HapGroup.tsv
    gene_features_04Gg010100.tsv      start/end/type (UTR, CDS); introns implied

Outputs (default prefix figures/Fig4D_gene_haplotype):
    *_PA.png / *_PA.pdf               binary states (Fig. 4D)
    *_DOSAGE.png                      optional dosage classes A / P/A / P

Usage:
    python scripts/07_gene_haplotype_figure.py
    python scripts/07_gene_haplotype_figure.py --dosage
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

TOKEN2DOS = {"A/A": 0, "A/G": 1, "G/G": 2, "G/A": 1}


def load_tokens(path):
    fh = open(path)
    header = fh.readline().rstrip("\n").split("\t")
    samples = header[3:]
    pos, dos = [], []
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        pos.append(int(p[1]))
        dos.append([TOKEN2DOS.get(t, np.nan) for t in p[3:]])
    fh.close()
    dos = np.array(dos, float)

    # minus strand: 5' -> 3' means descending coordinate
    order = np.argsort(pos)[::-1]
    pos = [pos[i] for i in order]
    dos = dos[order, :]
    return samples, pos, dos


def load_hapgroup(path):
    hap_of = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                hap_of[p[2]] = p[0]
    return hap_of


def load_features(path):
    feats = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                feats.append((int(p[0]), int(p[1]), p[2]))
    return feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tokens", default=os.path.join(DATA, "sample_genotype_tokens_bonf.tsv"))
    ap.add_argument("--hapgroup", default=os.path.join(DATA, "hapgroup_bonf.tsv"))
    ap.add_argument("--features", default=os.path.join(DATA, "gene_features_04Gg010100.tsv"))
    ap.add_argument("--gene", default="SoffiXsponR570.04Gg010100")
    ap.add_argument("--out-prefix", default=os.path.join(ROOT, "figures", "Fig4D_gene_haplotype"))
    ap.add_argument("--region-label", default="Chr4G:3213203-3231114")
    ap.add_argument("--dosage", action="store_true",
                    help="also write the A / P/A / P version")
    args = ap.parse_args()

    samples, pos, dos = load_tokens(args.tokens)
    hap_of = load_hapgroup(args.hapgroup)
    feats = load_features(args.features)

    haps = sorted(set(hap_of.values()))
    states, counts = [], []
    for h in haps:
        idx = [i for i, s in enumerate(samples) if hap_of.get(s) == h]
        counts.append(len(idx))
        row = []
        for j in range(len(pos)):
            v = dos[j, idx]
            v = v[~np.isnan(v)]
            if v.size == 0:
                row.append(np.nan)
            else:
                u, c = np.unique(v, return_counts=True)
                row.append(u[np.argmax(c)])
        states.append(row)
    states = np.array(states, float)

    o = np.lexsort((np.array(haps), -np.array(counts)))
    haps = [haps[i] for i in o]
    states = states[o]
    counts = [counts[i] for i in o]

    binary = np.where(np.isnan(states), np.nan, (states > 0).astype(float))

    gmin = min(f[0] for f in feats)
    gmax = max(f[1] for f in feats)
    pad = int((gmax - gmin) * 0.03)
    nr, nc = states.shape

    fig = plt.figure(figsize=(9.2, 6.0))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 3.2], hspace=0.22,
                           left=0.11, right=0.72, top=0.83, bottom=0.17)
    axg = fig.add_subplot(gs[0])
    axh = fig.add_subplot(gs[1])

    # ---------- gene model ----------
    axg.set_xlim(gmax + pad, gmin - pad)      # high coordinate on the left = 5'
    axg.set_ylim(0, 1)
    axg.axis("off")
    ybar, h_cds, h_utr = 0.42, 0.30, 0.225
    axg.plot([gmin, gmax], [ybar, ybar], color="#333333", lw=1.0,
             solid_capstyle="butt")
    for a, b, kind in feats:
        hh = h_cds if kind == "CDS" else h_utr
        col = "#1A1A1A" if kind == "CDS" else "#BFBFBF"
        axg.add_patch(Rectangle((a, ybar - hh / 2), b - a, hh,
                                facecolor=col, edgecolor="none"))
    axg.text((gmin + gmax) / 2.0, 0.90, args.gene, ha="center", va="center",
             fontsize=10)
    axg.text(gmax + pad, 0.42, "5'", ha="left", va="center", fontsize=10)
    axg.text(gmin - pad, 0.42, "3'", ha="right", va="center", fontsize=10)
    axg.legend(handles=[Patch(facecolor="#BFBFBF", label="UTR"),
                        Patch(facecolor="#1A1A1A", label="CDS"),
                        Line2D([0], [0], color="#333333", lw=1.0, label="Intron")],
               loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3,
               frameon=False, fontsize=8.5, handlelength=1.6)

    COL_AB = {0.0: "#A9D18E", 1.0: "#9DC3E6"}
    LAB_AB = {0.0: "A", 1.0: "P"}
    COL_DOS = {0.0: "#A9D18E", 1.0: "#F4B183", 2.0: "#2E75B6"}
    LAB_DOS = {0.0: "A", 1.0: "P/A", 2.0: "P"}

    def draw(mat, colmap, labmap, fname, items, suffix, pdf=False):
        axh.cla()
        fig.lines = []
        fig.suptitle("Gene structure and haplotypes: %s%s" % (args.region_label, suffix),
                     fontsize=12, y=0.97)
        for i in range(nr):
            for j in range(nc):
                v = mat[i, j]
                if np.isnan(v):
                    axh.add_patch(Rectangle((j, i), 1, 1, facecolor="white",
                                            edgecolor="#BFBFBF", lw=0.6))
                    continue
                axh.add_patch(Rectangle((j, i), 1, 1, facecolor=colmap[float(v)],
                                        edgecolor="#7F7F7F", lw=0.8))
                axh.text(j + 0.5, i + 0.5, labmap[float(v)], ha="center",
                         va="center", fontsize=10, color="#1A1A1A")
        axh.set_xlim(0, nc)
        axh.set_ylim(nr, 0)
        axh.set_xticks(np.arange(nc) + 0.5)
        axh.set_xticklabels([str(p) for p in pos], rotation=90, fontsize=8.5)
        axh.set_yticks(np.arange(nr) + 0.5)
        axh.set_yticklabels(haps, fontsize=9)
        axh.set_xlabel("Variant position", fontsize=10)
        axh.set_ylabel("Haplotype class", fontsize=10)
        for s in ("top", "right"):
            axh.spines[s].set_visible(False)
        axh.tick_params(length=0)
        for k, n in enumerate(counts):
            axh.text(nc + 0.15, k + 0.5, str(n), va="center", fontsize=8.5,
                     clip_on=False)
        axh.legend(handles=[Patch(facecolor=colmap[v], edgecolor="#7F7F7F",
                                  label=labmap[v] + "  " + t) for v, t in items],
                   title="State", loc="upper left", bbox_to_anchor=(1.12, 1.0),
                   frameon=False, fontsize=8.5)

        # connectors: heatmap column -> genomic position
        fig.canvas.draw()
        inv = fig.transFigure.inverted()
        for j, p in enumerate(pos):
            x0, y0 = inv.transform(axh.transData.transform((j + 0.5, 0)))
            x1, y1 = inv.transform(axg.transData.transform((p, ybar)))
            fig.lines.append(Line2D([x0, x1], [y0 - 0.005, y1], color="#9E9E9E",
                                    lw=0.6, transform=fig.transFigure, zorder=0))
        fig.savefig(fname, dpi=300)
        print("saved: " + fname)
        if pdf:
            out_pdf = fname[:-4] + ".pdf"
            fig.savefig(out_pdf)
            print("saved: " + out_pdf)

    draw(binary, COL_AB, LAB_AB, args.out_prefix + "_PA.png",
         [(0.0, "(absent, 0 copies)"), (1.0, "(present, >=1 copy)")],
         " (P/A)", pdf=True)
    if args.dosage:
        draw(states, COL_DOS, LAB_DOS, args.out_prefix + "_DOSAGE.png",
             [(0.0, "(absent, 0 copies)"), (1.0, "(1 copy)"), (2.0, "(2 copies)")],
             " (P/A with dosage)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
