#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步骤 3：从 GEMMA 全基因组关联结果中筛出目标区间内的显著 k-mer，
并按 Bonferroni 校正标注核心位点。

输入：
    ${GEMMA_SIG}  列：chr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald
输出：
    ${WORKDIR}/sig_region.tsv    区间内名义显著 k-mer：kmer / pos / af / beta / p_wald
    ${WORKDIR}/bonf_region.tsv   其中通过 Bonferroni 校正的 k-mer
    ${WORKDIR}/marker_list.txt   供步骤 4 使用，可用 --set 选择

环境变量（见 config.sh）：
    WORKDIR, GEMMA_SIG, CHROM, REGION_START, REGION_END, N_TESTS, ALPHA
"""
import os
import sys


def main() -> int:
    work = os.environ["WORKDIR"]
    sig_file = os.environ["GEMMA_SIG"]
    chrom = os.environ["CHROM"]
    start = int(os.environ["REGION_START"])
    end = int(os.environ["REGION_END"])
    n_tests = int(os.environ["N_TESTS"])
    alpha = float(os.environ.get("ALPHA", 0.05))

    bonf = alpha / n_tests
    print(f"[step3] Bonferroni 阈值 = {alpha} / {n_tests} = {bonf:.3e}")

    header = None
    rows = []
    for i, line in enumerate(open(sig_file)):
        p = line.rstrip("\n").split("\t")
        if i == 0:
            header = [x.strip() for x in p]
            continue
        if len(p) < len(header):
            continue
        d = dict(zip(header, p))
        try:
            pos = int(float(d["ps"]))
        except (KeyError, ValueError):
            continue
        if d.get("chr", "").strip() != chrom:
            continue
        if not (start <= pos <= end):
            continue
        rows.append({
            "kmer": d["rs"].strip(),
            "pos": pos,
            "af": d.get("af", "NA").strip(),
            "beta": d.get("beta", "NA").strip(),
            "p": float(d["p_wald"]),
        })

    rows.sort(key=lambda r: r["pos"])
    sig_rows = rows
    bonf_rows = [r for r in rows if r["p"] < bonf]

    with open(os.path.join(work, "sig_region.tsv"), "w") as w:
        w.write("kmer\tpos\taf\tbeta\tp_wald\n")
        for r in sig_rows:
            w.write(f"{r['kmer']}\t{r['pos']}\t{r['af']}\t{r['beta']}\t{r['p']:.6e}\n")

    with open(os.path.join(work, "bonf_region.tsv"), "w") as w:
        w.write("kmer\tpos\taf\tbeta\tp_wald\n")
        for r in bonf_rows:
            w.write(f"{r['kmer']}\t{r['pos']}\t{r['af']}\t{r['beta']}\t{r['p']:.6e}\n")

    # 供步骤 4 直接使用的标记清单
    for name, rr in (("region", rows), ("sig", sig_rows), ("bonf", bonf_rows)):
        with open(os.path.join(work, f"markers_{name}.txt"), "w") as w:
            w.write("\n".join(r["kmer"] for r in rr) + ("\n" if rr else ""))

    if rows:
        span = rows[-1]["pos"] - rows[0]["pos"] + 1
        print(f"[step3] 区间内名义显著 k-mer: {len(sig_rows)}")
        print(f"[step3]   P 范围: {min(r['p'] for r in sig_rows):.2e} ~ "
              f"{max(r['p'] for r in sig_rows):.2e}")
        print(f"[step3]   起点跨度: {rows[0]['pos']} - {rows[-1]['pos']}  ({span} bp)")
        print(f"[step3] Bonferroni 显著: {len(bonf_rows)}")
    else:
        print("[step3] 警告：区间内没有显著 k-mer")

    return 0


if __name__ == "__main__":
    sys.exit(main())
