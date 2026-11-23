#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步骤 4：将 k-mer 剂量矩阵转换为 EasyHap 可读的伪 VCF。

说明：
    每条 k-mer 视为一个双等位标记：
        CHROM = 比对染色体，POS = 比对起点，ID = k-mer 序列
        REF = A（占位），ALT = G（占位）
        FORMAT = GT:DS
        GT 由剂量取整得到（0 = 缺失 / 1 = 单拷贝 / 2 = 双拷贝）
        DS 保留原始剂量（连续值）
    注意：REF/ALT 为占位符，无生物学含义，仅用于满足 VCF 规范。

用法：
    python 04_build_vcf.py --markers bonf --ploidy 2 --out region_bonf.vcf
    python 04_build_vcf.py --markers sig  --ploidy 4 --out region_sig_p4.vcf

环境变量（见 config.sh）：
    WORKDIR, SAMPLE_ORDER, CHROM
"""
import argparse
import collections
import os
import sys


def load_positions(path):
    pos = {}
    with open(path) as fh:
        fh.readline()                      # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                pos[p[0]] = int(p[2])
    return pos


def load_dosages(path, keep):
    rows = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split(", ")
            k = p[0].strip()
            if k not in keep:
                continue
            vals = []
            for x in p[3:]:
                vals.append(float("nan") if x in ("NA", "", "nan", "NaN") else float(x))
            rows[k] = vals
    return rows


def gt_string(copies, ploidy):
    return "/".join(["0"] * (ploidy - copies) + ["1"] * copies)


def main() -> int:
    ap = argparse.ArgumentParser(description="k-mer dosage -> pseudo VCF (EasyHap input)")
    ap.add_argument("--markers", choices=["region", "sig", "bonf"], default="bonf",
                    help="使用哪一组标记（对应 step3 生成的 markers_*.txt）")
    ap.add_argument("--marker-file", default=None, help="自定义标记清单（每行一条 k-mer）")
    ap.add_argument("--ploidy", type=int, choices=[2, 4], default=2,
                    help="2 = 0/1/2 二体编码；4 = 0-4 四倍体编码（剂量 x2 后取整）")
    ap.add_argument("--out", required=True, help="输出 VCF 路径")
    args = ap.parse_args()

    work = os.environ["WORKDIR"]
    chrom = os.environ["CHROM"]

    samples = [l.strip() for l in open(os.environ["SAMPLE_ORDER"]) if l.strip()]
    pos = load_positions(os.path.join(work, "region_kmers.tsv"))

    if args.marker_file:
        keep = [l.strip() for l in open(args.marker_file) if l.strip()]
    else:
        keep = [l.strip() for l in open(os.path.join(work, f"markers_{args.markers}.txt"))
                if l.strip()]
    keep = list(dict.fromkeys(keep))
    if not keep:
        sys.stderr.write("[ERROR] 标记清单为空\n")
        return 1

    dos = load_dosages(os.path.join(work, "region_dosage.bimbam"), set(keep))

    scale = 2.0 if args.ploidy == 4 else 1.0
    order = sorted((k for k in keep if k in dos), key=lambda k: pos.get(k, 1 << 60))

    hist = collections.Counter()
    with open(args.out, "w") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write(f"##contig=<ID={chrom}>\n")
        w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (k-mer dosage)">\n')
        w.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Raw k-mer dosage">\n')
        w.write('##INFO=<ID=KMER_LEN,Number=1,Type=Integer,Description="k-mer length">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")

        for k in order:
            gts = []
            for v in dos[k]:
                if v != v:                       # NaN
                    gts.append(":".join(["/".join(["."] * args.ploidy), "."]))
                    continue
                c = int(round(v * scale))
                c = max(0, min(args.ploidy, c))
                hist[c] += 1
                gts.append(f"{gt_string(c, args.ploidy)}:{v:.4f}")
            w.write(f"{chrom}\t{pos.get(k, 1)}\t{k}\tA\tG\t.\t.\t.\tGT:DS\t"
                    + "\t".join(gts) + "\n")

    print(f"[step4] 位点数: {len(order)}   样品数: {len(samples)}   ploidy: {args.ploidy}")
    print(f"[step4] 剂量分布(取整后拷贝数 → 单元数): {dict(sorted(hist.items()))}")
    print(f"[step4] 写出: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
