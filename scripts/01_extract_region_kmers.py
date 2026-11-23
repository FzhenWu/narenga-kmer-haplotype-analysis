#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步骤 1：从 k-mer 坐标索引（sqlite）中提取目标区间内 MAPQ>20 的 k-mer。

输入：
    ${KMER_IDX_DIR}/kmer_mapq_gt20.sqlite   （表 positions(kmer, chr, pos)）
输出：
    ${WORKDIR}/region_kmers.tsv   列：kmer / chr / start / end
    ${WORKDIR}/patterns.txt       每行一条 k-mer 序列，供步骤 2 使用

环境变量（见 config.sh）：
    KMER_IDX_DIR, WORKDIR, CHROM, REGION_START, REGION_END, KMER_LEN
"""
import os
import sqlite3
import sys


def main() -> int:
    idx = os.environ["KMER_IDX_DIR"]
    work = os.environ["WORKDIR"]
    chrom = os.environ["CHROM"]
    start = int(os.environ["REGION_START"])
    end = int(os.environ["REGION_END"])
    klen = int(os.environ.get("KMER_LEN", 31))

    db = os.path.join(idx, "kmer_mapq_gt20.sqlite")
    if not os.path.exists(db):
        sys.stderr.write(f"[ERROR] 找不到坐标索引: {db}\n")
        return 1

    con = sqlite3.connect(db)
    cur = con.cursor()
    rows = cur.execute(
        "select kmer, chr, cast(pos as integer) "
        "from positions "
        "where chr = ? and cast(pos as integer) between ? and ? "
        "order by 3",
        (chrom, start, end),
    ).fetchall()
    con.close()

    if not rows:
        sys.stderr.write(
            f"[ERROR] 区间 {chrom}:{start}-{end} 内没有 k-mer，请检查染色体命名。\n")
        return 1

    tsv = os.path.join(work, "region_kmers.tsv")
    with open(tsv, "w") as w:
        w.write("kmer\tchr\tstart\tend\n")
        for k, c, p in rows:
            w.write(f"{k}\t{c}\t{p}\t{p + klen - 1}\n")

    pat = os.path.join(work, "patterns.txt")
    with open(pat, "w") as w:
        w.write("\n".join(k for k, _, _ in rows) + "\n")

    span = rows[-1][2] - rows[0][2] + 1
    print(f"[step1] 区间: {chrom}:{start}-{end}  (长度 {end - start + 1} bp)")
    print(f"[step1] k-mer 数: {len(rows)}")
    print(f"[step1] 坐标跨度: {rows[0][2]} - {rows[-1][2]}  ({span} bp)")
    print(f"[step1] 写出: {tsv}, {pat}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
