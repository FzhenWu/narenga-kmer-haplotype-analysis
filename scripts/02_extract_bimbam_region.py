#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步骤 2：从 BIMBAM 分块文件中提取步骤 1 得到的 k-mer 的剂量行。

BIMBAM 每行格式（逗号 + 空格分隔）：
    kmer_sequence, X, Y, sample1_dose, sample2_dose, ..., sample112_dose

输入：
    ${WORKDIR}/patterns.txt
    ${BIMBAM_DIR}/filtered_sample_k31.*.bimbam.gz
输出：
    ${WORKDIR}/region_dosage.bimbam   （按坐标顺序、每条 k-mer 只保留一行）

环境变量（见 config.sh）：
    WORKDIR, BIMBAM_DIR
"""
import glob
import gzip
import os
import sys


def main() -> int:
    work = os.environ["WORKDIR"]
    bimbam_dir = os.environ["BIMBAM_DIR"]

    pat = os.path.join(work, "patterns.txt")
    want = [l.strip() for l in open(pat) if l.strip()]
    wantset = set(want)

    files = sorted(glob.glob(os.path.join(bimbam_dir, "filtered_sample_k31.*.bimbam.gz")))
    if not files:
        sys.stderr.write(f"[ERROR] {bimbam_dir} 下没有 BIMBAM 分块文件\n")
        return 1

    seen = {}
    for i, f in enumerate(files, 1):
        n_before = len(seen)
        with gzip.open(f, "rt") as fh:
            for line in fh:
                k = line.split(",", 1)[0].strip()
                if k in wantset and k not in seen:
                    seen[k] = line if line.endswith("\n") else line + "\n"
        print(f"[step2] {i}/{len(files)} {os.path.basename(f)}  "
              f"命中累计 {len(seen)} (+{len(seen) - n_before})", flush=True)

    out = os.path.join(work, "region_dosage.bimbam")
    with open(out, "w") as w:
        for k in want:                    # 按坐标顺序输出
            if k in seen:
                w.write(seen[k])

    missing = [k for k in want if k not in seen]
    print(f"[step2] 查询 k-mer: {len(want)}  命中: {len(seen)}  缺失: {len(missing)}")
    if missing:
        print(f"[step2] 缺失示例: {missing[:5]}")
    print(f"[step2] 写出: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
