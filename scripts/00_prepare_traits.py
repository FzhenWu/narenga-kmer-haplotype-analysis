#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 0: Build traits.tsv (EasyHap input) from the phenotype file and check
the sample order.

EasyHap requires that the first column of the trait file contains sample IDs
matching the VCF sample columns. This script writes the file in the BIMBAM
column order (i.e. sample_order.txt), guaranteeing exact alignment.

Environment variables (see config.sh):
    WORKDIR, PHENO_FILE, SAMPLE_ORDER
"""
import os
import sys


def main() -> int:
    work = os.environ["WORKDIR"]
    pheno_file = os.environ["PHENO_FILE"]
    sample_file = os.environ["SAMPLE_ORDER"]

    samples = [l.strip() for l in open(sample_file) if l.strip()]

    pheno = {}
    pheno_order = []
    for line in open(pheno_file):
        p = line.split()
        if len(p) < 2:
            continue
        pheno[p[0]] = float(p[1])
        pheno_order.append(p[0])

    missing = [s for s in samples if s not in pheno]
    if missing:
        sys.stderr.write(
            f"[ERROR] {len(missing)} samples are missing from the phenotype file, "
            f"e.g. {missing[:5]}\n")
        return 1

    out = os.path.join(work, "traits.tsv")
    with open(out, "w") as w:
        w.write("Sample\tPlant_height\n")
        for s in samples:                 # strictly follow the BIMBAM column order
            w.write(f"{s}\t{pheno[s]:.6f}\n")

    print(f"[step0] traits.tsv written: {out}")
    print(f"[step0] samples: {len(samples)}")
    print(f"[step0] phenotype row order == BIMBAM column order: {pheno_order == samples}")

    if pheno_order != samples:
        print("[step0] NOTE: orders differ; the output was reordered to match "
              "sample_order.txt and is therefore aligned.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
