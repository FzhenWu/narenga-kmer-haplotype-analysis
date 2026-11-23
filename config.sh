#!/usr/bin/env bash

# ============================================================================
# Global configuration
#
# All paths below are placeholders. Edit them to match your own environment,
# then load this file before running any step:
#
#     source config.sh
# ============================================================================

# ---------- Project layout ----------
# Root of the KMERIA output directory (one sub-directory per stage).
export KMER_ROOT="/path/to/your/project/kmer_gwas"
export ASSOC_DIR="${KMER_ROOT}/05_association"
export KMER_IDX_DIR="${ASSOC_DIR}/02.minibwa/R570"   # sqlite index / positions / GEMMA output
export BIMBAM_DIR="${KMER_ROOT}/04_bimbam"           # filtered_sample_k31.*.bimbam.gz

# ---------- Input files ----------
export SAMPLE_ORDER="${KMER_ROOT}/01_kmer_counts/sample_order.txt"
export PHENO_FILE="${KMER_ROOT}/phenotype_Plant_height.txt"
export REF_FA="/path/to/reference/SofficinarumxspontaneumR570_771_v2.0.fa"
export GFF3="/path/to/reference/R570.v2023.gff3"
export GEMMA_SIG="${KMER_IDX_DIR}/Plant_height.significant.MAPQ_gt20.txt"
export GEMMA_ALL="${KMER_IDX_DIR}/Plant_height.merged_assoc.MAPQ_gt20.txt"

# ---------- Target region ----------
export CHROM="Chr4G"
export REGION_START=3213203          # gene +/- 5 kb
export REGION_END=3231114
export GENE_START=3218203
export GENE_END=3226114
export GENE_ID="SoffiXsponR570.04Gg010100"
export KMER_LEN=31

# ---------- Statistical parameters ----------
export N_TESTS=34786919              # distinct 31-mers tested per trait (Bonferroni family size)
export ALPHA=0.05
export PLOIDY=2                      # main analysis: 0/1/2 coding; use 4 for the sensitivity check

# ---------- Output directory ----------
export WORKDIR="${KMER_IDX_DIR}/hap04Gg010100"
export TRAITS="${WORKDIR}/traits.tsv"

mkdir -p "${WORKDIR}"

echo "[config] WORKDIR = ${WORKDIR}"
