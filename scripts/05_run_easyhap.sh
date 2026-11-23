#!/usr/bin/env bash
# ============================================================================
# 步骤 5：用 EasyHap 做区域单倍型分析并生成图件
#
# 用法：
#   bash scripts/05_run_easyhap.sh <input.vcf> <outdir> [额外 EasyHap 参数...]
#
# 例：
#   bash scripts/05_run_easyhap.sh region_bonf.vcf easyhap_bonf
#   bash scripts/05_run_easyhap.sh region_sig_p4.vcf easyhap_sig_p4 --cluster-threshold 0.3
# ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

if [[ $# -lt 2 ]]; then
    echo "用法: bash $0 <input.vcf> <outdir> [extra args...]" >&2
    exit 1
fi

VCF="$1"
OUTDIR="$2"
shift 2

if [[ ! -f "$VCF" ]]; then
    echo "[ERROR] 找不到 VCF: $VCF" >&2
    exit 1
fi

# ---------- 截取目标染色体的 GFF 子集（若尚未生成） ----------
GFF_SUB="${WORKDIR}/${CHROM}.gff3"
if [[ ! -f "$GFF_SUB" ]]; then
    echo "[step5] 生成 GFF 子集: $GFF_SUB"
    { grep -E '^#' "$GFF3"; awk -F'\t' -v c="$CHROM" '$1==c' "$GFF3"; } > "$GFF_SUB"
fi

GFF_ARG=()
if [[ -s "$GFF_SUB" ]]; then
    GFF_ARG=(--gff "$GFF_SUB" --gene-feature exon)
else
    echo "[step5] 警告：GFF 子集为空，跳过基因结构注释"
fi

echo "[step5] easyhap analyze  ${VCF}  ->  ${OUTDIR}"

easyhap analyze \
    --vcf "$VCF" \
    --region "${CHROM}:${REGION_START}-${REGION_END}" \
    --mode genotype \
    --traits "${TRAITS}" \
    --trait-cols Plant_height \
    ${GFF_ARG[@]+"${GFF_ARG[@]}"} \
    --plot --plot-format png \
    --outdir "$OUTDIR" \
    "$@"

echo "[step5] 完成: $OUTDIR"
