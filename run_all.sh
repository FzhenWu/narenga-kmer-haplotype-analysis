#!/usr/bin/env bash
# ============================================================================
# End-to-end run of the main analysis (Bonferroni set as primary; the
# region-significant set is also produced for comparison).
#
#     bash run_all.sh
#
# Prerequisite: EasyHap (>= 1.4.0) installed and paths in config.sh edited.
# ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/config.sh"

PY="${PYTHON:-python}"
cd "${WORKDIR}"

echo "==================== step 0: traits ===================="
"$PY" "${HERE}/scripts/00_prepare_traits.py"

echo "==================== step 1: region k-mers ============="
"$PY" "${HERE}/scripts/01_extract_region_kmers.py"

echo "==================== step 2: BIMBAM dosages ============"
"$PY" "${HERE}/scripts/02_extract_bimbam_region.py"

echo "==================== step 3: significant k-mers ========"
"$PY" "${HERE}/scripts/03_select_significant_kmers.py"

echo "==================== step 4: build VCF ================="
"$PY" "${HERE}/scripts/04_build_vcf.py" --markers bonf   --ploidy 2 --out region_bonf.vcf
"$PY" "${HERE}/scripts/04_build_vcf.py" --markers sig    --ploidy 2 --out region_sig.vcf
"$PY" "${HERE}/scripts/04_build_vcf.py" --markers region --ploidy 2 --out region_all.vcf
"$PY" "${HERE}/scripts/04_build_vcf.py" --markers sig    --ploidy 4 --out region_sig_p4.vcf

echo "==================== step 5: EasyHap ==================="
bash "${HERE}/scripts/05_run_easyhap.sh" region_bonf.vcf easyhap_bonf
bash "${HERE}/scripts/05_run_easyhap.sh" region_sig.vcf  easyhap_sig

echo "==================== step 6: dosage plot ==============="
"$PY" "${HERE}/scripts/06_plot_dosage.py" \
    --vcf region_bonf.vcf \
    --hapgroup "easyhap_bonf/${CHROM}_${REGION_START}_${REGION_END}.HapGroup.tsv" \
    --out dosage_vs_height_bonf.png

echo "==================== step 7: Fig. 4D =================="
"$PY" "${HERE}/scripts/07_gene_haplotype_figure.py" \
    --tokens "easyhap_bonf/${CHROM}_${REGION_START}_${REGION_END}.SampleGenotypeTokens.tsv" \
    --hapgroup "easyhap_bonf/${CHROM}_${REGION_START}_${REGION_END}.HapGroup.tsv" \
    --features "${HERE}/data/gene_features_04Gg010100.tsv" \
    --out-prefix "${HERE}/figures/Fig4D_gene_haplotype"

echo "==================== step 8: Table S16 ================="
"$PY" "${HERE}/scripts/08_table_s16.py" \
    --input "${WORKDIR}/dosage_by_accession.tsv" \
    --out-prefix "${HERE}/supplementary/Table_S16"

echo
echo "Done. Results are in: ${WORKDIR}"
ls -1d easyhap_* dosage_vs_height_bonf.png dosage_by_accession.tsv 2>/dev/null || true
