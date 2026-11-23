# K-mer haplotype analysis of a plant-height candidate locus in *Narenga porphyrocoma*

Regional haplotype analysis of a plant-height candidate gene
**SoffiXsponR570.04Gg010100** (Chr4G:3,218,203–3,226,114, minus strand,
Kinesin motor domain) in 112 *Narenga porphyrocoma* accessions
(allotetraploid, 2n = 4x = 30), based on k-mer GWAS results.
Figures are produced with [EasyHap](https://github.com/Guangqi-He/EasyHap).

This repository contains the complete code path from the k-mer coordinate index
and BIMBAM dosage matrix to the pseudo-VCF, haplotype statistics and
publication-quality figures. All results reported in the manuscript can be
reproduced with the scripts below.

---

## 1. Workflow

```
KMC 31-mer counting
   └─ filtering → 03_filtered_matrices/*.bin
        └─ kmeria m2b → 04_bimbam/*.bimbam.gz       (k-mer x sample dosage matrix)
                        │
k-mer alignment to the reference (minimap2/minibwa, MAPQ > 20)
   └─ positions.tsv → kmer_mapq_gt20.sqlite           (k-mer → chr/pos)
                        │
        ┌───────────────┴────────────────┐
        │  Step 1: SQL query of the      │
        │          target interval       │
        └───────────────┬────────────────┘
                        │ region_kmers.tsv / patterns.txt
        ┌───────────────┴────────────────┐
        │  Step 2: extract dosage rows   │
        │          from BIMBAM shards    │
        └───────────────┬────────────────┘
                        │ region_dosage.bimbam
        ┌───────────────┴────────────────┐
        │  Step 3: select GEMMA-significant│
        │          k-mers + Bonferroni    │
        └───────────────┬────────────────┘
                        │ sig_region.tsv / bonf_region.tsv
        ┌───────────────┴────────────────┐
        │  Step 4: dosage → pseudo-VCF   │
        └───────────────┬────────────────┘
                        │ region_*.vcf
        ┌───────────────┴────────────────┐
        │  Step 5: EasyHap haplotype     │
        │          analysis and figures  │
        └───────────────┬────────────────┘
                        │ HapGroup / TraitAssociationTests / *.png
        ┌───────────────┴────────────────┐
        │  Step 6: dosage–phenotype plot │
        └────────────────────────────────┘
```

> **Important note.** A k-mer is a 31-bp sequence, not a genomic variant.
> In the generated VCF, `REF`/`ALT` (`A`/`G`) are placeholders only, and
> `FORMAT/GT` encodes the estimated copy number of that k-mer in each sample
> (0 = absent, 1 = one copy, 2 = two copies). The VCF is a format adapter that
> allows EasyHap to read a dosage matrix; it is not a catalogue of real variants.

---

## 2. Requirements

| Software | Version | Purpose |
|---|---|---|
| Python | ≥ 3.9 | main pipeline |
| numpy / scipy / matplotlib | ≥ 1.24 / 1.10 / 3.6 | statistics and plotting |
| EasyHap | ≥ 1.4.0 | haplotype analysis and figures |
| pysam / cyvcf2 | any | EasyHap dependencies |

```bash
pip install -r requirements.txt

# EasyHap (a dedicated environment is recommended)
git clone https://github.com/Guangqi-He/EasyHap.git
cd EasyHap && pip install -e . --no-deps
```

---

## 3. Quick start

### 3.0 Configure paths

Edit `config.sh` (all paths are placeholders), then load it:

```bash
source config.sh
```

### 3.1 Run step by step

```bash
# 0) Build traits.tsv for EasyHap and check the sample order
python scripts/00_prepare_traits.py

# 1) Extract k-mers located in the target interval (Chr4G:3213203-3231114)
python scripts/01_extract_region_kmers.py

# 2) Extract the dosage rows of these k-mers from the 50 BIMBAM shards
python scripts/02_extract_bimbam_region.py

# 3) Select GEMMA-significant k-mers (nominal and Bonferroni)
python scripts/03_select_significant_kmers.py

# 4) Build VCFs (different marker sets and ploidy codings)
python scripts/04_build_vcf.py --markers bonf --ploidy 2 --out region_bonf.vcf
python scripts/04_build_vcf.py --markers sig  --ploidy 2 --out region_sig.vcf
python scripts/04_build_vcf.py --markers sig  --ploidy 4 --out region_sig_p4.vcf

# 5) Subset the GFF for Chr4G and run EasyHap
bash scripts/05_run_easyhap.sh region_bonf.vcf easyhap_bonf

# 6) Dosage–phenotype figure (letters, OLS fit, 95% CI) and the
#    per-accession table used for Table S16
python scripts/06_plot_dosage.py --vcf region_bonf.vcf \
    --hapgroup easyhap_bonf/Chr4G_3213203_3231114.HapGroup.tsv \
    --out dosage_vs_height_bonf.png
```

Or run the main analysis (Bonferroni set) end to end:

```bash
bash run_all.sh
```

### 3.2 Figure 4D and supplementary tables (no server needed)

The archived EasyHap tables in `data/` are enough to rebuild the published
figure panel and tables:

```bash
# Fig. 4D: gene model + P/A haplotype heatmap (add --dosage for A / P/A / P)
python scripts/07_gene_haplotype_figure.py --dosage

# Table S16: per-accession class, 4-dp mean dosage and BLUP (docx + xlsx)
python scripts/08_table_s16.py --format both
```

### 3.3 Outputs

| File | Content |
|---|---|
| `region_kmers.tsv` | MAPQ > 20 k-mers in the interval with coordinates |
| `region_dosage.bimbam` | dosage matrix (112 samples) of those k-mers |
| `sig_region.tsv` | nominally significant k-mers in the interval (β, P) |
| `bonf_region.tsv` | subset passing Bonferroni correction |
| `region_*.vcf` | pseudo-VCF (GT = dosage coding) |
| `easyhap_*/` | all EasyHap result tables and figures |
| `dosage_vs_height_bonf.png` | two-panel dosage–phenotype figure (letters, OLS fit and 95% CI band) |
| `dosage_by_accession.tsv` | per-accession class, 4-dp mean dosage and phenotype (source of Table S16) |
| `figures/Fig4D_gene_haplotype_PA.png/.pdf` | Fig. 4D, presence/absence (A/P) representation |
| `figures/Fig4D_gene_haplotype_DOSAGE.png` | the same heatmap with A / P/A / P dosage classes |
| `figures/Fig4F_G_dosage_vs_height.png/.pdf` | Fig. 4F–4G, boxes with letters and dosage regression |
| `supplementary/Table_S16.docx/.xlsx` | Table S16 in three-line layout |

---

## 4. Summary of results

Main analysis: **7 Bonferroni-significant k-mers**, which define **11 genotype
classes** (9 of them contain at least two accessions), with the Bonferroni
threshold α = 0.05 / 34,786,919 = 1.44 × 10⁻⁹.

| Analysis | Groups | Overall KW test | Pairwise result |
|---|---|---|---|
| 9 classes with n ≥ 2 (110 accessions; Table S15) | 9 | H = 18.16, df = 8, **P = 0.0200** | Hap001 vs Hap003: raw P = 1.53×10⁻³, BH q = 0.0552; none of the 36 pairs passed |
| Fig. 4F: Hap001–Hap003 + pooled "Other" (112 accessions; Table S16) | 4 | H = 10.81, df = 3, **P = 0.0128** | Hap001 vs Hap003: raw P = 1.53×10⁻³, **BH q = 0.0092**; only this pair passed |
| Region-significant set (17 markers) | 9 | **P = 0.013** | Hap001 vs Hap003: P = 9.1×10⁻⁴, BH q = 0.033 |
| Region set, 0–4 coding (17 markers) | 9 | **P = 0.038** | — |

Letter display for Fig. 4F: Hap001 = a, Hap002 = ab, Hap003 = b, Other = ab.

Plant height by haplotype class (Bonferroni set):

| Class | n | Mean height | Relative effect |
|---|---|---|---|
| Hap001 | 59 | 149.0 | −3.7% |
| Hap002 | 20 | 156.6 | +1.3% |
| Hap003 | 12 | 173.5 | **+12.2%** |
| Hap009 | 2 | 115.0 | −25.6% |

Continuous mean dosage (arithmetic mean of the 7 k-mers, stored to four
decimal places before analysis) versus plant-height BLUP:
Spearman ρ = −0.330, P = 3.8 × 10⁻⁴ (n = 112);
OLS fit Y = 179.6 − 37.3X, R² = 0.149, slope P = 2.6 × 10⁻⁵,
95% CI of the slope −54.2 to −20.5.

> The Spearman statistic must be computed on the 7-marker mean dosage
> (`--vcf region_bonf.vcf`). A value of ρ = −0.316 refers to the single lead
> k-mer and is **not** interchangeable with the 7-marker result.

**Conclusion.** k-mer copy number at this locus is negatively correlated with
plant height. The tall haplotype (Hap003) is rare in the population (10.7%),
which leaves room for selection in breeding.

---

## 5. Repository layout

```
.
├── README.md
├── config.sh                          # all paths and parameters
├── run_all.sh                         # end-to-end main analysis
├── requirements.txt
├── data/
│   ├── dosage_by_accession_bonf.tsv   # 112 accessions, Table S16 source
│   ├── sample_genotype_tokens_bonf.tsv# EasyHap genotype tokens (7 markers)
│   ├── hapgroup_bonf.tsv              # EasyHap class assignment
│   ├── haplotype_state_bonf_pa.tsv    # A / P state matrix behind Fig. 4D
│   ├── gene_features_04Gg010100.tsv   # UTR / CDS segments of the gene model
│   ├── table_s15_pairwise_9class.tsv  # Table S15 (9-class KW + pairwise tests)
│   ├── table_s15_pairwise_sig17.tsv   # same tests for the 17-marker set
│   └── hap_summary_sig17.tsv          # 17-marker haplotype summary (supplementary)
├── figures/                           # Fig. 4D, 4F, 4G as PNG and PDF
├── supplementary/                     # Table S16 (docx, xlsx)
└── scripts/
    ├── 00_prepare_traits.py           # phenotype → traits.tsv, sample-order check
    ├── 01_extract_region_kmers.py     # SQL query of the interval
    ├── 02_extract_bimbam_region.py    # dosage rows from BIMBAM shards
    ├── 03_select_significant_kmers.py # significant k-mers + Bonferroni
    ├── 04_build_vcf.py                # dosage → pseudo-VCF
    ├── 05_run_easyhap.sh              # EasyHap analysis and figures
    ├── 06_plot_dosage.py              # Fig. 4F-4G + KW/MWU/BH letters + regression
    ├── 07_gene_haplotype_figure.py    # Fig. 4D (A/P heatmap + gene model)
    └── 08_table_s16.py                # Table S16 (docx / xlsx)
```

---

## 6. Methodological notes

1. **k-mer mapping.** Only unique alignments with MAPQ > 20 are retained to avoid
   paralogue/repeat interference.
2. **Ploidy.** *N. porphyrocoma* is an allotetraploid with disomic inheritance
   within subgenomes. The 31-mers used here take only 0/1/2 dosage states in the
   population, so the main analysis uses diploid-like coding (`--ploidy 2`).
   A 0–4 coding (`--ploidy 4`) is provided as a sensitivity check; conclusions
   are unchanged.
3. **Phase-free analysis.** k-mer dosages carry no phase information; EasyHap is
   run in `--mode genotype`, so the classes are multilocus genotype combinations
   rather than chromosome-phased haplotypes.
4. **Marker redundancy.** All significant k-mers in the interval start within a
   24-bp window and represent redundant observations of a single variant; they
   should not be described as independent loci or as a long haplotype block.
5. **Mean dosage and rounding.** For each accession, the BIMBAM dosages of the
   selected markers are averaged ignoring missing values, and the mean is
   rounded to four decimal places before any correlation, regression or
   supplementary table is produced (`DECIMALS = 4` in `06_plot_dosage.py`).
   Rounding at this precision does not change the reported statistics
   (ρ = −0.3298 with or without rounding).
6. **Pairwise tests.** Mann-Whitney U tests are always two-sided and share one
   Benjamini-Hochberg family within the displayed panel (six comparisons for
   Fig. 4F, 36 for the nine-class analysis). The `alternative` argument is
   passed explicitly because its default differs between SciPy versions;
   omitting it on SciPy 1.2.x yields one-sided P values and would halve the
   reported P for Hap001 vs Hap003 (0.00077 instead of 0.00153).

---

## 7. Figure 4 legend and supplementary tables

Provenance of the panels and tables:

| Item | Produced by | Stored as |
|---|---|---|
| Fig. 4D | `scripts/07_gene_haplotype_figure.py` | `figures/Fig4D_gene_haplotype_PA.png/.pdf` |
| Fig. 4F, 4G | `scripts/06_plot_dosage.py` | `figures/Fig4F_G_dosage_vs_height.png/.pdf` |
| Table S15 | EasyHap `TraitAssociationTests.tsv` (step 5) | `data/table_s15_pairwise_9class.tsv` |
| Table S16 | `scripts/08_table_s16.py` | `supplementary/Table_S16.docx/.xlsx` |

### 7.1 Panel D legend (presence/absence wording)

(D) Eleven unphased genotype classes defined by seven Bonferroni-significant
k-mers with alignment starts at Chr4G:3,219,932–3,219,955. Rows show classes and
their sample sizes (right). Columns show markers in genomic order. Green
indicates A (k-mer absent, GT 0/0) and blue indicates P (k-mer present, GT 0/1
or 1/1). Because k-mer markers are presence/absence-type, the REF and ALT
alleles in the pseudo-VCF are placeholder alleles and do not represent actual
nucleotide substitutions; the two states therefore reflect k-mer absence and
presence rather than nucleotide identity. In the gene model, grey denotes
untranslated regions, black denotes coding sequence and horizontal lines denote
introns.

中文对应表述：绿色示 A（k-mer 缺失，GT 0/0），蓝色示 P（k-mer 存在，GT 0/1 或
1/1）；伪 VCF 中的 REF/ALT 仅为占位等位基因，两种状态表示 k-mer 缺失与存在，
而非碱基种类。

### 7.2 Full Fig. 4 legend

**Fig. 4** k-mer associations with six quantitative traits, functional
enrichment of candidate genes and multimarker genotyping of a plant-height
candidate region. (A) Manhattan plot for the six traits. The dashed line
indicates the Bonferroni threshold (raw P = 0.05/34,786,919). (B, C) KEGG
pathways (B) and GO terms (C) for candidate genes associated with plant height,
stalk diameter and internode length. Lines link genes to pathways or terms.
Bubble size indicates the number of genes assigned to a pathway or term, and
colour indicates −log₁₀ of the raw P value. (D) Eleven unphased genotype classes
defined by seven Bonferroni-significant k-mers with alignment starts at
Chr4G:3,219,932–3,219,955. Rows show classes and their sample sizes (right).
Columns show markers in genomic order. Green indicates A (k-mer absent, GT 0/0)
and blue indicates P (k-mer present, GT 0/1 or 1/1). Because k-mer markers are
presence/absence-type, the REF and ALT alleles in the pseudo-VCF are placeholder
alleles and do not represent actual nucleotide substitutions; the two states
therefore reflect k-mer absence and presence rather than nucleotide identity. In
the gene model, grey denotes untranslated regions, black denotes coding sequence
and horizontal lines denote introns. (E) Frequencies of the 11 genotype classes
among the 112 accessions. (F) Plant height BLUP of the EasyHap haplotype classes.
Boxes show the median and interquartile range, whiskers extend to the most
extreme values within 1.5×IQR, and grey dots represent individual accessions
(n given below each class; Hap004–Hap011, 21 accessions in total, were pooled as
"Other"). Overall differences were assessed by Kruskal–Wallis test (H = 10.81,
df = 3, P = 0.0128). Different lowercase letters denote significant differences
based on two-sided Mann–Whitney U tests with Benjamini–Hochberg correction
(P < 0.05); classes sharing a letter do not differ significantly. The 1-copy
(GT 0/1) and 2-copy (GT 1/1) classes are both shown as P in (D); their
quantitative difference is resolved by the continuous dosage analysis in (G).
(G) Relationship between the mean relative dosage of the seven
Bonferroni-significant k-mers (X) and plant height BLUP (Y) across 112
accessions. The red line is the ordinary least-squares fit
(Y = 179.6 − 37.3X, R² = 0.149, P = 2.6×10⁻⁵; 95% CI of the slope, −54.2 to
−20.5) and the shaded band is the 95% confidence band for the fitted line;
Spearman's ρ = −0.330, P = 3.8×10⁻⁴. Both analyses used accession means stored
to four decimal places (Table S16).

---

## 8. References

- Yu et al. (2022) *New Phytologist* 233:1953–1965 — allotetraploid nature of *N. porphyrocoma*.
- He & Shou (2026) *Crop Design* — EasyHap.

## 9. License

MIT License (adjust to your project requirements).
