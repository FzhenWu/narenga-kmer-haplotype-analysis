#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 8: Table S16 -- per-accession haplotype class, 4-dp mean k-mer dosage and
plant-height BLUP, laid out as a three-line table.

Input : data/dosage_by_accession_bonf.tsv (written by 06_plot_dosage.py or
        provided as archived data)
Output: <prefix>.docx and/or <prefix>.xlsx (English headers and notes)

Usage:
    python scripts/08_table_s16.py --format both
"""
import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

MAIN = ["Hap001", "Hap002", "Hap003"]

TITLE = ("Table S16  Haplotype class, mean relative dosage of the seven "
         "Bonferroni-significant k-mers, and plant height BLUP for the 112 "
         "Narenga porphyrocoma accessions")
HEADERS = ["Sample ID", "Haplotype class", "Mean relative dosage",
           "Plant height BLUP"]
NOTES = [
    ("Note: ", True),
    ("N = 112 accessions. Hap001-Hap003 are the three major haplotype classes; "
     "the remaining eight classes (Hap004-Hap011, 21 accessions) were pooled as "
     "\"Other\".", False),
    ("1. Mean relative dosage: for each accession, the relative dosages of the "
     "seven Bonferroni-significant k-mers were read from the BIMBAM matrix, "
     "missing values were ignored, and the arithmetic mean was stored to four "
     "decimal places before any statistic was calculated.", False),
    ("2. Haplotype class: dosages were rounded to 0/1/2 copies to build the GT "
     "field (0/0, 0/1, 1/1) and classes were obtained with EasyHap 1.4.0 in "
     "genotype mode. k-mer data carry no phase information, so classes are "
     "phase-agnostic.", False),
    ("3. Plant height BLUP: breeding values from a mixed linear model; used as "
     "the phenotype in all association tests and figures.", False),
    ("4. Kruskal-Wallis (four displayed groups): H = 10.81, df = 3, P = 0.0128. "
     "Two-sided Mann-Whitney U tests with Benjamini-Hochberg correction across "
     "six comparisons: only Hap001 vs Hap003 was significant (raw P = 0.00153; "
     "BH-adjusted P = 0.0092). Letter display: Hap001 = a, Hap002 = ab, "
     "Hap003 = b, Other = ab.", False),
    ("5. Environment: Python with SciPy (scipy.stats.kruskal, "
     "scipy.stats.mannwhitneyu).", False),
]


def read_rows(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    data = []
    for r in rows:
        data.append((r[0], r[1] if r[1] in MAIN else "Other",
                     float(r[2]), float(r[3])))
    return data


def write_docx(data, out):
    from docx import Document
    from docx.shared import Pt, Cm
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml.ns import qn
    from docx.oxml import OxmlElement

    def cell_font(cell, text, bold=False):
        cell.text = ""
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(text)
        r.bold = bold
        r.font.size = Pt(9)
        r.font.name = "Times New Roman"

    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.2)
    sec.top_margin = sec.bottom_margin = Cm(2.0)
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(10.5)

    p = doc.add_paragraph()
    r = p.add_run(TITLE)
    r.bold = True
    r.font.size = Pt(10.5)
    r.font.name = "Times New Roman"

    widths = [Cm(4.0), Cm(3.6), Cm(4.0), Cm(4.0)]
    table = doc.add_table(rows=1, cols=4)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    for i, h in enumerate(HEADERS):
        c = table.rows[0].cells[i]
        c.width = widths[i]
        cell_font(c, h, bold=True)
    for s, h, d, ph in data:
        cells = table.add_row().cells
        vals = [s, h, "%.4f" % d, "%.2f" % ph]
        for i, v in enumerate(vals):
            cells[i].width = widths[i]
            cell_font(cells[i], v)

    tblPr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "bottom"):
        el = OxmlElement("w:%s" % edge)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "12")
        el.set(qn("w:color"), "000000")
        borders.append(el)
    for edge in ("left", "right", "insideH", "insideV"):
        el = OxmlElement("w:%s" % edge)
        el.set(qn("w:val"), "none")
        borders.append(el)
    tblPr.append(borders)
    for cell in table.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        b = OxmlElement("w:tcBorders")
        el = OxmlElement("w:bottom")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "6")
        el.set(qn("w:color"), "000000")
        b.append(el)
        tcPr.append(b)

    doc.add_paragraph()
    for txt, bold in NOTES:
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(txt)
        r.bold = bold
        r.font.size = Pt(8.5)
        r.font.name = "Times New Roman"
    doc.save(out)
    print("saved: " + out)


def write_xlsx(data, out):
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment, Border, Side

    wb = Workbook()
    ws = wb.active
    ws.title = "Table S16"
    thin = Side(style="thin", color="000000")
    thick = Side(style="medium", color="000000")

    ws.merge_cells("A1:D1")
    c = ws["A1"]
    c.value = TITLE
    c.font = Font(name="Times New Roman", size=11, bold=True)
    c.alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    ws.row_dimensions[1].height = 46

    for j, h in enumerate(HEADERS, 1):
        cell = ws.cell(row=2, column=j, value=h)
        cell.font = Font(name="Times New Roman", size=10, bold=True)
        cell.alignment = Alignment(horizontal="center", vertical="center",
                                   wrap_text=True)
        cell.border = Border(top=thick, bottom=thin)
    ws.row_dimensions[2].height = 30

    r0 = 3
    for i, (s, h, d, ph) in enumerate(data):
        r = r0 + i
        ws.cell(row=r, column=1, value=s)
        ws.cell(row=r, column=2, value=h)
        ws.cell(row=r, column=3, value=d).number_format = "0.0000"
        ws.cell(row=r, column=4, value=round(ph, 2)).number_format = "0.00"
        for j in range(1, 5):
            cell = ws.cell(row=r, column=j)
            cell.font = Font(name="Times New Roman", size=10)
            cell.alignment = Alignment(horizontal="center", vertical="center")

    last = r0 + len(data) - 1
    for j in range(1, 5):
        ws.cell(row=last, column=j).border = Border(bottom=thick)

    nr = last + 2
    for txt, bold in NOTES:
        ws.merge_cells(start_row=nr, start_column=1, end_row=nr, end_column=4)
        cell = ws.cell(row=nr, column=1, value=txt)
        cell.font = Font(name="Times New Roman", size=9, bold=bold)
        cell.alignment = Alignment(horizontal="left", vertical="top",
                                   wrap_text=True)
        ws.row_dimensions[nr].height = 14 * (1 + len(txt) // 95)
        nr += 1

    for col, w in zip("ABCD", (16, 16, 21, 19)):
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A3"
    ws.sheet_view.showGridLines = False
    wb.save(out)
    print("saved: " + out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=os.path.join(ROOT, "data", "dosage_by_accession_bonf.tsv"))
    ap.add_argument("--out-prefix", default=os.path.join(ROOT, "supplementary", "Table_S16"))
    ap.add_argument("--format", choices=["docx", "xlsx", "both"], default="both")
    args = ap.parse_args()

    data = read_rows(args.input)
    if len(data) != 112:
        print("[warn] expected 112 accessions, found %d" % len(data))
    if args.format in ("docx", "both"):
        write_docx(data, args.out_prefix + ".docx")
    if args.format in ("xlsx", "both"):
        write_xlsx(data, args.out_prefix + ".xlsx")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
