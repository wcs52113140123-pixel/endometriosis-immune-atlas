# The unit of replication in single-cell studies of the endometriotic immune microenvironment

[Chinese version](README.zh-CN.md)

Analysis code and derived results for a patient-level re-analysis of eight public single-cell,
spatial and bulk datasets of endometriosis, together with a cross-disease calibration study of the
same statistical question.

Authors: X. Wang, K. Liang, D. Yang, B. Xiong, J. Liao, Z. Wang —
Department of Obstetrics and Gynecology, Affiliated Hospital of Zunyi Medical University, Zunyi,
Guizhou, China.

---

## The question

Case-control single-cell studies of endometriosis almost always pool cells from all donors in a
group and compare the pooled cell-level distributions. Cells from one donor are not independent
observations. This work asks what changes when the **patient**, rather than the **cell**, is the
unit of replication.

Three results:

1. **Treg enrichment reproduces when cells are pooled and disappears when patients are the unit.**
   Per cohort at patient level, *P* = 0.20, 0.63 and 0.46; across all 81 samples with cohort as a
   stratum, β = 0.064 on the log scale (95% CI −0.56 to 0.69, *P* = 0.84). Per-sample M2,
   CD8-exhaustion and NK-cytotoxicity scores are likewise null.
2. **One finding survives a patient-level test.** A nine-gene M2/immunosuppression programme is
   higher in lesion macrophage segments than in the same woman's eutopic endometrium in all five
   paired patients (paired *P* = 0.031), and only in macrophage segments (stroma *P* = 0.41,
   epithelium *P* = 0.97). Two further cohorts testing adjacent contrasts do not support it, so the
   claim is restricted to macrophage segments of peritoneal lesions.
3. **The calibration failure is general.** Across 523 public case-control settings in 124 diseases
   (6,111 cell-type comparisons), the pooled-cell test rejects a true null at a median rate of 0.84
   against 0.04 for the donor-level test, and 69% of pooled-significant results are not upheld at
   the donor level.

---

## Layout

```
.
├── scripts/                     analysis and figure code
├── results/
│   ├── figures/                 rendered figures
│   └── tables/                  every table behind the numbers above
├── data/
│   ├── README.md                dataset manifest: accessions, role, retrieval
│   ├── interim/                 derived intermediates
│   └── bulk/                    bulk expression and survival matrices
├── UNIT_OF_REPLICATION.md       documentation for the reusable diagnostic
├── environment.yml              conda environment (python 3.11)
└── requirements.txt             pip equivalent
```

## Running it

```bash
conda env create -f environment.yml && conda activate endo
python scripts/00_download_data.py           # retrieve the public matrices into data/raw
python scripts/01_load_qc_integrate.py       # QC, Harmony integration, lineage annotation
python scripts/02_epithelial_infercnv.py     # epithelial subset, inferred copy number
python scripts/05_statistics.py --part all   # patient-level tests, power, method comparison
python scripts/build_fig2_print.py           # figures, one script per figure
```

`scripts/run_all.py --list` prints the scripted steps. `scripts/figure_layout.py` holds the shared
figure layout pass and records the arguments used for each figure.

## The diagnostic

`scripts/unit_of_replication.py` applies the same test to any dataset. Give it a donor x cell-type
table and a group label; it runs the pooled-cell and donor-level tests, calibrates both by permuting
the group label across donors, and reports which pooled-significant results survive:

```python
from unit_of_replication import diagnose, risk_note
res = diagnose(df, donor="donor_id", group="group", cell_type="cell_type", case="case")
print(risk_note(res))
```

`UNIT_OF_REPLICATION.md` documents the arguments and the exact Mann-Whitney null used for small
arms.

## Data

All primary data are public. `data/README.md` lists every accession with its role in the analysis
and the retrieval command; `scripts/00_download_data.py` fetches them from the original
repositories.

## Licence and citation

Code is released under the MIT licence (`LICENSE`). The third-party datasets remain under the terms
set by their original depositors. If you use this code, please cite `CITATION.cff`.

## Fonts

Figures are typeset in Arial at an 11/10/9 pt ladder. Arial is licensed separately and is not
included here; place `ARIAL.TTF`, `ARIALBD.TTF`, `ARIALI.TTF` and `ARIALBI.TTF` in `assets/fonts/`
to reproduce the typography exactly.
