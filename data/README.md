# Dataset manifest

Every dataset used is public and was deposited by other groups. Nothing in `data/raw/` is
redistributed in this repository; `scripts/00_download_data.py` fetches it from the original
repositories (~3.3 GB). Each accession is also cited in full in the manuscript reference list under
the `[dataset]` tag, as the journal's data citation policy requires.

## Single-cell and single-nucleus RNA-seq

| Accession | Content | Role in this study |
|---|---|---|
| GSE179640 | Endometriosis, multi-site (eutopic endometrium, peritoneum, ovary) with controls | Atlas; per-sample case–control; ligand–receptor inference |
| GSE203191 | Menstrual endometrial tissue | Atlas |
| GSE213216 | Endometriosis atlas including ovarian endometrioma | Per-sample Treg comparison (third cohort) |
| GSE214411 | Minimal/mild endometriosis endometrium with controls | Atlas; per-sample case–control |
| GSE247695 | Eutopic endometrium and ectopic lesion **paired within patient** | Independent paired validation of the macrophage programme |
| GSE224333 | Ovarian clear cell carcinoma, mono- and co-culture | Malignant epithelium (EAOC extension) |
| GSE224334 | Ovarian clear cell carcinoma, single-nucleus | Malignant epithelium (EAOC extension) |
| GSE291389 | Endometriosis-associated ovarian cancer | Malignant epithelium (EAOC extension) |

## Cell-type-segmented spatial profiling (GeoMx DSP)

| Accession | Content | Role |
|---|---|---|
| GSE263897 | Superficial peritoneal endometriotic lesion vs eutopic endometrium, macrophage / stromal / epithelial segments, 5 patients paired | **Primary spatial finding**; patient-level paired test |
| GSE303150 | Adenomyosis lesion vs eutopic endometrium, CD45 / panCK / stroma / CD31 segments, 10 patients | Independent test of an adjacent contrast — **did not support** the finding |

Raw GeoMx data are per-AOI `.dcc` count files plus a `.pkc` probe definition; the gene-level matrix
is assembled by summing probes to their target gene (see the spatial section of the Methods).

## Bulk transcriptomics

| Accession | Content | Role |
|---|---|---|
| GSE130435 | Sorted M1/M2 endometrial macrophages, endometriosis vs control | Independent test of an adjacent contrast — **did not support** the finding |
| GSE65986 | Ovarian clear cell and endometrioid carcinoma with progression-free survival | Histology-matched survival test (null, 15 events in 55 patients) |
| TCGA-OV | Ovarian carcinoma, predominantly high-grade serous | Survival test in a histology-mismatched cohort (null) |

## Cross-disease comparison

| Source | Content | Role |
|---|---|---|
| CZ CELLxGENE Census, release 2025-11-08 | Harmonised cell metadata (`donor_id`, `disease`, `cell_type`, `tissue`, `assay`) for 159 million human cells | 523 case–control settings, 124 diseases, 6111 cell-type comparisons (Supplementary Fig. S1) |

Only **cell metadata** is required — no expression matrices are downloaded for this part, which is
why the whole cross-disease permutation runs in seconds. Access requires the `cellxgene-census`
package and network access to the Census S3 bucket.

## Derived files kept in the repository

`data/interim/` holds small derived tables (per-sample statistics, benchmark outputs, design tables,
gene-coordinate maps). Large binary intermediates (`*.pkl`, `*.npy`, `*.h5ad`) are excluded from git
and rebuilt by the scripts; the manifest of what each script writes is in `scripts/README.md`.
