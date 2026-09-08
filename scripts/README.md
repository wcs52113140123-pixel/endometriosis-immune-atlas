# Scripts

Run from the repository root with the `endo` environment active.

## Pipeline

| Script | Does | Writes |
|---|---|---|
| `00_download_data.py` | Fetches every public matrix listed in `data/README.md` | `data/raw/` (~3.3 GB, git-ignored) |
| `01_load_qc_integrate.py` | QC, normalisation, HVG, PCA, Harmony integration by sample, Leiden clustering, lineage and immune-state annotation | `checkpoints/full_atlas.h5ad`, `checkpoints/immune_annotated.h5ad` |
| `02_epithelial_infercnv.py` | Epithelial subset, gene coordinates from the Ensembl GTF, inferred copy number, per-cell burden | `checkpoints/comb_epi_final.h5ad` |
| `03_trajectory_dnb.py` | Malignancy axis along copy-number burden, transition-state index, and the two-population mixing null that refuted it | `data/interim/ci_expanded.csv`, `p4_obs_vs_mixing.csv`, `p4_mixnull.csv` |
| `utils.py` | Harmony wrapper (fixes the `Z_corr` orientation), marker sets, shared helpers | — |

## Figures

One script per figure. `*_print.py` renders each figure at full-page width with the
11/10/9 pt type ladder; the non-print variants are the earlier wide-canvas versions kept for the
record.

| Script | Figure | Content |
|---|---|---|
| `build_fig1_print.py` | Figure 1 | Atlas: lineage and cohort embeddings, composition by site and disease status, scIB integration benchmark |
| `build_fig2_print.py` | Figure 2 | Unit of replication: per-sample Treg fractions, pooled vs patient-mean estimates, cell-yield relationship, effect-size forest, permutation false-positive rate |
| `build_fig3_print.py` | Figure 3 | Spatial compartment specificity in 8 panels, including the two cohorts that do not support the finding |
| `build_fig4_print.py` | Figure 4 | Cell–cell communication around Treg: circle plot, chord diagram, ligand–receptor bubble, information flow, signalling role, expression dot plot |
| `build_fig5_print.py` | Figure 5 | EAOC extension: gain/loss frequency track, per-chromosome effect, burden on the embedding, per-sample violins, transition index against the mixing null |
| `build_fig6_print.py` | Supplementary Fig. S1 | Cross-disease calibration (the script name is unchanged for provenance) |

`unit_of_replication.py` is the reusable diagnostic; `UNIT_OF_REPLICATION.md` at the repository
root documents it.

`figure_layout.py` holds the shared layout pass — panel letter and title on one line, titles level
within a row, row spacing measured from rendered ink — and its `LAYOUT` table records the arguments
used for each figure.
