"""
02_epithelial_infercnv.py — subset epithelium, run inferCNV (normal+benign as
reference), compute per-cell CNV burden. Writes checkpoints/comb_epi_final.h5ad.
"""
import os, urllib.request
import scanpy as sc, numpy as np, pandas as pd, infercnvpy as cnv
from utils import gene_positions_from_gtf, cnv_burden, STAGE_MAP

CKPT = "checkpoints"
GTF = "data/interim/genes.gtf.gz"
GTF_URL = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

adata = sc.read_h5ad(f"{CKPT}/all_integrated.h5ad")
epi = adata[adata.obs.celltype == "Epithelial"].copy()
epi.obs["stage"] = epi.obs.tissue.map(STAGE_MAP).astype("category")
epi = epi[epi.obs.stage.isin(["0_normal","1_benign_lesion","2_malignant"])].copy()

# ---- gene genomic positions (Ensembl GTF) ----
if not os.path.exists(GTF):
    os.makedirs("data/interim", exist_ok=True); urllib.request.urlretrieve(GTF_URL, GTF)
ann = gene_positions_from_gtf(GTF).reindex(epi.var_names)
epi.var["chromosome"] = ("chr" + ann.chromosome.astype(str)).values
epi.var["start"] = ann.start.values; epi.var["end"] = ann.end.values
epi = epi[:, epi.var.chromosome.notna() & epi.var.start.notna()].copy()
epi.var["start"] = epi.var.start.astype(int); epi.var["end"] = epi.var.end.astype(int)

# ---- inferCNV: reference = normal + benign epithelium ----
epi.obs["cnv_ref"] = np.where(epi.obs.stage == "2_malignant", "tumor", "normal")
cnv.tl.infercnv(epi, reference_key="cnv_ref", reference_cat=["normal"], window_size=100)
epi.obs["cnv_burden"] = cnv_burden(epi)
print(epi.obs.groupby("stage")["cnv_burden"].agg(["count","median"]).round(5))
epi.write(f"{CKPT}/comb_epi_final.h5ad")
print("wrote", f"{CKPT}/comb_epi_final.h5ad", epi.shape)
