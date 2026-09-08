"""
01_load_qc_integrate.py — load 8 GEO datasets, QC, normalize, integrate (Harmony),
cluster, annotate lineages. Writes checkpoints/all_integrated.h5ad.

Assumes data/raw/<GSE>/ holds either per-GSM folders with 10x mtx triplets, or
GSM*.h5 (cellranger filtered). Tissue labels per GSM from data/interim/*tissue*.json.
"""
import os, glob, json, warnings
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from utils import run_harmony, annotate_lineage
warnings.simplefilter("ignore"); sc.settings.verbosity = 1

RAW, CKPT = "data/raw", "checkpoints"; os.makedirs(CKPT, exist_ok=True)
# per-GSM tissue map assembled during the study (GEO characteristics 'tissue' tag)
TISSUE = {}
for j in glob.glob("data/interim/*tissue*.json") + glob.glob("data/interim/*supp*.json"):
    try:
        d = json.load(open(j))
        for acc, v in d.items():
            if isinstance(v, dict):
                for gsm, t in v.items():
                    if isinstance(t, str): TISSUE[gsm] = t
    except Exception: pass

def load_dataset(acc):
    out = []
    for d in sorted(glob.glob(f"{RAW}/{acc}/GSM*")):
        if os.path.isdir(d):
            a = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True)
        elif d.endswith(".h5"):
            a = sc.read_10x_h5(d); a.var_names_make_unique()
        else: continue
        gsm = os.path.basename(d).replace(".h5", "")
        a.obs["dataset"] = acc; a.obs["gsm"] = gsm; a.obs["tissue"] = TISSUE.get(gsm, "NA")
        a.obs_names = [f"{gsm}_{b}" for b in a.obs_names]; out.append(a)
    for h in sorted(glob.glob(f"{RAW}/{acc}/*.h5")):
        gsm = os.path.basename(h)[:-3]
        a = sc.read_10x_h5(h); a.var_names_make_unique()
        a.obs["dataset"] = acc; a.obs["gsm"] = gsm; a.obs["tissue"] = TISSUE.get(gsm, "NA")
        a.obs_names = [f"{gsm}_{b}" for b in a.obs_names]; out.append(a)
    return out

accs = [os.path.basename(p) for p in glob.glob(f"{RAW}/GSE*") if os.path.isdir(p)]
ads = [a for acc in accs for a in load_dataset(acc)]
adata = ad.concat(ads, join="outer", index_unique=None); del ads
print("combined:", adata.shape, "| samples:", adata.obs.gsm.nunique())

# ---- QC ----
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[(adata.obs.n_genes_by_counts >= 300) & (adata.obs.n_genes_by_counts < 8000)
              & (adata.obs.pct_counts_mt < 20)].copy()
sc.pp.filter_genes(adata, min_cells=10)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata); adata.raw = adata

# ---- HVG -> PCA -> Harmony(gsm) -> cluster ----
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="dataset")
aw = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(aw, max_value=10); sc.tl.pca(aw, n_comps=50)
run_harmony(aw, batch_key="gsm")                      # -> aw.obsm['X_pca_harmony']
adata.obsm["X_pca_harmony"] = aw.obsm["X_pca_harmony"]
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata, min_dist=0.3)
annotate_lineage(adata)
print(adata.obs.celltype.value_counts())
adata.write(f"{CKPT}/all_integrated.h5ad")
print("wrote", f"{CKPT}/all_integrated.h5ad")
