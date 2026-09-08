"""Shared helpers for the EAOC-scRNA pipeline."""
import numpy as np, pandas as pd, scipy.sparse as sp

# ---- lineage marker sets (per-cluster mean-score annotation) ----
MARKERS = {
    "Epithelial":  ["EPCAM","KRT8","KRT18","KRT19","CDH1","PAX8","KRT17","MUC1"],
    "Stromal_Fib": ["DCN","COL1A1","COL1A2","PDGFRA","LUM","THY1"],
    "Endothelial": ["PECAM1","VWF","CLDN5","CD34"],
    "T_NK":        ["CD3D","CD3E","CD8A","IL7R","NKG7","GNLY"],
    "Myeloid":     ["CD68","LYZ","CD14","AIF1","C1QA"],
    "B_Plasma":    ["MS4A1","CD79A","MZB1","IGHG1"],
    "SmoothMuscle":["ACTA2","MYH11","TAGLN","DES"],
}

def run_harmony(adata, batch_key="gsm", rep="X_pca", out="X_pca_harmony", max_iter=20):
    """Harmony via harmonypy directly. The scanpy wrapper is incompatible with the
    installed harmonypy (mis-assigns Z_corr shape); this handles either orientation."""
    import harmonypy
    ho = harmonypy.run_harmony(adata.obsm[rep], adata.obs, [batch_key], max_iter_harmony=max_iter)
    Zc = np.asarray(ho.Z_corr)
    adata.obsm[out] = Zc if Zc.shape[0] == adata.n_obs else Zc.T
    return adata

def annotate_lineage(adata, cluster_key="leiden", use_raw=True):
    """Assign each cluster to the lineage with the highest mean marker score."""
    import scanpy as sc
    ref = adata.raw.var_names if use_raw else adata.var_names
    for k, v in MARKERS.items():
        sc.tl.score_genes(adata, [g for g in v if g in ref], score_name=f"sc_{k}", use_raw=use_raw)
    cols = [f"sc_{k}" for k in MARKERS]
    assign = adata.obs.groupby(cluster_key)[cols].mean().idxmax(axis=1).str[3:]
    adata.obs["celltype"] = adata.obs[cluster_key].map(assign).astype("category")
    return adata

def gene_positions_from_gtf(gtf_gz):
    """Parse an Ensembl GTF (.gtf.gz) -> DataFrame indexed by gene symbol with
    chromosome/start/end (autosomes + X/Y). Used to place genes for inferCNV."""
    import gzip, re
    pat = re.compile(r'gene_name "([^"]+)"'); rows = []
    with gzip.open(gtf_gz, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            c = line.split("\t")
            if len(c) > 8 and c[2] == "gene":
                m = pat.search(c[8])
                if m: rows.append((m.group(1), c[0], int(c[3]), int(c[4])))
    ann = pd.DataFrame(rows, columns=["gene","chromosome","start","end"])
    chroms = [str(i) for i in range(1,23)] + ["X","Y"]
    return ann[ann.chromosome.isin(chroms)].drop_duplicates("gene").set_index("gene")

def cnv_burden(adata, key="X_cnv"):
    """Per-cell CNV burden = mean(|inferCNV signal|)."""
    X = adata.obsm[key]; X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return np.mean(np.abs(X), axis=1)

def dnb_index(expr_df):
    """DNB-style critical index for cells (rows) x module genes (cols):
    CI = mean_gene_SD * mean_abs_pairwise_correlation. Peaks at a critical transition."""
    sd = expr_df.std(0).mean()
    C = np.corrcoef(expr_df.values.T)
    iu = np.triu_indices_from(C, k=1)
    return sd * np.nanmean(np.abs(C[iu]))

# ---- tissue -> stage mapping (normal / benign lesion / malignant) ----
STAGE_MAP = {
    # normal / eutopic / ovary
    "Control":"0_normal","Eutopic":"0_normal","Eutopic Endometrium":"0_normal","Eutopic endometrium":"0_normal",
    "No endometriosis detected":"0_normal","Ovary_normal":"0_normal","Ovary":"0_normal",
    "Normal endometrium":"0_normal","Endometrium":"0_normal","Endometirum":"0_normal",
    "MEtissue":"0_normal","wholeME":"0_normal","wholeME/Metissue":"0_normal",
    # benign lesion
    "Endometriosis":"1_benign_lesion","Ectopic_lesion":"1_benign_lesion","Ectopic endometrium":"1_benign_lesion",
    "Endometrioma":"1_benign_lesion","Peritoneum":"1_benign_lesion","endometriosis":"1_benign_lesion",
    "Patient-Derived Organoid":"1_benign_lesion",
    # malignant
    "OCCC tumor":"2_malignant","endometriosis associated ovarian cancer":"2_malignant",
}
