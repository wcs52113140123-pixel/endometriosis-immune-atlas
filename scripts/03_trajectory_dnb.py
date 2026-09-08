"""
03_trajectory_dnb.py — malignancy gradient (CNV burden) as progression axis;
DNB critical-transition index along it with a 1000-permutation null.
Writes results/tables/ci_expanded.csv.

Note: benign and malignant epithelium form two blobs (large CNV/expression gap and
no dedicated atypical/precursor samples), so a global diffusion pseudotime is
dominated by within-compartment heterogeneity. We therefore order cells along the
externally-grounded CNV-burden axis (malignancy replicated across 3 EAOC cohorts).
"""
import os
import scanpy as sc, numpy as np, pandas as pd, scipy.sparse as sp
from utils import dnb_index

CKPT, TAB = "checkpoints", "results/tables"; os.makedirs(TAB, exist_ok=True)
NBIN, NPERM, MODULE = 20, 1000, 300

epi = sc.read_h5ad(f"{CKPT}/comb_epi_final.h5ad")
# balanced set: all benign/normal + subsample malignant so the axis is not tumor-dominated
rng = np.random.default_rng(0)
mal = epi.obs_names[epi.obs.stage == "2_malignant"].to_numpy()
mal_sub = rng.choice(mal, min(20000, len(mal)), replace=False)
ben = epi.obs_names[epi.obs.stage != "2_malignant"].to_numpy()
b = epi[np.concatenate([ben, mal_sub])].copy()
b.X = b.layers["counts"].copy(); sc.pp.normalize_total(b, 1e4); sc.pp.log1p(b)
sc.pp.highly_variable_genes(b, n_top_genes=2000)
hvg = b.var_names[b.var.highly_variable][:MODULE]
X = b[:, hvg].X; X = X.toarray() if sp.issparse(X) else np.asarray(X)
X = pd.DataFrame(X, index=b.obs_names, columns=hvg)

b.obs["cbin"] = pd.qcut(b.obs.cnv_burden, NBIN, labels=False, duplicates="drop")
rows = []
for bb in sorted(b.obs.cbin.dropna().unique()):
    idx = b.obs_names[b.obs.cbin == bb]
    rows.append((int(bb), len(idx), b.obs.loc[idx, "cnv_burden"].mean(),
                 (b.obs.loc[idx, "stage"] == "2_malignant").mean(), dnb_index(X.loc[idx])))
ci = pd.DataFrame(rows, columns=["bin","n","cnv","frac_mal","CI"])

# permutation null: shuffle bin labels, record max CI
rng = np.random.default_rng(1); nullmax = []
for _ in range(NPERM):
    perm = b.obs.cbin.sample(frac=1, random_state=int(rng.integers(1e9))).values
    cis = [dnb_index(X.loc[b.obs_names[perm == bb]]) for bb in sorted(pd.unique(perm[~pd.isna(perm)]))
           if (perm == bb).sum() >= 50]
    nullmax.append(max(cis))
nullmax = np.array(nullmax); thr = np.percentile(nullmax, 95)
pk = ci.loc[ci.CI.idxmax()]
p = max((nullmax >= pk.CI).mean(), 1/NPERM)
print(ci.round(4).to_string(index=False))
print(f"\nnull95(n={NPERM})={thr:.4f} | peak bin{int(pk.bin)} CI={pk.CI:.4f} "
      f"frac_mal={pk.frac_mal:.2f} p<={p:.3f}")
ci.to_csv(f"{TAB}/ci_expanded.csv", index=False)
np.save("data/interim/dnb_nullmax.npy", nullmax)
print("wrote", f"{TAB}/ci_expanded.csv")
