"""
04_figures.py — main pilot figure (4 panels) from the epithelial checkpoint + DNB table.
Writes results/figures/pilot_v3_EAOC_transition.png.
"""
import os
import scanpy as sc, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

CKPT, FIG, TAB = "checkpoints", "results/figures", "results/tables"; os.makedirs(FIG, exist_ok=True)
epi = sc.read_h5ad(f"{CKPT}/comb_epi_final.h5ad")
ci = pd.read_csv(f"{TAB}/ci_expanded.csv")
nullmax = np.load("data/interim/dnb_nullmax.npy"); thr = np.percentile(nullmax, 95)
plt.rcParams.update({"figure.dpi":120,"font.size":10,"axes.spines.top":False,"axes.spines.right":False})

fig, ax = plt.subplots(2, 2, figsize=(13, 10.5))
um = epi.obsm["X_umap"]; scol = {"0_normal":"#2c7fb8","1_benign_lesion":"#41ae76","2_malignant":"#e34a33"}
a = ax[0,0]
for s,c in scol.items():
    m = (epi.obs.stage==s).values
    a.scatter(um[m,0],um[m,1],s=2,c=c,label=f"{s} ({m.sum()})",linewidths=0,rasterized=True)
a.set_title("A. Combined epithelium (8 datasets)",loc="left",fontweight="bold")
a.set_xlabel("UMAP1"); a.set_ylabel("UMAP2"); a.legend(markerscale=5,fontsize=8,frameon=False)

a = ax[0,1]
groups = [("normal",epi.obs.stage=="0_normal","#2c7fb8"),("benign",epi.obs.stage=="1_benign_lesion","#41ae76")]
for ds,c in [("GSE224334","#e34a33"),("GSE224333","#b30000"),("GSE291389","#7f0000")]:
    groups.append((ds, (epi.obs.dataset==ds)&(epi.obs.stage=="2_malignant"), c))
bp = a.boxplot([epi.obs.cnv_burden[m].values for _,m,_ in groups],
               tick_labels=[g[0] for g in groups], showfliers=False, patch_artist=True, widths=0.6)
for patch,(_,_,c) in zip(bp["boxes"],groups): patch.set_facecolor(to_rgba(c,0.6))
a.set_ylabel("inferCNV burden"); a.set_title("B. Malignancy across 3 EAOC cohorts",loc="left",fontweight="bold")
a.tick_params(axis="x",labelsize=8)

a = ax[1,0]
a.plot(ci.cnv, ci.frac_mal, "-o", color="#d95f0e", ms=4)
tz = ci[(ci.frac_mal>0.1)&(ci.frac_mal<0.9)]
a.axvspan(tz.cnv.min(),tz.cnv.max(),color="orange",alpha=0.15)
a.set_xlabel("CNV burden (malignancy axis)"); a.set_ylabel("fraction malignant")
a.set_title("C. Smooth benign->malignant gradient",loc="left",fontweight="bold")

a = ax[1,1]
a.plot(ci.cnv, ci.CI, "-o", color="#762a83", ms=4, label="DNB critical index")
a.axhline(thr, ls="--", color="grey", label=f"perm 95% null={thr:.4f}")
pk = ci.loc[ci.CI.idxmax()]
a.scatter([pk.cnv],[pk.CI],s=140,facecolors="none",edgecolors="red",linewidths=2,zorder=5,
          label=f"critical peak (frac_mal={pk.frac_mal:.2f})")
a.axvspan(tz.cnv.min(),tz.cnv.max(),color="orange",alpha=0.15,label="transition zone")
a.set_xlabel("CNV burden (malignancy axis)"); a.set_ylabel("critical index (SD x |PCC|)")
a.set_title("D. Critical-transition signal in transition zone",loc="left",fontweight="bold")
a.legend(fontsize=7.5,frameon=False)
fig.suptitle("EAOC epithelial malignant transition | 8 datasets, 3 malignant cohorts",
             fontsize=11, fontweight="bold", y=1.01)
fig.tight_layout(); fig.savefig(f"{FIG}/pilot_v3_EAOC_transition.png", bbox_inches="tight", dpi=150)
print("wrote", f"{FIG}/pilot_v3_EAOC_transition.png")
