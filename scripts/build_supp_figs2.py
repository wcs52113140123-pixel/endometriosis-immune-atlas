"""Build Supplementary Figures S8-S10 at journal print width (180 mm = 7.09 in).

Each shows data that exists in the analysis but appears in no other figure:

  S8   the cross-disease calibration result broken down by tissue and by disease,
       which Supplementary Figure S3 deliberately summarised in three panels
  S9   the full eight-lineage composition of the atlas by tissue, where the main
       atlas figure shows only the immune compartment in detail
  S10  the design and technical behaviour of the segmented spatial cohort:
       how many regions each patient contributed, how much signal each carries,
       and how well the two replicate regions of a patient agree

Plot types follow the conventions already established for this figure set.
"""

TIER = (11.0, 10.0, 9.0)
C_POOL, C_PAT = "#C0392B", "#2E6DA4"
MG = "0.42"
LINC = {"Stromal_Fib": "#7BA7C7", "Epithelial": "#C0392B", "T_NK": "#5B8C5A",
        "Myeloid": "#B8860B", "Endothelial": "#7D6B9E", "SmoothMuscle": "#4E9A9A",
        "B_Plasma": "#C77BA7", "Mast": "#8C6B5B"}
LINLAB = {"Stromal_Fib": "Stromal / fibroblast", "Epithelial": "Epithelial", "T_NK": "T / NK",
          "Myeloid": "Myeloid", "Endothelial": "Endothelial", "SmoothMuscle": "Smooth muscle",
          "B_Plasma": "B / plasma", "Mast": "Mast"}
TISLAB = {"Endometrium": "Endometrium",
          "MEtissue": "Menstrual effluent", "Normal endometrium": "Normal endometrium",
          "OCCC tumor": "Clear-cell carcinoma", "Ovary": "Ovary",
          "Patient-Derived Organoid": "Organoid", "Peritoneum": "Peritoneal lesion",
          "endometriosis": "Endometriotic lesion",
          "endometriosis associated ovarian cancer": "Cancer arising in endometriosis",
          "wholeME": "Whole menstrual effluent", "wholeME/Metissue": "Menstrual effluent, mixed"}


def _letters(fig, ax, keys, dx=0.075, dy=0.030):
    fig.canvas.draw()
    for k in keys:
        p = ax[k].get_position()
        fig.text(p.x0 - dx, p.y1 + dy, k, fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()


def real_overlaps(fig, mpl):
    """Only count two upright texts overlapping by more than two pixels each way."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible() and not t.get_rotation()]
    return [(a.get_text()[:24], b.get_text()[:24]) for i, (a, ba) in enumerate(tx)
            for b, bb in tx[i + 1:]
            if min(ba.x1, bb.x1) - max(ba.x0, bb.x0) > 2 and min(ba.y1, bb.y1) - max(ba.y0, bb.y0) > 2]


def fig_benchmark_breakdown(plt, np, pd, bt, bd, min_tests=150):
    """S8: the calibration failure is not carried by one tissue or one disease."""
    bt = bt[bt.tests >= min_tests].sort_values("fpr_pooled")
    bd = bd[bd.tests >= 100].sort_values("fpr_pooled")
    fig = plt.figure(figsize=(6.35, 7.4))
    gs = fig.add_gridspec(3, 1, hspace=0.55, left=0.30, right=0.955, top=0.955, bottom=0.075)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1]); axC = fig.add_subplot(gs[2])
    AX = {"A": axA, "B": axB, "C": axC}

    y = np.arange(len(bt))
    axA.hlines(y, 0, bt.fpr_pooled.values, color="0.85", lw=1.1, zorder=1)
    axA.scatter(bt.fpr_pooled.values, y, s=30, c=C_POOL, zorder=3, ec="w", lw=.3)
    axA.axvline(0.05, color="0.35", lw=1.0, ls="--", zorder=2)
    axA.text(0.05, len(bt) - 0.2, " nominal", fontsize=TIER[2], color=MG, va="top")
    axA.set_yticks(y); axA.set_yticklabels(bt.tissue.values)
    axA.set_xlim(0, 1.0); axA.set_ylim(-0.7, len(bt) - 0.3)
    axA.set_xlabel("False-positive rate of the pooled-cell test")
    axA.set_title("Every tissue between 0.77 and 0.88", loc="center")
    for yy, n in zip(y, bt.tests.values):
        axA.text(1.0, yy, f"  {n}", fontsize=TIER[2], color=MG, va="center")

    axB.hlines(y, 0, bt.pct_not_upheld.values, color="0.85", lw=1.1, zorder=1)
    axB.scatter(bt.pct_not_upheld.values, y, s=30, c=C_PAT, zorder=3, ec="w", lw=.3)
    axB.set_yticks(y); axB.set_yticklabels(bt.tissue.values)
    axB.set_xlim(0, 100); axB.set_ylim(-0.7, len(bt) - 0.3)
    axB.set_xlabel("Pooled-significant results not upheld at donor level (%)")
    axB.set_title("Most pooled hits fail in seven of the eight tissues", loc="center")

    yd = np.arange(len(bd))
    axC.hlines(yd, 0, bd.fpr_pooled.values, color="0.85", lw=1.1, zorder=1)
    axC.scatter(bd.fpr_pooled.values, yd, s=30, c=C_POOL, zorder=3, ec="w", lw=.3)
    axC.axvline(0.05, color="0.35", lw=1.0, ls="--", zorder=2)
    axC.set_yticks(yd); axC.set_yticklabels(bd.disease.values)
    axC.set_xlim(0, 1.0); axC.set_ylim(-0.7, len(bd) - 0.3)
    axC.set_xlabel("False-positive rate of the pooled-cell test")
    axC.set_title("Same disease by disease, dementia excepted", loc="center")
    for yy, n in zip(yd, bd.tests.values):
        axC.text(1.0, yy, f"  {n}", fontsize=TIER[2], color=MG, va="center")
    for a in (axA, axB, axC):
        for t in a.get_yticklabels(): t.set_fontsize(TIER[2])
    return fig, AX


def fig_lineage_composition(plt, np, pd, obs):
    """S9: all eight lineages by tissue, and how each cohort contributes.

    'Endometirum' is a spelling of 'Endometrium' in one cohort's GEO metadata,
    not a second tissue; the two are merged rather than shown as separate rows,
    which is what an earlier version did and no reader could interpret.
    """
    obs = obs.copy()
    obs["tissue"] = obs.tissue.astype(str).replace({"Endometirum": "Endometrium"})
    keep = [t for t in obs.tissue.unique() if (obs.tissue == t).sum() >= 500]
    comp = (pd.crosstab(obs.loc[obs.tissue.isin(keep), "tissue"],
                        obs.loc[obs.tissue.isin(keep), "celltype"], normalize="index") * 100)
    order = comp.index[np.argsort(-comp.get("Epithelial", pd.Series(0, index=comp.index)).values)]
    comp = comp.loc[order]
    ncell = obs.tissue.value_counts()
    nsamp = obs.groupby("tissue", observed=True).gsm.nunique()

    fig = plt.figure(figsize=(6.55, 6.3))
    gs = fig.add_gridspec(2, 1, hspace=0.60, height_ratios=[1.25, 1.0],
                          left=0.30, right=0.88, top=0.945, bottom=0.145)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])
    AX = {"A": axA, "B": axB}

    y = np.arange(len(comp)); left = np.zeros(len(comp))
    for l in LINC:
        if l not in comp.columns: continue
        axA.barh(y, comp[l].values, left=left, color=LINC[l], height=.72, label=LINLAB[l])
        left += comp[l].values
    axA.set_yticks(y); axA.set_yticklabels([TISLAB.get(t, t) for t in comp.index])
    axA.set_xlim(0, 100); axA.set_ylim(-0.7, len(comp) - 0.3)
    axA.set_xlabel("Share of cells in the tissue (%)")
    axA.set_title("Lineage composition differs by tissue", loc="center")
    for yy, t in zip(y, comp.index):
        axA.text(101.5, yy, f"{ncell[t] / 1000:.0f}k, {nsamp[t]}", fontsize=TIER[2], color=MG, va="center")
    fig.legend(handles=[plt.Line2D([], [], marker="s", ls="", ms=6, mfc=LINC[l], mec="none",
                                   label=LINLAB[l]) for l in LINC],
               frameon=False, fontsize=TIER[2], loc="lower center", bbox_to_anchor=(0.5, -0.005),
               ncol=4, handletextpad=.35, columnspacing=1.4)

    ct = pd.crosstab(obs.dataset, obs.celltype, normalize="index") * 100
    y2 = np.arange(len(ct)); left = np.zeros(len(ct))
    for l in LINC:
        if l not in ct.columns: continue
        axB.barh(y2, ct[l].values, left=left, color=LINC[l], height=.66)
        left += ct[l].values
    axB.set_yticks(y2); axB.set_yticklabels([d.replace("GSE", "GSE ") for d in ct.index])
    axB.set_xlim(0, 100); axB.set_ylim(-0.7, len(ct) - 0.3)
    axB.set_xlabel("Share of cells in the cohort (%)")
    axB.set_title("and by cohort, which is why cohort is a covariate", loc="center")
    for a in (axA, axB):
        for t in a.get_yticklabels(): t.set_fontsize(TIER[2])
    return fig, AX


def fig_spatial_design(plt, np, pd, mn, dm, sig):
    """S10: design and technical behaviour of the segmented spatial cohort."""
    SEG = ["Macrophages", "Stroma", "Epithelium"]
    SEGC = {"Macrophages": "#B8860B", "Stroma": "#7BA7C7", "Epithelium": "#C0392B"}
    TISC = {"eutopic endometrium": "#7BA7C7", "endometriotic lesion": "#C0392B"}
    fig = plt.figure(figsize=(7.09, 6.6))
    gs = fig.add_gridspec(3, 1, hspace=0.62, left=0.20, right=0.955, top=0.945, bottom=0.075)
    axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1]); axC = fig.add_subplot(gs[2])
    AX = {"A": axA, "B": axB, "C": axC}

    pats = sorted(dm.sampleid.unique())
    cnt = np.zeros((len(SEG) * 2, len(pats)), int)
    rowlab = []
    for i, seg in enumerate(SEG):
        for j, tis in enumerate(["eutopic endometrium", "endometriotic lesion"]):
            rowlab.append(f"{seg}, {'eutopic' if j == 0 else 'lesion'}")
            for k, p in enumerate(pats):
                cnt[i * 2 + j, k] = ((dm["cell type"] == seg) & (dm.tissue == tis)
                                     & (dm.sampleid == p)).sum()
    axA.imshow(np.ones_like(cnt), cmap="Blues", vmin=0, vmax=6, aspect="auto")
    for a in range(cnt.shape[0]):
        for b in range(cnt.shape[1]):
            axA.text(b, a, cnt[a, b], ha="center", va="center", fontsize=TIER[2], color="0.15")
    axA.set_xticks(range(len(pats))); axA.set_xticklabels(pats)
    axA.set_yticks(range(len(rowlab))); axA.set_yticklabels(rowlab)
    axA.set_xlabel("Patient"); axA.set_title("Two regions per patient in every segment and tissue", loc="center")
    for t in axA.get_yticklabels(): t.set_fontsize(TIER[2])

    det = (mn > mn.values.mean()).sum(1)
    rng = np.random.default_rng(3)
    for i, seg in enumerate(SEG):
        for j, tis in enumerate(["eutopic endometrium", "endometriotic lesion"]):
            m = (dm["cell type"] == seg) & (dm.tissue == tis)
            v = det[m.values].values; x = i + (j - 0.5) * 0.32
            axB.scatter(x + rng.uniform(-0.055, 0.055, len(v)), v, s=24, c=TISC[tis],
                        alpha=.9, ec="w", lw=.3, zorder=3)
            axB.plot([x - 0.1, x + 0.1], [np.median(v)] * 2, color="k", lw=1.3, zorder=4)
    axB.set_xticks(range(len(SEG))); axB.set_xticklabels(SEG)
    axB.set_ylabel("Genes above the array mean")
    axB.set_xlim(-0.6, len(SEG) - 0.4)
    axB.set_title("Signal depth is comparable across segments and tissues", loc="center")
    axB.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="none",
                                   label="Eutopic endometrium" if k.startswith("eutopic") else "Lesion")
                        for k, c in TISC.items()],
               frameon=False, fontsize=TIER[2], loc="lower right", handletextpad=.3)

    z = ((mn[sig] - mn[sig].mean()) / mn[sig].std().replace(0, 1)).mean(1)
    pairs = []
    for seg in SEG:
        for tis in ["eutopic endometrium", "endometriotic lesion"]:
            for p in pats:
                m = (dm["cell type"] == seg) & (dm.tissue == tis) & (dm.sampleid == p)
                v = z[m.values].values
                if len(v) == 2: pairs.append((seg, v[0], v[1]))
    PD_ = pd.DataFrame(pairs, columns=["segment", "r1", "r2"])
    for seg in SEG:
        s = PD_[PD_.segment == seg]
        axC.scatter(s.r1, s.r2, s=34, c=SEGC[seg], label=seg, ec="w", lw=.35, zorder=3)
    lim = [PD_[["r1", "r2"]].values.min() - .25, PD_[["r1", "r2"]].values.max() + .25]
    axC.plot(lim, lim, color="0.6", lw=.9, ls="--", zorder=1)
    rho = float(np.corrcoef(PD_.r1, PD_.r2)[0, 1])
    axC.set_xlim(lim); axC.set_ylim(lim)
    axC.set_xlabel("First region of the pair"); axC.set_ylabel("Second region")
    axC.set_title(f"Replicate regions of a patient agree (r = {rho:.2f}, {len(PD_)} pairs)", loc="center")
    axC.legend(frameon=False, fontsize=TIER[2], loc="upper left", handletextpad=.3)
    return fig, AX
