"""Build Supplementary Figures S3-S5 at journal print width (180 mm = 7.09 in).

Each figure closes a specific hole a reviewer would probe:

  S3  the ten-state patient-level composition behind Supplementary Table S5 -
      the direct evidence that no immune state, not only Treg, separates cases
      from controls once the patient is the unit
  S4  the sensitivity of the one surviving positive finding - leaving out any
      patient, any gene, and swapping in published marker sets
  S5  the transcriptional cycle-phase estimator: that it tracks the reported
      phase in the samples that have one, and that it does not differ by disease

Same conventions as the main figures: Arial at an 11/10/9 pt ladder, panel
letters bold outside the axes, axes frames levelled within a row, and blue for
control / red for endometriosis threaded across every panel.
"""

C_CTRL, C_ENDO = "#7BA7C7", "#C0392B"
C_POS, C_NEG = "#C0392B", "#9AA5AD"
TIER = (11.0, 10.0, 9.0)
MG = "0.35"


def _letters(fig, ax, keys, dy=0.028, dx=0.070):
    for a in ax.values():
        for t in list(a.texts):
            if t.get_text() in keys and t.get_fontweight() in ("bold", 700):
                t.remove()
    for t in list(fig.texts):
        if t.get_text() in keys and t.get_fontweight() in ("bold", 700):
            t.remove()
    fig.canvas.draw()
    for k in keys:
        p = ax[k].get_position()
        fig.text(p.x0 - dx, p.y1 + dy, k, fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()


def _level_titles(fig, ax, rows):
    """Recompute each title's pad so a row's titles share one baseline.

    title.set_position() does not hold: matplotlib recomputes the title from the
    pad on every draw.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
    fh = fig.get_size_inches()[1]
    for rw in rows:
        bot = {k: T.transform(ax[k].title.get_window_extent(r))[0][1] for k in rw}
        target = max(bot.values())
        for k in rw:
            ax[k].set_title(ax[k].get_title(), loc="center",
                            pad=(target - ax[k].get_position().y1) * fh * 72)
        fig.canvas.draw(); r = fig.canvas.get_renderer()


def fig_composition(plt, np, pd, comp, scan, n_ctrl, n_endo):
    """S3: per-patient composition of ten immune states, and the test on each."""
    STATES = list(scan.immune_state)                      # ordered by combined P
    fig = plt.figure(figsize=(7.09, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.42,
                          left=0.14, right=0.97, top=0.885, bottom=0.10)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    AX = {"A": axA, "B": axB}

    rng = np.random.default_rng(3)
    d = comp[comp.condition.isin(["Control", "Endometriosis"])]
    for i, st in enumerate(STATES):
        for off, (cond, col) in zip((0.18, -0.18), [("Control", C_CTRL), ("Endometriosis", C_ENDO)]):
            v = d.loc[d.condition == cond, st].values
            y = len(STATES) - 1 - i + off
            axA.scatter(v, y + rng.uniform(-0.07, 0.07, len(v)), s=15, c=col,
                        alpha=.85, edgecolors="w", linewidths=.3, zorder=3)
            axA.plot([np.median(v)] * 2, [y - 0.13, y + 0.13], color="k", lw=1.4, zorder=4)
    axA.set_yticks(range(len(STATES))); axA.set_yticklabels(STATES[::-1])
    axA.set_xlabel("Per-patient share of immune cells (%)")
    axA.set_title("Every state overlaps between groups", loc="center")
    axA.set_xscale("symlog", linthresh=1); axA.set_xticks([0, 1, 10, 50])
    axA.set_xticklabels(["0", "1", "10", "50"])
    axA.set_ylim(-0.7, len(STATES) - 0.3)
    # figure-level key: inside panel A it sat on the CD4T points, and above panel A
    # it reached the title
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="none", label=l)
                        for l, c in [(f"Control (n = {n_ctrl} patients)", C_CTRL),
                                     (f"Endometriosis (n = {n_endo} patients)", C_ENDO)]],
               frameon=False, fontsize=TIER[2], loc="upper center", bbox_to_anchor=(0.5, 1.0),
               ncol=2, handletextpad=.3, columnspacing=1.6)

    nlp = -np.log10(scan.P_combined.values)
    y = np.arange(len(STATES))[::-1]
    axB.hlines(y, 0, nlp, color="0.80", lw=1.0, zorder=1)
    axB.scatter(nlp, y, s=34, c=C_NEG, edgecolors="w", linewidths=.4, zorder=3)
    axB.axvline(-np.log10(0.05), color=C_POS, lw=1.0, ls="--", zorder=2)
    axB.text(-np.log10(0.05), len(STATES) - 0.45, " P = 0.05", color=C_POS,
             fontsize=TIER[2], va="center", ha="left")
    for yy, q in zip(y, scan.q_combined_BH.values):
        axB.text(nlp[len(STATES) - 1 - yy] + 0.045, yy, f"q = {q:.2f}", va="center",
                 fontsize=TIER[2], color=MG)
    axB.set_yticks(range(len(STATES))); axB.set_yticklabels([])
    axB.set_xlabel("Combined evidence, $-$log$_{10}$ P")
    axB.set_title("None survives correction", loc="center")
    axB.set_xlim(0, max(nlp) + 0.62); axB.set_ylim(-0.7, len(STATES) - 0.3)
    return fig, AX


def fig_sensitivity(plt, np, pd, lop, log, gs_):
    """S4: how much the spatial finding depends on one patient, one gene, one marker set."""
    fig = plt.figure(figsize=(7.09, 7.4))
    g = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.35, 1.0], hspace=0.92,
                         left=0.30, right=0.90, top=0.94, bottom=0.075)
    axA = fig.add_subplot(g[0]); axB = fig.add_subplot(g[1]); axC = fig.add_subplot(g[2])
    AX = {"A": axA, "B": axB, "C": axC}

    def rows(ax, labels, diffs, ps, ref, title, xlab, ups=None, italic=False):
        y = np.arange(len(labels))[::-1]
        ax.axvline(0, color="0.75", lw=0.8, zorder=1)
        ax.axvline(ref, color="0.55", lw=0.9, ls=":", zorder=1)
        ax.hlines(y, 0, diffs, color="0.85", lw=1.0, zorder=2)
        ax.scatter(diffs, y, s=32, c=[C_POS if p < 0.05 else C_NEG for p in ps],
                   edgecolors="w", linewidths=.4, zorder=3)
        for i_, (yy, dd, pp) in enumerate(zip(y, diffs, ps)):
            txt = f"P = {pp:.3f}".rstrip("0").rstrip(".")
            if ups is not None: txt = f"{ups[i_]} up   " + txt
            ax.text(dd + 0.045 * max(diffs), yy, txt, va="center", fontsize=TIER[2], color=MG)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels[::-1], style="italic" if italic else "normal")
        ax.set_xlabel(xlab); ax.set_title(title, loc="center")
        ax.set_xlim(min(0, min(diffs)) - 0.05, max(diffs) * 1.34)
        ax.set_ylim(-0.7, len(labels) - 0.3)

    rows(axA, [f"without patient {p}" for p in lop.patient_left_out],
         lop.mean_diff.values, lop.P.values, 0.785,
         "Every four-patient subset still moves one way", "Lesion $-$ eutopic score",
         ups=list(lop.up))
    axA.text(0.5, -0.42, "4 of 4 patients up in every subset; with four pairs the smallest attainable "
             "one-sided P is 0.062,\nso these are the strongest values the design allows. "
             "Dotted line: all five patients (0.785).",
             transform=axA.transAxes, ha="center", va="top", fontsize=TIER[2], color=MG)
    rows(axB, list(log.gene_removed),
         log.mean_diff.values, log.P.values, 0.785,
         "Dropping any one gene keeps it too", "Lesion $-$ eutopic score", italic=True)

    lab = ["Nine-gene composite\n(this study)", "WikiPathways\nmacrophage markers",
           "Ad hoc disjoint\nalternative activation", "Published M2 signature\n(188 genes)"]
    order = [0, 2, 3, 1]
    dd = gs_.macrophage_diff.values[order]; pp = gs_.macrophage_P.values[order]
    ng = gs_.n_genes_measured.values[order]; ov = gs_.overlap_with_composite.values[order]
    y = np.arange(len(lab))[::-1]
    axC.axvline(0, color="0.75", lw=0.8, zorder=1)
    axC.hlines(y, 0, dd, color="0.85", lw=1.0, zorder=2)
    axC.scatter(dd, y, s=34, c=[C_POS if p < 0.05 else C_NEG for p in pp],
                edgecolors="w", linewidths=.4, zorder=3)
    for yy, d_, p_, n_, o_ in zip(y, dd, pp, ng, ov):
        axC.text(max(dd) * 0.06 + max(d_, 0.02), yy,
                 f"P = {p_:.3f}".rstrip("0").rstrip(".") + f"   {n_} genes, {o_} shared",
                 va="center", fontsize=TIER[2], color=MG)
    axC.set_yticks(range(len(lab))); axC.set_yticklabels(lab[::-1])
    axC.set_xlabel("Lesion $-$ eutopic score")
    axC.set_title("But it depends on which markers define the programme", loc="center")
    axC.set_xlim(min(dd) - 0.10, max(dd) * 2.05); axC.set_ylim(-0.7, len(lab) - 0.3)
    axC.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="none", label=l)
                        for l, c in [("P < 0.05", C_POS), ("not significant", C_NEG)]],
               frameon=False, fontsize=TIER[2], loc="lower right", bbox_to_anchor=(1.0, 0.10),
               handletextpad=.3)
    return fig, AX


def fig_phase(plt, np, pd, ph, auc, p_cond):
    """S5: the cycle-phase estimator - validated where a phase is reported, flat by disease."""
    fig = plt.figure(figsize=(7.09, 3.5))
    g = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.34,
                         left=0.10, right=0.98, top=0.86, bottom=0.17)
    axA = fig.add_subplot(g[0]); axB = fig.add_subplot(g[1])
    AX = {"A": axA, "B": axB}
    rng = np.random.default_rng(5); nlab = {}

    lab = ph.dropna(subset=["reported_phase"])
    for i, (phase, col) in enumerate([("proliferative", "#8E9AAF"), ("secretory", "#4B6584")]):
        v = lab.loc[lab.reported_phase == phase, "phase_score"].values
        axA.scatter(i + rng.uniform(-0.10, 0.10, len(v)), v, s=26, c=col,
                    edgecolors="w", linewidths=.4, zorder=3)
        axA.plot([i - 0.20, i + 0.20], [np.median(v)] * 2, color="k", lw=1.5, zorder=4)
        nlab[i] = len(v)
    axA.set_xticks([0, 1])
    axA.set_xticklabels([f"Proliferative\nn = {nlab[0]}", f"Secretory\nn = {nlab[1]}"])
    axA.set_ylabel("Cycle-phase score")
    axA.set_title(f"Tracks the reported phase (AUC {auc:.2f})", loc="center")
    axA.set_xlim(-0.6, 1.6); axA.set_yticks([-2, -1, 0, 1, 2])

    cohorts = list(dict.fromkeys(ph.dataset))
    for i, ds in enumerate(cohorts):
        for off, (cond, col) in zip((-0.16, 0.16), [("Control", C_CTRL), ("Endometriosis", C_ENDO)]):
            v = ph.loc[(ph.dataset == ds) & (ph.condition == cond), "phase_score"].values
            if not len(v): continue
            axB.scatter(i + off + rng.uniform(-0.055, 0.055, len(v)), v, s=22, c=col,
                        alpha=.9, edgecolors="w", linewidths=.3, zorder=3)
            axB.plot([i + off - 0.11, i + off + 0.11], [np.median(v)] * 2, color="k", lw=1.3, zorder=4)
    axB.set_xticks(range(len(cohorts)))
    axB.set_xticklabels([c.replace("GSE", "GSE ") for c in cohorts])
    axB.set_ylabel("Cycle-phase score")
    axB.set_title(f"Does not differ by disease (P = {p_cond:.2f})", loc="center")
    axB.set_xlim(-0.6, len(cohorts) - 0.4); axB.set_yticks([-2, -1, 0, 1, 2])
    axB.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="none", label=l)
                        for l, c in [("Control", C_CTRL), ("Endometriosis", C_ENDO)]],
               frameon=False, fontsize=TIER[2], loc="lower left", handletextpad=.3)
    return fig, AX


def fig_markers(plt, np, pd, dot, cnt, MK):
    """S6: the marker evidence behind the ten immune states."""
    states = list(cnt.index); genes = list(dict.fromkeys(dot.gene))
    piv_m = dot.pivot(index="state", columns="gene", values="mean_expr").loc[states, genes]
    piv_f = dot.pivot(index="state", columns="gene", values="frac_pos").loc[states, genes]
    z = piv_m.apply(lambda c: (c - c.mean()) / (c.std() if c.std() else 1), axis=0)

    fig = plt.figure(figsize=(6.55, 4.3))
    gs = fig.add_gridspec(1, 1, left=0.155, right=0.845, top=0.90, bottom=0.30)
    ax = fig.add_subplot(gs[0, 0])
    vmax = float(np.nanmax(np.abs(z.values)))
    X, Y = np.meshgrid(np.arange(len(genes)), np.arange(len(states)))
    sc_ = ax.scatter(X.ravel(), Y.ravel(), s=piv_f.values.ravel() * 62,
                     c=z.values.ravel(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                     edgecolors="0.55", linewidths=.25, zorder=3)
    # block boundaries: each gene group marks one state
    gstate = {g: st for st in MK for g in MK[st]}
    bounds = [i for i in range(1, len(genes)) if gstate[genes[i]] != gstate[genes[i - 1]]]
    for b in bounds:
        ax.axvline(b - 0.5, color="0.88", lw=0.6, zorder=1)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, ha="center", style="italic")
    ax.set_yticks(range(len(states)))
    ax.set_yticklabels([f"{st}  ({cnt[st]:,})" for st in states])
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(-0.6, len(states) - 0.4)
    ax.invert_yaxis()
    ax.set_title("Each state expresses its own canonical markers", loc="center", pad=10)
    cax = fig.add_axes([0.872, 0.53, 0.020, 0.30])
    cb = fig.colorbar(sc_, cax=cax, ticks=[-2, 0, 2])
    cb.ax.tick_params(labelsize=TIER[2], length=2, pad=1)
    cb.set_label("Mean expression\n(z across states)", fontsize=TIER[2], labelpad=6)
    for frac, lab in [(0.25, "25%"), (0.75, "75%")]:
        ax.scatter([], [], s=frac * 62, c="0.75", ec="0.5", lw=.25, label=lab)
    # size key stacked under the colour bar: above the panel it reached the title
    ax.legend(frameon=False, fontsize=TIER[2], loc="upper left", bbox_to_anchor=(1.035, 0.42),
              handletextpad=.35, labelspacing=.7, title="cells\nexpressing",
              title_fontsize=TIER[2])
    return fig, {"A": ax}


def fig_qc(plt, np, pd, obs):
    """S7: what quality control kept, per cohort and per sample."""
    ds = list(dict.fromkeys(obs.dataset))
    fig = plt.figure(figsize=(7.09, 5.6))
    gs = fig.add_gridspec(2, 1, hspace=0.72, left=0.16, right=0.97, top=0.94, bottom=0.11)
    top = gs[0].subgridspec(1, 2, wspace=0.34)
    axA = fig.add_subplot(top[0]); axB = fig.add_subplot(top[1]); axC = fig.add_subplot(gs[1])
    AX = {"A": axA, "B": axB, "C": axC}

    def box(ax, col, ylab, title, thresh=None):
        data = [obs.loc[obs.dataset == d, col].values for d in ds]
        bp = ax.boxplot(data, positions=range(len(ds)), widths=.55, showfliers=False,
                        patch_artist=True, medianprops=dict(color="k", lw=1.2))
        for p in bp["boxes"]:
            p.set_facecolor("#C9D6E0"); p.set_edgecolor("0.45"); p.set_linewidth(.7)
        if thresh is not None:
            ax.axhline(thresh, color=C_POS, lw=1.0, ls="--", zorder=1)
        ax.set_xticks(range(len(ds)))
        ax.set_xticklabels([d.replace("GSE", "") for d in ds], rotation=90)
        ax.set_ylabel(ylab); ax.set_title(title, loc="center")
        ax.set_xlim(-0.7, len(ds) - 0.3)

    box(axA, "n_genes_by_counts", "Genes per cell", "Detected genes", thresh=None)
    for t_ in (300, 8000):
        axA.axhline(t_, color=C_POS, lw=1.0, ls="--", zorder=1)
    axA.text(len(ds) - 0.35, 8000, " cut-offs", color=C_POS, fontsize=TIER[2], va="center", ha="left")
    axA.set_xlabel("Cohort (GSE)")
    box(axB, "pct_counts_mt", "Mitochondrial reads (%)", "Mitochondrial fraction", thresh=20)
    axB.text(len(ds) - 0.35, 20, " cut-off", color=C_POS, fontsize=TIER[2], va="center", ha="left")
    axB.set_xlabel("Cohort (GSE)")

    per = obs.groupby("gsm", observed=True).agg(n=("dataset", "size"),
                                                dataset=("dataset", "first")).sort_values("n")
    col = {d: c for d, c in zip(ds, ["#7BA7C7", "#C0392B", "#5B8C5A", "#B8860B", "#7D6B9E"])}
    axC.bar(range(len(per)), per.n.values, color=[col[d] for d in per.dataset], width=.82)
    axC.set_yscale("log"); axC.set_ylabel("Cells retained")
    axC.set_xlabel(f"Sample, ordered by yield (n = {len(per)})")
    axC.set_xticks([]); axC.set_xlim(-1, len(per))
    axC.set_title("Cells passing quality control in each sample", loc="center")
    axC.legend(handles=[plt.Line2D([], [], marker="s", ls="", ms=6, mfc=col[d], mec="none",
                                   label=d.replace("GSE", "GSE ")) for d in ds],
               frameon=False, fontsize=TIER[2], loc="upper left", ncol=2, handletextpad=.3,
               columnspacing=1.2)
    return fig, AX
