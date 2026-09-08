"""Build Figure 2 (unit of replication, 5 panels).

Plot types follow a caption census of 77 spatial/single-cell immune papers
(629 figures): box/violin 52%, scatter 39%, bar 35%, forest 12%, while paired
connected-line plots appear in 1%.  The previous version of this figure used a
connected-line panel for the spatial result; that result is now shown in
Figure 3, where it belongs, and this figure is only about the unit of
replication.

  A  per-sample Treg fraction by cohort and condition, box plus points
  B  the same difference under the two estimators, grouped bars
  C  Treg fraction against immune cells recovered per sample, with Spearman rho
  D  forest plot of per-cohort, combined and sensitivity effects
  E  permutation null: rejection rate of the pooled-cell test versus the
     patient-level test when condition labels are shuffled

Panel letters are placed at a constant figure-space offset (see place_letters),
so each row of letters sits on one baseline regardless of panel height.
"""

C_CTRL, C_ENDO = "#7BA7C7", "#C0392B"
C_POOL, C_PAT = "#C0392B", "#5B6B7B"
TIER = (11.0, 10.0, 9.0)   # unified ladder across Figures 1-4
DSS = ["GSE179640", "GSE213216", "GSE214411"]


def place_letters(fig, ax, dy=0.026, dx=0.052):
    for a in fig.axes:
        for t in list(a.texts):
            if t.get_text() in list("ABCDE") and t.get_fontweight() in ("bold", 700):
                t.remove()
    for t in list(fig.texts):
        if t.get_text() in list("ABCDE") and t.get_fontweight() in ("bold", 700):
            t.remove()
    for lab in "ABCDE":
        p = ax[lab].get_position()
        fig.text(max(p.x0 - dx, 0.004), min(p.y1 + dy, 0.995), lab,
                 fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()
    return {lab: min(ax[lab].get_position().y1 + dy, 0.995) for lab in "ABCDE"}


def build(plt, np, pd, ctx):
    g = ctx
    allps, T, PM, FO, FPR, rho = (g["allps"], g["T"], g["PM"], g["FO"], g["FPR"], g["rho"])
    # wider at the 11/10/9 ladder: at 9/8/7 the forest row labels fitted inside
    # panel D, at the larger sizes they reached into panel E
    # PRINT-SIZE VARIANT: 180 mm journal column = 7.09 in, so the 11/10/9 ladder
    # lands near full size on the page instead of being scaled to ~50%
    # 2 columns x 3 rows.  The 2 x 5 arrangement of the wide version does not
    # survive at 180 mm: each panel would be 1.1 in across, too narrow for the
    # cohort tick labels and the forest row labels.  Panel titles are noun
    # phrases here, as in Figure 3, with the message left to the caption.
    fig = plt.figure(figsize=(7.09, 8.4))
    gs = fig.add_gridspec(3, 2, hspace=0.52, wspace=0.42,
                          left=0.115, right=0.985, top=0.945, bottom=0.055)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, :])   # spans both columns: a bar chart reads better wide

    # ---- A: box + points, per-sample Treg by cohort x condition ----
    xt = []
    for i, dsn in enumerate(DSS):
        gg = allps[allps.dataset == dsn]
        for j, (cd, col) in enumerate([("Control", C_CTRL), ("Endometriosis", C_ENDO)]):
            v = gg.loc[gg.cond == cd, "Treg_pct"].values
            x = i * 2.6 + j * 0.9
            axA.boxplot([v], positions=[x], widths=.66, patch_artist=True, showfliers=False,
                        medianprops=dict(color="0.15", lw=1.3),
                        boxprops=dict(fc=col, ec="0.35", lw=.8, alpha=.55),
                        whiskerprops=dict(color="0.5", lw=.8), capprops=dict(color="0.5", lw=.8))
            rngj = np.random.default_rng(11 + i * 2 + j)
            axA.scatter(x + rngj.uniform(-.17, .17, len(v)), v, s=20, c=col, ec="w", lw=.35, zorder=3)
        p = float(T.loc[T.dataset == dsn, "p_onesided"].iloc[0])
        axA.text(i * 2.6 + .45, 13.2, f"P = {p:.2f}", ha="center", fontsize=TIER[2], color="0.25")
        n0 = int((gg.cond == "Control").sum()); n1 = int((gg.cond == "Endometriosis").sum())
        axA.text(i * 2.6 + .45, -1.75, f"{n0} vs {n1}", ha="center", fontsize=TIER[2], color="0.45")
        xt.append(i * 2.6 + .45)
    axA.set_xticks(xt); axA.set_xticklabels(DSS)
    axA.set_ylabel("Treg (% of immune cells)"); axA.set_ylim(-2.4, 14.2); axA.margins(x=0.05)
    axA.set_title("Per-sample Treg fractions", loc="center", pad=10)
    axA.legend(handles=[plt.Line2D([], [], marker="s", ls="", ms=6, mfc=c, mec="none", label=l)
                        for l, c in [("Control", C_CTRL), ("Endometriosis", C_ENDO)]],
               frameon=False, fontsize=TIER[2], ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, -0.14), handletextpad=.3, columnspacing=1.6)

    # ---- B: grouped bars, the two estimators ----
    w = 0.34
    for i, dsn in enumerate(DSS):
        s_ = PM[PM.dataset == dsn].set_index("cond")
        d_pool = s_.loc["Endometriosis", "pooled"] - s_.loc["Control", "pooled"]
        d_pat = s_.loc["Endometriosis", "patient_mean"] - s_.loc["Control", "patient_mean"]
        axB.bar(i - w / 2, d_pool, w, color=C_POOL, alpha=.85,
                label="pooled cells" if i == 0 else None)
        axB.bar(i + w / 2, d_pat, w, color=C_PAT, alpha=.85,
                label="patient mean" if i == 0 else None)
        axB.text(i, max(d_pool, d_pat) + 0.07, f"Δ {d_pool - d_pat:+.2f}", ha="center",
                 va="bottom", fontsize=TIER[2], color="0.35")
    axB.axhline(0, color="0.6", lw=.8)
    axB.set_xticks(range(3)); axB.set_xticklabels(DSS)
    axB.set_ylabel("Endometriosis − control\n(Treg percentage points)")
    axB.set_ylim(0, 1.95); axB.margins(x=0.14)
    axB.set_title("Pooled vs patient-mean estimate", loc="center", pad=10)
    axB.legend(frameon=False, fontsize=TIER[2], ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, -0.14), handletextpad=.4, columnspacing=1.6)

    # ---- C: scatter with Spearman ----
    for cd, col in [("Control", C_CTRL), ("Endometriosis", C_ENDO)]:
        gg = allps[allps.cond == cd]
        axC.scatter(gg.n_immune, gg.Treg_pct, s=22, c=col, alpha=.85, ec="w", lw=.35,
                    zorder=3, label=cd)
    top2 = allps.nlargest(2, "Treg_pct")
    axC.scatter(top2.n_immune, top2.Treg_pct, s=54, facecolors="none", edgecolors="0.25",
                lw=1.0, zorder=4)
    axC.annotate("2 of 70 case samples,\n14% of all case Treg cells",
                 xy=(float(top2.n_immune.mean()), float(top2.Treg_pct.min())),
                 xytext=(2.1e3, 7.6), fontsize=TIER[2], color="0.25",
                 arrowprops=dict(arrowstyle="-", lw=.6, color="0.6"))
    axC.set_xscale("log"); axC.set_xlabel("Immune cells recovered per sample")
    axC.set_ylabel("Treg (% of immune cells)"); axC.margins(y=0.12)
    # the Spearman statistic moved to the caption: in the panel it sat on the
    # topmost point, and its place is better used by the condition key
    axC.set_title("Treg fraction vs cells recovered", loc="center", pad=10)
    axC.legend(frameon=False, fontsize=TIER[2], ncol=1, loc="upper left",
               bbox_to_anchor=(-0.01, 1.02), handletextpad=.3, labelspacing=.3)

    # ---- D: forest ----
    yy = np.arange(len(FO))[::-1]
    for y, (_, r) in zip(yy, FO.iterrows()):
        comb = r.label.startswith("Combined")
        col = "#1F2933" if comb else "#5B6B7B"
        axD.plot([r.lo, r.hi], [y, y], color=col, lw=2.2 if comb else 1.4,
                 solid_capstyle="butt", zorder=2)
        axD.scatter([r.beta], [y], s=64 if comb else 40, marker="D" if comb else "o",
                    c=col, zorder=3, ec="w", lw=.6)
        axD.text(-1.46, y, r.label.replace("Combined (stratified)", "Combined").replace("Combined + ", "+ "),
                 fontsize=TIER[2], va="center", ha="left", color="0.25", zorder=5,
                 bbox=dict(fc="white", ec="none", pad=1.0))
        axD.text(1.58, y, f"{r.n}  P={r.p:.2f}", fontsize=TIER[2], va="center", ha="right",
                 color="0.35", zorder=5, bbox=dict(fc="white", ec="none", pad=1.0))
    axD.axvline(0, color=C_ENDO, lw=.9, ls="--", zorder=1)
    axD.set_yticks([])
    axD.set_xlabel("Effect on log(Treg % + 0.1),\nendometriosis − control")
    axD.set_xlim(-1.5, 1.62); axD.margins(y=0.16)
    axD.set_title("Patient-level effect estimates", loc="center", pad=10)

    # ---- E: permutation null, rejection rate of the two tests ----
    w = 0.34
    for i, dsn in enumerate(DSS):
        f = FPR[dsn]
        axE.bar(i - w / 2, f["fpr_pooled"] * 100, w, color=C_POOL, alpha=.85,
                label="pooled-cell test" if i == 0 else None)
        axE.bar(i + w / 2, f["fpr_patient"] * 100, w, color=C_PAT, alpha=.85,
                label="patient-level test" if i == 0 else None)
        axE.text(i - w / 2, f["fpr_pooled"] * 100 + 1.6, f"{f['fpr_pooled']*100:.0f}%",
                 ha="center", fontsize=TIER[2], color="0.3")
        axE.text(i + w / 2, f["fpr_patient"] * 100 + 1.6, f"{f['fpr_patient']*100:.0f}%",
                 ha="center", fontsize=TIER[2], color="0.3")
    axE.axhline(5, color="0.35", lw=1.0, ls="--")
    axE.text(1.02, 5 / 96, "nominal 5%", transform=axE.transAxes, fontsize=TIER[2],
             color="0.35", ha="left", va="center", clip_on=False)
    axE.set_xticks(range(3)); axE.set_xticklabels(DSS)
    axE.set_ylabel("False positives at α = 0.05\n(%, 2,000 label permutations)")
    axE.set_ylim(0, 96); axE.margins(x=0.14)
    axE.set_title("False positives under permutation", loc="center", pad=10)
    # explicit in-range ticks: an off-range tick label stays a live Text object
    # that never renders but collides in the overlap audit
    _y0, _y1 = axE.get_ylim()
    axE.set_yticks([v for v in np.arange(0, 101, 20) if _y0 <= v <= _y1])
    axE.legend(frameon=False, fontsize=TIER[2], ncol=1, loc="upper center",
               bbox_to_anchor=(0.5, -0.16), handletextpad=.4)

    ax = dict(A=axA, B=axB, C=axC, D=axD, E=axE)
    g["normalize"](fig)
    return fig, ax


def audit(fig, mpl):
    """Text-overlap / typography audit. Draws first: measuring straight after a
    savefig(bbox_inches="tight") returns stale extents and reports phantom overlaps.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    t = [(x, x.get_window_extent(r)) for x in fig.findobj(mpl.text.Text)
         if x.get_text().strip() and x.get_visible()]
    ov = [(a.get_text()[:24], b.get_text()[:24])
          for i, (a, ba) in enumerate(t) for b, bb in t[i + 1:] if ba.overlaps(bb)]
    bold = [w for ax in fig.axes for lbl, w in
            [(ax.xaxis.label, "x"), (ax.yaxis.label, "y"), (ax.title, "title")]
            if lbl.get_text().strip() and lbl.get_fontweight() not in ("normal", 400)]
    return dict(n_overlap=len(ov), overlaps=ov, bold_axis_titles=bold,
                fonts=sorted({x.get_fontname() for x, _ in t}),
                sizes=sorted({round(x.get_fontsize(), 1) for x, _ in t}))


def align_row_bottoms(fig, axes_list, passes=2):
    """Make a row of panels end on one line.

    Panels in a gridspec row share an axes bottom, but what hangs below differs
    (a legend under one, a two-line axis title under another), so the visual
    bottoms do not match.  Measure each panel's lowest rendered descendant and
    shift the axes so those lowest elements coincide.
    """
    def lowest(fig, ax, r, T):
        ys = [T.transform(ax.get_window_extent(r))[0][1]]
        for t in ax.get_xticklabels():
            if t.get_text().strip():
                ys.append(T.transform(t.get_window_extent(r))[0][1])
        if ax.xaxis.label.get_text().strip():
            ys.append(T.transform(ax.xaxis.label.get_window_extent(r))[0][1])
        lg = ax.get_legend()
        if lg is not None:
            ys.append(T.transform(lg.get_window_extent(r))[0][1])
        return min(ys)

    for _ in range(passes):
        fig.canvas.draw()
        r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
        over = {id(a): a.get_position().y0 - lowest(fig, a, r, T) for a in axes_list}
        target = min(lowest(fig, a, r, T) for a in axes_list)
        top = max(a.get_position().y1 for a in axes_list)
        for a in axes_list:
            p = a.get_position()
            y0 = target + over[id(a)]
            a.set_position([p.x0, y0, p.width, top - y0])
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
    lows = [lowest(fig, a, r, T) for a in axes_list]
    return dict(lows=[round(x, 4) for x in lows], spread=round(max(lows) - min(lows), 5))
