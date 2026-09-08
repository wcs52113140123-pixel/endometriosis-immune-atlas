"""Build Figure 3 (compartment specificity of the spatial finding, 6 panels).

Plot types follow the convention of segmented digital-spatial-profiling papers
(e.g. the ESCC DSP atlas in Genome Medicine 2024, whose PanCK / CD45 / stroma
segmentation matches our epithelium / macrophage / stroma AOIs): an embedding of
all AOIs annotated by segment, a gene x AOI heatmap with segment and tissue
annotation bars, box plots with the paired test, a volcano, and a dot plot.
A caption census of 77 spatial-immune papers (629 figures) found heatmaps in
73%, box/violin in 52%, dot plots in 49% and volcanoes in 47%, while paired
connected-line plots appeared in 1%.  The previous version of this figure was
built entirely from connected-line plots.

A  PCA of all 60 AOIs, coloured by segment, open/filled by tissue
B  nine-gene heatmap across AOIs with segment and tissue annotation bars
C  composite score by segment and tissue, box plus per-patient points
D  volcano of lesion vs eutopic within macrophage AOIs (paired, 5 patients)
E  dot plot of the nine genes by segment: colour = log2FC, size = -log10 P
F  independent paired single-cell cohort, same nine-gene score
"""

C_EUT, C_LES = "#7BA7C7", "#C0392B"
C_UP, C_DOWN = "#C0392B", "#9AA5AD"
TIER = (11.0, 10.0, 9.0)   # unified ladder across Figures 1-4
# one legend baseline per row, so the bottom edge of each row is a single line
LEG_Y1, LEG_Y2 = -0.19, -0.15
SEGS = ["Macrophages", "Stroma", "Epithelium"]
SEGC = {"Macrophages": "#E1A730", "Stroma": "#7F8FA6", "Epithelium": "#59A14F"}
SEGSHORT = {"Macrophages": "Macrophage", "Stroma": "Stroma", "Epithelium": "Epithelium"}
TISC = {"eutopic endometrium": C_EUT, "endometriotic lesion": C_LES}
SEGLAB = {"Macrophages": "Macrophage segments", "Stroma": "Stromal segments",
          "Epithelium": "Epithelial segments"}


def place_letters(fig, ax, dy=0.030, dx=0.055):
    for a in fig.axes:
        for t in list(a.texts):
            if t.get_text() in list("ABCDEF") and t.get_fontweight() in ("bold", 700):
                t.remove()
    for t in list(fig.texts):
        if t.get_text() in list("ABCDEF") and t.get_fontweight() in ("bold", 700):
            t.remove()
    for lab in [k for k in ax if k in "ABCDEF"]:
        p = ax[lab].get_position()
        fig.text(max(p.x0 - dx, 0.004), min(p.y1 + dy, 0.995), lab,
                 fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()
    return {lab: min(ax[lab].get_position().y1 + dy, 0.995)
            for lab in [k for k in ax if k in "ABCDEF"]}


def _paired(ax, tab, c1, c2, xlabels, title, note, ylab=None, ylim=None):
    # one colour per patient line, set by direction; colouring the dots by tissue
    # instead would put a second, conflicting meaning on the same red
    for _, r in tab.iterrows():
        col = C_UP if r[c2] > r[c1] else C_DOWN
        ax.plot([0, 1], [r[c1], r[c2]], color=col, lw=1.4, alpha=.9, zorder=2)
        ax.scatter([0, 1], [r[c1], r[c2]], s=42, c=col, ec="w", lw=.5, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(xlabels)
    ax.set_xlim(-.36, 1.36)
    if ylim:
        ax.set_ylim(*ylim)
    else:
        ax.margins(y=0.20)
    ax.set_title(title, loc="left")
    ax.text(.5, .02, note, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=TIER[2], style="italic", color="0.25")
    if ylab:
        ax.set_ylabel(ylab)


def build(plt, np, pd, ctx):
    import matplotlib as mpl
    from matplotlib.colors import TwoSlopeNorm
    g = ctx
    S, evr = g["S"], g["evr"]
    meta, comp, dot, vol, g247 = g["meta"], g["comp"], g["dot"], g["vol"], g["g247"]
    Mz, tests, SIG = g["Mz"], g["tests"], g["SIG"]

    # 2 columns x 3 rows (portrait).  A/B, C/D, E/F.
    # Compact 2 x 3.  Axes filled only 55% of the previous canvas; the space went
    # to three copies of the tissue legend, sentence-length panel titles and
    # legend strips hanging below panels.  Legends now sit inside their panels or
    # once at figure level, titles are noun phrases, and the gaps shrink to match.
    # Equal columns.  Colourbars are NOT attached to their panels (that steals
    # width and leaves the column edges ragged); they go in the gutters, placed by
    # align_columns() after layout, so A/C/E share one left and right edge and
    # B/D/F share another.
    # wider at the 11/10/9 ladder: the gutter that holds E's colour bar was sized
    # for 7 pt tick labels and its numbers reached into F's axis title
    fig = plt.figure(figsize=(13.4, 9.4))
    gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.30,
                          width_ratios=[1.0, 1.0],
                          left=0.055, right=0.935, top=0.92, bottom=0.055)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, 0]); axF = fig.add_subplot(gs[2, 1])

    # ---- A: PCA of AOIs ----
    for seg in SEGS:
        for tis, filled in [("eutopic endometrium", False), ("endometriotic lesion", True)]:
            k = ((meta["cell type"] == seg) & (meta["tissue"] == tis)).values
            axA.scatter(S[k, 0], S[k, 1], s=34,
                        facecolors=SEGC[seg] if filled else "none",
                        edgecolors=SEGC[seg], linewidths=1.1, alpha=.9, zorder=3)
    axA.set_xlabel(f"PC1 ({evr[0]:.0f}%)"); axA.set_ylabel(f"PC2 ({evr[1]:.0f}%)")
    # explicit ticks inside the data range: matplotlib keeps off-range tick labels
    # as live Text objects, which collide in the overlap audit though they never render
    x0, x1 = axA.get_xlim(); y0, y1 = axA.get_ylim()
    axA.set_xticks([t for t in np.arange(-60, 41, 20) if x0 < t < x1])
    axA.set_yticks([t for t in np.arange(-20, 31, 10) if y0 < t < y1])
    axA.set_title("AOI structure", loc="center", pad=34)
    h = [plt.Line2D([], [], marker="o", ls="", mfc=SEGC[s_], mec=SEGC[s_], ms=6, label=SEGSHORT[s_])
         for s_ in SEGS]
    axA.legend(handles=h, frameon=False, fontsize=TIER[2], ncol=1, loc="upper left",
               bbox_to_anchor=(-0.015, 1.02), handletextpad=.25, labelspacing=.28)

    # ---- B: gene x AOI heatmap with annotation bars ----
    order = []
    for seg in SEGS:
        for tis in ["eutopic endometrium", "endometriotic lesion"]:
            order += list(meta.index[(meta["cell type"] == seg) & (meta["tissue"] == tis)])
    H = Mz.loc[order, SIG].T.values
    im = axB.imshow(H, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, interpolation="nearest")
    axB.set_yticks(range(len(SIG))); axB.set_yticklabels(SIG, fontsize=TIER[2])
    axB.set_xticks([])
    for i in range(1, 3):
        axB.axvline(i * 20 - .5, color="w", lw=2.0)
    axB.set_title("Nine-gene programme, all 60 AOIs", loc="center", pad=34)
    tr = axB.get_xaxis_transform()
    for i, seg in enumerate(SEGS):
        axB.add_patch(mpl.patches.Rectangle((i * 20 - .5, 1.02), 20, 0.055, transform=tr,
                                            clip_on=False, fc=SEGC[seg], ec="none"))
        axB.text(i * 20 + 9.5, 1.135, SEGSHORT[seg], transform=tr, ha="center",
                 fontsize=TIER[2], color="0.25")
        for j, tis in enumerate(["eutopic endometrium", "endometriotic lesion"]):
            axB.add_patch(mpl.patches.Rectangle((i * 20 + j * 10 - .5, 0.955), 10, 0.055,
                                                transform=tr, clip_on=False, fc=TISC[tis], ec="none"))
    axB.text(-0.155, 1.045, "segment", transform=axB.transAxes, ha="right", fontsize=TIER[2], color="0.4")
    axB.text(-0.155, 0.975, "tissue", transform=axB.transAxes, ha="right", fontsize=TIER[2], color="0.4")
    cbB = fig.add_axes([0.95, 0.75, 0.010, 0.15])   # repositioned in align_columns
    cb = fig.colorbar(im, cax=cbB)
    cb.set_label("z-score", fontsize=TIER[2]); cb.ax.tick_params(labelsize=TIER[2])
    # ---- C: composite by segment and tissue ----
    pos, xt = [], []
    for i, seg in enumerate(SEGS):
        for j, tis in enumerate(["eutopic endometrium", "endometriotic lesion"]):
            x = i * 2.6 + j * 0.9
            v = comp.loc[(comp["cell type"] == seg) & (comp["tissue"] == tis), "score"].values
            bp = axC.boxplot([v], positions=[x], widths=.66, patch_artist=True, showfliers=False,
                             medianprops=dict(color="0.15", lw=1.3),
                             boxprops=dict(fc=TISC[tis], ec="0.35", lw=.8, alpha=.55),
                             whiskerprops=dict(color="0.5", lw=.8), capprops=dict(color="0.5", lw=.8))
            rngj = np.random.default_rng(4 + i * 2 + j)
            axC.scatter(x + rngj.uniform(-.16, .16, len(v)), v, s=16, c=TISC[tis],
                        ec="w", lw=.35, zorder=3)
            pos.append(x)
        xt.append(i * 2.6 + 0.45)
        t = tests[seg]
        axC.text(i * 2.6 + 0.45, 2.02, f"P = {t['p']:.2f}", ha="center", fontsize=TIER[2], color="0.25")
    axC.set_xticks(xt); axC.set_xticklabels([SEGSHORT[s_] for s_ in SEGS])
    axC.set_ylabel("M2 / immunosuppression score\n(per AOI)")
    axC.set_ylim(-1.6, 2.35); axC.margins(x=0.06)
    axC.set_title("Composite score by compartment", loc="center", pad=10)


    # ---- D: volcano inside macrophage AOIs ----
    sig = (vol.p < 0.05) & (vol.lfc.abs() > 0.5)
    axD.scatter(vol.loc[~sig, "lfc"], vol.loc[~sig, "nlp"], s=3, c="#C8CED4", lw=0, rasterized=True)
    axD.scatter(vol.loc[sig & (vol.lfc > 0), "lfc"], vol.loc[sig & (vol.lfc > 0), "nlp"],
                s=4, c="#E2A6A0", lw=0, rasterized=True)
    axD.scatter(vol.loc[sig & (vol.lfc < 0), "lfc"], vol.loc[sig & (vol.lfc < 0), "nlp"],
                s=4, c="#A9C4D8", lw=0, rasterized=True)
    sub = vol[vol.gene.isin(SIG)]
    axD.scatter(sub.lfc, sub.nlp, s=30, c=C_LES, ec="w", lw=.6, zorder=4)
    lab = sub.nlargest(4, "lfc").reset_index(drop=True)
    ytop = vol.nlp.max() * 0.92
    for k, r in lab.iterrows():
        axD.annotate(r.gene, xy=(r.lfc, r.nlp), xytext=(2.05, ytop - k * 0.42),
                     fontsize=TIER[2], color="0.2", ha="left", va="center",
                     arrowprops=dict(arrowstyle="-", lw=.6, color="0.6", shrinkA=1, shrinkB=2))
    axD.axvline(0, color="0.75", lw=.8); axD.axhline(-np.log10(.05), color="0.75", lw=.8, ls="--")
    axD.set_xlabel("log$_2$FC, lesion − eutopic"); axD.set_ylabel("−log$_{10}$ P (paired)")
    axD.set_title("Macrophage AOIs, 18,675 genes", loc="center", pad=10)
    axD.margins(x=0.16, y=0.10)

    # ---- E: dot plot ----
    piv_l = dot.pivot(index="gene", columns="segment", values="lfc")[SEGS].loc[SIG]
    piv_p = dot.pivot(index="gene", columns="segment", values="p")[SEGS].loc[SIG]
    nrm = TwoSlopeNorm(vmin=-1.3, vcenter=0, vmax=1.3)
    for j, seg in enumerate(SEGS):
        for i, gn in enumerate(SIG):
            axE.scatter(j, i, s=18 + 110 * min(-np.log10(max(piv_p.loc[gn, seg], 1e-4)) / 2.0, 1.0),
                        c=[plt.get_cmap("RdBu_r")(nrm(piv_l.loc[gn, seg]))], ec="0.4", lw=.4)
    axE.set_xticks(range(3)); axE.set_xticklabels([SEGSHORT[s_] for s_ in SEGS], rotation=30, ha="right")
    for t, s_ in zip(axE.get_xticklabels(), SEGS):
        t.set_color(SEGC[s_])
    axE.set_yticks(range(len(SIG))); axE.set_yticklabels(SIG, fontsize=TIER[2])
    axE.set_xlim(-.5, 2.5); axE.set_ylim(-.8, len(SIG) + 1.0)   # top band reserved for the size key
    axE.set_title("Per-gene effect by compartment", loc="center", pad=10)
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=nrm)
    cbE = fig.add_axes([0.46, 0.10, 0.010, 0.15])   # repositioned in align_columns
    cb2 = fig.colorbar(sm, cax=cbE)
    # label above the bar, not beside it: a rotated side label reaches into panel F
    cb2.ax.set_title("log$_2$FC", fontsize=TIER[2], pad=4); cb2.ax.tick_params(labelsize=TIER[2])
    for s_, lab in [(18, "P = 1"), (128, "P < 0.01")]:
        axE.scatter([], [], s=s_, c="0.7", ec="0.4", lw=.4, label=lab)
    # size key in the gutter under the colourbar: keeping it inside E forced a
    # wide empty right margin (xlim ran to 3.5 for three columns at 0/1/2)
    # key inside the panel, in a reserved band above the top gene row: below the
    # panel it hung off the figure, in the gutter it hit F's axis title, and at the
    # right it forced the three columns into the left half
    lgE = axE.legend(frameon=False, fontsize=TIER[2], loc="upper right", bbox_to_anchor=(1.0, 1.0),
                     ncol=2, handletextpad=.4, columnspacing=1.2)

    # ---- F: independent cohort ----
    t247 = tests["GSE247695"]
    for j, (cn, tis) in enumerate([("Eutopic", "eutopic endometrium"), ("Lesion", "endometriotic lesion")]):
        v = g247[cn].values
        axF.boxplot([v], positions=[j], widths=.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="0.15", lw=1.3),
                    boxprops=dict(fc=TISC[tis], ec="0.35", lw=.8, alpha=.55),
                    whiskerprops=dict(color="0.5", lw=.8), capprops=dict(color="0.5", lw=.8))
        rngj = np.random.default_rng(9 + j)
        axF.scatter(j + rngj.uniform(-.13, .13, len(v)), v, s=22, c=TISC[tis], ec="w", lw=.4, zorder=3)
    axF.set_xticks([0, 1]); axF.set_xticklabels(["Eutopic", "Lesion"])
    axF.set_xlim(-.62, 1.62)
    axF.set_ylabel("M2 score (myeloid cells)")
    axF.margins(y=0.20)
    axF.set_title("Independent scRNA cohort (GSE247695)", loc="center", pad=10)
    axF.text(.5, .02, f"{t247['up']}/4 patients up · P = {t247['p']:.2f}", transform=axF.transAxes,
             ha="center", va="bottom", fontsize=TIER[2], style="italic", color="0.25")

    # tissue key under B, where the colours are defined (B's tissue annotation bar)
    # and reused by C and F.  It was at the figure's top-left corner, which reads as
    # belonging to A - but A encodes tissue by open/filled symbols, not by fill colour.
    axB.legend(handles=[plt.Line2D([], [], marker="s", ls="", ms=6, mfc=TISC[t], mec="none",
                                   label="eutopic endometrium" if "eutopic" in t else "endometriotic lesion")
                        for t in ["eutopic endometrium", "endometriotic lesion"]],
               frameon=False, fontsize=TIER[1], ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, -0.02), handletextpad=.35, columnspacing=2.0)

    ax = dict(A=axA, B=axB, C=axC, D=axD, E=axE, F=axF, cbB=cbB, cbE=cbE, lgE=lgE)
    g["normalize"](fig)
    for a in (axB, axE):
        for t in a.get_yticklabels():
            t.set_fontsize(TIER[2])
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

def align_columns(fig, ax, gap=0.012, cbw=0.010):
    """Give every panel in a column the same left and right edge, and park the
    colourbars in the gutters.

    Attaching a colourbar with ax= shrinks that panel only, so B was narrower
    than D and F and E narrower than A and C.  Here each column is squared off
    first, then cbE goes in the inter-column gutter and cbB in the right margin.
    """
    fig.canvas.draw()
    for keys in (["A", "C", "E"], ["B", "D", "F"]):
        x0 = min(ax[k].get_position().x0 for k in keys)
        x1 = max(ax[k].get_position().x1 for k in keys)
        for k in keys:
            p = ax[k].get_position()
            ax[k].set_position([x0, p.y0, x1 - x0, p.height])
    pe = ax["E"].get_position(); pb = ax["B"].get_position()
    ax["cbE"].set_position([pe.x1 + gap, pe.y0, cbw, pe.height])
    ax["cbB"].set_position([pb.x1 + gap, pb.y0, cbw, pb.height])
    fig.canvas.draw()
    return dict(left_edges=[round(ax[k].get_position().x0, 5) for k in "ACE"],
                left_right=[round(ax[k].get_position().x1, 5) for k in "ACE"],
                right_edges=[round(ax[k].get_position().x0, 5) for k in "BDF"],
                right_right=[round(ax[k].get_position().x1, 5) for k in "BDF"])

def align_row_frames(fig, axes_list):
    """Give a row of panels an identical axes rectangle.

    align_row_bottoms() equalises each panel's LOWEST RENDERED ELEMENT, so a
    panel carrying a legend or rotated tick labels ends up with a higher frame
    than its neighbour - the plot boxes then visibly disagree even though the
    ink stops on one line.  Here the frames themselves are made to coincide and
    the hanging furniture is allowed to extend past them.
    """
    fig.canvas.draw()
    y0 = max(a.get_position().y0 for a in axes_list)
    y1 = min(a.get_position().y1 for a in axes_list)
    for a in axes_list:
        p = a.get_position()
        a.set_position([p.x0, y0, p.width, y1 - y0])
    fig.canvas.draw()
    return dict(y0=[round(a.get_position().y0, 5) for a in axes_list],
                y1=[round(a.get_position().y1, 5) for a in axes_list])
