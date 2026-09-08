"""Build Figure 1 (atlas overview, 5 panels).

Layout is a nested gridspec: an outer 1x2 split whose wspace is set wide enough
that panel C's long row labels never reach panel B, and a 2x2 inner grid for
A/B/D/E.  After the first draw the left block is repositioned so its lowest
rendered element (D's axis title / E's legend) lines up with panel C's lowest
element (its lineage legend).

Expects these names in the caller's namespace: plt, mpl, np, pd, obs, um, sub,
LIN, LINLAB, LINC, DSC, INLINE, OFF, cent, xr, yr, compC, nfo, ROWS_T, ROWS_C,
ds, sb, cols, bio, META_GREY, PL, normalize.
"""

PT = 1.5            # UMAP marker size
TIER = (11.0, 10.0, 9.0)   # unified ladder across Figures 1-4


def build(plt, mpl, np, ctx):
    g = ctx
    # wider at the 11/10/9 ladder: at 9/8/7 the cohort legend fitted under panel B
    # and D's x ticks cleared E's rotated metric names; both collided at the
    # larger sizes
    # PRINT-SIZE VARIANT: 180 mm journal column = 7.09 in.  The wide version puts
    # C in a full-height right column; at this width its long row labels would
    # leave the bars about an inch across, so C spans the middle row instead.
    fig = plt.figure(figsize=(7.09, 9.3))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.25, 0.95], hspace=0.66,
                             left=0.225, right=0.985, top=0.955, bottom=0.055)
    top = outer[0].subgridspec(1, 2, wspace=0.24)
    bot_gs = outer[2].subgridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.52)
    axA = fig.add_subplot(top[0, 0]); axB = fig.add_subplot(top[0, 1])
    axD = fig.add_subplot(bot_gs[0, 0]); axE = fig.add_subplot(bot_gs[0, 1])
    axC = fig.add_subplot(outer[1])

    # ---- A: lineages ----
    for l in g["LIN"]:
        k = g["sub"][g["obs"].celltype.values[g["sub"]] == l]
        axA.scatter(g["um"][k, 0], g["um"][k, 1], s=PT, c=g["LINC"][l], lw=0, alpha=.60, rasterized=True)
    counts = g["obs"].celltype.value_counts()
    for l in g["INLINE"]:
        axA.text(*g["cent"][l], g["LINLAB"][l], fontsize=TIER[1], ha="center", va="center",
                 bbox=dict(fc="white", ec="none", alpha=.75, pad=1.2))
    for l, (dx, dy) in g["OFF"].items():
        x0, y0 = g["cent"][l]
        # counts dropped from the leader labels at this width; they are in the caption
        axA.annotate(g["LINLAB"][l], xy=(x0, y0),
                     xytext=(x0 + dx * g["xr"], y0 + dy * g["yr"]), fontsize=TIER[1],
                     ha="center", va="center", bbox=dict(fc="white", ec="none", alpha=.75, pad=1.2),
                     arrowprops=dict(arrowstyle="-", lw=.7, color="0.45", shrinkA=0, shrinkB=2))

    # ---- B: cohorts ----
    for d, c in g["DSC"].items():
        k = g["sub"][g["obs"].dataset.values[g["sub"]] == d]
        axB.scatter(g["um"][k, 0], g["um"][k, 1], s=PT, c=c, lw=0, alpha=.48, rasterized=True)
    for ax, ttl in [(axA, "Eight lineages, 344,510 cells"),
                    (axB, "Cohorts intermix after integration")]:
        ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.01)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(ttl, loc="center", pad=10)
    MG = g["META_GREY"]
    axA.annotate("", xy=(0.12, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=.8, color=MG))
    axA.annotate("", xy=(0.02, 0.12), xytext=(0.02, 0.02), xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=.8, color=MG))
    axA.text(0.135, 0.014, "UMAP 1", transform=axA.transAxes, fontsize=TIER[2], color=MG)
    axA.text(0.008, 0.135, "UMAP 2", transform=axA.transAxes, fontsize=TIER[2], color=MG, rotation=90)
    lgB = axB.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=5, mfc=c, mec="none", label=d)
                              for d, c in g["DSC"].items()],
                     frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.015), ncol=3,
                     fontsize=TIER[1], handletextpad=.3, columnspacing=1.1, labelspacing=.35)

    # ---- C: composition ----
    rows_t, rows_c = g["ROWS_T"], g["ROWS_C"]
    ypos = list(range(len(rows_t))) + [len(rows_t) + 0.9, len(rows_t) + 1.9]
    labels = rows_t + rows_c
    bot = np.zeros(len(labels))
    for l in g["LIN"]:
        axC.barh(ypos, g["compC"][l].values, left=bot, color=g["LINC"][l], height=.72, label=g["LINLAB"][l])
        bot += g["compC"][l].values
    for y, lab in zip(ypos, labels):
        n = g["nfo"].loc[lab]
        axC.text(101.5, y, f"{n.cells/1000:.0f}k · {n.samples}", va="center", fontsize=TIER[2], color="0.35")
    axC.axhline(len(rows_t) - 0.05, color="0.72", lw=.8, ls="--")
    axC.text(50, len(rows_t) + 0.32, "culture model — composition set by protocol",
             ha="center", va="center", fontsize=TIER[2], color=MG, style="italic")
    axC.set_yticks(ypos)
    # single-line row labels: two-line labels at 9 pt touch when ten rows share
    # 2.4 in, so the site and status are joined with a comma instead
    _SHORTROW = {"Endometrium — control": "Endometrium, control",
                 "Eutopic endometrium — endometriosis": "Eutopic endo., endometriosis",
                 "Endometrium — condition not stated": "Endometrium, not stated",
                 "Menstrual endometrium — not stated": "Menstrual endo., not stated",
                 "Peritoneal lesion — endometriosis": "Peritoneal lesion, endo.",
                 "Ovary — endometriosis": "Ovary, endometriosis",
                 "OCCC — malignant": "OCCC, malignant", "EAOC — malignant": "EAOC, malignant",
                 "Organoid — control": "Organoid, control",
                 "Organoid — endometriosis": "Organoid, endometriosis"}
    axC.set_yticklabels([_SHORTROW.get(l, l) for l in labels], fontsize=TIER[2])
    axC.invert_yaxis(); axC.set_xlim(0, 100); axC.set_ylim(len(rows_t) + 2.6, -0.8)
    axC.set_xlabel("Cells (%)")
    axC.set_title("Composition by site and disease status", loc="center", pad=8)
    axC.text(1.005, 1.02, "cells · samples", transform=axC.transAxes, fontsize=TIER[2],
             color="0.35", ha="left")
    # four columns on one line: the two-column block of the wide version is as
    # tall as a panel here and reached into the bottom row
    lgC = axC.legend(frameon=False, fontsize=TIER[2], ncol=4, loc="upper center",
                     bbox_to_anchor=(0.5, -0.155), handlelength=.9, columnspacing=1.0,
                     handletextpad=.4, labelspacing=.3)

    # ---- D: cohort scale ----
    ds = g["ds"]
    y = np.arange(len(ds))
    axD.barh(y, ds.cells.values / 1000, color=[g["DSC"][d] for d in ds.index], height=.66)
    for i, (cc, ss) in enumerate(zip(ds.cells.values, ds.samples.values)):
        axD.text(cc / 1000 + 3, i, f"{cc/1000:.0f}k · {ss}", va="center", fontsize=TIER[1], color="0.30")
    axD.set_yticks(y); axD.set_yticklabels(ds.index)
    axD.set_xlabel("Cells (thousands)"); axD.set_xlim(0, 172); axD.margins(y=0.08)
    # same pad as E so the two bottom-row titles sit on one line
    axD.set_title("Cohort contributions", loc="center", pad=14)

    # ---- E: scIB ----
    sb, cols, bio = g["sb"], g["cols"], g["bio"]
    x = np.arange(len(cols)); w = .38
    axE.bar(x - w / 2, sb.loc["Unintegrated", cols].values, w, color="#9AA5AD", label="Unintegrated PCA")
    axE.bar(x + w / 2, sb.loc["Harmony", cols].values, w, color="#1F6FB2", label="Harmony")
    axE.axvline(len(bio) - .5, color="0.78", lw=.8, ls="--")
    axE.set_xticks(x); axE.set_xticklabels(cols, rotation=38, ha="right")
    axE.set_ylabel("scIB score"); axE.set_ylim(0, 1.06); axE.margins(x=0.02)
    axE.set_title("scIB integration benchmark", loc="center", pad=14)
    _y0, _y1 = axE.get_ylim()
    axE.set_yticks([v for v in np.arange(0, 1.01, 0.2) if _y0 <= v <= _y1])
    axE.text((len(bio) - 1) / 2, 1.008, "biological conservation", ha="center", fontsize=TIER[2],
             color=MG, transform=axE.get_xaxis_transform())
    axE.text(len(bio) + 1.5, 1.008, "batch correction", ha="center", fontsize=TIER[2],
             color=MG, transform=axE.get_xaxis_transform())
    lgE = axE.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.44), ncol=2,
                     fontsize=TIER[1], handlelength=1.1, columnspacing=1.6)

    for ax, l in [(axA, "A"), (axB, "B"), (axC, "C"), (axD, "D"), (axE, "E")]:
        g["PL"](ax, l)
    g["normalize"](fig)
    for t in axC.get_yticklabels():
        t.set_fontsize(TIER[2])
    return fig, dict(A=axA, B=axB, C=axC, D=axD, E=axE, lgB=lgB, lgC=lgC, lgE=lgE)


def place_letters(fig, ax, dy=0.030, dx=0.046):
    """Put every panel letter at a constant figure-space offset from its axes.

    Placing letters in axes-fraction coordinates makes them drift: the same
    fraction is a larger absolute offset on a taller panel, so C (full height)
    floated above A and B.  A fixed figure-space offset puts all letters of a
    row on one baseline.
    """
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
    ys = {lab: min(ax[lab].get_position().y1 + dy, 0.995) for lab in "ABCDE"}
    return ys


def align_bottoms(fig, ax, plt):
    """Extend the left block down so its lowest rendered element meets C's."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    T = fig.transFigure.inverted()

    def bb(a):
        return T.transform(a.get_window_extent(r))

    c_bot = bb(ax["lgC"])[0][1]
    gap = (bb(ax["lgB"])[1][1] - bb(ax["lgB"])[0][1]) + 0.055
    top = ax["A"].get_position().y1
    low_now = min(bb(ax["D"].xaxis.label)[0][1], bb(ax["lgE"])[0][1])
    overhang = bb(ax["D"])[0][1] - low_now
    target = c_bot + overhang
    h = (top - target - gap) / 2
    for k in ("A", "B"):
        p = ax[k].get_position(); ax[k].set_position([p.x0, top - h, p.width, h])
    for k in ("D", "E"):
        p = ax[k].get_position(); ax[k].set_position([p.x0, target, p.width, h])
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    low = min(bb(ax["D"].xaxis.label)[0][1], bb(ax["lgE"])[0][1])
    return dict(mismatch=abs(low - c_bot), panel_height=h,
                B_legend_right=bb(ax["lgB"])[1][0], C_left=bb(ax["C"])[0][0],
                C_label_left=min(T.transform(t.get_window_extent(r))[0][0]
                                 for t in ax["C"].get_yticklabels()),
                B_right=bb(ax["B"])[1][0])


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
