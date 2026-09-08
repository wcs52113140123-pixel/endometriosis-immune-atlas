"""Build Figure 6 (is the pooling effect peculiar to endometriosis?) at print width.

Three panels only.  This figure exists to show that the calibration failure
behind the endometriosis result is general, so it carries the calibration, the
design question a reader will ask next (does a larger cohort fix it), and the
consequence — and deliberately not the tissue-by-tissue detail, which belongs to
a methodological paper rather than a reproductive-medicine one.

Panel types follow the same convention as Figures 2-5: violin with median and
interquartile range, line plot with a nominal reference, and a bar.
"""

import numpy as np

TIER = (11.0, 10.0, 9.0)
C_POOL, C_DON = "#C0392B", "#2E6DA4"
MG = "#5A6570"


def build(plt, np, pd, ctx):
    B = ctx["B"]
    fig = plt.figure(figsize=(7.09, 2.9))
    gs = fig.add_gridspec(1, 3, wspace=0.52, left=0.095, right=0.985, top=0.84, bottom=0.235)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])

    # ---- A: calibration of the two tests ----
    for i, (col, c) in enumerate([("fpr_pooled", C_POOL), ("fpr_donor", C_DON)]):
        v = B[col].values
        parts = axA.violinplot([v], positions=[i], widths=.82, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(c); b.set_alpha(.45); b.set_edgecolor("0.4"); b.set_linewidth(.7)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        axA.vlines(i, q1, q3, color="0.2", lw=3.0, zorder=4)
        axA.scatter([i], [med], s=15, c="w", ec="0.2", lw=.8, zorder=5)
        axA.text(i, 1.02, f"{med:.2f}", ha="center", va="bottom", fontsize=TIER[2], color=c)
    axA.axhline(0.05, color="0.35", ls="--", lw=.9)
    axA.text(-0.42, 0.075, "nominal 0.05", fontsize=TIER[2], color=MG, va="bottom", ha="left")
    axA.set_xticks([0, 1]); axA.set_xticklabels(["Pooled\ncells", "Donor\nlevel"])
    axA.set_xlim(-.6, 1.6); axA.set_ylim(0, 1.12)
    axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axA.set_ylabel("False-positive rate")
    axA.set_title("Calibration", loc="center", pad=8)

    # ---- B: does a larger cohort fix it ----
    bins = [0, 10, 20, 40, 80, 400]
    g = pd.cut(B.n_donors, bins)
    m = B.groupby(g, observed=True)[["fpr_pooled", "fpr_donor"]].median()
    centres = [iv.mid for iv in m.index]
    for col, c, lab in [("fpr_pooled", C_POOL, "pooled cells"),
                        ("fpr_donor", C_DON, "donor level")]:
        axB.plot(centres, m[col].values, "-o", ms=3.5, lw=1.5, color=c, label=lab)
    axB.axhline(0.05, color="0.35", ls="--", lw=.9)
    axB.set_xscale("log"); axB.set_xlim(3, 300); axB.set_ylim(0, 1)
    axB.set_xlabel("Donors in the study"); axB.set_ylabel("Median false-positive rate")
    axB.legend(frameon=False, loc="center right", handletextpad=.3, handlelength=1.4)
    axB.set_title("Cohort size does not help", loc="center", pad=8)

    # ---- C: consequence on observed data ----
    frac = [(B.p_pooled < .05).mean(), (B.p_donor < .05).mean()]
    upheld = ((B.p_pooled < .05) & (B.p_donor < .05)).sum() / (B.p_pooled < .05).sum()
    axC.bar([0, 1], frac, color=[C_POOL, C_DON], width=.62, alpha=.85)
    for i, v in enumerate(frac):
        axC.text(i, v + .025, f"{v*100:.0f}%", ha="center", fontsize=TIER[2])
    axC.set_xticks([0, 1]); axC.set_xticklabels(["Pooled\ncells", "Donor\nlevel"])
    axC.set_xlim(-.6, 1.6); axC.set_ylim(0, 1.0)
    axC.set_ylabel("Comparisons called significant")
    axC.set_title("Observed data", loc="center", pad=8)
    # over the shorter bar: centred it sat on top of the taller one
    axC.text(0.98, 0.62, f"{(1-upheld)*100:.0f}% of pooled\nhits not upheld",
             transform=axC.transAxes, ha="right", va="top", fontsize=TIER[2],
             style="italic", color="0.25")

    ax = dict(A=axA, B=axB, C=axC)
    for a in ax.values():
        a.xaxis.label.set_fontsize(TIER[0]); a.yaxis.label.set_fontsize(TIER[0])
        a.title.set_fontsize(TIER[0])
        for t in a.get_xticklabels() + a.get_yticklabels():
            t.set_fontsize(TIER[2])
    return fig, ax


def align_row_frames(fig, axes_list):
    fig.canvas.draw()
    y0 = max(a.get_position().y0 for a in axes_list)
    y1 = min(a.get_position().y1 for a in axes_list)
    for a in axes_list:
        p = a.get_position()
        a.set_position([p.x0, y0, p.width, y1 - y0])
    fig.canvas.draw()


def level_titles(fig, ax):
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
    fh = fig.get_size_inches()[1]
    bot = {k: T.transform(ax[k].title.get_window_extent(r))[0][1] for k in ax}
    tgt = max(bot.values())
    for k in ax:
        ax[k].set_title(ax[k].get_title(), loc="center",
                        pad=(tgt - ax[k].get_position().y1) * fh * 72)
    fig.canvas.draw()


def place_letters(fig, ax, dy=0.034, dx=0.055):
    for t in list(fig.texts):
        if t.get_text() in list("ABC") and t.get_fontweight() in ("bold", 700):
            t.remove()
    fig.canvas.draw()
    for k in "ABC":
        p = ax[k].get_position()
        fig.text(max(p.x0 - dx, 0.004), min(p.y1 + dy, 0.995), k,
                 fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()


def real_overlaps(fig, mpl):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible() and abs(t.get_rotation()) < 1]
    out = []
    for i, (a, ba) in enumerate(tx):
        for b, bb in tx[i + 1:]:
            if min(ba.x1, bb.x1) - max(ba.x0, bb.x0) > 2 and min(ba.y1, bb.y1) - max(ba.y0, bb.y0) > 2:
                out.append((a.get_text()[:22], b.get_text()[:22]))
    return out
