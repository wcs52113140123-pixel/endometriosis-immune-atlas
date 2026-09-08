"""Build Figure 4 (cell-cell communication around Treg, 6 panels).

PLOT TYPES follow a caption census of the cell-communication literature, not a
general single-cell one: 218 communication-related figure captions from 74
papers published 2023-2026 in Nat Commun, J Immunother Cancer, Cell Rep, Genome
Biol, Sci Adv, Genome Med, Cell Rep Med, Gut, Nat Cancer, Immunity, J Clin
Invest, eLife and Hepatology.  Panel types used by those papers:

    bubble / dot plot of LR pairs   47% of papers
    bar chart (information flow)    32%
    circle plot (node-edge)         30%
    chord / circos                  24%
    scatter (incoming vs outgoing)  24%
    ranked point / lollipop          1%   <- what the previous version used

The previous version of this figure was two lollipops, two heatmaps and two dot
plots, which is both the rarest type in this literature and a repeat of Figure
3's grammar.  This version uses the five communication-specific types above and
keeps one expression dot plot to ground the inference in measured expression.

  A  circle plot of significant interactions among all cell types, Treg marked
  B  chord ribbons of what reaches Treg, by sending cell type
  C  bubble plot of ligand-receptor pairs directed at Treg, by sender
  D  information flow: interaction strength each sender directs at Treg
  E  signalling role: outgoing versus incoming strength per cell type
  F  expression of the ligands and receptors involved, by immune state

Layout follows the standard settled on Figure 3: 2 columns x 3 rows, equal
column widths, axes FRAMES levelled per row, titles centred and level, colour
bars parked as free axes beneath their panel.
"""

import numpy as _np
from adjustText import adjust_text

TIER = (11.0, 10.0, 9.0)   # unified ladder across Figures 1-4
MG = "0.35"
SHORT = {"Stromal_Fib": "Stromal", "SmoothMuscle": "Smooth mus.", "Endothelial": "Endothel.",
         "Epithelial": "Epithel.", "Macrophage": "Macro."}
CTC = {"Treg": "#C0392B", "CD4T": "#2E86AB", "CD8T": "#4C7BA6", "NK": "#59A14F",
       "Macrophage": "#E1A730", "Mono": "#B07A5A", "Stromal_Fib": "#7F8FA6",
       "SmoothMuscle": "#9C6DAB", "Endothelial": "#1F6FB2", "Epithelial": "#8C564B",
       "B": "#7B9E89", "Plasma": "#C49A6C", "Mast": "#D08770", "pDC": "#6C7A89"}


def _ring(n, start=90.0):
    """Angles (degrees) for n nodes on a circle, clockwise from the top."""
    return [start - 360.0 * i / n for i in range(n)]


def _bezier(p0, p1, curve=0.55, npts=60):
    """Quadratic Bezier from p0 to p1 bending toward the origin."""
    mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
    c = (mx * (1 - curve), my * (1 - curve))
    t = _np.linspace(0, 1, npts)[:, None]
    p0, p1, c = _np.array(p0), _np.array(p1), _np.array(c)
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * c + t ** 2 * p1


def _circle_plot(ax, net, focus, plt, np, top_edges=60):
    cts = list(net.index)
    ang = {c: a for c, a in zip(cts, _ring(len(cts)))}
    pos = {c: (np.cos(np.deg2rad(a)), np.sin(np.deg2rad(a))) for c, a in ang.items()}
    edges = [(s, t, int(net.loc[s, t])) for s in cts for t in net.columns
             if s != t and net.loc[s, t] > 0]
    edges.sort(key=lambda e: e[2])
    edges = edges[-top_edges:]
    wmax = max(e[2] for e in edges)
    for s, t, w in edges:
        hit = focus in (s, t)
        col = CTC.get(s, "0.6")
        b = _bezier(pos[s], pos[t])
        ax.plot(b[:, 0], b[:, 1], lw=0.25 + 2.4 * w / wmax, color=col,
                alpha=.85 if hit else .18, zorder=3 if hit else 2,
                solid_capstyle="round")
    # only the single most-vertical node in each hemisphere gets a horizontal
    # label; making every near-vertical one horizontal made neighbours collide
    _sv = {c: _np.sin(_np.deg2rad(ang[c])) for c in cts}
    _hz = {max((c for c in cts if _sv[c] > 0), key=lambda c: _sv[c]),
           min((c for c in cts if _sv[c] < 0), key=lambda c: _sv[c])}
    tot = net.sum(1) + net.sum(0)
    for _i, c in enumerate(cts):
        x, y = pos[c]
        ax.scatter([x], [y], s=26 + 150 * float(tot[c] / tot.max()),
                   c=CTC.get(c, "0.6"), ec="w", lw=.7, zorder=5)
        a = ang[c]
        # a node near the top or bottom gets a horizontal label: rotated to the
        # radius it would stand upright and run into the panel title
        near_v = c in _hz
        kw = (dict(rotation=0, ha="center", va="bottom" if y > 0 else "top")
              if near_v else
              dict(rotation=a if -90 < a < 90 else a - 180, rotation_mode="anchor",
                   ha="left" if -90 < a < 90 else "right", va="center"))
        # alternate the label radius: at 180 mm the panels are narrow enough that
        # neighbouring radial labels in the crowded lower-left arc collide
        _rad = 1.13 if _i % 2 == 0 else 1.30
        ax.text(x * _rad, y * (1.22 if near_v else _rad), SHORT.get(c, c), fontsize=TIER[1],
                color="#C0392B" if c == focus else "0.25",
                fontweight="bold" if c == focus else "normal", **kw)
    # adjustable="datalim": with the default box adjustment matplotlib SHRINKS the
    # axes to the data aspect, so the circle came out small and the panel no longer
    # shared its left and right edges with the rectangular panels below
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _chord(ax, ser, focus, plt, np):
    """Ribbons from each sender to the focus node, width proportional to count."""
    send = [c for c in ser.index if c != focus and ser[c] > 0]
    vals = np.array([ser[c] for c in send], float)
    total = vals.sum() + vals.sum()          # senders occupy half the ring, focus the other
    gap = 2.0
    sizes = vals / total * (360 - gap * (len(send) + 1))
    # rotated 40 degrees: with the focus arc starting at 90 the small sender
    # segments land at the bottom of the ring, where their near-vertical labels set
    # how far the panel's ink reaches below its axes and hold the whole row apart
    a = 130 - gap / 2
    spans = {}
    for c, s_ in zip(send, sizes):
        spans[c] = (a - s_, a); a -= s_ + gap
    focus_span = (a - vals.sum() / total * (360 - gap * (len(send) + 1)), a)

    def arc(a0, a1, r0=1.0, r1=1.09, col="0.6", **kw):
        th = np.linspace(np.deg2rad(a0), np.deg2rad(a1), 40)
        ax.fill(np.r_[r0 * np.cos(th), (r1 * np.cos(th))[::-1]],
                np.r_[r0 * np.sin(th), (r1 * np.sin(th))[::-1]], color=col, lw=0, **kw)

    _mids = {c: (spans[c][0] + spans[c][1]) / 2 for c in send}
    _sv = {c: np.sin(np.deg2rad(_mids[c])) for c in send}
    _pos = [c for c in send if _sv[c] > 0]; _neg = [c for c in send if _sv[c] < 0]
    _hz = set()
    if _pos: _hz.add(max(_pos, key=lambda c: _sv[c]))
    if _neg: _hz.add(min(_neg, key=lambda c: _sv[c]))
    fa0, fa1 = focus_span
    arc(fa0, fa1, col=CTC[focus])
    fw = (fa1 - fa0)
    cur = fa1
    for _i, c in enumerate(send):
        a0, a1 = spans[c]
        arc(a0, a1, col=CTC.get(c, "0.6"))
        share = ser[c] / vals.sum() * fw
        f1, f0 = cur, cur - share
        cur = f0
        th_s = np.deg2rad(np.linspace(a0, a1, 24)); th_f = np.deg2rad(np.linspace(f0, f1, 24))
        outer = np.c_[np.cos(th_s), np.sin(th_s)]
        inner = np.c_[np.cos(th_f), np.sin(th_f)]
        b1 = _bezier(outer[-1], inner[0], curve=0.75, npts=28)
        b2 = _bezier(inner[-1], outer[0], curve=0.75, npts=28)
        poly = np.vstack([outer, b1, inner, b2])
        ax.fill(poly[:, 0], poly[:, 1], color=CTC.get(c, "0.6"), alpha=.42, lw=0)
        mid = (a0 + a1) / 2
        # every chord label stays rotated: the horizontal treatment collided with
        # its neighbour in the crowded small-sender cluster at the bottom
        near_v = False
        kw = (dict(rotation=0, ha="center", va="bottom" if np.sin(np.deg2rad(mid)) > 0 else "top")
              if near_v else
              dict(rotation=mid if -90 < mid < 90 else mid - 180, rotation_mode="anchor",
                   ha="left" if -90 < mid < 90 else "right", va="center"))
        _rad = 1.13 if _i % 2 == 0 else 1.30
        ax.text(_rad * np.cos(np.deg2rad(mid)),
                (1.26 if near_v else _rad) * np.sin(np.deg2rad(mid)),
                SHORT.get(c, c), fontsize=TIER[1], color="0.25", **kw)
    mid = (fa0 + fa1) / 2
    ax.text(1.14 * np.cos(np.deg2rad(mid)), 1.14 * np.sin(np.deg2rad(mid)), focus,
            fontsize=TIER[0], fontweight="bold", color=CTC[focus],
            ha="left" if -90 < mid < 90 else "right", va="center",
            rotation=mid if -90 < mid < 90 else mid - 180, rotation_mode="anchor")
    # adjustable="datalim": with the default box adjustment matplotlib SHRINKS the
    # axes to the data aspect, so the circle came out small and the panel no longer
    # shared its left and right edges with the rectangular panels below
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def build(plt, np, pd, ctx):
    g = ctx
    net, bub, inflow, role = g["net"], g["bub"], g["inflow"], g["role"]
    me, fe = g["expr_mean"], g["expr_frac"]

    # PRINT-SIZE VARIANT: 180 mm journal column = 7.09 in
    # canvas set so the exported width is ~198 mm, the median of the five main
    # figures: a figure submitted wider is scaled further by the publisher and its
    # type prints smaller, so equal widths give equal printed type
    fig = plt.figure(figsize=(7.00, 9.6))
    # row 1 holds two equal-aspect circular panels, which need less height than the
    # rectangular panels below; without this the row leaves a band of white
    gs = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.72, width_ratios=[1.0, 1.0],
                          height_ratios=[0.92, 1.0, 1.0],
                          left=0.075, right=0.885, top=0.945, bottom=0.045)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, 0]); axF = fig.add_subplot(gs[2, 1])

    # ---- A: circle plot of the whole network ----
    _circle_plot(axA, net, "Treg", plt, np)
    axA.set_title("Significant interactions among all cell types", loc="center", pad=22)

    # ---- B: chord of what reaches Treg ----
    _chord(axB, net["Treg"], "Treg", plt, np)
    axB.set_title("What reaches Treg, by sending cell type", loc="center", pad=22)

    # ---- C: bubble plot of LR pairs directed at Treg ----
    pairs = (bub.groupby("pair").lr_means.max().sort_values().index.tolist())
    senders = [c for c in net.index if c in set(bub.source)]
    px = {p: i for i, p in enumerate(pairs)}; sx = {s: i for i, s in enumerate(senders)}
    scC = axC.scatter([sx[r.source] for r in bub.itertuples()],
                      [px[r.pair] for r in bub.itertuples()],
                      s=[6 + 14 * r.nlp for r in bub.itertuples()],
                      c=[r.lr_means for r in bub.itertuples()],
                      cmap="Reds", ec="0.4", lw=.3, vmin=1.4, vmax=3.0)
    axC.set_xticks(range(len(senders)))
    axC.set_xticklabels([SHORT.get(s, s) for s in senders], rotation=42, ha="right")
    axC.set_yticks(range(len(pairs))); axC.set_yticklabels(pairs)
    axC.set_xlim(-.7, len(senders) - .3); axC.set_ylim(-.8, len(pairs) - .2)
    axC.grid(axis="both", color="0.92", lw=.5, zorder=0); axC.set_axisbelow(True)
    axC.set_title("Ligand–receptor pairs directed at Treg", loc="center", pad=10)
    for s_, lab in [(6 + 14 * 2, "2"), (6 + 14 * 5, "5")]:
        axC.scatter([], [], s=s_, c="0.75", ec="0.4", lw=.3, label=lab)
    # keys ABOVE the axes, not inside: at 180 mm a reserved in-panel band costs
    # about two of the fourteen data rows and the keys still sat on the top row
    # keys to the RIGHT of C, matching F; the column gutter is widened to clear
    # D's row labels
    axC.legend(frameon=False, fontsize=TIER[2], loc="upper left", bbox_to_anchor=(1.015, 1.02),
               ncol=1, handletextpad=.5, labelspacing=.9,
               title="−log$_{10}$\nspecificity rank", title_fontsize=TIER[2])
    cbC = axC.inset_axes([1.055, 0.02, 0.045, 0.36])
    cb = fig.colorbar(scC, cax=cbC, orientation="vertical", ticks=[1.5, 2.0, 2.5, 3.0])
    # both label and ticks on the OUTER side of the bar: on the inner flank the
    # vertical label sat against the panel it belongs to
    cb.set_label("expression mean", fontsize=TIER[2], rotation=90, labelpad=2)
    cb.ax.yaxis.set_ticks_position("right"); cb.ax.yaxis.set_label_position("right")
    cb.ax.tick_params(labelsize=TIER[2], length=2, pad=1.5)

    # ---- D: information flow toward Treg ----
    inf = inflow[inflow.strength > 0]
    y = np.arange(len(inf))
    axD.barh(y, inf.strength.values, color=[CTC.get(c, "0.6") for c in inf.index], height=.68)
    for i, (st, n) in enumerate(zip(inf.strength.values, inf.n.values)):
        axD.text(st + inf.strength.max() * 0.015, i, f"{int(n)}", va="center",
                 fontsize=TIER[2], color="0.30")
    axD.set_yticks(y); axD.set_yticklabels([SHORT.get(c, c) for c in inf.index])
    axD.set_xlabel("Summed expression mean of significant pairs")
    axD.set_xlim(0, inf.strength.max() * 1.14); axD.margins(y=0.05)
    axD.set_title("Information flow directed at Treg", loc="center", pad=10)

    # ---- E: signalling role ----
    for c in role.index:
        axE.scatter(role.loc[c, "outgoing"], role.loc[c, "incoming"],
                    s=26 + 150 * float((role.loc[c, "n_out"] + role.loc[c, "n_in"]) /
                                       (role.n_out + role.n_in).max()),
                    c=CTC.get(c, "0.6"), ec="w", lw=.7, zorder=3)
    # extra headroom so the repelled labels have somewhere to go
    lim = float(max(role.outgoing.max(), role.incoming.max())) * 1.14
    axE.plot([0, lim], [0, lim], ls="--", lw=.7, color="0.72", zorder=1)
    # adjustText repels each label from BOTH the other labels and the markers,
    # then draws a short leader where it had to travel.  A fixed per-point offset
    # cannot work here: seven of the fourteen cell types sit within about twenty
    # units of a neighbour, so any constant offset puts a label on someone's disc.
    # only the cell types the text discusses are labelled.  Fourteen labels in a
    # 2.3 in scatter cannot be placed without either covering a disc or colliding
    # with each other: a label is about 95 data units wide at this size (median of
    # the seven drawn), while 10 of the 14 cell types lie closer than that to their
    # nearest neighbour (nearest-neighbour distances 15.5, 15.5, 29.9, 29.9, 53.8,
    # 53.8, 54.9, 55.6, 56.6, 69.1, 101, 122, 122, 226).  A sweep of repulsion
    # strengths from 1.18 to 1.90 never reached zero disc coverage with all
    # fourteen.  The rest are identified by colour, which is shared with D.
    LAB_E = ["Mono", "Macrophage", "Endothelial", "Stromal_Fib", "Treg", "CD8T", "NK"]
    annsE = [axE.text(role.loc[c, "outgoing"], role.loc[c, "incoming"], SHORT.get(c, c),
                      fontsize=TIER[2], ha="center", va="center",
                      color="#C0392B" if c == "Treg" else "0.3",
                      fontweight="bold" if c == "Treg" else "normal",
                      bbox=dict(fc="white", ec="none", alpha=.75, pad=0.6))
             for c in LAB_E]
    axE.set_xlim(0, lim); axE.set_ylim(0, lim)
    _ticks = [v for v in np.arange(0, lim + 1, 100) if v <= lim]
    axE.set_xticks(_ticks); axE.set_yticks(_ticks)
    adjust_text(annsE, x=role.outgoing.values, y=role.incoming.values, ax=axE,
                expand=(1.45, 1.75), force_text=(0.6, 1.0), force_static=(0.5, 0.9),
                only_move={"text": "xy", "static": "xy", "explode": "xy"},
                arrowprops=dict(arrowstyle="-", color="0.55", lw=.6, shrinkA=1, shrinkB=3))
    axE.set_xlabel("Outgoing strength"); axE.set_ylabel("Incoming strength")
    axE.set_title("Signalling role of each cell type", loc="center", pad=10)

    # ---- F: expression of the ligands and receptors ----
    states = list(me.index); genes = list(me.columns)
    xs, ys, ss, cs = [], [], [], []
    # transposed relative to the wide version: eighteen genes as rows would give
    # 0.12 in per row at 180 mm, less than the 9 pt label height
    for i, gn in enumerate(genes):
        for j, st in enumerate(states):
            xs.append(i); ys.append(len(states) - 1 - j)
            ss.append(float(fe.loc[st, gn]) * 40); cs.append(float(me.loc[st, gn]))
    scF = axF.scatter(xs, ys, s=ss, c=cs, cmap="Reds", ec="0.45", lw=.25, vmin=0, vmax=4)
    # vertical gene labels: eighteen names at 42 degrees crowd a 2.6 in axis,
    # upright ones each stay inside their own 0.145 in column
    axF.set_xticks(range(len(genes)))
    axF.set_xticklabels(genes, rotation=90, ha="center", va="top", style="italic")
    axF.set_yticks(range(len(states)))
    axF.set_yticklabels([SHORT.get(c, c) for c in states][::-1])
    axF.set_xlim(-.8, len(genes) - .2); axF.set_ylim(-.8, len(states) - .2)
    # a centred title wider than its panel spills left onto the panel letter;
    # shortened so it fits inside the axes width
    axF.set_title("Expression of the ligands and receptors", loc="center", pad=10)
    # nudged off centre: the title is wider than the axes, and centred it sat
    # flush against the panel letter with no gap
    axF.title.set_x(0.565)
    # keys to the RIGHT of F rather than above it: the top strip pushed F's title
    # off its row and the right margin was otherwise unused
    for s_, lab in [(0.25 * 40, "25%"), (0.75 * 40, "75%")]:
        axF.scatter([], [], s=s_, c="0.75", ec="0.45", lw=.25, label=lab)
    axF.legend(frameon=False, fontsize=TIER[2], loc="upper left", bbox_to_anchor=(1.015, 1.02),
               ncol=1, handletextpad=.5, labelspacing=.9, title="cells\nexpressing",
               title_fontsize=TIER[2])
    cbF = axF.inset_axes([1.055, 0.02, 0.045, 0.36])
    cb = fig.colorbar(scF, cax=cbF, orientation="vertical", ticks=[0, 1, 2, 3, 4])
    # label parallel to the bar on its inner side, ticks on the outer side: the
    # two then sit on opposite flanks instead of competing for the same edge
    cb.set_label("mean expression", fontsize=TIER[2], rotation=90, labelpad=2)
    cb.ax.yaxis.set_ticks_position("right"); cb.ax.yaxis.set_label_position("right")
    cb.ax.tick_params(labelsize=TIER[2], length=2, pad=1.5)

    ax = dict(A=axA, B=axB, C=axC, D=axD, E=axE, F=axF, cbC=cbC, cbF=cbF, annsE=annsE)
    g["normalize"](fig)
    for a in (axC, axD, axE, axF):
        for t in a.get_yticklabels() + a.get_xticklabels():
            t.set_fontsize(TIER[2])
    return fig, ax


def close_top_gap(fig, ax, gap=0.030):
    """Pull the circular row down until it sits `gap` above the next row's title.

    The two circular panels need much less height than the rectangular ones, so
    a single gridspec hspace that clears the colour bar under row 2 leaves a
    band of white under row 1.  Shifting the row costs nothing: the white it
    leaves above is removed by the tight bounding box on save.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
    t2 = max(T.transform(ax[k].title.get_window_extent(r))[1][1] for k in "CD")
    y0 = min(ax[k].get_position().y0 for k in "AB")
    delta = y0 - (t2 + gap)
    if delta > 0:
        for k in "AB":
            p = ax[k].get_position()
            ax[k].set_position([p.x0, p.y0 - delta, p.width, p.height])
    fig.canvas.draw()
    return round(delta, 5)


def fit_circles(fig, ax, keys="AB", r=1.82):
    """Give each circular panel limits that match its frame, so it fills it.

    With aspect="equal" and the default box adjustment, matplotlib SHRINKS the
    axes to the data aspect: the drawn circle came out small and its left and
    right edges no longer lined up with the rectangular panels below.  Setting
    the limits from the frame's own aspect keeps aspect="equal" true while the
    panel keeps the rectangle the grid gave it.
    """
    fig.canvas.draw()
    for k in keys:
        ax[k].set_ylim(-r, r); ax[k].set_xlim(-r, r)   # datalim aspect widens x to the box
    fig.canvas.draw()
    return {k: [round(v, 3) for v in ax[k].get_xlim()] for k in keys}


def close_row_gaps(fig, ax, rows, gap=0.030):
    """Pull each row up until it sits `gap` above the next row's title.

    One gridspec hspace has to clear the deepest furniture in the figure, which
    leaves the other rows with a band of white beneath them.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()

    def low(rw):
        v = []
        for k in rw:
            a = ax[k]; ys = [T.transform(a.get_window_extent(r))[0][1]]
            for t in a.get_xticklabels():
                if t.get_text().strip():
                    ys.append(T.transform(t.get_window_extent(r))[0][1])
            if a.xaxis.label.get_text().strip():
                ys.append(T.transform(a.xaxis.label.get_window_extent(r))[0][1])
            # the circular panels' node labels are Text children that hang BELOW
            # the axes box; without them the gap is measured too small and a
            # label lands in the next row's title band
            for t in a.texts:
                if t.get_text().strip() and t.get_visible():
                    ys.append(T.transform(t.get_window_extent(r))[0][1])
            v.append(min(ys))
        return min(v)

    def top(rw):
        return max(T.transform(ax[k].title.get_window_extent(r))[1][1] for k in rw)

    moved = []
    for i in range(len(rows) - 1):
        d = low(rows[i]) - (top(rows[i + 1]) + gap)
        if d > 0:
            for k in rows[i]:
                p = ax[k].get_position()
                ax[k].set_position([p.x0, p.y0 - d, p.width, p.height])
            fig.canvas.draw(); r = fig.canvas.get_renderer()
        moved.append(round(d, 5))
    return moved


def decollide(fig, ax, anns, step=3.0, iters=60):
    """Push overlapping point labels apart along y.

    Fixed per-point offsets are not enough where two cell types sit almost on
    top of each other in the scatter (smooth muscle and epithelium differ by
    about 6 and 9 units on the two axes).
    """
    for _ in range(iters):
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        box = ax.get_window_extent(r)
        bb = [a.get_window_extent(r) for a in anns]
        moved = False
        # keep every label inside the axes: unbounded pushing sent the lowest
        # one onto the y tick labels
        for k, b in enumerate(bb):
            over = max(0.0, box.y0 - b.y0) - max(0.0, b.y1 - box.y1)
            if abs(over) > 0.5:
                x, y = anns[k].get_position()
                anns[k].set_position((x, y + over / fig.dpi * 72))
                moved = True
        if moved:
            fig.canvas.draw(); r = fig.canvas.get_renderer()
            bb = [a.get_window_extent(r) for a in anns]
        for i in range(len(anns)):
            for j in range(i + 1, len(anns)):
                if bb[i].overlaps(bb[j]):
                    moved = True
                    hi, lo = (i, j) if bb[i].y0 >= bb[j].y0 else (j, i)
                    for k, sgn in ((hi, 1), (lo, -1)):
                        x, y = anns[k].get_position()
                        anns[k].set_position((x, y + sgn * step))
        if not moved:
            break
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = [a.get_window_extent(r) for a in anns]
    return sum(1 for i in range(len(anns)) for j in range(i + 1, len(anns)) if bb[i].overlaps(bb[j]))


def place_letters(fig, ax, dy=0.012, dx=0.058):
    for a in fig.axes:
        for t in list(a.texts):
            if t.get_text() in list("ABCDEF") and t.get_fontweight() in ("bold", 700):
                t.remove()
    for t in list(fig.texts):
        if t.get_text() in list("ABCDEF") and t.get_fontweight() in ("bold", 700):
            t.remove()
    fig.canvas.draw()
    for k in "ABCDEF":
        p = ax[k].get_position()
        fig.text(p.x0 - dx, p.y1 + dy, k, fontsize=11, fontweight="bold", va="bottom", ha="left")
    fig.canvas.draw()


def audit(fig, mpl):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:26], b.get_text()[:26]) for i, (a, ba) in enumerate(tx)
          for b, bb in tx[i + 1:] if ba.overlaps(bb)]
    return dict(n_overlap=len(ov), overlaps=ov)


def level_titles(fig, ax, rows):
    # NOTE: this levels a row onto its HIGHEST title. That is only wanted when the
    # titles genuinely differ in height. When C and F carried key strips above their
    # axes their pads were 42 against 10 for D and E, so levelling lifted the whole
    # row and left every title floating a centimetre above its panel with the panel
    # letter stranded in between. The keys now sit to the right of those panels and
    # all pads are equal, so this is a no-op that only guards against drift.
    """Put both titles in a row on one baseline.

    A panel carrying a key strip above its axes needs a larger title pad than its
    neighbour, so equal pads leave the two titles at different heights.  The pad
    is recomputed from the axes top so every title in a row starts at the same
    figure height; setting title.set_position() does not work, because matplotlib
    recomputes the title position from the pad on every draw.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); T = fig.transFigure.inverted()
    fh = fig.get_size_inches()[1]
    out = {}
    for rw in rows:
        bot = {k: T.transform(ax[k].title.get_window_extent(r))[0][1] for k in rw}
        target = max(bot.values())
        for k in rw:
            ax[k].set_title(ax[k].get_title(), loc="center",
                            pad=(target - ax[k].get_position().y1) * fh * 72)
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        out[tuple(rw)] = [round(T.transform(ax[k].title.get_window_extent(r))[0][1], 5) for k in rw]
    return out

def title_block(fig, ax, keys, pads=None, gap_pt=5.0, size=11.0):
    """Set each panel's letter and title on one left-aligned line above its axes.

    A CENTRED title over a narrow panel spills past the axes on both sides: on the
    left it runs into the panel letter, which sits outside the axes box, and the
    letter then reads as detached from the words it belongs to. Left-aligning the
    title at the axes' left edge and putting the letter immediately before it, on
    the same baseline, fixes both. It also makes a row's titles level as soon as
    the axes frames are level, so no pad has to be reverse-computed from the
    tallest title -- the trick that previously lifted a whole row of titles a
    centimetre clear of their panels.

    `pads` is per key, in points: a panel whose own content reaches above its axes
    (a circle plot with node labels on its rim) needs more clearance than a
    rectangular one.
    """
    pads = pads or {}
    for a in fig.axes:
        for t in list(a.texts):
            if t.get_text() in keys and t.get_fontweight() in ("bold", 700):
                t.remove()
    for t in list(fig.texts):
        if t.get_text() in keys and t.get_fontweight() in ("bold", 700):
            t.remove()
    titles = {}
    for k in keys:
        # an Axes carries THREE independent title objects (left, centre, right).
        # set_title(loc="left") fills the left one and leaves the centre one with
        # its old text, so both render and each title overprints itself.
        # read whichever slot still holds text, so a second call after the rows
        # have been respaced does not blank every title
        txt = next((ax[k].get_title(loc=w) for w in ("center", "left", "right")
                    if ax[k].get_title(loc=w).strip()), "")
        ax[k].set_title("", loc="center")
        titles[k] = ax[k].set_title(txt, loc="left", pad=pads.get(k, 4.0))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    out = {}
    for k in keys:
        tb = inv.transform(titles[k].get_window_extent(r))
        # the letter's right edge sits gap_pt before the title's first glyph, and
        # both share a baseline, so the pair reads as one line
        x = tb[0][0] - gap_pt / 72 / fig.get_size_inches()[0]
        fig.text(x, tb[0][1], k, fontsize=size, fontweight="bold", va="baseline", ha="right")
        out[k] = (round(tb[0][0], 4), round(tb[0][1], 4))
    fig.canvas.draw()
    return out

def space_rows(fig, ax, rows, gap_mm=5.0, iters=3):
    """Give every pair of adjacent rows the same visual gap.

    Spacing computed from axes boxes is wrong whenever a panel's own content
    reaches outside its box: a circle plot's rim labels hang below it and rotated
    tick labels hang below a bubble plot, so equal box gaps produce very unequal
    visual ones -- here 6.4 mm above one row and 1.1 mm above the next, close
    enough to touching. This measures each row's lowest RENDERED element and the
    next row's title top, then shifts the lower rows until every gap matches.
    """
    def _low(ks):
        fig.canvas.draw()
        r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
        ys = []
        for k in ks:
            a = ax[k]
            ys.append(inv.transform(a.get_window_extent(r))[0][1])
            objs = list(a.texts) + list(a.get_xticklabels()) + list(a.get_yticklabels()) + [a.xaxis.label]
            lg = a.get_legend()
            if lg is not None:
                objs.append(lg)
            # a colour bar created with inset_axes is a CHILD axes, invisible to a
            # scan of the parent's own artists: leaving it out put Figure 5 panel
            # C's colour-bar label 3.7 px from the next row's title
            for c in a.child_axes:
                objs.append(c)
                objs.extend(list(c.get_xticklabels()) + list(c.get_yticklabels())
                            + [c.xaxis.label, c.yaxis.label])
            for t in objs:
                try:
                    if hasattr(t, "get_text") and not t.get_text().strip():
                        continue
                    ys.append(inv.transform(t.get_window_extent(r))[0][1])
                except Exception:
                    pass
        return min(ys)

    def _top(ks):
        fig.canvas.draw()
        r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
        return max(inv.transform(ax[k]._left_title.get_window_extent(r))[1][1] for k in ks)

    want = gap_mm / 25.4 / fig.get_size_inches()[1]
    for _ in range(iters):
        for i in range(len(rows) - 1):
            shift = (_low(rows[i]) - want) - _top(rows[i + 1])
            if abs(shift) < 1e-4:
                continue
            for k in rows[i + 1]:
                p = ax[k].get_position()
                ax[k].set_position([p.x0, p.y0 + shift, p.width, p.height])
            for j in range(i + 2, len(rows)):
                for k in rows[j]:
                    p = ax[k].get_position()
                    ax[k].set_position([p.x0, p.y0 + shift, p.width, p.height])
    fig.canvas.draw()
    return [round((_low(rows[i]) - _top(rows[i + 1])) * fig.get_size_inches()[1] * 25.4, 1)
            for i in range(len(rows) - 1)]
