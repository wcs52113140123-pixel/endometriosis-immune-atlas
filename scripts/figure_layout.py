#!/usr/bin/env python
"""One layout pass for the multi-panel figures: letters, titles, row spacing.

The rule, arrived at by measurement rather than preference:

  * A panel letter and its title sit on ONE left-aligned line above the axes.
    A centred title over a narrow panel spills past the axes on both sides -
    into the letter on the left, towards the neighbour on the right - and reads
    as detached from the letter it belongs to.
  * A row's titles are level because the row shares one pad and its axes frames
    are aligned. No pad is reverse-computed from the tallest title; doing that
    once lifted a whole row of titles a centimetre clear of their panels.
  * Row spacing is measured from RENDERED INK, not from axes boxes. A circle
    plot's rim labels hang below its box, rotated tick labels hang below a
    bubble plot, and a colour bar made with inset_axes is a child axes invisible
    to a scan of its parent. Each of those, left out, produced a visibly uneven
    gap that box arithmetic called even.

LAYOUT records the arguments used for each published figure, so the pass is
reproducible rather than a set of numbers typed into a session.
"""

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


def move_aux(fig, ax, pairs, before):
    """Shift a colour bar by the same delta its parent panel moved."""
    for parent, aux in pairs:
        o = ax[aux]
        if not hasattr(o, "get_position"):
            continue
        p0 = before[parent]; p1 = ax[parent].get_position(); q = o.get_position()
        o.set_position([q.x0 + (p1.x0 - p0.x0), q.y0 + (p1.y0 - p0.y0), q.width, q.height])
    fig.canvas.draw()


LAYOUT = {
    # figure: rows, per-panel title pads in points, row gap in mm, colour bars to
    # carry with their parent. A pad above the default clears something the panel
    # itself draws above its axes: annotation strips on Figure 3 B, the stacked
    # tracks on Figure 5 A. Every panel in a row shares one pad or the row's
    # titles will not be level.
    "Figure2": dict(rows=[["A", "B"], ["C", "D"], ["E"]], pads={}, gap_mm=3.0, aux=[]),
    "Figure3": dict(rows=[["A", "B"], ["C", "D"], ["E", "F"], ["G", "H"]],
                    pads={"A": 30, "B": 30}, gap_mm=3.0, aux=[("B", "cbB")]),
    "Figure4": dict(rows=[["A", "B"], ["C", "D"], ["E", "F"]],
                    pads={"A": 20, "B": 20}, gap_mm=3.0, aux=[]),
    "Figure5": dict(rows=[["A", "B"], ["C", "D"], ["E", "F"]],
                    pads={"A": 14, "B": 14}, gap_mm=3.0, aux=[]),
}


def _post_figure3(fig, ax):
    """Re-park Figure 3's two colour bars after the rows move.

    cbB spanned panel B's full height, which put its lowest tick level with B's
    legend row; cbE sat in the band BELOW panel E, where respacing dropped it onto
    panel H's axis label. cbB is shortened at the bottom and cbE moved into the
    gutter to the right of E.
    """
    pB = ax["B"].get_position(); qb = ax["cbB"].get_position()
    ax["cbB"].set_position([qb.x0, pB.y0 + 0.06 * pB.height, qb.width, 0.94 * pB.height])
    pE = ax["E"].get_position()
    ax["cbE"].set_position([pE.x1 + 0.016, pE.y0 + 0.18 * pE.height, 0.013, 0.56 * pE.height])
    ax["cbE"].yaxis.set_ticks_position("right")
    ax["cbE"].yaxis.set_label_position("right")
    fig.canvas.draw()


POST = {"Figure3": _post_figure3}


def apply_layout(fig, ax, figure_name, align_row_frames):
    """Run the whole pass for one figure. Returns the measured row gaps."""
    cfg = LAYOUT[figure_name]
    keys = "".join(k for row in cfg["rows"] for k in row)
    for row in cfg["rows"]:
        if len(row) > 1:
            align_row_frames(fig, [ax[k] for k in row])
    title_block(fig, ax, keys, pads=cfg["pads"])
    before = {k: ax[k].get_position() for k in keys}
    gaps = space_rows(fig, ax, cfg["rows"], gap_mm=cfg["gap_mm"])
    if cfg["aux"]:
        move_aux(fig, ax, cfg["aux"], before)
    if figure_name in POST:
        POST[figure_name](fig, ax)
    title_block(fig, ax, keys, pads=cfg["pads"])
    return gaps
