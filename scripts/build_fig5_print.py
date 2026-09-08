"""Build Figure 5 (EAOC extension, 6 panels) at journal print width.

PLOT TYPES follow a caption census of the copy-number single-cell literature in
CNS-family journals, not a general single-cell one: 600 figure captions from 54
papers (2021-2026) in Nature (16), Nat Genet (11), Cell (6), Cancer Discov (6),
Cell Rep Med (7), Nat Cell Biol, Cancer Cell, Nat Med, Nat Cancer and Science,
retrieved from PubMed Central.  Of the 38 papers whose captions concern copy
number:

    heatmap of cells x genomic bins   66% of papers   <- the canonical panel
    UMAP / t-SNE embedding            58%
    phylogeny / clone tree            37%   (needs scDNA; not available here)
    box plot                          39%
    scatter                           34%
    bar chart                         34%
    density / histogram               32%
    violin                            24%

and among the 29 papers that treat copy-number burden as a per-cell quantity,
box plots and violins are used equally often (31% each).  Of the 23 papers with
a null-model or permutation panel, box plots (52%) and density/histograms (30%)
dominate, with the observed value marked against the null.

Panel order follows the same convention those papers use: the raw inferred
signal first, then the derived quantity, then the test against a null.

A  inferred copy-number heatmap, cells x genomic bins, grouped by stage
B  genome-wide mean inferred signal per stage
C  epithelial embedding coloured by copy-number burden
D  per-sample burden by stage, with the three malignant cohorts distinguished
E  transition index along the burden axis against the two-population mixing null
F  excess over the null, which is the panel that refutes the tipping point

The unit of replication in D is the sample, not the cell, for the reason the
paper as a whole argues.
"""

TIER = (11.0, 10.0, 9.0)
STAGES = ["0_normal", "1_benign_lesion", "2_malignant"]
STLAB = {"0_normal": "Normal endometrium", "1_benign_lesion": "Endometriotic lesion",
         "2_malignant": "Malignant (EAOC/OCCC)"}
STSHORT = {"0_normal": "Normal", "1_benign_lesion": "Lesion", "2_malignant": "Malignant"}
STC = {"0_normal": "#7BA7C7", "1_benign_lesion": "#E8A33D", "2_malignant": "#C0392B"}
COHC = {"GSE224334": "#C0392B", "GSE224333": "#8E44AD", "GSE291389": "#16A085"}
MG = "0.35"


def build(plt, np, pd, ctx):
    import matplotlib as mpl
    from matplotlib.colors import TwoSlopeNorm

    FR = ctx["FR"]; CH = ctx["CH"]; chrpos = ctx["chrpos"]
    um = ctx["um"]; ps = ctx["ps"]; obs_null = ctx["obs_null"]; mixnull = ctx["mixnull"]
    meta = ctx["meta"]; normalize = ctx["normalize"]

    # widened to the ~198 mm common width (see build_fig4_print.py)
    fig = plt.figure(figsize=(8.28, 8.9))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.40,
                          left=0.085, right=0.905, top=0.945, bottom=0.055)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, 0]); axF = fig.add_subplot(gs[2, 1])
    ax = dict(A=axA, B=axB, C=axC, D=axD, E=axE, F=axF)

    chr_order = sorted(chrpos, key=lambda c: chrpos[c])
    bounds = [chrpos[c] for c in chr_order] + [len(FR)]
    mids = [(bounds[i] + bounds[i + 1]) / 2 for i in range(len(chr_order))]
    keep = [i for i, c in enumerate(chr_order) if c.replace("chr", "") in
            ("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "21")]

    # ---- A: gain / loss frequency across the genome ----
    # replaces a cells-by-bins heatmap: the same information (which chromosomes
    # are altered, in what fraction of cells) in the mirrored-frequency form used
    # by 11% of the surveyed papers, and a heatmap already carries Figure 3B
    for k, st in enumerate(STAGES):
        off = -k * 210.0
        axA.fill_between(np.arange(len(FR)), off, off + FR[f"{st}_gain"].values,
                         color="#C0392B", lw=0)
        axA.fill_between(np.arange(len(FR)), off, off - FR[f"{st}_loss"].values,
                         color="#2E6DA4", lw=0)
        axA.axhline(off, color="0.45", lw=.5)
        # inside the track: outside the axes the three names ran into the tick labels
        axA.text(0.012, off + 95, STSHORT[st], transform=axA.get_yaxis_transform(),
                 ha="left", va="top", fontsize=TIER[2], color=STC[st])
    for bnd in bounds[1:-1]:
        axA.axvline(bnd - .5, color="0.88", lw=.35, zorder=0)
    axA.set_xlim(0, len(FR)); axA.set_ylim(-2 * 210 - 105, 105)
    axA.set_yticks([100, 50, 0, -50, -100]); axA.set_yticklabels(["100", "50", "0", "50", "100"])
    axA.set_ylabel("Cells altered (%)")
    axA.set_xticks([mids[i] for i in keep])
    axA.set_xticklabels([chr_order[i].replace("chr", "") for i in keep], fontsize=TIER[2])
    axA.set_xlabel("Chromosome")
    axA.spines["left"].set_bounds(-105, 105)
    axA.set_title("Gain and loss frequency", loc="center", pad=8)
    axA.text(0.995, 0.985, "gain", transform=axA.transAxes, ha="right", va="top",
             fontsize=TIER[2], color="#C0392B")
    axA.text(0.995, 0.905, "loss", transform=axA.transAxes, ha="right", va="top",
             fontsize=TIER[2], color="#2E6DA4")

    # ---- B: per-chromosome effect, lollipop ----
    # a second genome-axis line plot would repeat A; this asks the next question
    # (which chromosomes carry the difference) at the sample level.  Vertical
    # stems, because 22 chromosome labels down a 2 in axis collide at 9 pt
    x = np.arange(len(CH))
    axB.vlines(x, 0, CH.delta.values, color="0.72", lw=1.0, zorder=1)
    sig = CH.p.values < 0.05
    axB.scatter(x[sig], CH.delta.values[sig], s=24, c="#C0392B", ec="w", lw=.4, zorder=3,
                label="P < 0.05")
    axB.scatter(x[~sig], CH.delta.values[~sig], s=24, c="0.72", ec="w", lw=.4, zorder=3,
                label="not significant")
    axB.axhline(0, color="0.35", lw=.8)
    # upright labels: 22 two-digit numbers across a 2.4 in axis touch at 9 pt
    axB.set_xticks(x)
    axB.set_xticklabels(CH.chrom.values, fontsize=TIER[2], rotation=90,
                        ha="center", va="top")
    axB.set_xlim(-0.8, len(CH) - 0.2)
    axB.set_xlabel("Chromosome (ordered by effect)")
    axB.set_ylabel("Excess mean |signal|\n(malignant − rest, per sample)")
    axB.set_title("Chromosomes carrying the difference", loc="center", pad=8)
    axB.legend(frameon=False, fontsize=TIER[2], loc="upper right", handletextpad=.2,
               labelspacing=.25, borderaxespad=.3)

    # ---- C: embedding coloured by burden ----
    o = np.argsort(um.burden.values)
    scC = axC.scatter(um.x.values[o], um.y.values[o], c=um.burden.values[o], s=.6,
                      cmap="magma_r", vmin=0, vmax=float(np.percentile(um.burden, 99)),
                      lw=0, rasterized=True)
    axC.set_xticks([]); axC.set_yticks([])
    for sp_ in axC.spines.values():
        sp_.set_visible(False)
    axC.set_title("Burden on the embedding", loc="center", pad=8)
    axC.text(0.01, -0.02, "UMAP 1", transform=axC.transAxes, fontsize=TIER[2], color=MG,
             va="top")
    axC.text(-0.02, 0.01, "UMAP 2", transform=axC.transAxes, fontsize=TIER[2], color=MG,
             rotation=90)
    # horizontal, below the panel: inside it the bar covered cells, and the gutter
    # to its right belongs to D's axis title.  The embedding has no x axis, so the
    # strip under it is free
    cbC = axC.inset_axes([0.20, -0.115, 0.60, 0.035])
    cb = fig.colorbar(scC, cax=cbC, orientation="horizontal")
    cb.set_label("Copy-number burden", fontsize=TIER[2], labelpad=2)
    cb.ax.xaxis.set_ticks_position("bottom"); cb.ax.xaxis.set_label_position("bottom")
    cb.ax.tick_params(labelsize=TIER[2], length=2, pad=1.5)

    # ---- D: per-sample burden, violin ----
    # box plots already carry Figure 2A and Figure 3C; among the surveyed papers
    # that show copy-number burden as a distribution, violins and boxes are used
    # equally often (31% each)
    rng = np.random.default_rng(1); nsamp = {}
    for i, st in enumerate(STAGES):
        g = ps[ps.stage == st]; nsamp[st] = len(g)
        vp = axD.violinplot([g.burden.values], positions=[i], widths=.78,
                            showextrema=False, showmedians=False)
        for bdy in vp["bodies"]:
            bdy.set_facecolor(STC[st]); bdy.set_alpha(.30)
            bdy.set_edgecolor("0.45"); bdy.set_linewidth(.7)
        q1, med, q3 = np.percentile(g.burden.values, [25, 50, 75])
        axD.vlines(i, q1, q3, color="0.25", lw=3.2, zorder=4)
        axD.scatter([i], [med], s=16, c="w", ec="0.25", lw=.8, zorder=5)
        cols = [COHC.get(d, STC[st]) for d in g.dataset]
        axD.scatter(i + rng.uniform(-.16, .16, len(g)), g.burden.values, s=13, c=cols,
                    alpha=.85, ec="w", lw=.25, zorder=3)
    axD.set_xticks(range(3))
    axD.set_xticklabels([f"{STSHORT[s]}\n{nsamp[s]}" for s in STAGES], fontsize=TIER[2])
    for t, st in zip(axD.get_xticklabels(), STAGES):
        t.set_color(STC[st])
    axD.set_xlim(-.6, 2.6)
    # headroom so the key sits above every sample rather than over the top ones
    axD.set_ylim(top=float(ps.burden.max()) * 1.30)
    axD.set_ylabel("Burden (sample mean)")
    axD.set_title("Burden per sample", loc="center", pad=8)
    # the test result is in the caption
    axD.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=4.5, mfc=c, mec="w",
                                   mew=.3, label=k) for k, c in COHC.items()],
               frameon=False, fontsize=TIER[2], loc="upper left", bbox_to_anchor=(-0.005, 1.005),
               handletextpad=.3, labelspacing=.28, title="malignant cohort",
               title_fontsize=TIER[2], alignment="left")

    # ---- E: transition index against the mixing null ----
    axE.plot(mixnull.frac_mal * 100, mixnull.CI_mix_mean, color="0.45", lw=1.0,
             label="two-population mixing null")
    axE.fill_between(mixnull.frac_mal * 100, mixnull.CI_mix_mean - 2 * mixnull.CI_mix_sd,
                     mixnull.CI_mix_mean + 2 * mixnull.CI_mix_sd, color="0.75", alpha=.45,
                     lw=0, label="null ± 2 s.d.")
    axE.plot(obs_null.frac_mal * 100, obs_null.CI, "o-", ms=3.2, lw=1.0, color="#C0392B",
             label="observed")
    axE.set_xlabel("Malignant cells in bin (%)")
    axE.set_ylabel("Transition index")
    axE.set_title("Transition index against the null", loc="center", pad=8)
    # lower right, inside the panel: the curve rises on the left and the corner
    # under its descending limb is empty
    axE.legend(frameon=False, fontsize=TIER[2], loc="lower right", bbox_to_anchor=(1.01, -0.01),
               ncol=1, handlelength=1.3, labelspacing=.28)

    # ---- F: excess over the null ----
    ex = obs_null.excess_CI.values
    axF.axhline(0, color="0.35", lw=.8)
    # width from the smallest gap between bins: the bins are unevenly spaced and a
    # fixed width made neighbouring bars overlap near the low end
    wF = 0.8 * float(np.median(np.diff(np.sort(obs_null.frac_mal.values * 100))))
    axF.bar(obs_null.frac_mal * 100, ex, width=wF,
            color=["#C0392B" if e > 0 else "#7BA7C7" for e in ex],
            edgecolor="0.35", lw=.3)
    axF.set_xlabel("Malignant cells in bin (%)")
    axF.set_ylabel("Excess over null")
    axF.set_title("Excess over the null", loc="center", pad=8)
    # explicit in-range ticks: out-of-view tick labels stay as live text objects
    # that never render but still register in an overlap audit
    axF.set_xticks([0, 20, 40, 60, 80, 100])
    axF.set_xlim(-4, 100)
    axF.set_yticks([-0.0015, -0.0010, -0.0005, 0.0000])
    # the two numbers that were printed here are in the caption

    normalize(fig)
    for a in (axA, axB, axC, axD, axE, axF):
        a.title.set_fontsize(TIER[0])
    fig.canvas.draw()
    return fig, ax


def align_row_frames(fig, axes_list):
    fig.canvas.draw()
    y0 = max(a.get_position().y0 for a in axes_list)
    y1 = min(a.get_position().y1 for a in axes_list)
    for a in axes_list:
        p = a.get_position()
        a.set_position([p.x0, y0, p.width, y1 - y0])
    fig.canvas.draw()
    return dict(y0=[round(a.get_position().y0, 5) for a in axes_list])


def level_titles(fig, ax, rows):
    """Recompute each title's pad from the axes top so a row shares one baseline.

    title.set_position() does not survive: matplotlib recomputes the title
    position from the pad on every draw.
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


def place_letters(fig, ax, dy=0.028, dx=0.078):
    for a in ax.values():
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
