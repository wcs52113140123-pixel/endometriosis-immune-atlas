"""Pooled-cell versus donor-level testing of cell-type composition.

Given a donor x cell-type count matrix and a group label per donor, this runs
the two tests that a single-cell case-control study can choose between:

  pooled  - cells are treated as independent observations.  Cells of the type
            under test are counted against all other cells in each arm and the
            2x2 table is tested with a chi-square test.  This is what pooling
            cells across donors does, and it ignores between-donor variance.
  donor   - one proportion per donor, compared between arms with a two-sided
            Mann-Whitney U test.  The donor is the unit of replication.

`permutation_fpr` shuffles the group label ACROSS DONORS (never across cells,
which would destroy the donor structure the pooled test ignores) and reports
how often each test rejects at alpha.  Under the null both should reject at
alpha; the pooled test does not.

The U distribution is enumerated exactly by dynamic programming when the arms
are small, which is the regime where a normal approximation is unreliable and
also the regime most public case-control studies sit in.
"""

from functools import lru_cache

import numpy as np
from scipy.stats import chi2


# ---------------------------------------------------------------- exact MWU

@lru_cache(maxsize=4096)
def _u_null(n1: int, n2: int):
    """Exact null distribution of U over 0 .. n1*n2.

    Uses the standard recurrence for the Gaussian binomial coefficient,
    f(i, j, u) = f(i-1, j, u-j) + f(i, j-1, u), which counts the arrangements
    of i first-sample and j second-sample items giving rank-sum displacement u.
    """
    umax = n1 * n2
    prev = [np.zeros(umax + 1) for _ in range(n2 + 1)]
    for j in range(n2 + 1):
        prev[j][0] = 1.0                       # i = 0: only U = 0
    for i in range(1, n1 + 1):
        cur = [np.zeros(umax + 1) for _ in range(n2 + 1)]
        cur[0][0] = 1.0                        # j = 0: only U = 0
        for j in range(1, n2 + 1):
            shifted = np.zeros(umax + 1)
            if j <= umax:
                shifted[j:] = prev[j][: umax + 1 - j]
            cur[j] = shifted + cur[j - 1]
        prev = cur
    dp = prev[n2]
    return dp / dp.sum()


def _mwu_p_exact(u, n1, n2):
    """Two-sided p from the exact null, vectorised over an array of U."""
    pmf = _u_null(n1, n2)
    cdf = np.cumsum(pmf)
    sf = 1.0 - np.concatenate([[0.0], cdf[:-1]])          # P(U >= u)
    u = np.clip(np.round(u).astype(int), 0, n1 * n2)
    lower = cdf[u]
    upper = sf[u]
    return np.minimum(1.0, 2.0 * np.minimum(lower, upper))


def _mwu_p_normal(u, n1, n2, tie_term=0.0):
    mu = n1 * n2 / 2.0
    n = n1 + n2
    sd = np.sqrt(n1 * n2 * (n + 1 - tie_term) / 12.0)
    if sd == 0:
        return np.ones_like(np.asarray(u, dtype=float))
    from scipy.stats import norm
    z = (np.abs(u - mu) - 0.5) / sd
    return 2.0 * norm.sf(z)


EXACT_MAX = 12          # exact null when both arms are at most this size


def mwu_p(u, n1, n2, tie_term=0.0):
    if max(n1, n2) <= EXACT_MAX:
        return _mwu_p_exact(u, n1, n2)
    return _mwu_p_normal(u, n1, n2, tie_term)


# ------------------------------------------------------------- the two tests

def donor_u(prop, is_case):
    """U statistic of the case arm for every column of `prop` (donors x types)."""
    order = np.argsort(prop, axis=0, kind="mergesort")
    ranks = np.empty_like(prop)
    d = prop.shape[0]
    idx = np.arange(1, d + 1)[:, None]
    np.put_along_axis(ranks, order, np.broadcast_to(idx, prop.shape).astype(float), axis=0)
    # average ranks within ties, column by column (ties are rare in proportions
    # but common when a cell type is absent from several donors)
    for k in range(prop.shape[1]):
        col = prop[:, k]
        uniq, inv, cnts = np.unique(col, return_inverse=True, return_counts=True)
        if (cnts > 1).any():
            sums = np.zeros(len(uniq))
            np.add.at(sums, inv, ranks[:, k])
            ranks[:, k] = (sums / cnts)[inv]
    n1 = int(is_case.sum())
    n2 = prop.shape[0] - n1
    r1 = is_case.astype(float) @ ranks
    return r1 - n1 * (n1 + 1) / 2.0, n1, n2


def donor_p(prop, is_case):
    u, n1, n2 = donor_u(prop, is_case)
    return mwu_p(u, n1, n2), n1, n2


def pooled_p(counts, is_case):
    """Chi-square on the 2x2 table [type / other] x [case / control] per type."""
    tot = counts.sum(axis=1)
    a = is_case.astype(float) @ counts                 # case cells of each type
    b = counts.sum(axis=0) - a                         # control cells of each type
    ca = float(is_case @ tot)                          # all case cells
    cb = float(tot.sum() - ca)
    c = ca - a
    d = cb - b
    n = ca + cb
    with np.errstate(divide="ignore", invalid="ignore"):
        num = n * (a * d - b * c) ** 2
        den = (a + b) * (c + d) * (a + c) * (b + d)
        stat = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    return chi2.sf(stat, 1)


# ------------------------------------------------------------- calibration

def permutation_fpr(counts, is_case, n_perm=1000, alpha=0.05, seed=0):
    """Reject rate of both tests when the group label is shuffled across donors."""
    rng = np.random.default_rng(seed)
    d = counts.shape[0]
    n1 = int(is_case.sum())
    prop = counts / counts.sum(axis=1, keepdims=True)
    tot = counts.sum(axis=1)
    grand = counts.sum(axis=0)
    n_cells = float(tot.sum())

    # ranks do not depend on the labels, so they are computed once
    _, _, _ = donor_u(prop, is_case)
    order = np.argsort(prop, axis=0, kind="mergesort")
    ranks = np.empty_like(prop)
    idx = np.arange(1, d + 1)[:, None]
    np.put_along_axis(ranks, order, np.broadcast_to(idx, prop.shape).astype(float), axis=0)
    for k in range(prop.shape[1]):
        col = prop[:, k]
        uniq, inv, cnts = np.unique(col, return_inverse=True, return_counts=True)
        if (cnts > 1).any():
            sums = np.zeros(len(uniq))
            np.add.at(sums, inv, ranks[:, k])
            ranks[:, k] = (sums / cnts)[inv]

    Pi = np.zeros((n_perm, d))
    for b in range(n_perm):
        Pi[b, rng.choice(d, n1, replace=False)] = 1.0

    r1 = Pi @ ranks
    u = r1 - n1 * (n1 + 1) / 2.0
    p_donor = mwu_p(u, n1, d - n1)          # vectorised over the whole B x K block

    A = Pi @ counts
    CA = (Pi @ tot)[:, None]
    B = grand[None, :] - A
    C = CA - A
    Dm = (n_cells - CA) - B
    with np.errstate(divide="ignore", invalid="ignore"):
        num = n_cells * (A * Dm - B * C) ** 2
        den = (A + B) * (C + Dm) * (A + C) * (B + Dm)
        stat = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    p_pool = chi2.sf(stat, 1)

    return (p_pool < alpha).mean(axis=0), (p_donor < alpha).mean(axis=0)


# ------------------------------------------------------------------ frontend

def diagnose(df, donor="donor_id", group="group", cell_type="cell_type",
             case=None, n_perm=1000, alpha=0.05, min_cells=50, min_donors=3, seed=0):
    """Compare the two tests on one case-control single-cell dataset.

    `df` is one row per cell (or a pre-aggregated table with a `n` column) with
    a donor column, a two-level group column and a cell-type column.  Returns
    one row per cell type: the p value each test gives, the permutation false
    positive rate of each test on this dataset, and whether a pooled-significant
    result is upheld when the donor is the unit.

    The permutation shuffles the group label across donors, so it preserves each
    donor's own cell counts and composition and asks only what the tests do when
    the labels carry no information.
    """
    import pandas as pd

    d = df.copy()
    if "n" not in d.columns:
        d = d.groupby([donor, group, cell_type], observed=True).size().rename("n").reset_index()
    piv = d.pivot_table(index=[donor, group], columns=cell_type, values="n",
                        aggfunc="sum", observed=True).fillna(0)
    groups = [g for _, g in piv.index]
    levels = sorted(set(groups))
    if len(levels) != 2:
        raise ValueError(f"group column must have exactly two levels, found {levels}")
    case = case if case is not None else levels[1]
    is_case = np.array([g == case for g in groups])
    if is_case.sum() < min_donors or (~is_case).sum() < min_donors:
        raise ValueError(f"need at least {min_donors} donors per arm, "
                         f"have {int(is_case.sum())} and {int((~is_case).sum())}")

    N = piv.values.astype(float)
    keep = (N.sum(0) >= min_cells) & ((N > 0).sum(0) >= min_donors)
    N = N[:, keep]
    types = list(np.asarray(piv.columns)[keep])
    if not types:
        raise ValueError("no cell type passes min_cells / min_donors")

    prop = N / N.sum(axis=1, keepdims=True)
    p_pool = pooled_p(N, is_case)
    p_don, n1, n2 = donor_p(prop, is_case)
    fpr_pool, fpr_don = permutation_fpr(N, is_case, n_perm=n_perm, alpha=alpha, seed=seed)

    out = pd.DataFrame({
        "cell_type": types,
        "n_cells": N.sum(0).astype(int),
        "n_case_donors": n1, "n_control_donors": n2,
        "prop_case": prop[is_case].mean(0), "prop_control": prop[~is_case].mean(0),
        "p_pooled": p_pool, "p_donor": p_don,
        "fpr_pooled": fpr_pool, "fpr_donor": fpr_don,
    })
    out["upheld_at_donor_level"] = np.where(
        out.p_pooled < alpha, out.p_donor < alpha, pd.NA)
    return out.sort_values("p_donor").reset_index(drop=True)


def risk_note(res, alpha=0.05):
    """One-line summary of what the diagnosis found."""
    sig = res.p_pooled < alpha
    upheld = int((sig & (res.p_donor < alpha)).sum())
    return (f"pooled-cell test rejects {int(sig.sum())}/{len(res)} cell types; "
            f"{upheld} of those hold when the donor is the unit; "
            f"median permutation false-positive rate {res.fpr_pooled.median():.2f} (pooled) "
            f"vs {res.fpr_donor.median():.2f} (donor-level) against a nominal {alpha}")
