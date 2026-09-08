# `unit_of_replication.py` — is your cell-type composition result driven by donors or by cells?

A single-cell case–control study that compares cell-type composition has to choose a unit
of replication. Pooling cells across donors and testing the resulting contingency table
treats every cell as an independent observation. It is not: cells from one donor share
that donor's biology, dissociation, sampling depth and batch.

This module runs both tests on your data and calibrates them by permutation, so you can
see what each one does on **your** dataset before you report a result.

## What it computes

| | how the test works | what it ignores |
|---|---|---|
| `pooled` | 2×2 chi-square, cells of the type under test against all other cells, case against control | between-donor variance |
| `donor` | one proportion per donor, Mann–Whitney between arms | nothing structural; loses power at small donor counts |

`permutation_fpr` shuffles the group label **across donors**, never across cells — each
donor keeps its own cells and its own composition, and only the labels are randomised.
Under that null both tests should reject at α. In a survey of 6,111 tests spanning 523
public case–control settings, 124 diseases and 244 datasets, the donor-level test rejected
at a median rate of 0.042 and the pooled test at 0.836.

## Minimal use

```python
import pandas as pd, unit_of_replication as uor

# one row per cell (or a pre-aggregated table with an `n` column)
cells = pd.DataFrame({
    "donor_id":  [...],   # patient / sample identifier
    "group":     [...],   # two levels, e.g. "case" / "control"
    "cell_type": [...],   # your annotation
})

res = uor.diagnose(cells, case="case", n_perm=1000)
print(res)
print(uor.risk_note(res))
```

`diagnose` returns one row per cell type:

| column | meaning |
|---|---|
| `p_pooled`, `p_donor` | what each test says on your data |
| `fpr_pooled`, `fpr_donor` | how often each test rejects when the labels are shuffled across your donors |
| `upheld_at_donor_level` | for pooled-significant cell types, whether the donor-level test agrees |

`risk_note(res)` condenses that into one sentence for a methods paragraph or a reviewer reply.

## Reading the output

- **`fpr_pooled` far above α.** Expected — it was above 0.05 in 97% of the settings we
  surveyed and above 0.50 in 91%. The pooled p value on your data is not interpretable at
  face value.
- **`upheld_at_donor_level` is `False`.** The pooled test rejected and the donor-level test
  did not. Across the survey this happened for 69% of pooled-significant results.
- **`fpr_donor` above α.** Uncommon (18% of tests, and usually mildly). It occurs when a
  cell type is absent from several donors, which creates ties in the proportions.
- **A donor-level null with few donors is not evidence of absence.** Report the smallest
  attainable p value at your donor count alongside the observed one; with 3 versus 3 donors
  the two-sided minimum is 0.1.

## What it does not do

Composition only. Differential expression has the same pseudo-replication problem, but the
fix there is pseudobulk aggregation per donor, not this function. Nothing here corrects for
covariates; if age, sex or batch differ between arms, a donor-level regression on the
proportions is the right next step and this module is only the diagnosis.

## Requirements

`numpy`, `scipy`, `pandas`. The exact Mann–Whitney null is enumerated by dynamic
programming when both arms have at most 12 donors, and a tie-corrected normal
approximation is used above that.
