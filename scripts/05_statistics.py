#!/usr/bin/env python
"""The statistical core of the manuscript, in four parts.

    python scripts/05_statistics.py --part all

    composition  per-sample immune composition; Treg case-control per cohort and
                 combined; the patient-level scan over all ten immune states
                 -> results/tables/per_sample_all_immune_states.csv
                    results/tables/patient_level_composition_scan.csv
    power        simulated power and minimum detectable effect for the
                 patient-level Treg comparison, per cohort and combined
                 -> results/tables/power_minimum_detectable_effect.csv
    methods      the same donor-by-cell-type counts under scCODA, and exported
                 for propeller (scripts/propeller.R)
                 -> results/tables/method_comparison_sccoda.csv
                    data/interim/donor_celltype_counts.csv
    phase        transcriptional menstrual-cycle phase, validated against the
                 samples with a reported phase, then used as a covariate
                 -> results/tables/cycle_phase_scores.csv

Inputs are the annotated objects written by 01_load_qc_integrate.py
(checkpoints/immune_annotated.h5ad, checkpoints/full_atlas.h5ad) and the sample
metadata in data/interim/gsm_full_meta.json. Every printed number is one that
appears in the manuscript.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
INTERIM = ROOT / "data" / "interim"
CKPT = ROOT / "checkpoints"

MIN_IMMUNE = 100          # a sample enters the analysis at this many immune cells
STATES = ["B", "CD4T", "CD8T", "Macrophage", "Mast", "Mono", "NK", "Plasma", "Treg", "pDC"]
SECRETORY = ["PAEP", "GPX3", "SPP1", "CXCL14", "DPP4", "MAOA", "SCGB2A1", "DEFB1", "GLUL", "MT1G"]
PROLIFERATIVE = ["MKI67", "TOP2A", "PCNA", "CCNB1", "SFRP1", "HMGB2", "TYMS", "RRM2", "MCM2", "CDK1"]


def sample_meta():
    return json.loads((INTERIM / "gsm_full_meta.json").read_text())


def condition_of(gsm, meta):
    d = meta.get(gsm, {})
    cond = (d.get("condition") or d.get("disease") or "").lower()
    tissue = (d.get("tissue") or "").lower()
    if "control" in cond or "control" in tissue or "normal" in tissue:
        return "Control"
    if "endometriosis" in cond or "endometriosis" in tissue:
        return "Endometriosis"
    return "NA"


def reported_phase(gsm, meta):
    d = meta.get(gsm, {})
    for k, v in (d.items() if isinstance(d, dict) else []):
        if any(x in k.lower() for x in ("phase", "cycle")):
            return str(v)
    return None


def _immune_obs():
    import scanpy as sc
    imm = sc.read_h5ad(CKPT / "immune_annotated.h5ad")
    meta = sample_meta()
    imm.obs["condition"] = [condition_of(g, meta) for g in imm.obs.gsm]
    return imm.obs[imm.obs.condition.isin(["Control", "Endometriosis"])]


def donor_counts(obs):
    """Donor x cell-type counts, the input every composition method here shares."""
    cnt = pd.crosstab(obs.gsm, obs.immune_type)
    meta = obs.groupby("gsm", observed=True).agg(dataset=("dataset", "first"),
                                                 condition=("condition", "first"))
    cnt = cnt.join(meta)
    return cnt[cnt[STATES].sum(axis=1) >= MIN_IMMUNE]


# ------------------------------------------------------------------ composition

def composition():
    obs = _immune_obs()
    prop = (pd.crosstab([obs.dataset, obs.gsm, obs.condition], obs.immune_type,
                        normalize="index") * 100).reset_index()
    prop.columns.name = None
    n = obs.groupby(["dataset", "gsm", "condition"], observed=True).size().rename("n_immune").reset_index()
    comp = prop.merge(n, on=["dataset", "gsm", "condition"])
    comp = comp[comp.n_immune >= MIN_IMMUNE]
    comp.to_csv(TABLES / "per_sample_all_immune_states.csv", index=False)

    cohorts = [c for c in comp.dataset.unique()
               if (comp.loc[comp.dataset == c, "condition"] == "Control").sum() >= 3]
    rows = []
    for state in STATES:
        zs, per = [], {}
        for c in cohorts:
            d = comp[comp.dataset == c]
            case = d.loc[d.condition == "Endometriosis", state].values
            ctrl = d.loc[d.condition == "Control", state].values
            p = mannwhitneyu(case, ctrl, alternative="two-sided").pvalue
            per[c] = p
            sign = np.sign(np.median(case) - np.median(ctrl)) or 1.0
            zs.append(sign * norm.isf(p / 2))
        z = np.sum(zs) / np.sqrt(len(zs))
        rows.append({"immune_state": state,
                     **{f"P_{c}": round(per[c], 4) for c in cohorts},
                     "z_combined": round(z, 3), "P_combined": 2 * norm.sf(abs(z))})
    scan = pd.DataFrame(rows)
    scan["q_combined_BH"] = multipletests(scan.P_combined, method="fdr_bh")[1]
    scan = scan.sort_values("P_combined").round({"P_combined": 4, "q_combined_BH": 3})
    scan.to_csv(TABLES / "patient_level_composition_scan.csv", index=False)
    # these are the cohorts carrying the full ten-state annotation AND at least
    # three controls; GSE213216 has four controls but no ten-state breakdown
    print(f"cohorts scanned across all ten states: {cohorts}")
    print(scan.to_string(index=False))
    print(f"states significant at q < 0.05: {int((scan.q_combined_BH < 0.05).sum())} of {len(scan)}")


# ------------------------------------------------------------------------ power

def _sim_power(ctrl, case, delta, n_ctrl=None, n_case=None, n_sim=20000, alpha=0.05, seed=0):
    """Power of the one-sided donor-level Mann-Whitney to detect `delta`.

    n_sim is large because the reported values are quoted in the manuscript and
    must not move between runs; at 20,000 draws the Monte Carlo error is under
    0.005.

    Both arms are drawn from a common null built by centring the observed case
    values on the control mean, so the simulation inherits the real
    between-sample dispersion rather than assuming a distribution.
    """
    rng = np.random.default_rng(seed)
    n_ctrl = n_ctrl or len(ctrl)
    n_case = n_case or len(case)
    base = np.concatenate([ctrl, case - np.mean(case) + np.mean(ctrl)])
    hits = 0
    for _ in range(n_sim):
        c = rng.choice(base, n_ctrl, replace=True)
        e = rng.choice(base, n_case, replace=True) + delta
        hits += mannwhitneyu(e, c, alternative="greater").pvalue < alpha
    return hits / n_sim


def _sim_power_combined(packs, n_sim=8000, alpha=0.05, seed=1):
    """Power of the dataset-stratified model that is the manuscript's primary test."""
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_sim):
        frames = []
        for name, base, n_ctrl, n_case, delta in packs:
            c = rng.choice(base, n_ctrl, replace=True)
            e = rng.choice(base, n_case, replace=True) + delta
            frames.append(pd.DataFrame({
                "y": np.log(np.clip(np.concatenate([c, e]), 0, None) + 0.1),
                "cond": ["Control"] * n_ctrl + ["Endometriosis"] * n_case,
                "ds": name}))
        sim = pd.concat(frames)
        m = smf.ols("y ~ C(cond) + C(ds)", data=sim).fit()
        hits += (m.pvalues.get("C(cond)[T.Endometriosis]", 1) < alpha and
                 m.params.get("C(cond)[T.Endometriosis]", 0) > 0)
    return hits / n_sim


def power():
    allps = pd.read_csv(TABLES / "per_sample_treg_all3.csv")
    rows, packs = [], []
    for ds, d in allps.groupby("dataset"):
        ctrl = d.loc[d.cond == "Control", "Treg_pct"].values
        case = d.loc[d.cond == "Endometriosis", "Treg_pct"].values
        # the effect the pooled-cell analysis implies, i.e. cell-weighted means
        pooled_c = np.average(ctrl, weights=d.loc[d.cond == "Control", "n_immune"])
        pooled_e = np.average(case, weights=d.loc[d.cond == "Endometriosis", "n_immune"])
        delta = pooled_e - pooled_c
        sd = np.sqrt(((len(ctrl) - 1) * np.var(ctrl, ddof=1) +
                      (len(case) - 1) * np.var(case, ddof=1)) / (len(ctrl) + len(case) - 2))
        mde_d = TTestIndPower().solve_power(effect_size=None, nobs1=len(ctrl), alpha=0.05,
                                            power=0.8, ratio=len(case) / len(ctrl),
                                            alternative="larger")
        need = ">60"
        for n_ctrl in range(4, 61):
            if _sim_power(ctrl, case, delta, n_ctrl=n_ctrl, n_sim=20000) >= 0.8:
                need = str(n_ctrl)
                break
        rows.append({"cohort": ds, "n_ctrl": len(ctrl), "n_case": len(case),
                     "pooled_ctrl_pct": round(pooled_c, 2), "pooled_case_pct": round(pooled_e, 2),
                     "pooled_implied_delta": round(delta, 2), "between_sample_SD": round(sd, 2),
                     "MDE_80pct_power": round(mde_d * sd, 2),
                     "power_for_pooled_delta": round(_sim_power(ctrl, case, delta), 3),
                     "controls_needed_80pct": need})
        base = np.concatenate([ctrl, case - np.mean(case) + np.mean(ctrl)])
        packs.append((ds, base, len(ctrl), len(case), delta))
    out = pd.DataFrame(rows)
    out["power_combined_81samples"] = round(_sim_power_combined(packs), 3)
    out.to_csv(TABLES / "power_minimum_detectable_effect.csv", index=False)
    print(out.to_string(index=False))


# ---------------------------------------------------------------------- methods

def methods():
    cnt = donor_counts(_immune_obs())
    cnt.reset_index().to_csv(INTERIM / "donor_celltype_counts.csv", index=False)
    print(f"donor x cell-type counts: {cnt.shape}; run scripts/propeller.R for the propeller half")

    from sccoda.util import cell_composition_data as scdat
    from sccoda.util import comp_ana as mod
    out = []
    for ds in ["GSE179640", "GSE214411"]:
        d = cnt[cnt.dataset == ds]
        df = d[STATES].reset_index(drop=True)
        df["condition"] = d.condition.values
        data = scdat.from_pandas(df, covariate_columns=["condition"])
        model = mod.CompositionalAnalysis(data, formula="condition",
                                          reference_cell_type="automatic")
        res = model.sample_hmc(num_results=8000, num_burnin=2000, verbose=False)
        eff = res.effect_df.reset_index()
        eff["dataset"] = ds
        out.append(eff)
    sc_res = pd.concat(out)
    cols = [c for c in ["dataset", "Cell Type", "Final Parameter",
                        "Inclusion probability", "log2-fold change"] if c in sc_res.columns]
    sc_res[cols].to_csv(TABLES / "method_comparison_sccoda.csv", index=False)
    print(sc_res.loc[sc_res["Cell Type"] == "Treg", cols].to_string(index=False))
    credible = sc_res[sc_res["Inclusion probability"] > 0.8]
    print(f"cell types above the 0.8 inclusion threshold: {len(credible)}")


# ------------------------------------------------------------------------ phase

def phase():
    import scanpy as sc
    meta = sample_meta()
    allps = pd.read_csv(TABLES / "per_sample_treg_all3.csv")
    wanted = set(allps.gsm)

    adata = sc.read_h5ad(CKPT / "full_atlas.h5ad", backed="r")
    secretory = [g for g in SECRETORY if g in adata.var_names]
    proliferative = [g for g in PROLIFERATIVE if g in adata.var_names]
    keep = (adata.obs.celltype.isin(["Epithelial", "Stromal_Fib"]).values &
            adata.obs.gsm.isin(wanted).values)
    sub = adata[keep, secretory + proliferative].to_memory()
    sub.X = sub.X.astype("float32")
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    expr = pd.DataFrame(sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X,
                        columns=secretory + proliferative)
    expr["gsm"] = sub.obs.gsm.values
    expr["dataset"] = sub.obs.dataset.values
    pb = expr.groupby(["dataset", "gsm"], observed=True).mean(numeric_only=True)
    z = (pb - pb.mean()) / pb.std().replace(0, 1)
    pb["phase_score"] = z[secretory].mean(axis=1) - z[proliferative].mean(axis=1)
    pb = pb.reset_index()
    pb["reported_phase"] = [reported_phase(g, meta) for g in pb.gsm]
    pb["condition"] = [condition_of(g, meta) for g in pb.gsm]

    lab = pb[pb.reported_phase.notna()]
    s = lab.loc[lab.reported_phase == "secretory", "phase_score"].values
    p = lab.loc[lab.reported_phase == "proliferative", "phase_score"].values
    u = mannwhitneyu(s, p, alternative="greater")
    print(f"estimator validation on {len(lab)} labelled samples: "
          f"AUC = {u.statistic / (len(s) * len(p)):.2f}, one-sided P = {u.pvalue:.3f}")

    cc = pb[pb.condition.isin(["Control", "Endometriosis"])].copy()
    for ds, d in cc.groupby("dataset"):
        case = d.loc[d.condition == "Endometriosis", "phase_score"]
        ctrl = d.loc[d.condition == "Control", "phase_score"]
        if len(ctrl) >= 3:
            print(f"  {ds}: inferred phase case vs control P = "
                  f"{mannwhitneyu(case, ctrl).pvalue:.3f}")
    cc["Treg_pct"] = cc.gsm.map(allps.set_index("gsm").Treg_pct)
    cc = cc.dropna(subset=["Treg_pct"])
    for label, formula in [("without phase", "np.log(Treg_pct+0.1) ~ C(condition) + C(dataset)"),
                           ("with phase score", "np.log(Treg_pct+0.1) ~ C(condition) + C(dataset) + phase_score")]:
        m = smf.ols(formula, data=cc).fit()
        beta = m.params.get("C(condition)[T.Endometriosis]")
        ci = m.conf_int().loc["C(condition)[T.Endometriosis]"]
        print(f"  {label:17s} n={int(m.nobs)} beta={beta:+.3f} "
              f"95%CI [{ci[0]:.2f},{ci[1]:.2f}] P={m.pvalues.get('C(condition)[T.Endometriosis]'):.3f}")
    pb.to_csv(TABLES / "cycle_phase_scores.csv", index=False)


# ---------------------------------------------------------------- primary test

def primary():
    """Tabulate the manuscript's headline combined test.

    The dataset-stratified model over every sample, and the Stouffer combination
    of the per-cohort tests, are the two numbers the paper's central claim rests
    on. They were previously only printed, which left the most important figures
    in the manuscript with no table a reader could check them against.

    The sensitivity fit adds the three GSE214411 samples that carry no condition
    field, treated as controls.
    """
    rows = []
    for label, src in [("primary", "per_sample_treg_all3.csv"),
                       ("sensitivity: unstated-condition samples as controls",
                        "per_sample_treg_sensitivity.csv")]:
        d = pd.read_csv(TABLES / src)
        y = "lt" if "lt" in d.columns else "y"
        m = smf.ols(f"{y} ~ C(cond) + C(dataset)", data=d).fit()
        term = "C(cond)[T.Endometriosis]"
        ci = m.conf_int().loc[term]
        rows.append({"analysis": label, "n_samples": int(m.nobs),
                     "n_control": int((d.cond == "Control").sum()),
                     "n_endometriosis": int((d.cond == "Endometriosis").sum()),
                     "beta_log_scale": round(m.params[term], 3),
                     "CI95_low": round(ci[0], 2), "CI95_high": round(ci[1], 2),
                     "P_two_sided": round(m.pvalues[term], 2)})

    per = pd.read_csv(TABLES / "treg_samplelevel_tests.csv")
    z = np.sum([norm.isf(p) for p in per.p_onesided]) / np.sqrt(len(per))
    rows.append({"analysis": f"Stouffer combination of {len(per)} per-cohort one-sided tests",
                 "n_samples": int(per.n_ctrl.sum() + per.n_endo.sum()),
                 "n_control": int(per.n_ctrl.sum()), "n_endometriosis": int(per.n_endo.sum()),
                 "beta_log_scale": np.nan, "CI95_low": np.nan, "CI95_high": np.nan,
                 "P_two_sided": round(norm.sf(z), 2)})

    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "treg_primary_tests.csv", index=False)
    print(out.to_string(index=False))


# ------------------------------------------------------- atlas scale and burden

def atlas():
    """Tabulate the atlas size, its lineage composition and its cohort scale.

    The headline atlas count, the number of lineages and the three small-lineage
    counts quoted in the Figure 1 legend had no table behind them, so none of
    them could be checked against anything.
    """
    import scanpy as sc
    a = sc.read_h5ad(CKPT / "full_atlas.h5ad")
    comp = a.obs.celltype.value_counts().rename_axis("lineage").reset_index(name="n_cells")
    comp["pct"] = (comp.n_cells / a.n_obs * 100).round(2)
    comp.to_csv(TABLES / "atlas_composition.csv", index=False)

    per = a.obs.groupby(["dataset", "gsm"], observed=True).size().rename("n_cells").reset_index()
    per.to_csv(TABLES / "atlas_per_sample.csv", index=False)
    coh = (per.groupby("dataset", observed=True)
           .agg(n_samples=("gsm", "nunique"), n_cells=("n_cells", "sum")).reset_index())
    coh.to_csv(TABLES / "atlas_by_cohort.csv", index=False)
    print(comp.to_string(index=False))
    print(f"\n{comp.n_cells.sum():,} cells, {len(comp)} lineages, "
          f"{per.gsm.nunique()} samples, {len(coh)} cohorts")


def burden():
    """Copy-number burden by stage in BOTH units, per sample and per cell.

    The manuscript quoted a cell-level median beside a per-sample figure. The
    per-sample column is the one the text uses, because the sample is the unit
    of replication throughout; the per-cell column is kept so the earlier
    figures can still be traced.
    """
    import scanpy as sc
    e = sc.read_h5ad(CKPT / "comb_epi_final.h5ad")
    o = e.obs.copy()
    o["cnv_burden"] = pd.to_numeric(o.cnv_burden, errors="coerce")
    per_cell = o.groupby("stage", observed=True).cnv_burden.agg(
        n_cells="size", median_per_cell="median").round(4)
    samp = (o.groupby(["stage", "gsm"], observed=True)
            .agg(burden=("cnv_burden", "mean"), n=("cnv_burden", "size")).reset_index())
    samp = samp[samp.n >= 50]
    per_samp = samp.groupby("stage", observed=True).burden.agg(
        n="size", median_per_sample="median").round(4)
    out = per_samp.join(per_cell).reset_index()
    out.to_csv(TABLES / "cnv_burden_by_stage.csv", index=False)
    print(out.to_string(index=False))


PARTS = {"composition": composition, "primary": primary, "power": power,
         "methods": methods, "phase": phase, "atlas": atlas, "burden": burden}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--part", default="all", help="composition, power, methods, phase, or all")
    args = ap.parse_args()
    parts = list(PARTS) if args.part == "all" else args.part.split(",")
    unknown = [p for p in parts if p not in PARTS]
    if unknown:
        raise SystemExit(f"unknown part(s): {', '.join(unknown)}")
    for p in parts:
        print(f"\n=== {p} ===")
        PARTS[p]()


if __name__ == "__main__":
    main()
