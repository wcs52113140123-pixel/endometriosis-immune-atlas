#!/usr/bin/env python
"""Collect Supplementary Tables S1-S8 into one workbook for submission.

    python scripts/06_supplementary_workbook.py

Reads the per-analysis CSVs in results/tables/ and writes
supplementary/Supplementary_Tables_S1_S8.xlsx, one sheet per supplementary
table, in the order the manuscript cites them. A missing input is reported
rather than silently skipped, because a gap here means a table the manuscript
promises does not exist.
"""

from pathlib import Path

import json

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "results" / "tables"
OUT = ROOT / "supplementary" / "Supplementary_Tables_S1_S8.xlsx"

STATES = ["B", "CD4T", "CD8T", "Macrophage", "Mast", "Mono", "NK", "Plasma", "Treg", "pDC"]


def read(name):
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run the analysis that writes it first")
    return pd.read_csv(path)


def build():
    sheets = {}

    # S1  every sample that enters the case-control analysis, with the
    #     proportions and programme scores used anywhere in the manuscript
    s1 = read("per_sample_treg_all3.csv").rename(columns={"lt": "log_Treg_pct"})
    # Mac_pct, NK_pct and CD8_pct duplicate pct_Macrophage, pct_NK and pct_CD8T
    # exactly (verified equal); only the programme scores are kept from this file
    prog = read("per_sample_immune_stats.csv")[["gsm", "M2_in_mac", "CD8_exh", "NK_cyt"]]
    states = read("per_sample_all_immune_states.csv").drop(
        columns=["dataset", "condition", "n_immune"], errors="ignore")
    states = states.rename(columns={c: f"pct_{c}" for c in STATES})
    sheets["S1_per_sample_composition"] = s1.merge(prog, on="gsm", how="left").merge(
        states, on="gsm", how="left")

    sheets["S2_spatial_paired_segments"] = read("dsp_patientlevel_paired.csv")
    sheets["S3_scIB_benchmark"] = read("scib_benchmark.csv")

    # S4  the four sensitivity analyses of the spatial finding, stacked
    s4 = []
    for name, label in [("sens_leave_one_patient_out.csv", "leave-one-patient-out"),
                        ("sens_leave_one_gene_out.csv", "leave-one-gene-out"),
                        ("sens_per_gene.csv", "single gene"),
                        ("sens_geneset_comparison.csv", "external gene sets")]:
        block = read(name)
        block.insert(0, "analysis", label)
        s4.append(block)
    sheets["S4_spatial_sensitivity"] = pd.concat(s4, ignore_index=True)

    sheets["S5_composition_scan"] = read("patient_level_composition_scan.csv")
    sheets["S6_power_MDE"] = read("power_minimum_detectable_effect.csv")

    prop = read("method_comparison_propeller.csv"); prop.insert(0, "tool", "propeller")
    coda = read("method_comparison_sccoda.csv"); coda.insert(0, "tool", "scCODA")
    sheets["S7_method_comparison"] = pd.concat([prop, coda], ignore_index=True)

    phase = read("cycle_phase_scores.csv")
    sheets["S8_cycle_phase"] = phase[[c for c in
        ["dataset", "gsm", "condition", "reported_phase", "phase_score"] if c in phase.columns]]

    OUT.parent.mkdir(exist_ok=True)
    # The journal requires a title on the same page as each supplementary table,
    # so every sheet carries its title in A1 and the data starts two rows below.
    # Titles are read from the manuscript's own supplementary legends at build
    # time; the previous hardcoded copy went stale when a legend was edited and
    # the audit caught the workbook and the paper disagreeing.
    TITLES = json.load(open(ROOT / "data" / "interim" / "supp_legends.json"))["tab"]
    from openpyxl.styles import Font
    with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
        for name, frame in sheets.items():
            sheet = name[:31]
            key = name.split("_")[0]
            if key not in TITLES:
                raise SystemExit(f"no title for sheet {sheet}; add it to TITLES")
            frame.to_excel(writer, sheet_name=sheet, index=False, startrow=2)
            ws = writer.sheets[sheet]
            ws["A1"] = f"Supplementary Table {key}. {TITLES[key]}"
            ws["A1"].font = Font(bold=True)
    for name, frame in sheets.items():
        print(f"  {name:34s} {frame.shape[0]:4d} rows x {frame.shape[1]:2d} cols")
    print(f"written: {OUT} ({OUT.stat().st_size / 1e3:.0f} KB)")


if __name__ == "__main__":
    build()
