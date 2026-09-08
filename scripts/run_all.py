#!/usr/bin/env python
"""Run the scripted parts of the analysis.

    python scripts/run_all.py --step all
    python scripts/run_all.py --list

Run from the repository root with the `endo` environment active. Each step
writes only into `results/` and `data/interim/` and prints the numbers that
appear in the manuscript, so a reader can compare them line by line.

WHAT IS SCRIPTED HERE

    download    fetch the accessions listed in data/README.md into data/raw/
    atlas       QC, Harmony integration, lineage and immune-state annotation
                -> checkpoints/full_atlas.h5ad, checkpoints/immune_annotated.h5ad
    epithelial  epithelial subset and inferred copy number
    statistics  the statistical core, in four parts (see 05_statistics.py):
                per-sample composition and the ten-state scan (Tables S1, S5),
                power and minimum detectable effect (Table S6), scCODA
                (Table S7), cycle phase (Table S8)
    workbook    collect Tables S1-S8 into supplementary/*.xlsx

    propeller   NOT run from here: it needs the R environment.
                Rscript scripts/propeller.R  (the other half of Table S7)

WHAT IS NOT SCRIPTED YET

The GeoMx spatial tests and their sensitivity analyses (Tables S2, S4), the two
non-replicating validation cohorts, the EAOC extension with its survival curves
(Supplementary Figure S2), the cross-disease Census benchmark (Supplementary
Figure S1) and the scIB benchmark (Table S3) were run interactively. Their
outputs are the CSVs in `results/tables/`, which the workbook step consumes, but
the code that produced them has not been extracted into scripts in this
checkout. Treat those tables as data, not as regenerable output, until it has.

Figures are built by the `build_fig*_print.py` modules, which are imported
rather than run as programs; `docs/FIGURE_PROVENANCE.md` records the checksum of
every exported file.

`03_trajectory_dnb.py` is deliberately absent from the step list: it implements
the critical-transition analysis that a two-population mixing control refuted,
and is kept only for provenance.

Timing on 128 cores / 251 GB (no GPU): download ~20 min for 3.3 GB, atlas ~40
min, statistics ~5 min, everything else seconds.
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
R_DEFAULT = "Rscript"

STEPS = [
    ("download", "00_download_data.py"),
    ("atlas", "01_load_qc_integrate.py"),
    ("epithelial", "02_epithelial_infercnv.py"),
    ("statistics", "05_statistics.py"),
    ("workbook", "06_supplementary_workbook.py"),
]


def run(script: str, r_bin: str) -> None:
    """Run one step. A missing script is a bug, not something to skip past."""
    path = SCRIPTS / script
    if not path.exists():
        raise SystemExit(f"{path} is missing; this checkout is incomplete")
    print(f"\n=== {script} ===", flush=True)
    subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
    if script == "05_statistics.py":
        print(f"\n[note] the propeller half of Table S7 needs R: "
              f"{r_bin} scripts/propeller.R", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", default="all",
                    help="step name, comma-separated list, or 'all'")
    ap.add_argument("--r-bin", default=R_DEFAULT,
                    help="Rscript executable used by the propeller step")
    ap.add_argument("--list", action="store_true", help="list steps and exit")
    args = ap.parse_args()

    if args.list:
        for name, script in STEPS:
            print(f"  {name:12s} {script}")
        return

    wanted = [n for n, _ in STEPS] if args.step == "all" else args.step.split(",")
    unknown = [w for w in wanted if w not in dict(STEPS)]
    if unknown:
        sys.exit(f"unknown step(s): {', '.join(unknown)}")
    for name, script in STEPS:
        if name in wanted:
            run(script, args.r_bin)
    print("\ndone", flush=True)


if __name__ == "__main__":
    main()
