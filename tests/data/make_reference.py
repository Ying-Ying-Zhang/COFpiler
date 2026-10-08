"""Regenerate tests/data/reference/*.extxyz by running the original 2022 scripts in legacy/.

The scripts need the package versions they were written for, e.g.

    python==3.10 numpy==1.23.5 pandas==1.5.3 scipy==1.10.1 ase==3.22.1 openpyxl

Usage:  python tests/data/make_reference.py
"""

import json
import runpy
import shutil
import sys
import tempfile
from pathlib import Path

import ase.visualize
import numpy as np
import pandas as pd
from ase.io import read, write

DATA = Path(__file__).resolve().parent
ROOT = DATA.parents[1]
REFERENCE = DATA / "reference"

# The scripts open the ASE GUI at the end.
ase.visualize.view = lambda *args, **kwargs: None
# With pandas < 3, the stacking script rotates rows of the input table in place
# through a view, so the result depends on the pandas version. Copy-on-write gives
# the pandas >= 3 behaviour, which the package reproduces.
pd.set_option("mode.copy_on_write", True)


def run_legacy(script: str, seed: int, argv: list[str]) -> None:
    np.random.seed(seed)
    sys.argv = [script, *argv]
    runpy.run_path(str(ROOT / "legacy" / script), run_name="__main__")


def main() -> None:
    cases = json.loads((REFERENCE / "cases.json").read_text())

    for case in cases["stacking"]:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            shutil.copy(DATA / "stacking_aa_ab.xlsx", work)
            shutil.copy(DATA / "hexagonal_monolayer.extxyz", work)
            run_legacy(
                "COFpiler.py",
                case["seed"],
                [
                    "-data",
                    "stacking_aa_ab.xlsx",
                    "-i",
                    "hexagonal_monolayer.extxyz",
                    "-path",
                    f"{work}/",
                    "-T",
                    str(case["temperature"]),
                    "-s",
                    case["symmetry"],
                    "-m",
                    "1" if case["mirror"] else "",  # the script parses this with bool()
                    "-L",
                    str(case["layers"]),
                    "-M",
                    "1",
                    "-o",
                    "extxyz",
                ],
            )
            out = work / "statistical_structures" / f"{case['layers']}_0.extxyz"
            write(REFERENCE / f"{case['name']}.extxyz", read(out))

    for case in cases["intercalate"]:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for f in (ROOT / "example").glob("*"):
                shutil.copy(f, work)
            run_legacy(
                "COFpiler_inter.py",
                case["seed"],
                [
                    "-data",
                    "example_intercalate.xlsx",
                    "-path",
                    f"{work}/",
                    "-T",
                    str(case["temperature"]),
                    "-L",
                    str(case["period"]),
                    "--set",
                    str(case["sets"]),
                    "-M",
                    "1",
                    "-o",
                    "extxyz",
                ],
            )
            (out,) = work.glob(f"*_1_{case['period']}l.extxyz")
            atoms = read(out)
            atoms.info["final_mode"] = out.name.rsplit("_1_", 1)[0]
            write(REFERENCE / f"{case['name']}.extxyz", atoms)


if __name__ == "__main__":
    main()
