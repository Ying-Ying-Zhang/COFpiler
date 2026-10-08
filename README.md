# COFpiler

[![CI](https://github.com/Ying-Ying-Zhang/COFpiler/actions/workflows/ci.yml/badge.svg)](https://github.com/Ying-Ying-Zhang/COFpiler/actions/workflows/ci.yml)

Build statistical structure models of layered materials such as covalent organic frameworks (COFs). Each layer is stacked on the previous one with a shift drawn from the Boltzmann distribution of the stacking modes. The method is described in:

> Y. Zhang, M. Položij, T. Heine, *Statistical Representation of Stacking Disorder in Layered Covalent Organic Frameworks*, Chem. Mater. **34**, 2376–2381 (2022). [doi:10.1021/acs.chemmater.1c04365](https://doi.org/10.1021/acs.chemmater.1c04365)

## Installation

Python 3.10 or newer.

```
pip install git+https://github.com/Ying-Ying-Zhang/COFpiler
```

## Input data

A table of stacking modes, as the sheet `data` of an `.xlsx` file or as a `.csv` file:

| Column | Meaning | Default unit |
| --- | --- | --- |
| `stacking_type` | Name, e.g. `AA_s1_serr` or `AB_ecl`. Contains `AA` or `AB_`, and one of `ecl`, `s1`, `s2`. | |
| `Erel` | Relative energy of the stacking mode | kJ/mol |
| `x`, `y`, `z` | Shift between two neighbouring layers | Å |

Other units can be given with `--energy-unit` and `--length-unit`, for example `eV` or `nm`. Energies per stacking unit (eV) are converted to kJ/mol.

The table is checked when it is read. Units must match, shifts need three finite components with a positive `z`, and names must be unique.

## Command line

Stack a monolayer:

```
cofpiler -data stacking_aa_ab.xlsx -i hexagonal_monolayer.extxyz -path tests/data -T 293 -s C6 -L 10 -M 5 --seed 1
```

This writes the models to `tests/data/statistical_structures/10_0.cif` … `10_4.cif`, with `record_10.txt` listing the stacking type and shift drawn for every layer.

| Option | Meaning |
| --- | --- |
| `-data` | Table of stacking modes |
| `-i` | Monolayer structure, any format ASE reads |
| `-path` | Folder with the input files; output goes here too (default: current folder) |
| `-T` | Synthesis temperature in K (default: 293) |
| `-s` | Symmetry of the layer: `C3`, `C4` or `C6` (default: `C6`) |
| `-L` | Number of layers |
| `-M` | Number of models (default: 1) |
| `--mirror` / `--no-mirror` | Mirror `s2` shifts before rotating them (default: on) |
| `-mp X Y` | Direction of the mirror line (default: `0.001 1`) |
| `--ab-slip axis\|diag` | For `C4`: which AB slip the table describes |
| `-o` | Output format (default: `cif`) |
| `--seed` | Random seed, for reproducible models |
| `--view` | Open the last model in the ASE GUI |

Stack with an intercalated monomer:

```
cofpiler-intercalate -data example_intercalate.xlsx -path example -L 3 --set 2 -M 2 --seed 0
```

For every stacking type `t` in the table, the folder holds `t1.xyz` (layers below the monomer), `t2.xyz` (the monomer) and `t3.xyz` (layers above). Each set has `L + 1` layers, one of them the monomer.

## Python

```python
import numpy as np
from ase.io import read
from cofpiler import StackingTable, build_stacked_model, ureg

table = StackingTable.read("tests/data/stacking_aa_ab.xlsx")
print(table.probabilities(ureg.Quantity(20, "degC")))

structure, records = build_stacked_model(
    read("tests/data/hexagonal_monolayer.extxyz"),
    table,
    n_layers=10,
    temperature=293,
    rng=np.random.default_rng(1),
)
structure.write("model.cif")
```

## Development

```
pip install -e ".[dev]"
pytest
ruff check . && ruff format --check .
mypy
```

The tests compare the package with reference structures made by the original 2022 scripts, which are kept unchanged in `legacy/`. `tests/data/make_reference.py` regenerates them.

## Changes from the 2022 scripts

- `python COFpiler.py …` is now `cofpiler …`, and `python COFpiler_inter.py …` is `cofpiler-intercalate …`.
- With `--seed N`, the models are identical to the 2022 scripts run after `np.random.seed(N)` with pandas 3. With pandas < 3, the old stacking script rotated rows of the input table in place, so its results depended on the pandas version.
- `-m False` used to switch mirroring on. Use `--no-mirror`. `-mp` now takes two numbers.
- For `C4`, the AB slip direction is given with `--ab-slip` instead of an interactive prompt.
- The mirror image of a shift is computed in closed form instead of with `scipy.optimize.fsolve`. The results agree to 1e-9 Å.
- CSV tables are read, as the old README promised.
- Unknown stacking types stop with an error message instead of a `NameError`.
- The ASE GUI only opens with `--view`. The unused options `-i`, `-m` and `-mp` of the intercalation script are gone.

## Known issues

These come from the 2022 scripts and are kept so that results stay reproducible.

1. Stacking types named `ABC…` are not supported. The 2022 script had no code path for them and failed at the first ABC layer. `example/example.xlsx` contains ABC types, so the stacking builder rejects it.
2. For `AB_` types, the rotated lattice offset replaces the in-plane shift from the table. Only the `z` column of AB rows has an effect.
3. The random rotation step is π/3 for every symmetry. For `C3` and `C4` the symmetric step would be 2π/3 and π/2, which the intercalation builder uses.
4. For `C4`, choosing an AB slip direction sets the energy of the other direction's `AB_s2` mode to zero, which makes that mode the most likely one.
5. The Boltzmann factor uses 1/R = 120.37 K mol/kJ. The CODATA value is 120.27 K mol/kJ, a 0.08 % difference.
6. In the intercalation builder, a mirrored monomer keeps its unmirrored shift, and each monomer is rotated twice with only the second angle kept.

## Citation

If you use COFpiler, please cite the paper above. Citation metadata for the software is in `CITATION.cff`.

## License

MIT © Yingying Zhang
