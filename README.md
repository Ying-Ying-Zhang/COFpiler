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
| `stacking_type` | Name, e.g. `AA_s1_serr`, `AB_ecl` or `ABC_s1`. Contains `AA`, `AB_` or `ABC`, and one of `ecl`, `s1`, `s2`. | |
| `Erel` | Relative energy of the stacking mode | kJ/mol |
| `x`, `y`, `z` | Shift between two neighbouring layers | Å |

Other units can be given with `--energy-unit` and `--length-unit`, for example `eV` or `nm`. Energies per stacking unit (eV) are converted to kJ/mol.

The table is checked when it is read. Units must match, shifts need three finite components with a positive `z`, and names must be unique.

AB and ABC layers sit over the pore of the layer below: they shift by [1/3, 1/3] of the cell (hexagonal lattices) on top of their slip. AB alternates the direction of this offset from layer to layer, ABC keeps it. Since every layer draws a random symmetry-equivalent direction, both give the same set of interlayer shifts and differ in energy, as AB(C)<sub>stat</sub> in the paper. If the table has AB or ABC rows, `--ab-shift` says how to read their `x` and `y`:

- `--ab-shift slip`: `x`, `y` are the slip on top of the AB position, like the AA rows.
- `--ab-shift full`: `x`, `y` are the whole interlayer shift, including the offset.

## Command line

Stack a monolayer:

```
cofpiler -data stacking_aa_ab.xlsx -i hexagonal_monolayer.extxyz -path tests/data -T 293 -s C6 -L 10 -M 5 --ab-shift slip --seed 1
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
| `--ab-shift slip\|full` | For tables with AB or ABC rows: how to read their `x`, `y` (see above) |
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
    ab_shift="slip",
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

The tests compare the package with reference structures made by the original 2022 scripts, which are kept unchanged in `legacy/`. `tests/data/make_reference.py` regenerates them. For stacking, the comparison covers models of AA layers; AB and ABC layers are tested against the definitions above.

## Changes from the 2022 scripts

- `python COFpiler.py …` is now `cofpiler …`, and `python COFpiler_inter.py …` is `cofpiler-intercalate …`.
- AB layers keep the slip from the table and add the offset to the AB position. The 2022 script replaced the slip by the offset, so only the `z` of AB rows had an effect.
- ABC stacking types are supported for `C3` and `C6`. The 2022 script had no code path for them and failed at the first ABC layer.
- With `--seed N`, models of AA layers are identical to the 2022 scripts run after `np.random.seed(N)` with pandas 3. With pandas < 3, the old stacking script rotated rows of the input table in place, so its results depended on the pandas version.
- `-m False` used to switch mirroring on. Use `--no-mirror`. `-mp` now takes two numbers.
- For `C4`, the AB slip direction is given with `--ab-slip` instead of an interactive prompt.
- The mirror image of a shift is computed in closed form instead of with `scipy.optimize.fsolve`. The results agree to 1e-9 Å.
- CSV tables are read, as the old README promised.
- Unknown stacking types stop with an error message.
- The ASE GUI only opens with `--view`. The unused options `-i`, `-m` and `-mp` of the intercalation script are gone.

## Known issues

These come from the 2022 scripts and are kept so that results stay reproducible.

1. The random rotation step is π/3 for every symmetry. For `C3` and `C4` the symmetric step would be 2π/3 and π/2, which the intercalation builder uses.
2. For `C4`, choosing an AB slip direction sets the energy of the other direction's `AB_s2` mode to zero, which makes that mode the most likely one.
3. The Boltzmann factor uses 1/R = 120.37 K mol/kJ. The CODATA value is 120.27 K mol/kJ, a 0.08 % difference.
4. In the intercalation builder, a mirrored monomer keeps its unmirrored shift, and each monomer is rotated twice with only the second angle kept.
5. ABC stacking on square (`C4`) lattices, with an offset of [1/3, 0] of the cell, is not implemented.

## Citation

If you use COFpiler, please cite the paper above. Citation metadata for the software is in `CITATION.cff`.

## License

MIT © Yingying Zhang
