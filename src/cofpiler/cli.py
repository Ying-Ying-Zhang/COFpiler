"""Command-line entry points `cofpiler` and `cofpiler-intercalate`."""

from __future__ import annotations

import argparse
import re
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO, cast

import numpy as np
from ase import Atoms
from ase.io import read, write

from cofpiler.intercalate import build_intercalated_model
from cofpiler.models import ENERGY_UNIT, LENGTH_UNIT, StackingTable
from cofpiler.stacking import LayerRecord, build_stacked_model


def _common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--data",
        "-data",
        required=True,
        help="table of stacking types (.xlsx sheet 'data', or .csv) with columns "
        "stacking_type, Erel, x, y, z",
    )
    parser.add_argument(
        "--path",
        "-path",
        type=Path,
        default=Path("."),
        help="folder with the input files; output is written here too",
    )
    parser.add_argument(
        "-T", "--tem", type=float, default=293, help="synthesis temperature in K (default: 293)"
    )
    parser.add_argument(
        "-o",
        "--outstr_format",
        default="cif",
        help="output structure format understood by ASE (default: cif)",
    )
    parser.add_argument("-M", type=int, default=1, help="number of models to build (default: 1)")
    parser.add_argument("--seed", type=int, help="random seed, for reproducible models")
    parser.add_argument(
        "--energy-unit",
        default=ENERGY_UNIT,
        help=f"unit of the Erel column (default: {ENERGY_UNIT})",
    )
    parser.add_argument(
        "--length-unit",
        default=LENGTH_UNIT,
        help=f"unit of the x, y, z columns (default: {LENGTH_UNIT})",
    )
    parser.add_argument("--view", action="store_true", help="open the last model in the ASE GUI")


def _read_table(args: argparse.Namespace) -> StackingTable:
    return StackingTable.read(args.path / args.data, args.energy_unit, args.length_unit)


def _write_header(f: TextIO, table: StackingTable, temperature: float) -> None:
    f.write(f"temperature: {temperature} K\n")
    f.write("stacking_type\tErel/(kJ/mol)\tprobability\tx/A\ty/A\tz/A\n")
    for mode, energy, p, shift in zip(
        table.names,
        table.energies(),
        table.probabilities(temperature),
        table.shifts(),
        strict=True,
    ):
        f.write(f"{mode}\t{energy:.6f}\t{p:.6e}\t" + "\t".join(f"{v:.6f}" for v in shift) + "\n")


def _write_records(f: TextIO, model: int, records: list[LayerRecord]) -> None:
    f.write(f"\nmodel {model}\nlayer\tstacking_type\tshift\ttotal_shift\n")
    for r in records:
        f.write(f"{r.layer}\t{r.mode}\t{np.round(r.shift, 6)}\t{np.round(r.total_shift, 6)}\n")


def _view(atoms: Atoms) -> None:
    from ase.visualize import view

    view(atoms)


def stack_main(argv: Sequence[str] | None = None) -> None:
    """Build statistically stacked models from a monolayer."""
    parser = argparse.ArgumentParser(prog="cofpiler", description=stack_main.__doc__)
    _common_arguments(parser)
    parser.add_argument("-i", "--instr", required=True, help="the monolayer structure")
    parser.add_argument("-L", type=int, required=True, help="number of layers")
    parser.add_argument(
        "-s", "--symmetry", default="C6", help="symmetry of the layer: C3, C4 or C6 (default: C6)"
    )
    parser.add_argument(
        "-m",
        "--mirror",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="mirror s2 shifts before rotating them (default: on)",
    )
    parser.add_argument(
        "-mp",
        "--mplane",
        type=float,
        nargs=2,
        default=[0.001, 1.0],
        metavar=("X", "Y"),
        help="direction of the mirror line (default: 0.001 1)",
    )
    parser.add_argument(
        "--ab-slip", choices=["axis", "diag"], help="for C4 only: which AB slip the table describes"
    )
    parser.add_argument(
        "--ab-shift",
        choices=["slip", "full"],
        help="for tables with AB or ABC types: their x, y are the slip on top of the AB "
        "position (slip) or include the offset to it (full)",
    )
    args = parser.parse_args(argv)
    try:
        _stack(args)
    except (OSError, ValueError, NotImplementedError) as e:
        parser.error(str(e))


def _stack(args: argparse.Namespace) -> None:
    table = _read_table(args)
    monolayer = cast(Atoms, read(args.path / args.instr))
    symmetry = int(re.sub(r"\D", "", args.symmetry))
    rng = np.random.RandomState(args.seed)
    out = args.path / "statistical_structures"
    out.mkdir(exist_ok=True)

    with open(out / f"record_{args.L}.txt", "w") as f:
        _write_header(f, table, args.tem)
        for model in range(args.M):
            structure, records = build_stacked_model(
                monolayer,
                table,
                args.L,
                args.tem,
                symmetry=symmetry,
                mirror=args.mirror,
                mirror_plane=args.mplane,
                ab_slip=args.ab_slip,
                ab_shift=args.ab_shift,
                rng=rng,
            )
            filename = out / f"{args.L}_{model}.{args.outstr_format}"
            write(filename, structure, format=args.outstr_format)
            _write_records(f, model, records)
            print(f"wrote {filename}")
    if args.view:
        _view(structure)


def intercalate_main(argv: Sequence[str] | None = None) -> None:
    """Build statistical models with an intercalated monomer."""
    parser = argparse.ArgumentParser(
        prog="cofpiler-intercalate", description=intercalate_main.__doc__
    )
    _common_arguments(parser)
    parser.add_argument("-L", type=int, required=True, help="a monomer layer every L + 1 layers")
    parser.add_argument("--set", type=int, required=True, help="number of sets of L + 1 layers")
    args = parser.parse_args(argv)
    try:
        _intercalate(args)
    except (OSError, ValueError) as e:
        parser.error(str(e))


def _intercalate(args: argparse.Namespace) -> None:
    table = _read_table(args)
    rng = np.random.RandomState(args.seed)

    with open(args.path / f"record_{args.L}.txt", "w") as f:
        _write_header(f, table, args.tem)
        for model in range(1, args.M + 1):
            structure, mode, records = build_intercalated_model(
                args.path, table, args.L, args.set, args.tem, rng=rng
            )
            filename = args.path / f"{mode}_{model}_{args.L}l.{args.outstr_format}"
            write(filename, structure, format=args.outstr_format)
            _write_records(f, model, records)
            print(f"wrote {filename}")
    if args.view:
        _view(structure)
