import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
REFERENCE = DATA / "reference"
EXAMPLE = ROOT / "example"
CASES = json.loads((REFERENCE / "cases.json").read_text())


@pytest.fixture
def aa_ab_table_path() -> Path:
    return DATA / "stacking_aa_ab.xlsx"


@pytest.fixture
def monolayer_path() -> Path:
    return DATA / "hexagonal_monolayer.extxyz"
