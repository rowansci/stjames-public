"""Read elemental data from files."""

import json
from importlib import resources

data_dir = resources.files("stjames").joinpath("data")

with data_dir.joinpath("symbol_element.json").open() as f:
    SYMBOL_ELEMENT: dict[str, int] = json.loads(f.read())

ELEMENT_SYMBOL = {v: k for k, v in SYMBOL_ELEMENT.items()}

with data_dir.joinpath("isotopes.json").open() as f:
    d = json.loads(f.read())
    # Convert masses to floats
    ISOTOPES: dict[int, dict[float, float]] = {k: {float(kk): vv for kk, vv in v.items()} for k, v in d.items()}

with data_dir.joinpath("bragg_radii.json").open() as f:
    BRAGG_RADII: dict[int, float] = json.loads(f.read())
