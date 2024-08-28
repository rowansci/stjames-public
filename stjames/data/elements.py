"""Read elemental data from files."""

import json
from collections import namedtuple
from importlib import resources

data_dir = resources.files("stjames").joinpath("data")

with data_dir.joinpath("symbol_element.json").open() as f:
    SYMBOL_ELEMENT: dict[str, int] = json.loads(f.read())

ELEMENT_SYMBOL = {v: k for k, v in SYMBOL_ELEMENT.items()}

Isotope = namedtuple("Isotope", ["relative_atomic_mass", "isotopic_composition", "standard_atomic_weight"])
with data_dir.joinpath("nist_isotopes.json").open() as f:
    d = json.loads(f.read())

    ISOTOPES: dict[int, dict[int, Isotope]] = {
        int(k): {
            int(kk): Isotope(*vv)
            for kk, vv in v.items()  # stay open
        }
        for k, v in d.items()
    }

with data_dir.joinpath("bragg_radii.json").open() as f:
    BRAGG_RADII: dict[int, float] = json.loads(f.read())
