"""
This code populates ELEMENT_DICTIONARY and ISOTOPE_DICTIONARY from a static datafile.

Adapted from peregrine
"""

from importlib.resources import files

from stjames import data

# Index -> Symbol
ELEMENT_DICTIONARY = {}
ISOTOPE_DICTIONARY = {}

with files(data).joinpath("isotopes.csv").open() as isotope_file:
    prev_number = 0
    current_dict = {}
    for line in isotope_file:
        symbol, number_str, mass, abundance = line.split(",")

        if symbol == "Symbol":
            continue

        number = int(number_str)
        ELEMENT_DICTIONARY[number] = symbol

        if number == prev_number:
            current_dict[float(mass)] = float(abundance.rstrip())
        else:
            if prev_number:
                ISOTOPE_DICTIONARY[prev_number] = current_dict
            current_dict = {}
            current_dict[float(mass)] = float(abundance.rstrip())

        prev_number = number

    ISOTOPE_DICTIONARY[prev_number] = current_dict
    ELEMENT_DICTIONARY[0] = "Bq"  # like Banquo, get it?

# Symbol -> Index
INV_ELEMENT_DICTIONARY = {v: int(k) for k, v in ELEMENT_DICTIONARY.items()}


def get_atomic_symbol(atomic_number: int) -> str:
    """Gets element symbol from a given atomic number."""
    if atomic_number in ELEMENT_DICTIONARY:
        return ELEMENT_DICTIONARY[atomic_number]
    else:
        raise ValueError(f"Unknown atomic number: ``{atomic_number}``")


def get_atomic_number(atomic_symbol: str) -> int:
    """Gets atomic number from a given element symbol (converted to titlecase using ``string.title()``)."""
    if atomic_symbol.title() in INV_ELEMENT_DICTIONARY:
        return INV_ELEMENT_DICTIONARY[atomic_symbol.title()]
    else:
        raise ValueError(f"Unknown atomic symbol: ``{atomic_symbol}``")


# values in .csv from https://github.com/psi4/psi4/blob/70d1c1448ea0388505b826b24d7f6dd2db723c9b/psi4/src/psi4/libfock/cubature.cc#L112

BRAGG_DICTIONARY = {}
with files(data).joinpath("bragg_radii.csv").open() as isotope_file:
    for line in isotope_file:
        number_str, radius = line.split(",")
        BRAGG_DICTIONARY[int(number_str)] = float(radius)


def get_bragg_radius(z: int) -> float:
    if z in BRAGG_DICTIONARY:
        return BRAGG_DICTIONARY[z]
    else:
        raise ValueError(f"No Bragg radius defined for Z=``{z}``")
