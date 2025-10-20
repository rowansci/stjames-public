"""
Read the NIST isotopes data file and write it to a JSON file.

NIST Isotopes data from:
https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2
"""

import json
from collections import defaultdict
from importlib import resources
from typing import Callable, TypeVar

data_dir = resources.files("stjames").joinpath("data")

_T = TypeVar("_T")


def process_line(line: str, fmt: Callable[[str], _T] = str) -> _T:  # type: ignore[assignment]
    """
    Process a line from the NIST data file.

    :param line: line to process
    :param fmt: function to format the value
    >>> process_line("Atomic Number = 1", int)
    1
    """
    return fmt(line.split("=")[-1].strip())


def fmt_float(val: str) -> float:
    """
    Format a float from the NIST data file.

    >>> fmt_float(" 1.00784(7)")
    1.00784
    """
    return float(val.strip().split("(")[0])


def fmt_maybe_list(val: str) -> float:
    """
    Format a float or list of floats from the NIST data file.

    Only the first value is returned.

    >>> fmt_maybe_list("1.00784(7)")
    1.00784
    >>> fmt_maybe_list(" [1.00784,1.00811]")
    1.00784
    >>> fmt_maybe_list(" [98]")
    98.0
    """
    val = val.strip()
    if val.startswith("["):
        val = val[1:-1].split(",")[0]
    return fmt_float(val)


def process_chunk(chunk: str) -> tuple[int, int, tuple[float, float, float]]:
    r"""
    Atomic Number, Mass Number, (Relative Atomic Mass, Isotopic Composition, Standard Atomic Weight)

    >>> process_chunk('''\
    ... Atomic Number = 1
    ... Atomic Symbol = H
    ... Mass Number = 1
    ... Relative Atomic Mass = 1.00784(7)
    ... Isotopic Composition = 0.999885(70)
    ... Standard Atomic Weight = [1.00784,1.00811]
    ... Notes = m
    ... ''')
    (1, 1, (1.00784, 0.999885, 1.00784))
    """
    lines = chunk.splitlines()

    atomic_number = process_line(lines[0], int)
    _atomic_symbol = process_line(lines[1], str)
    mass_number = process_line(lines[2], int)
    relative_atomic_mass = process_line(lines[3], fmt_float)
    try:
        isotopic_composition = process_line(lines[4], fmt_float)
    except ValueError:
        isotopic_composition = 0
    try:
        standard_atomic_weight = process_line(lines[5], fmt_maybe_list)
    except ValueError:
        standard_atomic_weight = relative_atomic_mass

    return atomic_number, mass_number, (relative_atomic_mass, isotopic_composition, standard_atomic_weight)


def read_nist_isotopes() -> dict[int, dict[int, tuple[float, float, float]]]:
    """
    Read the NIST data file and write it to a JSON file.

    {Atomic Number: {Mass Number, (Relative Atomic Mass, Isotopic Composition, Standard Atomic Weight)}}
    """
    with data_dir.joinpath("nist_isotopes.txt").open() as f:
        next(f), next(f)  # Skip the first two lines
        nist_isotopes = f.read()

    isotopes: dict[int, dict[int, tuple[float, float, float]]] = defaultdict(dict)
    for chunk in nist_isotopes.split("\n\n"):
        atomic_number, mass_number, values = process_chunk(chunk)
        isotopes[atomic_number][mass_number] = values

    with open("nist_isotopes.json", "w") as f:
        json.dump(isotopes, f)

    return isotopes


if __name__ == "__main__":
    from pprint import pprint

    pprint(read_nist_isotopes())
