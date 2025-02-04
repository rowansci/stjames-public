# mypy: ignore-errors
# fmt: off
"""Contains various file handling helper functions."""

import builtins
import gzip
from typing import Any

from requests import get

from .data import File, data_dict_to_file
from .mmcif import mmcif_dict_to_data_dict, mmcif_string_to_mmcif_dict
from .mmtf import mmtf_bytes_to_mmtf_dict, mmtf_dict_to_data_dict
from .pdb import pdb_dict_to_data_dict, pdb_string_to_pdb_dict


def open(path, *args, **kwargs) -> (File | dict[str, Any] | Any | dict):
    """Opens a file at a given path, works out what filetype it is, and parses
    it accordingly.

    For example:
    open('/path/to/file.pdb', data_dict=True)

    This will parse file.pdb as a .pdb file, but only go as far as converting it
    to an atomium data dictionary.

    If the file extension is .gz, the file will be unzipped first.

    :param str path: the location of the file.
    :param bool file_dict: if ``True``, parsing will stop at the file ``dict``.
    :param bool data_dict: if ``True``, parsing will stop at the data ``dict``.
    :rtype: ``File``"""

    if str(path)[-3:] == ".gz":
        try:
            with gzip.open(path) as f:
                filestring = f.read().decode()
        except Exception:
            with gzip.open(path, "rt") as f:
                filestring = f.read()
        return parse_string(filestring, path[:-3], *args, **kwargs)
    else:
        try:
            with builtins.open(path) as f:
                filestring = f.read()
        except Exception:
            with builtins.open(path, "rb") as f:
                filestring = f.read()
        return parse_string(filestring, path, *args, **kwargs)


def fetch(code, *args, **kwargs) -> (File | dict[str, Any] | Any | dict):
    """Fetches a file from a remote location via HTTP.

    If a PDB code is given, the .cif form of that struture will be fetched from
    the RCSB servers. If that code is given an extension, that file format will
    be obtained instead of .cif. If a URL is given, the function will simply
    look in that location.

    For example:
    fetch('1lol.mmtf', file_dict=True)

    This will get the .mmtf version of structure 1LOL, but only go as far as
    converting it to an atomium file dictionary.

    :param str code: the file to fetch.
    :param bool file_dict: if ``True``, parsing will stop at the file ``dict``.
    :param bool data_dict: if ``True``, parsing will stop at the data ``dict``.
    :raises ValueError: if no file is found.
    :rtype: ``File``"""

    if code.startswith("http"):
        url = code
    elif code.endswith(".mmtf"):
        url = "https://mmtf.rcsb.org/v1.0/full/{}".format(code[:-5].lower())
    else:
        if "." not in code:
            code += ".cif"
        url = "https://files.rcsb.org/view/" + code.lower()
    response = get(url, stream=True)
    if response.status_code == 200:
        text = response.content if code.endswith(".mmtf") else response.text
        return parse_string(text, code, *args, **kwargs)
    raise ValueError("Could not find anything at {}".format(url))


def parse_string(filestring, path, file_dict=False, data_dict=False):
    """Takes a filestring and parses it in the appropriate way. You must provide
    the string to parse itself, and some other string that ends in either .cif,
    .mmtf, or .cif - that will determine how the file is parsed.

    (If this cannot be inferred from the path string, atomium will guess based
    on the filestring contents.)

    :param str filestring: the contents of some file.
    :param str path: the filename of the file of origin.
    :param bool file_dict: if ``True``, parsing will stop at the file ``dict``.
    :param bool data_dict: if ``True``, parsing will stop at the data ``dict``.
    :rtype: ``File``"""

    file_func, data_func = get_parse_functions(filestring, path)
    parsed = file_func(filestring)
    if not file_dict:
        parsed = data_func(parsed)
        if not data_dict:
            filetype = data_func.__name__.split("_")[0].replace("mmc", "c")
            parsed = data_dict_to_file(parsed, filetype)
    return parsed


def get_parse_functions(filestring, path):
    """Works out which parsing functions to use for a given filestring and
    returns them.

    (If this cannot be inferred from the path string, atomium will guess based
    on the filestring contents.)

    :param str filestring: the filestring to inspect.
    :param str path: the path to inspect.
    :rtype: ``tuple``"""

    if "." in path:
        ending = path.split(".")[-1]
        if ending in ("mmtf", "cif", "pdb"):
            return {
                "cif": (mmcif_string_to_mmcif_dict, mmcif_dict_to_data_dict),
                "mmtf": (mmtf_bytes_to_mmtf_dict, mmtf_dict_to_data_dict),
                "pdb": (pdb_string_to_pdb_dict, pdb_dict_to_data_dict),
            }[ending]
    if isinstance(filestring, bytes):
        return (mmtf_bytes_to_mmtf_dict, mmtf_dict_to_data_dict)
    elif "_atom_sites" in filestring:
        return (mmcif_string_to_mmcif_dict, mmcif_dict_to_data_dict)
    else:
        return (pdb_string_to_pdb_dict, pdb_dict_to_data_dict)
