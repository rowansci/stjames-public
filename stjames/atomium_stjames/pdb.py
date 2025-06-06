"""Contains functions for dealing with the .pdb file format."""

import re
from datetime import datetime
from itertools import chain, groupby
from typing import Any, Callable, TypedDict

from .data import CODES
from .mmcif import add_secondary_structure_to_polymers


def pdb_string_to_pdb_dict(filestring: str) -> dict[str, Any]:
    """Takes a .pdb filestring and turns into a ``dict`` which represents its
    record structure. Only lines which aren't empty are used.

    The resultant dictionary has line types as the keys, which point to the
    lines as its value. So ``{"TITLE": ["TITLE line 1", "TITLE line 2"]}`` etc.

    The exceptions are the REMARK records, where there is a sub-dictionary with
    REMARK numbers as keys, and the structure records themselves which are just
    arranged into lists - one for each model.

    :param str filestring: the .pdb filestring to process.
    :rtype: ``dict``"""

    pdb_dict: dict[str, Any] = {}
    lines_1 = list(filter(lambda l: bool(l.strip()), filestring.split("\n")))
    lines: list[list[str]] = [[line[:6].rstrip(), line.rstrip()] for line in lines_1]
    model_recs = ("ATOM", "HETATM", "ANISOU", "MODEL", "TER", "ENDMDL")
    for head, line in lines:
        if head == "REMARK":
            if "REMARK" not in pdb_dict:
                pdb_dict["REMARK"] = {}
            number = line.lstrip().split()[1]
            update_dict(pdb_dict["REMARK"], number, line)
        elif head in model_recs:
            if "MODEL" not in pdb_dict:
                pdb_dict["MODEL"] = [[]]
            if head == "ENDMDL":
                pdb_dict["MODEL"].append([])
            elif head != "MODEL":
                pdb_dict["MODEL"][-1].append(line)
        else:
            update_dict(pdb_dict, head, line)
    if "MODEL" in pdb_dict and not pdb_dict["MODEL"][-1]:
        pdb_dict["MODEL"].pop()
    return pdb_dict


def update_dict(d: dict[str, Any], key: str, value: str) -> None:
    """Takes a dictionary where the values are lists, and adds a value to one of
    the lists at the specific key. If the list doesn't exist, it creates it
    first.

    The dictionary is changed in place.

    :param dict d: the dictionary to update.
    :param str key: the location of the list.
    :param str value: the value to add to the list."""

    try:
        d[key].append(value)
    except Exception:
        d[key] = [value]


def pdb_dict_to_data_dict(pdb_dict: dict[str, Any]) -> dict[str, Any]:
    """Converts an .pdb dictionary into an atomium data dictionary, with the
    same standard layout that the other file formats get converted into.

    :param dict pdb_dict: the .pdb dictionary.
    :rtype: ``dict``"""

    data_dict = {
        "description": {"code": None, "title": None, "deposition_date": None, "classification": None, "keywords": [], "authors": []},
        "experiment": {"technique": None, "source_organism": None, "expression_system": None, "missing_residues": []},
        "quality": {"resolution": None, "rvalue": None, "rfree": None},
        "geometry": {"assemblies": [], "crystallography": {}},
        "models": [],
    }
    update_description_dict(pdb_dict, data_dict)
    update_experiment_dict(pdb_dict, data_dict)
    update_quality_dict(pdb_dict, data_dict)
    update_geometry_dict(pdb_dict, data_dict)
    update_models_list(pdb_dict, data_dict)
    return data_dict


def update_description_dict(pdb_dict: dict[str, Any], data_dict: dict[str, Any]) -> None:
    """Creates the description component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    extract_header(pdb_dict, data_dict["description"])
    extract_title(pdb_dict, data_dict["description"])
    extract_keywords(pdb_dict, data_dict["description"])
    extract_authors(pdb_dict, data_dict["description"])


def update_experiment_dict(pdb_dict: dict[str, Any], data_dict: dict[str, Any]) -> None:
    """Creates the experiment component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    extract_technique(pdb_dict, data_dict["experiment"])
    extract_source(pdb_dict, data_dict["experiment"])
    extract_missing_residues(pdb_dict, data_dict["experiment"])


def update_quality_dict(pdb_dict: dict[str, Any], data_dict: dict[str, Any]) -> None:
    """Creates the quality component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    extract_resolution_remark(pdb_dict, data_dict["quality"])
    extract_rvalue_remark(pdb_dict, data_dict["quality"])


def update_geometry_dict(pdb_dict: dict[str, Any], data_dict: dict[str, Any]) -> None:
    """Creates the geometry component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    extract_assembly_remark(pdb_dict, data_dict["geometry"])
    extract_crystallography(pdb_dict, data_dict["geometry"])


def update_models_list(pdb_dict: dict[str, Any], data_dict: dict[str, Any]) -> None:
    """Creates model dictionaries in a data dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    sequences = make_sequences(pdb_dict)
    secondary_structure = make_secondary_structure(pdb_dict)
    full_names = get_full_names(pdb_dict)
    for model_lines in pdb_dict["MODEL"]:
        aniso = make_aniso(model_lines)
        last_ter = get_last_ter_line(model_lines)
        model: dict[str, Any] = {"polymer": {}, "non_polymer": {}, "water": {}}
        for index, line in enumerate(model_lines):
            if line[:6] in ["ATOM  ", "HETATM"]:
                chain_id = line[21] if index < last_ter else id_from_line(line)
                res_id = id_from_line(line)
                if index < last_ter:
                    add_atom_to_polymer(line, model, chain_id, res_id, aniso, full_names)
                else:
                    add_atom_to_non_polymer(line, model, res_id, aniso, full_names)

            for chain_id, _chain in model["polymer"].items():
                _chain["sequence"] = sequences.get(chain_id, "")
        add_secondary_structure_to_polymers(model, secondary_structure)
        data_dict["models"].append(model)


def extract_header(pdb_dict: dict[str, Any], description_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds header information to it by parsing the HEADER
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("HEADER"):
        line = pdb_dict["HEADER"][0]
        if line[50:59].strip():
            description_dict["deposition_date"] = datetime.strptime(line[50:59], "%d-%b-%y").date()
        if line[62:66].strip():
            description_dict["code"] = line[62:66]
        if line[10:50].strip():
            description_dict["classification"] = line[10:50].strip()


def extract_title(pdb_dict: dict[str, Any], description_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds header information to it by parsing the TITLE
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("TITLE"):
        description_dict["title"] = merge_lines(pdb_dict["TITLE"], 10)


def extract_keywords(pdb_dict: dict[str, Any], description_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds header information to it by parsing the KEYWDS
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("KEYWDS"):
        text = merge_lines(pdb_dict["KEYWDS"], 10)
        description_dict["keywords"] = [w.strip() for w in text.split(",")]


def extract_authors(pdb_dict: dict[str, Any], description_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds header information to it by parsing the AUTHOR
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("AUTHOR"):
        text = merge_lines(pdb_dict["AUTHOR"], 10)
        description_dict["authors"] = [w.strip() for w in text.split(",")]


def extract_technique(pdb_dict: dict[str, Any], experiment_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds technique information to it by parsing EXPDTA
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("EXPDTA"):
        if pdb_dict["EXPDTA"][0].strip():
            experiment_dict["technique"] = pdb_dict["EXPDTA"][0][6:].strip()


def extract_source(pdb_dict: dict[str, Any], experiment_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds source information to it by parsing SOURCE
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("SOURCE"):
        data = merge_lines(pdb_dict["SOURCE"], 10)
        patterns = {"source_organism": r"ORGANISM_SCIENTIFIC\: (.+?);", "expression_system": r"EXPRESSION_SYSTEM\: (.+?);"}
        for attribute, pattern in patterns.items():
            matches = re.findall(pattern, data)
            if matches:
                experiment_dict[attribute] = matches[0]


def extract_missing_residues(pdb_dict: dict[str, Any], experiment_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds missing residue information to it by parsing
    REMARK 465 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    for line in pdb_dict.get("REMARK", {}).get("465", []):
        chunks = line.strip().split()
        if len(chunks) == 5:
            experiment_dict["missing_residues"].append({"name": chunks[2], "id": f"{chunks[3]}.{chunks[4]}"})


def extract_resolution_remark(pdb_dict: dict[str, Any], quality_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds resolution information to it by parsing REMARK
    2 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict quality_dict: the ``dict`` to update."""

    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("2"):
        for remark in pdb_dict["REMARK"]["2"]:
            try:
                quality_dict["resolution"] = float(remark[10:].strip().split()[1])
                break
            except Exception:
                pass


def extract_rvalue_remark(pdb_dict: dict[str, Any], quality_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds resolution information to it by parsing REMARK
    3 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict quality_dict: the ``dict`` to update."""

    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("3"):
        patterns = {
            "rvalue": r"R VALUE.+WORKING.+?: (.+)",
            "rfree": r"FREE R VALUE[ ]{2,}: (.+)",
        }
        for attribute, pattern in patterns.items():
            for remark in pdb_dict["REMARK"]["3"]:
                matches = re.findall(pattern, remark.strip())
                if matches:
                    try:
                        quality_dict[attribute] = float(matches[0].strip())
                    except Exception:
                        pass
                    break


def extract_assembly_remark(pdb_dict: dict[str, Any], geometry_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds assembly information to it by parsing REMARK
    350 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict geometry_dict: the ``dict`` to update."""
    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("350"):
        groups = [list(g) for k, g in groupby(pdb_dict["REMARK"]["350"], lambda x: "ECULE:" in x)][1:]
        assemblies = [list(chain(*a)) for a in zip(groups[::2], groups[1::2])]
        for a in assemblies:
            geometry_dict["assemblies"].append(assembly_lines_to_assembly_dict(a))


def assembly_lines_to_assembly_dict(lines: list[str]) -> dict[str, Any]:
    """Takes the lines representing a single biological assembly and turns
    them into an assembly dictionary.

    :param list lines: The REMARK lines to read.
    :rtype: ``dict``"""

    assembly: dict[str, Any] = {"transformations": [], "software": None, "buried_surface_area": None, "surface_area": None, "delta_energy": None, "id": 0}
    patterns: list[tuple[str, str, Callable[[str], Any]]] = [
        (r"(.+)SOFTWARE USED: (.+)", "software", lambda x: x),
        (r"(.+)BIOMOLECULE: (.+)", "id", int),
        (r"(.+)SURFACE AREA: (.+) [A-Z]", "buried_surface_area", float),
        (r"(.+)AREA OF THE COMPLEX: (.+) [A-Z]", "surface_area", float),
        (r"(.+)FREE ENERGY: (.+) [A-Z]", "delta_energy", float),
    ]
    t = None
    for line in lines:
        for pattern, key, converter in patterns:
            matches = re.findall(pattern, line)
            if matches:
                assembly[key] = converter(matches[0][1].strip())
        if "APPLY THE FOLLOWING" in line:
            if t:
                assembly["transformations"].append(t)
            t = {"chains": [], "matrix": [], "vector": []}
        if "CHAINS:" in line and t:
            t["chains"] += [c.strip() for c in line.split(":")[-1].strip().split(",") if c.strip()]
        if "BIOMT" in line and t:
            values = [float(x) for x in line.split()[4:]]
            if len(t["matrix"]) == 3:
                assembly["transformations"].append(t)
                t = {"chains": t["chains"], "matrix": [], "vector": []}
            t["matrix"].append(values[:3])
            t["vector"].append(values[-1])
    if t:
        assembly["transformations"].append(t)
    return assembly


def extract_crystallography(pdb_dict: dict[str, Any], geometry_dict: dict[str, Any]) -> None:
    """Takes a ``dict`` and adds assembly information to it by parsing the
    CRYST1 record.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict geometry_dict: the ``dict`` to update."""

    if pdb_dict.get("CRYST1"):
        line = pdb_dict["CRYST1"][0]
        values = line.split()
        geometry_dict["crystallography"]["space_group"] = line[55:66].strip()
        geometry_dict["crystallography"]["unit_cell"] = [float(val) for val in values[1:7]] if len(values) >= 6 else []


def make_sequences(pdb_dict: dict[str, Any]) -> dict[str, str]:
    """Creates a mapping of chain IDs to sequences, by parsing SEQRES records.

    :param dict pdb_dict: the .pdb dictionary to read.
    :rtype: ``dict``"""

    seq: dict[str, Any] = {}
    if pdb_dict.get("SEQRES"):
        for line in pdb_dict["SEQRES"]:
            chain, residues = line[11], line[19:].strip().split()
            if chain not in seq:
                seq[chain] = []
            seq[chain] += residues
    return {k: "".join([CODES.get(r, "X") for r in v]) for k, v in seq.items()}


def inverse_make_sequences(seq: str, chain_id: str) -> list[str]:
    """Converts a mapping of chain IDs to sequences back into SEQRES format.

    :param dict seq_dict: A dictionary mapping chain IDs to sequences.
    """
    # Reverse CODES dictionary
    REVERSE_CODES = {v: k for k, v in CODES.items()}

    seqres_lines = []
    residues = [REVERSE_CODES.get(aa, "UNK") for aa in seq]
    # SEQRES records are typically formatted into lines of up to 13 residues
    for i in range(0, len(residues), 13):
        seqres_lines.append(f"SEQRES {i // 13 + 1:>3} {chain_id} {len(seq):>4}  " + " ".join(residues[i : i + 13]))

    return seqres_lines


def make_secondary_structure(pdb_dict: dict[str, Any]) -> dict[str, Any]:
    """Creates a dictionary of helices and strands, with each having a list of
    start and end residues.

    :param pdb_dict: the .pdb dict to read.
    :rtype: ``dict``"""

    helices, strands = [], []
    for helix in pdb_dict.get("HELIX", []):
        helices.append(
            [
                f"{helix[19]}.{helix[21:25].strip()}{helix[25].strip()}",
                f"{helix[31]}.{helix[33:37].strip()}{helix[37].strip() if len(helix) > 37 else ''}",
            ]
        )
    for strand in pdb_dict.get("SHEET", []):
        strands.append(
            [
                f"{strand[21]}.{strand[22:26].strip()}{strand[26].strip()}",
                f"{strand[32]}.{strand[33:37].strip()}{strand[37].strip() if len(strand) > 37 else ''}",
            ]
        )
    return {"helices": helices, "strands": strands}


def get_full_names(pdb_dict: dict[str, Any]) -> dict[str, Any]:
    """Creates a mapping of het names to full English names.

    :param pdb_dict: the .pdb dict to read.
    :rtype: ``dict``"""

    full_names: dict[str, Any] = {}
    for line in pdb_dict.get("HETNAM", []):
        try:
            full_names[line[11:14].strip()] += line[15:].strip()
        except Exception:
            full_names[line[11:14].strip()] = line[15:].strip()

    return full_names


def make_aniso(model_lines: list[str]) -> dict[int, list[float]]:
    """Creates a mapping of chain IDs to anisotropy, by parsing ANISOU records.

    :param dict pdb_dict: the .pdb dictionary to read.
    :rtype: ``dict``"""

    return {int(line[6:11].strip()): [int(line[n * 7 + 28 : n * 7 + 35]) / 10000 for n in range(6)] for line in model_lines if line[:6] == "ANISOU"}


def get_last_ter_line(model_lines: list[str]) -> int:
    """Gets the index of the last TER record in a list of records. 0 will be
    returned if there are none.

    :param list model_lines: the lines to search.
    :rtype: ``int``"""

    last_ter = 0
    for index, line in enumerate(model_lines[::-1]):
        if line[:3] == "TER":
            last_ter = len(model_lines) - index - 1
            break
    return last_ter


def id_from_line(line: str) -> str:
    """Creates a residue ID from an atom line.

    :param str line: the ATOM or HETATM line record.
    :rtype: ``str``"""

    return "{}.{}{}".format(line[21], line[22:26].strip(), line[26].strip())


def add_atom_to_polymer(line: str, model: dict[Any, Any], chain_id: str, res_id: str, aniso_dict: dict[Any, Any], full_names: dict[Any, Any]) -> None:
    """Takes an .pdb ATOM or HETATM record, converts it, and adds it to a
    polymer dictionary.

    :param dict line: the line to read.
    :param dict model: the model to update.
    :param str chain_id: the chain ID to add to.
    :param str res_id: the molecule ID to add to.
    :param dict aniso_dict: lookup dictionary for anisotropy information."""

    atom = atom_line_to_dict(line, aniso_dict)

    try:
        model["polymer"][chain_id]["residues"][res_id]["atoms"][int(line[6:11])] = atom
    except Exception:
        name = line[17:20].strip()
        try:
            model["polymer"][chain_id]["residues"][res_id] = {
                "name": name,
                "full_name": full_names.get(name),
                "atoms": {int(line[6:11]): atom},
                "number": len(model["polymer"][chain_id]["residues"]) + 1,
            }
        except Exception:
            model["polymer"][chain_id] = {
                "internal_id": chain_id,
                "helices": [],
                "strands": [],
                "residues": {
                    res_id: {
                        "name": line[17:20].strip(),
                        "atoms": {int(line[6:11]): atom},
                        "number": 1,
                        "full_name": None,
                    }
                },
            }


def add_atom_to_non_polymer(line: str, model: dict[Any, Any], res_id: str, aniso_dict: dict[Any, Any], full_names: dict[Any, Any]) -> None:
    """Takes an .pdb ATOM or HETATM record, converts it, and adds it to a
    non-polymer dictionary.

    :param dict line: the line to read.
    :param dict model: the model to update.
    :param str res_id: the molecule ID to add to.
    :param dict aniso_dict: lookup dictionary for anisotropy information."""
    atom = atom_line_to_dict(line, aniso_dict)

    key = "water" if line[17:20] in ["HOH", "DOD"] else "non_polymer"
    try:
        model[key][res_id]["atoms"][int(line[6:11])] = atom
    except Exception:
        name = line[17:20].strip()
        model[key][res_id] = {
            "name": name,
            "full_name": full_names.get(name),
            "internal_id": line[21],
            "polymer": line[21],
            "atoms": {int(line[6:11]): atom},
        }


def guess_element_from_name(atom_name: str) -> str | None:
    atom_name = atom_name.strip()
    if not atom_name:
        return None

    # Case 1: Atom name starts with a digit (e.g. '1HG1') → element is second character
    if atom_name[0].isdigit() and len(atom_name) > 1:
        return atom_name[1].upper()

    # # Case 2: Atom name starts with a letter
    # if len(atom_name) >= 2 and atom_name[:2].isalpha():
    #     possible = atom_name[:2].strip().capitalize()
    #     # Check for common two-letter elements
    #     if possible in {"Cl", "Br", "Fe", "Mg", "Zn", "Ca", "Na", "Cu", "Mn", "Co", "Ni"}:
    #         return possible
    # Fallback to first letter
    return atom_name[0].upper()


class AtomDict(TypedDict, total=False):
    """A dictionary representing an atom in a PDB file."""

    occupancy: float | None
    bvalue: float | None
    charge: int | None
    anisotropy: float | None
    is_hetatm: bool | None
    name: str | None
    alt_loc: str | None
    x: float
    y: float
    z: float
    element: str | None


def atom_line_to_dict(line: str, aniso_dict: dict[Any, Any]) -> AtomDict:
    """
    Converts an ATOM or HETATM record to an atom dictionary.

    :param str line: the record to convert.
    :param dict aniso_dict: the anisotropy dictionary to use.
    :return: atom dictionary
    """

    a: AtomDict = {"occupancy": 1, "bvalue": None, "charge": 0, "anisotropy": aniso_dict.get(int(line[6:11].strip()), None)}
    a["is_hetatm"] = line[:6] == "HETATM"
    a["name"] = line[12:16].strip() or None
    a["alt_loc"] = line[16].strip() or None
    a["x"] = float(line[30:38].strip())
    a["y"] = float(line[38:46].strip())
    a["z"] = float(line[46:54].strip())
    if line[54:60].strip():
        a["occupancy"] = float(line[54:60].strip())
    if line[60:66].strip():
        a["bvalue"] = float(line[60:66].strip())
    a["element"] = line[76:78].strip() or None
    if not a["element"]:
        if not a["name"]:
            raise ValueError("Cannot guess element from empty name.")
        assert isinstance(a["name"], str)
        a["element"] = guess_element_from_name(a["name"])
    if line[78:80].strip():
        try:
            a["charge"] = int(line[78:80].strip())
        except Exception:
            a["charge"] = int(line[78:80][::-1].strip())

    if a["charge"] == 0:
        a["charge"] = None
    if not a["is_hetatm"]:
        a["is_hetatm"] = None
    if not a["alt_loc"]:
        a["alt_loc"] = None
    if a["occupancy"] == 1:
        a["occupancy"] = None
    if a["name"] == a["element"]:
        a["name"] = None

    return a


def merge_lines(lines: list[str], start: int, join: str = " ") -> str:
    """Gets a single continuous string from a sequence of lines.

    :param list lines: The lines to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :rtype: ``str``"""

    string = join.join([line[start:].strip() for line in lines])
    return string
