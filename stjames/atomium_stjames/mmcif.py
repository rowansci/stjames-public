# mypy: ignore-errors

"""Contains functions for dealing with the .cif file format."""

import re
from collections import deque
from datetime import datetime
from typing import Any

import numpy as np

from .data import CODES


def mmcif_string_to_mmcif_dict(filestring) -> dict:
    """Takes a .cif filestring and turns into a ``dict`` which represents its
    table structure. Only lines which aren't empty and which don't begin with
    ``#`` are used.

    Multi-line strings are consolidated onto one line, and the whole thing is
    then split into the blocks that will become table lists. At the end, quote
    marks are removed from any string which retains them.

    :param str filestring: the .cif filestring to process.
    :rtype: ``dict``"""

    lines = deque(filter(lambda l: l and l[0] != "#", filestring.split("\n")))
    lines = consolidate_strings(lines)
    blocks = mmcif_lines_to_mmcif_blocks(lines)
    mmcif_dict = {}
    for block in blocks:
        if block["lines"][0] == "loop_":
            mmcif_dict[block["category"]] = loop_block_to_list(block)
        else:
            mmcif_dict[block["category"]] = non_loop_block_to_list(block)
    strip_quotes(mmcif_dict)
    return mmcif_dict


def consolidate_strings(lines):
    """Generally, .cif files have a one file line to one table row
    correspondence. Sometimes however, a string cell is given a line of its own,
    breaking the row over several lines. This function takes the lines of a .cif
    file and puts all table rows on a single line.

    :param deque lines: the .cif file lines.
    :rtype: ``deque``"""

    new_lines = deque()
    while lines:
        line = lines.popleft()
        if line.startswith(";"):
            string = [line[1:].strip()]
            while not lines[0].startswith(";"):
                string.append(lines.popleft())
            lines.popleft()
            new_lines[-1] += ' "{}"'.format(" ".join(string).replace('"', "\x1a").replace("'", "\x1b"))
        else:
            new_lines.append(line)
    return new_lines


def mmcif_lines_to_mmcif_blocks(lines):
    """A .cif file is ultimately a list of tables. This function takes a list of
    .cif file lines and splits them into these table blocks. Each block will be
    a ``dict`` containing a category name and a list of lines.

    :param deque lines: the .cif file lines.
    :rtype: ``list``"""

    category = None
    block, blocks = [], []
    while lines:
        line = lines.popleft()
        if line.startswith("data_"):
            continue
        if line.startswith("_"):
            line_category = line.split(".")[0]
            if line_category != category:
                if category:
                    blocks.append({"category": category[1:], "lines": block})
                category = line_category
                block = []
        if line.startswith("loop_"):
            if category:
                blocks.append({"category": category[1:], "lines": block})
            category = lines[0].split(".")[0]
            block = []
        block.append(line)
    if block:
        blocks.append({"category": category[1:], "lines": block})
    return blocks


def non_loop_block_to_list(block):
    """Takes a simple block ``dict`` with no loop and turns it into a table
    ``list``.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

    d = {}
    for index in range(len(block["lines"]) - 1):
        if block["lines"][index + 1][0] != "_":
            block["lines"][index] += " " + block["lines"][index + 1]
    block["lines"] = [l for l in block["lines"] if l[0] == "_"]
    for line in block["lines"]:
        name = line.split(".")[1].split()[0]
        value = line
        if line.startswith("_"):
            value = " ".join(line.split()[1:])
        d[name] = value
    return [d]


def loop_block_to_list(block):
    """Takes a loop block ``dict`` where the initial lines are table headers and
    turns it into a table ``list``. Sometimes a row is broken over several lines
    so this function deals with that too.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

    names, lines, _ = [], [], True
    body_start = 0
    for index, line in enumerate(block["lines"][1:], start=1):
        if not line.startswith("_" + block["category"]):
            body_start = index
            break
    names = [l.split(".")[1].rstrip() for l in block["lines"][1:body_start]]
    lines = [split_values(l) for l in block["lines"][body_start:]]
    l = []
    for n in range(len(lines) - 1):
        while n < len(lines) - 1 and len(lines[n]) + len(lines[n + 1]) <= len(names):
            lines[n] += lines[n + 1]
            lines.pop(n + 1)
    for line in lines:
        l.append({name: value for name, value in zip(names, line)})
    return l


def split_values(line):
    """The body of a .cif table is a series of lines, with each cell divided by
    whitespace. This function takes a string line and breaks it into cells.

    There are a few peculiarities to handle. Sometimes a cell is a string
    enclosed in quote marks, and spaces within this string obviously shouldn't
    be used to break the line. This function handles all of that.

    :param str line: the .cif line to split.
    :rtype: ``list``"""

    if not re.search("['\"]", line):
        return line.split()
    chars = deque(line.strip())
    values, value, in_string = [], [], False
    while chars:
        char = chars.popleft()
        if char == " " and not in_string:
            values.append(value)
            value = []
        elif char in "'\"":
            if in_string and chars and chars[0] != " ":
                value.append(char)
            else:
                in_string = not in_string
        else:
            value.append(char)
    values.append(value)
    return ["".join(v) for v in values if v]


def strip_quotes(mmcif_dict):
    """Goes through each table in the mmcif ``dict`` and removes any unneeded
    quote marks from the cells.

    :param dict mmcif_dict: the almost finished .mmcif dictionary to clean."""

    for name, table in mmcif_dict.items():
        for row in table:
            for k, value in row.items():
                for char in "'\"":
                    if value[0] == char and value[-1] == char:
                        row[k] = value[1:-1]
                    row[k] = row[k].replace("\x1a", '"').replace("\x1b", "'")


def mmcif_dict_to_data_dict(mmcif_dict) -> dict[str, Any]:
    """Converts an .mmcif dictionary into an atomium data dictionary, with the
    same standard layout that the other file formats get converted into.

    :param dict mmcif_dict: the .mmcif dictionary.
    :rtype: ``dict``"""

    data_dict = {
        "description": {"code": None, "title": None, "deposition_date": None, "classification": None, "keywords": [], "authors": []},
        "experiment": {"technique": None, "source_organism": None, "expression_system": None, "missing_residues": []},
        "quality": {"resolution": None, "rvalue": None, "rfree": None},
        "geometry": {"assemblies": [], "crystallography": {}},
        "models": [],
    }
    update_description_dict(mmcif_dict, data_dict)
    update_experiment_dict(mmcif_dict, data_dict)
    update_quality_dict(mmcif_dict, data_dict)
    update_geometry_dict(mmcif_dict, data_dict)
    update_models_list(mmcif_dict, data_dict)
    return data_dict


def update_description_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its description sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "code", "entry", "id")
    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "title", "struct", "title")
    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "deposition_date", "pdbx_database_status", "recvd_initial_deposition_date", date=True)
    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "classification", "struct_keywords", "pdbx_keywords")
    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "keywords", "struct_keywords", "text", split=True)
    mmcif_to_data_transfer(mmcif_dict, data_dict, "description", "authors", "audit_author", "name", multi=True)


def update_experiment_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its experiment sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict, "experiment", "technique", "exptl", "method")
    for cat, key in [
        ["entity_src_nat", "pdbx_organism_scientific"],
        ["entity_src_gen", "pdbx_gene_src_scientific_name"],
        ["pdbx_entity_src_syn", "organism_scientific"],
    ]:
        mmcif_to_data_transfer(mmcif_dict, data_dict, "experiment", "source_organism", cat, key)
        if data_dict["experiment"]["source_organism"] not in [None, "?"]:
            break
    mmcif_to_data_transfer(mmcif_dict, data_dict, "experiment", "expression_system", "entity_src_gen", "pdbx_host_org_scientific_name")
    for r in mmcif_dict.get("pdbx_unobs_or_zero_occ_residues", []):
        insert = "" if r["PDB_ins_code"] in "?." else r["PDB_ins_code"]
        data_dict["experiment"]["missing_residues"].append({"id": f"{r['auth_asym_id']}.{r['auth_seq_id']}{insert}", "name": r["auth_comp_id"]})


def update_quality_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its quality sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict, "quality", "resolution", "reflns", "d_resolution_high", func=float)
    if not data_dict["quality"]["resolution"]:
        mmcif_to_data_transfer(mmcif_dict, data_dict, "quality", "resolution", "refine", "ls_d_res_high", func=float)
    mmcif_to_data_transfer(mmcif_dict, data_dict, "quality", "rvalue", "refine", "ls_R_factor_R_work", func=float)
    if not data_dict["quality"]["rvalue"]:
        mmcif_to_data_transfer(mmcif_dict, data_dict, "quality", "rvalue", "refine", "ls_R_factor_obs", func=float)
    mmcif_to_data_transfer(mmcif_dict, data_dict, "quality", "rfree", "refine", "ls_R_factor_R_free", func=float)


def update_geometry_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its geometry sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    data_dict["geometry"]["assemblies"] = [
        {
            "id": int(a["id"]),
            "software": a.get("method_details", None),
            "delta_energy": None,
            "buried_surface_area": None,
            "surface_area": None,
            "transformations": [],
        }
        for a in mmcif_dict.get("pdbx_struct_assembly", [])
    ]
    operations = {
        o["id"]: [[float(o["matrix[{}][{}]".format(r, c)]) for c in [1, 2, 3]] + [float(o["vector[{}]".format(r)])] for r in [1, 2, 3]] + [[0, 0, 0, 1]]
        for o in mmcif_dict.get("pdbx_struct_oper_list", [])
    }
    for assembly in data_dict["geometry"]["assemblies"]:
        if assembly["software"] == "?":
            assembly["software"] = None
        assign_metrics_to_assembly(mmcif_dict, assembly)
        assign_transformations_to_assembly(mmcif_dict, operations, assembly)
    update_crystallography_dict(mmcif_dict, data_dict)


def assign_metrics_to_assembly(mmcif_dict, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant energy etc. information to update it with.

    :param dict mmcif_dict: The dictionary to read.
    :param dict assembly: The assembly to update."""

    for a in mmcif_dict.get("pdbx_struct_assembly_prop", []):
        if a["biol_id"] == str(assembly["id"]):
            if a["type"] == "MORE":
                assembly["delta_energy"] = float(a["value"].split("/")[0])
            elif a["type"] == "SSA (A^2)":
                assembly["surface_area"] = float(a["value"].split("/")[0])
            elif a["type"] == "ABSA (A^2)":
                assembly["buried_surface_area"] = float(a["value"].split("/")[0])


def assign_transformations_to_assembly(mmcif_dict, operations, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant transformation information to update it with.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict operations: the processed operations matrices.
    :param dict assembly: the assembly to update."""

    for gen in mmcif_dict.get("pdbx_struct_assembly_gen", []):
        if gen["assembly_id"] == str(assembly["id"]):
            op_ids_groups = get_operation_id_groups(gen["oper_expression"])
            ops = operation_id_groups_to_operations(operations, op_ids_groups)
            for operation in ops:
                assembly["transformations"].append(
                    {"chains": gen["asym_id_list"].split(","), "matrix": [row[:3] for row in operation[:3]], "vector": [row[-1] for row in operation[:3]]}
                )


def get_operation_id_groups(expression):
    """Takes an operator expression from an .mmcif transformation dict, and
    works out what transformation IDs it is referring to. For example, (1,2,3)
    becomes [[1, 2, 3]], (1-3)(8-11,17) becomes [[1, 2, 3], [8, 9, 10, 11, 17]],
    and so on.

    :param str expression: The expression to parse.
    :rtype: ``list``"""

    if expression[0] != "(":
        expression = "({})".format(expression)
    groups = re.findall(r"\((.+?)\)", expression)
    group_ids = []
    for group in groups:
        ids = []
        elements = group.split(",")
        for element in elements:
            if "-" in element:
                bounds = [int(x) for x in element.split("-")]
                ids += [str(n) for n in list(range(bounds[0], bounds[1] + 1))]
            else:
                ids.append(element)
        group_ids.append(ids)
    return group_ids


def update_crystallography_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its crystallography
    sub-sub-dictionary with information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    if mmcif_dict.get("cell"):
        mmcif_to_data_transfer(mmcif_dict, data_dict["geometry"], "crystallography", "space_group", "symmetry", "space_group_name_H-M")
        data_dict["geometry"]["crystallography"]["unit_cell"] = [
            float(mmcif_dict["cell"][0][key].replace("?", "0")) for key in ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"]
        ]
    if data_dict["geometry"]["crystallography"].get("space_group") == "NA":
        data_dict["geometry"]["crystallography"] = {}


def operation_id_groups_to_operations(operations, operation_id_groups):
    """Creates a list of operation matrices for an assembly, from a list of
    operation IDs - cross multiplying as required.

    :param dict operations: the parsed .mmcif operations.
    :param list operation_id_groups: the operation IDs."""

    operation_groups = [[operations[i] for i in ids] for ids in operation_id_groups]
    while len(operation_groups) and len(operation_groups) != 1:
        operations = []
        for op1 in operation_groups[0]:
            for op2 in operation_groups[1]:
                operations.append(np.matmul(op1, op2))
        operation_groups[0] = operations
        operation_groups.pop(1)
    return operation_groups[0]


def update_models_list(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its models list with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    data_dict["models"] = []
    types = {e["id"]: e["type"] for e in mmcif_dict.get("entity", {})}
    names = {e["id"]: e["name"] for e in mmcif_dict.get("chem_comp", {}) if e["mon_nstd_flag"] != "y"}
    entities = {m["id"]: m["entity_id"] for m in mmcif_dict.get("struct_asym", [])}
    # sequences = make_sequences(mmcif_dict)
    secondary_structure = make_secondary_structure(mmcif_dict)
    aniso = make_aniso(mmcif_dict)
    model = {"polymer": {}, "non-polymer": {}, "water": {}, "branched": {}}
    model_num = mmcif_dict["atom_site"][0]["pdbx_PDB_model_num"]
    for atom in mmcif_dict["atom_site"]:
        if atom["pdbx_PDB_model_num"] != model_num:
            data_dict["models"].append(model)
            model = {"polymer": {}, "non-polymer": {}, "water": {}, "branched": {}}
            model_num = atom["pdbx_PDB_model_num"]
        mol_type = types[entities[atom["label_asym_id"]]]
        if mol_type == "polymer" or mol_type == "branched":
            add_atom_to_polymer(atom, aniso, model, names)
        else:
            add_atom_to_non_polymer(atom, aniso, model, mol_type, names)
    data_dict["models"].append(model)
    for model in data_dict["models"]:
        add_sequences_to_polymers(model, mmcif_dict, entities)
        add_secondary_structure_to_polymers(model, secondary_structure)


def make_aniso(mmcif_dict):
    """Makes a mapping of atom IDs to anisotropy information.

    :param mmcif_dict: the .mmcif dict to read.
    :rtype: ``dict``"""

    return {
        int(a["id"]): [float(a["U[{}][{}]".format(x, y)]) for x, y in ["11", "22", "33", "12", "13", "23"]] for a in mmcif_dict.get("atom_site_anisotrop", [])
    }


def make_secondary_structure(mmcif_dict):
    """Creates a dictionary of helices and strands, with each having a list of
    start and end residues.

    :param mmcif_dict: the .mmcif dict to read.
    :rtype: ``dict``"""

    helices, strands = [], []
    for helix in mmcif_dict.get("struct_conf", []):
        helices.append(
            [
                "{}.{}{}".format(
                    helix[f"{x}_auth_asym_id"],
                    helix[f"{x}_auth_seq_id"],
                    helix[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
                )
                for x in ["beg", "end"]
            ]
        )
    for strand in mmcif_dict.get("struct_sheet_range", []):
        strands.append(
            [
                "{}.{}{}".format(
                    strand[f"{x}_auth_asym_id"],
                    strand[f"{x}_auth_seq_id"],
                    strand[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
                )
                for x in ["beg", "end"]
            ]
        )
    return {"helices": helices, "strands": strands}


def add_atom_to_polymer(atom, aniso, model, names):
    """Takes an MMCIF atom dictionary, converts it, and adds it to a polymer
    dictionary.

    :param dict atom: the .mmcif dictionary to read.
    :param dict aniso: lookup dictionary for anisotropy information.
    :param dict model: the model to update.
    :param dict names: the lookup dictionary for full name information."""

    mol_id = atom["auth_asym_id"]
    res_id = make_residue_id(atom)
    try:
        model["polymer"][mol_id]["residues"][res_id]["atoms"][int(atom["id"])] = atom_dict_to_atom_dict(atom, aniso)
    except Exception:
        name = atom["auth_comp_id"]
        try:
            model["polymer"][mol_id]["residues"][res_id] = {
                "name": name,
                "full_name": names.get(name),
                "atoms": {int(atom["id"]): atom_dict_to_atom_dict(atom, aniso)},
                "number": len(model["polymer"][mol_id]["residues"]) + 1,
            }
        except Exception:
            model["polymer"][mol_id] = {
                "internal_id": atom["label_asym_id"],
                "helices": [],
                "strands": [],
                "residues": {
                    res_id: {
                        "name": name,
                        "atoms": {int(atom["id"]): atom_dict_to_atom_dict(atom, aniso)},
                        "number": 1,
                        "full_name": names.get(name),
                    }
                },
            }


def add_atom_to_non_polymer(atom, aniso, model, mol_type, names):
    """Takes an MMCIF atom dictionary, converts it, and adds it to a non-polymer
    dictionary.

    :param dict atom: the .mmcif dictionary to read.
    :param dict aniso: lookup dictionary for anisotropy information.
    :param dict model: the model to update.
    :param str mol_type: non-polymer or water.
    :param dict names: the lookup dictionary for full name information."""

    mol_id = make_residue_id(atom)
    try:
        model[mol_type][mol_id]["atoms"][int(atom["id"])] = atom_dict_to_atom_dict(atom, aniso)
    except Exception:
        name = atom["auth_comp_id"]
        model[mol_type][mol_id] = {
            "name": name,
            "full_name": names.get(name),
            "internal_id": atom["label_asym_id"],
            "polymer": atom["auth_asym_id"],
            "atoms": {int(atom["id"]): atom_dict_to_atom_dict(atom, aniso)},
        }


def make_residue_id(d):
    """Generates a residue ID for an atom.

    :param dict d: the atom dictionary to read.
    :rtype: ``str``"""

    insert = "" if d["pdbx_PDB_ins_code"] in "?." else d["pdbx_PDB_ins_code"]
    return "{}.{}{}".format(d["auth_asym_id"], d["auth_seq_id"], insert)


def add_sequences_to_polymers(model, mmcif_dict, entities):
    """Takes a pre-populated mapping of chain IDs to entity IDs, and uses them
    to add sequence information to a model.

    :param dict model: the model to update.
    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict entities: a mapping of chain IDs to entity IDs."""

    sequences = make_sequences(mmcif_dict)
    for polymer in model["polymer"].values():
        polymer["sequence"] = sequences.get(entities.get(polymer["internal_id"], ""), "")


def add_secondary_structure_to_polymers(model, ss_dict) -> None:
    """Updates polymer dictionaries with secondary structure information, from
    a previously created mapping.

    :param dict model: the model to update.
    :param dict ss_dict: the mapping to read."""

    for ss in ("helices", "strands"):
        for segment in ss_dict[ss]:
            chain = model["polymer"].get(segment[0].split(".")[0])
            if chain:
                in_segment = False
                chain[ss].append([])
                for residue_id in chain["residues"].keys():
                    if residue_id == segment[0]:
                        in_segment = True
                    if in_segment:
                        chain[ss][-1].append(residue_id)
                    if residue_id == segment[1]:
                        break


def make_sequences(mmcif_dict):
    """Creates a mapping of entity IDs to sequences.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :rtype: ``dict``"""

    return {
        e["id"]: "".join([CODES.get(res["mon_id"], "X") for res in mmcif_dict.get("entity_poly_seq", []) if res["entity_id"] == e["id"]])
        for e in mmcif_dict.get("entity", [])
        if e["type"] == "polymer"
    }


def atom_dict_to_atom_dict(d, aniso_dict):
    """Turns an .mmcif atom dictionary into an atomium atom data dictionary.

    :param dict d: the .mmcif atom dictionary.
    :param dict d: the mapping of atom IDs to anisotropy.
    :rtype: ``dict``"""

    charge = "pdbx_formal_charge"
    atom = {
        "x": d["Cartn_x"],
        "y": d["Cartn_y"],
        "z": d["Cartn_z"],
        "element": d["type_symbol"],
        "name": d.get("label_atom_id"),
        "occupancy": d.get("occupancy", 1),
        "bvalue": d.get("B_iso_or_equiv"),
        "charge": d.get(charge, 0) if d.get(charge) != "?" else 0,
        "alt_loc": d.get("label_alt_id") if d.get("label_alt_id") != "." else None,
        "anisotropy": aniso_dict.get(int(d["id"]), [0, 0, 0, 0, 0, 0]),
        "is_hetatm": d.get("group_PDB", "ATOM") == "HETATM",
    }
    for key in ["x", "y", "z", "charge", "bvalue", "occupancy"]:
        if atom[key] is not None:
            atom[key] = float(atom[key])
    return atom


def mmcif_to_data_transfer(mmcif_dict, data_dict, d_cat, d_key, m_table, m_key, date=False, split=False, multi=False, func=None):
    """A function for transfering a bit of data from a .mmcif dictionary to a
    data dictionary, or doing nothing if the data doesn't exist.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update.
    :param str d_cat: the top-level key in the data dictionary.
    :param str d_key: the data dictionary field to update.
    :param str m_table: the name of the .mmcif table to look in.
    :param str m_key: the .mmcif field to read.
    :param bool date: if True, the value will be converted to a date.
    :param bool split: if True, the value will be split on commas.
    :param bool multi: if True, every row in the table will be read.
    :param function func: if given, this will be applied to the value."""

    try:
        if multi:
            value = [row[m_key] for row in mmcif_dict[m_table]]
        else:
            value = mmcif_dict[m_table][0][m_key]
        if date:
            value = datetime.strptime(value, "%Y-%m-%d").date()
        if split:
            value = value.replace(", ", ",").split(",")
        if func:
            value = func(value)
        data_dict[d_cat][d_key] = None if value == "?" else value
    except Exception:
        pass


def get_atom_name(atom):
    """Formats an atom name for packing in .cif.

    :param Atom atom: the atom to read.
    :rtype: ``str``"""

    return '"{}"'.format(atom._name) if "'" in atom._name else atom._name


def split_residue_id(atom):
    """Takes an atom and splits its het ID into components.

    :param Atom atom: the atom to read.
    :rtype: ``tuple``"""

    if atom.het:
        id = atom.het.id.split(".")[-1]
        num = "".join([c for c in id if c.isdigit()])
        insert = "".join([c for c in id if c.isalpha()]) or "?"
        return num, insert
    return ".."
