import re
from datetime import date, datetime
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

import stjames.atomium_stjames as astj
from stjames.atomium_stjames.mmcif import mmcif_dict_to_data_dict, mmcif_string_to_mmcif_dict
from stjames.atomium_stjames.pdb import inverse_make_sequences, pdb_dict_to_data_dict, pdb_string_to_pdb_dict
from stjames.types import Matrix3x3, Vector3D

# Mostly for testing purposes
EXTRA: Literal["allow", "ignore", "forbid"] = "allow"


class PDBAtom(BaseModel):
    """An atom within a residue."""

    model_config = ConfigDict(extra=EXTRA)

    x: float
    y: float
    z: float
    element: str
    name: str | None = None
    charge: float | None = None
    occupancy: float | None = None
    alt_loc: str | None = None
    anisotropy: list[float] | None = None
    bvalue: float
    is_hetatm: bool | None = None


class PDBWater(BaseModel):
    """A water molecule."""

    model_config = ConfigDict(extra=EXTRA)

    name: str | None = None
    full_name: str | None = None
    atoms: dict[int, PDBAtom] = {}
    internal_id: str | None = None
    polymer: str


class PDBResidue(BaseModel):
    """A structure."""

    model_config = ConfigDict(extra=EXTRA)

    name: str | None = None
    full_name: str | None = None
    atoms: dict[int, PDBAtom] = {}
    number: int


class PDBPolymer(BaseModel):
    """A polymer chain."""

    model_config = ConfigDict(extra=EXTRA)

    internal_id: str
    helices: list[list[str]] = []
    residues: dict[str, PDBResidue] = {}
    sequence: str | None = None
    strands: list[list[str]] = []


class PDBNonPolymer(BaseModel):
    """Non-polymeric molecules/atoms (e.g. ions and ligands)."""

    model_config = ConfigDict(extra=EXTRA)

    name: str
    full_name: str | None = None
    atoms: dict[int, PDBAtom]
    internal_id: str
    polymer: str


class PDBModel(BaseModel):
    """Structure data."""

    model_config = ConfigDict(extra=EXTRA)

    polymer: dict[str, PDBPolymer] = {}
    non_polymer: dict[str, PDBNonPolymer] = {}
    branched: dict[str, Any] = {}
    water: dict[str, PDBWater] = {}


class PDBTransformations(BaseModel):
    """Transformations applied to the structure."""

    model_config = ConfigDict(extra=EXTRA)

    chains: list[str]
    matrix: Matrix3x3
    vector: Vector3D


class PDBAssembly(BaseModel):
    """How the structure was assembled."""

    model_config = ConfigDict(extra=EXTRA)

    transformations: list[PDBTransformations]
    software: str | None = None
    buried_surface_area: float | None = None
    surface_area: float | None = None
    delta_energy: float | None = None
    id: int


class PDBCrystallography(BaseModel):
    """Crystallography related information."""

    model_config = ConfigDict(extra=EXTRA)

    space_group: str | None = None
    unit_cell: list[float] | None = None


class PDBGeometry(BaseModel):
    """Details of the geometry."""

    model_config = ConfigDict(extra=EXTRA)

    assemblies: list[PDBAssembly] = []
    crystallography: PDBCrystallography = Field(default_factory=PDBCrystallography)


class PDBQuality(BaseModel):
    """Quality metrics."""

    model_config = ConfigDict(extra=EXTRA)

    resolution: float | None = None
    rfree: float | None = None
    rvalue: float | None = None


class PDBMissingResidue(BaseModel):
    model_config = ConfigDict(extra=EXTRA)

    name: str
    id: str


class PDBExperiment(BaseModel):
    """Details of the experiment."""

    model_config = ConfigDict(extra=EXTRA)

    expression_system: str | None = None
    missing_residues: list[PDBMissingResidue] = []
    source_organism: str | None = None
    technique: str | None = None


class PDBDescription(BaseModel):
    """A description of the molecule."""

    model_config = ConfigDict(extra=EXTRA)

    code: str | None = None
    title: str | None = None
    authors: list[str] = []
    classification: str | None = None
    deposition_date: str | None = None
    keywords: list[str] = []

    @field_validator("deposition_date", mode="before")
    @classmethod
    def date_to_string(cls, v: str | date | None) -> str | None:
        if v is None:
            return v

        if isinstance(v, date):
            return v.isoformat()

        return str(v)


class PDB(BaseModel):
    """A PDB formatted file."""

    model_config = ConfigDict(extra=EXTRA)

    description: PDBDescription
    experiment: PDBExperiment
    geometry: PDBGeometry
    models: list[PDBModel] = []
    quality: PDBQuality = Field(default_factory=PDBQuality)


def read_pdb(path: Path | str) -> PDB:
    """Read a pdb located at path."""
    return PDB.model_validate(astj.open(str(path), data_dict=True))


def fetch_pdb(code: str) -> PDB:
    """Fetch a pdb from the Protein Data Bank."""
    return PDB.model_validate(astj.fetch(code, data_dict=True))


def fetch_pdb_from_mmcif(code: str) -> PDB:
    """Fetch a pdb from the Protein Data Bank."""
    code += ".cif"
    return PDB.model_validate(astj.fetch(code, data_dict=True))


def pdb_from_pdb_filestring(pdb: str) -> PDB:
    """Read a PDB from a string."""
    return PDB.model_validate(pdb_dict_to_data_dict(pdb_string_to_pdb_dict(pdb)))


def pdb_from_mmcif_filestring(pdb: str) -> PDB:
    """Read a PDB from a string."""
    return PDB.model_validate(mmcif_dict_to_data_dict(mmcif_string_to_mmcif_dict(pdb)))


def pdb_object_to_pdb_filestring(
    pdb: PDB,
    header: bool = False,
    source: bool = False,
    keyword: bool = False,
    seqres: bool = True,
    hetnam: bool = True,
    remark: bool = False,
    crystallography: bool = False,
) -> str:
    pdb_lines: list[str] = []
    chains: list[str] = []

    if header:
        pdb_lines.extend(_build_header_section(pdb))

    if source:
        pdb_lines.extend(_build_source_section(pdb))

    if keyword:
        pdb_lines.extend(_build_keyword_section(pdb))

    full_name_dict: dict[str, str] = {}
    seqres_lines, chains = _build_secondary_structure_and_seqres(pdb, full_name_dict)

    if seqres:
        pdb_lines.extend(seqres_lines)

    if hetnam:
        pdb_lines.extend(_build_hetname_section(full_name_dict))

    if remark:
        pdb_lines.extend(_build_remark_section(pdb, chains))

    if crystallography:
        pdb_lines.extend(_build_crystallography_section(pdb))

    for model_index, model in enumerate(pdb.models, start=1):
        # If more than one model, add MODEL line
        if len(pdb.models) > 1:
            pdb_lines.append(f"MODEL     {model_index:>4}")

        # === 1) Polymers (protein, DNA, etc.) ===
        for chain_id, polymer in model.polymer.items():
            # Use polymer's internal_id if you want that as the chain ID
            # otherwise just use the dictionary key
            this_chain_id = polymer.internal_id or chain_id

            for _residue_id, residue in polymer.residues.items():
                assert residue.name is not None
                for _atom_id, atom in residue.atoms.items():
                    line = _format_atom_line(
                        serial=_atom_id,
                        atom=atom,
                        chain_id=this_chain_id,
                        res_name=residue.name,
                        res_num=_residue_id[2:],
                        alt_loc=atom.alt_loc or "",
                    )
                    pdb_lines.append(line)
                    if atom.anisotropy and atom.anisotropy != [0, 0, 0, 0, 0, 0]:
                        line = _format_anisou_line(
                            serial=_atom_id,
                            atom=atom,
                            chain_id=this_chain_id,
                            res_name=residue.name,
                            res_num=_residue_id[2:],
                            alt_loc=atom.alt_loc or "",
                        )
                        pdb_lines.append(line)

            pdb_lines.append(f"TER   {_atom_id + 1:>5}      {residue.name:>3} {this_chain_id}{_residue_id[2:]:>4}")

        # === 2) Non-polymers (e.g. ligands, ions) ===
        for _np_id, nonpoly in model.non_polymer.items():
            # We'll treat each non-polymer as if it had a chain ID = nonpoly.polymer (or "Z")
            chain_id_for_np = nonpoly.polymer or "Z"

            # For residue name, we can just use nonpoly.name or a 3-letter variant
            # There's no standard "residue number" for these, so pick something
            # or let the user define it in the original model. We'll just use 1 for demonstration.
            # If you prefer incremental numbering, keep a separate counter.
            for _atom_id, atom in nonpoly.atoms.items():
                line = _format_atom_line(
                    serial=_atom_id,
                    atom=atom,
                    chain_id=chain_id_for_np,
                    res_name=nonpoly.name,
                    res_num=_np_id[2:],
                )
                pdb_lines.append(line)
                if atom.anisotropy and atom.anisotropy != [0, 0, 0, 0, 0, 0]:
                    line = _format_anisou_line(
                        serial=_atom_id,
                        atom=atom,
                        chain_id=chain_id_for_np,
                        res_name=nonpoly.name,
                        res_num=_np_id[2:],
                    )
                    pdb_lines.append(line)

        # === 3) Water ===
        for _w_id, water in model.water.items():
            # Water is typically "HOH" in PDB
            for _atom_id, atom in water.atoms.items():
                line = _format_atom_line(
                    serial=_atom_id,
                    atom=atom,
                    chain_id=_w_id[0],  # Or you can use water.polymer if set
                    res_name="HOH",
                    res_num=_w_id[2:],  # or an incrementing value
                )
                pdb_lines.append(line)
                if atom.anisotropy and atom.anisotropy != [0, 0, 0, 0, 0, 0]:
                    line = _format_anisou_line(
                        serial=_atom_id,
                        atom=atom,
                        chain_id=_w_id[0],
                        res_name="HOH",
                        res_num=_w_id[2:],
                    )
                    pdb_lines.append(line)

        # === 4) Branched ===
        # If your structure has branched molecules (glycans, etc.),
        # adapt similarly. For now, let's demonstrate if there's anything in branched
        for _b_id, branched_obj in model.branched.items():
            # "branched_obj" could be a custom structure. We'll assume it
            # mirrors the format of non_polymer or something similar.
            # If it has `.atoms`, we do the same:
            if isinstance(branched_obj, dict) and "atoms" in branched_obj:
                for _atom_id, atom in branched_obj["atoms"].items():
                    line = _format_atom_line(
                        serial=_atom_id,
                        atom=atom,
                        chain_id="B",
                        res_name="BRN",  # or branched_obj.get("name", "BRN")
                        res_num="1",
                    )
                    pdb_lines.append(line)
                    if atom.anisotropy and atom.anisotropy != [0, 0, 0, 0, 0, 0]:
                        line = _format_anisou_line(
                            serial=_atom_id,
                            atom=atom,
                            chain_id="B",
                            res_name="BRN",
                            res_num="1",
                        )
                        pdb_lines.append(line)

        if len(pdb.models) > 1:
            pdb_lines.append("ENDMDL")

    # Finally, the PDB standard ends with an END record
    pdb_lines.append("END")

    resulting_string = _create_filestring(pdb_lines)
    return resulting_string


def _create_filestring(lines: list[str]) -> str:
    # Join the lines with newline characters and add a newline at the end if desired
    filestring = "\n".join(lines) + "\n"
    return filestring


def _format_date(date_str: str | None) -> str | None:
    """
    Formats a date string from "YYYY-MM-DD" to "DD-MMM-YY".

    Args:
        date_str (str): Date string in "YYYY-MM-DD" format.

    Returns:
        str: Formatted date string in "DD-MMM-YY" format.
    """
    if date_str is None:
        return None
    date_obj = datetime.strptime(date_str, "%Y-%m-%d").date()
    return date_obj.strftime("%d-%b-%y").upper()


def _format_atom_line(
    serial: int,
    atom: PDBAtom,
    chain_id: str,
    res_name: str,
    res_num: str | None,
    alt_loc: str = "",
) -> str:
    """
    Return a single PDB ATOM/HETATM record line as a string, using standard
    column alignment conventions:

    See https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf for details
    """
    record_type = "HETATM" if atom.is_hetatm else "ATOM  "

    # Columns are typically strict. We'll use Python formatting with fixed widths.
    # Some fields might need defaults if missing.
    alt_loc_char = alt_loc if alt_loc else " "
    residue_name = (res_name or "UNK")[:3]  # limit to 3 chars
    chain_char = (chain_id or "A")[:1]  # PDB chain ID is 1 char
    residue_num_str = "1"
    insertion_code = " "
    if res_num:
        match = re.match(r"(\d+)([a-zA-Z]*)", res_num)
        if match:
            residue_num_str, insertion_code = match.groups()
            insertion_code = insertion_code if insertion_code != "" else " "

    residue_num = int(residue_num_str)

    # Format charge: PDB uses e.g. " 2-", " 1+" in columns 79-80
    # If your model stores charges differently, adapt as needed.
    # For simplicity, let's store integer/float charges as strings, e.g. " 0", " 2", etc.
    # Or we can leave it blank if zero.
    chg = ""
    if atom.charge and abs(atom.charge) > 0:
        # E.g., +1.0 -> " +1", -2.0 -> " -2"
        # Convert to integer if it's always integral
        chg_val = int(atom.charge) if float(atom.charge).is_integer() else atom.charge
        chg = f"{chg_val:2}"
    else:
        chg = "  "

    atom_name = atom.name if atom.name else atom.element
    occupancy = atom.occupancy if atom.occupancy else 1.0

    # Construct the line.
    # Use exact spacing & field widths to match PDB guidelines.
    line = (
        f"{record_type}"
        f"{serial:5d} "  # atom serial number (columns 7-11)
        f"{atom_name:<4}"  # atom name (columns 13-16, left-justified in this snippet)
        f"{alt_loc_char}"  # altLoc (column 17)
        f"{residue_name:>3}"  # residue name (columns 18-20)
        f" {chain_char}"  # chain ID (column 22)
        f"{residue_num:4d}"  # residue sequence number (columns 23-26)
        f"{insertion_code}"
        f"   "  # columns 27-30 (spacing)
        f"{atom.x:8.3f}"  # x (columns 31-38)
        f"{atom.y:8.3f}"  # y (columns 39-46)
        f"{atom.z:8.3f}"  # z (columns 47-54)
        f"{occupancy:6.2f}"  # occupancy (columns 55-60)
        f"{atom.bvalue:6.2f}"  # temp factor (columns 61-66)
        f"          "  # columns 67-76 (padding)
        f"{atom.element:>2}"  # element (columns 77-78)
        f"{chg:>2}"  # charge (columns 79-80)
    )
    return line


def _format_anisou_line(
    serial: int,
    atom: PDBAtom,
    chain_id: str,
    res_name: str,
    res_num: str | None,
    alt_loc: str = "",
) -> str:
    """
    Return a single PDB ANISOU record line as a string, using standard
    column alignment conventions:

    See https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf for details
    """
    record_type = "ANISOU"

    # Columns are typically strict. We'll use Python formatting with fixed widths.
    # Some fields might need defaults if missing.
    alt_loc_char = alt_loc if alt_loc else " "
    residue_name = (res_name or "UNK")[:3]  # limit to 3 chars
    chain_char = (chain_id or "A")[:1]  # PDB chain ID is 1 char
    residue_num_str = "1"
    insertion_code = " "
    if res_num:
        match = re.match(r"(\d+)([a-zA-Z]*)", res_num)
        if match:
            residue_num_str, insertion_code = match.groups()
            insertion_code = insertion_code if insertion_code != "" else " "

    residue_num = int(residue_num_str)

    chg = ""
    if atom.charge and abs(atom.charge) > 0:
        # E.g., +1.0 -> " +1", -2.0 -> " -2"
        # Convert to integer if it's always integral
        chg_val = int(atom.charge) if float(atom.charge).is_integer() else atom.charge
        chg = f"{chg_val:2}"
    else:
        chg = "  "

    atom_name = atom.name if atom.name else atom.element

    if atom.anisotropy:
        aniso_lines = (
            f"{_float_to_pdb_string(atom.anisotropy[0]):>7}"  # x (columns 29-35)
            f"{_float_to_pdb_string(atom.anisotropy[1]):>7}"  # x (columns 36-42)
            f"{_float_to_pdb_string(atom.anisotropy[2]):>7}"  # x (columns 43-49)
            f"{_float_to_pdb_string(atom.anisotropy[3]):>7}"  # x (columns 50-56)
            f"{_float_to_pdb_string(atom.anisotropy[4]):>7}"  # x (columns 57-63)
            f"{_float_to_pdb_string(atom.anisotropy[5]):>7}"
        )
    else:
        space = " "
        aniso_lines = (
            f"{space:>7}"  # x (columns 29-35)
            f"{space:>7}"  # x (columns 36-42)
            f"{space:>7}"  # x (columns 43-49)
            f"{space:>7}"  # x (columns 50-56)
            f"{space:>7}"  # x (columns 57-63)
            f"{space:>7}"
        )

    # Construct the line.
    # Use exact spacing & field widths to match PDB guidelines.
    line = (
        f"{record_type}"
        f"{serial:5d} "  # atom serial number (columns 7-11)
        f"{atom_name:<4}"  # atom name (columns 13-16, left-justified in this snippet)
        f"{alt_loc_char}"  # altLoc (column 17)
        f"{residue_name:>3}"  # residue name (columns 18-20)
        f" {chain_char}"  # chain ID (column 22)
        f"{residue_num:4d}"  # residue sequence number (columns 23-26)
        f"{insertion_code}"
        f" "  # columns 27-28 (plus spacing)
        f"{aniso_lines}"
        f"      "  # columns 70-76 (padding)
        f"{atom.element:>2}"  # element (columns 77-78)
        f"{chg:>2}"  # charge (columns 79-80)
    )
    return line


# chat code
def _float_to_pdb_string(x: float) -> str:
    # Determine the sign
    sign = "-" if x < 0 else ""
    a = abs(x)

    if a < 1:
        # Format with exactly 4 decimals, e.g. 0.0044 -> "0.0044"
        s = f"{a:.4f}"
        # Remove the "0." and then remove any leading zeros.
        significant = s[2:].lstrip("0")
        return sign + significant
    else:
        # Format with exactly 4 decimals. For example, 1.131 -> "1.1310"
        s = f"{a:.4f}"
        # Split into integer and fractional parts.
        integer_part, fractional_part = s.split(".")
        # We want a total of 5 digits. So, the number of fractional digits we need is:
        needed = 5 - len(integer_part)
        # Use the needed number of digits from the fractional part.
        result = integer_part + fractional_part[:needed]
        return sign + result


def _helix_list_to_pdb_helix(polymer_dict: dict[str, PDBPolymer], helices: list[list[str]]) -> list[str]:
    helix_lines = []
    for i, helix in enumerate(helices, start=1):
        start_aa_name = polymer_dict[helix[0][0]].residues[helix[0]].name
        end_aa_name = polymer_dict[helix[-1][0]].residues[helix[-1]].name
        helix_line = f"HELIX  {i:>3} {i:>3} {start_aa_name} {helix[0][0]} {helix[0][2:]:>4}  {end_aa_name} {helix[-1][0]} {helix[-1][2:]:>4}  1{len(helix):>36}"
        helix_lines.append(helix_line)
    return helix_lines


def _strand_list_to_pdb_sheets(polymer_dict: dict[str, PDBPolymer], strands: list[list[str]]) -> list[str]:
    strand_lines = []
    for i, strand in enumerate(strands, start=1):
        start_aa_name = polymer_dict[strand[0][0]].residues[strand[0]].name
        end_aa_name = polymer_dict[strand[-1][0]].residues[strand[-1]].name
        helix_line = (
            f"SHEET  {i:>3} {strand[0][0]:>3}{len(strands):>2} {start_aa_name} {strand[0][0]}{strand[0][2:]:>4}  "
            f"{end_aa_name} {strand[-1][0]}{strand[-1][2:]:>4} {-1 if i != 1 else 0:>2}"
        )
        strand_lines.append(helix_line)
    return strand_lines


def _build_header_section(pdb: PDB) -> list[str]:
    header = f"HEADER    {pdb.description.classification or '':<40}{_format_date(pdb.description.deposition_date) or '':<10}  {pdb.description.code or '':<5}"
    title = f"TITLE     {pdb.description.title or '':<70}"
    exp_dta = f"EXPDTA    {pdb.experiment.technique or '':<69}"
    authors = f"AUTHOR    {','.join(pdb.description.authors).upper():<69}"

    return [header, title, exp_dta, authors]


def _build_source_section(pdb: PDB) -> list[str]:
    """Builds the source organism and expression system lines."""
    organism_line = f"SOURCE    ORGANISM_SCIENTIFIC: {(pdb.experiment.source_organism + ';') if pdb.experiment.source_organism else '':<69}"
    expression_line = f"SOURCE    EXPRESSION_SYSTEM: {(pdb.experiment.expression_system + ';') if pdb.experiment.expression_system else '':<69}"
    return [organism_line, expression_line]


def _build_keyword_section(pdb: PDB) -> list[str]:
    """Builds the keyword (KEYWDS) lines."""
    lines = []
    for i, keyword in enumerate(pdb.description.keywords):
        if i == len(pdb.description.keywords) - 1:
            lines.append(f"KEYWDS    {keyword:<79}")
        else:
            lines.append(f"KEYWDS    {keyword + ',':<79}")
    return lines


def _build_secondary_structure_and_seqres(pdb: PDB, full_name_dict: dict[str, str]) -> tuple[list[str], list[str]]:
    """
    Iterates over models and polymers to build secondary structure lines (e.g. sheets, helices)
    and sequence records (SEQRES). Also collects full names for heterogen records.
    Returns a tuple: (list of seqres (and secondary structure) lines, list of chain IDs).
    """
    seqres_lines = []
    chains = []

    for model in pdb.models:
        for chain_id, polymer in model.polymer.items():
            chains.append(chain_id)
            # Add sheet and helix records (if available)
            for strand_line in _strand_list_to_pdb_sheets(model.polymer, polymer.strands):
                seqres_lines.append(strand_line)
            for helix_line in _helix_list_to_pdb_helix(model.polymer, polymer.helices):
                seqres_lines.append(helix_line)
            # Add SEQRES lines from the polymer’s sequence
            if polymer.sequence:
                seqres_lines.extend(inverse_make_sequences(polymer.sequence, chain_id))
            # Collect full names from each residue
            for _, residue in polymer.residues.items():
                if residue.full_name and residue.name:
                    full_name_dict[residue.name] = residue.full_name
        # Also collect full names for non-polymer molecules
        for _, non_polymer in model.non_polymer.items():
            if non_polymer.full_name and non_polymer.name:
                full_name_dict[non_polymer.name] = non_polymer.full_name

    return seqres_lines, chains


def _build_hetname_section(full_name_dict: dict[str, str]) -> list[str]:
    """Builds the HETNAM lines for non-polymer molecules."""
    lines = []
    for name, full_name in full_name_dict.items():
        if len(full_name) > 55:
            for i in range(0, len(full_name), 55):
                lines.append(f"HETNAM  {int(i / 55):>2} {name:<3} {full_name[i : i + 55]:<55}")
        else:
            lines.append(f"HETNAM     {name:<3} {full_name:<55}")
    return lines


def _build_remark_section(pdb: PDB, chains: list[str]) -> list[str]:
    """Builds REMARK lines (resolution, R factors, biomolecule and missing residues)."""
    lines = []
    lines.append(f"REMARK   2 RESOLUTION. {pdb.quality.resolution or '':>7} ANGSTROMS.")
    if pdb.quality.rfree:
        lines.append(f"REMARK   3   FREE R VALUE                     : {pdb.quality.rfree or ''}")
    if pdb.quality.rvalue:
        lines.append(f"REMARK   3   R VALUE            (WORKING SET) : {pdb.quality.rvalue or ''}")

    # REMARK 350: Biomolecule details
    lines.append("REMARK 350")
    lines.append("REMARK 350 COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN")
    lines.append("REMARK 350 BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE")
    lines.append("REMARK 350 MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS")
    lines.append("REMARK 350 GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND")
    lines.append("REMARK 350 CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN.")
    lines.append("REMARK 350")
    lines.append("REMARK 350 BIOMOLECULE: 1")
    lines.append("REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC")
    lines.append(f"REMARK 350 APPLY THE FOLLOWING TO CHAINS: {', '.join(chains)}")
    lines.append("REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000")
    lines.append("REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000")
    lines.append("REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000")

    # REMARK 465: Missing residues
    lines.append("REMARK 465 MISSING RESIDUES")
    lines.append("REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE")
    lines.append("REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN")
    lines.append("REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)")
    lines.append("REMARK 465")
    lines.append("REMARK 465   M RES C SSSEQI")
    for missing_residue in pdb.experiment.missing_residues:
        lines.append(f"REMARK 465     {missing_residue.name} {missing_residue.id[0]}   {missing_residue.id[2:]}")
    return lines


def _build_crystallography_section(pdb: PDB) -> list[str]:
    """Builds the CRYST1 line if unit cell data is provided."""
    lines = []
    if pdb.geometry.crystallography.unit_cell:
        lines.append(
            f"CRYST1{pdb.geometry.crystallography.unit_cell[0]:>9}"
            f"{pdb.geometry.crystallography.unit_cell[1]:>9}"
            f"{pdb.geometry.crystallography.unit_cell[2]:>9}"
            f"{pdb.geometry.crystallography.unit_cell[3]:>7}"
            f"{pdb.geometry.crystallography.unit_cell[4]:>7}"
            f"{pdb.geometry.crystallography.unit_cell[5]:>7} "
            f"{pdb.geometry.crystallography.space_group or '':<11}"
        )
    return lines
