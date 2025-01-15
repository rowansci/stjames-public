from datetime import date
from pathlib import Path
from typing import Any, Literal

import atomium  # type: ignore [import-untyped]
from atomium.pdb import pdb_dict_to_data_dict, pdb_string_to_pdb_dict  # type: ignore [import-untyped]
from pydantic import BaseModel, ConfigDict, Field

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
    name: str
    charge: float
    occupancy: float
    alt_loc: str | None
    anisotropy: list[float]
    bvalue: float
    is_hetatm: bool


class PDBWater(BaseModel):
    """A water molecule."""

    model_config = ConfigDict(extra=EXTRA)

    name: str | None
    full_name: str | None
    atoms: dict[int, PDBAtom] = {}
    internal_id: str | None
    polymer: str


class PDBResidue(BaseModel):
    """A structure."""

    model_config = ConfigDict(extra=EXTRA)

    name: str | None
    full_name: str | None = None
    atoms: dict[int, PDBAtom] = {}
    number: int


class PDBPolymer(BaseModel):
    """A polymer chain."""

    model_config = ConfigDict(extra=EXTRA)

    internal_id: str
    helices: list[list[str]] = []
    residues: dict[str, PDBResidue] = {}
    sequence: str | None
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
    non_polymer: dict[str, PDBNonPolymer] = Field(alias="non-polymer", default_factory=dict)
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
    software: str | None
    buried_surface_area: float | None
    surface_area: float | None
    delta_energy: float | None
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

    expression_system: str | None
    missing_residues: list[PDBMissingResidue] = []
    source_organism: str | None
    technique: str | None


class PDBDescription(BaseModel):
    """A description of the molecule."""

    model_config = ConfigDict(extra=EXTRA)

    code: str | None
    title: str | None
    authors: list[str] = []
    classification: str | None
    deposition_date: date | None
    keywords: list[str] = []


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
    return PDB.model_validate(atomium.open(str(path), data_dict=True))


def fetch_pdb(code: str) -> PDB:
    """Fetch a pdb from the Protein Data Bank."""
    return PDB.model_validate(atomium.fetch(code, data_dict=True))


def pdb_from_string(pdb: str) -> PDB:
    """Read a PDB from a string."""
    return PDB.model_validate(pdb_dict_to_data_dict(pdb_string_to_pdb_dict(pdb)))
