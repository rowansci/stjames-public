import re
from pathlib import Path
from typing import Annotated, Any, Iterable, Optional, Self, TypeAlias

import pydantic
from pydantic import AfterValidator, NonNegativeInt, PositiveInt, ValidationError
from rdkit import Chem
from rdkit.Chem import AllChem

from .atom import Atom
from .base import Base, round_float, round_optional_float
from .periodic_cell import PeriodicCell
from .types import (
    FloatPerAtom,
    Matrix3x3,
    Vector3D,
    Vector3DPerAtom,
    round_optional_float_per_atom,
    round_optional_matrix3x3,
    round_optional_vector3d,
    round_optional_vector3d_per_atom,
    round_vector3d_per_atom,
)

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol


class MoleculeReadError(RuntimeError):
    pass


class VibrationalMode(Base):
    frequency: Annotated[float, AfterValidator(round_float(3))]  # in cm-1
    reduced_mass: Annotated[float, AfterValidator(round_float(3))]  # amu
    force_constant: Annotated[float, AfterValidator(round_float(3))]  # mDyne/Å
    displacements: Annotated[Vector3DPerAtom, AfterValidator(round_vector3d_per_atom(6))]  # Å
    ir_intensity: Annotated[Optional[float], AfterValidator(round_optional_float(3))] = None  # km/mol


class Molecule(Base):
    charge: int
    multiplicity: PositiveInt
    atoms: list[Atom]

    # for periodic boundary conditions
    cell: Optional[PeriodicCell] = None

    energy: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None  # in Hartree
    scf_iterations: Optional[NonNegativeInt] = None
    scf_completed: Optional[bool] = None
    elapsed: Annotated[Optional[float], AfterValidator(round_optional_float(3))] = None  # in seconds

    homo_lumo_gap: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None  # in eV

    gradient: Annotated[Optional[Vector3DPerAtom], AfterValidator(round_optional_vector3d_per_atom(6))] = None  # Hartree/Å
    stress: Annotated[Optional[Matrix3x3], AfterValidator(round_optional_matrix3x3(6))] = None  # Hartree/Å

    velocities: Annotated[Optional[Vector3DPerAtom], AfterValidator(round_optional_vector3d_per_atom(6))] = None  # Å/fs

    mulliken_charges: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    mulliken_spin_densities: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    dipole: Annotated[Optional[Vector3D], AfterValidator(round_optional_vector3d(6))] = None  # in Debye

    vibrational_modes: Optional[list[VibrationalMode]] = None

    zero_point_energy: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None
    thermal_energy_corr: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None
    thermal_enthalpy_corr: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None
    thermal_free_energy_corr: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None

    smiles: Optional[str] = None
    calculation_index: int | None = None

    def __len__(self) -> int:
        return len(self.atoms)

    def distance(self, atom1: PositiveInt, atom2: PositiveInt) -> float:
        r"""
        Get the distance between atoms.

        >>> mol = Molecule.from_xyz("H 0 1 0\nH 0 0 1")
        >>> mol.distance(1, 2)
        1.4142135623730951
        """
        return sum((q2 - q1) ** 2 for q1, q2 in zip(self.atoms[atom1 - 1].position, self.atoms[atom2 - 1].position)) ** 0.5  # type: ignore [no-any-return,unused-ignore]

    @property
    def coordinates(self) -> Vector3DPerAtom:
        return [a.position for a in self.atoms]

    def translated(self, vector: Vector3D) -> Self:
        r"""
        Translate the molecule by a vector.

        >>> mol = Molecule.from_xyz("H 0 0 0\nH 0 0 1")
        >>> print(mol.translated((1, 0, 0)).to_xyz())
        2
        <BLANKLINE>
        H     1.0000000000    0.0000000000    0.0000000000
        H     1.0000000000    0.0000000000    1.0000000000
        """

        def translated(position: Vector3D) -> Vector3D:
            return tuple(q + v for q, v in zip(position, vector, strict=True))  # type: ignore [return-value]

        atoms = [atom.copy(update={"position": translated(atom.position)}) for atom in self.atoms]

        return self.copy(update={"atoms": atoms})

    @property
    def atomic_numbers(self) -> list[NonNegativeInt]:
        return [a.atomic_number for a in self.atoms]

    @property
    def sum_energy_zpe(self) -> Optional[float]:
        if (self.energy is None) or (self.zero_point_energy is None):
            return None
        return self.energy + self.zero_point_energy

    @property
    def sum_energy_thermal_corr(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_energy_corr is None):
            return None
        return self.energy + self.thermal_energy_corr

    @property
    def sum_energy_enthalpy(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_enthalpy_corr is None):
            return None
        return self.energy + self.thermal_enthalpy_corr

    @property
    def sum_energy_free_energy(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_free_energy_corr is None):
            return None
        return self.energy + self.thermal_free_energy_corr

    @pydantic.model_validator(mode="after")
    def check_electron_sanity(self) -> Self:
        num_electrons = sum(self.atomic_numbers) - self.charge
        num_unpaired_electrons = self.multiplicity - 1
        if (num_electrons - num_unpaired_electrons) % 2 != 0:
            raise ValueError(
                f"The combination of {num_electrons} electrons, charge {self.charge}, and multiplicity {self.multiplicity} is impossible. "
                "Double-check the charge and multiplicity values given and verify that they are correct."
            )

        return self

    @classmethod
    def from_file(cls: type[Self], filename: Path | str, format: str | None = None, charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        r"""
        Read a molecule from a file.

        >>> import tempfile
        >>> with tempfile.NamedTemporaryFile("w+", suffix=".xyz") as f:
        ...    _ = f.write("2\nComment\nH 0 0 0\nF 0 0 1")
        ...    _ = f.seek(0)
        ...    mol = Molecule.from_file(f.name)
        >>> print(mol.to_xyz())
        2
        <BLANKLINE>
        H     0.0000000000    0.0000000000    0.0000000000
        F     0.0000000000    0.0000000000    1.0000000000
        """
        filename = Path(filename)
        if not format:
            format = filename.suffix[1:]

        with open(filename) as f:
            match format:
                case "xyz":
                    return cls.from_xyz_lines(f.readlines(), charge=charge, multiplicity=multiplicity)
                case "extxyz":
                    return cls.from_extxyz_lines(f.readlines(), charge=charge, multiplicity=multiplicity)
                case _:
                    raise ValueError(f"Unsupported {format=}")

    @classmethod
    def from_xyz(cls: type[Self], xyz: str, charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        r"""
        Generate a Molecule from an XYZ string.

        Note: only supports single molecule inputs.

        >>> len(Molecule.from_xyz("2\nComment\nH 0 0 0\nH 0 0 1"))
        2
        """
        return cls.from_xyz_lines(xyz.strip().splitlines(), charge=charge, multiplicity=multiplicity)

    @classmethod
    def from_xyz_lines(cls: type[Self], lines: Iterable[str], charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        lines = list(lines)
        if len(lines[0].split()) == 1:
            natoms = lines[0].strip()
            if not natoms.isdigit() or (int(lines[0]) != len(lines) - 2):
                raise MoleculeReadError(f"First line of XYZ file should be the number of atoms, got: {lines[0]} != {len(lines) - 2}")
            lines = lines[2:]

        try:
            return cls(atoms=[Atom.from_xyz(line) for line in lines], charge=charge, multiplicity=multiplicity)
        except (ValueError, ValidationError) as e:
            raise MoleculeReadError("Error reading molecule from xyz") from e

    def to_xyz(self, comment: str = "", out_file: Path | str | None = None) -> str:
        r"""
        Generate an XYZ string.

        >>> mol = Molecule.from_xyz("2\nComment\nH 0 1 2\nF 1 2 3")
        >>> print(mol.to_xyz(comment="HF"))
        2
        HF
        H     0.0000000000    1.0000000000    2.0000000000
        F     1.0000000000    2.0000000000    3.0000000000
        >>> import tempfile
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     file = Path(directory) / "mol.xyz"
        ...     out = mol.to_xyz(comment="HF", out_file=file)
        ...     with file.open() as f:
        ...         Molecule.from_xyz(f.read()).to_xyz("HF") == out
        True
        """
        geom = "\n".join(map(str, self.atoms))
        out = f"{len(self)}\n{comment}\n{geom}"

        if out_file:
            with Path(out_file).open("w") as f:
                f.write(out)

        return out

    @classmethod
    def from_extxyz(cls: type[Self], extxyz: str, charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        r"""
        Generate a Molecule from a EXTXYZ string. Currently only supporting Lattice and Properties fields.

        >>> Molecule.from_extxyz('''
        ... 2
        ... Lattice="6.0 0.0 0.0 6.0 0.0 0.0 6.0 0.0 0.0"Properties=species:S:1:pos:R:3
        ... H 0 0 0
        ... H 0 0 1
        ... ''').cell.lattice_vectors
        ((6.0, 0.0, 0.0), (6.0, 0.0, 0.0), (6.0, 0.0, 0.0))
        """

        return cls.from_extxyz_lines(extxyz.strip().splitlines(), charge=charge, multiplicity=multiplicity)

    @classmethod
    def from_extxyz_lines(cls: type[Self], lines: Iterable[str], charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        """
        Parses an EXTXYZ file, extracting atom positions, forces (if present), and metadata.

        Supports:
        - Lattice vectors (cell information)
        - Properties field (species, positions, forces, etc.)
        - Other metadata like charge, multiplicity, energy, etc.
        """

        lines = list(lines)

        # Ensure first line contains number of atoms
        if len(lines[0].split()) == 1:
            natoms = lines[0].strip()
            if not natoms.isdigit() or (int(natoms) != len(lines) - 2):
                raise MoleculeReadError(f"First line should be number of atoms, got: {lines[0]} != {len(lines) - 2}")
            lines = lines[1:]
        else:
            raise MoleculeReadError(f"First line should be an integer denoting atom count. Got {lines[0].split()}")

        # Extract metadata from the second line
        metadata = parse_comment_line(lines[0])

        lines = lines[1:]  # Move to atomic data

        # Extract key metadata values
        cell = metadata.get("cell")
        charge = metadata.get("total_charge", charge)
        multiplicity = metadata.get("multiplicity", multiplicity)
        energy = metadata.get("energy", None)

        # Extract properties like species, positions, and forces
        properties_str = metadata.get("Properties", "")
        properties = properties_str.split(":")

        if len(properties) < 4 or properties[0].lower() != "species":
            raise MoleculeReadError(f"Invalid or missing 'Properties' field in EXTXYZ. Got: {properties_str}")

        # Identify column indices for position and force data
        pos_idx = None
        force_idx = None
        current_idx = 0  # Start after 'species:S'

        while current_idx < len(properties):
            if properties[current_idx].lower() == "pos" and properties[current_idx + 1].lower() == "r" and properties[current_idx + 2] == "3":
                pos_idx = current_idx
            elif properties[current_idx].lower() == "forces" and properties[current_idx + 1].lower() == "r" and properties[current_idx + 2] == "3":
                force_idx = current_idx
            current_idx += 3  # Move to next property group

        if pos_idx is None:
            raise MoleculeReadError("No position data ('pos:R:3') found in Properties field.")

        atoms = []
        gradients: list[tuple[float, float, float]] = []  # Will store gradients (-forces)

        for line in lines:
            values = line.split()
            species = values[0]
            position = tuple(map(float, values[1:4]))

            atom = Atom.from_xyz(f"{species} {position[0]} {position[1]} {position[2]}")
            atoms.append(atom)

            if force_idx is not None:
                force = tuple(map(float, values[4:7]))
                gradients.append(tuple(-f for f in force))  # type:ignore [arg-type]

        actual_gradients = gradients or None

        return cls(atoms=atoms, cell=cell, charge=charge, multiplicity=multiplicity, energy=energy, gradient=actual_gradients)

    @classmethod
    def from_rdkit(cls: type[Self], rdkm: RdkitMol, cid: int = 0) -> Self:
        if len(rdkm.GetConformers()) == 0:
            rdkm = _embed_rdkit_mol(rdkm)

        atoms = []
        atomic_numbers = [atom.GetAtomicNum() for atom in rdkm.GetAtoms()]  # type: ignore [no-untyped-call, unused-ignore]
        geom = rdkm.GetConformers()[cid].GetPositions()

        for i in range(len(atomic_numbers)):
            atoms.append(Atom(atomic_number=atomic_numbers[i], position=geom[i]))

        charge = Chem.GetFormalCharge(rdkm)
        multiplicity = 1

        return cls(atoms=atoms, charge=charge, multiplicity=multiplicity)

    @classmethod
    def from_smiles(cls: type[Self], smiles: str) -> Self:
        rdkm = Chem.MolFromSmiles(smiles)
        assert rdkm is not None
        return cls.from_rdkit(rdkm)


def _embed_rdkit_mol(rdkm: RdkitMol) -> RdkitMol:
    try:
        AllChem.SanitizeMol(rdkm)  # type: ignore [attr-defined, unused-ignore]
    except Exception as e:
        raise ValueError("Molecule could not be generated -- invalid chemistry!\n") from e

    rdkm = AllChem.AddHs(rdkm)  # type: ignore [attr-defined, unused-ignore]
    try:
        status1 = AllChem.EmbedMolecule(rdkm, maxAttempts=200)  # type: ignore [attr-defined, unused-ignore]
        assert status1 >= 0
    except Exception as e:
        status1 = AllChem.EmbedMolecule(rdkm, maxAttempts=200, useRandomCoords=True)  # type: ignore [attr-defined, unused-ignore]
        if status1 < 0:
            raise ValueError(f"Cannot embed molecule! Error: {e}")

    AllChem.MMFFOptimizeMolecule(rdkm, maxIters=200)  # type: ignore [attr-defined, call-arg, unused-ignore]

    return rdkm


def parse_comment_line(line: str) -> dict[str, Any]:
    """
    Parse the comment line of an EXTXYZ file, extracting lattice, properties, and metadata.

    Supports:
    - Lattice vectors (cell information)
    - Properties field (species, positions, forces, etc.)
    - Other metadata fields like charge, multiplicity, energy, etc.

    Returns:
        dict: Parsed properties including `cell`, `Properties`, and any other metadata.
    """

    # Regular expression to match key="value", key='value', or key=value
    pattern = r"(\S+?=(?:\".*?\"|\'.*?\'|\S+))"
    pairs = re.findall(pattern, line)

    prop_dict = {}
    cell = None

    for pair in pairs:
        key, value = pair.split("=", 1)
        key = key.lower().strip()
        value = value.strip("'\"")

        if key == "lattice":
            lattice_values = value.split()
            if len(lattice_values) != 9:
                raise MoleculeReadError(f"Lattice should have 9 entries, got {len(lattice_values)}")

            # Convert to a 3x3 tuple of floats
            try:
                cell = tuple(tuple(map(float, lattice_values[i : i + 3])) for i in range(0, 9, 3))
            except ValueError:
                raise MoleculeReadError(f"Lattice should be floats, got {lattice_values}")

        elif key == "properties":
            prop_dict["Properties"] = value  # Store properties as raw text

        elif key in {"total_charge", "multiplicity"}:
            prop_dict[key] = int(value)  # Store charge and multiplicity as integers

        elif key in {"energy"}:
            prop_dict[key] = float(value)  # Store energy as float

        else:
            prop_dict[key] = value  # Store any other metadata as a string

    if cell is not None:
        prop_dict["cell"] = PeriodicCell(lattice_vectors=cell)

    return prop_dict
