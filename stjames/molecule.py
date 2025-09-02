import re
from pathlib import Path
from typing import Annotated, Any, Iterable, Optional, Self, Sequence, TypeAlias, TypedDict, TypeVar

import numpy as np
import pydantic
from pydantic import AfterValidator, NonNegativeInt, PositiveInt, ValidationError
from rdkit import Chem
from rdkit.Chem import AllChem

from .atom import Atom
from .base import Base, round_float, round_optional_float
from .data import SYMBOL_ELEMENT
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

    def distance(self, i: PositiveInt, j: PositiveInt) -> float:
        r"""
        Calculate the distance between atoms.

        >>> mol = Molecule.from_xyz("H 0 1 0\nH 0 0 1")
        >>> mol.distance(1, 2)
        1.4142135623730951
        """
        return sum((q2 - q1) ** 2 for q1, q2 in zip(self.atoms[i - 1].position, self.atoms[j - 1].position)) ** 0.5  # type: ignore [no-any-return,unused-ignore]

    def angle(self, i: PositiveInt, j: PositiveInt, k: PositiveInt, degrees: bool = True) -> float:
        r"""
        Calculate the angle between three atoms.

        >>> Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1").angle(1, 2, 3)
        90.0
        """

        return angle(self.coordinates[i - 1], self.coordinates[j - 1], self.coordinates[k - 1], degrees=degrees)

    def dihedral(self, i: int, j: int, k: int, l: int, degrees: bool = True, positive_domain: bool = True) -> float:
        r"""
        Calculate the dihedral angle between four atoms.

        >>> Molecule.from_xyz("H 0 0 0\nO 0 0 1\nO 0 1 1\nH 1 1 1").dihedral(1, 2, 3, 4)
        270.0
        """
        return dihedral(
            self.coordinates[i - 1],
            self.coordinates[j - 1],
            self.coordinates[k - 1],
            self.coordinates[l - 1],
            degrees=degrees,
            positive_domain=positive_domain,
        )

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
    def enthalpy(self) -> Optional[float]:
        return self.sum_energy_enthalpy

    @property
    def sum_energy_free_energy(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_free_energy_corr is None):
            return None
        return self.energy + self.thermal_free_energy_corr

    @property
    def gibbs_free_energy(self) -> Optional[float]:
        return self.sum_energy_free_energy

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
    def from_extxyz_lines(
        cls: type[Self],
        lines: Iterable[str],
        charge: int | None = None,
        multiplicity: PositiveInt | None = None,
        cell: PeriodicCell | None = None,
    ) -> Self:
        """
        Parses an EXTXYZ file, extracting atom positions, forces (if present), and metadata.

        Supports:
        - Lattice vectors (cell information)
        - Properties field (species, positions, forces, etc.)
        - Other metadata like charge, multiplicity, energy, etc.

        :param lines: Iterable of lines from an EXTXYZ file
        :param charge: total charge of the molecule (default: 0 if not found)
        :param multiplicity: spin multiplicity of the molecule (default: 1 if not found)
        :param cell: PeriodicCell containing lattice vectors
        :return: Molecule
        :raises MoleculeReadError: if the file is not in the correct format
        """
        if not isinstance(lines, Sequence):
            lines = list(lines)

        # Ensure first line contains number of atoms
        if len(lines[0].split()) == 1:
            natoms = lines[0].strip()
            if not natoms.isdigit() or (int(natoms) != len(lines) - 2):
                raise MoleculeReadError(f"First line should be number of atoms, got: {lines[0]} != {len(lines) - 2}")
            data_line, *lines = lines[1:]
        else:
            raise MoleculeReadError(f"First line should be an integer denoting atom count. Got {lines[0].split()}")

        metadata = parse_extxyz_comment_line(data_line)

        T = TypeVar("T")

        def metadata_optional_get(key: str, value: T | None, default: T) -> T:
            """Set key to default if not found in metadata"""
            if value is None:
                return metadata.get(key, default)  # type: ignore [return-value]

            return value

        charge = metadata_optional_get("total_charge", charge, 0)
        multiplicity = metadata_optional_get("multiplicity", multiplicity, 1)
        cell = cell or metadata.get("cell")
        energy = metadata.get("energy", None)

        force_idx = None
        if properties := metadata.get("properties", "").split(":"):
            if properties[0].lower() != "species":
                raise MoleculeReadError(f"Invalid or missing 'Properties' field in EXTXYZ, got: {properties}")

            # Identify column indices for position and force data
            pos_idx = None
            current_idx = 0  # Start after 'species:S'

            while current_idx < len(properties):
                if properties[current_idx].lower() == "pos" and properties[current_idx + 1].lower() == "r" and properties[current_idx + 2] == "3":
                    pos_idx = current_idx
                elif properties[current_idx].lower() == "forces" and properties[current_idx + 1].lower() == "r" and properties[current_idx + 2] == "3":
                    force_idx = current_idx
                current_idx += 3

            if pos_idx is None:
                raise MoleculeReadError("No position data ('pos:R:3') found in Properties field.")

        def parse_line_atoms(line: str) -> Atom:
            symbol, sx, sy, sz, *_ = line.split()
            atomic_number = SYMBOL_ELEMENT[symbol.title()]
            x, y, z = map(float, (sx, sy, sz))

            return Atom(atomic_number=atomic_number, position=(x, y, z))

        def parse_line_with_grad(line: str) -> tuple[Atom, Vector3D]:
            symbol, sx, sy, sz, sgx, sgy, sgz, *_ = line.split()
            atomic_number = SYMBOL_ELEMENT[symbol.title()]
            x, y, z = map(float, (sx, sy, sz))
            gx, gy, gz = map(float, (sgx, sgy, sgz))

            return (
                Atom(atomic_number=atomic_number, position=(x, y, z)),
                (-gx, -gy, -gz),
            )

        atoms: list[Atom]
        gradients: list[Vector3D] | None
        if force_idx is not None:
            atoms, gradients = zip(*map(parse_line_with_grad, lines), strict=True)  # type: ignore [assignment]
        else:
            atoms = [parse_line_atoms(line) for line in lines]
            gradients = None

        return cls(atoms=atoms, cell=cell, charge=charge, multiplicity=multiplicity, energy=energy, gradient=gradients)

    @classmethod
    def from_rdkit(cls: type[Self], rdkm: RdkitMol, cid: int = 0, multiplicity: int = 1) -> Self:
        if len(rdkm.GetConformers()) == 0:
            rdkm = _embed_rdkit_mol(rdkm)

        atomic_numbers = [atom.GetAtomicNum() for atom in rdkm.GetAtoms()]  # type: ignore [no-untyped-call, unused-ignore]
        atoms = [
            Atom(atomic_number=atom, position=xyz)  # keep open
            for atom, xyz in zip(atomic_numbers, rdkm.GetConformers()[cid].GetPositions(), strict=True)
        ]

        charge = Chem.GetFormalCharge(rdkm)

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


class EXTXYZMetadata(TypedDict, total=False):
    properties: Any
    total_charge: int
    multiplicity: int
    energy: float
    cell: PeriodicCell


def parse_extxyz_comment_line(line: str) -> EXTXYZMetadata:
    """
    Parse the comment line of an EXTXYZ file, extracting lattice, properties, and metadata.

    Supports:
    - Lattice vectors (cell information)
    - Properties field (species, positions, forces, etc.)
    - Other metadata fields like charge, multiplicity, energy, etc.

    :param line: comment line from an EXTXYZ file
    :return: parsed properties

    >>> parse_extxyz_comment_line('Lattice="6.0 0.0 0.0 6.0 0.0 0.0 6.0 0.0 0.0"Properties=species:S:1:pos:R:3')
    {'cell': PeriodicCell(lattice_vectors=((6.0, 0.0, 0.0), (6.0, 0.0, 0.0), (6.0, 0.0, 0.0)), is_periodic=(True, True, True), volume=0.0), 'properties': 'species:S:1:pos:R:3'}
    """  # noqa: E501

    # Regular expression to match key="value", key='value', or key=value
    pattern = r"(\S+?=(?:\".*?\"|\'.*?\'|\S+))"
    pairs = re.findall(pattern, line)

    prop_dict: EXTXYZMetadata = {}
    for pair in pairs:
        key, value = pair.split("=", 1)
        key = key.lower().strip()
        value = value.strip("'\"")

        if key == "lattice":
            lattice_values = value.split()
            if len(lattice_values) != 9:
                raise MoleculeReadError(f"Lattice should have 9 entries, got {len(lattice_values)}")

            try:
                cell = tuple(tuple(map(float, lattice_values[i : i + 3])) for i in range(0, 9, 3))
            except ValueError:
                raise MoleculeReadError(f"Lattice should be floats, got {lattice_values}")

            prop_dict["cell"] = PeriodicCell(lattice_vectors=cell)

        elif key == "properties":
            prop_dict["properties"] = value

        elif key == "total_charge":
            prop_dict["total_charge"] = int(value)
        elif key == "multiplicity":
            prop_dict["multiplicity"] = int(value)
        elif key == "energy":
            prop_dict["energy"] = float(value)
        else:
            prop_dict[key] = value  # type: ignore [literal-required]

    return prop_dict


def angle(p0: Sequence[float], p1: Sequence[float], p2: Sequence[float], degrees: bool = True) -> float:
    """
    Angle between three points.

    :param i, j, k: positions of points
    :param degrees: whether to return in degrees
    :return: angle in radians or degrees
    """
    a0, a1, a2 = map(np.asarray, (p0, p1, p2))
    u = a1 - a0
    v = a1 - a2

    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    cos_theta = np.dot(u, v) / (nu * nv)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    ang = np.arccos(cos_theta)
    if degrees:
        return float(np.degrees(ang))
    return float(ang)


def dihedral(p0: Sequence[float], p1: Sequence[float], p2: Sequence[float], p3: Sequence[float], degrees: bool = True, positive_domain: bool = True) -> float:
    """
    Dihedral angle between four points.

    :param p0, p1, p2, p3: points
    :param degrees: whether to return in degrees
    :param positive_domain: (0, 360] if True else (-180, 180]
    :return: angle in degrees or radians (or nan if collinearities detected)

    >>> a = [0, 0, 0]
    >>> b = [0, 0, 1]
    >>> c = [0, 1, 1]
    >>> d1 = [0, 1, 2]
    >>> d2 = [0.5, 1, 1.5]
    >>> dihedral(a, b, c, d1)
    180.0
    >>> dihedral(a, b, c, d2, positive_domain=False)
    -135.0
    """
    a0, a1, a2, a3 = map(np.asarray, (p0, p1, p2, p3))
    b0 = a1 - a0
    b1 = a2 - a1
    b2 = a3 - a2

    b1 = b1 / np.linalg.norm(b1)

    v = b1 * np.dot(b0, b1) - b0
    w = b2 - b1 * np.dot(b2, b1)

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    ang = np.arctan2(y, x)

    if positive_domain and ang < 0:
        ang += 2 * np.pi

    if degrees:
        return float(np.degrees(ang))
    return float(ang)
