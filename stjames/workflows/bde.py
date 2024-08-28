"""Bond Dissociation Energy (BDE) workflow."""

import itertools
from typing import Any, Iterable, Self

from pydantic import BaseModel, Field, PositiveInt, field_validator, model_validator

from ..mode import Mode
from ..molecule import Molecule
from ..types import UUID
from .multistage_opt import MultiStageOptMixin
from .workflow import Workflow

# the id of a mutable object may change, thus using object()
_sentinel_mso_mode = object()


class BDE(BaseModel):
    """
    Bond Dissociation Energy (BDE) result.

    energy => (E_{fragment1} + E_{fragment2}) - E_{starting molecule}

    :param fragment_idxs: indices of the atoms in the fragment that was dissociated (1-indexed)
    :param energy: BDE in kcal/mol
    :param fragment1_energy: energy of fragment 1
    :param fragment2_energy: energy of fragment 2
    :param calculations: calculation UUIDs
    """

    fragment_idxs: tuple[PositiveInt, ...]
    energy: float
    fragment1_energy: float
    fragment2_energy: float
    calculations: tuple[UUID, ...]

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """
        Return a string representation of the BDE result.

        >>> BDE(fragment_idxs=(1, 2), energy=1.0, fragment1_energy=50.0, fragment2_energy=50.0, calculations=[])
        <BDE (1, 2)  1.00>
        """
        return f"<{type(self).__name__} {self.fragment_idxs} {self.energy:>5.2f}>"


class BDEWorkflow(Workflow, MultiStageOptMixin):
    """
    Bond Dissociation Energy (BDE) workflow.

    Uses the modes from MultiStageOptSettings to compute BDEs.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Overridden:
    :param mso_mode: Mode for MultiStageOptSettings

    New:
    :param mode: Mode for workflow
    :param optimize_fragments: whether to optimize the fragments, or just the starting molecule (default depends on mode)
    :param atoms: atoms to dissociate (1-indexed)
    :param fragment_indices: fragments to dissociate (all fields feed into this, 1-indexed)
    :param all_CH: dissociate all C–H bonds
    :param all_CX: dissociate all C–X bonds (X ∈ {F, Cl, Br, I, At, Ts})
    :param opt_molecule: optimized starting molecule
    :param bdes: BDE results
    """

    mode: Mode
    mso_mode: Mode = _sentinel_mso_mode  # type: ignore [assignment]
    optimize_fragments: bool = None  # type: ignore [assignment]

    atoms: tuple[PositiveInt, ...] = Field(default_factory=tuple)
    fragment_indices: tuple[tuple[PositiveInt, ...], ...] = Field(default_factory=tuple)

    all_CH: bool = False
    all_CX: bool = False

    # Results
    opt_molecule: Molecule | None = None
    bdes: tuple[BDE] = Field(default_factory=tuple)

    def __str__(self) -> str:
        r"""
        Return a string representation of the BDE workflow.

        >>> print(BDEWorkflow(initial_molecule=Molecule.from_xyz("H 0 0 0\nF 0 0 1"), mode=Mode.METICULOUS, atoms=[1, 2]))
        BDEWorkflow METICULOUS
        (1,)
        (2,)
        """
        return f"{type(self).__name__} {self.mode.name}\n" + "\n".join(map(str, self.fragment_indices))

    def __repr__(self) -> str:
        """
        Return a string representation of the BDE workflow.

        >>> BDEWorkflow(initial_molecule=Molecule.from_xyz("He 0 0 0"), mode=Mode.METICULOUS, atoms=[])
        <BDEWorkflow METICULOUS>
        """
        return f"<{type(self).__name__} {self.mode.name}>"

    @property
    def energies(self) -> tuple[float, ...]:
        return tuple(bde.energy for bde in self.bdes)

    @field_validator("initial_molecule", mode="before")
    @classmethod
    def no_charge_or_spin(cls, mol: Molecule) -> Molecule:
        """Ensure the molecule has no charge or spin."""
        if mol.charge != 0 or mol.multiplicity != 1:
            raise ValueError("Charge and spin partitioning undefined for BDE, only neutral singlet molecules supported.")

        return mol

    @model_validator(mode="before")
    @classmethod
    def set_mso_mode(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Set the MultiStageOptSettings mode to match current SpinStates mode."""
        values["mso_mode"] = values["mode"]
        return values

    @model_validator(mode="after")
    def validate_and_build(self) -> Self:
        """Validate atomic numbers and build the atoms field."""
        self.atoms = tuple(self.atoms)
        self.fragment_indices = tuple(map(tuple, self.fragment_indices))

        match self.mode:
            case Mode.RECKLESS | Mode.RAPID:
                # Default off
                self.optimize_fragments = self.optimize_fragments or False
            case Mode.CAREFUL | Mode.METICULOUS:
                # Default on
                self.optimize_fragments = self.optimize_fragments or self.optimize_fragments is None
            case _:
                raise NotImplementedError(f"{self.mode} not implemented.")

        for atom in self.atoms:
            if atom > len(self.initial_molecule):
                raise ValueError(f"{atom=} is out of range.")
        for fragment_idxs in self.fragment_indices:
            if any(a > len(self.initial_molecule) for a in fragment_idxs):
                raise ValueError(f"{fragment_idxs=} contains atoms that are out of range.")

        if self.all_CH:
            self.atoms = self.atoms + tuple(H for _C, H in find_CH_bonds(self.initial_molecule))
        if self.all_CX:
            self.atoms = self.atoms + tuple(X for _C, X in find_CX_bonds(self.initial_molecule))

        # Combine atoms and fragments, remove duplicates, and sort
        self.fragment_indices = self.fragment_indices + tuple((a,) for a in self.atoms)
        assert isinstance(self.fragment_indices, tuple)
        self.fragment_indices = tuple(sorted(set(self.fragment_indices)))

        return self


def atomic_number_indices(molecule: Molecule, atomic_numbers: set[PositiveInt] | PositiveInt) -> tuple[PositiveInt, ...]:
    r"""
    Return the indices of the atoms with the given atomic numbers.

    :param molecule: Molecule of interest
    :param atomic_numbers: atomic number(s) of interest

    >>> H2O = Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")
    >>> atomic_number_indices(H2O, 1)
    (1, 3)
    >>> atomic_number_indices(H2O, {8, 1})
    (1, 2, 3)
    >>> atomic_number_indices(H2O, {6, 9})
    ()
    """
    if isinstance(atomic_numbers, int):
        return tuple(i for i, an in enumerate(molecule.atomic_numbers, start=1) if an == atomic_numbers)
    return tuple(i for i, an in enumerate(molecule.atomic_numbers, start=1) if an in atomic_numbers)


def find_CH_bonds(molecule: Molecule, distance_max: float = 1.2) -> Iterable[tuple[PositiveInt, PositiveInt]]:
    r"""
    Find all C–H bonds in the molecule.

    :param molecule: Molecule of interest
    :param distance_max: distance max for bond

    >>> H2O = Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")
    >>> CH4 = Molecule.from_xyz("C 0 0 0\nH 0 0 1\nH 0 1 0\nH 1 0 0\nH 0 0 -1")
    >>> ethane = Molecule.from_xyz("C -1 0 0\nC 1 0 0\nH -1 0 1\nH -1 0 -1\nH -1 1 0\nH 1 0 1\nH 1 0 -1\nH 1 1 0")
    >>> list(find_CH_bonds(H2O))
    []
    >>> list(find_CH_bonds(CH4))
    [(1, 2), (1, 3), (1, 4), (1, 5)]
    >>> list(find_CH_bonds(ethane))
    [(1, 3), (1, 4), (1, 5), (2, 6), (2, 7), (2, 8)]
    """
    yield from find_AB_bonds(molecule, 6, 1, distance_max)


def find_CX_bonds(molecule: Molecule) -> Iterable[tuple[PositiveInt, PositiveInt]]:
    r"""
    Find all C–X bonds in the molecule.

    :param molecule: Molecule of interest
    :param distance_max: distance max for bond

    >>> HCF = Molecule.from_xyz("H 0 0 0\nC 0 0 1\nF 0 1 1")
    >>> list(find_CX_bonds(HCF))
    [(2, 3)]
    """
    halogens = {9: 2.0, 17: 2.2, 35: 2.5, 53: 2.8, 85: 3.0, 117: 4.0}
    yield from itertools.chain.from_iterable(find_AB_bonds(molecule, 6, x, distance) for x, distance in halogens.items())


def find_AB_bonds(molecule: Molecule, a: int, b: int, distance_max: float) -> Iterable[tuple[PositiveInt, PositiveInt]]:
    r"""
    Find all A–B bonds in the molecule.

    :param molecule: Molecule of interest
    :param a: atomic number of atom A
    :param b: atomic number of atom B
    :param distance_max: distance max for bond

    >>> H2O = Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")
    >>> list(find_AB_bonds(H2O, 8, 1, 1.1))
    [(2, 1), (2, 3)]
    """

    def close(idxs: tuple[PositiveInt, PositiveInt]) -> bool:
        return molecule.distance(*idxs) < distance_max

    yield from filter(
        close,
        itertools.product(
            atomic_number_indices(molecule, a),
            atomic_number_indices(molecule, b),
        ),
    )