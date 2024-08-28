"""Bond Dissociation Energy (BDE) workflow."""

from typing import Any, Self

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
            self.atoms = self.atoms + atomic_number_indices(self.initial_molecule, 1)
        if self.all_CX:
            X = {9, 17, 35, 53, 85, 117}
            self.atoms = self.atoms + atomic_number_indices(self.initial_molecule, X)

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
