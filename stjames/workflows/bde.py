"""Bond Dissociation Energy (BDE) workflow."""

from typing import Any, Self, Sequence

from pydantic import BaseModel, Field, PositiveInt, model_validator
from pydantic import NonNegativeInt as NNInt

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

    :param fragments: two fragments after dissociation
    :param energy: BDE in kcal/mol
    """

    fragments: tuple[Molecule, Molecule]
    energy: float
    fragment1_energy: float
    fragment2_energy: float
    calculation: list[UUID]

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.energy:>5.2f}>"


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
    :param atoms: atoms to dissociate
    :param fragments: fragments to dissociate (all fields feed into this)
    :param all_CH: dissociate all C–H bonds
    :param all_CX: dissociate all C–X bonds (X ∈ {F, Cl, Br, I, At, Ts})
    :param opt_molecule: optimized starting molecule
    :param bdes: list of BDE results
    """

    mode: Mode
    mso_mode: Mode = _sentinel_mso_mode  # type: ignore [assignment]

    atoms: Sequence[NNInt] = Field(default_factory=tuple)
    fragments: Sequence[Sequence[NNInt]] = Field(default_factory=tuple)

    all_CH: bool = False
    all_CX: bool = False

    # Results
    opt_molecule: Molecule | None = None
    bdes: list[BDE] = Field(default_factory=list)

    @property
    def energies(self) -> list[float]:
        return [bde.energy for bde in self.bdes]

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
        self.fragments = tuple(map(tuple, self.fragments))

        for atom in self.atoms:
            if atom > len(self.initial_molecule):
                raise ValueError(f"{atom=} is out of range.")
        for fragment in self.fragments:
            if any(a > len(self.initial_molecule) for a in fragment):
                raise ValueError(f"{fragment=} contains atoms that are out of range.")

        if self.all_CH:
            self.atoms = self.atoms + atomic_number_indices(self.initial_molecule, 1)
        if self.all_CX:
            X = {9, 17, 35, 53, 85, 117}
            self.atoms = self.atoms + atomic_number_indices(self.initial_molecule, X)

        # Combine atoms and fragments, remove duplicates, and sort
        self.fragments = self.fragments + tuple((a,) for a in self.atoms)
        assert isinstance(self.fragments, tuple)
        self.fragments = tuple(sorted(set(self.fragments)))

        return self


def atomic_number_indices(molecule: Molecule, atomic_numbers: set[PositiveInt] | PositiveInt) -> tuple[PositiveInt, ...]:
    """
    Return the indices of the atoms with the given atomic numbers.

    :param molecule: Molecule of interest
    :param atomic_numbers: atomic number(s) of interest

    >>> from stjames.molecule import Atom, Molecule
    >>> H2O = Molecule(
    ...     charge=0, multiplicity=1,
    ...     atoms=[
    ...         Atom(atomic_number=1, position=[0, 0, 0]),
    ...         Atom(atomic_number=8, position=[0, 0, 1]),
    ...         Atom(atomic_number=1, position=[0, 1, 1])]
    ... )
    >>> atomic_number_indices(H2O, 1)
    (0, 2)
    >>> atomic_number_indices(H2O, {8, 1})
    (0, 1, 2)
    >>> atomic_number_indices(H2O, {6, 9})
    ()
    """
    if isinstance(atomic_numbers, int):
        return tuple(i for i, an in enumerate(molecule.atomic_numbers) if an == atomic_numbers)
    return tuple(i for i, an in enumerate(molecule.atomic_numbers) if an in atomic_numbers)
