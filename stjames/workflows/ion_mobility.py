"""Ion mobility workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, computed_field, model_validator

from ..base import Base, round_optional_float
from ..data import ELEMENT_SYMBOL
from ..types import UUID, round_list
from .workflow import MoleculeWorkflow


class IonMobilityForcefieldElement(Base):
    """
    A single atom specification for the ion-mobility forcefield.

    :param name: the name of the element (e.g. "Hydrogen")
    :param atomic_number: the element's atomic number
    :param mass: the mass of the element in Daltons (e.g. 1.00794)
    :param sigma: the sigma Lennard-Jones parameter, in Å
    :param epsilon: the epsilon Lennard-Jones parameter, in kcal/mol
    """

    name: str
    atomic_number: int
    mass: float
    sigma: float
    epsilon: float

    @computed_field  # type: ignore[misc, prop-decorator, unused-ignore]
    @property
    def symbol(self) -> str:
        return ELEMENT_SYMBOL[self.atomic_number]


class IonMobilityWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating hydrogen bond basicity.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param protonate: automatically protonate the molecule
    :param temperature: the temperature, in Kelvin
    :param do_csearch: whether to perform a conformational search
    :param do_optimization: whether to perform an optimization
    :param forcefield: the forcefield used to describe atom–gas interactions.
        if None, the default forcefield will be used.

    Results:
    :param conformers: the UUIDs of the conformers
    :param conformer_ccs: the collision cross section (Å**2) per conformer
    :param conformer_ccs_stdev: the uncertainty in the same
    :param conformer_weights: the Boltzmann weights at RT
    :param average_ccs: the Boltzmann-weighted CCS for the ensemble
    :param average_ccs_stdev: the uncertainty in the same
    """

    protonate: bool = False
    temperature: float = 300
    do_csearch: bool = True
    do_optimization: bool = True

    forcefield: list[IonMobilityForcefieldElement] | None = None

    conformers: list[UUID] = []

    conformer_ccs: Annotated[list[float], AfterValidator(round_list(3))] = []
    conformer_ccs_stdev: Annotated[list[float], AfterValidator(round_list(3))] = []
    boltzmann_weights: Annotated[list[float], AfterValidator(round_list(3))] = []

    average_ccs: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    average_ccs_stdev: Annotated[float | None, AfterValidator(round_optional_float(3))] = None

    @model_validator(mode="after")
    def check_supported_atoms(self) -> Self:
        """Validate that user-supplied forcefields have correct atoms."""
        if self.forcefield is not None:
            supported_atoms = set(e.atomic_number for e in self.forcefield)
            if not all(z in supported_atoms for z in self.initial_molecule.atomic_numbers):
                raise ValueError("Provided forcefield does not support all elements in input structure!")

        return self
