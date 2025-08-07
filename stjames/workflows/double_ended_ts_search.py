"""Double ended transition-state-finding workflow."""

from typing import Self

from pydantic import model_validator

from stjames.settings import Settings

from ..molecule import Molecule
from ..optimization.freezing_string_method import FSMSettings
from ..types import UUID
from .workflow import Workflow


class DoubleEndedTSSearchWorkflow(Workflow):
    """
    Settings for running a Double-Ended-Transition-State-Search Workflow.

    # Inputs
    :param reactant: reactant Molecule
    :param product: product Molecule
    :param calculation_settings: Settings for the calculations
    :param search_settings: Settings for the search method
    :param optimize_inputs: Whether to optimize the inupt reactant and product
    :param optimize_ts: Whether to optimize the guess transition state
    :param reactant_string_distances: distances along the reactant string

    # Results
    :param product_string_distances: distances along the string (starting from the reactant as 0)
    :param reactant_calculation_uuids: UUIDs of the calculations for the reactant string nodes
    :param product_calculation_uuids: UUIDs of the calculations for the product string nodes
    :param ts_guess_calculation_uuid: UUID of the calculation for the guess (or optimized) transition state
    """

    reactant: Molecule
    product: Molecule

    calculation_settings: Settings
    search_settings: FSMSettings
    optimize_inputs: bool = False
    optimize_ts: bool = True

    # Data
    reactant_string_distances: list[float] = []
    product_string_distances: list[float] = []

    reactant_calculation_uuids: list[UUID | None] = []
    product_calculation_uuids: list[UUID | None] = []

    ts_guess_calculation_uuid: UUID | None = None

    @model_validator(mode="after")
    def validate_input_molecules(self) -> Self:
        """Validate that reactants and products are compatible."""

        if self.reactant.charge != self.product.charge:
            raise ValueError("Reactants and products must have the same charge.")
        if self.reactant.multiplicity != self.product.multiplicity:
            raise ValueError("Reactants and products must have the same multiplicity.")

        if len(self.reactant) != len(self.product):
            raise ValueError("Reactants and products must have the same number of atoms.")

        for r_atom, p_atom in zip(self.reactant.atoms, self.product.atoms, strict=True):
            if r_atom.atomic_number != p_atom.atomic_number:
                raise ValueError("Reactants and products must have the same atom ordering.")

        return self

    @property
    def path_uuids(self) -> list[UUID | None]:
        """Return the path from reactant to product."""
        return self.reactant_calculation_uuids + self.product_calculation_uuids

    @property
    def distances(self) -> list[float]:
        """Return the path from reactant to product."""
        return self.reactant_string_distances
