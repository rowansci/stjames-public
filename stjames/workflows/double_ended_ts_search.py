"""Double ended transition-state-search workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, PositiveFloat, model_validator

from ..molecule import Molecule
from ..optimization.freezing_string_method import FSMSettings
from ..settings import Settings
from ..types import UUID, round_list
from .workflow import Workflow


class DoubleEndedTSSearchWorkflow(Workflow):
    """
    Settings for running a double-ended-transition-state-search Workflow.

    # Inputs
    :param reactant: reactant Molecule
    :param product: product Molecule
    :param calculation_settings: Settings for the calculations
    :param search_settings: Settings for the search
    :param optimize_inputs: Whether to optimize the input reactant and product
    :param optimize_ts: Whether to optimize the guess transition state
    :param reactant_string_distances: distances along the reactant string

    # Results
    :param reactant_string_distances: distances along the string
    :param product_string_distances: distances along the string (starting from the end of the reactant)
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

    # Results
    reactant_string_distances: Annotated[list[PositiveFloat], AfterValidator(round_list(5))] = []
    product_string_distances: Annotated[list[PositiveFloat], AfterValidator(round_list(5))] = []

    reactant_calculation_uuids: list[UUID | None] = []
    product_calculation_uuids: list[UUID | None] = []

    ts_guess_calculation_uuid: UUID | None = None

    @model_validator(mode="after")
    def validate_input_molecules(self) -> Self:
        """Validate that reactant and product are compatible."""

        if self.reactant.charge != self.product.charge:
            raise ValueError("Reactant and product must have the same charge.")
        if self.reactant.multiplicity != self.product.multiplicity:
            raise ValueError("Reactant and product must have the same multiplicity.")

        if len(self.reactant) != len(self.product):
            raise ValueError("Reactant and product must have the same number of atoms.")

        for r_atom, p_atom in zip(self.reactant.atoms, self.product.atoms, strict=True):
            if r_atom.atomic_number != p_atom.atomic_number:
                raise ValueError("Reactant and product must have the same atom ordering.")

        return self

    @property
    def path_uuids(self) -> list[UUID | None]:
        """Return the path from reactant to product."""
        return self.reactant_calculation_uuids + self.product_calculation_uuids

    @property
    def distances(self) -> list[float]:
        """Return the path from reactant to product."""
        return self.reactant_string_distances + self.product_string_distances[::-1]
