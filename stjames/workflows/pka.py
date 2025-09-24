"""pKa workflow."""

from typing import Annotated, Optional, Self

from pydantic import AfterValidator, model_validator

from ..base import Base, round_float, round_optional_float, LowercaseStrEnum
from ..mode import Mode
from ..solvent import Solvent
from .workflow import DBCalculation, MoleculeWorkflow, SMILESWorkflow

CHEMPROP_NEVOLIANUS2025_ALLOWED_SOLVENTS = {
    Solvent.WATER,
    Solvent.DIMETHYLSULFOXIDE,
    Solvent.DIMETHYLFORMAMIDE,
    Solvent.ACETONITRILE,
    Solvent.METHANOL,
    Solvent.ETHANOL,
    Solvent.ETHYLENE_GLYCOL,
    Solvent.N_METHYLPYRROLIDONE,
}


class microscopicpKaMethod(LowecaseStrEnum):
    AIMNET2_WAGEN2024 = "aimnet2_wagen2024"
    CHEMPROP_NEVOLIANUS2025 = "chemprop_nevolianus2025"


class pKaMicrostate(Base):
    """
    A microstate for pKa calculations.

    :param atom_index: index of the atom
    :param smiles: SMILES of the microstate
    :param structures: DBCalculation for the microstate
    :param deltaG: relative free energy
    :param pka: pKa
    :uncertainty: uncertainty
    """

    atom_index: int
    smiles = Optional[str]
    structures: list[DBCalculation] = []
    deltaG: Annotated[Optional[float, AfterValidator(round_float(3))]]
    pka: Annotated[float, AfterValidator(round_float(3))]
    uncertainty: Annotated[Optional[float, AfterValidator(round_float(3))]]


class pKaWorkflow(SMILESWorkflow, MoleculeWorkflow):
    """
    Workflow for calculating pKa.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    :param initial_smiles: SMILES of the molecule of interest

    New:
    :param microscopic_pKa_method: the method used for pka prediciton
    :param solvent: the solvent for the pka prediction
    :param pka_range: range of pKa values to consider
    :param deprotonate_elements: elements to deprotonate
    :param deprotonate_atoms: atoms to deprotonate
    :param protonate_elements: elements to protonate
    :param protonate_atoms: atom indices to protonate
    :param reasonableness_buffer: buffer for pKa reasonableness

    Results:
    :param structures: DBCalculations for the microstates
    :param conjugate_acids: conjugate acid microstates
    :param conjugate_bases: conjugate base microstates
    :param strongest_acid: pKa of the strongest acid
    :param strongest_base: pKa of the strongest base
    """

    mode: Mode = Mode.CAREFUL

    microscopic_pKa_method: microscopicpKaMethod = microscopicpKaMethod.CHEMPROP_NEVOLIANUS2025
    solvent: Solvent = Solvent.WATER
    pka_range: tuple[float, float] = (2, 12)
    deprotonate_elements: list[int] = [7, 8, 16]
    deprotonate_atoms: list[int] = []
    protonate_elements: list[int] = [7]
    protonate_atoms: list[int] = []

    reasonableness_buffer: float = 5

    structures: list[DBCalculation] = []
    conjugate_acids: list[pKaMicrostate] = []
    conjugate_bases: list[pKaMicrostate] = []
    strongest_acid: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    strongest_base: Annotated[float | None, AfterValidator(round_optional_float(3))] = None

    @model_validator(mode="after")
    def check_method_settings(self) -> Self:
        """Check that models with limited domain of applicability are predicting within correct domain."""
        match self.microscopic_pKa_method:
            case microscopicpKaMethod.AIMNET2_WAGEN2024:
                if self.solvent is not Solvent.WATER:
                    raise ValueError(f"Method `{self.microscopic_pKa_method}` can only predict microscopic pKa for water so `solvent` must be `water` only.")
                if len(self.protonate_atoms) > 0 or len(self.deprotonate_atoms) > 0:
                    raise ValueError(f"Method `{self.microscopic_pKa_method}` does not support selecting atoms by number.")
            case microscopicpKaMethod.CHEMPROP_NEVOLIANUS2025:
                if self.solvent not in CHEMPROP_NEVOLIANUS2025_ALLOWED_SOLVENTS:
                    raise ValueError(f"Solvent `{self.solvent}` is invalid for method `{self.microscopic_pKa_method}`.")
