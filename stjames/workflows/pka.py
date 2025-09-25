"""pKa workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, model_validator

from ..base import Base, LowercaseStrEnum, round_float, round_optional_float
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


class MicroscopicpKaMethod(LowercaseStrEnum):
    AIMNET2_WAGEN2024 = "aimnet2_wagen2024"
    CHEMPROP_NEVOLIANUS2025 = "chemprop_nevolianus2025"


class pKaMicrostate(Base):
    """
    A microstate for pKa calculations.

    :param atom_index: index of the atom
    :param smiles: SMILES of the microstate
    :param structures: DBCalculation for the microstate
    :param deltaG: relative free energy (where applicable)
    :param pka: pKa value associated with this microstate
    :uncertainty: uncertainty in the pKa prediction
    """

    atom_index: int
    smiles: str | None = None
    structures: list[DBCalculation] = []
    deltaG: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    pka: Annotated[float, AfterValidator(round_float(3))]
    uncertainty: Annotated[float | None, AfterValidator(round_optional_float(3))] = None


class pKaWorkflow(SMILESWorkflow, MoleculeWorkflow):
    """
    Workflow for calculating pKa.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    :param initial_smiles: SMILES of the molecule of interest

    New:
    :param microscopic_pka_method: method used for pka prediciton
    :param solvent: solvent for the pka prediction
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

    microscopic_pka_method: MicroscopicpKaMethod = MicroscopicpKaMethod.CHEMPROP_NEVOLIANUS2025
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
        match self.microscopic_pka_method:
            case MicroscopicpKaMethod.AIMNET2_WAGEN2024:
                if self.solvent is not Solvent.WATER:
                    raise ValueError(f"{self.microscopic_pka_method} only supports water")
            case MicroscopicpKaMethod.CHEMPROP_NEVOLIANUS2025:
                if self.solvent not in CHEMPROP_NEVOLIANUS2025_ALLOWED_SOLVENTS:
                    raise ValueError(f"Solvent `{self.solvent}` is not supported by method `{self.microscopic_pka_method}`.")
                if len(self.protonate_atoms) or len(self.deprotonate_atoms):
                    raise ValueError(f"Method `{self.microscopic_pka_method}` does not support selecting atoms by number.")
        return self

    @model_validator(mode="after")
    def validate_mol_input(self) -> Self:
        """Ensure that only one of initial_molecule or initial_smiles is set."""

        if not (bool(self.initial_smiles) ^ bool(self.initial_molecule)):
            raise ValueError("Can only set one of initial_molecule should and initial_smiles")

        return self
