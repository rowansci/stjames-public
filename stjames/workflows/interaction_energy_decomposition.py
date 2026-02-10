from typing import Annotated, Literal

from pydantic import AfterValidator, PositiveInt

from ..base import Base, LowercaseStrEnum, round_float
from ..basis_set import BasisSet
from .workflow import MoleculeWorkflow


class EnergyDecompositionMethod(LowercaseStrEnum):
    SAPT0 = "sapt0"


class EnergyDecompositionSettings(Base):
    """
    How energy-decomposition calculations should be performed.

    :param method: which energy decomposition method
    :param basis_set: which basis set to employ
    """

    method: EnergyDecompositionMethod = EnergyDecompositionMethod.SAPT0
    basis_set: BasisSet = BasisSet(name="jun-cc-pVDZ")

    def __str__(self) -> str:
        """
        >>> str(EnergyDecompositionSettings())
        'SAPT0(jun-cc-pVDZ)'
        """
        return f"{self.method.value.upper()}({self.basis_set.name})"


class SAPT0Result(Base):
    """
    Stores the result of a SAPT0 calculation.

    :param total_interaction_energy: total interaction energy, in kcal/mol
    :param electrostatic_interaction_energy: electrostatic interaction energy, in kcal/mol
    :param exchange_interaction_energy: exchange interaction energy, in kcal/mol
    :param dispersion_interaction_energy: dispersion interaction energy, in kcal/mol
    :param induction_interaction_energy: induction interaction energy, in kcal/mol
    :param result_type: to disambiguate from possible future alternate returns
    """

    total_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    electrostatic_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    exchange_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    dispersion_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    induction_interaction_energy: Annotated[float, AfterValidator(round_float(3))]

    result_type: Literal["sapt0"] = "sapt0"


class InteractionEnergyDecompositionWorkflow(MoleculeWorkflow):
    """
    Performs an interaction-energy-decomposition calculation.

    Inherited:
    :param initial_molecule: Molecule in question

    New:
    :param fragment1_indices: which atoms go to fragment #1 (fragment #2 takes the rest)
    :param energy_decomposition_settings: settings for SAPT calculations

    Results:
    :param energy_decomposition_result: results from the energy decomposition
    """

    fragment1_indices: list[PositiveInt]

    energy_decomposition_settings: EnergyDecompositionSettings = EnergyDecompositionSettings()

    energy_decomposition_result: SAPT0Result | None = None
