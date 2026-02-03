from typing import Annotated, Literal

from pydantic import AfterValidator

from ..base import Base, LowercaseStrEnum, round_float
from ..basis_set import BasisSet
from .workflow import MoleculeWorkflow


class SAPTMethod(LowercaseStrEnum):
    SAPT0 = "sapt0"


class SAPTSettings(Base):
    """
    How SAPT calculations should be performed.

    :param method: which SAPT method
    :param basis_set: which basis set to employ
    """

    method: SAPTMethod = SAPTMethod.SAPT0
    basis_set: BasisSet = BasisSet(name="jun-cc-pVDZ")

    def __str__(self) -> str:
        """
        >>> str(SAPTSettings())
        'SAPT0(jun-cc-pVDZ)'
        """
        return f"{self.method.value.upper()}({self.basis_set.name})"


class SAPT0Result(Base):
    """
    Stores the result of a SAPT0 calculation.

    :param electrostatic_interaction_energy: electrostatic interaction energy, in kcal/mol
    :param exchange_interaction_energy: exchange interaction energy, in kcal/mol
    :param dispersion_interaction_energy: dispersion interaction energy, in kcal/mol
    :param induction_interaction_energy: induction interaction energy, in kcal/mol
    :param result_type: to disambiguate from possible future alternate SAPT returns
    """

    electrostatic_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    exchange_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    dispersion_interaction_energy: Annotated[float, AfterValidator(round_float(3))]
    induction_interaction_energy: Annotated[float, AfterValidator(round_float(3))]

    result_type: Literal["sapt0"] = "sapt0"


class SymmetryAdaptedPerturbationTheoryWorkflow(MoleculeWorkflow):
    """
    Performs a SAPT calculation.

    Inherited:
    :param initial_molecule: Molecule in question

    New:
    :param fragment1_indices: which atoms go to fragment #1 (fragment #2 takes the rest)
    :param sapt_settings: settings for SAPT calculations

    Results:
    :param sapt_result: results from SAPT
    """

    fragment1_indices: list[int]

    sapt_settings: SAPTSettings = SAPTSettings()

    sapt_result: SAPT0Result | None = None
