from typing import Annotated, Literal

from pydantic import AfterValidator

from ..base import Base, LowercaseStrEnum, round_optional_float
from ..basis_set import BasisSet
from .workflow import MoleculeWorkflow


class SAPTMethod(LowercaseStrEnum):
    SAPT0 = "sapt_0"


class SAPTSettings(Base):
    """
    How SAPT calculations should be performed.

    :param method: which SAPT method
    :param basis_set: which basis set to employ
    """

    method: SAPTMethod = SAPTMethod.SAPT0
    basis_set: BasisSet = BasisSet(name="jun-cc-pVDZ")


class SAPT0Result(Base):
    """
    Stores the result of a SAPT0 calculation.

    :param electrostatic_interaction_energy: the electrostatic interaction energy (in kcal/mol)
    :param exchange_interaction_energy: the exchange interaction energy (in kcal/mol)
    :param dispersion_interaction_energy: the dispersion interaction energy (in kcal/mol)
    :param induction_interaction_energy: the induction interaction energy (in kcal/mol)
    :param result_type: to disambiguate from possible future alternate SAPT returns
    """

    electrostatic_interaction_energy: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    exchange_interaction_energy: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    dispersion_interaction_energy: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    induction_interaction_energy: Annotated[float | None, AfterValidator(round_optional_float(3))] = None

    result_type: Literal["sapt0"] = "sapt0"


class SymmetryAdaptedPerturbationTheoryWorkflow(MoleculeWorkflow):
    """
    Performs a SAPT calculation.

    Inherited:
    :param initial_molecule: the molecule in question
    :param fragment1_indices: which atoms go to fragment #1
    :param fragment2_indices: which atoms go to fragment #2
    :param sapt_settings: the settings for SAPT calculations

    Results:
    :param sapt_result: the results from SAPT
    """

    fragment1_indices: list[int]
    fragment2_indices: list[int]

    sapt_settings: SAPTSettings = SAPTSettings()

    sapt_result: SAPT0Result | None = None
