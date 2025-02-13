"""pKa workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import Base, round_float, round_optional_float
from ..mode import Mode
from .workflow import DBCalculation, MoleculeWorkflow


class pKaMicrostate(Base):
    """
    A microstate for pKa calculations.

    :param atom_index: index of the atom
    :param structures: DBCalculation for the microstate
    :param deltaG: relative free energy
    :param pka: pKa
    """

    atom_index: int
    structures: list[DBCalculation] = []
    deltaG: Annotated[float, AfterValidator(round_float(3))]
    pka: Annotated[float, AfterValidator(round_float(3))]


class pKaWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating pKa.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    New:
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
