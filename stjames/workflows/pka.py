from typing import Annotated

from pydantic import AfterValidator

from ..base import Base, round_float, round_optional_float
from ..mode import Mode
from .workflow import DBCalculation, Workflow


class pKaMicrostate(Base):
    atom_index: int
    structures: list[DBCalculation] = []
    deltaG: Annotated[float, AfterValidator(round_float(3))]
    pka: Annotated[float, AfterValidator(round_float(3))]


class pKaWorkflow(Workflow):
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
