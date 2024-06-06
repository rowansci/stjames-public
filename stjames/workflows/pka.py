from .workflow import Workflow


class pKaWorkflow(Workflow):
    pka_range: tuple[float, float] = (2, 12)
    deprotonate_elements: list[int] = [7, 8, 16]
    deprotonate_atoms: list[int] = []
    protonate_elements: list[int] = [7]
    protonate_atoms: list[int] = []

    reasonableness_buffer: float = 5
