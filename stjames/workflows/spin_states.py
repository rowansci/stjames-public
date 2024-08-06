from typing import Any, Self, Sequence

from pydantic import PositiveInt, model_validator

from ..constraint import Constraint
from ..mode import Mode
from ..solvent import Solvent
from .multistage_opt import MultiStageOpt
from .workflow import Workflow


class SpinStates(Workflow):
    """
    Workflow for computing spin states of molecules.

    Uses the modes from MultiStageOpt.

    Inherited
    :param initial_molecule: Molecule of interest

    :param states: multiplicities of the spin state targetted
    :param mode: Mode to use
    :param MultiStageOpt: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state

    >>> from stjames.molecule import Molecule, Atom
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> ss = SpinStates(initial_molecule=He, states=[1, 3, 5], mode=Mode.RAPID, solvent="water")
    >>> str(ss)
    '<SpinStates [1, 3, 5] r2scan_3c/cpcm(water)//gfn2_xtb>'
    """

    states: list[PositiveInt]
    mode: Mode | None = None
    multistage_opt: MultiStageOpt | None = None
    solvent: Solvent | None = None
    xtb_preopt: bool | None = None
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False

    def __str__(self) -> str:
        return f"<SpinStates {self.states} {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        assert self.multistage_opt
        return self.multistage_opt.level_of_theory

    @model_validator(mode="after")
    def validate_states(self) -> Self:
        """Confirm that all spin states are valid."""
        multiplicity = self.initial_molecule.multiplicity
        if not self.states:
            raise ValueError("Expected at least one spin state.")

        for state in self.states:
            # If of a different parity than molecule
            if (state - multiplicity) % 2:
                raise ValueError(f"Invalid multiplicity for molecule. Molecule started with {multiplicity=}, got {state=} of different parity")

        return self

    @model_validator(mode="before")
    def check_mode(cls, values: dict[str, Any]) -> dict[str, Any]:
        mso = values.get("multistage_opt")
        match mode := values.get("mode"):
            case None:
                values["mode"] = Mode.MANUAL
                if not mso:
                    raise ValueError("Must specify mode or multistage_opt")
            case Mode.MANUAL:
                if not mso:
                    raise ValueError("Must specify multistage_opt with MANUAL mode")
            case Mode.DEBUG:
                raise NotImplementedError(f"Unsupported mode: {mode}")
            case _:
                if mso:
                    raise ValueError(f"Cannot specify multistage_opt with {mode=}")

                if mode == mode.AUTO:
                    mode = mode.RAPID

                values["multistage_opt"] = MultiStageOpt(
                    initial_molecule=values["initial_molecule"],
                    mode=mode,
                    solvent=values.get("solvent"),
                    xtb_preopt=values.get("xtb_preopt"),
                    constraints=values.get("constraints", ()),
                    transition_state=values.get("transition_state", False),
                )

        return values
