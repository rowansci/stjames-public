from typing import Self, Sequence

from pydantic import BaseModel, Field, PositiveInt, field_validator, model_validator

from ..constraint import Constraint
from ..mode import Mode
from ..solvent import Solvent
from .multistage_opt import MultiStageOptWorkflow
from .workflow import UUID, Workflow

# the id of a mutable object may change, thus using object()
_sentinel_msow = object()


class SpinState(BaseModel):
    """
    The result of a SpinState calculation.

    :param multiplicity: multiplicity of the SpinState
    :param energy: energy of the optimized Molecule
    :param calculations: the UUIDs of the Calculations that produced this SpinState

    >>> from stjames.molecule import Atom, Molecule
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> SpinState(multiplicity=3, energy=-13.4291, calculation=['8a031a27-30d2-4ac7-8ade-efae9e9fc94a'])
    <SpinState 3 -13.429>
    """

    multiplicity: PositiveInt
    energy: float
    calculation: list[UUID]

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<SpinState {self.multiplicity} {self.energy:.3f}>"


class SpinStatesWorkflow(Workflow):
    """
    Workflow for computing spin states of molecules.

    Uses the modes from MultiStageOptWorkflow.

    Influenced by:
    [Performance of Quantum Chemistry Methods for Benchmark Set of Spin–State
    Energetics Derived from Experimental Data of 17 Transition Metal Complexes
    (SSE17)](https://chemrxiv.org/engage/chemrxiv/article-details/66a8b15cc9c6a5c07a792487)

    Inherited
    :param initial_molecule: Molecule of interest

    :param mode: Mode for workflow
    :param states: multiplicities of the spin state targetted
    :param mode: Mode to use
    :param multistage_opt_workflow: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state

    :param spin_states: resulting spin states data

    >>> from stjames.molecule import Atom, Molecule
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> ss = SpinStatesWorkflow(initial_molecule=He, states=[1, 3, 5], mode=Mode.RAPID, solvent="water")
    >>> str(ss)
    '<SpinStatesWorkflow [1, 3, 5] RAPID>'
    """

    mode: Mode
    states: list[PositiveInt]
    # Need to use a sentinel object to make both mypy and pydantic happy
    multistage_opt_workflow: MultiStageOptWorkflow = _sentinel_msow  # type: ignore [assignment]
    solvent: Solvent | None = None
    xtb_preopt: bool | None = None
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False

    spin_states: list[SpinState] = Field(default_factory=list)

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        if self.mode != Mode.MANUAL:
            return f"<SpinStatesWorkflow {self.states} {self.mode.name}>"

        return f"<SpinStatesWorkflow {self.states} {self.level_of_theory}>"

    def __len__(self) -> int:
        return len(self.states)

    @property
    def level_of_theory(self) -> str:
        return self.multistage_opt_workflow.level_of_theory

    @field_validator("states")
    @classmethod
    def validate_states(cls, states: list[PositiveInt]) -> list[PositiveInt]:
        """Confirm that all spin states are valid."""
        if not states:
            raise ValueError("Expected at least one spin state.")

        if any((states[0] - mult) % 2 for mult in states):
            raise ValueError(f"Inconsistent multiplicities found: {states}")

        return states

    @model_validator(mode="after")
    def set_mode_and_settings(self) -> Self:
        if self.mode == Mode.AUTO:
            self.mode = Mode.RAPID

        match self.mode, self.multistage_opt_workflow:
            case (Mode.DEBUG, _):
                raise NotImplementedError("Unsupported mode: DEBUG")

            case (Mode.MANUAL, msow) if msow is _sentinel_msow:
                raise ValueError("Must specify multistage_opt_workflow with MANUAL mode")
            case (Mode.MANUAL, _):
                pass

            case (mode, msow) if msow is not _sentinel_msow:
                raise ValueError(f"Cannot specify multistage_opt_workflow with {mode=}, {msow=}")

            case (mode, _):
                self.multistage_opt_workflow = MultiStageOptWorkflow(
                    initial_molecule=self.initial_molecule,
                    mode=self.mode,
                    solvent=self.solvent,
                    xtb_preopt=self.xtb_preopt,
                    constraints=self.constraints,
                    transition_state=self.transition_state,
                )

        return self

    def str_results(self) -> str:
        return "SpinStatesResults\n" + "\n".join(f"  {ss.multiplicity}: {ss.energy:>5.2f}" for ss in self.spin_states)

    @property
    def energies(self) -> list[float]:
        return [ss.energy for ss in self.spin_states]

    @field_validator("spin_states")
    @classmethod
    def validate_spin_states(cls, spin_states: list[SpinState]) -> list[SpinState]:
        """Ensure that all spin states have consistent results."""
        mults = [ss.multiplicity for ss in spin_states]

        if len(mults) != len(set(mults)):
            raise ValueError(f"Duplicate multiplicities found: {mults}")

        if any((mults[0] - mult) % 2 for mult in mults[1:]):
            raise ValueError(f"Inconsistent multiplicities found: {mults}")

        return spin_states
