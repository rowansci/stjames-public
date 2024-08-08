from typing import Self, Sequence

from pydantic import BaseModel, PositiveInt, field_validator, model_validator

from ..calculation import Calculation
from ..constraint import Constraint
from ..mode import Mode
from ..molecule import Molecule
from ..solvent import Solvent
from .multistage_opt import MultiStageOptInput
from .workflow import WorkflowInput, WorkflowResults

# the id of a mutable object may change, thus using object()
_sentinel_mso = object()


class SpinStatesInput(WorkflowInput):
    """
    Workflow for computing spin states of molecules.

    Uses the modes from MultiStageOptInput.

    Inherited
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    :param states: multiplicities of the spin state targetted
    :param mode: Mode to use
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state

    >>> from stjames.molecule import Atom
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> ss = SpinStatesInput(initial_molecule=He, states=[1, 3, 5], mode=Mode.RAPID, solvent="water")
    >>> str(ss)
    '<SpinStatesInput [1, 3, 5] RAPID>'
    """

    states: list[PositiveInt]
    # Need to use a sentinel object to make both mypy and pydantic happy
    multistage_opt_settings: MultiStageOptInput = _sentinel_mso  # type: ignore [assignment]
    solvent: Solvent | None = None
    xtb_preopt: bool | None = None
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False

    def __repr__(self) -> str:
        if self.mode != Mode.MANUAL:
            return f"<SpinStatesInput {self.states} {self.mode.name}>"

        return f"<SpinStatesInput {self.states} {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        return self.multistage_opt_settings.level_of_theory

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

        match self.mode, self.multistage_opt_settings:
            case (Mode.DEBUG, _):
                raise NotImplementedError("Unsupported mode: DEBUG")

            case (Mode.MANUAL, mso) if mso is _sentinel_mso:
                raise ValueError("Must specify multistage_opt_settings with MANUAL mode")
            case (Mode.MANUAL, _):
                pass

            case (mode, mso) if mso is not _sentinel_mso:
                raise ValueError(f"Cannot specify multistage_opt_settings with {mode=}, {mso=}")

            case (mode, _):
                self.multistage_opt_settings = MultiStageOptInput(
                    initial_molecule=self.initial_molecule,
                    mode=self.mode,
                    solvent=self.solvent,
                    xtb_preopt=self.xtb_preopt,
                    constraints=self.constraints,
                    transition_state=self.transition_state,
                )

        return self


class SpinState(BaseModel):
    """
    The result of a SpinState calculation.

    :param multiplicity: multiplicity of the SpinState
    :param relative_energy: energy of the optimized Molecule
    :param calculation: the Calculation that produced this SpinState

    >>> from stjames.molecule import Atom
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> calc = Calculation(molecules=[He])
    >>> SpinState(multiplicity=3, relative_energy=0, calculation=calc)
    <SpinState 3 0.0>
    """

    multiplicity: PositiveInt
    relative_energy: float
    calculation: Calculation

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<SpinState {self.multiplicity} {self.relative_energy:.1f}>"

    @property
    def molecule(self) -> Molecule:
        """The last molecule is always the one we want."""
        return self.calculation.molecules[-1]


class SpinStatesResults(WorkflowResults):
    """
    Results of a SpinStates workflow.

    :param spin_states: list of SpinStates

    >>> from stjames.molecule import Atom
    >>> mols = [
    ...     Molecule(charge=0, multiplicity=i, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    ...     for i in [1, 3, 5]
    ... ]
    >>> calcs = [Calculation(molecules=[mol]) for mol in mols]
    >>> relative_energies = [0, 1, 3]
    >>> spin_states = [
    ...     SpinState(molecule=mol, multiplicity=mol.multiplicity, relative_energy=re, calculation=calc)
    ...     for mol, calc, re in zip(mols, calcs, relative_energies)
    ... ]
    >>> ssr = SpinStatesResults(spin_states=spin_states)
    >>> ssr
    <SpinStatesResults [1, 3, 5]>
    >>> ssr.multiplicities
    [1, 3, 5]
    >>> ssr.relative_energies
    [0.0, 1.0, 3.0]
    >>> print(ssr)
    SpinStatesResults
      1:  0.00
      3:  1.00
      5:  3.00
    """

    spin_states: list[SpinState]

    def __repr__(self) -> str:
        return f"<SpinStatesResults {self.multiplicities}>"

    def __str__(self) -> str:
        return "SpinStatesResults\n" + "\n".join(f"  {ss.multiplicity}: {ss.relative_energy:>5.2f}" for ss in self.spin_states)

    def __len__(self) -> int:
        return len(self.spin_states)

    @property
    def molecules(self) -> list[Molecule]:
        return [ss.molecule for ss in self.spin_states]

    @property
    def multiplicities(self) -> list[PositiveInt]:
        return [ss.multiplicity for ss in self.spin_states]

    @property
    def relative_energies(self) -> list[float]:
        return [ss.relative_energy for ss in self.spin_states]

    @field_validator("spin_states")
    @classmethod
    def validate_spin_states(cls, spin_states: list[SpinState]) -> list[SpinState]:
        """Ensure that all spin states have consistent results."""
        mults = [ss.multiplicity for ss in spin_states]
        relative_energies = [ss.relative_energy for ss in spin_states]

        if len(mults) != len(set(mults)):
            raise ValueError(f"Duplicate multiplicities found: {mults}")

        if any((mults[0] - mult) % 2 for mult in mults[1:]):
            raise ValueError(f"Inconsistent multiplicities found: {mults}")

        if any(re < 0 for re in relative_energies):
            raise ValueError(f"All relative energies must be >= 0, got: {relative_energies}")

        if any(len(spin_states[0].molecule) != len(ss.molecule) for ss in spin_states[1:]):
            raise ValueError("Inconsistent molecules found")

        return spin_states
