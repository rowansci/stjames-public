"""Conformer search workflow."""

from collections import Counter
from typing import Annotated, Self

from pydantic import AfterValidator, Field, model_validator

from ..conformers import ConformerProperties, ConformerSearchMixin, iMTDSettings
from ..molecule import Molecule
from ..types import UUID, FloatPerAtom, round_float_per_atom
from .workflow import MoleculeWorkflow, SMILESWorkflow


class ConformerSearchWorkflow(ConformerSearchMixin, SMILESWorkflow, MoleculeWorkflow):
    """
    ConformerSearch Workflow.

    This workflow supports both SMILES and 3D molecular input. Some conformer generation settings
    support both methods; others (like CREST) require 3D information. Only one should be supplied.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param initial_smiles: SMILES of the molecule of interest
    :param conf_gen_settings: settings for conformer generation
    :param multistage_opt_settings: multi-stage optimization settings
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    New:
    :param initial_conformers: input conformers (if no conformer-generation is requested)
    :param conformer_uuids: list of UUIDs of the Molecules generated
    :param energies: energies of the molecules
    :param conformer_properties: each conformer's properties
    :param ensemble_properties: overall ensemble's properties
    """

    initial_smiles: str = ""
    initial_molecule: Molecule | None = None  # type: ignore [assignment]
    initial_conformers: list[Molecule] = []

    # Results
    conformer_uuids: list[list[UUID | None]] = Field(default_factory=list)
    energies: Annotated[FloatPerAtom, AfterValidator(round_float_per_atom(6))] = Field(default_factory=list)

    conformer_properties: list[ConformerProperties] = []
    ensemble_properties: ConformerProperties | None = None

    @model_validator(mode="after")
    def validate_mol_input(self) -> Self:
        """Ensure that a valid combination of input types is set."""
        if (not self.initial_conformers) and (not (bool(self.initial_smiles) ^ bool(self.initial_molecule))):
            raise ValueError("Can only set one of initial_molecule and initial_smiles")

        if isinstance(self.conf_gen_settings, iMTDSettings) and (self.initial_molecule is None):
            raise ValueError("iMTDSettings requires initial_molecule to be set")

        if self.conf_gen_settings is None and not self.initial_conformers:
            raise ValueError("Need initial_conformers to be set without a conformer-generation method")

        if len(self.initial_conformers) > 1:
            initial_count = Counter(self.initial_conformers[0].atomic_numbers)
            for conformer in self.initial_conformers[1:]:
                if Counter(conformer.atomic_numbers) != initial_count:
                    raise ValueError("Not all molecules in initial_conformers have the same atomic formula")

        return self
