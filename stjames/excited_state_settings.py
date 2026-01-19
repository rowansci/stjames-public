from abc import ABC
from typing import Annotated, Literal, Self

from pydantic import BaseModel, Field, PositiveInt, model_validator


class ExcitedStateSettings(BaseModel, ABC):
    """Settings for excited state calculations."""


class TDDFTSettings(ExcitedStateSettings):
    """
    Settings for TDDFT calculations.

    New:
    :param tda: use Tamm-Dancoff approximation
    :param num_excitations: number of excitations to calculate
    :param target_root: root to target (for gradient/optimization)
    """

    tda: bool = True
    num_excitations: PositiveInt = 5
    target_root: PositiveInt | None = None

    settings_type: Literal["tddft_settings"] = "tddft_settings"

    @model_validator(mode="after")
    def validate_setup(self) -> Self:
        if self.target_root and (self.target_root > self.num_excitations):
            raise ValueError("target_root cannot be greater than num_excitations.")

        return self


ExcitedStateSettingsUnion = Annotated[TDDFTSettings, Field(discriminator="settings_type")]
