import pydantic
from enum import Enum
import numpy as np
from typing import Any


class Base(pydantic.BaseModel):
    @pydantic.field_validator("*", mode="before")
    @classmethod
    def coerce_numpy(cls, val: Any):
        if isinstance(val, np.ndarray):
            return val.tolist()
        else:
            return val


class LowercaseStrEnum(str, Enum):
    """Enum where hyphens and case are ignored."""

    @classmethod
    def _missing_(cls, value: str):
        for member in cls:
            if member.lower().replace("-", "") == value.lower().replace("-", ""):
                return member
        return None
