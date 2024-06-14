from enum import Enum
from typing import Annotated, Hashable, TypeVar

import numpy as np
import pydantic

T = TypeVar("T")


class Base(pydantic.BaseModel):
    @pydantic.field_validator("*", mode="before")
    @classmethod
    def coerce_numpy(cls, val: T) -> T | list:
        if isinstance(val, np.ndarray):
            return val.tolist()
        else:
            return val


class LowercaseStrEnum(str, Enum):
    """Enum where hyphens and case are ignored."""

    @classmethod
    def _missing_(cls, value: str) -> str | None:  # type: ignore
        # Type note: technically breaking Liskov, value: object in Enum
        for member in cls:
            if member.lower().replace("-", "") == value.lower().replace("-", ""):
                return member
        return None


# cf. https://github.com/pydantic/pydantic-core/pull/820#issuecomment-1670475909
H = TypeVar("H", bound=Hashable)


def _validate_unique_list(v: list[H]) -> list[H]:
    if len(v) != len(set(v)):
        raise ValueError("this list must be unique, and isn't!")
    return v


UniqueList = Annotated[list[H], pydantic.AfterValidator(_validate_unique_list)]
