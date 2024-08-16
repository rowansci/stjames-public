from enum import Enum
from typing import Annotated, Any, Hashable, TypeVar

import numpy as np
import pydantic

_T = TypeVar("_T")


class Base(pydantic.BaseModel):
    @pydantic.field_validator("*", mode="before")
    @classmethod
    def coerce_numpy(cls, val: _T) -> _T | list[Any]:
        if isinstance(val, np.ndarray):
            return val.tolist()  # type: ignore [no-any-return, unused-ignore]

        return val


class LowercaseStrEnum(str, Enum):
    """Enum where hyphens, underscores, and case are ignored."""

    @classmethod
    def _missing_(cls, value: object) -> str | None:
        for member in cls:
            if isinstance(value, str):
                if member.lower().replace("-", "").replace("_", "") == value.lower().replace("-", "").replace("_", ""):
                    return member
        return None


# cf. https://github.com/pydantic/pydantic-core/pull/820#issuecomment-1670475909
_H = TypeVar("_H", bound=Hashable)


def _validate_unique_list(v: list[_H]) -> list[_H]:
    if len(v) != len(set(v)):
        raise ValueError("this list must be unique, and isn't!")
    return v


UniqueList = Annotated[list[_H], pydantic.AfterValidator(_validate_unique_list)]
