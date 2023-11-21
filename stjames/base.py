import pydantic
from enum import Enum
import numpy as np
from typing import Any, Annotated, TypeVar, Hashable


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


# cf. https://github.com/pydantic/pydantic-core/pull/820#issuecomment-1670475909
# there are python details here I don't really grasp
T = TypeVar("T", bound=Hashable)


def _validate_unique_list(v: list[T]) -> list[T]:
    if len(v) != len(set(v)):
        raise ValueError("this list must be unique, and isn't!")
    return v


UniqueList = Annotated[list[T], pydantic.AfterValidator(_validate_unique_list)]
