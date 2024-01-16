import pydantic

from .base import Base, LowercaseStrEnum


class ERIStrategy(LowercaseStrEnum):
    # direct SCF, as per Almlof/Faegri/Korsell
    DIRECT = "direct"

    # Raffinetti-style supermatrix
    SUPERMATRIX = "supermatrix"

    # let the software choose between direct and supermatrix based on size
    AUTO = "auto"

    # resolution of the identity for J and K
    RIJK = "rijk"


class IntSettings(Base):
    strategy: ERIStrategy = ERIStrategy.AUTO

    # these will get overwritten by ``mode`` anyway, for the most part
    eri_threshold: pydantic.PositiveFloat = 1e-9
    csam_multiplier: pydantic.PositiveFloat = pydantic.Field(default=1, ge=1)
    pair_overlap_threshold: pydantic.PositiveFloat = 1e-10

    @property
    def resolution_of_the_identity(self) -> bool:
        """Abstracting in case we add RIJCOSX, etc later"""
        return self.strategy == ERIStrategy.RIJK
