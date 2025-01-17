from .base import Base, LowercaseStrEnum


class ComputeType(LowercaseStrEnum):
    CPU = "cpu"
    GPU = "gpu"


class ComputeSettings(Base):
    compute_type: ComputeType = ComputeType.CPU
