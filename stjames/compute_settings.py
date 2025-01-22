from .base import Base, LowercaseStrEnum


class ComputeType(LowercaseStrEnum):
    CPU = "cpu"
    GPU = "gpu"
    AUTO = "auto"


class ComputeSettings(Base):
    requested_compute_type: ComputeType = ComputeType.AUTO
    compute_type_used: ComputeType | None = None
