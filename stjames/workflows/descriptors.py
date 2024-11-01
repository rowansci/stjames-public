from ..types import UUID
from .workflow import Workflow

Descriptors = dict[str, dict[str, float] | tuple[float | None, ...] | float]


class DescriptorsWorkflow(Workflow):
    # UUID of optimization
    optimization: UUID | None = None

    descriptors: Descriptors | None = None
