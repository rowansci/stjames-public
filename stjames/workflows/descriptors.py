from .workflow import Workflow

Descriptors = dict[str, dict[str, float] | tuple[float | None, ...] | float]


class DescriptorsWorkflow(Workflow):
    # uuid of optimization
    optimization: str | None = None

    descriptors: Descriptors | None = None
