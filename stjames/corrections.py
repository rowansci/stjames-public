from .base import LowercaseStrEnum


class Correction(LowercaseStrEnum):
    """Various post hoc corrections."""

    # Grimme's D3 dispersion correction, with Becke–Johnson damping
    D3BJ = "d3bj"

    # Grimme's geometric counterpoise correction
    GCP = "gcp"
