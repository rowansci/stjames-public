from .base import LowercaseStrEnum


class Correction(LowercaseStrEnum):
    """Various post hoc corrections."""

    # Grimme's D3 dispersion correction, with Becke–Johnson damping
    D3BJ = "d3bj"

    # Grimme's D3 dispersion correction, *without* Becke–Johnson damping
    D3 = "d3"

    # Grimme's D4 dispersion correction
    D4 = "d4"

    # Grimme's geometric counterpoise correction
    GCP = "gcp"
