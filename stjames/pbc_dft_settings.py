from pydantic import BaseModel, PositiveFloat, PositiveInt

from .base import LowercaseStrEnum


class PBCDFTSmearing(LowercaseStrEnum):
    """Smearing types for occupations in PBC DFT calculations."""

    MV = "mv"  # Marzari–Vanderbilt cold smearing
    MP = "mp"  # Methfessel–Paxton
    FD = "fd"  # Fermi–Dirac
    GAUSSIAN = "gaussian"  # Gaussian smearing


class PBCDFTSettings(BaseModel):
    """
    PBC DFT settings.

    :param pw_cutoff: plane-wave kinetic-energy cutoff (Hartree)
    :param kpoints: Monkhorst–Pack k-point-grid dimensions
    :param smearing: occupations smearing type
    :param degauss: smearing width, if relevant (Hartree)
    """

    pw_cutoff: PositiveFloat
    kpoints: tuple[PositiveInt, PositiveInt, PositiveInt]

    smearing: PBCDFTSmearing | None = None
    degauss: PositiveFloat = 0.005
