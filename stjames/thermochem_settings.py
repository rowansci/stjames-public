import pydantic

from .base import Base


class ThermochemistrySettings(Base):
    """
    Thermochemistry calculation settings.

    :param cutoff_frequency: Cramer/Truhlar quasi-harmonic cutoff, in cm^-1
    :param temperature: temperature for thermochemistry, in K
    :param scaling_factor: frequency scaling factor
    :param concentration: concentration, in M (defaults to 1 atm)
    """

    # Cramer/Truhlar cutoff freq (cm^-1)
    cutoff_frequency: pydantic.NonNegativeFloat = 100

    # temp (K)
    temperature: pydantic.NonNegativeFloat = 298

    # scaling factor, for frequencies
    scaling_factor: pydantic.NonNegativeFloat = 1.0

    # conc (M), defaults to 1 atm
    # (cf. https://github.com/patonlab/GoodVibes/blob/master/goodvibes/GoodVibes.py#L403)
    concentration: pydantic.NonNegativeFloat = 0.0408740470708
