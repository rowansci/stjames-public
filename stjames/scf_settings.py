from .base import Base, LowercaseStrEnum


class UseSOSCF(LowercaseStrEnum):
    ALWAYS = "always"
    UPON_FAILURE = "upon_failure"
    NEVER = "never"


class SCFSettings(Base):
    """
    Settings for SCF convergence.

    :param max_iters: the maximum number of SCF iterations to permit
    :param soscf: whether or not to use SOSCF (second-order SCF).
    """

    max_iters: int = 250
    soscf: UseSOSCF = UseSOSCF.UPON_FAILURE
