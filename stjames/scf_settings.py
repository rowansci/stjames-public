from pydantic import BaseModel


class SCFSettings(BaseModel):
    """
    Settings for SCF convergence.

    :param max_iters: the maximum number of SCF iterations to permit
    :param soscf: whether or not to use SOSCF (second-order SCF).
        SOSCF can converge tricky cases, but is also slower than regular SCF algorithms.
        If `None`, SOSCF will automatically be attempted where regular (DIIS) convergence fails.
    """

    max_iters: int = 250
    soscf: bool | None = None
