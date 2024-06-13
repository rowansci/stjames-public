from .base import LowercaseStrEnum


class Mode(LowercaseStrEnum):
    # choose based on job type
    AUTO = "auto"

    # very rapid, possibly unreliable
    RECKLESS = "reckless"

    # standard for single-point energies
    RAPID = "rapid"

    # standard for gradient / geometry optimizations
    CAREFUL = "careful"

    # very careful!
    METICULOUS = "meticulous"

    # for debugging purposes
    DEBUG = "debug"

    # set parameters manually, prevent overwriting
    MANUAL = "manual"
