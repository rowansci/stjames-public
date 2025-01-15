from pytest import mark

from stjames.pdb import fetch_pdb, pdb_from_string


def test_1ema() -> None:
    """Green fluorescent protein."""
    fetch_pdb("1EMA")


def test_read_pdb() -> None:
    """Rest reading of a pdb string."""
    with open("tests/data/1ema.pdb") as f:
        data = f.read()
    pdb = pdb_from_string(data)
    assert pdb.description.code == "1EMA"


# fmt: off
@mark.regression
@mark.parametrize(
    "code",
    [
        # Codes from molecules of the month August 2024–January 2025
        "7S6B", "8UCS", "8UPL", "7CGO", "2ZVY", "1F4V", "6YKM", "6E10", "6E11",
        "3VCM", "2X0B", "6OS0", "1N9U", "1O86", "2V0Z", "6KI1", "6KI2", "7YYO",
        "6TJV", "7EGL", "7CYF", "7EGK", "7ZCG", "3FRT", "6AP1",
        # Codes created by o1
        "1CRN", "1MBN", "4HHB", "1HHO", "1BNA", "1CAG", "2JHO", "1EVV", "3ZOJ",
        "4AGG", "2Y69", "6R1V", "6ND2", "7NZ6", "1S72", "3G5U", "7DFT", "6AI0",
        "6NG2", "1A0I", "1B7C", "1C8R", "1D4T", "1E7O", "1F9J", "1G5K", "1H8L",
        "1I2M", "1J3N", "1K4P", "1M7R", "1N8S", "1O9T", "1P0U", "1Q1V", "1R2W",
        "1S3X", "1T4Y", "1U5Z", "1V6A", "1W7B", "1X8C", "1Y9D", "1Z0E", "2A1F",
        "2B2G", "2C3H", "2D4I", "2E5J", "2F6K", "2G7L", "2H8M", "2I9N", "2J0O",
        "2K1P", "2L2Q", "2N4S", "2O5T", "2P6U", "2Q7V", "2R8W", "2V2A", "2W3B",
        "2X4C", "2Y5D", "2Z6E", "3A7F", "3B8G", "3C9H", "3D0I", "3E1J", "3F2K",
        "3G3L", "3H4M", "3I5N", "3J6O", "3K7P", "3L8Q", "3N0S", "3O1T", "3P2U",
        "3Q3V", "3S5X", "3T6Y", "3U7Z", "3V8A", "3W9B", "3X0C", "4A3F", "4B4G",
        "4C5H", "4D6I", "4E7J", "4F8K", "4G9L", "4H0M", "4I1N", "4J2O", "4K3P",
        "4L4Q", "4M5R", "4N6S", "4O7T", "4P8U", "4Q9V", "4R0W", "4S1X",
    ]
)  # fmt: on
def test_pdb(code: str) -> None:
    fetch_pdb(code)
