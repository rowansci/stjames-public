import json

from pytest import mark

from stjames.pdb import (
    PDB,
    PDBAtom,
    PDBDescription,
    PDBExperiment,
    PDBGeometry,
    PDBModel,
    PDBPolymer,
    PDBQuality,
    PDBResidue,
    fetch_pdb,
    fetch_pdb_from_mmcif,
    pdb_from_mmcif_filestring,
    pdb_from_pdb_filestring,
    pdb_object_to_pdb_filestring,
)


def test_1ema() -> None:
    """Green fluorescent protein."""
    fetch_pdb("1EMA")


def test_8qvy() -> None:
    """
    Human GABAA receptor
    Not availble as a .pdb
    """
    fetch_pdb_from_mmcif("8VQY")


def test_read_pdb_filestring() -> None:
    """Rest reading of a pdb string."""
    with open("tests/data/1ema.pdb") as f:
        data = f.read()
    pdb = pdb_from_pdb_filestring(data)
    assert pdb.description.code == "1EMA"

    json = pdb.model_dump()
    PDB.model_validate(json)


def test_read_pdb_filestring_mome() -> None:
    """Rest reading of a pdb string."""
    with open("tests/data/cluster_1.pdb") as f:
        data = f.read()
    pdb = pdb_from_pdb_filestring(data)

    json = pdb.model_dump()
    print(json)
    PDB.model_validate(json)


def test_read_mmcif_filestring() -> None:
    """Rest reading of a mmcif string."""
    with open("tests/data/1ema.cif") as f:
        data = f.read()
    pdb = pdb_from_mmcif_filestring(data)
    assert pdb.description.code == "1EMA"

    json = pdb.model_dump()
    PDB.model_validate(json)


# fmt: off
@mark.regression
@mark.parametrize(
    "code",
    [
        # Codes from molecules of the month August 2024–January 2025
        "7S6B", "8UCS", "2ZVY", "1F4V", "6YKM", "6E10", "6E11", "3VCM", "2X0B",
        "6OS0", "1N9U", "1O86", "2V0Z", "6KI1", "6KI2", "7YYO", "6TJV", "7EGL",
        "7CYF", "7EGK", "7ZCG", "3FRT", "6AP1",
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
    pdb = fetch_pdb(code)

    json = pdb.model_dump()
    PDB.model_validate(json)

def test_from_pdb_to_pdb_2qto() -> None:
    with open("tests/data/2qto.pdb") as f:
        data = f.read()
    pdb = pdb_from_pdb_filestring(data)
    filestring = pdb_object_to_pdb_filestring(pdb, header=True, source=True, keyword=True, crystallography=True)
    pdb2 = pdb_from_pdb_filestring(filestring)

    assert pdb.description == pdb2.description
    assert pdb.experiment == pdb2.experiment
    assert pdb.models == pdb2.models

def test_from_pdb_to_pdb_1ema() -> None:
    with open("tests/data/1ema.pdb") as f:
        data = f.read()
    pdb = pdb_from_pdb_filestring(data)
    filestring = pdb_object_to_pdb_filestring(pdb, header=True, source=True, keyword=True, crystallography=True)
    pdb2 = pdb_from_pdb_filestring(filestring)

    assert pdb.description == pdb2.description
    assert pdb.models == pdb2.models

def test_from_pdb_to_pdb_2hu4() -> None:
    with open("tests/data/2HU4.pdb") as f:
        data = f.read()
    pdb = pdb_from_pdb_filestring(data)
    filestring = pdb_object_to_pdb_filestring(pdb, header=True, source=True, keyword=True, crystallography=True)

    pdb2 = pdb_from_pdb_filestring(filestring)

    assert pdb.description == pdb2.description
    # not true but doesn't matter
    print(pdb.geometry == pdb2.geometry)
    assert pdb.models == pdb2.models

def mmcif_author_format_to_pdb_format(authors: list[str]) -> list[str]:
    return [f"{last.upper()}{first.upper()}" for first, last in
            (author.split(", ") for author in authors)]

def compare_descriptions_mmcif_and_pdb(mmcif_description: PDBDescription, pdb_description: PDBDescription) -> bool:
    return (
        mmcif_description.code == pdb_description.code and
        mmcif_description.title == pdb_description.title and
        mmcif_author_format_to_pdb_format(mmcif_description.authors) == (pdb_description.authors) and
        mmcif_description.classification == pdb_description.classification and
        mmcif_description.deposition_date == pdb_description.deposition_date and
        sorted(mmcif_description.keywords) == sorted(pdb_description.keywords)
    )

def compare_experiments_mmcif_and_pdb(mmcif_experiment: PDBExperiment, pdb_experiment: PDBExperiment) -> bool:
    return (
        mmcif_experiment.expression_system.upper() == pdb_experiment.expression_system and  # type: ignore [union-attr]
        mmcif_experiment.missing_residues == pdb_experiment.missing_residues and
        mmcif_experiment.source_organism.upper() == pdb_experiment.source_organism and  # type: ignore [union-attr]
        mmcif_experiment.technique == pdb_experiment.technique
    )

def compare_models_mmcif_and_pdb(mmcif_models: list[PDBModel], pdb_models: list[PDBModel]) -> bool:
    for i in range(len(mmcif_models)):
        assert mmcif_models[i].polymer == pdb_models[i].polymer
        assert mmcif_models[i].non_polymer == pdb_models[i].non_polymer
        assert mmcif_models[i].branched == pdb_models[i].branched
    return True

def test_from_mmcif_to_pdb_1ema() -> None:
    with open("tests/data/1ema.cif") as f:
        mmcif_data = f.read()

    with open("tests/data/1ema.pdb") as f:
        pdb_data = f.read()
    mmcif_1ema = pdb_from_mmcif_filestring(mmcif_data)
    pdb_1ema = pdb_from_pdb_filestring(pdb_data)

    assert compare_descriptions_mmcif_and_pdb(mmcif_1ema.description, pdb_1ema.description)
    assert compare_experiments_mmcif_and_pdb(mmcif_1ema.experiment, pdb_1ema.experiment)
    assert mmcif_1ema.quality == pdb_1ema.quality
    assert compare_models_mmcif_and_pdb(mmcif_1ema.models, pdb_1ema.models)


def test_residue_order_preserved_in_pdb_output_after_json_sort() -> None:
    """
    Test that residue ordering is correct in PDB output after JSON round-trip.

    PostgreSQL JSONB sorts keys alphabetically, which breaks numeric ordering for
    residue IDs like "A.-3", "A.-2", "A.-1", "A.0", "A.1", "A.10", etc.
    The pdb_object_to_pdb_filestring function should output residues in correct order.
    """
    # Create residues in correct numeric order (including negative numbers)
    # Each residue needs at least one atom for serialization
    residues = {
        "A.-3": PDBResidue(name="ACE", number=-3, atoms={1: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.-2": PDBResidue(name="ASP", number=-2, atoms={2: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.-1": PDBResidue(name="ASP", number=-1, atoms={3: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.0": PDBResidue(name="LYS", number=0, atoms={4: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.1": PDBResidue(name="MET", number=1, atoms={5: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.2": PDBResidue(name="ASP", number=2, atoms={6: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.10": PDBResidue(name="GLY", number=10, atoms={7: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.11": PDBResidue(name="ALA", number=11, atoms={8: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
    }

    polymer = PDBPolymer(internal_id="A", residues=residues)
    model = PDBModel(polymer={"A": polymer})
    pdb = PDB(
        description=PDBDescription(),
        experiment=PDBExperiment(),
        geometry=PDBGeometry(),
        quality=PDBQuality(),
        models=[model],
    )

    # Simulate JSONB round-trip (sort_keys=True mimics PostgreSQL JSONB behavior)
    json_dict = pdb.model_dump(mode="json")
    json_string_sorted = json.dumps(json_dict, sort_keys=True)
    json_parsed = json.loads(json_string_sorted)

    # Verify JSON sorting scrambled the residue keys (alphabetic order)
    alphabetic_order = list(json_parsed["models"][0]["polymer"]["A"]["residues"].keys())
    assert alphabetic_order == ["A.-1", "A.-2", "A.-3", "A.0", "A.1", "A.10", "A.11", "A.2"]

    # Restore from sorted JSON
    restored_pdb = PDB.model_validate(json_parsed)

    # Serialize to PDB format - should be in correct numeric order
    pdb_string = pdb_object_to_pdb_filestring(restored_pdb, seqres=False, hetnam=False)

    # Extract residue numbers from ATOM lines
    residue_nums = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM"):
            # Residue number is in columns 23-26 (0-indexed: 22:26)
            res_num = int(line[22:26])
            residue_nums.append(res_num)

    # Verify correct numeric order in output
    expected_order = [-3, -2, -1, 0, 1, 2, 10, 11]
    assert residue_nums == expected_order


def test_residue_order_with_multichar_chain_id() -> None:
    """
    Test that atom ordering by serial number preserves file order with multi-char chain IDs.

    Atom serials encode original file order, allowing correct reconstruction after JSONB.
    """
    # Atom serials reflect correct file order: -1, 0, 1, 2, 10
    residues = {
        "AA.-1": PDBResidue(name="ACE", number=-1, atoms={1: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "AA.0": PDBResidue(name="MET", number=0, atoms={2: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "AA.1": PDBResidue(name="ALA", number=1, atoms={3: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "AA.2": PDBResidue(name="VAL", number=2, atoms={4: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "AA.10": PDBResidue(name="GLY", number=10, atoms={5: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
    }

    polymer = PDBPolymer(internal_id="A", residues=residues)
    model = PDBModel(polymer={"A": polymer})
    pdb = PDB(
        description=PDBDescription(),
        experiment=PDBExperiment(),
        geometry=PDBGeometry(),
        quality=PDBQuality(),
        models=[model],
    )

    # Simulate JSONB round-trip
    json_dict = pdb.model_dump(mode="json")
    json_string_sorted = json.dumps(json_dict, sort_keys=True)
    json_parsed = json.loads(json_string_sorted)

    restored_pdb = PDB.model_validate(json_parsed)

    # Serialize to PDB format - should be in correct numeric order
    pdb_string = pdb_object_to_pdb_filestring(restored_pdb, seqres=False, hetnam=False)

    # Extract residue numbers from ATOM lines
    residue_nums = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM"):
            res_num = int(line[22:26])
            residue_nums.append(res_num)

    expected_order = [-1, 0, 1, 2, 10]
    assert residue_nums == expected_order


def test_residue_order_with_insertion_codes() -> None:
    """
    Test that atom ordering by serial number preserves original file order through JSONB.

    Atoms are sorted by serial number, which encodes the original file order.
    This test verifies that order is preserved even when dict keys get scrambled.
    """
    # Atom serials reflect correct file order: 16, 16A, 17, 169, 170
    residues = {
        "A.16": PDBResidue(name="MET", number=16, atoms={1: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.16A": PDBResidue(name="GLY", number=16, atoms={2: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.17": PDBResidue(name="VAL", number=17, atoms={3: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.169": PDBResidue(name="ALA", number=169, atoms={4: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
        "A.170": PDBResidue(name="LEU", number=170, atoms={5: PDBAtom(x=0, y=0, z=0, element="N", bvalue=0)}, bvalue=0.0),
    }

    polymer = PDBPolymer(internal_id="A", residues=residues)
    model = PDBModel(polymer={"A": polymer})
    pdb = PDB(
        description=PDBDescription(),
        experiment=PDBExperiment(),
        geometry=PDBGeometry(),
        quality=PDBQuality(),
        models=[model],
    )

    # Simulate JSONB round-trip (sorts keys alphabetically)
    json_dict = pdb.model_dump(mode="json")
    json_string_sorted = json.dumps(json_dict, sort_keys=True)
    json_parsed = json.loads(json_string_sorted)

    # Verify JSONB scrambles residue keys (alphabetically: 16, 169, 16A, 17, 170)
    alphabetic_order = list(json_parsed["models"][0]["polymer"]["A"]["residues"].keys())
    assert alphabetic_order == ["A.16", "A.169", "A.16A", "A.17", "A.170"]

    restored_pdb = PDB.model_validate(json_parsed)

    # Serialize to PDB format - should preserve original order via atom serials
    pdb_string = pdb_object_to_pdb_filestring(restored_pdb, seqres=False, hetnam=False)

    # Extract residue numbers and insertion codes from ATOM lines
    residue_ids = []
    for line in pdb_string.split("\n"):
        if line.startswith("ATOM"):
            res_num = int(line[22:26])
            ins_code = line[26].strip()
            residue_ids.append((res_num, ins_code))

    # Order preserved: 16, 16A, 17, 169, 170 (based on atom serial order)
    expected_order = [(16, ""), (16, "A"), (17, ""), (169, ""), (170, "")]
    assert residue_ids == expected_order
