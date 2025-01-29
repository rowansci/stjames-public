from typing import Annotated, Any, Optional

from pydantic import AfterValidator

from ..base import Base, round_float, round_optional_float
from ..constraint import Constraint
from ..method import Method
from ..mode import Mode
from ..solvent import Solvent
from .workflow import Workflow


class ConformerSettings(Base):
    num_confs_considered: int = 100
    num_confs_taken: int = 50
    rmsd_cutoff: float = 0.1

    transition_state: bool = False
    final_method: Method = Method.AIMNET2_WB97MD3
    solvent: Optional[Solvent] = Solvent.WATER
    max_energy: float = 5

    constraints: list[Constraint] = []


class RdkitConformerSettings(ConformerSettings):
    num_initial_confs: int = 100
    max_mmff_energy: float = 10
    max_mmff_iterations: int = 500


class CrestConformerSettings(ConformerSettings):
    flags: str = "--quick --ewin 10"
    gfn: int | str = "ff"


class Conformer(Base):
    energy: Annotated[float, AfterValidator(round_float(6))]
    weight: Annotated[float | None, AfterValidator(round_optional_float(6))] = None

    # uuid, optionally
    uuid: Optional[str] = None


class ConformerWorkflow(Workflow):
    mode: Mode = Mode.RAPID
    settings: ConformerSettings = ConformerSettings()
    conformers: list[Conformer] = []

    def model_post_init(self, __context: Any) -> None:
        self.settings = csearch_settings_by_mode(self.mode, self.settings)


def csearch_settings_by_mode(mode: Mode, old_settings: Optional[ConformerSettings] = None) -> ConformerSettings:
    if mode == Mode.MANUAL:
        assert old_settings is not None
        return old_settings

    settings: ConformerSettings

    if mode == Mode.METICULOUS:
        settings = CrestConformerSettings(
            gfn=2,
            flags="--ewin 15 --noreftopo",
            max_energy=10,
            num_confs_considered=500,
            num_confs_taken=150,
        )

    elif mode == Mode.CAREFUL:
        settings = CrestConformerSettings(
            gfn="ff",
            flags="--quick --ewin 10 --noreftopo",
            num_confs_considered=150,
            num_confs_taken=50,
        )

    elif mode == Mode.RAPID or Mode.AUTO:
        settings = RdkitConformerSettings(
            num_initial_confs=300,
            max_mmff_energy=15,
            num_confs_considered=100,
            num_confs_taken=50,
        )

    elif mode == Mode.RECKLESS:
        settings = RdkitConformerSettings(
            num_initial_confs=200,
            max_mmff_energy=10,
            num_confs_considered=50,
            num_confs_taken=20,
            rmsd_cutoff=0.25,
        )

    else:
        raise ValueError(f"invalid mode ``{mode.value}`` for conformer settings")

    if old_settings is not None:
        settings.final_method = old_settings.final_method
        settings.solvent = old_settings.solvent
        settings.constraints = old_settings.constraints

    return settings
