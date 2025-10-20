from .base import Base, LowercaseStrEnum


class Solvent(LowercaseStrEnum):
    WATER = "water"
    NITROMETHANE = "nitromethane"
    NITROBENZENE = "nitrobenzene"
    TOLUENE = "toluene"
    BENZENE = "benzene"
    CHLOROBENZENE = "chlorobenzene"
    CARBONTETRACHLORIDE = "carbontetrachloride"
    DICHLOROETHANE = "dichloroethane"
    DICHLOROMETHANE = "dichloromethane"
    CHLOROFORM = "chloroform"
    DIETHYLETHER = "diethylether"
    DIISOPROPYLETHER = "diisopropylether"
    DIMETHYLSULFOXIDE = "dimethylsulfoxide"
    TETRAHYDROFURAN = "tetrahydrofuran"
    CYCLOHEXANE = "cyclohexane"
    OCTANE = "octane"
    ACETICACID = "aceticacid"
    HEXANE = "hexane"
    ETHYLACETATE = "ethylacetate"
    ACETONE = "acetone"
    ACETONITRILE = "acetonitrile"
    METHANOL = "methanol"
    ETHANOL = "ethanol"
    ISOPROPANOL = "isopropanol"
    DIMETHYLACETAMIDE = "dimethylacetamide"
    DIMETHYLFORMAMIDE = "dimethylformamide"
    N_METHYLPYRROLIDONE = "n_methylpyrrolidone"
    ETHYLENE_GLYCOL = "ethylene_glycol"


class SolventModel(LowercaseStrEnum):
    PCM = "pcm"
    CPCM = "cpcm"
    ALPB = "alpb"
    COSMO = "cosmo"
    GBSA = "gbsa"
    CPCMX = "cpcmx"


class SolventSettings(Base):
    solvent: Solvent
    model: SolventModel
