from typing import TypedDict

ConfigType = TypedDict(
    "ConfigType",
    {
        "doMP2": bool,
        "doMP3": bool,
        "doCCSD": bool,
        "doTriples": bool,
        "doSFX": bool,
        "toFreeze": bool,
        "doChengBasis": bool,
        "basisLibKey": str,
        "toPrintDensity": bool,
    },
)

GEN_CONFIG: ConfigType = dict(
    doMP2=True,
    doMP3=False,
    doCCSD=True,
    doTriples=True,
    doSFX=True,
    toFreeze=True,
    doChengBasis=False,
    basisLibKey="6-311+G(3df)",
    toPrintDensity=False,
)

Q_CONFIG: ConfigType = dict(**GEN_CONFIG)
Q_CONFIG["doMP3"] = False

PENTUPLE_CONFIG: ConfigType = dict(**GEN_CONFIG)
PENTUPLE_CONFIG.update(
    {
        "doMP3": False,
        "doCCSD": False,
        "doTriples": False,
    }
)
TEST_CONFIG: ConfigType = dict(
    doMP2=False,
    doMP3=False,
    doCCSD=False,
    doTriples=False,
    doSFX=True,
    toFreeze=True,
    doChengBasis=False,
    basisLibKey="6-311+G(3df)",
    toPrintDensity=False,
)
LOCAL_CONFIG: ConfigType = dict(
    doMP2=True,
    doMP3=False,
    doCCSD=False,
    doTriples=False,
    doSFX=True,
    toFreeze=True,
    doChengBasis=False,
    basisLibKey="6-311+G(3df)",
    toPrintDensity=False,
)

MP3_CONFIG: ConfigType = dict(
    doMP2=True,
    doMP3=True,
    doCCSD=False,
    doTriples=False,
    doSFX=True,
    toFreeze=True,
    doChengBasis=False,
    basisLibKey="6-311+G(3df)",
    toPrintDensity=True,
)
