"""Combined-layer configuration for ohc_map_emitter.

Same shallowest-first `n_fac` combination as the GCOS / OHCA-OHU emitters, but this product combines
on the **grid** (not integrated series), so only `n_fac` and the contributor tags are used here — the
`dz` column (volume weight) is irrelevant to a per-area map and is ignored.

NOTE: this is the *third* verbatim copy of this table (also in ohc_gcos_emitter and
ohc_ohca_ohu_emitter). Three consumers trips the rule-of-three, so a shared module is now earned —
but the emitters are separate git repos, so extracting it is its own task. Kept duplicated for the
deadline; if you edit one, edit all three.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Contributor:
    tag: str      # mapped-layer tag "top_bottom", matches ohc_ingest/publish `mapped_layer`
    n_fac: int    # multiplier (thin-layer proxy scaling)
    dz: float     # contributor thickness in m; unused by the map emitter (kept for table parity)


@dataclass(frozen=True)
class Level:
    name: str
    contributors: tuple

    @property
    def low(self):
        return int(self.name.split("_")[0])

    @property
    def high(self):
        return int(self.name.split("_")[1])


# Verbatim with the GCOS / OHCA-OHU tables (the third copy). `--levels` picks the subset a run emits.
LEVELS = [
    Level("0_300", (
        Contributor("15_20", 3, 5),
        Contributor("15_300", 1, 285),
    )),
    Level("0_700", (
        Contributor("15_20", 3, 5),
        Contributor("15_300", 1, 285),
        Contributor("300_700", 1, 400),
    )),
    Level("0_1000", (                       # needs the 700_1000 mapped layer
        Contributor("15_20", 3, 5),
        Contributor("15_300", 1, 285),
        Contributor("300_700", 1, 400),
        Contributor("700_1000", 1, 300),
    )),
    Level("700_2000", (
        Contributor("700_1850", 1, 1150),
        Contributor("1800_1850", 3, 50),
    )),
    Level("0_2000", (
        Contributor("15_20", 3, 5),
        Contributor("15_300", 1, 285),
        Contributor("300_700", 1, 400),
        Contributor("700_1850", 1, 1150),
        Contributor("1800_1850", 3, 50),
    )),
]

LEVELS_BY_NAME = {lv.name: lv for lv in LEVELS}


def select_levels(names=None):
    if names is None:
        return list(LEVELS)
    out = []
    for n in names:
        if n not in LEVELS_BY_NAME:
            raise SystemExit("unknown level %r; known: %s" % (n, list(LEVELS_BY_NAME)))
        out.append(LEVELS_BY_NAME[n])
    return out


def required_tags(levels):
    tags = []
    for lv in levels:
        for c in lv.contributors:
            if c.tag not in tags:
                tags.append(c.tag)
    return tags
