#!/usr/bin/env python3
"""Combine mapped-layer gridded OHC (publish OHC_/OHCENS_ output) into combined-layer anomaly maps.

    python combine.py OHC_*.nc --tag "OHC-maps-2026-OP20260127b" --code-version URL \
        [--time-window 2005:2024] [--levels 0_700,0_1000,0_2000] [--collaborators STR] [--out DIR]

Each OHC_*.nc is one mapped layer's publish submission; the OHCENS_ sibling (from
`publish.py --ensemble`) beside it supplies the per-cell spread. Reads the grid directly — this
product never integrates, so it does NOT go through ohc_derive. One .nc is written per combined
level, ME4OH layout, with `-999` fill.

`--time-window` sets the per-cell anomaly baseline (default: whole record), like the series emitters.
The whole provenance chain (the `localgp_*` blocks on the submissions, plus this step's own) is folded
into one `config_record` attribute, so each file stays in HDF5 compact attribute storage.

Requires: numpy, xarray>=2024.10, netCDF4.
"""
import argparse
import json
import os

import numpy as np

import layers as layers_mod
import aggregate
import emit

STAGE = "ohc_map_emitter"
_PROV_SUFFIXES = ("_run_config", "_run_facts", "_code_version")


def _sanitize_tag(tag):
    """Strip all whitespace from a provenance tag; never lowercase or otherwise munge it — it must
    match the provenance record char-for-char."""
    return "".join(tag.split())


def _parse_window(s):
    """YEAR0:YEAR1 -> (int, int); None/empty -> None. Separator `:`, `-`, or `_` (so the filename token
    2004_2025 works too)."""
    if not s:
        return None
    y0, y1 = (int(x) for x in s.replace("-", ":").replace("_", ":").split(":"))
    return (y0, y1)


def _window_token(args, sample_field):
    """Filename year-range token: the `--time-window` years if given, else the record span (from any
    contributor's days-since-1900 time axis)."""
    if args.time_window:
        return "%d_%d" % args.time_window
    t = np.round(sample_field["time"].values).astype("timedelta64[D]") + np.datetime64("1900-01-01")
    yrs = t.astype("datetime64[Y]").astype(int) + 1970
    return "%d_%d" % (int(yrs.min()), int(yrs.max()))


def _compact(obj):
    """One-line JSON — reads as a single clean line in `ncdump -h`."""
    return json.dumps(obj, separators=(",", ":"), default=str)


def _maybe_json(v):
    try:
        return json.loads(v)
    except (TypeError, ValueError):
        return v


def _shared_and_per(group_map):
    groups = list(group_map)
    common = set(group_map[groups[0]])
    for g in groups[1:]:
        common &= set(group_map[g])
    shared = {}
    for k in sorted(common):
        vals = [group_map[g][k] for g in groups]
        if all(v == vals[0] for v in vals):
            shared[k] = vals[0]
    per = {g: {k: v for k, v in group_map[g].items() if k not in shared} for g in groups}
    return shared, per


def _compact_block(block, axis):
    values = list(block.values())
    if all(isinstance(v, dict) for v in values):
        shared, per = _shared_and_per(block)
        if not any(per.values()):
            return shared
        return {"shared": shared, "per_" + axis: per}
    if all(v == values[0] for v in values):
        return values[0]
    return {"per_" + axis: dict(block)}


def _stamp_config_record(ds, level, by_tag, args):
    """Consolidate the whole chain into one `config_record` attribute (keeps the file in HDF5 compact
    storage). Forwards each contributor's `localgp_*` blocks grouped by constituent and DRY'd, then adds
    this step's own block. `provenance_tag` / `provenance_link` stay separate as the run's identity."""
    ds.attrs["provenance_tag"] = args.tag
    if args.provenance_link is not None:
        ds.attrs["provenance_link"] = args.provenance_link

    record = {}                                              # stage -> part -> {constituent: block}
    for c in level.contributors:
        for k, v in by_tag[c.tag]["src"].attrs.items():
            for suffix in _PROV_SUFFIXES:
                if k.endswith(suffix):
                    record.setdefault(k[:-len(suffix)], {}).setdefault(suffix[1:], {})[c.tag] = _maybe_json(v)
                    break
    roster = set(c.tag for c in level.contributors)          # this step supplies its own roster
    for parts in record.values():
        for name, val in list(parts.items()):
            if isinstance(val, dict) and len(val) > 1 and set(val) <= roster:
                parts[name] = _compact_block(val, "constituent")
    record[STAGE] = {
        "run_config": vars(args),
        "run_facts": {
            "level": level.name,
            "contributors": [c.tag for c in level.contributors],
            "n_fac": {c.tag: c.n_fac for c in level.contributors},
            "time_window": ("%d-%d" % args.time_window) if args.time_window else "all",
            "ensemble": "ohca_std" in ds.data_vars,
            "source_submissions": [os.path.abspath(by_tag[c.tag]["path"]) for c in level.contributors],
        },
        "code_version": args.code_version,
    }
    ds.attrs["config_record"] = _compact(record)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("submissions", nargs="+", help="OHC_ submissions, one per mapped layer")
    ap.add_argument("--tag", required=True, help='run tag, e.g. "OHC-maps-2026-OP20260127b"')
    ap.add_argument("--levels", default=None,
                    help="comma list of combined levels to emit (default: all in layers.py)")
    ap.add_argument("--time-window", default=None,
                    help="YEAR0:YEAR1 anomaly baseline window (default: whole record); separator "
                         "`:`, `-`, or `_`")
    ap.add_argument("--collaborators", default="LocalGP by Giglio, Sukianto, Kuusela, Mills")
    ap.add_argument("--provenance-link", default=None,
                    help="URL/path to the provenance record; written to the provenance_link header attr")
    ap.add_argument("--code-version", required=True,
                    help="URL to the exact ohc_map_emitter code (commit/release); stamped in config_record")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    args.tag = _sanitize_tag(args.tag)
    args.time_window = _parse_window(args.time_window)

    names = None if not args.levels else [s.strip() for s in args.levels.split(",")]
    levels = layers_mod.select_levels(names)

    by_tag = {}
    for nc in args.submissions:
        layer = aggregate.read_layer(nc)
        tag = layer["tag"]
        if tag in by_tag:                                       # one file per native level (no last-wins)
            raise SystemExit("two submissions map to native level %s:\n  %s\n  %s\n"
                             "the pool must hold exactly one file per native level."
                             % (tag, by_tag[tag]["path"], nc))
        layer["path"] = nc
        by_tag[tag] = layer

    need = layers_mod.required_tags(levels)
    absent = [t for t in need if t not in by_tag]
    if absent:
        raise SystemExit("missing contributor submission(s): %s (needed by %s)"
                         % (absent, [lv.name for lv in levels]))

    spread = aggregate.uncertainty_available(need, by_tag)      # raises on a partial ensemble mix

    window_token = _window_token(args, by_tag[next(iter(need))]["field"])
    baseline = ("%d-%d" % args.time_window) if args.time_window else "all-time"

    os.makedirs(args.out, exist_ok=True)
    written = []
    for lv in levels:
        cl = aggregate.combine_level_maps(lv, by_tag, window=args.time_window)
        ds = emit.build_level_dataset(cl, baseline)
        _stamp_config_record(ds, lv, by_tag, args)
        path = os.path.join(args.out, emit.filename(cl, args.tag, window_token))
        enc = {v: {"_FillValue": -999.0} for v in ds.data_vars}   # target fill (NaN -> -999)
        ds.to_netcdf(path, engine="netcdf4", encoding=enc)
        written.append(path)

    print("wrote %d map(s):" % len(written))
    for p in written:
        print("  ", p)
    print("spread:", "on — ohca_std written" if spread else "off (no OHCENS_ siblings)")


if __name__ == "__main__":
    main()
