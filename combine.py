#!/usr/bin/env python3
"""Combine mapped-layer gridded OHC (publish OHC_/OHCENS_ output) into combined-layer anomaly maps.

    python combine.py OHC_*.nc --tag "OHC-maps-2026-OP20260127b" \
        [--levels 0_700,0_1000,0_2000] [--collaborators STR] [--out DIR]

Each OHC_*.nc is one mapped layer's publish submission; the OHCENS_ sibling (from
`publish.py --ensemble`) beside it supplies the per-cell spread. Reads the grid directly — this
product never integrates, so it does NOT go through ohc_derive. One .nc is written per combined
level, ME4OH layout, with `-999` fill.

Requires: numpy, xarray>=2024.10, netCDF4.
"""
import argparse
import os

import layers as layers_mod
import aggregate
import emit


def _sanitize_tag(tag):
    """Strip all whitespace from a provenance tag; never lowercase or otherwise munge it — it must
    match the provenance record char-for-char."""
    return "".join(tag.split())


def _add_provenance(ds, args):
    """Stamp provenance attrs: provenance_tag is the required --tag (also the filename run token)."""
    ds.attrs["provenance_tag"] = args.tag
    if args.provenance_link is not None:
        ds.attrs["provenance_link"] = args.provenance_link


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("submissions", nargs="+", help="OHC_ submissions, one per mapped layer")
    ap.add_argument("--tag", required=True, help='run tag, e.g. "OHC-maps-2026-OP20260127b"')
    ap.add_argument("--levels", default=None,
                    help="comma list of combined levels to emit (default: all in layers.py)")
    ap.add_argument("--collaborators", default="LocalGP by Giglio, Sukianto, Kuusela, Mills")
    ap.add_argument("--provenance-link", default=None,
                    help="URL/path to the provenance record; written to the provenance_link header attr")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    args.tag = _sanitize_tag(args.tag)

    names = None if not args.levels else [s.strip() for s in args.levels.split(",")]
    levels = layers_mod.select_levels(names)

    by_tag = {}
    for nc in args.submissions:
        layer = aggregate.read_layer(nc)
        by_tag[layer["tag"]] = layer

    need = layers_mod.required_tags(levels)
    absent = [t for t in need if t not in by_tag]
    if absent:
        raise SystemExit("missing contributor submission(s): %s (needed by %s)"
                         % (absent, [lv.name for lv in levels]))

    spread = aggregate.uncertainty_available(need, by_tag)      # raises on a partial ensemble mix

    os.makedirs(args.out, exist_ok=True)
    written = []
    for lv in levels:
        cl = aggregate.combine_level_maps(lv, by_tag)
        ds = emit.build_level_dataset(cl, args.tag, args.collaborators)
        _add_provenance(ds, args)
        path = os.path.join(args.out, emit.filename(cl, args.tag))
        enc = {v: {"_FillValue": -999.0} for v in ds.data_vars}   # target fill (NaN -> -999)
        ds.to_netcdf(path, engine="netcdf4", encoding=enc)
        written.append(path)

    print("wrote %d map(s):" % len(written))
    for p in written:
        print("  ", p)
    print("spread:", "on — ohca_std written" if spread else "off (no OHCENS_ siblings)")


if __name__ == "__main__":
    main()
