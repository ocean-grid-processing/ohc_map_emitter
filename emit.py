#!/usr/bin/env python3
"""OHCA map packaging: one ohc_derive `map` blob -> the gridded per-cell deliverable.

The factory has already done the analysis — the cross-layer mask, the n_fac combine over constituents,
and the per-cell anomaly against the baseline window. Its blob carries `map` (with `map_sd` when the
ensemble is on) as the combined per-cell OHC density anomaly in TJ/m^2, on the native monthly grid,
plus `area_m2`, the `level`, and the `time_window` it was built with. This step is only the packaging:
scale TJ/m^2 -> J/m^2, lay the grid out in the ME4OH order, relabel, and write one file per level.

Output: `ohca(LONGITUDE, LATITUDE, TIME)` in J/m^2 (the baseline-referenced anomaly), plus `ohca_std`
(per-cell ensemble 1-sigma) when the members were present. The grid keeps its native monthly TIME axis,
re-encoded to days-since-1900 on write (the loader decodes it, so we pin the units back on).

Target: the gridded map has no Zenodo reference of its own; its per-cell referencing is the same
baseline-mean removal the OHCA series uses, applied before the integral instead of after — so
area-weighting the map reproduces the OHCA series (postflight_integral_check.py).
"""
import argparse
import json
import os

import xarray as xr

TERA = 1e12                 # TJ/m^2 -> J/m^2

# This step's identity, used to key its block inside the consolidated `config_record`. Each output file
# is built from one derive blob (one level), so this step is a 1-in-1-out courier: it rolls the blob's
# whole provenance chain forward and folds its own block in, emitting the lot as one `config_record`
# attribute (one attribute keeps the file in HDF5 compact storage — see `stamp_config_record`).
STAGE = "ohc_map_emitter"
_PROV_SUFFIXES = ("_run_config", "_run_facts", "_code_version")


def _compact(obj):
    """One-line JSON — reads as a single clean line in `ncdump -h`."""
    return json.dumps(obj, separators=(",", ":"), default=str)


def _shared_and_per(group_map):
    """{group: block} -> (shared, per): keys present in every group with an equal value go to `shared`;
    everything else stays per group. Lossless — block[g] == {**shared, **per[g]}."""
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
    """Factor one fan-out `{group: value}`: object values -> shared + per_<axis> (a fully-shared block
    collapses to the bare shared object); scalar values -> the bare value if all agree, else per_<axis>."""
    values = list(block.values())
    if all(isinstance(v, dict) for v in values):
        shared, per = _shared_and_per(block)
        if not any(per.values()):
            return shared
        return {"shared": shared, "per_" + axis: per}
    if all(v == values[0] for v in values):
        return values[0]
    return {"per_" + axis: dict(block)}


def _maybe_json(v):
    """Parse a forwarded block back to JSON so it nests as a real object; leave non-JSON (a bare
    code_version URL) as-is."""
    try:
        return json.loads(v)
    except (TypeError, ValueError):
        return v


def _dry_constituent_fanouts(record):
    """Within the assembled record, factor any per-constituent fan-out into shared + per_constituent
    (lossless). Driven by the constituents roster in `ohc_derive.run_facts`, so only genuine fan-outs
    are touched — a value like `n_fac`, nested inside a non-fanned block, is never a candidate."""
    facts = record.get("ohc_derive", {}).get("run_facts")
    roster = (set(facts["constituents"]) if isinstance(facts, dict)
              and isinstance(facts.get("constituents"), list) else None)
    if not roster:
        return record
    for parts in record.values():
        if not isinstance(parts, dict):
            continue
        for name, val in list(parts.items()):
            if isinstance(val, dict) and len(val) > 1 and set(val) <= roster:
                parts[name] = _compact_block(val, "constituent")
    return record


def stamp_config_record(out, blob, cfg, source_path):
    """Assemble the whole provenance chain into ONE `config_record` attribute, keyed by stage and DRY'd
    per constituent. A single attribute keeps the file at <=8 global attributes, i.e. HDF5 *compact*
    attribute storage — which every reader handles. Emitting a dozen separate `*_run_config` etc. tips
    HDF5 into dense (fractal-heap) storage, whose exact layout some netcdf builds mis-read."""
    record = {}
    # the forwarded chain: group the blob's *_run_config/_run_facts/_code_version by stage
    for k, v in blob.attrs.items():
        for suffix in _PROV_SUFFIXES:
            if k.endswith(suffix):
                record.setdefault(k[:-len(suffix)], {})[suffix[1:]] = _maybe_json(v)
                break
    # this step's own block
    record[STAGE] = {
        "run_config": vars(cfg),
        "run_facts": {
            "level": blob.attrs.get("level"),
            "time_window": blob.attrs.get("time_window", "all"),
            "quantities_present": [q for q in ("map",) if q in blob],
            "ensemble": "map_sd" in blob,
            "source_blob": os.path.abspath(source_path),
        },
        "code_version": cfg.code_version,
    }
    out.attrs["config_record"] = _compact(_dry_constituent_fanouts(record))


def _me4oh(da):
    """(time, lat, lon) -> (LONGITUDE, LATITUDE, TIME), the ME4OH submission layout."""
    return da.transpose("lon", "lat", "time").rename(
        {"lon": "LONGITUDE", "lat": "LATITUDE", "time": "TIME"})


def build_dataset(blob, tag, provenance_link):
    """A derive `map` blob -> the gridded OHCA anomaly deliverable Dataset, ME4OH layout."""
    window = blob.attrs.get("time_window", "all")
    baseline = "all-time mean" if window == "all" else "%s mean" % window

    anom = _me4oh(blob["map"] * TERA)                        # TJ/m^2 -> J/m^2, per-cell anomaly
    anom.attrs = {"units": "J/m2",
                  "long_name": "ocean heat content anomaly (%s removed)" % baseline}
    dv = {"ohca": anom}
    if "map_sd" in blob:
        sd = _me4oh(blob["map_sd"] * TERA)
        sd.attrs = {"units": "J/m2", "long_name": "ohca ensemble standard deviation (1-sigma)"}
        dv["ohca_std"] = sd

    out = xr.Dataset(dv)
    for coord, src in (("LONGITUDE", "lon"), ("LATITUDE", "lat")):   # carry coord attrs the loader kept
        if src in blob.coords and blob[src].attrs:
            out[coord].attrs = dict(blob[src].attrs)

    out.attrs["level"] = blob.attrs["level"]
    out.attrs["time_window"] = window
    out.attrs["provenance_tag"] = tag
    if provenance_link is not None:
        out.attrs["provenance_link"] = provenance_link
    return out


def _data_span(blob):
    """`YYYY_YYYY` for the years the blob's own axis spans (`year` annual, `time` monthly)."""
    if "year" in blob.coords:
        yrs = blob["year"].values.astype(int)
        return "%d_%d" % (int(yrs.min()), int(yrs.max()))
    if "time" in blob.coords:
        yrs = blob["time"].values.astype("datetime64[Y]").astype(int) + 1970
        return "%d_%d" % (int(yrs.min()), int(yrs.max()))
    return "all"


def _file_token(blob):
    """Combined filename token `<data>_tw<baseline>`: the blob's own data span, then its baseline
    window (defaulting to the whole data span when the derive run was windowless) — matching the derive
    filename, so runs differing only in baseline don't collide here either."""
    data = _data_span(blob)
    win = blob.attrs.get("time_window", "all")
    baseline = data if (not win or win == "all") else win.replace("-", "_")
    return "%s_tw%s" % (data, baseline)


def filename(level, tag, token):
    """Per-level map name: ohca_map_<tag>_<lo>_<hi>_dbar_<data>_tw<baseline>.nc (tag leads, after the step)."""
    low, high = level.split("_")
    return "ohca_map_%s_%s_%s_dbar_%s.nc" % (tag, low, high, token)


def main():
    ap = argparse.ArgumentParser(description="OHCA map packaging: ohc_derive `map` blob -> gridded deliverable")
    ap.add_argument("blobs", nargs="+", help="ohc_derive output NetCDFs built with --quantities map")
    ap.add_argument("--tag", required=True, help="provenance tag: filename token + provenance_tag attr")
    ap.add_argument("--provenance-link", default=None, help="URL/path to the provenance record")
    ap.add_argument("--code-version", required=True,
                    help="URL to the exact ohc_map_emitter code (commit/release); stamped in config_record")
    ap.add_argument("--out", default=".")
    cfg = ap.parse_args()
    os.makedirs(cfg.out, exist_ok=True)
    for path in cfg.blobs:
        blob = xr.open_dataset(path)
        if "map" not in blob:
            raise SystemExit("%s carries no map; run ohc_derive with --quantities map" % path)
        dest = os.path.join(cfg.out, filename(blob.attrs["level"], cfg.tag, _file_token(blob)))
        out = build_dataset(blob, cfg.tag, cfg.provenance_link)
        stamp_config_record(out, blob, cfg, path)                  # whole chain -> one config_record attr
        enc = {v: {"_FillValue": -999.0} for v in out.data_vars}   # target fill (NaN -> -999)
        enc["TIME"] = {"units": "days since 1900-01-01", "calendar": "proleptic_gregorian"}
        out.to_netcdf(dest, engine="netcdf4", encoding=enc)
        print("wrote", dest)


if __name__ == "__main__":
    main()
