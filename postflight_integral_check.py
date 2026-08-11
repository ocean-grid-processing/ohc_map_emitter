#!/usr/bin/env python3
"""Postflight: prove the gridded map integrates back to the Zenodo-validated OHCA series.

Ties this product to a known-good path. The map's per-cell anomaly and the OHCA emitter's
`integral_anom` share one referencing convention (subtract the whole-record mean); because
demeaning and the area integral are both linear and the mask is time-constant, they commute, so:

    Σ_x  ohca_map(x, t) · A(x)   ==   1e12 · integral_anom(t)            [per month]
    yearly-mean of that / area_m2  ==   ohca_series(year)               [J/m^2, the OHCA file]

So area-weighting the map, annualizing (calendar-year mean, as the OHCA emitter does), and dividing
by the OHCA file's own `area_m2` must reproduce that file's `ohca` to machine precision. Any real
gap means the two products drifted apart (different layers/n_fac, mask, cell area, or time base) —
this script localizes that.

`A(x)` is regenerated from the map's grid with the SAME pure formula ohc_derive/loader uses
(_cell_area, R = 6_371_000 m), so the check needs nothing but the two emitters' output files.

    python postflight_integral_check.py --map ohca_map_*.nc --ohca ohca_ohu_*.nc [--rtol 1e-8]

Files are paired by their `level` attr; every level present in both is checked. Exits non-zero if
any level/year exceeds tolerance.

Requires: numpy, xarray>=2024.10, netCDF4.
"""
import argparse
import glob
import sys

import numpy as np
import xarray as xr

EARTH_RADIUS_M = 6_371_000.0   # identical to ohc_derive/loader._cell_area


def _cell_area(lat, lon):
    """Spherical cell area [lat, lon] in m^2 from coordinate spacing — verbatim with derive's loader."""
    lat = np.asarray(lat, dtype="float64")
    lon = np.asarray(lon, dtype="float64")
    dlat = abs(lat[1] - lat[0])
    dlon = abs(lon[1] - lon[0])
    s_hi = np.sin(np.deg2rad(lat + dlat / 2))
    s_lo = np.sin(np.deg2rad(lat - dlat / 2))
    area_lat = EARTH_RADIUS_M ** 2 * np.deg2rad(dlon) * (s_hi - s_lo)
    grid = np.repeat(area_lat[:, None], len(lon), axis=1)
    return xr.DataArray(grid, dims=("LATITUDE", "LONGITUDE"),
                        coords={"LATITUDE": lat, "LONGITUDE": lon})


def _by_level(paths):
    out = {}
    for p in paths:
        ds = xr.open_dataset(p, decode_times=True)
        lvl = ds.attrs.get("level")
        if lvl is None:
            raise SystemExit("%s has no `level` attr — not an emitter output?" % p)
        out[str(lvl)] = (p, ds)
    return out


def _years_of(coord):
    """Calendar years for a time coord, whether xarray decoded it to datetime or left it as
    `days since <epoch>` floats."""
    if np.issubdtype(np.asarray(coord.values).dtype, np.datetime64):
        return coord.dt.year.values.astype(int)
    units = coord.attrs.get("units", "")
    if "since" not in units:
        raise SystemExit("time coord %r is numeric but has no 'days since' units to decode" % coord.name)
    epoch = np.datetime64(units.split("since", 1)[1].strip().split()[0])
    stamps = epoch + coord.values.astype("timedelta64[D]")
    return stamps.astype("datetime64[Y]").astype(int) + 1970


def _map_integral_yearly(mds):
    """Area-weighted spatial integral of the map's ohca, as a calendar-year-mean series [J]."""
    ohca = mds["ohca"].astype("float64")                       # (LON, LAT, TIME) J/m^2
    area = _cell_area(mds["LATITUDE"].values, mds["LONGITUDE"].values)
    integ = (ohca * area).sum(("LATITUDE", "LONGITUDE"))       # skipna: NaN cells drop; already J
    yrs = xr.DataArray(_years_of(mds["TIME"]), dims=("TIME",), coords={"TIME": integ["TIME"]})
    # Calendar-year mean, matching ohc_ohca_ohu_emitter.emit._yearly (groupby("time.year").mean()).
    return integ.groupby(yrs.rename("year")).mean("TIME")      # (year,) J


def _ohca_years(ods):
    """OHCA series as (year -> J/m^2), decoding time_ohca (labelled at each year's 1-June)."""
    ohca = ods["ohca"].astype("float64")
    years = _years_of(ods["time_ohca"])
    area_m2 = float(ohca.attrs["area_m2"])
    return {int(y): float(v) for y, v in zip(years, ohca.values)}, area_m2


def check_level(level, mpath, mds, opath, ods, rtol, atol):
    map_int = _map_integral_yearly(mds)                         # (year,) J
    ohca_by_year, area_m2 = _ohca_years(ods)                    # {year: J/m^2}, m^2

    rows, worst = [], 0.0
    for y in map_int["year"].values.astype(int):
        if y not in ohca_by_year:
            continue
        implied = float(map_int.sel(year=y)) / area_m2          # map -> J/m^2
        target = ohca_by_year[y]
        denom = max(abs(target), 1e-300)
        rel = abs(implied - target) / denom
        worst = max(worst, rel)
        rows.append((y, implied, target, implied - target, rel))

    ok = all(abs(i - t) <= atol + rtol * abs(t) for _, i, t, _, _ in rows)
    print("\nlevel %s   map=%s   ohca=%s" % (level, mpath.split("/")[-1], opath.split("/")[-1]))
    print("  area_m2 = %.6e   years compared = %d   worst rel = %.2e   %s"
          % (area_m2, len(rows), worst, "PASS" if ok else "FAIL"))
    print("    year   map∫/area (J/m^2)     ohca (J/m^2)         abs diff       rel diff")
    for y, i, t, d, r in rows:
        print("    %4d   %+.9e   %+.9e   %+.3e   %.3e" % (y, i, t, d, r))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", nargs="+", required=True, help="ohca_map_*.nc (one or more, any level)")
    ap.add_argument("--ohca", nargs="+", required=True, help="ohca_ohu_*.nc (one or more, any level)")
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--atol", type=float, default=0.0, help="J/m^2 absolute floor (default 0)")
    args = ap.parse_args()

    maps = _by_level(sum([glob.glob(g) for g in args.map], []) or args.map)
    ohcas = _by_level(sum([glob.glob(g) for g in args.ohca], []) or args.ohca)

    shared = [lv for lv in maps if lv in ohcas]
    if not shared:
        raise SystemExit("no shared level between map files %s and ohca files %s"
                         % (list(maps), list(ohcas)))

    all_ok = True
    for lv in shared:
        mpath, mds = maps[lv]
        opath, ods = ohcas[lv]
        all_ok &= check_level(lv, mpath, mds, opath, ods, args.rtol, args.atol)

    print("\n%s — %d level(s) checked" % ("ALL PASS" if all_ok else "FAILURES PRESENT", len(shared)))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
