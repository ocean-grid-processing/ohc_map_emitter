"""Write one combined-layer anomaly map as an ME4OH-layout NetCDF.

`ohca(LONGITUDE, LATITUDE, TIME)` in J/m^2 (baseline-referenced anomaly), plus `ohca_std` (per-cell
ensemble 1-sigma) when the ensemble was present. Coordinates and their attrs are carried from the
source OHC_ submission, so the grid and time axis match the ME4OH submissions exactly.

Global attributes are kept minimal (`level`, `period`) so that with the run's `provenance_tag` /
`provenance_link` and the one consolidated `config_record`, the file stays at <=8 global attributes,
i.e. HDF5 compact attribute storage — every reader handles it. The mask/preset/cp0/rho0/source details
that used to sit loose here now live inside `config_record` (carried in the forwarded `localgp_*`
blocks), so nothing is lost.
"""
import xarray as xr


def _me4oh(da):
    """(time, lat, lon) -> (LONGITUDE, LATITUDE, TIME), the ME4OH submission layout."""
    return da.transpose("lon", "lat", "time").rename(
        {"lon": "LONGITUDE", "lat": "LATITUDE", "time": "TIME"})


def build_level_dataset(cl, baseline):
    """`baseline` is the anomaly-baseline label for the long_name (e.g. "all-time" or "2005-2024")."""
    src = cl["src"]

    anom = _me4oh(cl["anom"])
    anom.attrs = {"units": "J/m2",
                  "long_name": "ocean heat content anomaly (%s mean removed)" % baseline}
    dv = {"ohca": anom}
    if cl["sd"] is not None:
        sd = _me4oh(cl["sd"])
        sd.attrs = {"units": "J/m2", "long_name": "ohca ensemble standard deviation (1-sigma)"}
        dv["ohca_std"] = sd

    out = xr.Dataset(dv)
    for c in ("LONGITUDE", "LATITUDE", "TIME"):          # carry coord attrs (units, axis, calendar)
        if c in src.coords:
            out[c].attrs = dict(src[c].attrs)

    out.attrs["level"] = cl["name"]
    if "period" in src.attrs:                            # the data span (kept top-level as discoverable)
        out.attrs["period"] = src.attrs["period"]
    return out


def filename(cl, tag, window):
    # `tag` is already whitespace-sanitized by the CLI; used verbatim (no lowercasing/munging).
    return "ohca_map_%d_%d_dbar_%s_%s.nc" % (cl["low"], cl["high"], window, tag)
