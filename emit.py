"""Write one combined-layer anomaly map as an ME4OH-layout NetCDF.

`ohca(LONGITUDE, LATITUDE, TIME)` in J/m^2 (whole-record anomaly), plus `ohca_std` (per-cell ensemble
1-sigma) when the ensemble was present. Coordinates and their attrs are carried from the source OHC_
submission, so the grid and time axis match the ME4OH submissions exactly.
"""
import xarray as xr


def _me4oh(da):
    """(time, lat, lon) -> (LONGITUDE, LATITUDE, TIME), the ME4OH submission layout."""
    return da.transpose("lon", "lat", "time").rename(
        {"lon": "LONGITUDE", "lat": "LATITUDE", "time": "TIME"})


def build_level_dataset(cl, tag, collaborators):
    src = cl["src"]

    anom = _me4oh(cl["anom"])
    anom.attrs = {"units": "J/m2", "long_name": "ocean heat content anomaly (all-time mean removed)"}
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
    out.attrs["description"] = "%s, %s" % (tag, collaborators)
    for k in ("mask_preset", "mask_applied", "cp0", "rho0", "source", "product", "period"):
        if k in src.attrs:
            out.attrs[k] = src.attrs[k]
    return out


def filename(cl, tag):
    return "ohca_map_%d_%d_dbar_%s.nc" % (cl["low"], cl["high"], tag.lower().replace(" ", ""))
