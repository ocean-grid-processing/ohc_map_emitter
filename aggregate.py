"""Read gridded OHC (the publish OHC_/OHCENS_ output) and combine mapped layers into combined-layer
anomaly maps + per-cell ensemble spread.

Unlike the series emitters, this stays on the grid, so it reads the **publish** output directly
(derive integrates the grid away). Per mapped layer:

    OHC_...nc     DATA(LONGITUDE, LATITUDE, TIME)          posterior-mean OHC density, TJ/m^2
    OHCENS_...nc  DATA(MEMBER, LONGITUDE, LATITUDE, TIME)  the 100 members (for the spread)

both already masked (NaN) by the publish preset (`wmo` = mean∪members footprint when the ensemble
is on). Combination is the shallowest-first `n_fac` sum, on the grid, with the agreed rule:

    combined(x, t) = nansum_i [ n_fac_i * field_i(x, t) ]     # missing contributor -> 0
    then NaN where EVERY contributor is NaN (so all-missing cells don't read as 0).

This mirrors the series combine's skipna sum and lands on the GCOS area/volume footprint, with no
shallowest-reference policy. The **anomaly** subtracts each cell's whole-record mean; the **spread**
is the ensemble std of the ABSOLUTE combined field (offset included, matching the OHCA emitter's
data_yearly_std choice), not of the demeaned anomaly.
"""
import os

import xarray as xr

TERA = 1e12   # TJ/m^2 -> J/m^2


def _to_tlatlon(da):
    return da.transpose("TIME", "LATITUDE", "LONGITUDE").rename(
        {"TIME": "time", "LATITUDE": "lat", "LONGITUDE": "lon"})


def _ohcens_path(oc_nc):
    base = os.path.basename(oc_nc)
    if not base.startswith("OHC_"):
        return None
    return os.path.join(os.path.dirname(oc_nc), "OHCENS_" + base[len("OHC_"):])


def read_layer(nc):
    """Read one mapped layer's gridded mean field (+ its OHCENS_ ensemble sibling, if present).

    Returns tag, field(time, lat, lon) [TJ/m^2], ens(member, time, lat, lon) or None, and `src`
    (the opened OHC_ Dataset, for its LON/LAT/TIME coords to carry to the output).
    """
    ds = xr.open_dataset(nc, decode_times=False)
    if "DATA" not in ds.data_vars:
        raise SystemExit("%s has no DATA variable — expected an OHC_ submission" % nc)
    tag = ds.attrs.get("mapped_layer") or ds.attrs.get("layer_m")
    if not tag or "_" not in str(tag):
        raise SystemExit("%s has no usable mapped_layer/layer_m attr (got %r)" % (nc, tag))

    field = _to_tlatlon(ds["DATA"]).astype("float64")            # (time, lat, lon) TJ/m^2

    ens = None
    ens_nc = _ohcens_path(nc)
    if ens_nc and os.path.exists(ens_nc):
        eda = xr.open_dataset(ens_nc, decode_times=False)["DATA"]
        ens = (eda.transpose("MEMBER", "TIME", "LATITUDE", "LONGITUDE")
                  .rename({"MEMBER": "member", "TIME": "time", "LATITUDE": "lat", "LONGITUDE": "lon"})
                  .astype("float64"))                            # (member, time, lat, lon)

    return {"tag": str(tag), "field": field, "ens": ens, "src": ds}


def _nansum(das, nfacs):
    """nansum of n_fac*da over contributors; NaN where EVERY contributor is NaN.

    Broadcasts over any leading member axis, so the same call combines the mean field and the
    ensemble. Missing (NaN) contributors count as 0; a cell with no valid contributor is restored
    to NaN so all-missing cells don't read as a spurious 0.
    """
    total = None
    all_missing = None
    for da, n in zip(das, nfacs):
        total = (n * da).fillna(0.0) if total is None else total + (n * da).fillna(0.0)
        miss = da.isnull()
        all_missing = miss if all_missing is None else (all_missing & miss)
    return total.where(~all_missing)


def combine_level_maps(level, by_tag):
    """Combine one Level's contributors on the grid.

    Returns {name, low, high, anom(time, lat, lon) [J/m^2], sd(time, lat, lon) [J/m^2] or None}.
    """
    missing = [c.tag for c in level.contributors if c.tag not in by_tag]
    if missing:
        raise SystemExit("level %s: missing contributor input(s) %s" % (level.name, missing))
    contribs = [by_tag[c.tag] for c in level.contributors]
    nfacs = [c.n_fac for c in level.contributors]

    combined = _nansum([c["field"] for c in contribs], nfacs)        # (time, lat, lon) TJ/m^2
    anom = (combined - combined.mean("time")) * TERA                 # J/m^2, whole-record anomaly

    sd = None
    if all(c["ens"] is not None for c in contribs):
        # NOTE (memory): accumulates the combined (member, time, lat, lon) cube — ~7-14 GB. Fine on a
        # big node; if tight, stream member-by-member (Welford) instead of materializing the cube.
        combined_ens = _nansum([c["ens"] for c in contribs], nfacs)  # (member, time, lat, lon)
        sd = combined_ens.std("member", ddof=1) * TERA              # (time, lat, lon) J/m^2, absolute

    return {"name": level.name, "low": level.low, "high": level.high, "anom": anom, "sd": sd,
            "src": contribs[0]["src"]}                              # coords from any contributor


def uncertainty_available(needed, by_tag):
    """True if every needed layer carries an ensemble, False if none do; raise on a partial mix."""
    have = [t for t in needed if by_tag[t]["ens"] is not None]
    if not have:
        return False
    if len(have) != len(needed):
        missing = [t for t in needed if by_tag[t]["ens"] is None]
        raise SystemExit("ensemble present for %s but missing for %s — publish those layers with "
                         "--ensemble, or none of the maps get a spread." % (have, missing))
    return True
