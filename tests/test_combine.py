"""Gridded combine: n_fac nansum, nan-if-all-nan, whole-record/windowed anomaly, absolute-member spread."""
import json
import types

import numpy as np
import xarray as xr

import aggregate
import combine
from layers import Contributor, Level
from conftest import make_layer

TERA = 1e12
L = Level("0_700", (Contributor("15_20", 3, 5),
                    Contributor("15_300", 1, 285),
                    Contributor("300_700", 1, 400)))


def _uniform_ramp(nt, scale):
    """(nt, 1, 2) field: scale*(1..nt), constant in space, ramping in time."""
    r = scale * (1.0 + np.arange(nt, dtype="float64"))
    return r[:, None, None] * np.ones((nt, 1, 2))


def test_nansum_weighting_and_anomaly():
    nt = 3
    by = {
        "15_20":  make_layer("15_20",  _uniform_ramp(nt, 1.0)),
        "15_300": make_layer("15_300", _uniform_ramp(nt, 2.0)),
        "300_700": make_layer("300_700", _uniform_ramp(nt, 4.0)),
    }
    out = aggregate.combine_level_maps(L, by)
    combined = 3 * _uniform_ramp(nt, 1.0) + _uniform_ramp(nt, 2.0) + _uniform_ramp(nt, 4.0)
    expect = (combined - combined.mean(axis=0)) * TERA          # whole-record anomaly, J/m^2
    assert out["anom"].dims == ("time", "lat", "lon")
    assert np.allclose(out["anom"].values, expect)
    assert out["sd"] is None                                    # no ensemble supplied


def test_missing_contributor_counts_as_zero():
    # 300_700 missing (NaN) at lon=1 for all time -> nansum treats it as 0 there
    nt = 3
    c = _uniform_ramp(nt, 4.0).copy()
    c[:, 0, 1] = np.nan
    by = {
        "15_20":  make_layer("15_20",  _uniform_ramp(nt, 1.0)),
        "15_300": make_layer("15_300", _uniform_ramp(nt, 2.0)),
        "300_700": make_layer("300_700", c),
    }
    out = aggregate.combine_level_maps(L, by)
    combined = 3 * _uniform_ramp(nt, 1.0) + _uniform_ramp(nt, 2.0)   # lon=1: 300_700 dropped
    combined[:, 0, 0] += _uniform_ramp(nt, 4.0)[:, 0, 0]             # lon=0: 300_700 present
    expect = (combined - combined.mean(axis=0)) * TERA
    assert np.allclose(out["anom"].values, expect)                  # both cells finite, deep=0 at lon=1
    assert np.all(np.isfinite(out["anom"].values))


def test_nan_where_all_contributors_missing():
    nt = 2
    def with_nan_col0(scale):
        a = _uniform_ramp(nt, scale)
        a[:, 0, 0] = np.nan          # lon=0 missing in every contributor
        return a
    by = {
        "15_20":  make_layer("15_20",  with_nan_col0(1.0)),
        "15_300": make_layer("15_300", with_nan_col0(2.0)),
        "300_700": make_layer("300_700", with_nan_col0(4.0)),
    }
    out = aggregate.combine_level_maps(L, by)
    assert np.all(np.isnan(out["anom"].values[:, 0, 0]))   # all-missing cell -> NaN (not 0)
    assert np.all(np.isfinite(out["anom"].values[:, 0, 1]))


def test_windowed_baseline_vs_whole_record():
    ref = np.datetime64("1900-01-01")
    days = np.array([float((np.datetime64("%d-07-01" % y) - ref) / np.timedelta64(1, "D"))
                     for y in (2004, 2005, 2006)])
    field = xr.DataArray((np.array([10.0, 2.0, 3.0])[:, None, None]) * np.ones((3, 1, 1)),
                         dims=("time", "lat", "lon"), coords={"time": ("time", days)})
    lv = Level("0_20", (Contributor("15_20", 1, 5),))
    by = {"15_20": {"tag": "15_20", "field": field, "ens": None, "src": xr.Dataset()}}

    win = aggregate.combine_level_maps(lv, by, window=(2005, 2006))    # baseline = mean(2,3) = 2.5
    assert np.allclose(win["anom"].values[:, 0, 0], np.array([7.5, -0.5, 0.5]) * TERA)
    whole = aggregate.combine_level_maps(lv, by, window=None)          # baseline = mean(10,2,3) = 5
    assert np.allclose(whole["anom"].values[:, 0, 0], np.array([5.0, -3.0, -2.0]) * TERA)


def test_config_record_consolidates_chain(tmp_path):
    lv = Level("0_700", (Contributor("15_20", 3, 5), Contributor("15_300", 1, 285)))

    def _src(tag):
        d = xr.Dataset()
        d.attrs.update({
            "localgp_ingest_run_config": json.dumps({"var_name": "pt", "cp0": 3989.0, "dir": "/" + tag}),
            "localgp_publish_code_version": "https://x/commit/p",   # agrees across constituents
        })
        return d

    by = {"15_20":  {"tag": "15_20",  "src": _src("15_20"),  "path": "/a.nc"},
          "15_300": {"tag": "15_300", "src": _src("15_300"), "path": "/b.nc"}}
    ds = xr.Dataset({"ohca": (("t",), [1.0])})
    args = types.SimpleNamespace(tag="M", provenance_link="http://m", code_version="https://x/commit/m",
                                 time_window=None, levels=None, collaborators="me", out=str(tmp_path),
                                 submissions=["/a.nc", "/b.nc"])
    combine._stamp_config_record(ds, lv, by, args)

    assert "config_record" in ds.attrs
    assert not any(k.endswith(("_run_config", "_run_facts", "_code_version")) for k in ds.attrs)
    rec = json.loads(ds.attrs["config_record"])

    ig = rec["localgp_ingest"]["run_config"]                          # DRY'd per constituent
    assert ig["shared"] == {"cp0": 3989.0, "var_name": "pt"}
    assert set(ig["per_constituent"]) == {"15_20", "15_300"}
    assert rec["localgp_publish"]["code_version"] == "https://x/commit/p"   # agreeing scalar -> bare
    assert rec["ohc_map_emitter"]["code_version"].endswith("/m")            # this step's own
    assert rec["ohc_map_emitter"]["run_facts"]["contributors"] == ["15_20", "15_300"]
    assert rec["ohc_map_emitter"]["run_facts"]["n_fac"] == {"15_20": 3, "15_300": 1}


def test_spread_is_member_std_of_absolute_combined():
    nt, nm = 2, 4
    offs = np.array([-1.0, 0.0, 1.0, 2.0])[:, None, None, None]
    ens_1520 = np.full((nm, nt, 1, 1), 10.0) + offs               # members differ by a constant offset
    zeros = np.zeros((nm, nt, 1, 1))
    L1 = Level("0_700", (Contributor("15_20", 3, 5),))            # single contributor for a clean check
    by = {"15_20": make_layer("15_20", np.full((nt, 1, 1), 10.0), ens=ens_1520)}
    out = aggregate.combine_level_maps(L1, by)
    # combined member = 3 * (10 + offs); std over members = 3 * std(offs) (the level cancels)
    expect = 3.0 * np.std([-1.0, 0.0, 1.0, 2.0], ddof=1) * TERA
    assert np.allclose(out["sd"].values, expect)
