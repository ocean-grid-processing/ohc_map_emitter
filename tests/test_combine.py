"""Gridded combine: n_fac nansum, nan-if-all-nan, whole-record anomaly, absolute-member spread."""
import numpy as np

import aggregate
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
