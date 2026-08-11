import os
import sys

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def make_layer(tag, field, ens=None):
    """read_layer()-shaped dict. `field` is (time, lat, lon); `ens` (member, time, lat, lon) or None."""
    f = xr.DataArray(np.asarray(field, dtype="float64"), dims=("time", "lat", "lon"))
    e = None if ens is None else xr.DataArray(np.asarray(ens, dtype="float64"),
                                              dims=("member", "time", "lat", "lon"))
    return {"tag": tag, "field": f, "ens": e, "src": xr.Dataset()}
