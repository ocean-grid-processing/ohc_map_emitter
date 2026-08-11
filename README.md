# ohc_map_emitter

`ohc_map_emitter` combines mapped-layer gridded OHC into **combined-layer OHC anomaly maps** — one
NetCDF per level (`0_700`, `0_1000`, `0_2000`), ME4OH layout, monthly, per grid cell.

```
ohc_ingest ─▶ publish (--preset wmo --ensemble) ─▶ ohc_map_emitter ─▶ per-level map .nc
```

Unlike the GCOS / OHCA-OHU emitters, this one stays **on the grid** — it never integrates — so it
reads the publish `OHC_`/`OHCENS_` output **directly** rather than going through `ohc_derive` (derive
collapses the grid to series).

## What it computes

Per combined layer, on the grid:

```
combined(x,t) = Σᵢ n_facᵢ · fieldᵢ(x,t)     nansum: a missing (NaN) contributor counts as 0,
                                            and a cell where EVERY contributor is NaN → NaN
ohca(x,t)     = combined(x,t) − mean_t combined(x)         [J/m², whole-record anomaly]
ohca_std(x,t) = std over members of the combined ABSOLUTE field   [J/m², 1-sigma]
```

- **`nansum`, not `sum`.** A plain sum would NaN-propagate and collapse the footprint to the deepest
  layer (bottom-must-be-wet); `nansum` treats a missing deep contributor as zero, mirroring the
  series combine's skipna sum and landing on the GCOS area/volume footprint. The one masking step is
  restoring NaN where *all* contributors are missing (else `nansum` reads land as a spurious 0). No
  shallowest-reference policy — the union footprint equals the shallowest by construction.
- **Anomaly is the whole-record mean** removed per cell (`OHC − OHC_time_mean`), not a windowed
  baseline.
- **Spread is of the absolute field**, offset included — same choice as the OHCA emitter
  (`data_yearly_std`), not the demeaned anomaly. So (as there) the value is the anomaly and the
  `_std` is the spread of the absolute; faithful, not obviously symmetric.
- **Units J/m²** (density) — the map is per-area; multiplying out by cell area is left to the
  consumer.

Masking is inherited from the `OHC_`/`OHCENS_` files: publish with `--preset wmo --ensemble` gives
the mean∪members footprint (the `ensemble_incomplete` bit active), so the maps' NaNs match the mask
the GCOS file uses for areas/volumes.

## Combined layers (config)

[`layers.py`](layers.py) — same table as the other emitters (the third copy; see the note there).
Only `n_fac` and the contributor tags are used (a map has no volume, so `dz` is ignored). Bands:
`0_300`, `0_700`, `0_1000`, `700_2000`, `0_2000`. `--levels` selects a subset.

## Usage

### Test

Basic unit tests run locally in a container:
   
```bash
docker image build -t ohc_map_emitter:test .
docker container run -v $(pwd):/app ohc_map_emitter:test pytest
```

### Run
```bash
python combine.py OHC_*.nc --tag "OHC-maps 2026 <run>" [--levels ...] [--out DIR]
```
The `OHC_*.nc` are publish submissions built with `--preset wmo --ensemble`; the `OHCENS_` siblings
must sit beside them for the spread.

## Memory

The spread materializes the combined `(member, time, lat, lon)` cube — ~7–14 GB per level (f32/f64).
Fine on a big node; if tight, stream member-by-member (Welford) in `aggregate.combine_level_maps`
instead of building the cube. Compute is I/O-bound (reading the `OHCENS_` files), tens of minutes,
not the compute-week an early estimate wrongly implied.

## Open / first-pass details

- **Variable names** (`ohca`, `ohca_std`) and the **filename** are our choice — no target file to
  match, unlike OHCA/OHU. Reconcile if a map submission spec turns up.
- **"J/m²"** confirmed with the boss (density, not per-cell integrated Joules).
- **Depth-partial cells** (surface data missing but a deep layer present) get a value under plain
  union-restore; that pathological case is deliberately left un-editorialized (net-new, deferred).
