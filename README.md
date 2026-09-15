# ohc_map_emitter

`ohc_map_emitter` packages one `ohc_derive` blob into a **combined-layer OHC anomaly map** — one
NetCDF per level, `ohca_map_<tag>_<lo>_<hi>_dbar_<data>_tw<baseline>.nc`, ME4OH layout, monthly, per grid cell.

```
ohc_ingest ─▶ publish ─▶ ohc_derive (--quantities map) ─▶ ohc_map_emitter ─▶ per-level map .nc
```

The analysis is all upstream. `ohc_derive` does the cross-layer mask, the `n_fac` combine over
constituents, the per-cell anomaly against the baseline window, and the ensemble → SD collapse — the
map is just another derive deliverable (`--quantities map`), the gridded sibling of `ohca`: where
`ohca` area-integrates the masked field, `map` carries it through on the grid. Its blob hands over
`map` (with `map_sd` when the ensemble was on) as the combined per-cell OHC density anomaly in TJ/m²,
on the native monthly grid, plus `area_m2`, the `level`, and the `time_window` it was built with. This
emitter is only the packaging: scale to J/m², lay the grid out in ME4OH order, relabel, and write.

> **Units:** **`ohca` in J/m²** (density; the map is per-area, so multiplying out by cell area is left
> to the consumer). Conversion is `ohca = map[TJ/m²] × 1e12`; the `_std` scales the same way.

## What it computes

Per blob (one synthetic level):

```
ohca(x, t)     = blob.map    × 1e12      # J/m², the baseline-referenced per-cell anomaly
ohca_std(x, t) = blob.map_sd × 1e12      # J/m², present when the ensemble was on
```

The grid keeps its **native monthly TIME axis** from the submissions (re-encoded to days-since-1900 on
write — the loader decodes it, so the units are pinned back). `LONGITUDE`/`LATITUDE` and their attrs
carry straight from the blob. The `low`/`high` in the filename come from the blob's `level` attr; the
anomaly baseline label comes from its `time_window`.

**Spread convention.** The `_std` is derive's ensemble spread of the map deliverable — the member
spread of the **anomaly**, `n_fac` worst-case summed across constituents — i.e. the same collapse that
produces `ohca_sd` in the series. So the map's uncertainty matches the OHCA series' by construction; the
value is the anomaly and the `_std` is its ensemble 1-sigma.

## Building the input

`ohc_map_emitter` consumes one `ohc_derive` blob per synthetic level, built with `--quantities map`:

```bash
python ../ohc_derive/run.py OHC_<constituents>.nc \
    --level 0_2000 --bathy etopo60.nc --quantities map \
    --tag <tag> --code-version URL --out <dir>
```

`--quantities map` can also ride along with the series quantities in one run (`--quantities
ohca,ohu,map`) so the mask pass is shared — the same blob then feeds both this emitter and
`ohc_ohca_ohu_emitter`, each reading only the variables it packages. A dedicated `map` run keeps the
gridded member cube (≈ 7–14 GB) isolated from the series jobs; either way works.

The anomaly baseline is set **at the derive step** by `--time-window` (default: whole record, matching
the OHCA series; pass `2005:2024` to match the GCOS convention). Each blob **must** carry:

- data var **`map`** (per-cell, monthly, TJ/m²), and its **`map_sd`** companion when the derive run
  kept the ensemble. The emitter errors if `map` is absent, and skips `ohca_std` when `map_sd` isn't
  present.
- attrs **`level`** and **`time_window`** (which becomes the anomaly baseline label).

## Usage

### Test
```bash
docker image build -t ohc_map_emitter:test .
docker container run -v $(pwd):/app ohc_map_emitter:test pytest
```

### Run
```bash
python emit.py derive_<tag>_<data>_tw<baseline>_<level>.nc [more levels …] --tag <tag> --code-version URL \
    --project LocalGP --author Giglio_etal2026 --citation "…" [--provenance-link URL] [--out DIR]
```

`--project` / `--author` become the filename's trailing pair (`…_<project>_<author>.nc`) and are recorded
in `config_record`; `--citation` is written to a standalone top-level `citation` attribute.

One blob in, one map out, per level. `-999` fill on every data variable on write (NaN off-footprint
lands as `-999` on disk and decodes back to NaN on read).

**Provenance chain.** Each map is built from one derive blob, so this step is a 1-in-1-out courier: it
rolls that blob's whole provenance chain forward (every `*_run_config` / `*_run_facts` /
`*_code_version` — the grouped `localgp_ingest_*` / `localgp_publish_*` and the `ohc_derive_*` blocks)
and folds in its own block, emitting the lot as **one** `config_record` attribute keyed by stage:

```
config_record = {
  "localgp_ingest":  {"run_config": {…}, "run_facts": {…}, "code_version": "…"},
  "localgp_publish": {…},
  "ohc_derive":      {…},
  "ohc_map_emitter": {"run_config": {resolved args}, "run_facts": {level, window, …}, "code_version": "…"}
}
```

*Why one attribute:* a dozen separate global attributes tips HDF5 into **dense (fractal-heap) attribute
storage**, whose exact layout some netcdf builds mis-read; a single attribute keeps the file at ≤ 8
global attributes (`level`, `time_window`, `provenance_tag`, `provenance_link`, `config_record`), i.e.
**compact** storage, which every reader handles. The global `provenance_tag` / `provenance_link` stay
separate (the run's discoverable identity).

The per-constituent `localgp_*` blocks are **DRY'd**: each `{15_20:{…}, 15_300:{…}, …}` fan-out is
factored into `{"shared": {common config}, "per_constituent": {only what differs}}`, and a fan-out
whose entries fully agree collapses to a bare value. Lossless and reversible, driven by the
`constituents` roster in `ohc_derive.run_facts` — so only genuine fan-outs are touched and a value like
`n_fac`, nested inside a non-fanned block, is never mistaken for one.

#### emit.py options

| option | default | effect |
|---|---|---|
| `derive_*.nc` (positional, 1+) | *(required)* | `ohc_derive` blobs, one per synthetic level (`derive_<tag>_<data>_tw<baseline>_<level>.nc`). Each must carry `map`. |
| `--tag` | *(required)* | run token in the filename (`ohca_map_<tag>_<lo>_<hi>_dbar_<data>_tw<baseline>_<project>_<author>.nc`) and the `provenance_tag` attr. Used verbatim; should match the tag the blob was derived under. |
| `--provenance-link` | *(none)* | URL/path to the provenance record; written to the `provenance_link` attr. |
| `--code-version` | *(required)* | URL to the exact ohc_map_emitter code (commit/release); stamped inside `config_record`. |
| `--project` | *(required)* | project string; first of the filename's trailing pair (whitespace-stripped, case preserved), a standalone top-level `project` attr, and recorded in `config_record`. |
| `--author` | *(required)* | author string; last of the filename's trailing pair (e.g. `Giglio_etal2026`) and recorded in `config_record`. |
| `--citation` | *(required)* | citation sentence; written to the standalone top-level `citation` attr (kept out of `config_record` so it isn't duplicated). |
| `--out` | `.` | output directory (created if absent). |

## Validation

The gridded map has no Zenodo reference of its own, but it ties back to one that does. Its per-cell
referencing is the same baseline-mean removal as the OHCA series (validated against Zenodo 14720478),
and — because it and the series come off the **same masked field** in derive (`integral` is the
area-weighted sum of the `map` primitive, and demeaning/annualizing are linear over that sum) —
area-weighting the map and annualizing reproduces the OHCA file's `ohca` series **by construction**:

```bash
python postflight_integral_check.py --map ohca_map_*.nc --ohca ohca_ohu_*.nc
```

Match the baselines (both whole-record, or both the same `--time-window`) when cross-checking. The
check is a regression guard rather than a discrepancy hunt now — the footprint can't drift, because
there's a single mask upstream instead of a second copy here.

## Memory

With `--quantities map`, derive materializes the combined `(realization, time, lat, lon)` cube — ≈ 7–14
GB per level (f32/f64). Fine on a big node; run `map` as its own derive job to keep it off the series
jobs, or fold it into the series run to share the mask pass at the cost of holding the grid alongside.

## Open / first-pass details

- **Variable names** (`ohca`, `ohca_std`) and the **filename** are our choice — no target map file to
  match, unlike OHCA/OHU. Reconcile if a map submission spec turns up.
- **"J/m²"** confirmed (density, not per-cell integrated Joules).
