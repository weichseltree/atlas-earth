"""Make the three arrays this repository publishes (land, residents, heights) from their sources, on one
latitude-longitude grid of 2048 x 4096 cells: level 3 of the tile pyramid, about 9.8 km at the equator.

It runs once, on the machine that holds the sources (SirBase, in MAP-L1's download cache of weichseltree/chronon),
with a Python that has rasterio, netCDF4, pyogrio and shapely. The Action never runs it: it builds from results/.

    python tools/prepare.py [<folder with the sources>]     # writes results/earth/*.npy and results/earth/summary.json

| array      | source                                                   | each cell holds                            |
|------------|----------------------------------------------------------|--------------------------------------------|
| land       | OpenStreetMap land polygons of 2026-01-03, EPSG:3857     | the share of its area that is not ocean    |
| people     | GHS-POP R2023A, epoch 2025, 30 arc seconds               | residents per km2                          |
| height     | ETOPO 2022, 60 arc seconds, surface                      | mean height above sea level in metres      |

Land: the polygons are drawn onto a Mercator grid of 16384 x 16384 samples (a sample is land when its centre is
inside a polygon) and every lat-lon cell takes the share of the samples whose centres fall into it. The polygons end at
85.05 degrees; north of that there is only sea (share 0), south of it only Antarctica (share 1).
Residents: each 30" cell's count goes to the lat-lon cell holding its centre, so the world's total is kept; the sum is
then divided by the cell's area on a sphere of the Earth's mean radius.
Heights: the mean of the 60" cells whose centres fall into the lat-lon cell.
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

try:  # progress on SirBase's board, when run through exprun; nothing anywhere else
    import expmetrics
except ImportError:
    expmetrics = None

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "earth"
H, W = 2048, 4096
R_EARTH = 6_371_008.8
HALF = 20037508.342789244  # half the side of the Web Mercator square, metres
MERC = 16384  # Mercator samples per side for the land polygons
BATCH = 100_000  # polygons drawn at a time (all at once needs ~12 GB)

SOURCES = {
    "land": "osm-e8c4bde2f90a7312/land-polygons-split-3857/land_polygons.shp",
    "people": "ghs-cec59a706722493b/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif",
    "height": "ETOPO_2022_v1_60s_N90W180_surface.nc",
}
ARCHIVES = {  # the downloads the sources were unpacked from, as MAP-L1 pinned them
    "land": "land-polygons-split-3857.zip",
    "people": "GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip",
    "height": "ETOPO_2022_v1_60s_N90W180_surface.nc",
}
PEOPLE_TOTAL = 8_191_988_536  # the source's own total, as MAP-L1 measured it


def report(**kw):
    if expmetrics is not None:
        expmetrics.report(**kw)


def say(*a):
    print(*a, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def rows_of(lat):
    """The grid row of each latitude (degrees); row 0 is the north."""
    return np.clip(np.floor((90.0 - np.asarray(lat, dtype=np.float64)) / 180.0 * H).astype(np.int64), 0, H - 1)


def cols_of(lon):
    """The grid column of each longitude (degrees); column 0 starts at 180 W, and 180 E is 180 W."""
    return np.floor((np.asarray(lon, dtype=np.float64) + 180.0) / 360.0 * W).astype(np.int64) % W


def row_area_km2():
    """The area of one cell of each row, on a sphere of the Earth's mean radius."""
    edges = np.radians(90.0 - np.arange(H + 1) * 180.0 / H)
    return R_EARTH**2 * (2 * math.pi / W) * (np.sin(edges[:-1]) - np.sin(edges[1:])) / 1e6


def land(path):
    import pyogrio
    import rasterio.features
    import shapely
    from rasterio.transform import from_bounds

    n = int(pyogrio.read_info(path)["features"])
    say("land: %d polygons onto %d x %d Mercator samples" % (n, MERC, MERC))
    transform = from_bounds(-HALF, -HALF, HALF, HALF, MERC, MERC)
    mask = np.zeros((MERC, MERC), dtype=bool)
    for s0 in range(0, n, BATCH):
        geoms = shapely.from_wkb(pyogrio.raw.read(path, columns=[], skip_features=s0, max_features=BATCH)[2])
        mask |= rasterio.features.rasterize(((g, 1) for g in geoms), out_shape=(MERC, MERC), transform=transform,
                                            fill=0, all_touched=False, dtype="uint8").astype(bool)
        del geoms
        say("land: polygons %d of %d" % (min(s0 + BATCH, n), n))
        report(step=min(s0 + BATCH, n), total=n, phase="land")
    # columns: the two grids cut longitude alike, MERC / W samples per cell
    per = MERC // W
    across = mask.reshape(MERC, W, per).sum(axis=2, dtype=np.uint32)
    del mask
    # rows: the latitude of each sample row's centre decides its cell row
    y = (np.arange(MERC) + 0.5) / MERC
    lat = np.degrees(np.arctan(np.sinh(math.pi * (1.0 - 2.0 * y))))
    r = rows_of(lat)
    starts = np.flatnonzero(np.r_[True, r[1:] != r[:-1]])
    sums = np.add.reduceat(across, starts, axis=0)
    counts = np.diff(np.r_[starts, MERC]) * per
    share = np.empty((H, W), dtype=np.float64)
    share[: r[0]] = 0.0  # north of the polygons: the Arctic Ocean
    share[r[-1] + 1:] = 1.0  # south of them: Antarctica
    share[r[starts]] = sums / counts[:, None]
    assert (np.diff(r[starts]) == 1).all(), "every row between the polygons' ends got samples"
    return share


def people(path):
    import rasterio
    from rasterio.windows import Window

    sums = np.zeros(H * W, dtype=np.float64)
    with rasterio.open(path) as src:
        tr, nodata = src.transform, src.nodata
        cols = cols_of(tr.c + (np.arange(src.width) + 0.5) * tr.a)
        say("people: %d x %d cells of %.5f degrees" % (src.width, src.height, tr.a))
        step = 512
        for r0 in range(0, src.height, step):
            h = min(step, src.height - r0)
            a = src.read(1, window=Window(0, r0, src.width, h)).astype(np.float64)
            if nodata is not None:
                a[a == nodata] = 0.0
            a[~np.isfinite(a) | (a < 0)] = 0.0
            rows = rows_of(tr.f + (r0 + np.arange(h) + 0.5) * tr.e)
            sums += np.bincount((rows[:, None] * W + cols[None, :]).ravel(), weights=a.ravel(), minlength=H * W)
            report(step=r0 + h, total=src.height, phase="people")
    total = float(sums.sum())
    assert abs(total - PEOPLE_TOTAL) < 1e-6 * PEOPLE_TOTAL, "people: %.0f, not the source's %d" % (total, PEOPLE_TOTAL)
    return sums.reshape(H, W) / row_area_km2()[:, None], total


def height(path):
    import netCDF4

    ds = netCDF4.Dataset(path)
    z = ds.variables["z"]
    lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
    cols = cols_of(np.asarray(ds.variables["lon"][:], dtype=np.float64))
    sums = np.zeros(H * W)
    counts = np.zeros(H * W)
    say("height: %d x %d cells" % (len(cols), len(lat)))
    for i0 in range(0, len(lat), 256):
        e = np.ma.filled(np.ma.asarray(z[i0:i0 + 256, :]).astype(np.float64), np.nan)
        rows = rows_of(lat[i0:i0 + 256])
        idx = (rows[:, None] * W + cols[None, :]).ravel()
        ok = np.isfinite(e.ravel())
        sums += np.bincount(idx[ok], weights=e.ravel()[ok], minlength=H * W)
        counts += np.bincount(idx[ok], minlength=H * W)
        report(step=i0 + 256, total=len(lat), phase="height")
    assert counts.min() > 0, "every cell got heights"
    return (sums / counts).reshape(H, W)


def main(work):
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {k: work / v for k, v in SOURCES.items()}
    for k, p in paths.items():
        assert p.is_file(), "%s: %s is missing" % (k, p)
    summary = {"grid": {"rows": H, "cols": W, "extent": "lon -180..180, lat 90..-90", "level": 3},
               "sources": {k: {"file": ARCHIVES[k], "sha256": sha256(work / ARCHIVES[k])} for k in ARCHIVES}}
    area = row_area_km2()[:, None] * np.ones((1, W))

    h = height(paths["height"])
    summary["height"] = {"min_m": float(h.min()), "max_m": float(h.max()), "mean_m": float((h * area).sum() / area.sum())}
    assert -11000 < h.min() < -8000 and 5000 < h.max() < 8848, summary["height"]
    np.save(OUT / "height.npy", h.astype(np.float16))
    say("height:", summary["height"])

    p, total = people(paths["people"])
    summary["people"] = {"total": round(total), "max_per_km2": float(p.max()), "cells_with_people": int((p > 0).sum())}
    np.save(OUT / "people.npy", p.astype(np.float32))  # the densest cells hold more than a float16 can (65504)
    say("people:", summary["people"])

    s = land(paths["land"])
    summary["land"] = {"km2": float((s * area).sum()), "share_of_sphere": float((s * area).sum() / area.sum())}
    assert 1.46e8 < summary["land"]["km2"] < 1.51e8, summary["land"]
    np.save(OUT / "land.npy", s.astype(np.float16))
    say("land:", summary["land"])

    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    say("written to", OUT)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / "weichseltree/chronon/data/maps/maplayers")
