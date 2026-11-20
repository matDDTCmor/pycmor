"""Check the OIFS reduced Gaussian -> regular remapping of a cmorized run.

Usage: python check_remap.py <exp> [<samples root>]   (root defaults to $ECE4_SAMPLES_ROOT)

For tas, psl and pr (monthly) compares the output with the source file
<root>/samples/<exp>/output/*/oifs/<exp>_atm_cmip6_1m_*.nc:
  - output latitudes against the Gaussian latitudes, no NaN;
  - global mean per month, source (row Gaussian weight / cells in the row) against
    output (cos(lat) weights), relative difference below TOL;
  - row means (each output latitude against the source cells of that latitude), relative to
    the standard deviation of the field, below ROW_TOL;
  - pattern: output against the nearest source cell found on the sphere with a KD-tree (no row
    logic, so it also catches shifted longitudes or a wrong seam), RMS relative to the standard
    deviation of the field, below NN_TOL.
Works for any reduced Gaussian source grid with Gaussian rows (TL linear, TCO octahedral).
Exit 1 when a check fails. Thresholds are 2 to 5 times the worst value seen (abxi, abx2, abx3).
Two months (January, July) and no check of the polar rows beyond the row means.
"""
import glob
import os
import sys

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

TOL = 1e-3  # relative difference of the global means (observed: 1e-4 or less)
ROW_TOL = 0.01  # row means, in standard deviations of the field
NN_TOL = 0.15  # pattern RMS against the nearest source cell, in standard deviations

exp = sys.argv[1]
root = sys.argv[2] if len(sys.argv) > 2 else os.environ["ECE4_SAMPLES_ROOT"]
src = xr.open_dataset(glob.glob(f"{root}/samples/{exp}/output/*/oifs/{exp}_atm_cmip6_1m_*.nc")[0])

# source weights: each row (same latitude) gets its Gaussian weight, shared by its cells
ul, inv, cnt = np.unique(np.round(src["lat"].values, 4), return_inverse=True, return_counts=True)
sw = np.polynomial.legendre.leggauss(ul.size)[1][np.argsort(np.argsort(ul))][inv] / cnt[inv]


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


tree = cKDTree(xyz(src["lat"].values, src["lon"].values))

failed = False
for var in ("tas", "psl", "pr"):
    files = sorted(glob.glob(f"{root}/pycmor-out/{exp}/**/mon/{var}/tavg*/**/{var}_*.nc", recursive=True))
    if not files:
        print(f"{var}: no output")
        failed = True
        continue
    v = xr.open_dataset(files[-1])[var]
    lat, lon = v["lat"].values, v["lon"].values
    gauss = np.degrees(np.arcsin(np.polynomial.legendre.leggauss(lat.size)[0]))
    dlat = np.abs(np.sort(lat) - gauss).max()
    w = np.cos(np.radians(lat))[:, None] * np.ones(lon.size)
    out = np.array([(v.isel(time=t).values * w).sum() / w.sum() for t in range(v.sizes["time"])])
    ref = np.array([(src[var].isel(time_counter=t).values * sw).sum() / sw.sum() for t in range(src.sizes["time_counter"])])
    rel = np.abs(out - ref) / np.abs(ref)
    # row means and nearest-cell pattern, first month and July (index 6)
    row_of = np.abs(np.unique(np.round(src["lat"].values, 4))[None, :] - lat[:, None]).argmin(axis=1)
    _, nn = tree.query(xyz(lat[:, None] * np.ones(lon.size), np.ones(lat.size)[:, None] * lon).reshape(-1, 3))
    row_err = nn_err = 0.0
    for t in (0, 6):
        sv = src[var].isel(time_counter=t).values
        std = sv.std()
        src_rows = np.bincount(inv, weights=sv) / cnt
        row_err = max(row_err, np.abs(v.isel(time=t).values.mean(axis=1) - src_rows[row_of]).max() / std)
        nn_err = max(nn_err, np.sqrt(((v.isel(time=t).values.reshape(-1) - sv[nn]) ** 2).mean()) / std)
    ok = dlat < 1e-3 and not np.isnan(v.values).any() and rel.max() < TOL and row_err < ROW_TOL and nn_err < NN_TOL
    failed |= not ok
    print(f"{var} {v.shape} lat vs Gaussian {dlat:.1e} deg, lon step {np.diff(lon).mean():.4f}, "
          f"global mean {rel.max():.1e}, row means {row_err:.1e} sd, pattern {nn_err:.1e} sd: {'ok' if ok else 'FAIL'}")
sys.exit(int(failed))
