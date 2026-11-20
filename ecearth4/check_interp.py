"""Check reduced_to_regular on the real source grid of a run with analytic fields.

Usage: python check_interp.py <exp> [<samples root>]   (root defaults to $ECE4_SAMPLES_ROOT)

Applies the interpolation weights of examples/custom_steps.py to smooth fields of zonal
wavenumber 1, 3 and 6 on the cells of <root>/samples/<exp>/output/*/oifs/<exp>_atm_cmip6_1m_*.nc
and compares with the exact field on the regular grid. Linear interpolation is second order:
the error must grow about as the square of the wavenumber, and a constant must be exact.
Works for any reduced Gaussian grid (TL linear, TCO octahedral). Exit 1 on failure.
"""
import glob
import importlib.util
import os
import sys

import numpy as np
import xarray as xr

MAX_ERR = {1: 2e-3, 3: 2e-2, 6: 8e-2}  # about 4 times the errors seen on TL63 and TCO79

exp = sys.argv[1]
root = sys.argv[2] if len(sys.argv) > 2 else os.environ["ECE4_SAMPLES_ROOT"]
spec = importlib.util.spec_from_file_location(
    "custom_steps", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples", "custom_steps.py"))
cs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cs)

src = xr.open_dataset(glob.glob(f"{root}/samples/{exp}/output/*/oifs/{exp}_atm_cmip6_1m_*.nc")[0])
lat, lon = src["lat"].values, src["lon"].values
rlat, rlon, i0, i1, w = cs._reduced_gaussian_index(lat, lon)


def field(la, lo, k):
    la, lo = np.radians(la), np.radians(lo)
    return np.cos(la) * np.cos(k * lo) + 0.5 * np.sin(la) ** 2


def interp(x):
    return np.where(w == 0, x[i0], x[i0] + w * (x[i1] - x[i0]))


failed = np.abs(interp(np.ones(lat.size)) - 1).max() > 1e-12
print(f"{exp} grid {rlat.size}x{rlon.size}: constant {'exact' if not failed else 'NOT exact'}")
for k, tol in MAX_ERR.items():
    err = np.abs(interp(field(lat, lon, k)) - field(rlat[:, None], rlon[None, :], k).ravel())
    failed |= err.max() > tol
    print(f"  k={k}: max error {err.max():.2e}, rms {np.sqrt((err**2).mean()):.2e} (limit {tol:.0e})")
sys.exit(int(failed))
