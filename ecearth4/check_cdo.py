"""Compare reduced_to_regular with `cdo setgridtype,regular` on the monthly tas of a run.

Usage: python check_cdo.py <exp> [<samples root>]   (root defaults to $ECE4_SAMPLES_ROOT)
Needs cdo on PATH (run.sh and env_mn5.sh use <PYCMOR_WORK>/pixi-global/bin).

CDO reads the XIOS files as unstructured grids, so the reduced Gaussian structure is given with a
grid description built from the row counts of the source. The result must be on the same grid as
the cmorized output <root>/pycmor-out/<exp> and differ by less than REL_TOL of the field's
standard deviation (observed: below 1e-4). Exit 1 on failure.
"""
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import xarray as xr

REL_TOL = 1e-3

exp = sys.argv[1]
root = sys.argv[2] if len(sys.argv) > 2 else os.environ["ECE4_SAMPLES_ROOT"]
srcfile = glob.glob(f"{root}/samples/{exp}/output/*/oifs/{exp}_atm_cmip6_1m_*.nc")[0]
lat = xr.open_dataset(srcfile)["lat"].values
cnt = np.diff(np.concatenate([[0], np.flatnonzero(np.abs(np.diff(lat)) > 1e-6) + 1, [lat.size]]))

with tempfile.TemporaryDirectory() as tmp:
    grid = os.path.join(tmp, "grid")
    with open(grid, "w") as f:
        f.write(f"gridtype = gaussian_reduced\ngridsize = {lat.size}\nxsize = 0\nysize = {cnt.size}\n"
                f"np = {cnt.size // 2}\nreducedPoints = {' '.join(map(str, cnt))}\n")
    out = os.path.join(tmp, "reg.nc")
    subprocess.run(f"cdo -s setgridtype,regular -setgrid,{grid} -selname,tas {srcfile} {out}",
                   shell=True, check=True, stderr=subprocess.DEVNULL)
    cdo = xr.open_dataset(out)["tas"]
    cv = cdo.values.reshape(cdo.shape[-3:])
    if cdo[cdo.dims[-2]].values[0] < cdo[cdo.dims[-2]].values[-1]:  # CDO writes south to north
        cv = cv[:, ::-1]
    files = sorted(glob.glob(f"{root}/pycmor-out/{exp}/**/mon/tas/tavg-h2m*/**/tas_*mon*.nc", recursive=True))
    ov = xr.open_dataset(files[-1])["tas"]
    ovals = ov.values if ov["lat"].values[0] > ov["lat"].values[-1] else ov.values[:, ::-1]
    if cv.shape != ovals.shape:
        print(f"{exp}: different grids, cdo {cv.shape}, pycmor {ovals.shape}")
        sys.exit(1)
    dd = np.abs(cv - ovals)
    rel = dd.max() / ovals.std()
    print(f"{exp} {ovals.shape}: max |cdo - pycmor| {dd.max():.2e} K, rms {np.sqrt((dd**2).mean()):.2e} K, "
          f"{rel:.1e} of the field sd: {'ok' if rel < REL_TOL else 'FAIL'}")
    sys.exit(int(rel >= REL_TOL))
