"""One block per netCDF file: data variables with shape, units and statistics,
and the global attributes that identify the dataset. Used by run.sh to keep
results in the log instead of tracking output files."""
import sys

import netCDF4 as nc
import numpy as np

GLOBALS = ("variable_id", "branded_variable", "frequency", "grid_label",
           "source_id", "experiment_id", "variant_label", "activity_id")


def stats(var):
    """min, max, mean, nan count over a variable, one slice of the first axis at a time
    (whole variables do not fit in memory on a login node); masked values count as nan."""
    lo, hi, tot, n, nan = np.inf, -np.inf, 0.0, 0, 0
    for i in range(var.shape[0]):
        x = np.ma.filled(np.ma.asarray(var[i]).astype("f8"), np.nan)
        ok = ~np.isnan(x)
        nan += x.size - int(ok.sum())
        if ok.any():
            lo, hi = min(lo, x[ok].min()), max(hi, x[ok].max())
            tot += float(x[ok].sum())
            n += int(ok.sum())
    return (lo, hi, tot / n, nan) if n else (np.nan, np.nan, np.nan, nan)


for path in sys.argv[1:]:
    print(f"file {path}")
    with nc.Dataset(path) as ds:
        print("  " + " | ".join(f"{k}={ds.getncattr(k)}" for k in GLOBALS if k in ds.ncattrs()))
        for name, var in ds.variables.items():
            if name in ds.dimensions or "bnds" in name or var.ndim < 2:
                continue
            lo, hi, mean, nan = stats(var)
            print(f"  {name}{tuple(var.dimensions)} {var.shape} {getattr(var, 'units', '')} "
                  f"min={lo:.6g} max={hi:.6g} mean={mean:.6g} nan={nan}")
        if "time" in ds.variables and ds["time"].size:
            t = ds["time"]
            print(f"  time {t[0]:g}..{t[-1]:g} {t.units} {getattr(t, 'calendar', '')} n={t.size}")
