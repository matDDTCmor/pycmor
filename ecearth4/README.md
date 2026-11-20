# EC-Earth4 CMORization cases (MN5)

`run.sh plan|run <case> [<exp>]` renders `<case>.yaml.in` and runs `pycmor process` on it.
Cases: `nemo_core`, `oifs_core`, `lpjg_carbon`, `lpjg_yearly`. Settings: `settings.example`
(copy it to `ecearth4/ecearth4.env`, which git ignores). Run `run.sh` on a compute or interactive node
(heavy work is killed on the login nodes).

## 0. Prepare the pixi env (once)

`env_mn5.sh` builds everything under `PYCMOR_WORK` (pixi, caches, DReq venv, `metadata.json`,
`cdo`); the pycmor env itself goes in `<checkout>/.pixi/envs/default`, from `pixi.lock`.

1. `cp ecearth4/settings.example ecearth4/ecearth4.env` (git-ignored) and
   set `PYCMOR_WORK` (a directory with room for the env and caches, not `$HOME`) and
   `ECE4_SAMPLES_ROOT` (the samples root, see section 1). The other entries are commented out;
   uncomment one only to override its default: `CMIP7_CMOR_COMMON` (local cmip7-cmor-common
   checkout, 0.3.0 or later), `RUN_DIR`, `EXP_OUTPUT`, `NEMO_DOMAIN`.
2. Check what it will do: `ecearth4/env_mn5.sh plan`
3. Install on a login node with internet access (it downloads pixi, conda-forge packages,
   the submodules and the DReq API):
   `ecearth4/env_mn5.sh install`
   It installs the locked env (`pixi install --frozen`), cmip7-cmor-common, pycdo, cdo, and exports
   the DReq `metadata.json` for `DREQ_VERSION` (default v1.2.2.4). It ends by running `check`.
4. Check any time, offline: `ecearth4/env_mn5.sh check`. Every module must say `ok`, and
   none `OUTSIDE ENV`. Rerun `install` if it reports an older cmip7-cmor-common or a missing
   pycdo or cdo.

Compute nodes have no internet: install first, then run.

## 1. Run abxi (eORCA1, N128 atmosphere)

Sample data layout: `<ECE4_SAMPLES_ROOT>/samples/abxi/output/<leg>/{nemo,oifs}/abxi_*.nc`.

1. Prepare the env (section 0).
2. `ECE4_SAMPLES_ROOT` is set in `ecearth4.env`. Leave `NEMO_DOMAIN` commented out: `run.sh` uses
   `samples/abxi/domain_cfg.nc` (eORCA1), so the file must be there.
3. Check the paths: `ecearth4/run.sh plan nemo_core abxi`
4. Run the NEMO test: `ecearth4/run.sh run nemo_core abxi` (or submit it with `hpc submit`).
5. Run the OIFS test: `ecearth4/run.sh run oifs_core abxi` (same options; no domain file needed).
   The input is `*/oifs/abxi_atm_cmip6_{1m,1d}_*.nc` (reduced Gaussian, regridded to a regular
   256 x 512 grid by `reduced_to_regular`).
6. Check the remapping: `pixi run python ecearth4/check_remap.py abxi <ECE4_SAMPLES_ROOT>`
   (as in section 0, with `--frozen --manifest-path pyproject.toml`). It compares the monthly
   `tas`, `psl`, `pr` output with the source (grid latitudes, NaNs, global means, row means and the
   pattern against the nearest source cell) and exits 1 on a failure. It works for any reduced
   Gaussian source grid (TL linear or TCO octahedral). Same command for abx2 and abx3.
   `reduced_to_regular` is linear interpolation along each latitude row (periodic, no meridional
   interpolation, output longitudes = the longest row), so it needs no change for TCO grids. Two more
   checks, same arguments (`<exp> <ECE4_SAMPLES_ROOT>`), both exit 1 on a failure:
   - `ecearth4/check_interp.py`: analytic fields on the real source grid; the error must grow as the
     square of the wavenumber (second order) and a constant must be exact;
   - `ecearth4/check_cdo.py`: against `cdo setgridtype,regular` on `tas` (cdo on `PATH`, as
     `<PYCMOR_WORK>/pixi-global/bin`): same grid, differences below 1e-3 of the field sd. CDO reads the
     XIOS files as unstructured, so the script gives it a `gaussian_reduced` grid description built
     from the row counts.
   All three scripts pass for abxi, abx2 and abx3.
7. Run the LPJ-GUESS tests: `ecearth4/run.sh run lpjg_carbon abxi` and
   `ecearth4/run.sh run lpjg_yearly abxi` (same options). The input is
   `*/lpjg/<var>_{monthly,yearly}.out[.gz]` (`.out.gz` is read as is). The land points are put on the
   reduced Gaussian grid of the OIFS files (`*/oifs/*_1m_*.nc`, so the OIFS output must exist) and
   remapped conservatively to the regular N128 grid with `cdo` (installed by `env_mn5.sh`).
   - `lpjg_carbon` (monthly): `cLand`, `netAtmosLandCO2Flux`, `nbp` (land means, missing over
     ocean), `fco2antt`, `fco2nat` (grid-cell means, 0 over ocean).
   - `lpjg_yearly`: `vegFrac`, `treeFrac` (fractions of the grid cell, 0 over ocean).
   Both ran for abxi (1990, rc 0); the values were not checked.

Output: `<ECE4_SAMPLES_ROOT>/pycmor-out/abxi`, log: `<ECE4_SAMPLES_ROOT>/logs/abxi/`.

## 2. Run abx2 (ORCA2, TL63 atmosphere)

Same steps, with these differences:

- **Domain file**: `samples/abx2/domain_cfg.nc` (ORCA2), found automatically. Keep
  `NEMO_DOMAIN` commented out in `ecearth4.env`: it would override it.
- **Experiment name**: the file names must start with the experiment name, in the `nemo/` and
  the `oifs/` folders (`abx2_oce_1m_T_*.nc`, `abx2_atm_cmip6_1m_*.nc`), because the patterns in
  `nemo_core.yaml.in` and `oifs_core.yaml.in` use `@EXP@`. A mismatch gives
  `OSError: no files to open` (the abx2 OIFS files must be renamed from `abxi_atm_*` too).
- **OIFS grid**: TL63 is regridded to a regular 64 x 128 grid (N32, not 256 x 512).
  `oifs_core.yaml.in` still sets the abxi values (`grid_label` g114, "N128", `nominal_resolution`
  "100 km"), so these labels are also provisional. `g210` in the EMD (regular Gaussian,
  128 x 64) is the matching grid cell, but is not used yet.
- **Grid labels (provisional)**: the EMD has no ORCA2 grid, so `nemo_core.yaml.in` still sets
  the eORCA1 values (`grid_label` g147, g148 for `tauuo`, `nominal_resolution` "100 km", `grid`
  text "eORCA1 ... 360 x 331"). The abx2 files carry these labels in their attributes and in the
  output path. Treat the output as test data until an ORCA2 entry exists.
- **Queue**: `gp_debug` allows one job per user; another debug job blocks it.

Results (1990, rc 0), for abxi and abx2:
- `nemo_core`: `areacello`, `tos`, `thetao`, `tauuo`, `siconc` (monthly).
- `oifs_core`: `tas` (daily `tavg`, `tmax`, `tmin`; monthly `tavg`, `tmaxavg`, `tminavg`), `pr`, `psl`
  (monthly). The reduced-to-regular remapping preserves the global means (checked on `tas`, `psl`, `pr`).

## 3. abx3 (TCO atmosphere)

Octahedral reduced Gaussian grid (TCO79): 28480 cells, 160 latitude rows of 20 to 336 points. The OIFS
test (`ecearth4/run.sh run oifs_core abx3`) regrids it to a regular 160 x 336 grid (the longitudes
follow the longest row, so not N80's 320).

- **Layout**: `samples/abx3/output/<leg>/{nemo,oifs}/abx3_*.nc` (`001`), with the OIFS files named
  `abx3_atm_cmip6_*.nc`, as for abx2.
- **Result** (1990, rc 0): the same 8 `oifs_core` files as abxi. `check_remap.py abx3` passes
  (latitudes at the Gaussian values, global means within 1.1e-4, pattern within 0.08 sd).
  `check_interp.py` and `check_cdo.py` pass too.
- **Labels (provisional)**: `oifs_core.yaml.in` has the abxi values (`g114`, "N128", "100 km"),
  which do not match this grid.
