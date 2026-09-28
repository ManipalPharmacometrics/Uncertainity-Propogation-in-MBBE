# Corrected re-run (2026-09-28)

The grids stored in `../batches_*.jls` and `../runs/` were generated with a bug in
`sample_from_vcov` (`../common.jl`): the covariance from `vcov(fit)` is on the natural
scale of the parameters, but it was used as a log-scale covariance. For V that gave a
log-SD of 2.1 instead of about 0.04, so drawn volumes spanned about 1–2,100 L. Under the
additive error (σ = 0.2 mg/L), large-V trials are mostly noise, the NCA ratios collapse
towards 1, and TOST passes at T/R 0.70 and 1.30. That produced the 10–23% false
acceptance in the stored VCOV, Bootstrap and SIR results, rising with n.
**Those stored results should not be used.** The Fixed arm did not call the function.

`../common.jl` now converts the covariance to the log scale by the delta method,
`Cov(log θ) ≈ D⁻¹ V D⁻¹` with `D = diag(θ̂)`, before sampling.

## What this folder contains

| File | Purpose |
|---|---|
| `00_scale_error_diagnostic.jl` | Reproduces the failure: original stored draws vs corrected draws, 600 trials per cell, per-trial drawn V |
| `rerun_common.jl` | Loads `../common.jl` and sets the AUC endpoint to AUC0–t (`auclast`) |
| `01_uncertainty.jl` | Refits the reference study (seed 142, reproduces θ̂) and builds the parameter sets: VCOV 1,000 draws; Bootstrap 500 refit vectors; SIR 2,000 proposals / 1,000 resamples (empirical vectors, not an MVN) |
| `02_grid.jl` | 4 methods × n = 12/24/30/40/50/60/80 × T/R 0.70/0.90/1.00/1.10/1.30, 20 batches × 100 trials per cell; resumable per cell |
| `03_summarise.jl` | Operating characteristics, differences from Fixed, n for 80% power |
| `04_cv.jl`, `04_fig3_scale_error.jl`, `05_gmr_check.jl` | Endpoint CVs, scale-error figure, realised GMRs at T/R 0.70 / 1.30 |
| `run_all.sh` | Runs 01 then 02 (set `JULIA` to your Pumas invocation) |
| `output/` | Summary CSVs, per-cell CSVs, `tables.md`, run log |

## Results (140 cells, 280,000 trials, none excluded)

- False acceptance ≤ 0.2% at T/R 0.70 and 1.0–2.2% at T/R 1.30 for all four methods, with no trend in n.
- Power at T/R 1.00 differs from Fixed by at most 3.4 points (n ≤ 40). n for 80% power:
  Fixed 21.7, VCOV 22.1, Bootstrap 21.5, SIR 22.2 per arm.
- Changing the AUC endpoint from AUC0–∞ to AUC0–t alone moved n for 80% power from about 33 to 22.
- The additive error compresses observed ratios towards 1 for both endpoints (AUC0–t more than
  Cmax). At a true T/R of 1.30 the median AUC0–t ratio is 1.20, so false acceptance stays near 1%
  only because Cmax is co-primary (`output/gmr_check.txt`).

`../01_model_fit.jl` and the `02_run_*.jl` scripts still use AUC0–∞ and the original
replicate counts; re-running them with the fixed `common.jl` does not reproduce this folder.
