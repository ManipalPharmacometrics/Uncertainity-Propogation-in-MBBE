# Uncertainty Propagation in MBBE — Summary

## Purpose

This is a Quarto document (Julia + Pumas) that investigates **uncertainty propagation in Model-Based Bioequivalence (MBBE)**. It asks: *how does parameter uncertainty from PK model fitting affect bioequivalence (BE) decisions?*

## Workflow

### 1. PK Model Definition
A one-compartment oral PK model with depot absorption (CL, V, Ka), log-normal inter-individual variability (IIV), and a treatment covariate (`TRT`) that scales absorption by a bioavailability parameter (`tvbio`).

### 2. Simulate a BE Trial
`simulate_trial()` generates a parallel-design dataset (Reference vs Test, 100 mg dose, 30 subjects/arm, observations 0–72h).

### 3. Fit the Model
The simulated data is fit back using `Pumas.FOCE()` to recover parameter estimates.

### 4. Three Uncertainty Methods
Parameter uncertainty around the fitted estimates is quantified via:
- **VCOV** — Variance-covariance matrix from the fit (extracts SE of `tvbio`)
- **Bootstrap** — 200 bootstrap resamples
- **SIR** (Sampling Importance Resampling) — 200 samples / 200 resamples

### 5. NCA-based BE Decision
Runs NCA on the simulated data, computes AUC/Cmax, then calls `pumas_be()` to get a standard BE result (90% CI within 80–125%).

### 6. Virtual BE Simulation Loop
`virtual_be()` is the core experiment: it repeatedly draws parameter sets from the uncertainty distributions, simulates new trials, computes AUC, and runs TOST to get a BE pass/fail — estimating **power** (GMR=1.0) and **type I error** (GMR=0.75).

### 7. Summary Comparison
Compares power across the three uncertainty methods plus a "fixed parameter" (no uncertainty) baseline.

## Current State

The file is a **work in progress**. The early sections (model, simulation, fitting, uncertainty estimation, NCA-based BE) are fairly complete. The later sections (`virtual_be` function and the summary) have inconsistencies:

- References undefined variables: `samples`, `boot_mat`, `sir_mat`, `θ̂`, `compute_auc`, `tost_decision`
- The function signature of `virtual_be` doesn't match how it's called
- Several empty code cells at the end
