# Uncertainty Propagation in Model-Based Bioequivalence (MBBE)

A simulation framework for comparing how different parameter uncertainty methods affect the operating characteristics (power, Type I error) of model-based bioequivalence (MBBE) studies using a parallel-group design.


> **2026-09-28 — sampler bug fixed.** `sample_from_vcov` in `mbbe-uncertainty-v3-split/common.jl` used the
> natural-scale covariance from Pumas as a log-scale covariance. The VCOV, Bootstrap and SIR results stored in
> `batches_*.jls` and `runs/` (10–23% false acceptance, 13–16 pp power loss) are artefacts of that bug and
> should not be used. The fix and the corrected re-run are in
> [`mbbe-uncertainty-v3-split/rerun-corrected/`](mbbe-uncertainty-v3-split/rerun-corrected/README.md).

## Background

In MBBE, a pharmacokinetic (PK) model is fitted to clinical data, and the estimated parameters are used to simulate virtual bioequivalence trials. A key question is: **how should uncertainty in the estimated parameters be propagated into these simulations?**

This project compares four approaches:

| Method | Description |
|--------|-------------|
| **Fixed (plug-in)** | Uses point estimates only; no parameter uncertainty |
| **Asymptotic Variance-Covariance (VCOV)** | Samples parameters from the asymptotic multivariate normal distribution derived from the Fisher information matrix |
| **Bootstrap** | Re-fits the model on bootstrapped datasets to obtain an empirical parameter distribution |
| **SIR (Sampling Importance Resampling)** | Uses importance-weighted resampling from a proposal distribution to approximate the posterior |

For each method, virtual BE trials are simulated across a grid of:
- **Sample sizes**: 12, 24, 30, 40, 50, 60, 80 subjects per arm
- **Test/Reference (T/R) ratios**: 0.70, 0.90, 1.00, 1.10, 1.30

Bioequivalence is assessed via the Two One-Sided Tests (TOST) procedure on AUC and Cmax derived from non-compartmental analysis (NCA), using the standard 80-125% acceptance limits.

## PK Model

A one-compartment oral absorption model with first-order elimination:
- **Parameters**: CL (clearance), V (volume), Ka (absorption rate), bioavailability (controls T/R ratio)
- **Random effects**: Log-normal IIV on CL, V, Ka (~30% CV)
- **Residual error**: Additive normal

## Project Structure

```
mbbe-uncertainty-v3-split/
├── common.jl                  # Shared code: model, helpers, simulation engine, constants
├── 01_model_fit.jl            # Step 1: Fit PK model, extract uncertainty (VCOV/Bootstrap/SIR)
├── 02_run_fixed.jl            # Step 2a: Run virtual BE sims — Fixed method
├── 02_run_vcov.jl             # Step 2b: Run virtual BE sims — VCOV method
├── 02_run_bootstrap.jl        # Step 2c: Run virtual BE sims — Bootstrap method
├── 02_run_sir.jl              # Step 2d: Run virtual BE sims — SIR method
├── 03_results_tables.qmd      # Step 3: Quarto doc — absolute & relative results tables
├── 04_visualizations.qmd      # Step 4: Quarto doc — power/Type I error plots
├── model_fit_results.jls      # Serialized model fit + parameter samples
├── batches_*.jls              # Serialized batch simulation results per method
└── runs/                      # Parallelized job runner for HPC/multi-process execution
    ├── _runner.jl             # Shared job engine (one method x one T/R ratio per job)
    ├── check_status.jl        # Monitor job completion status
    ├── combine_results.jl     # Merge per-job outputs into final batches_*.jls files
    ├── job_01_fixed_tr070.jl  # Individual job files (20 total: 4 methods x 5 T/R ratios)
    ├── ...
    ├── job_20_sir_tr130.jl
    └── output/                # Per-job results (.jls) and completion markers (.done)
```

## How to Run

### Prerequisites

- [Julia](https://julialang.org/) (tested with Julia 1.10+)
- Julia packages: `Pumas`, `Bioequivalence`, `CairoMakie`, `DataFrames`, `Distributions`, `Serialization`
- [Quarto](https://quarto.org/) (for rendering results documents)

### Sequential Execution

```bash
cd mbbe-uncertainty-v3-split

# Step 1: Fit model and extract parameter uncertainty
julia 01_model_fit.jl

# Step 2: Run simulations for each method (each takes significant time)
julia 02_run_fixed.jl
julia 02_run_vcov.jl
julia 02_run_bootstrap.jl
julia 02_run_sir.jl

# Step 3-4: Generate results
quarto render 03_results_tables.qmd
quarto render 04_visualizations.qmd
```

### Parallel Execution (recommended)

The `runs/` directory splits the simulation grid into 20 independent jobs (one per method-ratio combination), suitable for HPC clusters or parallel local execution:

```bash
# Run individual jobs (can be parallelized across processes/nodes)
julia runs/job_01_fixed_tr070.jl
julia runs/job_02_fixed_tr090.jl
# ... or use the generic runner:
julia runs/_runner.jl Fixed 0.90

# Monitor progress
julia runs/check_status.jl

# Combine into final result files
julia runs/combine_results.jl
```

Each job runs 5 sample sizes x 200 batches x 200 trials = **200,000 simulations**.

## Simulation Design

The simulation uses a **macro-replication** strategy:

1. **200 batches** (macro-replications) per grid cell
2. **200 trials** per batch — each trial draws parameters from the method's uncertainty distribution, simulates a full BE trial, runs NCA, and applies TOST
3. Pass rates are computed per batch, then aggregated with 95% CIs across batches

This design provides both point estimates and uncertainty intervals for operating characteristics.

## Outputs

- **Results tables** (`03_results_tables.html`): Absolute power/Type I error at each grid point, plus differences and ratios relative to the Fixed method
- **Visualizations** (`04_visualizations.html`): Pass rate curves with CI bands, method comparison overlays, and relative-to-Fixed difference/ratio plots

## Key Dependencies

- [Bioequivalence.jl](https://docs.pumas.ai/) — BE-specific analysis tools
- [CairoMakie](https://docs.makie.org/stable/) — Publication-quality plotting
- [Quarto](https://quarto.org/) — Reproducible scientific documents

