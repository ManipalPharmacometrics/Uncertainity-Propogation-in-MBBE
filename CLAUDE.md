# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Simulation framework comparing four parameter uncertainty propagation methods (Fixed, VCOV, Bootstrap, SIR) in model-based bioequivalence (MBBE) studies. Written in Julia, using the Pumas/Bioequivalence ecosystem for PK modeling, NCA, and TOST analysis.

## Commands

### Sequential execution (from `mbbe-uncertainty-v3-split/`)
```bash
julia 01_model_fit.jl          # Fit PK model, extract uncertainty samples
julia 02_run_fixed.jl          # Simulate Fixed method
julia 02_run_vcov.jl           # Simulate VCOV method
julia 02_run_bootstrap.jl      # Simulate Bootstrap method
julia 02_run_sir.jl            # Simulate SIR method
quarto render 03_results_tables.qmd   # Generate results tables HTML
quarto render 04_visualizations.qmd   # Generate visualization plots HTML
```

### Parallel execution (20 independent jobs: 4 methods × 5 T/R ratios)
```bash
julia runs/job_01_fixed_tr070.jl   # Run a single job
julia runs/_runner.jl Fixed 0.90   # Or use the generic runner with CLI args
julia runs/check_status.jl         # Monitor completion (.done markers)
julia runs/combine_results.jl      # Aggregate job outputs into batches_*.jls
```

The runner also accepts env vars: `JOB_METHOD=Fixed JOB_TR_RATIO=0.90 julia runs/_runner.jl`

## Architecture

### Data flow (4 phases)
1. **Model fitting** (`01_model_fit.jl`) — Simulates a reference trial with true parameters, fits via FOCE, extracts parameter samples using VCOV/Bootstrap/SIR → saves `model_fit_results.jls`
2. **Virtual BE simulation** (`02_run_*.jl` or `runs/` jobs) — For each (method, T/R ratio, sample size), runs 200 batches × 200 trials: draw params → simulate PK → NCA → TOST → pass/fail → saves `batches_*.jls`
3. **Results aggregation** (`runs/combine_results.jl`) — Merges 20 per-job outputs into 4 method-level `batches_*.jls` files (only needed for parallel execution)
4. **Analysis** (`03_results_tables.qmd`, `04_visualizations.qmd`) — Loads all batches, computes summary statistics and relative metrics vs Fixed, renders HTML

### Key file: `common.jl`
Central shared library included by all scripts. Contains:
- PK model definition (one-compartment oral, first-order elimination)
- True parameter values and simulation grid constants (sample sizes, T/R ratios, batch/trial counts)
- Core functions: `create_be_population`, `simulate_be_trial`, `sample_from_vcov`, `tost_test_detailed`, `run_be_on_sim_detailed`, `virtual_be_sim_batched`, `summarize_batches`, `load_all_method_results`
- Plotting constants (colors, shapes, offsets)

### `runs/_runner.jl`
Universal job engine with three parameter resolution modes (priority order): constants from job file → CLI args → env vars. When invoked with no args, runs all 20 combinations. Uses a reduced set of sample sizes (5 instead of 7) in parallel mode.

## Important Patterns

- **VCOV sampling is done on log-scale** then exponentiated back, ensuring all sampled parameters remain positive
- **Deterministic seeding**: batch seeds use `hash((method, tr_ratio, n_arm, batch_id)) % typemax(Int32)` for reproducibility and independence
- **NCA column detection is dynamic** — tries multiple candidate names (`aucinf_obs`, `aucinf`, `auclast`, etc.) to stay software-agnostic
- **BE pass requires both AUC and Cmax** to fall within 80–125% limits (AND logic)
- **Error handling**: NCA failures are caught and excluded from pass rate (`n_pass / n_valid`), not from total count
- All intermediate results use Julia's `Serialization` format (`.jls` files)

## Dependencies

Julia 1.10+ with: Pumas, Bioequivalence, CairoMakie, DataFrames, Distributions, LinearAlgebra, Statistics, Random, Serialization. Quarto required for rendering `.qmd` reports.
