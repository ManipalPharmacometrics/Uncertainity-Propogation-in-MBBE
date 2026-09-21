# ──────────────────────────────────────────────────────────────────────────────
# common.jl — Shared libraries, model, helpers, and configuration
#
# This file is `include()`d by all other scripts in this folder.
# It provides:
#   - Package imports
#   - PK model definition and true parameter values
#   - Simulation helpers (population creation, trial simulation)
#   - VCOV sampling, TOST, NCA helpers
#   - %CV estimation
#   - Batch simulation engine and summary aggregation
#   - Formatting helpers for display
#   - Simulation and plotting constants
# ──────────────────────────────────────────────────────────────────────────────

using Pumas
using Random
using Distributions
using LinearAlgebra
using Statistics
using DataFrames
using Bioequivalence
using CairoMakie
using Makie
using Logging
using Serialization

# ── Global NCA column names (set by 01_model_fit.jl, loaded by other scripts) ──
AUC_COL  = :aucinf_obs   # default placeholder — overridden after NCA detection
CMAX_COL = :cmax          # default placeholder — overridden after NCA detection

# ── PK Model ─────────────────────────────────────────────────────────────────
# `tvbio` controls relative bioavailability of Test vs Reference via @dosecontrol,
# so AUC_T/AUC_R = tvbio exactly.

pkmodel = @model begin
    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvv   ∈ RealDomain(lower=0)
        tvka  ∈ RealDomain(lower=0)
        tvbio ∈ RealDomain(lower=0)
        Ω     ∈ PDiagDomain(3)
        σ     ∈ RealDomain(lower=0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates TRT

    @pre begin
        CL  = tvcl * exp(η[1])
        V   = tvv  * exp(η[2])
        Ka  = tvka * exp(η[3])
    end

    @dosecontrol begin
        bioav = (Depot = TRT == "T" ? tvbio : 1.0,)
    end

    @dynamics begin
        Depot'   = -Ka * Depot
        Central' =  Ka * Depot - (CL / V) * Central
    end

    @derived begin
        cp = @. Central / V
        dv ~ @. Normal(cp, σ)
    end
end

# ── True Parameter Values ────────────────────────────────────────────────────
θ_true = (
    tvcl  = 5.0,
    tvv   = 50.0,
    tvka  = 1.2,
    tvbio = 1.0,
    Ω     = Diagonal([0.09, 0.09, 0.09]),   # ~30% CV IIV
    σ     = 0.2
)

# ── Simulation Configuration ─────────────────────────────────────────────────
# Adjust N_BATCHES and N_SIMS_PER_BATCH for faster testing
# (e.g., 10 batches × 50 trials).
const SAMPLE_SIZES     = [12, 24, 30, 40, 50, 60, 80]
const TR_RATIOS        = [0.7, 0.9, 1.0, 1.1, 1.3]
const N_BATCHES        = 200   # macro-replications
const N_SIMS_PER_BATCH = 200   # trials per macro-replication

# ── Plotting Constants ───────────────────────────────────────────────────────
const METHODS_LIST   = ["Fixed", "VCOV", "Bootstrap", "SIR"]
const METHOD_COLORS  = [:steelblue, :crimson, :seagreen, :darkorange]
const METHOD_SHAPES  = [:circle, :diamond, :utriangle, :rect]
const METHOD_OFFSETS = [-1.5, -0.5, 0.5, 1.5]
const SIZE_COLORS    = [:indigo, :dodgerblue, :steelblue, :seagreen, :orange, :crimson, :darkred]
const REL_METHODS    = ["VCOV", "Bootstrap", "SIR"]
const REL_COLORS     = [:crimson, :seagreen, :darkorange]
const REL_SHAPES     = [:diamond, :utriangle, :rect]
const REL_OFFSETS    = [-1.0, 0.0, 1.0]

# ── Simulation Helpers ───────────────────────────────────────────────────────

function create_be_population(n_per_arm)
    subjects = Subject[]
    for trt in ["R", "T"], i in 1:n_per_arm
        push!(subjects, Subject(
            id     = length(subjects) + 1,
            events = DosageRegimen(100, time = 0),
            covariates = (TRT = trt,)
        ))
    end
    return Population(subjects)
end

function simulate_be_trial(params; n_per_arm = 30)
    pop = create_be_population(n_per_arm)
    return simobs(pkmodel, pop, params; obstimes = 0:1:72)
end

# ── VCOV Sampling Helper ────────────────────────────────────────────────────

"""
Sample `n` parameter NamedTuples from a MvNormal on the estimation (log) scale.
"""
function sample_from_vcov(mean_params, V, n)
    μ = [
        log(mean_params.tvcl), log(mean_params.tvv),
        log(mean_params.tvka), log(mean_params.tvbio),
        log(mean_params.Ω[1,1]), log(mean_params.Ω[2,2]), log(mean_params.Ω[3,3]),
        log(mean_params.σ)
    ]
    dist = MvNormal(μ, Symmetric(V))
    return [begin
        s = rand(dist)
        (tvcl  = exp(s[1]), tvv  = exp(s[2]),
         tvka  = exp(s[3]), tvbio = exp(s[4]),
         Ω     = Diagonal([exp(s[5]), exp(s[6]), exp(s[7])]),
         σ     = exp(s[8]))
    end for _ in 1:n]
end

# ── TOST for Parallel-Design BE ──────────────────────────────────────────────

"""
Two One-Sided Tests on log-scale data for two independent groups.
Returns GMR, 90% CI bounds, and pass/fail status.
"""
function tost_test_detailed(log_R, log_T; alpha = 0.05, θ_L = 0.80, θ_U = 1.25)
    n_R, n_T = length(log_R), length(log_T)
    (n_R < 2 || n_T < 2) && return (gmr=NaN, ci_lo=NaN, ci_hi=NaN, pass=false)
    Δ   = mean(log_T) - mean(log_R)
    sp² = ((n_R - 1) * var(log_R) + (n_T - 1) * var(log_T)) / (n_R + n_T - 2)
    se  = sqrt(sp² * (1 / n_R + 1 / n_T))
    df  = n_R + n_T - 2
    tc  = quantile(TDist(df), 1 - alpha)
    gmr   = exp(Δ)
    ci_lo = exp(Δ - tc * se)
    ci_hi = exp(Δ + tc * se)
    pass  = ci_lo >= θ_L && ci_hi <= θ_U
    return (gmr=gmr, ci_lo=ci_lo, ci_hi=ci_hi, pass=pass)
end

"""Find first matching column from `candidates` in DataFrame `df`."""
function find_col(df, candidates)
    ns = names(df)
    for c in candidates
        if string(c) in ns
            return Symbol(c)
        end
    end
    error("None of $candidates found in columns: $ns")
end

"""Extract non-missing, positive values from a column as Float64 vector."""
function clean_pk_col(df, col)
    vals = df[!, col]
    return Float64[v for v in vals if !ismissing(v) && v > 0]
end

"""
Run NCA on a simulated trial, then apply TOST to AUC and Cmax.
Returns `nothing` if insufficient data.
"""
function run_be_on_sim_detailed(sim)
    simdf = DataFrame(sim)
    simdf.route .= "ev"
    simdf.dv .= max.(simdf.dv, 0.0)

    tbl = with_logger(NullLogger()) do
        ncapop = read_nca(simdf;
            id = :id, time = :time, amt = :amt,
            observations = :dv, route = :route, group = [:TRT]
        )
        nca_out = run_nca(ncapop)
        nca_out.reportdf
    end

    R = tbl[tbl.TRT .== "R", :]
    T = tbl[tbl.TRT .== "T", :]

    auc_R  = clean_pk_col(R, AUC_COL)
    auc_T  = clean_pk_col(T, AUC_COL)
    cmax_R = clean_pk_col(R, CMAX_COL)
    cmax_T = clean_pk_col(T, CMAX_COL)

    (length(auc_R) < 2 || length(auc_T) < 2)  && return nothing
    (length(cmax_R) < 2 || length(cmax_T) < 2) && return nothing

    auc_res  = tost_test_detailed(log.(auc_R),  log.(auc_T))
    cmax_res = tost_test_detailed(log.(cmax_R), log.(cmax_T))

    return (
        auc_gmr=auc_res.gmr, auc_lo=auc_res.ci_lo, auc_hi=auc_res.ci_hi, auc_pass=auc_res.pass,
        cmax_gmr=cmax_res.gmr, cmax_lo=cmax_res.ci_lo, cmax_hi=cmax_res.ci_hi, cmax_pass=cmax_res.pass,
        overall_pass=auc_res.pass && cmax_res.pass
    )
end

# ── %CV Estimation ───────────────────────────────────────────────────────────

"""
Estimate %CV for AUC and Cmax from simulated parallel-design trials.
"""
function estimate_cv(param_samples; n_per_arm=60, n_reps=50)
    auc_vars  = Float64[]
    cmax_vars = Float64[]

    for _ in 1:n_reps
        θ_draw = param_samples[rand(1:length(param_samples))]
        sim_params = merge(θ_draw, (tvbio = 1.0,))

        try
            sim = simulate_be_trial(sim_params; n_per_arm = n_per_arm)
            simdf = DataFrame(sim)
            simdf.route .= "ev"
            simdf.dv .= max.(simdf.dv, 0.0)

            tbl = with_logger(NullLogger()) do
                ncapop = read_nca(simdf;
                    id = :id, time = :time, amt = :amt,
                    observations = :dv, route = :route, group = [:TRT]
                )
                nca_out = run_nca(ncapop)
                nca_out.reportdf
            end

            R = tbl[tbl.TRT .== "R", :]
            T = tbl[tbl.TRT .== "T", :]

            auc_R  = clean_pk_col(R, AUC_COL)
            auc_T  = clean_pk_col(T, AUC_COL)
            cmax_R = clean_pk_col(R, CMAX_COL)
            cmax_T = clean_pk_col(T, CMAX_COL)

            if length(auc_R) >= 2 && length(auc_T) >= 2
                n_r, n_t = length(auc_R), length(auc_T)
                sp² = ((n_r - 1) * var(log.(auc_R)) + (n_t - 1) * var(log.(auc_T))) / (n_r + n_t - 2)
                push!(auc_vars, sp²)
            end
            if length(cmax_R) >= 2 && length(cmax_T) >= 2
                n_r, n_t = length(cmax_R), length(cmax_T)
                sp² = ((n_r - 1) * var(log.(cmax_R)) + (n_t - 1) * var(log.(cmax_T))) / (n_r + n_t - 2)
                push!(cmax_vars, sp²)
            end
        catch
            continue
        end
    end

    auc_cv  = 100 * sqrt(exp(mean(auc_vars)) - 1)
    cmax_cv = 100 * sqrt(exp(mean(cmax_vars)) - 1)
    return (auc_cv=auc_cv, cmax_cv=cmax_cv)
end

# ── Batch Simulation Engine ──────────────────────────────────────────────────

"""
Run batched virtual BE simulations. Returns a DataFrame with one row per
(method, tr_ratio, n_per_arm, batch_id).
"""
function virtual_be_sim_batched(param_samples, method_name;
                                tr_ratios  = TR_RATIOS,
                                n_per_arms = SAMPLE_SIZES,
                                n_batches  = N_BATCHES,
                                n_sims     = N_SIMS_PER_BATCH)

    batch_results = DataFrame(
        method    = String[],
        tr_ratio  = Float64[],
        n_per_arm = Int[],
        batch_id  = Int[],
        n_pass    = Int[],
        n_valid   = Int[],
        pass_rate = Float64[]
    )

    total_cells = length(n_per_arms) * length(tr_ratios)
    cell_count  = 0

    for n_arm in n_per_arms, tr in tr_ratios
        cell_count += 1

        for b in 1:n_batches
            # Seed each batch independently so macro-replications are truly
            # independent.  The seed is a deterministic function of the design
            # cell and batch id, which also makes results reproducible.
            batch_seed = abs(hash((method_name, tr, n_arm, b)) % typemax(Int32))
            Random.seed!(batch_seed)

            n_pass  = 0
            n_error = 0

            for i in 1:n_sims
                θ_draw     = param_samples[rand(1:length(param_samples))]
                sim_params = merge(θ_draw, (tvbio = tr,))

                try
                    sim = simulate_be_trial(sim_params; n_per_arm = n_arm)
                    res = run_be_on_sim_detailed(sim)
                    if res === nothing
                        n_error += 1
                    else
                        res.overall_pass && (n_pass += 1)
                    end
                catch
                    n_error += 1
                end
            end

            n_valid = n_sims - n_error
            rate    = n_valid > 0 ? (n_pass / n_valid) * 100 : NaN

            push!(batch_results, (
                method    = method_name,
                tr_ratio  = tr,
                n_per_arm = n_arm,
                batch_id  = b,
                n_pass    = n_pass,
                n_valid   = n_valid,
                pass_rate = rate
            ))
        end

        println("  [$method_name] ($cell_count/$total_cells) n=$n_arm, T/R=$tr — $n_batches batches done")
    end

    return batch_results
end

"""
Aggregate batch-level results into point estimates with CIs.
Returns a DataFrame with one row per (method, tr_ratio, n_per_arm).
"""
function summarize_batches(batch_df)
    summary = combine(groupby(batch_df, [:method, :tr_ratio, :n_per_arm]),
        :pass_rate => mean   => :mean_rate,
        :pass_rate => std    => :std_rate,
        :pass_rate => length => :n_batches,
        :pass_rate => (x -> quantile(x, 0.025)) => :pctl_025,
        :pass_rate => (x -> quantile(x, 0.975)) => :pctl_975
    )
    # 95% CI for the true rate (SE-based)
    summary.se    = summary.std_rate ./ sqrt.(summary.n_batches)
    summary.ci_lo = max.(summary.mean_rate .- 1.96 .* summary.se, 0.0)
    summary.ci_hi = min.(summary.mean_rate .+ 1.96 .* summary.se, 100.0)
    return summary
end

# ── Formatting Helpers ───────────────────────────────────────────────────────

fmt_ci(m, lo, hi) = "$(round(m, digits=1)) [$(round(lo, digits=1)), $(round(hi, digits=1))]"

function fmt_diff(d, lo, hi)
    s = round(d, digits=1) >= 0 ? "+" : ""
    "$(s)$(round(d, digits=1)) [$(round(lo, digits=1)), $(round(hi, digits=1))]"
end

function fmt_ratio(r, lo, hi)
    isnan(r) && return "—"
    "$(round(r, digits=3)) [$(round(lo, digits=3)), $(round(hi, digits=3))]"
end

# ── Load & Combine Per-Method Results ────────────────────────────────────────

"""
Load the four per-method batch result files, combine them, compute summaries
and relative-to-Fixed metrics. Returns (all_batches, all_summary, relative_df).
"""
function load_all_method_results(; dir = @__DIR__)
    method_files = [
        "batches_fixed.jls",
        "batches_vcov.jls",
        "batches_bootstrap.jls",
        "batches_sir.jls",
    ]

    all_batches = DataFrame()
    for f in method_files
        path = joinpath(dir, f)
        if !isfile(path)
            error("Missing result file: $path — run the corresponding 02_run_*.jl first")
        end
        batches = Serialization.deserialize(path)
        all_batches = nrow(all_batches) == 0 ? batches : vcat(all_batches, batches)
        println("  Loaded $f ($(nrow(batches)) rows)")
    end

    # Aggregate into summaries with CIs
    all_summary = summarize_batches(all_batches)

    # ── Relative-to-Fixed computation ────────────────────────────────────────
    fixed_ref = select(
        filter(r -> r.method == "Fixed", all_summary),
        :tr_ratio, :n_per_arm,
        :mean_rate => :fixed_mean,
        :se        => :fixed_se
    )

    relative_df = innerjoin(
        filter(r -> r.method != "Fixed", all_summary),
        fixed_ref,
        on = [:tr_ratio, :n_per_arm]
    )

    # Absolute difference (pp): method − Fixed
    relative_df.abs_diff       = relative_df.mean_rate .- relative_df.fixed_mean
    relative_df.abs_diff_se    = sqrt.(relative_df.se.^2 .+ relative_df.fixed_se.^2)
    relative_df.abs_diff_ci_lo = relative_df.abs_diff .- 1.96 .* relative_df.abs_diff_se
    relative_df.abs_diff_ci_hi = relative_df.abs_diff .+ 1.96 .* relative_df.abs_diff_se

    # Ratio: method / Fixed  (delta-method SE)
    relative_df.ratio = relative_df.mean_rate ./ max.(relative_df.fixed_mean, 1e-10)
    relative_df.ratio_se = relative_df.ratio .* sqrt.(
        (relative_df.se ./ max.(relative_df.mean_rate, 1e-10)).^2 .+
        (relative_df.fixed_se ./ max.(relative_df.fixed_mean, 1e-10)).^2
    )
    relative_df.ratio_ci_lo = relative_df.ratio .- 1.96 .* relative_df.ratio_se
    relative_df.ratio_ci_hi = relative_df.ratio .+ 1.96 .* relative_df.ratio_se

    # Guard: set ratio columns to NaN where Fixed rate ≈ 0
    for col in [:ratio, :ratio_se, :ratio_ci_lo, :ratio_ci_hi]
        relative_df[!, col] = ifelse.(relative_df.fixed_mean .< 0.1, NaN, relative_df[!, col])
    end

    return (all_batches = all_batches, all_summary = all_summary, relative_df = relative_df)
end
