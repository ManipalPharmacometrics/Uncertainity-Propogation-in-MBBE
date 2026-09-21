# ──────────────────────────────────────────────────────────────────────────────
# _base_runner.jl — Shared engine for distributed per-method runners
#
# This file is included by _runner_fixed.jl, _runner_vcov.jl, etc.
# The caller MUST define `const DIST_METHOD` before including this file.
#
# Resolves job parameters from (in priority order):
#   1. CLI arguments:
#        julia _runner_<method>.jl                     # all ratios, all batches
#        julia _runner_<method>.jl 0.90                # one ratio, all batches
#        julia _runner_<method>.jl 0.90 1 50           # one ratio, batches 1-50
#   2. ENV variables:
#        JOB_TR_RATIO=0.90 BATCH_START=1 BATCH_END=50 julia _runner_<method>.jl
#
# Key difference from runs/_runner.jl:
#   - Each method has its own runner file (no method resolution at runtime)
#   - Supports batch-range slicing for distributed execution
#   - Every individual trial is seeded independently:
#       seed = hash((method, tr_ratio, n_arm, batch_id, trial_index))
#     so results are reproducible regardless of how batches are split
#
# Output: distributed_runs/output/<method>_tr<ratio>_b<start>-<end>.jls
# ──────────────────────────────────────────────────────────────────────────────

if !isdefined(Main, :DIST_METHOD)
    error("DIST_METHOD must be defined before including _base_runner.jl")
end

# ── Resolve parameters ───────────────────────────────────────────────────────
_VALID_METHODS = ["Fixed", "VCOV", "Bootstrap", "SIR"]
_VALID_RATIOS  = [0.70, 0.90, 1.00, 1.10, 1.30]

if DIST_METHOD ∉ _VALID_METHODS
    error("Unknown method '$(DIST_METHOD)'. Must be one of: $(join(_VALID_METHODS, ", "))")
end

# Defaults
_run_ratios   = _VALID_RATIOS   # which T/R ratios to run
_batch_start  = 1
_batch_end    = nothing          # resolved after loading common.jl (uses N_BATCHES)

if length(ARGS) >= 3
    # Mode: julia _runner_<method>.jl <tr_ratio> <batch_start> <batch_end>
    _run_ratios  = [parse(Float64, ARGS[1])]
    _batch_start = parse(Int, ARGS[2])
    _batch_end   = parse(Int, ARGS[3])
elseif length(ARGS) == 2
    # Mode: julia _runner_<method>.jl <tr_ratio> <batch_start>
    #   (single batch — handy for per-batch distribution)
    _run_ratios  = [parse(Float64, ARGS[1])]
    _batch_start = parse(Int, ARGS[2])
    _batch_end   = parse(Int, ARGS[2])
elseif length(ARGS) == 1
    # Mode: julia _runner_<method>.jl <tr_ratio>
    _run_ratios = [parse(Float64, ARGS[1])]
else
    # Check env vars
    if haskey(ENV, "JOB_TR_RATIO")
        _run_ratios = [parse(Float64, ENV["JOB_TR_RATIO"])]
    end
    if haskey(ENV, "BATCH_START")
        _batch_start = parse(Int, ENV["BATCH_START"])
    end
    if haskey(ENV, "BATCH_END")
        _batch_end = parse(Int, ENV["BATCH_END"])
    end
end

# Validate ratios
for r in _run_ratios
    if r ∉ _VALID_RATIOS
        error("Unknown T/R ratio $(r). Must be one of: $(join(_VALID_RATIOS, ", "))")
    end
end

# ── Load shared code ─────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "..", "common.jl"))
using Dates

# Resolve batch_end now that N_BATCHES is available
if _batch_end === nothing
    _batch_end = N_BATCHES
end

if _batch_start < 1 || _batch_end < _batch_start
    error("Invalid batch range: $(_batch_start)-$(_batch_end)")
end
if _batch_end > N_BATCHES
    @warn "batch_end ($(_batch_end)) exceeds N_BATCHES ($N_BATCHES) — clamping."
    _batch_end = N_BATCHES
end

# ── Configuration ────────────────────────────────────────────────────────────
const DIST_SAMPLE_SIZES = SAMPLE_SIZES   # use full set from common.jl

# ── Load Model Fit Results ───────────────────────────────────────────────────
fit_data = Serialization.deserialize(joinpath(@__DIR__, "..", "model_fit_results.jls"))

global AUC_COL  = fit_data[:auc_col]
global CMAX_COL = fit_data[:cmax_col]

_param_key_map = Dict(
    "Fixed"     => :fixed_param_samples,
    "VCOV"      => :vcov_param_samples,
    "Bootstrap" => :boot_param_samples,
    "SIR"       => :sir_param_samples,
)

param_samples = fit_data[_param_key_map[DIST_METHOD]]

# ── Simulation Engine (per-trial seeding) ────────────────────────────────────

"""
Run batched virtual BE simulations with per-trial independent seeding.
Each individual trial gets its own seed derived from:
    hash((method, tr_ratio, n_per_arm, batch_id, trial_index))
This ensures full reproducibility regardless of batch splitting.
"""
function distributed_be_sim_batched(param_samples, method_name;
                                     tr_ratios   = _run_ratios,
                                     n_per_arms  = DIST_SAMPLE_SIZES,
                                     batch_start = _batch_start,
                                     batch_end   = _batch_end,
                                     n_sims      = N_SIMS_PER_BATCH)

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
    n_batches_run = batch_end - batch_start + 1

    for n_arm in n_per_arms, tr in tr_ratios
        cell_count += 1

        for b in batch_start:batch_end
            n_pass  = 0
            n_error = 0

            for i in 1:n_sims
                # Per-trial seed: every single simulation is independently seeded.
                # Use the full UInt64 hash space (~1.8e19 values) so that collisions
                # are statistically impossible across every (method, tr, n_arm,
                # batch, trial) tuple — no two sims ever share a seed.
                trial_seed = hash((method_name, tr, n_arm, b, i))  # UInt64
                Random.seed!(trial_seed)

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

        println("  [$method_name] ($cell_count/$total_cells) n=$n_arm, T/R=$tr — $n_batches_run batches done (b$batch_start-b$batch_end)")
    end

    return batch_results
end

# ── Helper: format T/R ratio tag ─────────────────────────────────────────────
_ratio_tag(tr) = "tr" * replace(string(Int(round(tr * 100))), "." => "")

# ── Helper: run and save one T/R ratio ───────────────────────────────────────
function _run_and_save(method, tr_ratio, param_samples, batch_start, batch_end)
    tag     = lowercase(method) * "_" * _ratio_tag(tr_ratio)
    brange  = "b$(lpad(batch_start, 3, '0'))-$(lpad(batch_end, 3, '0'))"
    outdir  = joinpath(@__DIR__, "output")
    outfile  = joinpath(outdir, "$(tag)_$(brange).jls")
    donefile = joinpath(outdir, "$(tag)_$(brange).done")

    n_batches_run = batch_end - batch_start + 1

    println("=" ^ 70)
    println("METHOD      : $method")
    println("T/R ratio   : $tr_ratio")
    println("Batches     : $batch_start - $batch_end ($n_batches_run batches)")
    println("Sample sizes: $(DIST_SAMPLE_SIZES)")
    println("Trials/batch: $N_SIMS_PER_BATCH")
    println("Total sims  : $(length(DIST_SAMPLE_SIZES) * n_batches_run * N_SIMS_PER_BATCH)")
    println("Output      : $outfile")
    println("AUC col     : $AUC_COL")
    println("Cmax col    : $CMAX_COL")
    println("-" ^ 70)

    t_start = time()

    batches = distributed_be_sim_batched(
        param_samples, method;
        tr_ratios   = [tr_ratio],
        n_per_arms  = DIST_SAMPLE_SIZES,
        batch_start = batch_start,
        batch_end   = batch_end,
        n_sims      = N_SIMS_PER_BATCH
    )

    elapsed = time() - t_start
    elapsed_min = round(elapsed / 60, digits=1)

    mkpath(outdir)
    Serialization.serialize(outfile, batches)
    println("\nSaved $(nrow(batches)) rows to: $outfile")

    open(donefile, "w") do io
        println(io, "method      : $method")
        println(io, "tr_ratio    : $tr_ratio")
        println(io, "batch_start : $batch_start")
        println(io, "batch_end   : $batch_end")
        println(io, "n_rows      : $(nrow(batches))")
        println(io, "elapsed_s   : $(round(elapsed, digits=1))")
        println(io, "elapsed_min : $elapsed_min")
        println(io, "completed   : $(Dates.now())")
    end

    println("Done! $method @ T/R=$tr_ratio (b$batch_start-b$batch_end) completed in $elapsed_min minutes.")
    println("=" ^ 70)

    return batches
end

# ── Execute ──────────────────────────────────────────────────────────────────
if length(_run_ratios) > 1
    println("\nRunning $(DIST_METHOD) for all $(length(_run_ratios)) T/R ratios, batches $(_batch_start)-$(_batch_end)\n")
end

for tr in _run_ratios
    _run_and_save(DIST_METHOD, tr, param_samples, _batch_start, _batch_end)
end

if length(_run_ratios) > 1
    println("\n✓ All $(length(_run_ratios)) T/R ratios completed for $(DIST_METHOD).")
end
