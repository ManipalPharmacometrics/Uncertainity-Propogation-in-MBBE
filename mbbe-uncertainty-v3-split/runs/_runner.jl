# ──────────────────────────────────────────────────────────────────────────────
# _runner.jl — Shared engine for individual job files
#
# Resolves job parameters from (in priority order):
#   1. Constants defined by job_XX_*.jl (JOB_METHOD, JOB_TR_RATIO, etc.)
#   2. CLI arguments:   julia runs/_runner.jl Fixed 0.90
#   3. ENV variables:   JOB_METHOD=Fixed JOB_TR_RATIO=0.90
#   4. No arguments:    runs ALL method × ratio combinations (JuliaHub mode)
#
# This script:
#   1. Loads common.jl and model_fit_results.jls
#   2. Runs virtual_be_sim_batched for the resolved (method, T/R ratio) pair(s)
#   3. Saves output to runs/output/job_XX_<tag>.jls
#   4. Writes a .done marker with timing info
# ──────────────────────────────────────────────────────────────────────────────

# ── Resolve job parameters ───────────────────────────────────────────────────
_VALID_METHODS = ["Fixed", "VCOV", "Bootstrap", "SIR"]
_VALID_RATIOS  = [0.70, 0.90, 1.00, 1.10, 1.30]
_RUN_ALL = false

if !isdefined(Main, :JOB_METHOD)
    if length(ARGS) >= 2
        # Mode 2: CLI arguments
        global JOB_METHOD   = ARGS[1]
        global JOB_TR_RATIO = parse(Float64, ARGS[2])
    elseif haskey(ENV, "JOB_METHOD") && haskey(ENV, "JOB_TR_RATIO")
        # Mode 3: Environment variables (JuliaHub batch jobs)
        global JOB_METHOD   = ENV["JOB_METHOD"]
        global JOB_TR_RATIO = parse(Float64, ENV["JOB_TR_RATIO"])
    else
        # Mode 4: No arguments — run all combinations
        _RUN_ALL = true
        println("No method/ratio specified — running ALL $(length(_VALID_METHODS) * length(_VALID_RATIOS)) combinations.")
    end

    if !_RUN_ALL
        if JOB_METHOD ∉ _VALID_METHODS
            error("Unknown method '$(JOB_METHOD)'. Must be one of: $(join(_VALID_METHODS, ", "))")
        end
        if JOB_TR_RATIO ∉ _VALID_RATIOS
            error("Unknown T/R ratio $(JOB_TR_RATIO). Must be one of: $(join(_VALID_RATIOS, ", "))")
        end

        # Auto-derive JOB_ID and JOB_TAG
        _method_offset = Dict("Fixed" => 0, "VCOV" => 5, "Bootstrap" => 10, "SIR" => 15)
        _ratio_index   = findfirst(==(JOB_TR_RATIO), _VALID_RATIOS)
        global JOB_ID  = _method_offset[JOB_METHOD] + _ratio_index
        global JOB_TAG = lowercase(JOB_METHOD) * "_tr" * replace(string(Int(round(JOB_TR_RATIO * 100))), "." => "")
    end
end

include(joinpath(@__DIR__, "..", "common.jl"))
using Dates

# ── Configuration ─────────────────────────────────────────────────────────────
# Reduced sample sizes: removed 12 and 50 (5 values instead of 7)
const JOB_SAMPLE_SIZES = [24, 30, 40, 60, 80]

# ── Load Model Fit Results ────────────────────────────────────────────────────
fit_data = Serialization.deserialize(joinpath(@__DIR__, "..", "model_fit_results.jls"))

global AUC_COL  = fit_data[:auc_col]
global CMAX_COL = fit_data[:cmax_col]

_param_key_map = Dict(
    "Fixed"     => :fixed_param_samples,
    "VCOV"      => :vcov_param_samples,
    "Bootstrap" => :boot_param_samples,
    "SIR"       => :sir_param_samples,
)

# ── Helper: run and save one (method, ratio) job ─────────────────────────────
function _run_single_job(method, tr_ratio, job_id, job_tag, param_samples)
    outdir   = joinpath(@__DIR__, "output")
    outfile  = joinpath(outdir, "job_$(lpad(job_id, 2, '0'))_$(job_tag).jls")
    donefile = joinpath(outdir, "job_$(lpad(job_id, 2, '0'))_$(job_tag).done")

    println("=" ^ 70)
    println("JOB $job_id: $method @ T/R = $tr_ratio")
    println("=" ^ 70)
    println("  Sample sizes : $(JOB_SAMPLE_SIZES)")
    println("  Batches      : $N_BATCHES")
    println("  Trials/batch : $N_SIMS_PER_BATCH")
    println("  Total sims   : $(length(JOB_SAMPLE_SIZES) * N_BATCHES * N_SIMS_PER_BATCH)")
    println("  Output       : $outfile")
    println("  AUC col      : $AUC_COL")
    println("  Cmax col     : $CMAX_COL")
    println("-" ^ 70)

    t_start = time()

    batches = virtual_be_sim_batched(
        param_samples, method;
        tr_ratios  = [tr_ratio],
        n_per_arms = JOB_SAMPLE_SIZES,
        n_batches  = N_BATCHES,
        n_sims     = N_SIMS_PER_BATCH
    )

    elapsed = time() - t_start
    elapsed_min = round(elapsed / 60, digits=1)

    mkpath(outdir)
    Serialization.serialize(outfile, batches)
    println("\nSaved $(nrow(batches)) rows to: $outfile")

    open(donefile, "w") do io
        println(io, "job_id     : $job_id")
        println(io, "method     : $method")
        println(io, "tr_ratio   : $tr_ratio")
        println(io, "n_rows     : $(nrow(batches))")
        println(io, "elapsed_s  : $(round(elapsed, digits=1))")
        println(io, "elapsed_min: $elapsed_min")
        println(io, "completed  : $(Dates.now())")
    end

    println("Done! Job $job_id completed in $(elapsed_min) minutes.")
    println("=" ^ 70)

    return batches
end

# ── Execute ──────────────────────────────────────────────────────────────────
if _RUN_ALL
    # Run every method × ratio combination
    _method_offset = Dict("Fixed" => 0, "VCOV" => 5, "Bootstrap" => 10, "SIR" => 15)
    total_jobs = length(_VALID_METHODS) * length(_VALID_RATIOS)
    job_count  = 0

    for method in _VALID_METHODS
        global param_samples = fit_data[_param_key_map[method]]
        for (ri, tr) in enumerate(_VALID_RATIOS)
            global job_count += 1
            job_id  = _method_offset[method] + ri
            job_tag = lowercase(method) * "_tr" * replace(string(Int(round(tr * 100))), "." => "")
            println("\n>>> Starting job $job_count / $total_jobs <<<\n")
            _run_single_job(method, tr, job_id, job_tag, param_samples)
        end
    end
    println("\n✓ All $total_jobs jobs completed.")
else
    # Run the single resolved job
    param_samples = fit_data[_param_key_map[JOB_METHOD]]
    _run_single_job(JOB_METHOD, JOB_TR_RATIO, JOB_ID, JOB_TAG, param_samples)
end
