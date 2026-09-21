# ──────────────────────────────────────────────────────────────────────────────
# check_status.jl — Show completion status of distributed jobs
#
# Scans distributed_runs/output/ for .done markers and reports:
#   - Which (method, ratio, batch-range) jobs are complete
#   - Total batch coverage per (method, ratio)
#   - Overall progress
#
# Run:  julia distributed_runs/check_status.jl
# ──────────────────────────────────────────────────────────────────────────────

using Dates

outdir = joinpath(@__DIR__, "output")

if !isdir(outdir)
    println("No output directory found at: $outdir")
    println("No jobs have been run yet.")
    exit(0)
end

# Parse .done files
done_files = filter(f -> endswith(f, ".done"), readdir(outdir))

if isempty(done_files)
    println("No completed jobs found in: $outdir")
    exit(0)
end

println("=" ^ 78)
println("  DISTRIBUTED JOB STATUS — $(Dates.now())")
println("=" ^ 78)
println()

# Collect info from .done files
struct JobInfo
    method::String
    tr_ratio::Float64
    batch_start::Int
    batch_end::Int
    elapsed_min::Float64
    completed::String
    filename::String
end

jobs = JobInfo[]

for f in sort(done_files)
    filepath = joinpath(outdir, f)
    info = Dict{String,String}()
    for line in readlines(filepath)
        if contains(line, ":")
            key, val = split(line, ":"; limit=2)
            info[strip(key)] = strip(val)
        end
    end

    push!(jobs, JobInfo(
        get(info, "method", "?"),
        parse(Float64, get(info, "tr_ratio", "0")),
        parse(Int, get(info, "batch_start", "0")),
        parse(Int, get(info, "batch_end", "0")),
        let v = tryparse(Float64, get(info, "elapsed_min", "0")); v === nothing ? 0.0 : v end,
        get(info, "completed", "?"),
        f
    ))
end

# Group by (method, ratio) and show coverage
methods = ["Fixed", "VCOV", "Bootstrap", "SIR"]
ratios  = [0.70, 0.90, 1.00, 1.10, 1.30]

total_elapsed = 0.0
n_done = 0

for method in methods
    method_jobs = filter(j -> j.method == method, jobs)
    if isempty(method_jobs)
        continue
    end

    println("  $(method)")
    println("  " * "-" ^ 74)

    for tr in ratios
        ratio_jobs = filter(j -> j.tr_ratio == tr, method_jobs)
        if isempty(ratio_jobs)
            println("    T/R=$tr  · PENDING")
            continue
        end

        # Compute batch coverage
        covered = Set{Int}()
        elapsed_sum = 0.0
        for j in ratio_jobs
            for b in j.batch_start:j.batch_end
                push!(covered, b)
            end
            elapsed_sum += j.elapsed_min
            total_elapsed += j.elapsed_min
            n_done += 1
        end

        n_covered = length(covered)
        ranges = join(["b$(j.batch_start)-$(j.batch_end)" for j in sort(ratio_jobs, by=j->j.batch_start)], ", ")
        status = n_covered >= 200 ? "✓ COMPLETE" : "… PARTIAL ($n_covered/200 batches)"

        println("    T/R=$tr  $status  [$ranges]  $(round(elapsed_sum, digits=1)) min")
    end
    println()
end

println("-" ^ 78)
println("  Total job files: $n_done")
if n_done > 0
    println("  Total compute time: $(round(total_elapsed, digits=1)) min ($(round(total_elapsed/60, digits=1)) hours)")
end
println("=" ^ 78)
