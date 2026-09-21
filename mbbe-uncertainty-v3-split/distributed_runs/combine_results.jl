# ──────────────────────────────────────────────────────────────────────────────
# combine_results.jl — Merge all distributed job outputs into final .jls files
#
# Scans distributed_runs/output/ for all .jls files, groups by method,
# deduplicates overlapping batch ranges, and saves:
#   batches_fixed.jls, batches_vcov.jls, batches_bootstrap.jls, batches_sir.jls
# to the parent directory (compatible with 03_results_tables.qmd / 04_visualizations.qmd).
#
# Run:  julia distributed_runs/combine_results.jl
# ──────────────────────────────────────────────────────────────────────────────

using Serialization
using DataFrames

outdir    = joinpath(@__DIR__, "output")
parentdir = joinpath(@__DIR__, "..")

if !isdir(outdir)
    error("No output directory found at: $outdir — run some jobs first.")
end

jls_files = filter(f -> endswith(f, ".jls"), readdir(outdir))

if isempty(jls_files)
    error("No .jls files found in: $outdir")
end

println("Combining distributed job outputs...")
println()

# Method → output filename
method_filenames = Dict(
    "Fixed"     => "batches_fixed.jls",
    "VCOV"      => "batches_vcov.jls",
    "Bootstrap" => "batches_bootstrap.jls",
    "SIR"       => "batches_sir.jls",
)

# Load all .jls files and concatenate
method_dfs = Dict{String, DataFrame}()
n_loaded = 0

for f in sort(jls_files)
    filepath = joinpath(outdir, f)
    df = Serialization.deserialize(filepath)

    if nrow(df) == 0
        println("  SKIP $f (empty)")
        continue
    end

    method = df.method[1]
    if haskey(method_dfs, method)
        method_dfs[method] = vcat(method_dfs[method], df)
    else
        method_dfs[method] = df
    end
    n_loaded += 1
    println("  Loaded $f — $(nrow(df)) rows ($(method))")
end

println()

if n_loaded == 0
    error("No valid result files found.")
end

# Deduplicate: if overlapping batch ranges were run, keep only unique
# (method, tr_ratio, n_per_arm, batch_id) rows. In case of duplicates,
# keep the last occurrence (most recent run).
println("Deduplicating overlapping batch ranges...")

for (method, df) in method_dfs
    before = nrow(df)
    # Keep last occurrence of each unique key
    df.row_idx = 1:nrow(df)
    deduped = combine(
        groupby(df, [:method, :tr_ratio, :n_per_arm, :batch_id]),
        last
    )
    select!(deduped, Not(:row_idx))
    method_dfs[method] = sort(deduped, [:tr_ratio, :n_per_arm, :batch_id])
    after = nrow(method_dfs[method])
    if before != after
        println("  $method: $before → $after rows ($(before - after) duplicates removed)")
    end
end

println()

# Save combined results
missing_methods = String[]

for (method, filename) in method_filenames
    if !haskey(method_dfs, method)
        push!(missing_methods, method)
        println("  SKIP $filename — no data for $method")
        continue
    end
    df = method_dfs[method]
    outpath = joinpath(parentdir, filename)
    Serialization.serialize(outpath, df)
    println("  Saved $filename — $(nrow(df)) rows")
end

println()

if !isempty(missing_methods)
    println("WARNING: Missing data for: $(join(missing_methods, ", "))")
    println("Run the corresponding _runner_*.jl scripts first.")
    println()
end

println("Done! Combined results saved to: $parentdir")
println("You can now run 03_results_tables.qmd and 04_visualizations.qmd.")
