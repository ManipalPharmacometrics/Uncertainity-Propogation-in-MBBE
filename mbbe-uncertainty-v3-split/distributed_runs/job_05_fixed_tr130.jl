# Job 05: Fixed method @ T/R = 1.30
# Run:  julia distributed_runs/job_05_fixed_tr130.jl
# Output: distributed_runs/output/fixed_tr130_b001-200.{jls,done}

push!(ARGS, "1.30")
include(joinpath(@__DIR__, "_runner_fixed.jl"))
