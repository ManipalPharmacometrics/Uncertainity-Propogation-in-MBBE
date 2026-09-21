# Job 15: Bootstrap method @ T/R = 1.30
# Run:  julia distributed_runs/job_15_bootstrap_tr130.jl
# Output: distributed_runs/output/bootstrap_tr130_b001-200.{jls,done}

push!(ARGS, "1.30")
include(joinpath(@__DIR__, "_runner_bootstrap.jl"))
