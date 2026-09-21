# Job 13: Bootstrap method @ T/R = 1.00
# Run:  julia distributed_runs/job_13_bootstrap_tr100.jl
# Output: distributed_runs/output/bootstrap_tr100_b001-200.{jls,done}

push!(ARGS, "1.00")
include(joinpath(@__DIR__, "_runner_bootstrap.jl"))
