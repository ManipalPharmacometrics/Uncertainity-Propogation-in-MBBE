# Job 14: Bootstrap method @ T/R = 1.10
# Run:  julia distributed_runs/job_14_bootstrap_tr110.jl
# Output: distributed_runs/output/bootstrap_tr110_b001-200.{jls,done}

push!(ARGS, "1.10")
include(joinpath(@__DIR__, "_runner_bootstrap.jl"))
