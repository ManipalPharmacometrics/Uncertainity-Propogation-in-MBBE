# Job 11: Bootstrap method @ T/R = 0.70
# Run:  julia distributed_runs/job_11_bootstrap_tr070.jl
# Output: distributed_runs/output/bootstrap_tr070_b001-200.{jls,done}

push!(ARGS, "0.70")
include(joinpath(@__DIR__, "_runner_bootstrap.jl"))
