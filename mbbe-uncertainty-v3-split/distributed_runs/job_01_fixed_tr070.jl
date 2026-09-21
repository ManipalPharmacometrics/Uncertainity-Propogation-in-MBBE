# Job 01: Fixed method @ T/R = 0.70
# Run:  julia distributed_runs/job_01_fixed_tr070.jl
# Output: distributed_runs/output/fixed_tr070_b001-200.{jls,done}

push!(ARGS, "0.70")
include(joinpath(@__DIR__, "_runner_fixed.jl"))
