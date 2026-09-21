# Job 03: Fixed method @ T/R = 1.00
# Run:  julia distributed_runs/job_03_fixed_tr100.jl
# Output: distributed_runs/output/fixed_tr100_b001-200.{jls,done}

push!(ARGS, "1.00")
include(joinpath(@__DIR__, "_runner_fixed.jl"))
