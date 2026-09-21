# Job 04: Fixed method @ T/R = 1.10
# Run:  julia distributed_runs/job_04_fixed_tr110.jl
# Output: distributed_runs/output/fixed_tr110_b001-200.{jls,done}

push!(ARGS, "1.10")
include(joinpath(@__DIR__, "_runner_fixed.jl"))
