# Job 02: Fixed method @ T/R = 0.90
# Run:  julia distributed_runs/job_02_fixed_tr090.jl
# Output: distributed_runs/output/fixed_tr090_b001-200.{jls,done}

push!(ARGS, "0.90")
include(joinpath(@__DIR__, "_runner_fixed.jl"))
