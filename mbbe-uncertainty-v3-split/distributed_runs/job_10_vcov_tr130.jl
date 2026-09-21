# Job 10: VCOV method @ T/R = 1.30
# Run:  julia distributed_runs/job_10_vcov_tr130.jl
# Output: distributed_runs/output/vcov_tr130_b001-200.{jls,done}

push!(ARGS, "1.30")
include(joinpath(@__DIR__, "_runner_vcov.jl"))
