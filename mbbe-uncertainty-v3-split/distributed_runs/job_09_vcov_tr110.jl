# Job 09: VCOV method @ T/R = 1.10
# Run:  julia distributed_runs/job_09_vcov_tr110.jl
# Output: distributed_runs/output/vcov_tr110_b001-200.{jls,done}

push!(ARGS, "1.10")
include(joinpath(@__DIR__, "_runner_vcov.jl"))
