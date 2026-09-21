# Job 08: VCOV method @ T/R = 1.00
# Run:  julia distributed_runs/job_08_vcov_tr100.jl
# Output: distributed_runs/output/vcov_tr100_b001-200.{jls,done}

push!(ARGS, "1.00")
include(joinpath(@__DIR__, "_runner_vcov.jl"))
