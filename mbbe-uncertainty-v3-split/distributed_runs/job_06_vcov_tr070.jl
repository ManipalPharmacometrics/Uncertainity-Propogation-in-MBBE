# Job 06: VCOV method @ T/R = 0.70
# Run:  julia distributed_runs/job_06_vcov_tr070.jl
# Output: distributed_runs/output/vcov_tr070_b001-200.{jls,done}

push!(ARGS, "0.70")
include(joinpath(@__DIR__, "_runner_vcov.jl"))
