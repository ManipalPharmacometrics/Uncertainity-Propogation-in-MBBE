# Job 07: VCOV method @ T/R = 0.90
# Run:  julia distributed_runs/job_07_vcov_tr090.jl
# Output: distributed_runs/output/vcov_tr090_b001-200.{jls,done}

push!(ARGS, "0.90")
include(joinpath(@__DIR__, "_runner_vcov.jl"))
