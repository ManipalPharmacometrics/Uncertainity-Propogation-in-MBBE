# Job 17: SIR method @ T/R = 0.90
# Run:  julia distributed_runs/job_17_sir_tr090.jl
# Output: distributed_runs/output/sir_tr090_b001-200.{jls,done}

push!(ARGS, "0.90")
include(joinpath(@__DIR__, "_runner_sir.jl"))
