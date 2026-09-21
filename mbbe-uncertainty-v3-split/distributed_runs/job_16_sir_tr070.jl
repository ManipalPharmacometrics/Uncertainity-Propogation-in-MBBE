# Job 16: SIR method @ T/R = 0.70
# Run:  julia distributed_runs/job_16_sir_tr070.jl
# Output: distributed_runs/output/sir_tr070_b001-200.{jls,done}

push!(ARGS, "0.70")
include(joinpath(@__DIR__, "_runner_sir.jl"))
