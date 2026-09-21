# Job 19: SIR method @ T/R = 1.10
# Run:  julia distributed_runs/job_19_sir_tr110.jl
# Output: distributed_runs/output/sir_tr110_b001-200.{jls,done}

push!(ARGS, "1.10")
include(joinpath(@__DIR__, "_runner_sir.jl"))
