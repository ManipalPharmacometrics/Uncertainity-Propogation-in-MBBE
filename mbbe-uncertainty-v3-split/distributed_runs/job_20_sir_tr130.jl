# Job 20: SIR method @ T/R = 1.30
# Run:  julia distributed_runs/job_20_sir_tr130.jl
# Output: distributed_runs/output/sir_tr130_b001-200.{jls,done}

push!(ARGS, "1.30")
include(joinpath(@__DIR__, "_runner_sir.jl"))
