# ──────────────────────────────────────────────────────────────────────────────
# _runner_fixed.jl — Distributed runner for the Fixed (plug-in) method
#
# Usage:
#   julia _runner_fixed.jl                     # all T/R ratios, all batches
#   julia _runner_fixed.jl 0.90                # T/R=0.90, all batches
#   julia _runner_fixed.jl 0.90 1 50           # T/R=0.90, batches 1–50
#   julia _runner_fixed.jl 0.90 25             # T/R=0.90, batch 25 only
#
# Environment variables (override defaults, overridden by CLI args):
#   JOB_TR_RATIO=0.90 BATCH_START=1 BATCH_END=50 julia _runner_fixed.jl
# ──────────────────────────────────────────────────────────────────────────────

const DIST_METHOD = "Fixed"
include(joinpath(@__DIR__, "_base_runner.jl"))
