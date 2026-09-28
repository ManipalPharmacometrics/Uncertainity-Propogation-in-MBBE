# 04 — between-subject %CV of AUC0–t and Cmax under the plug-in estimates (same settings as the original 01_model_fit.jl).
include(joinpath(@__DIR__, "rerun_common.jl"))
global AUC_COL  = :auclast
global CMAX_COL = :cmax
θ̂ = Serialization.deserialize(joinpath(@__DIR__, "output", "parameter_sets.jls")).θ̂
Random.seed!(20260928)
cv = estimate_cv([θ̂]; n_per_arm = 60, n_reps = 100)
println("AUC0-t %CV = $(round(cv.auc_cv, digits=1))%   Cmax %CV = $(round(cv.cmax_cv, digits=1))%")
