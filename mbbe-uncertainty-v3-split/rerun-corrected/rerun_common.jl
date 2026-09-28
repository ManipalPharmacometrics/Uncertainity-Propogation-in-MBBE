# Corrected re-run: the delta-method fix to sample_from_vcov now lives in ../common.jl.
include(joinpath(@__DIR__, "..", "common.jl"))
using Base.Threads, CSV

# FIX 2 — AUC endpoint is AUC0-t (auclast, 0–72 h); AUC extrapolated to infinity is not used.
const AUC_CANDIDATES = ["auclast", "aucobs", "auc_tlast", "AUClast"]
