# 05 — realised test/reference GMRs at T/R = 1.30 and 0.70 under Fixed: explains the flat false-acceptance rate at 1.30.
include(joinpath(@__DIR__, "rerun_common.jl"))
global AUC_COL  = :auclast
global CMAX_COL = :cmax
θ̂ = Serialization.deserialize(joinpath(@__DIR__, "output", "parameter_sets.jls")).θ̂
for tr in (0.7, 1.3), n in (24, 80)
    Random.seed!(hash((tr, n)) % typemax(Int32))
    rs = filter(!isnothing, [run_be_on_sim_detailed(simulate_be_trial(merge(θ̂, (tvbio = tr,)); n_per_arm = n)) for _ in 1:200])
    q(v) = round.(quantile(v, [0.5, 0.05, 0.95]), digits = 3)
    println("T/R=$tr n=$n  AUC0-t GMR med/5/95 = $(q([r.auc_gmr for r in rs]))  Cmax GMR = $(q([r.cmax_gmr for r in rs]))  ",
            "AUC pass $(round(100mean(r.auc_pass for r in rs),digits=1))%  Cmax pass $(round(100mean(r.cmax_pass for r in rs),digits=1))%  both $(round(100mean(r.overall_pass for r in rs),digits=1))%")
end
