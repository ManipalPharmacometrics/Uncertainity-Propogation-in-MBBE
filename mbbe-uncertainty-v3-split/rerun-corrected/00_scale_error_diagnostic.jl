# Diagnostic: is the type I error inflation caused by treating a natural-scale vcov as log-scale?
cd(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "common.jl"))
ENV["SCR"] = get(ENV, "SCR", joinpath(@__DIR__, "output"))
using Base.Threads
d = Serialization.deserialize("model_fit_results.jls")
global AUC_COL = d[:auc_col]; global CMAX_COL = d[:cmax_col]
println("NCA cols: ", AUC_COL, " ", CMAX_COL)

# 1. Re-fit exactly as 01_model_fit.jl did
Random.seed!(142)
ref_sim = simulate_be_trial(θ_true; n_per_arm = 30)
sim_pop = read_pumas(DataFrame(ref_sim); id=:id, time=:time, observations=[:dv], amt=:amt, evid=:evid, covariates=[:TRT], cmt=:cmt)
init = (tvcl=4.5, tvv=45.0, tvka=1.0, tvbio=1.0, Ω=Diagonal([0.1,0.1,0.1]), σ=0.3)
fit_res = fit(pkmodel, sim_pop, init, Pumas.FOCE())
θ̂ = coef(fit_res)
println("refit θ̂: ", θ̂)
println("stored θ̂: ", d[:θ̂])
V = vcov(fit_res)
inf = infer(fit_res)
println(coeftable(inf))
println("sqrt(diag(vcov)) = ", sqrt.(diag(V)))

# 2. Correct log-scale covariance by delta method
θvec = [θ̂.tvcl, θ̂.tvv, θ̂.tvka, θ̂.tvbio, θ̂.Ω[1,1], θ̂.Ω[2,2], θ̂.Ω[3,3], θ̂.σ]
Dinv = Diagonal(1 ./ θvec)
Vlog = Symmetric(Matrix(Dinv * V * Dinv))
Random.seed!(2026)
corr_samples = sample_from_vcov(θ̂, V, 200)      # common.jl now applies the delta method itself
bug_samples = d[:vcov_param_samples]              # as-run draws stored by the original 01_model_fit.jl
for (nm, s) in (("buggy", bug_samples), ("corrected", corr_samples))
  for p in (:tvcl, :tvv, :tvka, :σ)
    v = [x[p] for x in s]; println(nm, " ", p, " 2.5/50/97.5: ", round.(quantile(v,[0.025,0.5,0.975]), sigdigits=3))
  end
  v = [x.Ω[2,2] for x in s]; println(nm, " Ω22 ", round.(quantile(v,[0.025,0.5,0.975]), sigdigits=3))
end
Serialization.serialize(joinpath(ENV["SCR"], "corr_samples.jls"), (θ̂=θ̂, V=Matrix(V), Vlog=Matrix(Vlog), corr=corr_samples))

# 3. Mini-grid with per-trial logging
t0 = time(); simulate_be_trial(θ̂; n_per_arm=80) |> run_be_on_sim_detailed; println("1 trial n=80: ", time()-t0, " s")
NT = parse(Int, get(ENV, "NTRIALS", "600"))
rows = DataFrame(method=String[], tr=Float64[], n=Int[], tvv=Float64[], valid=Bool[], pass=Bool[], auc_gmr=Float64[], cmax_gmr=Float64[])
lk = ReentrantLock()
sets = [("Fixed", fill(θ̂, 1)), ("VCOV_asrun", bug_samples), ("VCOV_corrected", corr_samples)]
for (nm, s) in sets, tr in (0.7, 1.0, 1.3), n in (24, 80)
  @threads for i in 1:NT
    rng = Random.Xoshiro(hash((nm, tr, n, i)))
    θd = s[rand(rng, 1:length(s))]
    p = merge(θd, (tvbio = tr,))
    local r
    try
      Random.seed!(hash((nm,tr,n,i)) % typemax(Int32))
      r = run_be_on_sim_detailed(simulate_be_trial(p; n_per_arm=n))
    catch
      r = nothing
    end
    lock(lk) do
      if r === nothing
        push!(rows, (nm, tr, n, θd.tvv, false, false, NaN, NaN))
      else
        push!(rows, (nm, tr, n, θd.tvv, true, r.overall_pass, r.auc_gmr, r.cmax_gmr))
      end
    end
  end
  sub = rows[(rows.method .== nm) .& (rows.tr .== tr) .& (rows.n .== n) .& rows.valid, :]
  println(rpad(nm,15), " T/R=$tr n=$n  pass=", round(100*mean(sub.pass), digits=2), "%  (valid ", nrow(sub), ")  t=", round(time()-t0))
  flush(stdout)
end
CSV_PATH = joinpath(ENV["SCR"], "diag_trials.jls")
Serialization.serialize(CSV_PATH, rows)
println("DONE")
