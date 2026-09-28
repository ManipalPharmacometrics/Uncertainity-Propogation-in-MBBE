# 01 — refit the reference study and build the four parameter sets (corrected).
include(joinpath(@__DIR__, "rerun_common.jl"))
const OUT = joinpath(@__DIR__, "output")
nt(p) = (tvcl=p.tvcl, tvv=p.tvv, tvka=p.tvka, tvbio=p.tvbio, Ω=Diagonal(collect(diag(Matrix(p.Ω)))), σ=p.σ)

Random.seed!(142)                                   # identical to 01_model_fit.jl
ref_sim = simulate_be_trial(θ_true; n_per_arm = 30)
sim_pop = read_pumas(DataFrame(ref_sim); id=:id, time=:time, observations=[:dv], amt=:amt,
                     evid=:evid, covariates=[:TRT], cmt=:cmt)
init = (tvcl=4.5, tvv=45.0, tvka=1.0, tvbio=1.0, Ω=Diagonal([0.1,0.1,0.1]), σ=0.3)
fpm  = fit(pkmodel, sim_pop, init, Pumas.FOCE(); optim_options=(show_trace=false,))
θ̂    = nt(coef(fpm))
println("θ̂ = ", θ̂); flush(stdout)

rng = Random.Xoshiro(20260927)
vcov_set = sample_from_vcov(θ̂, vcov(fpm), 1000; rng = rng)          # corrected sampler
println("VCOV done"); flush(stdout)

# each piece is checkpointed so a later failure never repeats a finished step
fb = joinpath(OUT, "ckpt_bootstrap.jls")
if isfile(fb)
    boot_set = Serialization.deserialize(fb)
else
    t0 = time()
    bt = infer(fpm, Bootstrap(samples = 500, ensemblealg = EnsembleThreads()))
    boot_set = [nt(coef(f)) for f in bt.vcov.fits if f !== nothing]
    Serialization.serialize(fb, boot_set)
    println("Bootstrap: $(length(boot_set)) usable fits, $(round((time()-t0)/60, digits=1)) min"); flush(stdout)
end

t0 = time()
Random.seed!(20260928)
sir = infer(infer(fpm), SIR(samples = 2000, resamples = 1000))
smp = sir.vcov.samples; rs = sir.vcov.resamples
smp = eltype(smp) <: NamedTuple ? smp : last(smp)     # samples come nested one level deep
sir_set = [nt(p) for p in smp[last(rs)]]
Serialization.serialize(joinpath(OUT, "ckpt_sir.jls"), sir_set)
println("SIR: $(length(sir_set)) resamples, $(round((time()-t0)/60, digits=1)) min"); flush(stdout)

sets = Dict("Fixed" => [θ̂], "VCOV" => vcov_set, "Bootstrap" => boot_set, "SIR" => sir_set)
Serialization.serialize(joinpath(OUT, "parameter_sets.jls"), (θ̂ = θ̂, vcov = Matrix(vcov(fpm)), sets = sets))

# quantile summary per method, for the paper
rows = DataFrame(method=String[], parameter=String[], q025=Float64[], q50=Float64[], q975=Float64[])
getters = ["tvcl"=>p->p.tvcl, "tvv"=>p->p.tvv, "tvka"=>p->p.tvka, "tvbio"=>p->p.tvbio,
           "omega_cl"=>p->p.Ω[1,1], "omega_v"=>p->p.Ω[2,2], "omega_ka"=>p->p.Ω[3,3], "sigma"=>p->p.σ]
for m in ("VCOV", "Bootstrap", "SIR"), (nm, g) in getters
    v = g.(sets[m]); q = quantile(v, [0.025, 0.5, 0.975])
    push!(rows, (m, nm, q...))
end
CSV.write(joinpath(OUT, "parameter_set_quantiles.csv"), rows)
println(rows); println("01 DONE")
