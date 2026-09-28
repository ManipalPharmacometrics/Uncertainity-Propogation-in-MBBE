# 02 — virtual BE grid, all four methods, checkpointed per design cell (resumable).
include(joinpath(@__DIR__, "rerun_common.jl"))
const OUT  = joinpath(@__DIR__, "output")
const CELL = joinpath(OUT, "cells"); mkpath(CELL)
global AUC_COL  = :auclast
global CMAX_COL = :cmax

const NS   = [12, 24, 30, 40, 50, 60, 80]
const TRS  = [0.7, 0.9, 1.0, 1.1, 1.3]
const METH = ["Fixed", "VCOV", "Bootstrap", "SIR"]
const B    = parse(Int, get(ENV, "N_BATCHES", "20"))
const M    = parse(Int, get(ENV, "N_TRIALS",  "100"))

sets = Serialization.deserialize(joinpath(OUT, "parameter_sets.jls")).sets
t_start = time()
println("Grid: $(length(NS)) n × $(length(TRS)) T/R × $(length(METH)) methods, $B × $M trials/cell, $(nthreads()) threads"); flush(stdout)

for n in NS, tr in TRS, meth in METH
    f = joinpath(CELL, "$(meth)_tr$(tr)_n$(n).csv")
    isfile(f) && continue
    ps = sets[meth]
    res = DataFrame(method=String[], tr_ratio=Float64[], n_per_arm=Int[], batch_id=Int[],
                    n_pass=Int[], n_valid=Int[], pass_rate=Float64[])
    t0 = time()
    for b in 1:B
        pass = zeros(Bool, M); ok = zeros(Bool, M)
        @threads for i in 1:M
            s = hash((meth, tr, n, b, i))
            Random.seed!(s % typemax(Int32))
            θd = ps[rand(Random.Xoshiro(s), 1:length(ps))]
            try
                r = run_be_on_sim_detailed(simulate_be_trial(merge(θd, (tvbio = tr,)); n_per_arm = n))
                if r !== nothing
                    ok[i] = true; pass[i] = r.overall_pass
                end
            catch
            end
        end
        nv = sum(ok)
        push!(res, (meth, tr, n, b, sum(pass), nv, nv > 0 ? 100 * sum(pass) / nv : NaN))
    end
    CSV.write(f * ".tmp", res); mv(f * ".tmp", f; force = true)
    println(rpad(meth, 10), " n=$(lpad(n,2)) T/R=$tr  pass=$(round(mean(res.pass_rate), digits=2))%  ",
            "valid=$(sum(res.n_valid))/$(B*M)  cell $(round(time()-t0, digits=0)) s  total $(round((time()-t_start)/3600, digits=2)) h")
    flush(stdout)
end
println("02 DONE")
