# 03 — summarise finished cells into the paper's tables and figures. Safe to run on a partial grid.
using CSV, DataFrames, Statistics, CairoMakie
const OUT  = joinpath(@__DIR__, "output")
const FIG  = normpath(joinpath(@__DIR__, "..", "figures"))
const METH = ["Fixed", "VCOV", "Bootstrap", "SIR"]
const COLS = Dict("Fixed" => :black, "VCOV" => :crimson, "Bootstrap" => :seagreen, "SIR" => :darkorange)
const MRK  = Dict("Fixed" => :circle, "VCOV" => :diamond, "Bootstrap" => :utriangle, "SIR" => :rect)

files = filter(f -> endswith(f, ".csv"), readdir(joinpath(OUT, "cells"); join = true))
batches = reduce(vcat, [CSV.read(f, DataFrame) for f in files])
S = combine(groupby(batches, [:method, :tr_ratio, :n_per_arm]),
    :pass_rate => mean => :rate, :pass_rate => std => :sd, :pass_rate => length => :B,
    :n_valid => sum => :valid, :n_pass => sum => :passed)
S.se = S.sd ./ sqrt.(S.B)
S.lo = max.(S.rate .- 1.96 .* S.se, 0.0); S.hi = min.(S.rate .+ 1.96 .* S.se, 100.0)
sort!(S, [:tr_ratio, :n_per_arm, :method])
CSV.write(joinpath(OUT, "operating_characteristics_corrected.csv"), S)

# differences from Fixed (pp) with 95% CI
F = select(filter(:method => ==("Fixed"), S), :tr_ratio, :n_per_arm, :rate => :f_rate, :se => :f_se)
D = innerjoin(filter(:method => !=("Fixed"), S), F, on = [:tr_ratio, :n_per_arm])
D.diff = D.rate .- D.f_rate
D.dse  = sqrt.(D.se .^ 2 .+ D.f_se .^ 2)
D.dlo  = D.diff .- 1.96 .* D.dse; D.dhi = D.diff .+ 1.96 .* D.dse
CSV.write(joinpath(OUT, "difference_from_fixed.csv"), select(D, :method, :tr_ratio, :n_per_arm, :diff, :dlo, :dhi))

# sample size for 80% power at T/R = 1.0 (linear interpolation)
N80 = DataFrame(method = String[], n80 = Union{Missing,Float64}[])
for m in METH
    s = sort(filter(r -> r.method == m && r.tr_ratio == 1.0, S), :n_per_arm)
    k = findfirst(>=(80.0), s.rate)
    v = (k === nothing || k == 1) ? missing :
        s.n_per_arm[k-1] + (80 - s.rate[k-1]) * (s.n_per_arm[k] - s.n_per_arm[k-1]) / (s.rate[k] - s.rate[k-1])
    push!(N80, (m, v))
end
CSV.write(joinpath(OUT, "n_for_80pct_power.csv"), N80)

# markdown tables for the manuscript
fmt(r) = "$(round(r.rate, digits=1)) ($(round(r.lo, digits=1))–$(round(r.hi, digits=1)))"
open(joinpath(OUT, "tables.md"), "w") do io
    for (title, trs) in (("Power (%) — T/R within limits", [0.9, 1.0, 1.1]),
                         ("Probability (%) of concluding BE — T/R outside limits", [0.7, 1.3]))
        println(io, "### $title\n\n| T/R | n per arm | ", join(METH, " | "), " |\n|---|---|", repeat("---|", length(METH)))
        for tr in trs, n in sort(unique(S.n_per_arm))
            cells = [let r = filter(x -> x.method == m && x.tr_ratio == tr && x.n_per_arm == n, S)
                         nrow(r) == 1 ? fmt(r[1, :]) : "—" end for m in METH]
            println(io, "| $tr | $n | ", join(cells, " | "), " |")
        end
        println(io)
    end
    println(io, "### n per arm for 80% power at T/R = 1.00\n\n| Method | n |\n|---|---|")
    for r in eachrow(N80); println(io, "| $(r.method) | $(r.n80 === missing ? "—" : round(r.n80, digits=1)) |"); end
end

# figures
function panel!(ax, tr)
    for m in METH
        s = sort(filter(r -> r.method == m && r.tr_ratio == tr, S), :n_per_arm)
        nrow(s) == 0 && continue
        off = (findfirst(==(m), METH) - 2.5) * 0.6
        rangebars!(ax, s.n_per_arm .+ off, s.lo, s.hi, color = COLS[m], whiskerwidth = 6)
        scatterlines!(ax, s.n_per_arm .+ off, s.rate, color = COLS[m], marker = MRK[m], markersize = 10, label = m)
    end
end
fig = Figure(size = (1500, 520), fontsize = 17)
for (j, tr) in enumerate([0.9, 1.0, 1.1])
    ax = Axis(fig[1, j], title = "T/R = $tr", xlabel = "Subjects per arm", ylabel = j == 1 ? "Power (%)" : "",
              xticks = [12, 24, 30, 40, 50, 60, 80], limits = (nothing, (0, 101)))
    hlines!(ax, [80], color = :gray50, linestyle = :dot); panel!(ax, tr)
    j == 3 && try axislegend(ax, position = :rb) catch end
end
save(joinpath(FIG, "fig1_power.png"), fig, px_per_unit = 2)
fig = Figure(size = (1100, 520), fontsize = 17)
for (j, tr) in enumerate([0.7, 1.3])
    ax = Axis(fig[1, j], title = "T/R = $tr", xlabel = "Subjects per arm",
              ylabel = j == 1 ? "Probability of concluding BE (%)" : "", xticks = [12, 24, 30, 40, 50, 60, 80])
    hlines!(ax, [5], color = :gray50, linestyle = :dot); panel!(ax, tr)
    j == 2 && try axislegend(ax, position = :rt) catch end
end
save(joinpath(FIG, "fig2_false_acceptance.png"), fig, px_per_unit = 2)
println("cells summarised: $(nrow(S)) of 140")
