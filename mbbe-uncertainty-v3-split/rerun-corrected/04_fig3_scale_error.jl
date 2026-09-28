using Serialization, DataFrames, Statistics, CairoMakie
P = joinpath(@__DIR__, "output") * "/"   # diag_trials.jls from 00_scale_error_diagnostic.jl
r = deserialize(P * "diag_trials.jls"); r = r[r.valid, :]
fig = Figure(size=(1400, 540), fontsize=17)
ax1 = Axis(fig[1,1], xlabel="Drawn volume of distribution V (L, log scale)", ylabel="Number of trials",
           title="A. Drawn V in the trials at T/R = 0.70", xticks=(log10.([1,10,50,1000,10000]),["1","10","50","1,000","10,000"]))
for (m,c,l) in (("VCOV_asrun",(:crimson,0.55),"Original sampler"),("VCOV_corrected",(:steelblue,0.7),"Corrected sampler"))
    v = r[(r.method.==m) .& (r.tr.==0.7), :tvv]
    hist!(ax1, log10.(v), bins=range(-1.5, 4.5, length=61), color=c, label=l)
end
vlines!(ax1, [log10(50)], color=:black, linestyle=:dash); axislegend(ax1, position=:rt)
a = r[(r.method.=="VCOV_asrun") .& (r.tr.==0.7), :]
bins = [(0,20,"<20"),(20,100,"20–100"),(100,500,"100–500"),(500,Inf,"≥500")]
ax2 = Axis(fig[1,2], xlabel="Drawn V (L)", ylabel="Concluded BE at T/R = 0.70 (%)", title="B. Original sampler: false acceptance by drawn V",
           xticks=(1:4, [b[3] for b in bins]))
for (k,(n,c,off)) in enumerate(((24,:darkorange,-0.18),(80,:purple4,0.18)))
    y = [100*mean(a[(a.n.==n) .& (a.tvv .>= lo) .& (a.tvv .< hi), :pass]) for (lo,hi,_) in bins]
    barplot!(ax2, (1:4) .+ off, y, width=0.34, color=c, label="n = $n per arm")
end
axislegend(ax2, position=:lt)
save(P * "figures/fig3_scale_error.png", fig, px_per_unit=2); println("ok")
