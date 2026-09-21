# MBBE Uncertainty Propagation - Batch Job for JuliaHub
# Extracted from mbbe-uncertainity-v2.qmd

# ── Load Libraries ──────────────────────────────────────────────────
using Pumas
using Random
using Distributions
using LinearAlgebra
using Statistics
using DataFrames
using Bioequivalence
using CairoMakie
using Makie
using CSV

# ── Set Random Seed ─────────────────────────────────────────────────
Random.seed!(142)

# ── PK Model ────────────────────────────────────────────────────────
pkmodel = @model begin
    @param begin
        tvcl  ∈ RealDomain(lower=0)
        tvv   ∈ RealDomain(lower=0)
        tvka  ∈ RealDomain(lower=0)
        tvbio ∈ RealDomain(lower=0)
        Ω     ∈ PDiagDomain(3)
        σ     ∈ RealDomain(lower=0)
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates TRT

    @pre begin
        CL  = tvcl * exp(η[1])
        V   = tvv  * exp(η[2])
        Ka  = tvka * exp(η[3])
    end

    @dosecontrol begin
        bioav = (Depot = TRT == "T" ? tvbio : 1.0,)
    end

    @dynamics begin
        Depot'   = -Ka * Depot
        Central' =  Ka * Depot - (CL / V) * Central
    end

    @derived begin
        cp = @. Central / V
        dv ~ @. Normal(cp, σ)
    end
end

# ── True Parameter Values ───────────────────────────────────────────
θ_true = (
    tvcl  = 5.0,
    tvv   = 50.0,
    tvka  = 1.2,
    tvbio = 1.0,
    Ω     = Diagonal([0.09, 0.09, 0.09]),
    σ     = 0.2
)

# ── Simulation Helpers ──────────────────────────────────────────────
function create_be_population(n_per_arm)
    subjects = Subject[]
    for trt in ["R", "T"], i in 1:n_per_arm
        push!(subjects, Subject(
            id     = length(subjects) + 1,
            events = DosageRegimen(100, time = 0),
            covariates = (TRT = trt,)
        ))
    end
    return Population(subjects)
end

function simulate_be_trial(params; n_per_arm = 30)
    pop = create_be_population(n_per_arm)
    return simobs(pkmodel, pop, params; obstimes = 0:1:72)
end

# ── Simulate Reference Trial and Fit ───────────────────────────────
println("Simulating reference trial...")
ref_sim = simulate_be_trial(θ_true; n_per_arm = 30)

sim_pop = read_pumas(
    DataFrame(ref_sim);
    id = :id, time = :time, observations = [:dv],
    amt = :amt, evid = :evid, covariates = [:TRT], cmt = :cmt
)

init = (
    tvcl  = 4.5,
    tvv   = 45.0,
    tvka  = 1.0,
    tvbio = 1.0,
    Ω     = Diagonal([0.1, 0.1, 0.1]),
    σ     = 0.3
)

println("Fitting model...")
fit_res = fit(pkmodel, sim_pop, init, Pumas.FOCE())
θ̂ = coef(fit_res)

# ── Single-Trial BE (Reference Check) ──────────────────────────────
simdf_ref = DataFrame(ref_sim)
simdf_ref.route .= "ev"
simdf_ref.dv .= max.(simdf_ref.dv, 0.0)

ncapop_ref = read_nca(simdf_ref;
    id = :id, time = :time, amt = :amt,
    observations = :dv, route = :route, group = [:TRT]
)
nca_ref     = run_nca(ncapop_ref)
nca_tbl_ref = nca_ref.reportdf

println("NCA report columns: ", names(nca_tbl_ref))

# ── Uncertainty Extraction ──────────────────────────────────────────
function sample_from_vcov(mean_params, V, n)
    μ = [
        log(mean_params.tvcl), log(mean_params.tvv),
        log(mean_params.tvka), log(mean_params.tvbio),
        log(mean_params.Ω[1,1]), log(mean_params.Ω[2,2]), log(mean_params.Ω[3,3]),
        log(mean_params.σ)
    ]
    dist = MvNormal(μ, Symmetric(V))
    return [begin
        s = rand(dist)
        (tvcl  = exp(s[1]), tvv  = exp(s[2]),
         tvka  = exp(s[3]), tvbio = exp(s[4]),
         Ω     = Diagonal([exp(s[5]), exp(s[6]), exp(s[7])]),
         σ     = exp(s[8]))
    end for _ in 1:n]
end

# Method 1: VCOV
println("Sampling from VCOV...")
vcov_param_samples = sample_from_vcov(θ̂, vcov(fit_res), 20)

# Method 2: Bootstrap
println("Running Bootstrap inference (100 samples)...")
boot_settings = Bootstrap(
    samples = 100,
    stratify_by = nothing,
    ensemblealg = EnsembleThreads()
)
boot_inf = infer(fit_res, boot_settings; level = 0.95)
boot_param_samples = sample_from_vcov(coef(boot_inf), vcov(boot_inf), 20)

# Method 3: SIR
println("Running SIR inference...")
sir_settings = SIR(samples = 200, resamples = 100)
sir_inf = infer(fit_res, sir_settings; level = 0.95)
sir_param_samples = sample_from_vcov(coef(sir_inf), vcov(sir_inf), 20)

# ── TOST Helpers ────────────────────────────────────────────────────
function tost_test_detailed(log_R, log_T; alpha = 0.05, θ_L = 0.80, θ_U = 1.25)
    n_R, n_T = length(log_R), length(log_T)
    (n_R < 2 || n_T < 2) && return (gmr=NaN, ci_lo=NaN, ci_hi=NaN, pass=false)
    Δ   = mean(log_T) - mean(log_R)
    sp² = ((n_R - 1) * var(log_R) + (n_T - 1) * var(log_T)) / (n_R + n_T - 2)
    se  = sqrt(sp² * (1 / n_R + 1 / n_T))
    df  = n_R + n_T - 2
    tc  = quantile(TDist(df), 1 - alpha)
    gmr   = exp(Δ)
    ci_lo = exp(Δ - tc * se)
    ci_hi = exp(Δ + tc * se)
    pass  = ci_lo >= θ_L && ci_hi <= θ_U
    return (gmr=gmr, ci_lo=ci_lo, ci_hi=ci_hi, pass=pass)
end

function find_col(df, candidates)
    ns = names(df)
    for c in candidates
        if string(c) in ns
            return Symbol(c)
        end
    end
    error("None of $candidates found in columns: $ns")
end

const AUC_COL  = find_col(nca_tbl_ref, ["aucinf_obs", "aucinf", "auclast", "AUC", "auc"])
const CMAX_COL = find_col(nca_tbl_ref, ["cmax", "Cmax", "CMAX"])
println("Using NCA columns: AUC = $AUC_COL, Cmax = $CMAX_COL")

function clean_pk_col(df, col)
    vals = df[!, col]
    return Float64[v for v in vals if !ismissing(v) && v > 0]
end

function run_be_on_sim_detailed(sim)
    simdf = DataFrame(sim)
    simdf.route .= "ev"
    simdf.dv .= max.(simdf.dv, 0.0)

    ncapop = read_nca(simdf;
        id = :id, time = :time, amt = :amt,
        observations = :dv, route = :route, group = [:TRT]
    )
    nca_out = run_nca(ncapop)
    tbl     = nca_out.reportdf

    R = tbl[tbl.TRT .== "R", :]
    T = tbl[tbl.TRT .== "T", :]

    auc_R  = clean_pk_col(R, AUC_COL)
    auc_T  = clean_pk_col(T, AUC_COL)
    cmax_R = clean_pk_col(R, CMAX_COL)
    cmax_T = clean_pk_col(T, CMAX_COL)

    (length(auc_R) < 2 || length(auc_T) < 2)  && return nothing
    (length(cmax_R) < 2 || length(cmax_T) < 2) && return nothing

    auc_res  = tost_test_detailed(log.(auc_R),  log.(auc_T))
    cmax_res = tost_test_detailed(log.(cmax_R), log.(cmax_T))

    return (
        auc_gmr=auc_res.gmr, auc_lo=auc_res.ci_lo, auc_hi=auc_res.ci_hi, auc_pass=auc_res.pass,
        cmax_gmr=cmax_res.gmr, cmax_lo=cmax_res.ci_lo, cmax_hi=cmax_res.ci_hi, cmax_pass=cmax_res.pass,
        overall_pass=auc_res.pass && cmax_res.pass
    )
end

# ── %CV Estimation Helper ──────────────────────────────────────────
function estimate_cv(param_samples; n_per_arm=60, n_reps=50)
    auc_vars  = Float64[]
    cmax_vars = Float64[]

    for _ in 1:n_reps
        θ_draw = param_samples[rand(1:length(param_samples))]
        sim_params = merge(θ_draw, (tvbio = 1.0,))

        try
            sim = simulate_be_trial(sim_params; n_per_arm = n_per_arm)
            simdf = DataFrame(sim)
            simdf.route .= "ev"
            simdf.dv .= max.(simdf.dv, 0.0)

            ncapop = read_nca(simdf;
                id = :id, time = :time, amt = :amt,
                observations = :dv, route = :route, group = [:TRT]
            )
            nca_out = run_nca(ncapop)
            tbl = nca_out.reportdf

            R = tbl[tbl.TRT .== "R", :]
            T = tbl[tbl.TRT .== "T", :]

            auc_R  = clean_pk_col(R, AUC_COL)
            auc_T  = clean_pk_col(T, AUC_COL)
            cmax_R = clean_pk_col(R, CMAX_COL)
            cmax_T = clean_pk_col(T, CMAX_COL)

            if length(auc_R) >= 2 && length(auc_T) >= 2
                n_r, n_t = length(auc_R), length(auc_T)
                sp² = ((n_r - 1) * var(log.(auc_R)) + (n_t - 1) * var(log.(auc_T))) / (n_r + n_t - 2)
                push!(auc_vars, sp²)
            end
            if length(cmax_R) >= 2 && length(cmax_T) >= 2
                n_r, n_t = length(cmax_R), length(cmax_T)
                sp² = ((n_r - 1) * var(log.(cmax_R)) + (n_t - 1) * var(log.(cmax_T))) / (n_r + n_t - 2)
                push!(cmax_vars, sp²)
            end
        catch
            continue
        end
    end

    auc_cv  = 100 * sqrt(exp(mean(auc_vars)) - 1)
    cmax_cv = 100 * sqrt(exp(mean(cmax_vars)) - 1)
    return (auc_cv=auc_cv, cmax_cv=cmax_cv)
end

# ── Virtual BE Simulation Engine ────────────────────────────────────
const SAMPLE_SIZES = [12, 24, 30, 40, 50, 60, 80]
const TR_RATIOS    = [0.8, 0.9, 1.0, 1.1, 1.2]

function virtual_be_sim(param_samples, method_name;
                        tr_ratios  = TR_RATIOS,
                        n_per_arms = SAMPLE_SIZES,
                        n_sims     = 20)

    summary_results = DataFrame(
        method    = String[],
        tr_ratio  = Float64[],
        n_per_arm = Int[],
        n_pass    = Int[],
        n_valid   = Int[],
        n_error   = Int[],
        pass_rate = Float64[]
    )

    detailed_results = DataFrame(
        method       = String[],
        tr_ratio     = Float64[],
        n_per_arm    = Int[],
        sim_id       = Int[],
        auc_gmr      = Float64[],
        auc_lo       = Float64[],
        auc_hi       = Float64[],
        auc_pass     = Bool[],
        cmax_gmr     = Float64[],
        cmax_lo      = Float64[],
        cmax_hi      = Float64[],
        cmax_pass    = Bool[],
        overall_pass = Bool[]
    )

    total_combos = length(n_per_arms) * length(tr_ratios)
    combo = 0

    for n_arm in n_per_arms, tr in tr_ratios
        combo += 1
        n_pass  = 0
        n_error = 0

        for i in 1:n_sims
            θ_draw     = param_samples[rand(1:length(param_samples))]
            sim_params = merge(θ_draw, (tvbio = tr,))

            try
                sim = simulate_be_trial(sim_params; n_per_arm = n_arm)
                res = run_be_on_sim_detailed(sim)
                if res === nothing
                    n_error += 1
                else
                    res.overall_pass && (n_pass += 1)
                    push!(detailed_results, (
                        method       = method_name,
                        tr_ratio     = tr,
                        n_per_arm    = n_arm,
                        sim_id       = i,
                        auc_gmr      = res.auc_gmr,
                        auc_lo       = res.auc_lo,
                        auc_hi       = res.auc_hi,
                        auc_pass     = res.auc_pass,
                        cmax_gmr     = res.cmax_gmr,
                        cmax_lo      = res.cmax_lo,
                        cmax_hi      = res.cmax_hi,
                        cmax_pass    = res.cmax_pass,
                        overall_pass = res.overall_pass
                    ))
                end
            catch e
                n_error += 1
                if n_error == 1
                    @warn "virtual_be_sim error" method=method_name tr_ratio=tr n_per_arm=n_arm exception=(e, catch_backtrace())
                end
            end
        end

        n_valid = n_sims - n_error
        rate    = n_valid > 0 ? (n_pass / n_valid) * 100 : NaN

        push!(summary_results, (
            method    = method_name,
            tr_ratio  = tr,
            n_per_arm = n_arm,
            n_pass    = n_pass,
            n_valid   = n_valid,
            n_error   = n_error,
            pass_rate = rate
        ))

        println("[$method_name] ($combo/$total_combos) n=$n_arm, T/R=$tr => pass rate = $(round(rate, digits=1))%")
    end

    return (summary = summary_results, detailed = detailed_results)
end

# ── Run All Methods ─────────────────────────────────────────────────
println("\n" * "="^60)
println("Running Virtual BE Simulations")
println("="^60)

fixed_param_samples = fill(θ̂, 100)

println("\n--- Fixed ---")
res_fixed = virtual_be_sim(fixed_param_samples, "Fixed")
println("\n--- VCOV ---")
res_vcov  = virtual_be_sim(vcov_param_samples,  "VCOV")
println("\n--- Bootstrap ---")
res_boot  = virtual_be_sim(boot_param_samples,  "Bootstrap")
println("\n--- SIR ---")
res_sir   = virtual_be_sim(sir_param_samples,   "SIR")

# ── %CV Estimation ─────────────────────────────────────────────────
cv_estimates = estimate_cv(fixed_param_samples; n_per_arm = 60, n_reps = 100)

println("\n" * "=" ^ 55)
println("  Olanzapine %CV Estimates (Parallel Design)")
println("=" ^ 55)
println()
println("  NCA-derived (pooled from 100 simulated trials):")
println("    AUC  %CV = $(round(cv_estimates.auc_cv, digits=1))%")
println("    Cmax %CV = $(round(cv_estimates.cmax_cv, digits=1))%")
println()
println("  Model-derived IIV (between-subject) %CV:")
for (name, idx) in [("CL", 1), ("V", 2), ("Ka", 3)]
    ω² = θ̂.Ω[idx, idx]
    cv = 100 * sqrt(exp(ω²) - 1)
    println("    $name:  $(round(cv, digits=1))%  (omega² = $(round(ω², digits=4)))")
end
println("    Residual σ (additive) = $(round(θ̂.σ, digits=4))")
println("=" ^ 55)

# ── Combine Results ─────────────────────────────────────────────────
all_summary  = vcat(res_fixed.summary, res_vcov.summary, res_boot.summary, res_sir.summary)
all_detailed = vcat(res_fixed.detailed, res_vcov.detailed, res_boot.detailed, res_sir.detailed)

# ── Save CSV Results ────────────────────────────────────────────────
outdir = get(ENV, "RESULTS_DIR", "results")
mkpath(outdir)

CSV.write(joinpath(outdir, "mbbe_summary.csv"), all_summary)
CSV.write(joinpath(outdir, "mbbe_detailed.csv"), all_detailed)
println("\nCSV results saved to $outdir/")

# ── Print Summary Tables ───────────────────────────────────────────
for n in [24, 30, 48, 60]
    println("\n=== n per arm = $n ===")
    sub = filter(r -> r.n_per_arm == n, all_summary)
    wide = unstack(
        select(sub, :tr_ratio, :method, :pass_rate),
        :tr_ratio, :method, :pass_rate
    )
    display(wide)
end

# ── Power & Type I Error ───────────────────────────────────────────
power_df = filter(r -> r.tr_ratio in [0.9, 1.0, 1.1], all_summary)
power_summary = combine(
    groupby(power_df, [:method, :n_per_arm]),
    :pass_rate => mean => :mean_power
)

t1e_df = filter(r -> r.tr_ratio in [0.8, 1.2], all_summary)
t1e_summary = combine(
    groupby(t1e_df, [:method, :n_per_arm]),
    :pass_rate => mean => :mean_type1_error
)

println("\n--- Power (T/R in {0.9, 1.0, 1.1}) ---")
display(unstack(power_summary, :n_per_arm, :method, :mean_power))

println("\n--- Type I Error (T/R in {0.8, 1.2}) ---")
display(unstack(t1e_summary, :n_per_arm, :method, :mean_type1_error))

CSV.write(joinpath(outdir, "mbbe_power_summary.csv"), power_summary)
CSV.write(joinpath(outdir, "mbbe_type1error_summary.csv"), t1e_summary)

# ── Generate and Save Figures ───────────────────────────────────────
println("\nGenerating figures...")
methods_list = ["Fixed", "VCOV", "Bootstrap", "SIR"]
size_colors  = [:indigo, :dodgerblue, :steelblue, :seagreen, :orange, :crimson, :darkred]
method_colors = [:steelblue, :crimson, :seagreen, :darkorange]
method_shapes = [:circle, :diamond, :utriangle, :rect]

# Fig 1: Pass Rate by T/R Ratio and Sample Size
fig1 = Figure(size = (1400, 900))
for (idx, m) in enumerate(methods_list)
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1
    ax = Axis(fig1[row, col];
        xlabel = "T/R Ratio", ylabel = "BE Pass Rate (%)",
        title  = "$m Method", xticks = TR_RATIOS, yticks = 0:10:100
    )
    sub = filter(r -> r.method == m, all_summary)
    for (j, n) in enumerate(SAMPLE_SIZES)
        sub_n = sort(filter(r -> r.n_per_arm == n, sub), :tr_ratio)
        scatterlines!(ax, sub_n.tr_ratio, sub_n.pass_rate;
            label = "n=$n", color = size_colors[j], marker = :circle, markersize = 8, linewidth = 2)
    end
    hlines!(ax, [5.0]; color = :red, linestyle = :dash, linewidth = 1, label = "5% α")
    vlines!(ax, [0.80, 1.25]; color = :gray70, linestyle = :dot, linewidth = 1)
    idx == 2 && axislegend(ax; position = :rt, labelsize = 9)
end
Label(fig1[0, :], "BE Pass Rate by T/R Ratio and Sample Size"; fontsize = 20)
save(joinpath(outdir, "fig1_passrate_lines.png"), fig1)

# Fig 2: Heatmaps
fig2 = Figure(size = (1400, 900))
for (idx, m) in enumerate(methods_list)
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1
    sub = filter(r -> r.method == m, all_summary)
    mat = zeros(length(SAMPLE_SIZES), length(TR_RATIOS))
    for (i, n) in enumerate(SAMPLE_SIZES), (j, tr) in enumerate(TR_RATIOS)
        r = filter(r -> r.n_per_arm == n && r.tr_ratio == tr, sub)
        mat[i, j] = nrow(r) > 0 ? r.pass_rate[1] : NaN
    end
    ax = Axis(fig2[row, col];
        xlabel = "T/R Ratio", ylabel = "Subjects per Arm", title = "$m Method",
        xticks = (1:length(TR_RATIOS), string.(TR_RATIOS)),
        yticks = (1:length(SAMPLE_SIZES), string.(SAMPLE_SIZES))
    )
    hm = heatmap!(ax, 1:length(TR_RATIOS), 1:length(SAMPLE_SIZES), mat';
        colormap = :RdYlGn, colorrange = (0, 100))
    for (i, n) in enumerate(SAMPLE_SIZES), (j, tr) in enumerate(TR_RATIOS)
        txt = isnan(mat[i, j]) ? "—" : "$(round(Int, mat[i, j]))%"
        text!(ax, j, i; text = txt, align = (:center, :center),
            fontsize = 11, color = mat[i, j] > 50 ? :black : :white)
    end
    idx == 2 && Colorbar(fig2[row, col + 1], hm; label = "Pass Rate (%)")
end
Label(fig2[0, :], "BE Pass Rate Heatmap: Sample Size vs T/R Ratio"; fontsize = 20)
save(joinpath(outdir, "fig2_heatmaps.png"), fig2)

# Fig 3: Power Curves
fig3 = Figure(size = (1000, 600))
ax3 = Axis(fig3[1, 1];
    xlabel = "Subjects per Arm", ylabel = "Mean Power (%)",
    title = "Power (T/R in {0.9, 1.0, 1.1}) by Sample Size", xticks = SAMPLE_SIZES
)
for (i, m) in enumerate(methods_list)
    sub = sort(filter(r -> r.method == m, power_summary), :n_per_arm)
    scatterlines!(ax3, sub.n_per_arm, sub.mean_power;
        label = m, color = method_colors[i], marker = method_shapes[i], markersize = 12, linewidth = 2)
end
hlines!(ax3, [80.0]; color = :gray, linestyle = :dash, linewidth = 1, label = "80% power")
axislegend(ax3; position = :rb)
save(joinpath(outdir, "fig3_power_curves.png"), fig3)

# Fig 4: AUC GMR Histograms
fig4 = Figure(size = (1400, 500))
detail_sub = filter(r -> r.method == "Fixed" && r.n_per_arm == 30, all_detailed)
for (j, tr) in enumerate(TR_RATIOS)
    sub_tr = filter(r -> r.tr_ratio == tr, detail_sub)
    nrow(sub_tr) == 0 && continue
    ax = Axis(fig4[1, j]; xlabel = "AUC GMR", ylabel = j == 1 ? "Count" : "", title = "T/R = $tr")
    hist!(ax, sub_tr.auc_gmr; bins = 25, color = (:steelblue, 0.7), strokecolor = :black, strokewidth = 0.5)
    vlines!(ax, [0.80, 1.25]; color = :red, linestyle = :dash, linewidth = 2)
    vlines!(ax, [tr]; color = :black, linestyle = :solid, linewidth = 1.5)
end
Label(fig4[0, :], "AUC Geometric Mean Ratio Distribution (Fixed, n=30)"; fontsize = 18)
save(joinpath(outdir, "fig4_auc_gmr_hist.png"), fig4)

# Fig 5: Cmax GMR Histograms
fig5 = Figure(size = (1400, 500))
for (j, tr) in enumerate(TR_RATIOS)
    sub_tr = filter(r -> r.tr_ratio == tr, detail_sub)
    nrow(sub_tr) == 0 && continue
    ax = Axis(fig5[1, j]; xlabel = "Cmax GMR", ylabel = j == 1 ? "Count" : "", title = "T/R = $tr")
    hist!(ax, sub_tr.cmax_gmr; bins = 25, color = (:seagreen, 0.7), strokecolor = :black, strokewidth = 0.5)
    vlines!(ax, [0.80, 1.25]; color = :red, linestyle = :dash, linewidth = 2)
    vlines!(ax, [tr]; color = :black, linestyle = :solid, linewidth = 1.5)
end
Label(fig5[0, :], "Cmax Geometric Mean Ratio Distribution (Fixed, n=30)"; fontsize = 18)
save(joinpath(outdir, "fig5_cmax_gmr_hist.png"), fig5)

# Fig 6: CI Bounds
fig6 = Figure(size = (1400, 800))
for (j, tr) in enumerate(TR_RATIOS)
    sub_tr = filter(r -> r.tr_ratio == tr, detail_sub)
    nrow(sub_tr) == 0 && continue
    ax_lb = Axis(fig6[1, j]; xlabel = j == 3 ? "AUC 90% CI Lower Bound" : "", ylabel = j == 1 ? "Count" : "", title = "T/R = $tr")
    hist!(ax_lb, sub_tr.auc_lo; bins = 25, color = (:coral, 0.7), strokecolor = :black, strokewidth = 0.5)
    vlines!(ax_lb, [0.80]; color = :red, linestyle = :dash, linewidth = 2)
    ax_ub = Axis(fig6[2, j]; xlabel = j == 3 ? "AUC 90% CI Upper Bound" : "", ylabel = j == 1 ? "Count" : "")
    hist!(ax_ub, sub_tr.auc_hi; bins = 25, color = (:mediumpurple, 0.7), strokecolor = :black, strokewidth = 0.5)
    vlines!(ax_ub, [1.25]; color = :red, linestyle = :dash, linewidth = 2)
end
Label(fig6[0, :], "AUC 90% CI Bounds Distribution (Fixed, n=30)"; fontsize = 18)
Label(fig6[1, 0], "Lower Bound"; fontsize = 14, rotation = pi/2)
Label(fig6[2, 0], "Upper Bound"; fontsize = 14, rotation = pi/2)
save(joinpath(outdir, "fig6_ci_bounds.png"), fig6)

# Fig 7: CI Width vs Sample Size
fig7 = Figure(size = (1000, 600))
detail_tr1 = filter(r -> r.method == "Fixed" && r.tr_ratio == 1.0, all_detailed)
ci_width_summary = combine(
    groupby(detail_tr1, :n_per_arm),
    :auc_lo  => median => :median_lb,
    :auc_hi  => median => :median_ub,
    :auc_gmr => median => :median_gmr,
    [:auc_hi, :auc_lo] => ((hi, lo) -> median(hi .- lo)) => :median_ci_width
)
sort!(ci_width_summary, :n_per_arm)
ax7 = Axis(fig7[1, 1];
    xlabel = "Subjects per Arm", ylabel = "Value",
    title = "AUC: Median GMR and 90% CI Bounds vs Sample Size (T/R = 1.0, Fixed)", xticks = SAMPLE_SIZES
)
band!(ax7, ci_width_summary.n_per_arm, ci_width_summary.median_lb, ci_width_summary.median_ub;
    color = (:steelblue, 0.3), label = "Median 90% CI")
scatterlines!(ax7, ci_width_summary.n_per_arm, ci_width_summary.median_gmr;
    color = :steelblue, marker = :circle, markersize = 10, linewidth = 2, label = "Median GMR")
scatterlines!(ax7, ci_width_summary.n_per_arm, ci_width_summary.median_lb;
    color = :coral, marker = :dtriangle, markersize = 8, linewidth = 1.5, linestyle = :dash, label = "Median LB")
scatterlines!(ax7, ci_width_summary.n_per_arm, ci_width_summary.median_ub;
    color = :mediumpurple, marker = :utriangle, markersize = 8, linewidth = 1.5, linestyle = :dash, label = "Median UB")
hlines!(ax7, [0.80, 1.25]; color = :red, linestyle = :dot, linewidth = 1, label = "BE limits")
axislegend(ax7; position = :rt, labelsize = 10)
save(joinpath(outdir, "fig7_ci_width.png"), fig7)

# Fig 8: Method Comparison
fig8 = Figure(size = (900, 500))
ax8 = Axis(fig8[1, 1];
    xlabel = "T/R Ratio", ylabel = "BE Pass Rate (%)",
    title = "Impact of Parameter Uncertainty on BE Decisions (n = 30)", xticks = TR_RATIOS
)
for (i, m) in enumerate(methods_list)
    sub = sort(filter(r -> r.method == m && r.n_per_arm == 30, all_summary), :tr_ratio)
    scatterlines!(ax8, sub.tr_ratio, sub.pass_rate;
        label = m, color = method_colors[i], marker = method_shapes[i], markersize = 12, linewidth = 2)
end
hlines!(ax8, [5.0]; color = :gray, linestyle = :dash, linewidth = 1, label = "5% α")
vlines!(ax8, [0.80, 1.25]; color = :gray70, linestyle = :dot, linewidth = 1, label = "BE limits")
axislegend(ax8; position = :rt)
save(joinpath(outdir, "fig8_method_comparison.png"), fig8)

println("\nAll figures saved to $outdir/")
println("\n" * "="^60)
println("MBBE batch job complete!")
println("="^60)
