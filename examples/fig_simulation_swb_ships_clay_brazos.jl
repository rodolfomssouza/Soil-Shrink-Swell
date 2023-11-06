"""
Script to plot soil water balance results
"""

# %% Packages -----------------------------------------------------------------
using DataFrames
using KernelDensity
using Parquet
using Plots
using StatsBase
using StatsPlots


# %% Load data ----------------------------------------------------------------
# Paths
path = "examples/outputs/"

df_shrink = DataFrame(read_parquet(path * "Simulation_swb_shrink_soil.parquet"))
df_rigid = DataFrame(read_parquet(path * "Simulation_swb_rigid_soil.parquet"))


# %% Figures ------------------------------------------------------------------

# Figure size
w = 1200
h = w / 1.618

gr(
    size=(w, h),
    ytickfontsize=10,
    xtickfontsize=10,
    legendfontsize=10,
    xguidefontsize=10,
    yguidefontsize=10,
    bottom_margin=4Plots.mm,
    left_margin=4Plots.mm,
)


# Soil moisture
p1 = plot(
    df_shrink.Days,
    df_shrink.θf,
    ylims=(0.1, 0.5001),
    ylabel="Soil moisture",
    label="Shrink-swell soil",
)

p1 = plot!(
    df_rigid.Days,
    df_rigid.θf,
    label="Rigid soil",
)

# Crack area density
p2 = plot(
    df_shrink.Days,
    df_shrink.CrackDensity,
    ylims=(0.1, 0.5001),
    ylabel="Crack area density",
    label="Shrink-swell soil",
    legend=false
)

# Evapotranspiration
p3 = plot(
    df_shrink.Days,
    df_shrink.ETf,
    ylims=(0.0, 0.5001),
    xlabel="Time (days)",
    ylabel="Evapotranspiration",
    label="Shrink-swell soil",
    legend=false
)

p3 = plot!(
    df_rigid.Days,
    df_rigid.ET,
    label="Rigid soil",
)

# Leakage
p4 = plot(
    df_shrink.Days,
    df_shrink.Lkf,
    xlabel="Time (days)",
    ylabel="Leakage",
    label="Shrink-swell soil",
    legend=false
)

p4 = plot!(
    df_rigid.Days,
    df_rigid.Lk,
    label="Rigid soil",
)

# Render figure
plot(p1, p2, p3, p4, layout=(2, 2))

savefig(path * "Fig_water_balance.pdf")


# %% Soil moisture distribution -----------------------------------------------

θf_shrink = sort(df_shrink.θf)
θf_rigid = sort(df_rigid.θf)

d_shrink = kde(Float64.(θf_shrink))
d_rigid = kde(Float64.(θf_rigid))


w = 800
h = w / 1.618

gr(
    size=(w, h),
    ytickfontsize=10,
    xtickfontsize=10,
    legendfontsize=10,
    xguidefontsize=10,
    yguidefontsize=10,
    bottom_margin=4Plots.mm,
    left_margin=4Plots.mm,
)

p5 = plot(
    θf_shrink,
    pdf(d_shrink, θf_shrink),
    linewidth=2,
    xlabel="Soil moisture",
    ylabel="Probability density",
    label="Shrink-swell soil",
)

p5 = plot!(
    θf_rigid,
    pdf(d_rigid, θf_rigid),
    linewidth=2,
    label="Rigid soil",
)

savefig(path * "Fig_soil_moisture_probability_density.pdf")


# %% Total fluxes -------------------------------------------------------------
# Total fluxes
et = [sum(df_shrink.ETf) / sum(df_shrink.Rain) sum(df_rigid.ET) / sum(df_rigid.Rain)] .* 100
lk = [sum(df_shrink.Lkf) / sum(df_shrink.Rain) sum(df_rigid.Lk) / sum(df_rigid.Rain)] .* 100
q = [sum(df_shrink.Q) / sum(df_shrink.Rain) sum(df_rigid.Q) / sum(df_rigid.Rain)] .* 100

yf = [et lk q]
gf = ["1SS", "2RS", "1SS", "2RS", "1SS", "2RS"]
xf = ["ET", "ET", "L", "L", "Q", "Q"]

p6 = groupedbar(
    xf,
    yf,
    group=gf,
    label=["Shrink-swell soil" "Rigid soil"],
    ylabel="%Rain",
    xlabel="Water fluxes",
    ylims=(0, 110),
)

savefig(path * "Fig_water_fluxes_rainfall.pdf")
