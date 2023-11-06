module SoilShrinkSwell

# %% Packages ---------------------------------------------------------------
using Distributions
using DataFrames
using Roots
using LsqFit
using StatsBase
using LinearAlgebra


# %% Exported functions -----------------------------------------------------
export rainfall_poisson, rainfall_series_poisson
export evapotranspiration, ϕ_agg, porosity_θ
export leakage_ss, swb_ss, sol_swb_ss, pdf_soil_moisture, dt2daily
export pedon_scale_variables


# %% Functions --------------------------------------------------------------

"""
`rainfall_poisson(n, α, λ)`

Generates rainfall using Poisson process.

# Arguments
- `n::Integer`: number of events.
- `α::Float64`: mean rain depth.
- `λ::Float64`: rain event interval (1/d).

If rainfall needs to be similated in timescale smaller than daily, multiply `λ` by dt.

# Example
```
nevents = 30
α = 0.75
λ = 0.45
rain = rainfall_poisson(n, α, λ)
```
"""
function rainfall_poisson(nevents::Int64, α::Float64, λ::Float64)
    rain::Vector{Float64} = zeros(nevents)
    for i = 1:nevents
        r1 = rand(Uniform(0, 1))
        if r1 < λ
            r2 = rand(Uniform(0, 1))
            rain[i] = -α * log((1 - r2))
        else
            rain[i] = 0
        end
    end
    return rain
end


"""
Create rainfall series
"""
function rainfall_series_poisson(α, λ, dt, years)
    nr = years * 365 / dt
    rain = rainfall_poisson(Int(nr), α, λ * dt)

    days = Int.(repeat(1:1:(nr*dt), inner=(Int(1 / dt))))
    years = 1:years
    years = repeat(years, inner=Int(365 / dt))
    df = DataFrame(Years=years, Days=days, Rain=rain)
    return df

end


"""
Computes evapotranspiration losses

# Arguments
- `s`: soil moisture
- `θh`: soil moisture at hygroscopic point
- `θw`: soil moisture at wilting point
- `θstar`: soil moisture below field capacity
- `emax`: maximum evapotranspiration
- `ew`: evaporation

"""
function evapotranspiration(θ, θh, θw, θstar, emax, ew)
    if θ <= θh
        et = 0
    elseif θh <= θ <= θw
        et = ew * (θ - θh) / (θw - θh)
    elseif θw <= θ <= θstar
        et = ew + (emax - ew) * (θ - θw) / (θstar - θw)
    else
        et = emax
    end
    return et
end


"""
Computes soil porosity as a function of soil moisture based on Steawrt et al. (2015).
https://doi.org/10.2136/vzj2015.11.0146

:param s: soil moisture
:param ϕmax: maximum porosity
:param ϕmin: minimum porosity
:param ϵ: fitting parameter
:param q: fitting parameter
:return: soil porosity
"""
function ϕ_agg(θ, ϕmax, ϕmin, ϵ, q)
    Φ = θ / ϕmax
    ϕagg = (ϕmax - ϕmin) * ((ϵ + 1) / (ϵ + Φ^-q)) + ϕmin
    return ϕagg
end


"""
Computes matrix porosity as a function of volumetric soil moisture (θ).
This function is dependent of ``find_zero`` from Roots package.

# Arguments
- `θ::Float64`: volumetric soil moisture
- `ϕmax::Float64`: maximum porosity
- `ϕmin::Float64`: minimum porosity
- `ϵ::Float64`: parameter of form
- `q::Float64`: shape parameter

"""
function porosity_θ(θ, ϕmax, ϕmin, ϵ, q)
    y = find_zero(
        ϕ ->
            ϕ - ((ϕmax - ϕmin) * ((ϵ + 1) / (ϵ + ((θ / ϕmax) * (1 + ϕ - ϕmax))^-q)) + ϕmin),
        (ϕmin, ϕmax),
    )
    return y
end


"""
`leakage_ss(θ, ϕmax, ϕmax, ϕmin, ϵ, q, ksmax, b)`

Salvo's proposal
:param θ: volumetric soil moisture
:param ϕmax: maximum porosity
:param ϕmin: minimum porosity
:param ϵ: fitting parameter
:param q: fitting parameter
:param ksmax: maximum soil hydraulic conductivity
:param b: fitting parameter
:return: leakage
"""
function leakage_ss(θ, θfc, ϕmax, ϕmin, ϵ, q, ks, b)
    # ϕ = ϕ_agg(θ, ϕmax, ϕmin, ϵ, q)
    ϕ = porosity_θ(θ, ϕmax, ϕmin, ϵ, q)
    if θ <= θfc
        y = 0
    else
        y = (ks / ϕmax^2) * ϕ^2 * (θ / ϕ)^((4 - 2 * b) / (1 - b))
    end
    return y
end


"""
`swb_ss(rain, θ, θh, θw, θstar, ϕmax, ϕmin, ksmax, ksmin, ϵ, q, emax, ew, zr, dt)`

Soil water balance in shrink-swell soil.

:param rain: rainfall
:param θ: volumetric soil moisture
:param θh: soil moisture at hygroscopic point
:param θw: soil moisture at wilting point
:param θstar: soil moisture below field capacity
:param ϕmax: maximum porosity
:param ϕmin: minimum porosity
:param ksmax: maximum soil hydraulic conductivity
:param ksmin: minimum soil hydraulic conductivity
:param ϵ: fitting parameter
:param q: fitting parameter
:param emax: maximum evapotranspiration
:param ew: evaporation
:param zr: soil depth
:param dt: time step
:return: dictionary with soil water balance components
"""
function swb_ss(rain, θ, θh, θw, θstar, θfc, ϕmax, ϕmin, ksmax, b, ϵ, q, emax, ew, zr, dt)

    ϕ = porosity_θ(θ, ϕmax, ϕmin, ϵ, q)

    θ = θ + rain / zr

    if θ > ϕmax
        Q = (θ - ϕmax) * zr
        θ = ϕmax
    else
        Q = 0
    end
    lk = leakage_ss(θ, θfc, ϕmax, ϕmin, ϵ, q, ksmax, b) * dt
    θ = θ - lk / zr

    et = evapotranspiration(θ, θh, θw, θstar, emax, ew) * dt
    θ = θ - et / zr

    out = (rain=rain, θ=θ, Q=Q, lk=lk, et=et, ϕ=ϕ)
    return out
end


"""
`sol_swb_ss(rain, θ, params)`

Solve soil water balance in shrink-swell soil.
:param rain: rainfall
:param θ: volumetric soil moisture
:param params: dictionary with parameters
"""
function sol_swb_ss(rain, θ, params)
    # Unzip parameter
    θh = params[:θh]
    θw = params[:θw]
    θstar = params[:θstar]
    θfc = params[:θfc]
    ϕmax = params[:ϕmax]
    ϕmin = params[:ϕmin]
    ksmax = params[:ksmax]
    b = params[:b]
    ϵ = params[:ϵ]
    q = params[:q]
    emax = params[:emax]
    ew = params[:ew]
    zr = params[:zr]
    dt = params[:dt]

    nr = length(rain) + 1
    θout = zeros(nr)
    qout = zeros(nr)
    etout = zeros(nr)
    lkout = zeros(nr)
    nout = zeros(nr)
    θout[1] = θ

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner=Int(1 / dt))
    days = vcat(0, days)
    daysc = [(1+dt):dt:(ndays+1);] .- dt

    for i in eachindex(rain)
        sol = swb_ss(
            rain[i],
            θout[i],
            θh,
            θw,
            θstar,
            θfc,
            ϕmax,
            ϕmin,
            ksmax,
            b,
            ϵ,
            q,
            emax,
            ew,
            zr,
            dt,
        )
        θout[i+1] = sol.θ
        qout[i+1] = sol.Q
        lkout[i+1] = sol.lk
        etout[i+1] = sol.et
        nout[i+1] = sol.ϕ
    end
    df = DataFrame(
        Days=days,
        DaysC=vcat(0, daysc),
        Rain=vcat(0, rain),
        Q=qout,
        θ=θout,
        Lk=lkout,
        ET=etout,
        ϕ=nout,
    )
    return df
end


# Soil moisture pdf
"""
`pdf_soil_moisture(s, sh, sw, sstar, sfc, β, ks, n, zr, ew, emax, α, λ):`

Computes the pdf of soil moisture in steady-state conditions.
https://doi.org/10.1016/S0309-1708(01)00005-7
:param s: soil moisture
:param sh: soil moisture at hygroscopic point
:param sw: soil moisture at wilting point
:param sstar: soil moisture below field capacity
:param sfc: soil moisture at field capacity
:param β: leakage parameter
:param ks: saturated hydraulic conductivity
:param n: porosity
:param zr: soil depth
:param ew: evaporation
:param emax: maximum evapotranspiration
:param α: mean rain depth
:param λ: rain event interval (1/d)
:return: probability density function of soil moisture
"""
function pdf_soil_moisture(s, sh, sw, sstar, sfc, β, ks, n, zr, ew, emax, α, λ)

    η_w = ew / (n * zr)
    η = emax / (n * zr)
    m = ks / (n * zr * (exp(β * (1 - sfc)) - 1))
    γ = (n * zr) / α
    C = 1

    spdf = zeros(length(s))
    for i in 1:lastindex(s)
        if s[i] <= sh
            spdf[i] = 0
        elseif sh < s[i] <= sw
            p1 = (C / η_w)
            p2 = ((s[i] - sh) / (sw - sh))^((λ * (sw - sh) / η_w) - 1)
            p3 = exp(-γ * s[i])
            spdf[i] = p1 * p2 * p3
        elseif sw < s[i] <= sstar
            p1 = (C / η_w)
            p2 =
                (
                    1 + ((η / η_w) - 1) * ((s[i] - sw) / (sstar - sw))
                )^((λ * (sstar - sw)) / (η - η_w) - 1)
            p3 = exp(-γ * s[i])
            spdf[i] = p1 * p2 * p3
        elseif sstar < s[i] <= sfc
            p1 = (C / η)
            p2 = exp(-γ * s[i] + λ * (s[i] - sstar) / η)
            p3 = (η / η_w)^(λ * (sstar - sw) / (η - η_w))
            spdf[i] = p1 * p2 * p3
        else
            p1 = (C / η) * exp(-(β + γ) * s[i] + β * sfc)
            p2 =
                (
                    (η * exp(β * s[i])) / ((η - m) * exp(β * sfc) + m * exp(β * s[i]))
                )^((λ / (β * (η - m))) + 1)
            p3 = (η / η_w)^(λ * (sstar - sw) / (η - η_w))
            p4 = exp((λ / η) * (sfc - sstar))
            spdf[i] = p1 * p2 * p3 * p4
        end
    end
    ds = s[2] - s[1]
    ints = sum(spdf) * ds
    spdf2 = spdf ./ ints

    return spdf2
end


# Daily results
"""
`dt2daily(df)`

Return daily values of the soil water balance results.
:param df: dataframe with soil water balance results
"""
function dt2daily(df)
    df1 = combine(
        groupby(df, :Days),
        :Rain => sum,
        :Q => sum,
        :θ => mean,
        :Lk => sum,
        :ET => sum,
    )
    # Check if there is mulching layer
    if "ϕ" in names(df)
        df2 = combine(groupby(df, :Days), :ϕ => mean)
        df1 = rightjoin(df1, df2, on=:Days)
    end

    if "θf" in names(df)
        df3 = combine(groupby(df, :Days), :θf => mean)
        df1 = rightjoin(df1, df3, on=:Days)
    end

    if "ETf" in names(df)
        df4 = combine(groupby(df, :Days), :ETf => sum)
        df1 = rightjoin(df1, df4, on=:Days)
    end

    if "Lkf" in names(df)
        df5 = combine(groupby(df, :Days), :Lkf => sum)
        df1 = rightjoin(df1, df5, on=:Days)
    end

    # Get columns name and remove the operation
    cn = names(df1)
    cn2 = []
    for i in eachindex(cn)
        push!(cn2, split(cn[i], "_")[1])
    end

    rename!(df1, Symbol.(cn2))
    deleteat!(df1, 1)    # Delete first row (initial conditions)
    df1 = round.(df1, digits=3)
    df1[!, :Days] = Int.(df1.Days)
    return df1
end


"""
Calculate pedon scale variables

:param df: dataframe with soil water balance results
:param ϕmax: maximum porosity
"""
function pedon_scale_variables(df, ϕmax)
    df[!, :θf] = @. df[!, :θ] * (1 - (ϕmax - df[!, :ϕ]))
    df[!, :ETf] = @. df[!, :ET] * (1 - (ϕmax - df[!, :ϕ]))
    df[!, :Lkf] = @. df[!, :Lk] * (1 - (ϕmax - df[!, :ϕ])) + df[!, :Rain] * (ϕmax - df[!, :ϕ])
    df[!, :AggDensity] = @. 1 - (ϕmax - df[!, :ϕ])
    df[!, :CrackDensity] = @. (ϕmax - df[!, :ϕ])
    return df
end

end # module SoilShrinkSwell
