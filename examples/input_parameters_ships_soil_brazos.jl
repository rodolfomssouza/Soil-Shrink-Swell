"""
Soil parameters to be imported in files that runs simulations
"""


# Base parameters -------------------------------------------------------------
ϕmax = 0.56
ϕmin = 0.22
ksmax = 1.80
b = 7.72963
ϵ = 0.38    # Original - Ships Clay
q = 2.3     # Original - Ships Clay

# Shrink soil water retentation curve parameters - Fitted
θh_shrink = 0.172
θw_shrink = 0.220
θstar_shrink = 0.356
θfc_shrink = 0.360

# Rigid soil water retentation curve parameters - Fitted
θh_rigid = 0.163
θw_rigid = 0.208
θstar_rigid = 0.341
θfc_rigid = 0.345

emax = 0.5
ew = 0.05
zr = 30
dt = 1 / 24


# Shrink soil -----------------------------------------------------------------
param_shrink_soil = (
    θh=θh_shrink,
    θw=θw_shrink,
    θstar=θstar_shrink,
    θfc=θfc_shrink,
    emax=emax,
    ew=ew,
    ksmax=ksmax,
    b=b,
    ϕmax=ϕmax,
    ϕmin=ϕmin,
    ϵ=ϵ,
    q=q,
    zr=zr,
    dt=dt,
    nsave="",
)


# Rigid soil ------------------------------------------------------------------
param_rigid_soil = (
    θh=θh_rigid,
    θw=θw_rigid,
    θstar=θstar_rigid,
    θfc=θfc_rigid,
    emax=emax,
    ew=ew,
    ksmax=ksmax,
    b=b,
    ϕmax=ϕmax,
    ϕmin=ϕmin,
    ϵ=ϵ,
    q=0,
    zr=zr,
    dt=dt,
    nsave="",
)

# Rainfall parameters ---------------------------------------------------------
param_rain_chirps = (
    α=1.832,
    λ=0.175,
)

param_rain_station = (
    α=1.116,
    λ=0.257,
)
