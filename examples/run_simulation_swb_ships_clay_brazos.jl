"""
Run simulation of soil water balance for the Ships Clay soil in Brazos County,
Texas with rainfall simulated from the airport weather station parameters.

"""

# %% Packages and libraries ---------------------------------------------------
using Parquet
using Revise
using SoilShrinkSwell

# Project package
include("input_parameters_ships_soil_brazos.jl")


# %% Run simulation -----------------------------------------------------------
@time begin

    # Simulation set up -------------------------------------------------------
    path_in = "examples/outputs/"
    path_out = "examples/outputs/"

    # Rainfall series ---------------------------------------------------------
    # Rainfall parameters
    α = param_rain_station.α
    λ = param_rain_station.λ
    years = 3
    df_rain = rainfall_series_poisson(α, λ, dt, years)

    # Save rainfall series
    write_parquet(
        path_out * "rainfall_series_simulated.parquet",
        df_rain,
        compression_codec="GZIP")

    # # Read rainfall simulated
    # df_rain = DataFrame(
    #     read_parquet(
    #         path_out * "rainfall_series_simulated.parquet",
    #     ),
    # )

    # Create vector of rainfall
    rain = df_rain[!, :Rain]

    # Simulations -------------------------------------------------------------
    θ0 = 0.35    # Initial conditions

    # Simulation shrink soil
    sim_shrink_soil_dt = sol_swb_ss(rain, θ0, param_shrink_soil)
    sim_shrink_soil_dt[!, :θf] =
        @. sim_shrink_soil_dt.θ * (1 - (ϕmax - sim_shrink_soil_dt.ϕ))
    sim_shrink_soil_daily = dt2daily(sim_shrink_soil_dt)
    sim_shrink_soil_daily = pedon_scale_variables(sim_shrink_soil_daily, ϕmax)

    # Simulation rigid soil
    sim_rigid_soil_dt = sol_swb_ss(rain, θ0, param_rigid_soil)
    sim_rigid_soil_dt[!, :θf] = sim_rigid_soil_dt.θ
    sim_rigid_soil_daily = dt2daily(sim_rigid_soil_dt)


    # Save results ------------------------------------------------------------
    write_parquet(
        path_out * "Simulation_swb_shrink_soil.parquet",
        sim_shrink_soil_daily,
        compression_codec="GZIP",
    )

    write_parquet(
        path_out * "Simulation_swb_rigid_soil.parquet",
        sim_rigid_soil_daily,
        compression_codec="GZIP",
    )

end
