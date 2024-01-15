# Soil shrink-swell water balance simulations

This code runs simulation of soil water balance for both shrink-swell and rigid soils, complementary to the paper
**Dynamic coupling between soil properties and water content in shrink-swell soils: effects on surface hydrologic partitioning**.
Please read the paper for a detailed model explanation and parameterization.

This code is written in [Julia](https://julialang.org/) with the goal of having functions in a package, and imported later when designing a simulation.
The main file is [SoilShrinkSwell.jl](src/SoilShrinkSwell.jl).
Besides the functions individually, this module also contains functions to solve the model numerically using the [Euler method](https://en.wikipedia.org/wiki/Euler_method).

A simple example of how to use this code is shown in [run_simulations.jl](examples/run_simlation_swb.jl). This code imports the parameters from the file [input_parameters.jl](examples/input_parameters.jl), which are the same parameters used in the paper mentioned before. A stochastic rainfall series will be generated using the marked Poisson process. Then, the soil moisture model will be executed. The results will be saved in the folder [examples/outputs](examples/outputs) in the [parquet](https://en.wikipedia.org/wiki/Apache_Parquet) format. To visualize the results, run the code [figures_swb_simulation.jl](examples/figures_swb_simulation.jl).

## Requirements

This code was written using Julia v.1.9 and tested around December 2023. Updates in the language and dependencies might result in some errors.
The dependencies versions are listed in the file [Project.toml](Project.toml).

To run the code, install [Julia](https://julialang.org/downloads/) and add this package.

```{julia}
using Pkg
Pkg.add(url="https://github.com/rodolfomssouza/Soil-Shrink-Swell")
```

# Alternative

Additionally, you can use Docker to build an image and run the Soil Shrink Swell simulations

## Docker

To build the docker image run:

```bash
docker build -t julia-soil -f Dockerfile .
```

After building the image run:

```
docker run it --rm --name soil -v "$PWD":/env julia-soil
```

## Contacts

For any questions, comments or suggestions, please email us:

- Rodolfo Souza: rodolfo.souza@tamu.edu
- Salvatore Calabrese: salvatore.calabrese@ag.tamu.edu

