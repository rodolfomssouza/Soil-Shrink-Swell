# Soil shrink-swell water balance simulations

This code runs simulations of the soil water balance for both shrink-swell and rigid soils, following the model described in the paper **Dynamic coupling between soil properties and water content in shrink-swell soils: effects on surface hydrologic partitioning**.
Please read the paper for a detailed model explanation and parameterization.

This code is written in [Julia](https://julialang.org/) with the goal of creating a package, which can be readily incorporated into larger projects.
The main file is [SoilShrinkSwell.jl](src/SoilShrinkSwell.jl).
This module contains the equations of the model and the numerical solution using the [Euler method](https://en.wikipedia.org/wiki/Euler_method).

A simple example of how to use this code is shown in [run_simulations.jl](examples/run_simlation_swb.jl) and in [notebook simulation](examples/example_swb_simulation.ipynb).
This code imports the parameters from the file [input_parameters.jl](examples/input_parameters.jl), which are the same parameters used in the published paper.
A stochastic rainfall series will be generated using the marked Poisson process. Then, the soil moisture model will be executed.
The results will be saved in the folder examples/outputs in the parquet format.
To visualize the results, run the code [figures_swb_simulation.jl](examples/figures_swb_simulation.jl).


## Requirements
This code was written using Julia v.1.9 and last tested in December 2023
Updates in the language and dependencies might result in some errors.
The dependencies versions are listed in the file [Project.toml](Project.toml).

To run the code, install [Julia](https://julialang.org/downloads/) and add this package.

```{julia}
using Pkg
Pkg.add(url="https://github.com/rodolfomssouza/Soil-Shrink-Swell")
```

# Alternative

Additionally, you can use [Docker](https://www.docker.com/) to build an image and run the Soil Shrink Swell simulations.

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

