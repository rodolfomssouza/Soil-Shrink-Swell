# Soil shrink-swell water balance simulations

This code runs simulation of soil water balance for a shrink-swell soil.

## Docker

To build the docker image run:

```bash
docker build -t julia-soil -f Dockerfile .
```

After building the image run:

```
docker run it --rm --name soil -v "$PWD":/env julia-soil
```



