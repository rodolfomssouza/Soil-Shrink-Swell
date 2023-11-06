# Soil shrink-swell water balance simulations

## Docker

To build the docker image run:

```bash
docker build -t julia-soil -f Dockerfile .
```

After building the image run:

```
docker run it --rm --name soil -v "$PWD":/env julia-soil
```



