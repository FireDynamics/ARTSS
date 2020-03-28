# How to compile and run ARTSS

## Requirements

- Linux
- Docker-CE 19.03 or later
- [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) (Can be neglected for serial or multicore versions)


## Building the Docker image

In the directory where the Dockerfile is located run: `docker build -t artss_docker .`

## Building ARTSS inside the container

Navigate to the ARTSS main directory (where `compile.sh` is located).

Then run: `docker run --gpus all -it --rm -v $(pwd):/host_pwd -w /host_pwd artss_docker`

inside the container run: `nvidia-smi` to make sure you have access to all the GPUs needed

then run `compile.sh [OPTIONS]` as you normally would. `[OPTIONS]` could be for example `-c 10 -d -g`

## Building ARTSS inside the container without a GPU (serial, multicore)

Navigate to the ARTSS main directory (where `compile.sh` is located).

Then run: `docker run -it --rm -v $(pwd):/host_pwd -w /host_pwd artss_docker`

then run `compile.sh [OPTIONS]`. Don't use any GPU flags, since they will not work here.


## Running ARTSS on the GPU

after compiling the binaries should be able to run on all GPUs visible in `nvidia-smi`

example: `cd tests/burgers && bash run.sh ../../build/bin/artss_gpu`

[docker-build-output]:docker-build-output.png
