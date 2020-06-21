# How to compile and run ARTSS

## Requirements

- Linux
- Docker-CE 19.03 or later
- [nvidia-docker](https://github.com/NVIDIA/nvidia-docker) (Can be neglected for serial or multicore versions)


## Building the Docker image

In the directory where the Dockerfile is located run: `docker build -t artss_docker .`

It is not necessary to recreate the Docker image each time. It has to be created initially, and when a change is submitted to the Docker file.

## Building ARTSS inside the container

Navigate to the ARTSS main directory (where `compile.sh` is located).

Then run: `docker run --gpus all -it --rm -v $(pwd):/host_pwd -w /host_pwd artss_docker`

inside the container run: `nvidia-smi` to make sure you have access to all the GPUs needed

then run `compile.sh [OPTIONS]` as you normally would. `[OPTIONS]` could be for example `-c 10 -d -g`

## Building ARTSS inside the container without a GPU (serial, multicore)

Navigate to the ARTSS main directory (where `compile.sh` is located).

Then run: `docker run -it --rm -v $(pwd):/host_pwd -w /host_pwd artss_docker`

then run `compile.sh [OPTIONS]`. Don't use any GPU flags, since they will not work here.


## Running ARTSS

After compiling, ARTSS can be used as usual.


### Important Docker commands

- `exit` -> Leave container
- `docker container ls` -> List all running Docker containers
- `docker ps -a` -> List all Docker containers
- `docker image ls` -> List all Docker images
- `docker container rm [tag]` -> Removes Docker container with tag [tag]
- `docker image rm [tag]` -> Removes Docker image with tag [tag]

More Docker commands: https://docs.docker.com/engine/reference/commandline/cli/


## Using Docker with `compile.sh`

The Docker image can be built with `./compile.sh --docker-build`

To build ARTSS inside the container use `--docker-run` or `--docker-run-cpu`. To change the hostname of the Docker add the option `--docker-hostname`.

The Docker commands can be combined with each other.
