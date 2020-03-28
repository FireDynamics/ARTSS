# How to compile and run ARTSS

## Requirements

- Linux
- Docker-CE 19.03 or later
- [nvidia-docker](https://github.com/NVIDIA/nvidia-docker)

*it might work with older versions of Docker and nvidia-docker2 but it was not tested*

## Building the Docker image

In the directory where the Dockerfile is located run: `docker build .`

output should look like this:

![example build output][docker-build-output]

## Building ARTSS inside the container

navigate to the ARTSS main directory (where `compile.sh` is located).

Then run: `docker run --gpus all -it --rm -v $(pwd):/host_pwd -w /host_pwd {DOCKERIMAGEID}`

`{DOCKERIMAGEID}` has to be replaced by the image id from the docker build, in this example it would be `8a7bf3be21d0`

*for easier usage the Image should be tagged with a name, see [tagging](https://docs.docker.com/engine/reference/commandline/tag/)*

inside the container run: `nvidia-smi` to make sure you have access to all the GPUs needed

then run `compile.sh {params}` as you normally would, `{params}` could be for example `-c 10 -d -g`

## Running ARTSS on the GPU

after compiling the binaries should be able to run on all GPUs visible in `nvidia-smi`

example: `cd tests/burgers && bash run.sh ../../build/bin/artss_gpu`

[docker-build-output]:docker-build-output.png