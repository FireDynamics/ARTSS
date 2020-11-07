#!/usr/bin/env bash
COMPILE=""
GPU=1
JURECA=1
P100=1
GPU_CC="cc35"
CUDA_VERSION=10.1
BUILDTYPE=Release

YELLOW='\033[1;33m'
NC='\033[0;m'

DESCRIPTION="Script for compiling ARTSS. Select one or multiple executables for compiling (default: all executables). Specify CUDA version (default: $CUDA_VERSION) and GPU model (default: $GPU_CC). For a parallel execution choose --jobs, if there is no integer after the --jobs option, the number of processing units available to the current process will be used.
"
OPTIONS="
Available Options:

Load modules:
  ${YELLOW}--jureca${NC}                       \t load modules for JURECA
  ${YELLOW}--p100${NC}                         \t setup for P100

Executables:
  Production (with data output, visualization and analysis):
   ${YELLOW}-s${NC}
  ${YELLOW}--serial${NC} 
  ${YELLOW}--artss_serial${NC}                  \t Executable: artss

  ${YELLOW} -m${NC}
  ${YELLOW}--multicore${NC}
  ${YELLOW}--artss_multicore_cpu${NC}           \t Executable: artss_multicore_cpu

   ${YELLOW}-g${NC}
  ${YELLOW}--gpu${NC}
  ${YELLOW}--artss_gpu${NC}                      \t Executable: artss_gpu
----  
  Benchmarking (without output, visualization, analysis):
  ${YELLOW}--sb${NC}
  ${YELLOW}--serial_benchmark${NC}
  ${YELLOW}--artss_benchmark${NC}                  \t Executable artss_serial_benchmark

  ${YELLOW}--mb${NC}
  ${YELLOW}--multicore_benchmark${NC}
  ${YELLOW}--artss_multicore_cpu_benchmark${NC}    \t Executable artss_multicore_cpu_benchmark

  ${YELLOW}--gb${NC}
  ${YELLOW}--gpu_benchmark${NC}
  ${YELLOW}--artss_gpu_benchmark${NC}               \t Executable artss_gpu_benchmark

Other:
   ${YELLOW}-c${NC}
  ${YELLOW}--cudaversion${NC}                     \t set CUDA Version
  ${YELLOW}--cc${NC}
  ${YELLOW}--computecompatibility${NC}            \t set compute compability of the GPU (30|35|50|60|70|75)
   ${YELLOW}-d${NC}
   ${YELLOW}--debugmode${NC}                      \t set debug flag for build type (default: ${BUILDTYPE})

  ${YELLOW}--jobs${NC}                            \t set the number of recipes to execute at once (-j/--jobs flag in make)

  ${YELLOW}--gcc${NC}                             \t use gcc as compiler (optional: specify version)
  ${YELLOW}--pgi${NC}                             \t use pgcc ac compiler (optional: specify version)

Docker - ! cannot be combined with other commands ! (more information about docker commands in docker/README.md):
  ${YELLOW}--docker-build${NC}                    \t build docker image
  ${YELLOW}--docker-hostname${NC}                 \t set hostname for docker image (default: docker_\$(hostname))
  ${YELLOW}--docker-run${NC}                      \t run docker with gpu support
  ${YELLOW}--docker-run-cpu${NC}                  \t run docker without gpu support
"

HELP="$DESCRIPTION$OPTIONS"
COMPILE=""
DOCKERRUN=1
DOCKERBUILD=1
DOCKERHOST=docker_$(hostname)
DOCKERRUNCPU=1
PROCS=-1
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -c|--cuda)
      CUDA_VERSION="$2"
      shift
      shift
      ;;
    --cc|--computecompatibility)
      GPU_CC="cc$2"
      shift
      shift
      ;;
    -d|--debug|--debugmode)
      BUILDTYPE=Debug 
      shift
      ;;
    --docker-build)
      DOCKERBUILD=0
      shift
      ;;
    --docker-hostname)
      DOCKERHOST=$2
      shift
      shift
      ;;
    --docker-run)
      DOCKERRUN=0
      shift
      ;;
    --docker-run-cpu)
      DOCKERRUNCPU=0
      shift
      ;;
    -g|--gpu|--artss_gpu)
      COMPILE="$COMPILE artss_gpu" 
      GPU=0
      shift
      ;;
    --gb|--gpu_benchmark|--artss_gpu_benchmark)
      COMPILE="$COMPILE artss_gpu_benchmark "
      GPU=0
      shift
      ;;
    --gcc)
      COMPILER="GCC"
      if [[ $2 != -* ]]
      then
        GCC_VERSION="$2"
        shift
      fi
      shift
      ;;
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
    -j|--jobs)
      if [[ $2 =~ ^-?[0-9]+$ ]]
      then
        PROCS=$2
        shift
      else
        PROCS=$(nproc)
      fi
      shift
      ;;
    -m|--multicore|--artss_multicore_cpu)
      COMPILE="$COMPILE artss_multicore_cpu"
      GPU=0
      shift
      ;;
    --mb|--multicore_benchmark|--artss_multicore_cpu_benchmark)
      COMPILE="$COMPILE artss_multicore_cpu_benchmark"
      GPU=0
      shift
      ;;
    --pgi)
      COMPILER="PGI"
      if [[ $2 != -* ]]
      then
        PGI_VERSION="$2"
        shift
      fi
      shift
      ;;
    -s|--serial|--artss_serial)
      COMPILE="$COMPILE artss_serial"
      shift
      ;;
    --sb|--serial_benchmark|--artss_serial_benchmark)
      COMPILE="$COMPILE artss_serial_benchmark"
      shift
      ;;
    --jureca)
      JURECA=0
      shift
      ;;
    --p100)
      P100=0
      shift
      ;;
    *)
      POSITIONAL+=("$1")
      echo "unknown option: $1"
      shift
      ;;
  esac
done

if [ $DOCKERBUILD -eq 0 ]
then
  cd docker
  docker build -t artss_docker --no-cache .
  cd ..
fi

if [ $DOCKERRUN -eq 0 ]
then
  docker run --gpus all -it --rm --hostname=${DOCKERHOST} -v $(pwd):/host_pwd -w /host_pwd artss_docker
fi

if [ $DOCKERRUNCPU -eq 0 ]
then
  docker run -it --rm --hostname=${DOCKERHOST} -v $(pwd):/host_pwd -w /host_pwd artss_docker # /bin/bash -c "./compile.sh"
fi

if [ $DOCKERRUN -eq 0 -o $DOCKERRUNCPU -eq 0 -o $DOCKERBUILD -eq 0 ]
then
  exit
fi

if [[ $JURECA -eq 1 && $P100 -eq 1 ]]
then
  HOSTNAME=$(hostname)
  if [[ $HOSTNAME = jrl* ]]; then JURECA=0; fi
  if [ "$HOSTNAME" = "ias7139" ]; then P100=0; fi
fi

if [ "$COMPILE" = "" ]
then
  GPU=0
fi
if [ -z $COMPILER ]
then
  COMPILER="PGI"
fi

if [ $JURECA -eq 0 ]
then
  module load CMake
  module load PGI/19.10-GCC-8.3.0
  module load CUDA/10.1.105
  export CUDA_LIB=${CUDA_ROOT}/lib64/
  export CUDA_INC=${CUDA_ROOT}/include/
  CUDA_VERSION=10.1
  GPU_CC=cc35
  GPU=0
fi

if [ ${P100} -eq 0 ]
then
  if [ -z "${PGI_VERSION}" ]
  then
    PGI_VERSION=19.4
  fi
  if [ -z ${CUDA_VERSION} ]
  then
    CUDA_VERSION=10.1
  fi
  GPU_CC=cc60
fi

if [ $GPU -eq 0 ]
then
  export CUDA_LIB=$CUDA_ROOT/lib64
  export CUDA_INC=$CUDA_ROOT/include
fi

rm -rf build/
mkdir build
cd build || exit

if [[ $GPU -eq 0 ]] || [[ $COMPILER = "PGI" ]]
then
  CCOMPILER=pgcc
  CXXCOMPILER=pgc++
else
  CCOMPILER=gcc
  CXXCOMPILER=g++
fi

if [ -z ${CUDA_VERSION} ]
then
  cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE=${BUILDTYPE} -DCMAKE_C_COMPILER=${CCOMPILER} -DCMAKE_CXX_COMPILER=${CXXCOMPILER} ..
else
  cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE=${BUILDTYPE} -DCMAKE_C_COMPILER=${CCOMPILER} -DCMAKE_CXX_COMPILER=${CXXCOMPILER} -DGPU_CC=${GPU_CC} -DCUDA_VERSION=${CUDA_VERSION} ..
fi

if [ "$PROCS" -le 0 ]
then
  make $COMPILE
else
  echo "-- Parallel execution with $PROCS processing units"
  make $COMPILE -j"$PROCS"
fi
