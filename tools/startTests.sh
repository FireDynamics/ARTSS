HELP="
Run all test files in tests folder.

-g|--gpu\trun gpu version
-m|--multicore\trun multicore version
-M|--mpi\trun mpi only version
-x <integer>\tset number of meshes in x direction
-y <integer>\tset number of meshes in y direction
-z <integer>\tset number of meshes in z direction
"

for arg in "$@"; do
  shift
  case "$arg" in
    "--help")       set -- "$@" "-h" ;;
    "--gpu")        set -- "$@" "-g" ;;
    "--multicore")  set -- "$@" "-m" ;;
    "--mpi")        set -- "$@" "-M" ;;
    *)        set -- "$@" "$arg"
  esac
done

X=1
Y=1
Z=1
MPI=0
COMPILE="serial"
while getopts ":x:y:z:hgm" opt; do
  case $opt in
    x)
      X=$OPTARG
      MPI=1
      ;;
    y)
      Y=$OPTARG
      MPI=1
      ;;
    z)
      Z=$OPTARG
      MPI=1
      ;;
    h)
      echo -e "$HELP"
      ;;
    g)
      COMPILE="gpu"
      ;;
    m)
      COMPILE="multicore_cpu"
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

PROCS=$(($X*$Y*$Z))

cd ..
p=$(pwd)
build="$(pwd)/build/bin/artss_$COMPILE"

array=(advection burgers diffusion diffusion/hat diffusionTurb dissipation navierStokes/channelFlow navierStokes/cavityFlow navierStokes/mcDermott navierStokes/vortex navierStokesTemp/mms navierStokesTempTurb/mms navierStokesTurb/mcDermott navierStokesTurb/vortex pressure)

for i in ${array[@]}
do
  echo -e "\n$i"
  cd ${p}/tests/${i}
  if [ $MPI -eq 1 ]
  then
    mpiexec -n $PROCS $build Test_*.xml -x $X -y $Y -z $Z
  else
    $build Test_*.xml
  fi
  rm *.vtk *.log
done
