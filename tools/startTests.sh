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
COMPILE="mpi"
while getopts ":x:y:z:hgmM" opt; do
  case $opt in
    x)
      X=$OPTARG
      ;;
    y)
      Y=$OPTARG
      ;;
    z)
      Z=$OPTARG
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
    M)
      COMPILE="mpi"
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
  sed -i -e "s@<MESHX> 1 </MESHX>@<MESHX> $X </MESHX>@" Test_*.xml
  sed -i -e "s@<MESHY> 1 </MESHY>@<MESHY> $Y </MESHY>@" Test_*.xml
  sed -i -e "s@<MESHZ> 1 </MESHZ>@<MESHZ> $Z </MESHZ>@" Test_*.xml
  mpiexec -n $PROCS $build Test_*.xml
  rm *.vtk *.log
  sed -i -e "s@<MESHX> $X </MESHX>@<MESHX> 1 </MESHX>@" Test_*.xml
  sed -i -e "s@<MESHY> $Y </MESHY>@<MESHY> 1 </MESHY>@" Test_*.xml
  sed -i -e "s@<MESHZ> $Z </MESHZ>@<MESHZ> 1 </MESHZ>@" Test_*.xml
done
