HELP="
Run all test files in tests folder.

-g|--gpu\trun gpu version
-m|--multicore\trun multicore version
-s|--serial\trun serial version
"
COMPILE="serial"
while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo -e "$HELP"
      exit
      ;;
    -g|--gpu)
      COMPILE="gpu"
      shift
      ;;
    -m|--multicore)
      COMPILE="multicore_cpu"
      shift
      ;;
    -s|--serial)
      COMPILE="serial"
      shift
      ;;
  esac
done
cd ..
p=$(pwd)
build="$(pwd)/build/bin/artss_$COMPILE"

array=(advection burgers diffusion diffusion/hat diffusionTurb dissipation navierStokes/channelFlow navierStokes/cavityFlow navierStokes/mcDermott navierStokes/vortex navierStokesTemp/mms navierStokesTempTurb/mms navierStokesTurb/mcDermott navierStokesTurb/vortex pressure)

for i in ${array[@]}
do
  echo -e "\n$i"
  cd ${p}/tests/${i}
  $build Test_*.xml > /dev/null
#  rm *.vtk
done
