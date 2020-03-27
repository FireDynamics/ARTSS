HELP="
Run and compile all test files in tests folder and compare it with their corresponding reference files.

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
execRef=$(pwd)/tools/compareReferenceFiles.sh

array=(advection burgers diffusion diffusion/hat diffusionTurb dissipation navierStokes/channelFlow navierStokes/cavityFlow navierStokes/mcDermott navierStokes/vortex navierStokesTemp/mms navierStokesTempTurb/mms navierStokesTurb/mcDermott navierStokesTurb/vortex pressure)

for i in ${array[@]}
do
  echo -e "\n$i"
  cd ${p}/tests/${i}
  if [ -f u.dat ]; then rm u.dat; fi
  if [ -f v.dat ]; then rm v.dat; fi
  if [ -f w.dat ]; then rm w.dat; fi
  if [ -f p.dat ]; then rm p.dat; fi
  if [ -f T.dat ]; then rm T.dat; fi
  if [ -f rho.dat ]; then rm rho.dat; fi
  ${build} Test*.xml > /dev/null
  rm *.vtk
  ${execRef}
done
