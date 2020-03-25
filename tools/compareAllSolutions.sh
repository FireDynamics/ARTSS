cd ..
p=$(pwd)
build="$(pwd)/build/bin/artss_serial"
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
