cd ..
p=$(pwd)
build="$(pwd)/build/bin/artss_serial"

array=(advection burgers diffusion diffusion/hat diffusionTurb dissipation navierStokes/channelFlow navierStokes/cavityFlow navierStokes/mcDermott navierStokes/vortex navierStokesTemp/mms navierStokesTempTurb/mms navierStokesTurb/mcDermott navierStokesTurb/vortex pressure)

for i in ${array[@]}
do
  echo -e "\n$i"
  cd ${p}/tests/${i}
  $build Test_*.xml > /dev/null
  rm *.vtk
done
