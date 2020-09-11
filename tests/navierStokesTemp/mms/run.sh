$1 ./Test_NavierStokesTemp_BuoyancyMMS.xml

python3 ./verify.py

rm -f *.log
rm -f *.vtk
ls *.dat | grep -v '_ref.dat' | xargs rm -r

if [ $? -eq 0 ]
then
    exit 0
else
    exit 1
fi
