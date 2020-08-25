$1 ./Test_NavierStokesTempTurb_BuoyancyMMS.xml

python3 ./verify.py

if [ $? -eq 0 ]
then
    rm -f *.log
    rm -f *.vtk
    ls *.dat | grep -v '_ref.dat' | xargs rm -r

    exit 0
else
    exit 1
fi
