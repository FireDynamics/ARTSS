$1 ./Test_Pressure.xml

rm -rf *.vtk

python3 ./verify.py

if [ $? -eq 0 ]
then
    exit 0
else
    exit 1
fi
