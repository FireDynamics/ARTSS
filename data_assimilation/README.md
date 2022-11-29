# Data assimilation
Data assimilation is a technique which, generally speaking, combines observational data with numerical simulation. In the case of ARTSS, it is now possible to directly or indirectly modify the simulation data during runtime with the help of TCP. The data can be changed directly by rewriting the respective cells or indirectly by changing parameters.

## Configuration
`pip install -r requirements.txt`

## How to use
1. Compile target data_assimilation (`./compile.sh --mpi`)
1. Go to directory `ARTSS/data_assimilation/example`
2. Start ARTSS with either `mpirun` `or mpiexec`
`mpirun --np 2 path_to_build_dir/bin/artss_data_assimilation_serial Test_NavierStokesTempTurb_DataAssimilation.xml`
`mpiexec --np 2 path_to_build_dir/bin/artss_data_assimilation_serial Test_NavierStokesTempTurb_DataAssimilation.xml`
3. Run Python script from example directory
`python3 ../main.py`

