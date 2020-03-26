# ARTSS
Accelerated Based Real Time Smoke Simulation

ARTSS is a real-time and prognosis capable CFD code basis simulating buoyancy-driven turbulent smoke spread
based on finite differences and a large eddy simulation turbulence model, being performance portable on CPU
and GPU contained in one expandable, open source code, successfully verified using unit, analytical and
semi-analytical tests, and successfully validated with scenarios relevant for fire protection. 

## Getting Started

### Requirements
The serial CPU version of ARTSS can be compiled on Linux or OS X systems with very few tools, 
whereas the multicore and GPU version need an OpenACC capable compiler. 
Detailed requirements are listed in the table below (general requirements for serial version, specific for multicore and GPU version).

|          | Purpose                                             | Tool     | Version       |
|--------- | --------------------------------------------------- | -------- | --------------|
| General  | Version control system to obtain the source code    | git      |   >= 2.0      |
|          | Build processor using a compiler-independent method | CMake    |   >= 2.8      |
|          | Compiler fully supporting C++-11                    | gcc      |   >= 4.7      |
|          |                                                     | clang    |   >= 6.1      |
|          | Visualization of output                             | vtk      |   >= 5.8      |
|          |                                                     | Paraview/|   VisIT       |
|          | Testing for consistency of output while developing  | python   |   >= 2.7      |
| Specific | Compiler fully supporting C++-11 and OpenACC        | PGI      |   >=17.3/17.4 |
|          | Profiling (using NVTX)                              | CUDA     |   >= 8.0      |


### Arguments
The non PROFILING-version accepts folling parameter provided via command line

| Flag |                                      Purpose                                               | Argument |
| :--: | :----------------------------------------------------------------------------------------- | :------- |
|  -l  | determins the minmal loglevel you will only see messages greater that                      | loglevel |
|  -o  | determins the path for the logfile if not provided it uses './log.txt'. Use '-' for stdout | path     |


loglevel can have following arguments.

| Loglevel | Default |
| :------- | :-----: |
| trace    |         |
| debug    |         |
| info     |    X    |
| warning  |         |
| error    |         |
| critical |         |


#### Examples

Creates `log.txt` if not existing or appends on it. Shows only messages on level info or above
```
../build/bin/artss_serial XML.xml
```

Creates `myfile.txt` if not existing or appends on it. Shows only messages on level info or above
```
../build/bin/artss_serial -o myfile.txt XML.xml
```

Shows even debug messages on stdout
```
../build/bin/artss_serial -l debug -o - XML.xml
```
