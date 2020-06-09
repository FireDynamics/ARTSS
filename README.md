# ARTSS
Accelerated Based Real Time Smoke Simulation

ARTSS is a real-time and prognosis capable CFD code basis simulating buoyancy-driven turbulent smoke spread
based on finite differences and a large eddy simulation turbulence model. The open source code is portable on CPU and GPU and successfully verified using unit, analytical and
semi-analytical tests. It is also successfully validated with scenarios relevant for fire protection.
ARTSS is based on JuROr, which was originally developed within the [ORPHEUS](http://www.orpheus-projekt.de) project
(funded by [BMBF](https://www.bmbf.de/)) by Dr. Anne Küsters.

## Getting Started

### Requirements
The serial CPU version of ARTSS can be compiled on Linux or MacOS systems with very few tools,
whereas the multicore and GPU version needs an OpenACC capable compiler.
Detailed requirements are listed in the table below (general requirements for serial version, specific for multicore and GPU version).

|          | Purpose                                             | Tool     | Version       |
|--------- | --------------------------------------------------- | -------- | --------------|
| General  | Version control system (optional)                   | git      |   >= 2.0      |
|          | Build processor using a compiler-independent method | CMake    |   >= 2.8      |
|          | Compiler fully supporting C++-17 (gcc or clang)     | gcc      |   >= 7.0      |
|          |                                                     | clang    |   >= 8.0      |
|          | Visualization of output                             | vtk      |   >= 5.8      |
|          |                                                     | Paraview |   >= 5.50     |
|          | Testing for consistency of output while developing  | Python   |   >= 3.6      |
| Specific | Compiler fully supporting C++-17 and OpenACC        | PGI      |   >= 19.10    |

### Compiling the Code
Once the code has been checked out and all required software has been installed, JuROr
can be built from the terminal by first running cmake to configure the build, then
running make. The steps are summarized below.  

1. Clone
```
git clone --recurse-submodules  https://github.com/FireDynamics/ARTSS.git
cd ARTSS
```
2. Compiling the code

```
./compile.sh [OPTIONS]
```
*Note: Without options all executables will be compiled using the PGI compiler and CUDA 8.0.*

OPTIONS (selection; show all by using --help flag):
- '-s' -> Compile serial ARTSS version
- '-m' -> Compile multicore ARTSS version
- '-g' -> Compile GPU ARTSS version
- '--jobs [tasks]' -> Specifies  the  number of tasks to run simultaneously (default $(nproc))
- '--gcc' Use GCC as compiler

EXAMPLE:
- Compile serial version of ARTSS using 4 cores and GCC
`./compile.sh -s --jobs 4 --gcc`
- Compile multicore version of ARTSS
`./compile.sh -m`

Extra:
It is also possible to compile ARTSS using a docker file. Instructions and further information can be found in the folder `docker`.


### Code structure
```
ARTSS
│   compile.sh
│   README
│   CMakeLists.txt
│   LICENSE   
│
└───examples
│   # simple examples to demonstrate the capabilities of ARTSS
|
└───src
│   # source code
│
└───tests
│   # files to test individual modules of ARTSS
│  
└───tools
│   # different tools that make your life easier
```

## Authors
* [**Anne Küsters**](https://www.fz-juelich.de/SharedDocs/Personen/IAS/JSC/EN/staff/kuesters_a.html?nn=361682): *Initial work*
* [**Lukas Arnold**](https://www.fz-juelich.de/ias/ias-7/EN/AboutUs/Staff/Current/Arnold_Lukas/main.html): *Contributor* and *Supervisor*
* [**My Linh Würzburger**](https://www.fz-juelich.de/ias/ias-7/EN/AboutUs/Staff/Current/Wuerzburger_My_Linh/main.html?nn=2302136): *Contributor*
* [**Max Böhler**](https://www.fz-juelich.de/ias/ias-7/EN/AboutUs/Staff/Current/Boehler_Max/_node.html)
* Suryanarayana Maddu Kondaiah

<!--
### Arguments
The non PROFILING-version accepts following parameter provided via command line

| Flag  |                                      Purpose                                                | Argument |
| :---: | :------------------------------------------------------------------------------------------ | :------- |
|  -l   | determines the minimal loglevel you will only see messages greater that                     | loglevel |
|  -o   | determines the path for the logfile if not provided it uses './log.txt'. Use '-' for stdout | path     |


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
-->

## Development
If you want to participate in this project, this section is for you.

### Coding style
https://google.github.io/styleguide/cppguide.html

### Development Requirements
These additional tools can make your life a lot easier.

|     Tool   | Version |                Ref                 |
| :--------- | :-----: | :--------------------------------- |
| cpplint    |         | https://github.com/cpplint/cpplint |
| pre-commit |         | https://pre-commit.com/            |
| spdlog     |         | https://github.com/gabime/spdlog   |
