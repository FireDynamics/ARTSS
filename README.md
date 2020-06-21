# ARTSS
Accelerator-based Real Time Smoke Simulator

ARTSS is a real-time and prognosis capable CFD code basis simulating buoyancy-driven turbulent smoke spread
based on finite differences and a large eddy simulation turbulence model. The open source code is portable on CPU and GPU and successfully verified using analytical and
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
|          | Visualization of output                             | vtk      |   >= 5.8      |
|          | Testing for consistency of output while developing  | Python   |   >= 3.6      |
| Specific | Compiler fully supporting C++-17 and OpenACC        | PGI      |   >= 19.10    |

### Compiling the Code
Once the code has been checked out and all required software has been installed, JuROr
can be built from the terminal by first running cmake to configure the build, then
running make. The steps are summarized below.  

1. Clone
```
git clone https://github.com/FireDynamics/ARTSS.git
cd ARTSS
```
2. Compiling the code

```
./compile.sh [OPTIONS]
```
*Note: Without options all executables will be compiled using the PGI compiler and CUDA 10.1.*

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


### Code Structure
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
* [**My Linh Würzburger**](https://www.fz-juelich.de/ias/ias-7/EN/AboutUs/Staff/Current/Wuerzburger_My_Linh/main.html?nn=2302136): *Principal developer*
* [**Christian von Mach**]: *Developer*
* [**Max Böhler**](https://www.fz-juelich.de/ias/ias-7/EN/AboutUs/Staff/Current/Boehler_Max/_node.html): *Contributor*
* Suryanarayana Maddu Kondaiah
