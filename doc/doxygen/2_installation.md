\page installation Installation

* [Prerequisites](#Prerequisites)
* [Installation](#Installation)
    - [Quick installation](#Quick_installation)
        + [Linux](#Linux)
        + [Windows](#Windows)
            * [Generate project files](#Generate_project_files)
                - [ Using GUI](#Using_GUI)
                - [ Using command line](#Using_command_line)
            * [Build the project](#Build_the_project)
                - [With Visual Studio](#With_Visual_Studio)
                - [With NMake](#With_NMake)
    - [Customizing installation](#Customizing_installation)
        + [Project specific options and variables](#Project_specific_options_and_variables)
            * [Compiler options](#Compiler_options)
    - [Installation of the example 'refine'](#Installation_of_the_example_refine)
        + [Building refine during DiSCaMB compilation](#Building_refine_during_DiSCaMB_compilation)
        + [Building refine in separate step](#Building_refine_in_separate_step)


<a name="Prerequisites"></a>
## Prerequisites 

Required:

- C++ compiler (DiSCaMB has been tested with MSVC 9.0 (VS 2008), 14.0 (VS 2015) and 14.1 (VS 2017), gcc 5.4<br>
  Intel C++ compiler 14.0 and Clang 3.8)
- [Cmake](https://cmake.org/) (version 2.8 and above)
- build system - a tool managing build processes - i.e. compilation, linking, installation e.t.c. <br>
  it has to correspond to one of cmake generators, e.g. for Linux gnu make is popular choice, for <br>
  MS Visual Studio compiler the Visual Studio itself is the build system or alternatively nmake, <br>
  other choices are e.g. Ninja, CodeBlocks, CodeLite, Eclipse, KDevelop, Kate - run cmake --help to <br>
  list of available generators (and build systems)
      
Optionally:

- for compiling GPU part - NVIDIA Kepler (GK110 or GK110b, arch sm_35) GPU (or newer), CUDA Toolkit ( >= 6.0)
- for building some examples (slater_type_atomic_wfn, custom_model and type_assignment_info) compiler supporting C++11 standard
- for building documentation - [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html), latex and perl
- for building refinement example - [Alglib](http://www.alglib.net/) library (C++ version)

<a name="Installation"></a>
##Installation 

DiSCaMB uses [CMake](https://cmake.org/) to generate project files or makefiles for a particular development environment. <br/>
Cmake can be used from command line or via GUI. In order to use DiSCaMB it should be not only build but <br>
also installed (otherwise header files build during installation will be not placed in include directory).<br>
If Cmake is run more than once (e.g. to create project files with different than initially used options)<br>
it might be necessary to delete CMakeCache.txt file generated previously in build directory.

<a name="Quick_installation"></a>
### Quick installation 

<a name="Linux"></a>
#### Linux

Assuming that DiSCaMB distribution is located at /some_path/DiSCaMB:
```
cd /some_path 
mkdir build
cd build
cmake ../DiSCaMB
make
make install
```
The resulting files can be found in /some_path/build/build. The library location is<br>
/some_path/build/build/lib/libdiscamb.a. Header files are in /some_path/build/build/include/discamb.<br>
Tests can be run by calling:  
```
./build/tests/test_1 -all
```
(it can take few minuts), if they are all passed then 'ALL TESTS PASSES' messege wil be print.
To install call:
```
sudo make install
```
In this case DiSCaMB will be installed in default localization /usr/local. 
<a name="Windows"></a>
#### Windows

Installation with Microsoft Visual Studio is decribed here. 
<a name="Generate_project_files"></a>
#####Generate project files:
<a name="Using_GUI"></a>
###### Using GUI

Open CMake (cmake-gui). If you would like to use nmake generator (build from command line, probably <br>
faster) open CMake from developer command prompt (console with visual studio environment variables <br>
preset) by calling cmake-gui. Fill 'Where is the source code' and 'Where to build the binaries' filds. <br>
A window will pop up allowing for generator specification, on default some version of Visual Studio <br>
is chosen when present. Change to nmake if you prefere command line build instead of building from <br>
Visual Studio IDE (integrated developer environment - GUI for coding and program buiding). Click <br>
configure. CMake window will fill with variables which can be modified to customize build. Here you <br>
can customize build options. Click generate. The project files will be generated in the build directory.<br>
If Visual Studio was used clicking 'Open Project' will open the project in Visual Studio. See Build <br>
the project paragraph' for next steps. 
<a name="Using_command_line"></a>
###### Using command line 

Open command prompt (click WIN+r to open 'run' dialog box and call cmd to start the command line). Make<br>
 folder when you would like to build the library e.g.:
```
 mkdir build
 ```
 change to that folder 
 ```
 cd build
 ```
 run cmake 
 ```
 cmake path/to/DiSCaMB
 ```
 It should create project files in the build folder (for Visual Studio IDE). If you prefere building<br> 
 with nmake (building with Microsoft compiler from command line) use developer command prompt for Visual <br>
 Studio and run instead:<br>
 ```
 cmake path/to/DiSCaMB -G"NMake Makefiles" [options]
 ```
 where the additional options may e.g. -DCMAKE_BUILD_TYPE=release which will make release version of<br>
 the library.
<a name="Build_the_project"></a>
##### Build the project
<a name="With_Visual_Studio"></a>
###### With Visual Studio
Assuming that you use Visual Studio click on discamb.sln in the build directory (if the is not open <br>
in Visual Studio yet). On the menu bar, choose 'Build' and then 'Build Solution' ([see also here](https://msdn.microsoft.com/en-us/library/5tdasz7h.aspx))<br>
Resulting files are kept in build subdirectory of the build folder. For installation build <br>
CMPredefinedTargets/INSTALL. The default installation directory on Windows is the build subdirectory<br>
 of the build folder. Installation build is necessary for further use of the library.

<a name="With_NMake"></a>
###### With NMake

From developer command prompt for the compiler call nmake in the build directory. This should build <br>
the library (and other targes if chosen, e.g. most examples are build on default). Call nmake install<br>
to install DiSCaMB (necessary for further use of the library).
 
<a name="Customizing_installation"></a>
### Customizing installation

Command line: make a directory where DiSCaMB will be created. <br/>
Call:  <br/>
&nbsp;&nbsp;cmake "path/to/DiSCaMB/root/dir" [options]<br/>

the most useful **command line cmake options**:<br/>

&nbsp;&nbsp;-G - generator, specifies development environment for the build (e.g. unix makefiles,<br/>
 eclipse, visual studio). Call Cmake --help to see list of available generators. If not specified <br/>
 "Unix Makefiles" will be used on Unix like systems and Visual Studio on Windows if available,<br/>
 otherwise "NMake Makefiles"   
 
&nbsp;&nbsp;-T - toolset (for generators which supports toolsets), e.g. -T "v140_clang_3_7" to use<br>
 clang toolset on Visual studio 14,...  
  
&nbsp;&nbsp;-Dvariable_name=x sets variable in cmake script variable_name to x, e.g. CMAKE_CXX_FLAGS<br>
will set C++ compiler flags for project with no multiple configuration types (e.g. "Unix Makefiles"<br>
generator)  
 
<a name="Project_specific_options_and_variables"></a>
#### project specific options and variables

BUILD_DOCS  - Build the doxygen documentation (default OFF)  <br/>
BUILD_EXAMPLE_REFINER - Build example refinement program (requires alglib) (default OFF)<br/>
BUILD_EXAMPLES - Build the example programs (default ON)  <br/>
BUILD_FOR_GPU - Build version for GPU (default OFF)<br/>
BUILD_SHARED_LIBS  - Build shared libraries (default OFF)  <br/>
BUILD_TESTS - Build the test programs (default ON)  <br/>
GENERATE_INSTALL - Generate installation target (default ON)  <br/>
USE_CPP11 - Compile components requireing C++11 (some examples), disabled for GPU build (default ON)<br/>
<br/>
INSTALLATION_PREFIX - prefix for instalation, default is /usr/local on Linux 
                      and build folder in build directory on Windows <br/>
ALGLIB_INCLUDE_PATH - path to alglib header folder (necessary when building example refinement program) <br/>
ALGLIB_LIBRARY - path to alglib library (necessary when building example refinement program) <br/>
ALGLIB_LIBRARY_DEBUG - path to debug version of alglib library debug version (if not specified ALGLIB_LIBRARY is used) <br/>
ALGLIB_LIBRARY_OPTIMIZED - path to non debug version of alglib library debug version (if not specified ALGLIB_LIBRARY is used) <br/>

<a name="Compiler_options"></a>
##### Compiler 
Compiler can be set on command line with -DCMAKE_CXX_COMPILER=your_compiler.<br> 
Compiler options:  <br>
 * can be added with -DCMAKE_CXX_FLAGS="your flags" and -DCMAKE_C_FLAGS="your flags" for single<br>
 configuration generator like Makefile generator 
 * for multiple configuration generators use options for specific configuration (see <br>
 https://cmake.org/Wiki/CMake_Useful_Variables) e.g. -DCMAKE_CXX_FLAGS_RELEASE="your flags" for <br>
 release configuration  

<a name="Installation_of_the_example_refine"></a>
### Installation of the example 'refine'

The example is not built on default. One can either built it when compiling DiSCaMB or later in separate<br>
step. 

<a name="Building_refine_during_DiSCaMB_compilation"></a>
#### Building refine during DiSCaMB compilation.

To do that a compiled alglib library is needed. For some systems prebuild alglib is available e.g. for <br>
Debian and Ubuntu. It can be installed in this case with 'sudo apt-get install name_of_package', where<br>
name_of_package is a particular version of alglib (the name can be find out by calling 'apt-cache search alglib').<br> 
  
Alternatively one can build alglib. To do that copy alglib source files from /cpp/src folder in its<br>
 distrubution and compile into library, e.g. with g++ it would be:
```
g++ -c *.cpp
ar rvs alglib.a *.o 
```

Having alglib library ready we can build the refine example when compiling DiSCaMB. It is not done by<br>
default. Option BUILD_EXAMPLE_REFINER has to be set (in command line it would be -DBUILD_EXAMPLE_REFINER=ON).<br>
Two other variables have to be specified to cmake:<br>
ALGLIB_INCLUDE_PATH - the folder with header files have to specified. If alglib was compiled by user<br>
it can be the folder where compilation was performed (with header and cpp files from /cpp/src folder<br>
of alglib distribution). If alglib was preinstalled it might be e.g. /usr/include/libalglib but in <br>
general it may depend on particular package installed, installation setting e.t.c..<br>
ALGLIB_LIBRARY - path to alglib library. 

<a name="Building_refine_in_separate_step"></a>
#### Building refine in separate step.

Once DiSCaMB is built copy content of DiSCaMB/examples/refine and constent of folder /cpp/src from <br>
alglib distribution into one folder and compile together pointing to DiSCaMB (see also getting started <br>
page of help, point ['Compiling the simple program - linking against DiSCaMB'](#linking)). <br>
With gcc system the compilation command would be in general something like that:<br>
```
g++ refine.cpp Target.cpp alglib_cpp_files -Idiscamb -fopenmp
```
where alglib_cpp_files are all *.cpp alglib files. <br>
Below is example with pointing to DiSCaMB files directly (should not be necessary if DiSCaMB is<br>
 installed):<br>
``` 
g++ alglibinternal.cpp alglibmisc.cpp ap.cpp dataanalysis.cpp diffequations.cpp fasttransforms.cpp
 integration.cpp interpolation.cpp linalg.cpp optimization.cpp refine.cpp solvers.cpp 
 specialfunctions.cpp statistics.cpp Target.cpp -I/mnt/f/tmp/11/include /mnt/f/tmp/11/lib/libdiscamb.a -fopenmp
``` 
/mnt/f/tmp/11/ is here installaction location for DiSCaMB (exmaple for Ubuntu on Windows 10).<br>
