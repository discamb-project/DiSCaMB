\page getting_started Getting started 

In this section a simple program which uses DiSCaMB is introduced and compiled. Then another simple <br>
program for structure factor calculations is shown and finally a longer one is discussed fragment by<br>
fragment in order to overview DiSCaMB usage.

* [The first simple program](#simple_program_structural_data)
* [Compiling the simple program - linking against DiSCaMB](#linking)
* [Another simple program - structure factors calculation](#simple_program_structure_factors)
* [Longer program - overview of DiSCaMB functionalities](#longer_program)
    - [Crystal structure](#crystal_structure)
        + [defining atoms](#atoms)
        + [unit cell and coordinate transformations](#unit_cell_and_transform)
        + [handling symmetry](#handling_symmetry)
        + [defining crystal](#defining_crystal)
    - [Multipolar model](#multipolar_model)
        + [wave function related data](#wave_function_data)
        + [atom type specific data](#atom_type_data)
        + [putting it together](#multipolar_together)
    - [Structure factor calculation](#structure_factors_calculation)
        + [large molecules](#large_molecules)
        + [small molecules](#small_molecules)
    - [The code for the example](#long)        


## Simple program <a name="simple_program_structural_data"></a>
For very introductory example lets use the simple 'toy' program shown below (the next section shows <br>
how to build it):<br>

\includelineno example.cpp

It converts position in fractional coordinates into Cartesian ones and applies symmetry operation.<br>
 
The first three lines include discamb header files necessary for this example and standard library<br>
header (iostream) used here for printing out to console (see the next paragraph for the location of the <br>
header files. Discamb uses discamb namespace (thanks to line 5 there is no need to put the namespace<br> 
name in front of each DiSCaMB's function and structure). Vector3d is commonly used by discamb for handling<br>
real valued 3D vectors. Two such vectors are declared in line 10, one for fractional coordinates (r_f)<br>
and one for Cartesian ones. For the fractional to Cartesian coordinates transformation a unit cell <br>
has to be specified (line 12). Some value of the fractional coordinaes vector is assigned at line 13<br>
and then the vector is transformed into Cartesian one (line 17) and printed out (line 19).<br>

We also apply symmetry operation (defined at line 23) to the vector in fractional coorinates. It is <br>
transformed (line 26) to another symmetry equivalent vector (a vector transformed_r_f defined in line 24.<br>
Similarily as before the vector is transformed to Cartesian coordinates and then printed out.<br>

Let's see how to compile the example:<br>

<a name="linking"></a>
## Compiling the simple program - linking against DiSCaMB 

In order to compile the above example the compiler should be able to:
- find folder with discamb header files
- link the DiSCaMB library to the program

It can be as simple as calling<br>
```
g++ example.cpp -ldiscamb -std=c++0x -fopenmp
```
for default discamb installation on Linux, for GPU supporting version: 
```
g++ example.cpp -ldiscamb -ldiscamb_cuda -std=c++0x -fopenmp
``` 
The flag -std=c++0x might be required if C++11 standard features were used in example.cpp and -fopenmp<br>
DiSCaMB was compiled with support for openMP (which is the default).


In order to find folder with discamb header files discamb has to be installed i.e. target install has<br>
to be built. It can be done by calling 'make install' in build directory for typical linux command <br>
line setup, 'nmake install' in developer command prompt for Visual Studio command line or by building<br>
target install from GUI based code editors (e.g. Visual Studio, Eclipse itd). The headers folder will <br>
be located at installation_prefix/include. Here installation_prefix is the installation prefix used <br>
in discamb build. By default it is /usr/local/ on Linux and /build/include subfolder of build directory <br>
on Windows. If the headers folder is installed into directory which is searched by default by the <br>
 compiler then nothing more has to be done for the task. It can be the case for e.g. default linux <br>
 installation, which usually is searched by g++ compiler.  If the folder is not on compiler search path then <br>
the compiler has to be informed where the header files are (where is folder discamb with the header files). <br>
Let's assume that they are at path some_path. In the case of g++ compiler option -Isome_path should <br>
be used, when calling Microsoft c++ compiler from command line it would be /I some_path (see here for <br>
[more information for Visual Studio](https://msdn.microsoft.com/en-us/library/73f9s62w.aspx )).


The library is in installation_prefix/lib. Its name depends on operating system and on build type <br>
(debug or optimized). Optimized version name is libdiscamb.a on linux (and similar systems and environments)<br>
 and discamb.lib when compiled with Visual Studio on Windows. Debug version is libdiscambd.a or discambd.lib.<br>
To link against discamb library it is usually enough to specify te library when compiling from command line.<br>
Similarily as in the case of header files whole path is neee if the library containing is not at compiler<br>
search path for libraries.<br>

Assuming that the exmaple file is example.cpp and the compilation is done in the folder file, compilation<br>
command usually looks like that:<br>
```
compiler_name example.cpp optional_option_for_header_files_search path_to_discamb_lib
```
e.g. for Visual Studio  compiler called from command line:
```
cl example.cpp /I installation_prefix\include installation_prefix\lib\discamb.lib /MD /EHsc
```
/MD specifies version of compiler runtime library which works in this case (one of those /MD /MDd /MT /MTd <br>
should work ). E.g. if build directory is C:\\build, the command for default installation directory<br>
on Windows (in this case C:\\build\\build) would be:<br>

```
cl example.cpp /I C\build\build\include C\build\build\lib\discamb.lib /MD /EHsc
```

GUI based tools usually have its own specific means for defining external libraries and include <br>
directories so they are not discussed here.

<a name="simple_program_structure_factors"></a>
##  Another simple program - structure factors calculation

The next example shows how to calculate structure factors with XD files (xd.mas and xd.inp) used as <br>
input. DiSCaMB support for XD files is basic (files generated with lsdb based on UBDB are supported).<br>
This is the code for sf_ubdb example which is build on default:<br> 

\includelineno sf_ubdb.cpp

In the first line tools for reading XD files are included, in the second one tools for exception <br>
handling and in the third line structure factor calculation engine. 

Line 10 starts definition of function 'run' which is called from 'main' function. The main function <br>
contains only functionality for checking if number of input arguments is OK and for handling errors.<br>


HC_ModelParameters (line 12) is a container for hansen-Coppens multipolar model parameters. Crystal <br>
(line 13) is a structure for storing 'structural' data - unit cell and space group specification and <br>
atomic parameters including labels, positions ADPs and occupancy. Hkl indices are represented with<br>
3D integer vector structure Vector3i, local coordinate systems in XD like notation are represented by<br>
XdLocalCoordinateSystem class (line 16) but the structure factor calculation engine takes 3D matrices <br>
as the for locl coordinate specification (line 17). So localCoordinateSystems are used both to store<br>
information on local coordinates definition and to calculate it as 3D matrices for given positions of<br>
 atoms in crystal. <br>
 
 In the next lines (20-22) the Miller indices for which structure factors will be calculated are specified.<br>
Then multipolar model, crystal structure information and information about local coordinate systems<br>
 is read from XD files (line 25). <br>
 
 Next (lines 28-31) the local coordinate systems as 3D matrices are calculated for current geometry. <br>
 
 At line 35 engine for structure factor calculations is defined. Its constructor requires crystal structure<br>
 and multipolar model parameters as input argumnets.<br>  

The engine calculates structure factors (line 39) and finally they are printed out (l. 42-44).<br>

The example can be built in the same way as the previous one.<br>

<a name="longer_program"></a> 
## Longer program - overview of DiSCaMB functionalities

[Structure factor calculation](#structure_factor) require information on [crystal structure](#crystal_structure), parameterization<br>
 of [multipolar model](#multipolar_model) and [hkl vectors](#hkl_vectors) related information. <br>

 Various DiSCaMB functionalities are ilustrated with fragments of code, they can be combined into one<br>
 (the whole [source code](#long) is shown on the bottom of the page). The example uses urea structure <br>
 (XD files with UBDB parameterization for urea can be find in folder DiSCaMB/examples/data/urea). It<br>
 shows steps to be done for defining urea structure and multipolar model parameters manually (without <br>
 just setting it by read from XD files) and it shows how to calculate structure factors and related <br>
 derivatives in the case of both large and small molecule type calculations.<br> 
 
 
<a name="crystal_structure"></a> 
### Crystal structure 

Structural information is represented by a simple structure discamb::Crystal which contains information on:<br>
- atoms in crystal
- space group
- unit cell, and
- coordinate system conventions for atomic parameters (specifying e.g. that ADPs are represented by<br>
  \f$ U_{cif} \f$  or \f$U^*\f$)
  
<a name="atoms"></a>   
###### Defining atoms  
  
As an example we will use urea structure. Atoms are represented by discamb::AtomInCrystal structure, e.g. C <br>
atom in urea can be defined as follows: 

\includelineno atom_in_crystal.cpp

Coordinates are represented with discamb::Vector3d class (specialization of discamb::Vector3 type for type double). The<br>
coordinates can be either fractional or cartesian - discamb::Crystal stores information on the coordinate<br>
system in member xyzCoordinateSystem. Atomic displacement parameters are stored in container of type<br>
std::vector<double> and type of ADPs is deducted from size of the container. discamb::Crystal has a member<br>
adpConvention for indicating if the ADPs are Cartesian, \f$ U_{cif} \f$  or \f$U^*\f$. 

<a name="unit_cell_and_transform"></a>   
###### Unit cell and coordinate transformations

Unit cell is represented by class discamb::UnitCell which besides storing unit cell parameters allows<br>
for coordinate transformation between Cartesian and fractional coordinate systems. E.g. for Cartesian<br>
 coordinates calculation of C atom in urea can be done as follows:
\code{.cpp}
    discamb::UnitCell unitCell(5.565, 5.565, 4.684, 90, 90, 90);
    discamb::Vector3d cartesian;
    unitCell.fractionalToCartesian(c1.coordinates, cartesian);
    std::cout << c1.type << " " << cartesian[0] << " " << cartesian[1] << " " << cartesian[2] << std::endl;
\endcode

<a name="handling_symmetry"></a> 

###### Handling symmetry
 
Space group is defined by list of symmetry operations, either list similar to the one in CIF file or<br>
more typical to refinement programs list without operations generated by cell centering and inversion<br>
center (see also discamb::SpaceGroup). 

Space group operations (see also discamb::SpaceGroup) can be defined using either matrix-vector notation<br>
or strings as e.g. "x,y,z". For example:

\includelineno space_group.cpp

<a name="defining_crystal"></a> 
###### Defining crystal

Finally crystal structure for urea can be defined:
\code{.cpp}

    discamb::Crystal crystal;

    crystal.adpConvention = structural_parameters_convention::U_cif;
    crystal.atoms = { c1, o1, n1, h1, h2 }; // we have omited above declarations for atoms other than c1
    crystal.spaceGroup = p4bar21m;
    crystal.unitCell = unitCell;
    crystal.xyzCoordinateSystem = discamb::structural_parameters_convention::fractional;
        
\endcode

<a name="multipolar_model"></a>
### Multipolar model 

We will define multipolar model for urea, showing only parameterization for C atom. Multipole model<br>
parameters are divided into wave function related and atom type related. 

<a name="wave_function_data"></a>

###### Wave function related data

We will start with setting wave function reated data.

\includelineno multipole_model_wfn.cpp 

Class discamb::ClementiRoettiData was used for providing atomic wavefunction and default occupancy<br>
while discamb::DeformationValenceParameters for setting deformation valence term parameters \f$ n_l \f$ and \f$ \zeta \f$<br>

<a name="atom_type_data"></a>
###### Atom type specific data

Another part of model are type specific parameters - corresponding to parameters of types in banks of<br>
aspherical atomic densities. They include kappa parameters (\f$ \kappa \f$ and \f$ \kappa' \f$) and population parameters (\f$ P_{val} \f$ and \f$ P_{lm}\f$)

\includelineno multipole_model_type_data.cpp 
Definition of local coordinate system is also usually an element of atom type parameterization. In<br>
DiSCaMB it is however stored separately for greater flexibily of the coordinate system specification.<br>
Its usage is shown later.

<a name="multipolar_together"></a>
###### Putting it together

Finally we can put these to parts together into multipolar model parameterization:

\includelineno multipole_model_put_together.cpp 

Besides list of possible wave function related types and atom types a mapping between atoms in structure<br>
and the types should be set via members atom_to_wfn_map and atom_to_type_map of class HC_ModelParameters.<br>
The maps define index of wavefunction type in member wfn_parameters (and similarily for atom type)<br>
corresponding to the atom (indexaion starts from 0).

<a name="structure_factor">
##Structure factor calculation </a>

Having crystal structure and multipole model parameterization defined we can to calculate structure<br>
factors. discamb::HansenCoppensStructureFactorCalculator is designed for this task. We will also<br>
calculate derivatives. Here we have two options (1) derivatives of target function with respect to (w.r.t.)<br>
structural parameters and (2) derivatives of structure factor w.r.t. structural parameters. We will<br>
analyse the separately. Besides the data already defined the calculations will require infomation on<br>
local coordinate systems and Miller indices (hkl vectors). 

\includelineno lcs_and_hkl.cpp

<a name="large_molecules"> 
#### Macromolecular systems.

In this case derivatives of target function \f$ T \f$ w.r.t. structural parameters are more commonly used.<br>
For structure factor \f$ F_k = A_k + iB_k \f$ information on the following quantity: \f$ \frac{\partial T}{\partial A_k} + i \frac{\partial T}{\partial B_k}\f$ <br>
is needed (it let obtain the desired derivatives using chain rule). The code for such calculations <br>
can looks lke the one below:

\includelineno sf_macro.cpp

<a name="small_molecules"> 
#### 'Small molecules'.

In the case of 'small molecules' information on derivatives of structure factor w.r.t. structural<br>
parameters is calculated. Since it is calculated for each parameter and reflection there can be a large<br>
amount of data generated. Instead gathering all the data in one container the results are returned<br>
for each reflection. Implementation dedicated for calculations for one user provided hkl vector at time<br>
 is not provided (yet). Instead user have to provide and object of wiht the following member function:
\code{.cpp}
   void onPerHklCalculation(
           size_t hkl_idx, 
           std::complex<double> &structureFactor, 
           discamb::SfDerivativesAtHkl &derivatives);
\endcode
The results are passed for each hkl vectors as structure factor value and container with the derivatives<br>
at given hkl. An order of calculations does not have to correspond to the order of hkl vectors passed<br>
as an input argument (performing caculations for hkl indices ordered in specific way may result in<br>
more efficient code, the feature is not implemented yet). An example code for calculations for small<br>
molecules can looks as follows:

\includelineno sf_small.cpp 

This implementation of the object for collecting results is voversimplified, more realistic version <br>
would probably update some matrix for further use in least squares procedure. 

<a name="long"></a>

### The code for the example

The whole source code for the example is given below (it uses C++11):

\includelineno example2.cpp
