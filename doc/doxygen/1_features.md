\page features Main features of DiSCaMB  
## Features <a name="Features"></a>


### Main features

* Calculation of multipolar structure factors and their derivatives with respect to structural <br>
  parameters (position, occupancy and ADPs), various versions suitable for large and small<br>
  molecules. 
* Parallelized for use with multicore processors and graphics processing units. 
* Easy access to computational steps required to calculate the structure factors. 
* Algorithms aware of pseudoatom types.

####Detailed list

##### Structure factor calculations. 

* Access to atomic wavefunction data. 
* Convertion of the atomic wavefunction into corresponding spherically averaged density (linear <br>
  combination of Slaters ) for given electronic configuration.
* Calculation of the corresponding form factor. 
* Calculation of the whole deforamation density form factor for given Plms and radial function<br> parameters. 
* Calculation of individual terms in deformation density form factor. 
* Separate calculation of angular and radial parts of the form factor for each term (allows for<br>
  e.g. more complex radial functions). 
* Calculation of Fourier-Bessel transform of Slater type functions. 
* Local coordinate systems calculation.

##### Structural data handling. 
* Symmetry operations - multiplication, application to vectors, string and matrix notation. 
* Conversions between fractional and Cartesian coordinate system in direct and reciprocal space,<br>
 conversions of ADPs and parameters derivatives.

#####Utilities. 
* Math utilities - basic 3D algebra and calculations involving spherical harmonics.
* Utilities for string operations, error handling, performance measurement. 
* Reading Hansen-Coppens<br> model parameters from [lsdb](http://crystal.chem.uw.edu.pl/software.html) generated XD files. 
 
#####Project organization. 
* Build process managed with [CMake](http://www.cmake.org.) for enabling compilation environment of<br>
  user choice. 
* Code documentation including [Doxygen](www.doxygen.org) generated API documentation. 
* Extensive test set. 
* Examples illustrating usage.
