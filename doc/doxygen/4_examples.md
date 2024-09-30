\page examples_page Summary of examples

## custom_model

Structure factor calculations for deformation valence like electron density term:
\f[
   \rho(\mathbf{r}) = R_1(r) * Y_{10}(\theta,\phi) + R_2(r) * Y_{20}(\theta,\phi)
\f]
   
where \f$ R_i \f$ are radial functions and \f$ Y_lm \f$ are density normalized spherical harmonics.<br>
The radial functions are defined as follows:
\f[   
   R_1(r) = a_1 * r * e^{-b_1*r} + a_2 * r^2 * e^{-b_2*r} + a_3 * r^3 * e^{-b_3*r}
\f]
\f[
   R_2(r) = a_4 * r^3 * e^{-b_4*r}
\f]
The program does not need arguments to run.

## refine

Extremely basic program for refinement. Usage:
```
refine xd_mas xd_par hkl format [options_file] [-nd]
```
where:
- xd_mas is name of XD master file
- xd_par is name of XD parameter file
- hkl is name of file with intensity data, no headers, first three columns are h,k, l indices, columns<br>
  with observed intensities and their standard deviation are specified with the next argument - format
- format - it can be 5 or 6, in the case of 5 - intensities are in 4-th column and standard deviations <br>
  in 5-th column, if the argument takes value 6 then intensities are in 5-th column and standard deviations<br> 
  in 5-th column
- options_file (optional argument) (to be documented, specifies defines which variables should be optimized, <br>
  by default all are optimized except of occupancies if they are equal to 1.0)
- option -nd which switch off distorting crystal structure before refinement (the distortion is made<br>
  to show that refine program really refines)  

## sf_ubdb

Calculates structure factor for few Miller indices. On input takes  2 arguments: name of xd master file<br>
and xd parameter file.

##  slater_type_atomic_wfn

Parameters of wavefunction for carbon atom taken from E. Clementi and C. Roetti \cite clementi_roetti_1974 are used to construct<br>
spherically averaged atomic electron density and calculate the corresponding form factor. The parameters are<br>
also obtained from tables included in the discamb and the corresponding form factor are calculated for comparison.<br>
No input arguments needed.

## type assignment info

Reads model from XD files and prints information on grouping atoms into types. Expected 2 or 3 arguments:<br>
 [all] name of xd master file and xd parameter file. Produces more information when called with the option 'all'.
