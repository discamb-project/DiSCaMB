\page xd Differences in DiSCaMB and XD parameterization

# Differences in DiSCaMB and XD parameterization

## Electron configuration

Electron configuration including valence orbitals assignment is the same as in XD2006
with the exception of Cu which has no valence orbitals in XD in Clementi-Roetti wavefunction bank 
(but it has in the other ones) and in DiSCaMB 4s orbital is default valence orbital

## Exponents and powers in deformation valence terms

For Zn and Cu DiSCaMB has by default 4 as power of r and XD has 6.
Exponent for Al3+ - in DiSCaMB 5.13817 in XD2006 5.1544 (which is \f$ \xi_{3s} + \xi_{3p} \f$ )

## Types

Types present in XD Clementi-Roetti bank but not in DiSCaMB:
Mn4+, Ni3+, Si4+, Ti2+ 

Types present in DiSCaMB Clementi-Roetti wave functions set but not in XD:
S-, Se-, V4+

## Atomic orbital parameters 

from Clementi and Roetti publication (see \cite Clementi_Roetti_1974 ):

For S atom one of coefficients in the publication pdf is given as stars (********)
We have taken value which gives correct normalization (two values are possible, 
the one with lower absolute value is chosen). In XD the coefficient is -18.5018, 
in DiSCaMB -18.50193.

Ar  exponent for 3-rd STO of 2P and 3P:
in DiSCaMB  12.39970
in XD       12.39910  

Fe - coefficient for 11-th STO of 1S orbital
XD      -0.00008 
DiSCaMB -0.00006

Ga - exponent for 5-th STO of 1S, 2S, 3S, 4S
XD      12.54240 
DiSCaMB 12.54424

As - coefficient for 2-nd STO of 4P
XD      0.00735
DiSCaMB 0.00736

Ti4+ - coefficient for 5-th STO of 1S
XD      0.00236
DiSCaMB 0.00238

Fe3+ - exponent for 6-th STO of 2P and 3P
XD      2.82300 
DiSCaMB 2.82800


