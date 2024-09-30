\page algorithms_for_parameterization Algorithms for parameters assignment
   

## Default assignment to virtual orbitals

Default orbital assignent algorithm as core or valence (limited to periods 1-4):
- if nuber of electrons if one of the following - 0, 2, 10, 18, 36 - there are no valence electrons
- if cation with 28 electrons then there are no valence electrons
- if non of the above and s-block atom/ion - the outermost s orbitals are valence orbitals
- if non of the above and p-block atom/ion - the outermost s and p orbitals are valence orbitals
- if non of the above and d-block atom/ion with:
 -# 10 d electrons - no valence orbitals
 -# less than 10 d electrons - the corresponding d orbitals are valence orbitals


## Deformation valence parameters

Algorithms for assignment of parameters  
(\f$ n_l \f$ and \f$ \zeta \f$) in deformation-valence term
(\f$ \rho_{def-val}(\mathbf{r})) \f$ in Hansen-Coppens model:
\f[
\rho_{def-val}(\mathbf{r}) = \sum_l r^{n_l} e^{ - \zeta \kappa^{'} r} \sum_m P_{lm} Y_{lm} (\theta , \phi)
\f]

### Powers or r  - \f$ n_l \f$

Default powers are generated as follows:
- 1 - st period : powers[i] = i
- 2 - nd period : powers = { 2, 2, 2, 3, 4 }
- 3 - rd period elements and 4 period transition metals : powers[i] = 4
- 4 - rth period main group's elements: powers[i] = 6 
With the exception of Cu and Zn they are the same as in XD2006 for atoms types included in
DiSCaMB subset of Clementi-Roetti wave functions. 

### Exponent multiplier in radial function - \f$\zeta\f$

Default exponents (\f$\zeta\f$) are based on single zeta atomic orbitals exponents (\f$\xi_v\f$)
taken from Clementi and Raimondi publication \cite Clementi_Raimondi_1963.

If valence orbitals are present they are given by:
\f[
\zeta = 2 \frac{ \sum_{v} o_v \xi_v }{\sum_{v} o_v }
\f]
sum runs over occupied valence orbitals, \f$o_v\f$ are occupation factors 
(taken on defult from data Clementi and Roetti publication \cite Clementi_Roetti_1974, 
default assignment to valence orbitals is described above).

If there are no valence orbitals, then:
   - for S and D block atoms the exponent \f$ \zeta \f$ is equal to twice of the value of the exponent for 
     the outermost orbital of respectively s or d subshell of reference (neutral) atom
   - for P block atoms the exponent \f$ \zeta \f$ is equal to twice of the weighted sum of the
     exponents for 
     the outermost S and P subshells orbitals of reference (neutral) atom: 
     \f$ \zeta = 2 (2 \xi_s + 6 \xi_p )/8 \f$

