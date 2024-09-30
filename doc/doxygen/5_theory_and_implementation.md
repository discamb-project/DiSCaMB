\page theory Theory and Implemention

* [Electron density model](#ED_model)
    - [general form](#EC_general_form)
    - [core and valence spherical density](#ED_core_valence)
    - [deformation-valence density](#ED_def_val)
* [Structure factor](#Structure_factor)



# Hansen Coppens multipole model 
(For theoty see also original paper introducing the model \cite hansen_coppens_1978 and P. Coppens book \cite coppens1997).
## Electron density model <a name="ED_model"></a>

### General form <a name="EC_general_form"></a>

Atomic electron density in Hansen Coppens mulipole model is defined as follows:
<table border="0" cellpadding="20"><tr><td>
\f[
\rho_A(\mathbf{r}) =
         \rho_{core}(r) +  P_{val} \kappa^3 \rho_{val}(\kappa r) +
         \sum_l {\kappa '}^3 R_l( \kappa ' r) \sum_m P_{lm} d_{lm} (\theta_L , \phi_L)
\f]
</td><td> (hc.d.1) \anchor hcd1 </td></tr></table>

where \f$ \rho_{core} \f$ and \f$ \rho_{val} \f$ are core and valence (spherically averaged) <br>
densities of reference atom and the last term is so called deformation-valence density. The indices <br>
\f$ L \f$ in \f$ \theta_L , \phi_L \f$ in the last term indicates that values of the angular coordinates <br>
depends of local coordinate system for given atom.

In DiSCaMB parameters of atomic density are represented by three components:
- Parameters related to atomic wave function are represented by discamb::HC_WfnType. Such data are <br>
  commonly shared by multiple atoms in molecule (e.g. all carbon atoms). They include information on<br>
  atomic orbitals, assignment to core and valence orbitals, exponents ( \f$ \zeta \f$ ) and powers ( \f$ n_l \f$ ) <br>
  of r for radial functions of deformation valence term (see eqs. \ref hcd11 "hc.d.11-12") and anomalous <br>
  dispersion term.
- Parameters not related to atomic wave function including coefficients \f$ P_{val} \f$ , \f$ P_{l,m} \f$  <br>
  and parameters \f$ \kappa \f$ and \f$ \kappa ' \f$ . These data are represented by discamb::HC_AtomTypeParameters <br>
  These parameters are usually not shared by atoms in multipole model refinement except of atoms which<br>
  are assumed to be chemically equivalent. However they are shared by atoms of the same type <br>
  parameterized with help of banks of aspherical atom parameters used in Transferable Aspherical Atom <br>
  Model (TAAM).
- Local coordinate system represented as matrix type discamb::Matrix3d. A class for handling XD type<br>
  localc coordinate system is provided (discamb::XdLocalCoordinateSystem).  

Total atomic density is a sum of atomic contributions. Its parameterization related to Hansen-Coppens<br>
model is represented by structure discamb::HC_ModelParameters which contains sets of the above mentioned<br>
components of atomic density parameterization (except atomic local coordinate systems) of and mapping<br>
between atoms and the components.  

### Core and valence densities <a name="ED_core_valence"></a>

Core and valence densities are generated from the corresponding core and valence orbitals by
spherical averaging, e.g.: 
<table border="0" cellpadding="20"><tr><td>
\f[
\rho_{core}(\mathbf{r}) = \hat{S}_A \sum_{n,l,m} o_{n,l,m} \|\varphi_{n,l,m} (\mathbf{r}) \|^2
\f]
</td><td> (hc.d.2) </td></tr></table>
the summation runs over occupied core orbitals 
(indices \f$ n, l \f$ and \f$ m \f$ correspond to principal, azimuthal and magnetic <br>
 quantum numbers), 
 \f$ \hat{S}_A \f$ symbolizes spherical avraging \f$ o_{n,l,m} \f$ are occupancy factors
  and \f$ \varphi_{n,l,m} (\mathbf{r}) \f$ are atomic orbitals: 
<table border="0" cellpadding="20"><tr><td>
\f[
\varphi_{n,l,m} (\mathbf{r}) = \mathcal{R}_{n,l}(r) Y_{l,m}(\theta, \phi)
\f]
</td><td> (hc.d.3) </td></tr></table>
with radial function \f$ \mathcal{R}_{n,l}(r) \f$ and angular part given in terms of
('wavefuntion normalized') spherical harmonics \f$ Y_{l,m}(\theta, \phi) \f$. 

### Spherical averaging

The spherical averaging  is just a projection of the angular part of the averaged function (say 
\f$ f(r) h(\theta, \phi) \f$ onto 1D space <br>
 of spherically symmetric angular functions spanned by \f$ Y_{0,0} = (2\sqrt{\pi})^{-1} \f$ ) i.e.:
<table border="0" cellpadding="20"><tr><td>
\f[
\hat{S}_A f(r) h(\theta, \phi) = 
 f(r) Y_{0,0} \int \ Y_{0,0}(\theta, \phi) h(\theta, \phi) dS = 
 \frac{f(r)}{4\pi} \int h(\theta, \phi) dS
\f]
</td><td> (hc.d.3) </td></tr></table>
in the case of the square norms of atomic orbital \f$ \varphi_{n,l,m} \f$ the angular part is 
\f$ Y_{l,m}Y_{l,m}^{*} \f$ <br>
which integrates out to 1 due to orthonormality of the spherical harmonics. For example for 
\f$ \rho_{core}(\mathbf{r}) \f$ <br> 
this leads to the following expression:
<table border="0" cellpadding="20"><tr><td>
\f[
\rho_{core}(\mathbf{r}) = (4\pi)^{-1} \sum_{n,l} o_{n,l} \mathcal{R}_{n,l}^2(r)
\f]
</td><td> (hc.d.4) </td></tr></table>
the occupancy factors (already summed over \f$ m \f$ ) corresponds to electronic configuration and <br> 
the summation runs over the core orbitals.

### Radial functions of atomic orbitals and atomic densities.

The radial part \f$ \mathcal{R}_{n,l}(r) \f$ of atomic orbital \f$ \varphi_{n,l,m} \f$ is
 a linear combination of Slater type functions <br>
 (the indices \f$ l,m \f$ are temporarly skept for the sake of readability):
 <table border="0" cellpadding="20"><tr><td>
\f[
\mathcal{R}(r) = \sum_p c_p \chi_p(r) 
\f]
</td><td> (hc.d.5) \anchor hcd5 </td></tr></table>
 The Slater type functions (\f$ \chi_p(r) \f$) are defined as:
<table border="0" cellpadding="20"><tr><td>
\f[
\chi_p(r) = \mathcal{N} (\alpha_p,k_p) r^{k_p} e^{-\alpha_p r}
\f]
</td><td> (hc.d.6) </td></tr></table>
and the normalization factors \f$ \mathcal{N}(\alpha,n) \f$ are given by the following expression:
<table border="0" cellpadding="20"><tr><td>
\f[
\mathcal{N}(\alpha,k)  = \frac{(2\alpha)^{k+3/2}}{ \sqrt{[2(k+1)]!}}
\f]
</td><td> (hc.d.7) \anchor hcd7 </td></tr></table>
(implemented in discamb::sto_atomic_wfn::stoNormalizationFactor). Inserting explicitly  the radial<br>
 functions into expression for core density gives the following lenghty expression :
<table border="0" cellpadding="20"><tr><td>
\f[
\rho_{core}(\mathbf{r}) = (4\pi)^{-1} \sum_{n,l} o_{n,l} \sum_{p,q} 
                          \mathcal{N} (\alpha_{n,l,p},k_{n,l,p}) \mathcal{N} (\alpha_{n,l,q},k_{n,l,q})
                          r^{k_{n,l,p}+k_{n,l,q}} e^{-(\alpha_{n,l,p} + \alpha_{n,l,q}) r}
                        
\f]
</td><td> (hc.d.8) \anchor hcd8 </td></tr></table>
Which can be written in shorter form as:
<table border="0" cellpadding="20"><tr><td>
\f[
\rho_{core}(\mathbf{r}) = \sum_v d_v r^{n_v} e^{-\beta_v r}
\f]
</td><td> (hc.d.9) \anchor hcd9 </td></tr></table>


The information about orbital radial functions published by Clementi and Roetti \cite Clementi_Roetti_1974 can be<br>
accessed using class discamb::ClementiRoettiData which pass this information together with information<br>
on default orbital occupancy and assignement to core/valence orbitals via structure discamb::HC_WfnBankEntry.<br>
It can be used for calculation of the core/valence density parameters appearing in eq. hc.d.9 using<br>
 using discamb::sto_atomic_wfn::atomicStoWfnToSphericalDensity.

## Deformation valence term of atomic electron density <a name="ED_def_val"></a>

The deformation density \f$ \rho_{def-val}(\mathbf{r}) \f$ term in Hansen-Coppens model
is given by the following sum:

<table border="0" cellpadding="20"><tr><td>
\f[
\rho_{def-val}(\mathbf{r}) = \sum_l {\kappa '}^3 R_l( \kappa ' r)
                             \sum_{m=-l}^{l} P_{lm} d_{lm} (\theta , \phi)
\f]
</td><td> (hc.d.10) </td></tr></table>

where \f$ R_l(r) \f$ are radial functions, \f$ P_{lm} \f$ are coefficients in multipole expansion, <br>
\f$ \kappa ' \f$ is contraction-expansion factor and \f$ d_{lm}(\theta , \phi) \f$ are real valued
density normalized spherical harmonics. 

### Radial functions

The radial functions are density normalized slater type functions:

<table border="0" cellpadding="20"><tr><td>
\f[
R_l(r) = N(\zeta,n_l) r^{n_l} e^{-\kappa ' \zeta r}
\f]
</td><td> (hc.d.11) \anchor hcd11 </td></tr></table>
with the noramlization factor ( \f$ \int R_l(r) r^2 dr = 1 \f$ ) defined in the following way:
<table border="0" cellpadding="20"><tr><td>
\f[
N(\zeta,n) = \frac{\zeta^{n+3}}{(n+2)!}
\f]
</td><td> (hc.d.12) \anchor hcd12 </td></tr></table>
(implemented in discamb::sto_atomic_wfn::stoDensityNormalizationFactor). The parameters \f$ n_l \f$ and \f$ \zeta \f$<br>
can be obtained from discamb::DeformationValenceParameters.

### Spherical harmonics
Real spherical harmonics \f$ d_{lm} \f$ are identical to those [defined](https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form)
                                                                         and [tabularized](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics)
                      	    																 at Wikipedia<br>
except of normalization factors.
 
## Structure factor <a name="Structure_factor"></a>

Atomic electron density in Hansen-Coppens multipole model is a sum of multipolar terms  \f$ \rho_t(\mathbf{r}) \f$: <br>

<table border="0" cellpadding="20"><tr><td>
\f[
\rho_t(\mathbf{r}) = R(r)d_{lm}(\theta,\phi)
\f]
</td><td> (hc.s.1) \anchor hcs1 </td></tr></table>

scattering corresponding to such term is given by the following expression:
<table border="0" cellpadding="20"><tr><td>
\f[
\mathcal{F}[ \rho_t(\mathbf{r})](\mathbf{h}) =
              4 \pi i^l d_{lm}(\theta,\phi) \int_0^{\infty} R(r) j_l(2\pi h r) r^2 \mathrm{d}r

\f]
</td><td> (hc.s.2) \anchor hcs2 </td></tr></table>

where \f$ j_l \f$ is \f$l\f$-th order spherical Bessel function. We will refere to the integral factor<br>
as 'radial scattering' and to the factor \f$ 4 \pi i^l d_{lm}(\theta,\phi) \f$ as angular scattering. <br>
Replacing \f$d_{lm}\f$ by spherical harmonics in another convention would still lead to valid expression<br>
For Slater type radial functions the radial scattering functions \f$ g(l,n,h,z) \f$:
<table border="0" cellpadding="20"><tr><td>
\f[
g(l,n,h,z) = \int_0^{\infty} e^{-z r} r^n j_l(2\pi h r) \mathrm{d}r
\f]
</td><td> (hc.s.3) \anchor hcs3 </td></tr></table>
are tabularized, they can be also calculated recursively. Tabularized version is available in DiSCaMB<br>
in discamb::sto_scattering::gFunction(). Various variants of 'angular scattering' evaluation functions<br>
are available in discamb::multipole_scattering namespace. Spherical core and valence density scattering<br>
may be evaluated with function discamb::sto_scattering::scatteringSphericalDensity() which takes linear<br>
combination of Slater type functions (given by eg. \ref hcd9 "hc.d.9") as input. Evaluation of deformation<br>
valence density contribution to scattering is provided via discamb::multipole_scattering::hcDeformationValence()<br>
Scattering factor calculation for given crystallographic structure can be performed using dedicated<br>
class discamb::HansenCoppensStructureFactorCalculator. It allows also for scattering factor calculations<br>
for spherical atom model.



