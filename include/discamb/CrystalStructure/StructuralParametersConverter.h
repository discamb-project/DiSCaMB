
#ifndef _DISCAMB_CRYSTALSTRUCTURE_STRUCTURALPARAMETERSCONVERTER_H_
#define _DISCAMB_CRYSTALSTRUCTURE_STRUCTURALPARAMETERSCONVERTER_H_

//#include "config.h"

#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include <complex>
#include <vector>


namespace discamb {

/**
* \addtogroup CrystalStructure
* @{
*/


namespace structural_parameters_convention
{
    /**Indicates atomic displacement parameters (ADPs) convention (see e.g.
    R.W. Grosse-Kunstleve and P. D. Adams, J. Appl. Cryst. (2002). 35, 477-480    
    "On the handling of atomic anisotropic displacement parameters")
    */
    enum class AdpConvention { 
        /** \f$U^*\f$ */
        U_star,
        /** Cartesian atomic displacement parameters.*/
        U_cart,
        /** As defined for _atom_site_aniso_U_11 type entries in CIF file. */
        U_cif };
    /** Indicates coordinate system for atomic positions. */
    enum class XyzCoordinateSystem { 
        /** Indicates fractional coordinates.*/ fractional, 
        /** Indicates Cartesian coordinates.*/ cartesian };
}

/**
Converts between fractional and Cartesian coordinates, 
various conventions of anisotropic displacement parameters and
between derivatives of scattering factor with respect to structural parameters 
in various coordinate systems / conventions.
It is assumed that the a-axis is collinear with the x-axis, and the c*-axis parallel to z.
*/

class StructuralParametersConverter
{
public:
    /** \brief Constructs converter for unit cell 1, 1, 1, 90, 90, 90 (it can be changed with set(const UnitCell &)).*/
    StructuralParametersConverter();
    /** \brief Constructs converter for the provided unit cell.*/
    StructuralParametersConverter(const UnitCell &unitCell);
    ~StructuralParametersConverter();
    
    /** \brief Sets the unit sell defining coordinate system.*/
    void set(const UnitCell &unitCell);

    /** \brief Converts fractional coordinates into Cartesian coordinates.*/
    void xyzFractionalToCartesian(const Vector3d &fractional, Vector3d &cartesian) const;
    /** \brief Converts Cartesian coordinates into fractional coordinates.*/
    void xyzCartesianToFractional(const Vector3d &cartesian, Vector3d &fractional) const;


    /** \brief 
    Returns a matrix \f$\mathbf{M}\f$ which transforms vector \f$\mathbf{r}_f\f$ in fractional coordinates
    into vector \f$\mathbf{r}_c\f$ in cartesian coordinates:  \f$\mathbf{r}_c = \mathbf{M} \mathbf{r}_f \f$.
    */
    Matrix3d xyzFractionalToCartesianMatrix() const;
    /** \brief
    Returns a matrix \f$\mathbf{M}\f$ which transforms vector \f$\mathbf{r}_f\f$ in fractional coordinates
    into vector \f$\mathbf{r}_c\f$ in cartesian coordinates:  \f$\mathbf{r}_c = \mathbf{M} \mathbf{r}_f \f$.
    */

    Matrix3d xyzCartesianToFractionalMatrix() const;


    /** \brief
    Returns a matrix \f$\mathbf{M}\f$ which transforms vector of ADP components \f$\mathbf{u}_1\f$
    into to vector expressed in other convention \f$\mathbf{u}_2\f$: 
     \f$\mathbf{u}_2 = \mathbf{M} \mathbf{u}_1 \f$.
    */
    void getAdpConversionMatrix(structural_parameters_convention::AdpConvention adpConventionIn,
                                structural_parameters_convention::AdpConvention adpConventionOut,
                                std::vector<std::vector<double> > &m) const; 


    /** \brief Converts positional coordinates from system defined by coordinateSystemIn to system defined by coordinateSystemOut.*/
    void convertXyz(const Vector3d &xyzIn, Vector3d &xyzOut, 
                    structural_parameters_convention::XyzCoordinateSystem coordinateSystemIn,
                    structural_parameters_convention::XyzCoordinateSystem coordinateSystemOut) const;

    /** \brief Converts atomic displacement parameters in adpIn from parameterization defined by adpConventionIn to parameterization defined by
    adpConventionOut and saves them to adpOut.*/
    void convertADP(const std::vector<double> &adpIn, std::vector<double> &adpOut, 
                    structural_parameters_convention::AdpConvention adpConventionIn,
                    structural_parameters_convention::AdpConvention adpConventionOut) const;

    /** \brief Converts positional coordinates from system defined by coordinateSystemIn to system defined by coordinateSystemOut.*/
    void convertDerivativesXyz(const Vector3<std::complex<double> > &derivativesIn,
                               Vector3<std::complex<double> > &derivativesOut,
                               structural_parameters_convention::XyzCoordinateSystem coordinateSystemIn,
                               structural_parameters_convention::XyzCoordinateSystem coordinateSystemOut) const;

    /** \brief Converts ADPs from convention defined by adpConventionIn to the one defined by adpConventionOut.*/
    void convertDerivativesADP(const std::vector<std::complex<double> > &dAdpIn,
                               std::vector<std::complex<double> > &dAdpOut,
                               structural_parameters_convention::AdpConvention adpConventionIn,
                               structural_parameters_convention::AdpConvention adpConventionOut) const;

    /** \brief Converts \f$U_{cif}\f$ to \f$U^*\f$ */
    void uCifToU_star(const std::vector<double> &uCif, std::vector<double> &uStar) const;

    /** \brief Converts \f$U^*\f$ to \f$U_{cif}\f$ */
    void uStarToU_cif(const std::vector<double> &uStar, std::vector<double> &uCif) const;

    /** \brief Converts \f$U_{cif}\f$ to \f$U_{cart}\f$ */
    void uCifToU_cart(const std::vector<double> &uCif, std::vector<double> &uCart) const;

    /** \brief Converts \f$U_{cart}\f$ to \f$U_{cif}\f$ */
    void uCartToU_cif(const std::vector<double> &uCart, std::vector<double> &uCif) const;

    /** \brief Converts \f$U_{cart}\f$ to \f$U^*\f$ */
    void uCartToU_star(const std::vector<double> &uCart, std::vector<double> &uStar) const;

    /** \brief Converts \f$U^*\f$ to \f$U_{cart}\f$ */
    void uStarToU_cart(const std::vector<double> &uStar, std::vector<double> &uCart) const;

    
    /** \brief Converts derivatives of structure factor with respect to atomic coordinates in fractional coordinates to
    derivatives w.r.t. atomic coordinates in Cartesian coordinates.

    \f$ [ \frac{ \partial F(\mathbf{h}) } { \partial x_{frac} },
          \frac{ \partial F(\mathbf{h}) } { \partial y_{frac} },
          \frac{ \partial F(\mathbf{h}) } { \partial z_{frac} } ] \rightarrow
        [ \frac{ \partial F(\mathbf{h}) } { \partial x_{cart} },
          \frac{ \partial F(\mathbf{h}) } { \partial y_{cart} },
          \frac{ \partial F(\mathbf{h}) } { \partial z_{cart} } ]
    
    \f$ */
    void derivativesXyzFractionalToCartesian(const Vector3<std::complex<double> > &derivativesWrtXyzFractional,
                                  Vector3<std::complex<double> > &derivativesWrtXyzCartesian) const;

    /** \brief Converts derivatives of structure factor with respect to atomic coordinates in Cartesian coordinates to
    derivatives w.r.t. atomic coordinates in fractional coordinates.

    \f$ [ \frac{ \partial F(\mathbf{h}) } { \partial x_{cart} },
    \frac{ \partial F(\mathbf{h}) } { \partial y_{cart} },
    \frac{ \partial F(\mathbf{h}) } { \partial z_{cart} } ] \rightarrow
    [ \frac{ \partial F(\mathbf{h}) } { \partial x_{frac} },
    \frac{ \partial F(\mathbf{h}) } { \partial y_{frac} },
    \frac{ \partial F(\mathbf{h}) } { \partial z_{frac} } ]

    \f$ */

    void derivativesXyzCartesianToFractional(const Vector3<std::complex<double> > &derivativesWrtXyzCartesian,
                                  Vector3<std::complex<double> > &derivativesWrtXyzFractional) const;

    /** \brief Converts derivatives of structure factor with respect to ADPs in CIF convention to
    derivatives w.r.t. ADPs in \f$U^*\f$ convention.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cif})_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial U^*_{ij} } \}
    \f$ */


    void derivativesU_CifToU_Star(const std::vector<std::complex<double> > &derivativesWrtU_Cif,
                      std::vector<std::complex<double> > &derivativesWrtU_Star) const;

    /** \brief Converts derivatives of structure factor with respect to ADPs in CIF convention to
    derivatives w.r.t. Cartesian ADPs.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cif})_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cart})_{ij} } \}
    \f$ */

    void derivativesU_CifToU_Cartesian(const std::vector<std::complex<double> > &derivativeWrtU_Cif,
                           std::vector<std::complex<double> > &derivativeWrtU_Cartesian) const;

    /** \brief Converts derivatives of structure factor with respect to \f$U^*\f$ ADPs to
    derivatives w.r.t. ADPs in CIF convention.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial U^*_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cif})_{ij} } \}
    \f$ */

    void derivativesU_StarToU_Cif(const std::vector<std::complex<double> > &derivativeWrtU_Star,
                      std::vector<std::complex<double> > &derivativeWrtU_Cif) const;
    /** \brief Converts derivatives of structure factor with respect to \f$U^*\f$ ADPs to
    derivatives w.r.t. Cartesian ADPs.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial U^*_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cart})_{ij} } \}
    \f$ */


    void derivativesU_StarToU_Cartesian(const std::vector<std::complex<double> > &derivativeWrtU_Star,
                            std::vector<std::complex<double> > &derivativeWrtU_Cartesian) const;
    
    /** \brief Converts derivatives of structure factor with respect to Cartesian ADPs to
    derivatives w.r.t. ADPs in \f$U^*\f$ convention.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cart})_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial U^*_{ij} } \}
    \f$ */

    void derivativesU_CartesianToU_Star(const std::vector<std::complex<double> > &derivativeWrtU_Cartesian,
                            std::vector<std::complex<double> > &derivativeWrtU_Star) const;
    
    /** \brief Converts derivatives of structure factor with respect to Cartesian ADPs to 
    derivatives w.r.t. ADPs in CIF convention.

    \f$ \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cart})_{ij} } \} \rightarrow
    \{ \frac{ \partial F(\mathbf{h}) } { \partial (U_{cif})_{ij} } \}
    \f$ */

    void derivativesU_CartesianToU_Cif(const std::vector<std::complex<double> > &derivativeWrtU_Cartesian,
                           std::vector<std::complex<double> > &derivativeWrtU_Cif) const;
private:
    Matrix3d mD_XyzCartToFrac, mD_XyzFracToCart, mD_UcartToUcif, mD_UcartToUcifTransposed,
             mD_UcartToUstar, mD_UcartToUstarTransposed, mD_UstarToUcart, mD_UstarToUcartTransposed, 
             mD_UstarToUcif, mD_UstarToUcifTransposed, mD_UcifToUcart, mD_UcifToUcartTransposed,
             mD_UcifToUstar, mD_UcifToUstarTransposed,
             mFracToCart, mCartToFrac,
             mCifToCart,mCifToCartT,mCartToCif,mCartToCifT,
             mCartToStar,mCartToStarT,mStartToCart,mStarToCartT;
    double mRS_Edges[3]; // lenths of edges of reciprocal space unit cell
    double mRS_EdgesInv[3]; // inversions of lenths of edges of reciprocal space unit cell
    double mCifToStar[6],mStarToCif[6];
    mutable Matrix3<std::complex<double> > mDerivativeWrtAdpIn, mDerivativeWrtAdpOut;
    mutable Matrix3d mAdpIn,mAdpOut;

    void adpDerivativeTransform(const std::vector<std::complex<double> > &adpIn,
                                std::vector<std::complex<double> > &adpOut, 
                                const Matrix3d & trnasformMatrix,
                                const Matrix3d & trnasformMatrixTransposed) const;

    void adpTransform(const std::vector<double> &adpIn,
                      std::vector<double> &adpOut,
                      const Matrix3d & trnasformMatrix,
                      const Matrix3d & trnasformMatrixTransposed) const;

};

/**@}*/

} //namespace discamb

#endif /*_DISCAMB_CRYSTALSTRUCTURE_STRUCTURALPARAMETERSCONVERTER_H_*/


