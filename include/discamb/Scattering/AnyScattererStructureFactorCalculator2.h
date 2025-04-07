#ifndef _DISCAMB_SCATTERING_AnyScattererStructureFactorCalculator22
#define _DISCAMB_SCATTERING_AnyScattererStructureFactorCalculator22


#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/Scattering/SF_Engine_DataTypes.h"
#include "discamb/Scattering/AtomicFormFactorCalculationsManager.h"
#include "discamb/Scattering/SF_CalcDataTypes.h"

#include <memory>


namespace discamb{

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class AnyScattererStructureFactorCalculator2
    {
        AnyScattererStructureFactorCalculator2();
        Crystal mCrystal;
        ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
        StructuralParametersConverter mConverter;
		std::shared_ptr<AtomicFormFactorCalculationsManager> mManager;
        std::vector<sf_engine_data_types::SymmetryOperation> mSymOps;
        void makeCartesianHkl(const std::vector<Vector3i> &hklIndices,std::vector<Vector3d> &hklCartesian);
        std::vector<double> mAtomicMultiplicityWeight;
        bool mMultiManager=false;

        std::vector<Vector3d> mR_at;
        //[atomIdx][symmIdx]
        std::vector< std::vector<Vector3d> > mR_at_symm;
        std::vector<std::vector<double> > mAdp;
        std::vector<std::complex<double> > mAnomalous;
		/*a buffer for storing atomic form factors [atomIdx][symmetry operation index]*/
		std::vector< std::vector<std::complex<double> > > mAtomsFormFactors;
		/*a buffer for storing atomic form factors [atomIdx]*/
		std::vector<std::complex<double> > mAtomsFormFactor;
        //vector<REAL> mOccupancies;
        //vector<REAL> atomic_multiplicity_weight = _atomic_multiplicity_weight;
        void pre_atom_loop_sf_calc(
            //in:
            const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
            const Vector3<REAL> &hVector,
            //out:
            std::vector<Vector3<REAL> > &rotated_h,
            std::vector<REAL> &translation_factor,
            std::vector<std::vector<REAL> > &adp_multipliers);
        
        void calculate_line_phase_factors(
            //in:
            const std::vector<Vector3d> &hkl_cart,
            int hkl_idx,
            const std::vector<Vector3i>& hkl_line,
            int idx_in_line,
            int line_direction,
            const std::vector<Vector3<double> >& rotated_h,
            const std::vector< std::vector<std::complex<double > > >& phase_factor_multiplier,
            //std::vector<Vector3d> &rotated_h_2pi,
            //std::vector<double> &translation_factor_2pi,

            //out:
            std::vector< std::vector<std::complex<double > > >& line_phase_factor);

        void calculate_line_temperature_factors(
            //in:
            const std::vector<Vector3d>& hkl_cart,
            int hkl_idx,
            const std::vector<Vector3i>& hkl_line,
            int idx_in_line,
            int line_direction,
            const std::vector<Vector3<double> >& rotated_h,
            std::vector< std::vector<double> >& temp_factor_multiplier,
            const std::vector< std::vector<double> >& temp_factor_multiplier_multiplier,
            std::vector<std::vector<double> > &adpMultipliers,
            //std::vector<Vector3d> &rotated_h_2pi,
            //std::vector<double> &translation_factor_2pi,

            //out:
            std::vector< std::vector<double> >& line_temperature_factors);


        void updateStructuralData(const std::vector<AtomInCrystal> &atoms);
        structural_parameters_convention::AdpConvention mAdpConvention = structural_parameters_convention::AdpConvention::U_cif;
        structural_parameters_convention::XyzCoordinateSystem mXyzConvention = structural_parameters_convention::XyzCoordinateSystem::fractional;

        void convertDerivatives(std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam) const;
        void convertDerivatives(discamb::SfDerivativesAtHkl &derivatives) const;

    public:
        static int findPreferredHklOrderingDirection(
            const std::vector<Vector3i>& hkl, 
            std::vector<std::vector<Vector3i> > &orderedHklLines,
            std::vector<std::vector<int> > &mapToOriginalSetIndices);
        /** Calculate structure factors and their derivatives*/
        AnyScattererStructureFactorCalculator2(const Crystal &crystal);

        ~AnyScattererStructureFactorCalculator2();

        void setAtomicFormfactorManager(std::shared_ptr<AtomicFormFactorCalculationsManager> &manager);
        /** Sets model of the electron density (note that the model is also set in constructor).*/
        void update(const std::vector<AtomInCrystal> &atoms);

        void setAnoumalous(const std::vector<std::complex<double> > &anoumalous);

        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const {};
        /**
        Sets convention for derivatives calculation:
        \param dXyzConvention specifies if structure factor derivatives with respect to coordinates should be calculated for
        fractional or cartesian coordinates.
        \param dAdpConvention specifies if structure factor derivatives with respect to aniosotropic ADPs should be
        calculated for \f$U_{cart}\f$, \f$U_{cif}\f$ or \f$U^{*}\f$

        */

        void setParametersConvention(structural_parameters_convention::XyzCoordinateSystem dXyzConvention,
            structural_parameters_convention::AdpConvention dAdpConvention);

        /**
        Gets convention for derivatives calculation.
        \sa setDerivativesConvention
        */
        void getParametersConvention(structural_parameters_convention::XyzCoordinateSystem &dXyzConvention,
            structural_parameters_convention::AdpConvention &dAdpConvention) const;



        /** \brief Calculates structure factors.

        \param atoms carry atomic structural parameters
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param f the resulting structure factors
        \param countAtomContribution flags for counting atom contribution to structure factor calculations, size should match
        number of atoms
        */
        void calculateStructureFactors(
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution);

        // parallel

        void calculateStructureFactors2(
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution);

        void calculateStructureFactors_iso_aniso_split(
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution);

        void calculateStructureFactors_ordered_hkl_alg(
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution);


        /** \brief Calculates structure factors.

        \param atoms carry atomic structural parameters
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param f the resulting structure factors
        */

        void calculateStructureFactors(
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f);

        /** \brief Calculates structure factors and derivatives of target function with respet to structural parameters.

        \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param f the resulting structure factors
        \param dTarget_dparam derivatives of target function with respet to structural parameters
        \param dTarget_df derivatives of target function \f$ T \f$ with respect to structure factors \f$ F_k = A_k + iB_k \f$,
        \f$A_k, B_k \in \mathbb{R}\f$ given as:
        \f$ \frac{\partial T}{\partial A_k} + i \frac{\partial T}{\partial B_k}\f$

        \param countAtomContribution flags for counting atom contribution to structure factor calculations, size should match
        number of atoms
        */
        void calculateStructureFactorsAndDerivatives(
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution);


        void calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution);

        /**
        \brief Calculates structure factors and their derivatives with respect to structural parameters.

        For each hkl vector an user provided object \p singleHklResultsProcessor function is called
        which collects the results for given h vector.

        \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param singleHklResultsProcessor object of user defined type intended for collecting results,
        for calculation for each Miller index its member function:
        \code{.cpp}
        void onPerHklCalculation(
        int hkl_idx,
        std::complex<double> &structureFactor,
        discamb::SfDerivativesAtHkl &derivatives);
        \endcode
        is called (hkl_idx starts from 0). See e.g. Tests/correctness/test_1/GradientDataCollector.h for an example implementation.

        */

        template<typename T>
        void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Matrix3d> &localCoordinateSystems,
            const std::vector<Vector3i> &hkl,
            T &singleHklResultsProcessor) const;

        //void calculateStructureFactorsAndDerivatives(
          //  const std::vector<AtomInCrystal> &atoms,
          //  const std::vector<Matrix3d> &localCoordinateSystems,
          //  const std::vector<Vector3i> &hkl,
          //  T &singleHklResultsProcessor) const;


        /**
        \brief Calculates structure factors and their derivatives with respect to structural parameters.


        \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param singleHklResultsProcessor object of user defined type intended for collecting results,
        for calculation for each Miller index its member function:
        \code{.cpp}
        void onPerHklCalculation(
        int hkl_idx,
        std::complex<double> &structureFactor,
        discamb::SfDerivativesAtHkl &derivatives);
        \endcode
        is called (hkl_idx starts from 0). See e.g. Tests/correctness/test_1/GradientDataCollector.h for an example implementation.
        \param countAtomContribution (if present) flags for counting atom contribution to structure factor calculations, size should match
        number of atoms

        */


        template<typename T>
        void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Matrix3d> &localCoordinateSystems,
            const std::vector<Vector3i> &hkl,
            T &singleHklResultsProcessor,
            const std::vector<bool> &countAtomContribution) const;

		void calculateFormFactors(
			const Vector3i& hkl,
			std::vector<std::complex<double> >& formFactors,
			const std::vector<bool>& includeAtom) const;

        void calculateFormFactors(
            const std::vector<Vector3i>& hkl,
            std::vector< std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;


        //void calculateFormFactorsFrac(
        //    const Vector3d& hkl,
        //    std::vector<std::complex<double> >& formFactors,
        //    const std::vector<bool>& includeAtom) const;

        //void calculateFormFactorsFrac(
        //    const std::vector<Vector3d>& hkl,
        //    std::vector< std::vector<std::complex<double> > >& formFactors,
        //    const std::vector<bool>& includeAtom) const;

        //void calculateFormFactorsCart(
        //    const Vector3d& hkl,
        //    std::vector<std::complex<double> >& formFactors,
        //    const std::vector<bool>& includeAtom) const;

        //void calculateFormFactorsCart(
        //    const std::vector<Vector3d>& hkl,
        //    std::vector< std::vector<std::complex<double> > >& formFactors,
        //    const std::vector<bool>& includeAtom) const;


	};

    /** @}*/
    } // namespace discamb

#endif
