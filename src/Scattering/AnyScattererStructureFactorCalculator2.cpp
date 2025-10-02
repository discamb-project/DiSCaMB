#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/MathUtilities/math_utilities.h"

#include <set>
#include <algorithm>
#include <omp.h>

using namespace std;

namespace spc = discamb::structural_parameters_convention;

namespace discamb{

    AnyScattererStructureFactorCalculator2::AnyScattererStructureFactorCalculator2() 
    {
        mMultiManager = false;
        mAdpConvention = structural_parameters_convention::AdpConvention::U_cif;
        mXyzConvention = structural_parameters_convention::XyzCoordinateSystem::fractional;
    }

    
        /** Calculate structure factors and their derivatives*/
    AnyScattererStructureFactorCalculator2::AnyScattererStructureFactorCalculator2(
        const Crystal &crystal)
    {
        mMultiManager = false;
        mManager = shared_ptr<AtomicFormFactorCalculationsManager>(new IamFormFactorCalculationsManager(crystal));
        
        //

        mCrystal = crystal;
        mAnomalous.assign(crystal.atoms.size(), 0.0);
        mReciprocalLatticeUnitCell.set(crystal.unitCell);
        int nAtoms = crystal.atoms.size();
        vector<AtomInCrystal> &atoms = mCrystal.atoms;
        mConverter.set(crystal.unitCell);

        //setDerivativesConvention

        // space group 

        Matrix3d frac2cart = crystal.unitCell.getFractionalToCartesianMatrix();
        Matrix3d cart2frac = crystal.unitCell.getCartesianToFractionalMatrix();
        Matrix3d rotationFractional;
        Vector3d translationFractional;
        int symmOpIndex, nSymmOp;

        nSymmOp = mCrystal.spaceGroup.nSymmetryOperations();// .nSy nSymmetryOperationsInSubset();
        mSymOps.resize(nSymmOp);

        for (symmOpIndex = 0; symmOpIndex < nSymmOp; symmOpIndex++)
        {
            crystal.spaceGroup.getSpaceGroupOperation(symmOpIndex).get(rotationFractional, translationFractional);
            mSymOps[symmOpIndex].rotation = frac2cart * rotationFractional * cart2frac;
            mSymOps[symmOpIndex].translation = frac2cart * translationFractional;
        }

        // set mAtomicMultiplicityWeight


        double nSymmOpsAsDouble = double(crystal.spaceGroup.nSymmetryOperations());
        mAtomicMultiplicityWeight.resize(nAtoms);
        for (int i = 0; i<nAtoms; i++)
        {
            //mAtomicMultiplicityWeight[i] = atoms[i].multiplicity;
            mAtomicMultiplicityWeight[i] = atoms[i].multiplicity / nSymmOpsAsDouble;
        }

		// allocate buffers for atomic form factors

		mAtomsFormFactors.assign(nAtoms, vector<complex<double> >(nSymmOp, 0.0));
		mAtomsFormFactor.assign(nAtoms, 0.0);

        // allocate atomic positions and adps 

        mR_at.resize(nAtoms);
        mR_at_symm.resize(nAtoms, vector<Vector3d>(nSymmOp));


        mAdp.resize(nAtoms);
        int nAdpComponents;

        for (int atom_index = 0; atom_index < nAtoms; atom_index++)
        {

            nAdpComponents = atoms[atom_index].adp.size();
            mAdp[atom_index].resize(nAdpComponents);
        }

        // set input convetions

        mAdpConvention = crystal.adpConvention;
        mXyzConvention = crystal.xyzCoordinateSystem;

        // set structural data

        updateStructuralData(crystal.atoms);

    }

    AnyScattererStructureFactorCalculator2::~AnyScattererStructureFactorCalculator2() {}

    void AnyScattererStructureFactorCalculator2::updateStructuralData(
        const std::vector<AtomInCrystal> &atoms)
    {
        // set atomic positions
        int nAtoms = atoms.size();

        // reset occupancies
        for (int i = 0; i < nAtoms; i++)
            mCrystal.atoms[i].occupancy = atoms[i].occupancy;


        for (int i = 0; i < nAtoms; i++)
        {
            mConverter.convertXyz(atoms[i].coordinates, mR_at[i], mXyzConvention, spc::XyzCoordinateSystem::cartesian);
            for (int j = 0; j < mSymOps.size(); j++)
                mR_at_symm[i][j] = mSymOps[j].rotation * mR_at[i] + mSymOps[j].translation;
        }
        //  mCrystal.unitCell.fractionalToCartesian(atoms[i].coordinates, mR_at[i]);

        // set ADPs
        
        int nAdpComponents;
        double two_pi2 = 2 * M_PI*M_PI;
        for (int atom_index = 0; atom_index<nAtoms; atom_index++)
        {

            nAdpComponents = atoms[atom_index].adp.size();
            mAdp[atom_index].resize(nAdpComponents);
            if (nAdpComponents == 6)
            {
                //if (mInputAdpConvetion == spc::U_cart)
                //    mAdp[atom_index] = atoms[atom_index].adp;
                //else
                //    mConverter.convertADP(atoms[atom_index].adp, mAdp[atom_index], mCrystal.adpConvention, spc::U_cart);

                mConverter.convertADP(atoms[atom_index].adp, mAdp[atom_index], mAdpConvention, spc::AdpConvention:: U_cart);

                for (int i = 0; i<6; i++)
                    mAdp[atom_index][i] *= two_pi2;
            }
            else
            {
                if (nAdpComponents == 1)
                {
                    mAdp[atom_index] = atoms[atom_index].adp;
                    mAdp[atom_index][0] *= two_pi2;
                }
            }
        }

    }

    void AnyScattererStructureFactorCalculator2::makeCartesianHkl(
        const std::vector<Vector3i> &hklIndices,
        std::vector<Vector3d> &hklCartesian)
    {
        
        int nHkl = hklIndices.size();
        hklCartesian.resize(nHkl);
        for (int i = 0; i < nHkl; i++)
            mReciprocalLatticeUnitCell.fractionalToCartesian(hklIndices[i], hklCartesian[i]);

    }

    void AnyScattererStructureFactorCalculator2::setAtomicFormfactorManager(
		std::shared_ptr<AtomicFormFactorCalculationsManager> &manager)
    {
        mMultiManager = false;
        mManager = manager;
    }


    void AnyScattererStructureFactorCalculator2::update(
        const std::vector<AtomInCrystal> &atoms)
    {
        mManager->update(atoms);
        updateStructuralData(atoms);
    }
        /**
        Sets convention for derivatives calculation:
        \param dXyzConvention specifies if structure factor derivatives with respect to coordinates should be calculated for
        fractional or cartesian coordinates.
        \param dAdpConvention specifies if structure factor derivatives with respect to aniosotropic ADPs should be
        calculated for \f$U_{cart}\f$, \f$U_{cif}\f$ or \f$U^{*}\f$

        */

        void AnyScattererStructureFactorCalculator2::setParametersConvention(
            spc::XyzCoordinateSystem dXyzConvention,
            spc::AdpConvention dAdpConvention)
        {
            mAdpConvention = dAdpConvention;
            mXyzConvention = dXyzConvention;
        }

        /**
        Gets convention for derivatives calculation.
        \sa setDerivativesConvention
        */
        void AnyScattererStructureFactorCalculator2::getParametersConvention(
            structural_parameters_convention::XyzCoordinateSystem &dXyzConvention,
            structural_parameters_convention::AdpConvention &dAdpConvention) 
        const 
        {
            dAdpConvention = mAdpConvention;
            dXyzConvention = mXyzConvention;
        }


        void AnyScattererStructureFactorCalculator2::calculateStructureFactors(
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f)
        {
            std::vector<bool> countAtomContribution(mCrystal.atoms.size(), true);          
            calculateStructureFactors(hkl, f, countAtomContribution);
        }



        
        void AnyScattererStructureFactorCalculator2::setAnoumalous(
            const std::vector<std::complex<double> > &anoumalous)
        {
            mAnomalous = anoumalous;
        }

		void AnyScattererStructureFactorCalculator2::calculateFormFactors(
			const Vector3i& hkl,
			std::vector<std::complex<double> >& formFactors,
			const std::vector<bool>& includeAtom)
			const
		{
			mManager->calculateFrac(hkl, formFactors, includeAtom);
		}

        void AnyScattererStructureFactorCalculator2::calculateFormFactors(
            const std::vector<Vector3i>& hkl,
            std::vector< std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom)
            const
        {
            mManager->calculateFrac(hkl, formFactors, includeAtom);
        }


        //int  AnyScattererStructureFactorCalculator2::findPreferredHklOrderingDirection(
        //    const std::vector<Vector3i>& hkl,
        //    std::vector<std::vector<Vector3i> >& orderedHklLines,
        //    std::vector<std::vector<int> >& mapToOriginalSetIndices)
        //{
        //    orderedHklLines.clear();

        //    set<vector<int> > lines[3];
        //    for(auto const &h: hkl)
        //    {
        //        lines[0].insert({ h[1], h[2] });
        //        lines[1].insert({ h[0], h[2] });
        //        lines[2].insert({ h[0], h[1] });
        //    }
        //    
        //    int preferredDirection = 0;
        //    int nLines = lines[0].size();
        //    for (int i = 1; i < 3; i++)
        //    {
        //        if (nLines > lines[i].size())
        //        {
        //            preferredDirection = i;
        //            nLines = lines[i].size();
        //        }
        //    }
        //    map<vector<int>, int> line2idx;
        //    int idx = 0;
        //    for (auto& pq : lines[preferredDirection])
        //        line2idx[pq] = idx++;

        //    orderedHklLines.resize(nLines);
        //    vector<int> otherIndices;   
        //    if (preferredDirection == 0)
        //        otherIndices = { 1, 2 };
        //    else if (preferredDirection == 1)
        //        otherIndices = { 0, 2 };
        //    else
        //        otherIndices = { 0, 1 };

        //    for (auto const& h : hkl)
        //    {
        //        vector<int> pq = { h[otherIndices[0]], h[otherIndices[1]]};
        //        orderedHklLines[line2idx[pq]].push_back(h);
        //    }

        //    for (auto& line : orderedHklLines)
        //        std::sort(line.begin(), line.end());

        //    mapToOriginalSetIndices.clear();
        //    mapToOriginalSetIndices.resize(nLines);
        //    for (int i = 0; i < nLines; i++)
        //    {
        //        for (auto const& h : orderedHklLines[i])
        //        {
        //            auto it = find(hkl.begin(), hkl.end(), h);
        //            int idx = distance(hkl.begin(), it);
        //            mapToOriginalSetIndices[i].push_back(idx);
        //        }
        //    }

        //    return preferredDirection;
        //}

        void AnyScattererStructureFactorCalculator2::calculate_line_phase_factors(
            //in:
            const std::vector<Vector3d>& hkl_cart,
            int hkl_idx,
            const std::vector<Vector3i>& hkl_line,
            int idx_in_line,
            int line_direction,
            const std::vector<Vector3<double> >& rotated_h,
            const std::vector< std::vector<std::complex<double > > >& phase_factor_multiplier,
            //out:
            std::vector< std::vector<std::complex<double > > >& line_phase_factor)
        {
            int nAtoms = mR_at.size();
            int nSymmOps = mSymOps.size();
            double two_pi = 2 * M_PI;
            if (idx_in_line == 0)
            {
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        double phase_angle = two_pi * mR_at_symm[atomIdx][symOpIdx] * hkl_cart[hkl_idx];
                        line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                    }
                }
            }
            else
            {


                int lineBreak = hkl_line[idx_in_line][line_direction] - hkl_line[idx_in_line - 1][line_direction] - 1;
                if (lineBreak == 0)
                    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                            line_phase_factor[atomIdx][symOpIdx] *= phase_factor_multiplier[atomIdx][symOpIdx];
                else
                    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        {
                            double phase_angle = two_pi * mR_at_symm[atomIdx][symOpIdx] * hkl_cart[hkl_idx];
                            //double phase_angle = rotated_h_2pi[symOpIdx] * mR_at[atomIdx] + translation_factor_2pi[symOpIdx];
                            line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                        }

            }


        }

        void AnyScattererStructureFactorCalculator2::calculate_line_temperature_factors(
            //in:
            const std::vector<Vector3d>& hkl_cart,
            int hkl_idx,
            const std::vector<Vector3i>& hkl_line,
            int idx_in_line,
            int line_direction,
            const std::vector<Vector3<double> >& rotated_h,
            std::vector< std::vector<double> >& temp_factor_multiplier,
            const std::vector< std::vector<double> >& temp_factor_multiplier_multiplier,
            std::vector<std::vector<double> >& adpMultipliers,
            //out:
            std::vector< std::vector<double> >& line_temperature_factors)
        {
            int nAtoms = mR_at.size();
            int nSymmOps = mSymOps.size();
            double two_pi = 2 * M_PI;

            //bool calc_tf_from_scratch; // first point after break/start or second point after break/start, but no continuation
            //bool calc_tf_and_multiplier_from_scratch; // second point after break/start and there will be continuation
            //bool calc_tf_from_multiplier; // third or further points afger start/break

            int points_after_break = 0;
            int points_before_break = 0;
            vector<int> line;
            for (int i = 0; i < hkl_line.size(); i++)
                line.push_back(hkl_line[i][line_direction]);

            for (int i = idx_in_line + 1; i < line.size(); i++)
            {
                if (line[i] == line[i - 1] + 1)
                    points_before_break++;
                else
                    break;
            }
            for (int i = idx_in_line - 1; i >= 0; i--)
            {
                if (line[i] == line[i + 1] - 1)
                    points_after_break++;
                else
                    break;
            }
            double const* adps;
            double hVectorLength2 = hkl_cart[hkl_idx] * hkl_cart[hkl_idx];

            // first point after break/start or second point after break/start, but no continuation
            // just calculate temperature factors from scratch
            if (points_after_break == 0 || (points_after_break == 1 && points_before_break == 0))
            {
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {

                    int n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                    if (n_adp_components > 0)
                        adps = &mAdp[atomIdx][0];

                    if (n_adp_components == 6)
                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        {
                            double* multipliers = &adpMultipliers[symOpIdx][0];
                            line_temperature_factors[atomIdx][symOpIdx] = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                - multipliers[4] * adps[4] - multipliers[5] * adps[5]);
                        }
                    else if (n_adp_components == 1)
                        line_temperature_factors[atomIdx][0] = exp(-hVectorLength2 * (*adps));

                }
            }
            else if (points_after_break == 1 && points_before_break > 0)
            {
                // second point after break/start and there will be continuation
                
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {

                    int n_adp_components = mAdp[atomIdx].size();

                    if (n_adp_components > 0)
                        adps = &mAdp[atomIdx][0];

                    if (n_adp_components == 6)
                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        {
                            double* multipliers = &adpMultipliers[symOpIdx][0];
                            double previous_temp_factor = line_temperature_factors[atomIdx][symOpIdx];
                            line_temperature_factors[atomIdx][symOpIdx] = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                - multipliers[4] * adps[4] - multipliers[5] * adps[5]);
                            temp_factor_multiplier[atomIdx][symOpIdx] =
                                line_temperature_factors[atomIdx][symOpIdx] / previous_temp_factor;
                        }
                    else if (n_adp_components == 1)
                    {
                        double previous_temp_factor = line_temperature_factors[atomIdx][0];
                        line_temperature_factors[atomIdx][0] = exp(-hVectorLength2 * (*adps));
                        temp_factor_multiplier[atomIdx][0] =
                            line_temperature_factors[atomIdx][0] / previous_temp_factor;
                    }

                }
            }
            else
            {
                // third or further points after start/break
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    int n_adp_components = mAdp[atomIdx].size();
                    
                    if(n_adp_components == 6)
                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        {
                            temp_factor_multiplier[atomIdx][symOpIdx] *= temp_factor_multiplier_multiplier[atomIdx][symOpIdx];
                            line_temperature_factors[atomIdx][symOpIdx] *= temp_factor_multiplier[atomIdx][symOpIdx];
                        }
                    else if (n_adp_components == 1)
                    {
                        temp_factor_multiplier[atomIdx][0] *= temp_factor_multiplier_multiplier[atomIdx][0];
                        line_temperature_factors[atomIdx][0] *= temp_factor_multiplier[atomIdx][0];
                    }
                        
                }
            }
                
            // first point after break/start or second point after break/start, but no continuation
            // third or further points afger start/break



            //if (idx_in_line == 0)
            //{
            //    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            //    {
            //        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
            //        {
            //            line_temperature_factors[atomIdx][symOpIdx] = temp_factor_multiplier[atomIdx][symOpIdx];
            //        }
            //    }
            //}
            //else
            //{
            //    int lineBreak = hkl_line[idx_in_line][line_direction] - hkl_line[idx_in_line - 1][line_direction] - 1;
            //    if (lineBreak == 0)
            //        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            //            for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
            //                line_temperature_factors[atomIdx][symOpIdx] *= temp_factor_multiplier[atomIdx][symOpIdx];
            //    else
            //        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            //            for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
            //            {
            //                line_temperature_factors[atomIdx][symOpIdx] = temp_factor_multiplier_multiplier[atomIdx][symOpIdx];
            //            }
            //}
        }


        void AnyScattererStructureFactorCalculator2::calculateStructureFactors_ordered_hkl_alg(
            const std::vector<Vector3i>& _hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution)
        {
            std::vector<Vector3d> hkl;
            makeCartesianHkl(_hkl, hkl);



            const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));

            const REAL two_pi = 2 * REAL(M_PI);
            const REAL two_pi_sqare = 2 * REAL(M_PI * M_PI);
            REAL temperature_factor, hVectorLength, hVectorLength2;
            REAL const* adps;
            vector<Vector3<REAL> > rotated_h(mSymOps.size());

            vector<REAL> translation_factor(mSymOps.size());
            complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum;

            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
            REAL realFContrib, imagFContrib;


            //REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];

            complex<REAL> anomalousScattering, atom_ff, aux;


            

            //--

            int nSymmOps = mSymOps.size();
            int n_adp_components, nAtoms, nHklVectors = hkl.size();
            vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));
            nAtoms = mR_at.size();


            vector<double> atomWeights(nAtoms);
            for (int i = 0; i < nAtoms; i++)
                atomWeights[i] = mCrystal.atoms[i].occupancy * mAtomicMultiplicityWeight[i];

            f.assign(nHklVectors, 0.0);


            //

            vector < vector<complex<double> > > atomic_ff(nSymmOps, vector < complex<double> >(nAtoms));
            vector < vector<complex<double> > > atomic_ff_anom(nSymmOps, vector < complex<double> >(nAtoms));
            vector<double> phaseFactor(nAtoms * nSymmOps);
            vector<Vector3d> rotated_h_2pi(nSymmOps);
            vector<double> translation_factor_2pi(nSymmOps);

            vector<vector<Vector3i> > orderedHklLines;
            vector<vector<int> > mapToOriginalSetIndices;
            int lineDirection = scattering_utilities::findPreferredHklOrderingDirection(_hkl, orderedHklLines, mapToOriginalSetIndices);
            int nLines = orderedHklLines.size();
            //[atom][symmOp]
            vector<vector<complex<double> > > phase_factor_multiplier(nAtoms, vector<complex<double> >(nSymmOps));
            vector<vector<complex<double> > > line_phase_factor(nAtoms, vector<complex<double> >(nSymmOps));
            //vector<vector<double> > temp_factor_0(nAtoms, vector<double>(nSymmOps));
            //vector<vector<double> > temp_factor_1(nAtoms, vector<double>(nSymmOps)); 
            vector<vector<double> > temp_factor(nAtoms, vector<double>(nSymmOps));
            vector<vector<double> > temp_factor_multiplier_iter(nAtoms, vector<double>(nSymmOps));
            vector<vector<double> > temp_factor_multiplier_multiplier(nAtoms, vector<double>(nSymmOps));

            Vector3d step, step_frac(0.0,0.0,0.0);
            step_frac[lineDirection] = 1.0;
            mReciprocalLatticeUnitCell.fractionalToCartesian(step_frac, step);

            
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                {
                    double phase_angle = two_pi * mR_at_symm[atomIdx][symOpIdx] * step;
                    phase_factor_multiplier[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                }

            //temp_factor_multiplier_multiplier
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                int nADPComponents = mAdp[atomIdx].size();
                if(nADPComponents==1)
                    temp_factor_multiplier_multiplier[atomIdx][0] = exp(-2.0 * step * step * mAdp[atomIdx][0]);
                else if (nADPComponents == 6)
                {
                    
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        Vector3d step_rot = step * mSymOps[symOpIdx].rotation;
                        double exponent = step_rot[0] * step_rot[0] * mAdp[atomIdx][0] +
                                          step_rot[1] * step_rot[1] * mAdp[atomIdx][1] +
                                          step_rot[2] * step_rot[2] * mAdp[atomIdx][2] +
                                   2.0 * (step_rot[0] * step_rot[1] * mAdp[atomIdx][3] +
                                          step_rot[0] * step_rot[2] * mAdp[atomIdx][4] +
                                          step_rot[1] * step_rot[2] * mAdp[atomIdx][5]);
                        temp_factor_multiplier_multiplier[atomIdx][symOpIdx] = exp(-2.0 * exponent);
                    }
                }
            }
            
            for (int lineIdx = 0; lineIdx < nLines; lineIdx++)
            {
                auto& hklLine = orderedHklLines[lineIdx];
                auto& mapToOriginalSet = mapToOriginalSetIndices[lineIdx];
                int nHklLine = hklLine.size();

                /*
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        double phase_angle = two_pi * mR_at_symm[atomIdx][symOpIdx] * step;
                        line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                        
                    }
                */


                for (int idxInLine = 0; idxInLine < nHklLine; idxInLine++)
                {
                    int hklIndex = mapToOriginalSet[idxInLine];
                    hVectorLength2 = hkl[hklIndex] * hkl[hklIndex];


                    pre_atom_loop_sf_calc( /*  in: */  mSymOps, hkl[hklIndex],
                        /* out: */ rotated_h, translation_factor, adpMultipliers);


                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        mManager->calculateCart(rotated_h[symOpIdx], atomic_ff[symOpIdx], countAtomContribution);
                        rotated_h_2pi[symOpIdx] = rotated_h[symOpIdx] * two_pi;
                        translation_factor_2pi[symOpIdx] = translation_factor[symOpIdx] * two_pi;
                    }


                    
                    //if (idxInLine == 0)
                    //{
                    //    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    //    {
                    //        int nADPComponents = mAdp[atomIdx].size();
                    //        
                    //        if (nADPComponents == 1)
                    //            temp_factor[atomIdx][0] = exp(-hVectorLength2 * mAdp[atomIdx][0]);
                    //        if (nADPComponents == 6)
                    //        {
                    //            for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    //            {
                    //                double* multipliers = &adpMultipliers[symOpIdx][0];
                    //                temp_factor[atomIdx][0] = exp(-multipliers[0] * mAdp[atomIdx][0] - multipliers[1] * mAdp[atomIdx][1]
                    //                    - multipliers[2] * mAdp[atomIdx][2] - multipliers[3] * mAdp[atomIdx][3]
                    //                    - multipliers[4] * mAdp[atomIdx][4] - multipliers[5] * mAdp[atomIdx][5]);
                    //            }
                    //        }

                    //        //for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    //        //{
                    //        //    double phase_angle = two_pi * mR_at_symm[atomIdx][symOpIdx] * hkl[hklIndex];
                    //        //    line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                    //        //}
                    //    }
                    //}
                    //else
                    //{
                    //    if (idxInLine == 1)
                    //    {
                    //        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    //        {
                    //            int nADPComponents = mAdp[atomIdx].size();

                    //            if (nADPComponents == 1)
                    //            {
                    //                double new_temp_factor = exp(-hVectorLength2 * mAdp[atomIdx][0]);
                    //                temp_factor_multiplier_iter[atomIdx][0] = new_temp_factor / temp_factor[atomIdx][0];
                    //                temp_factor[atomIdx][0] = new_temp_factor;
                    //                temp_factor_multiplier_multiplier[atomIdx][0] = exp(-2.0 * step * step * mAdp[atomIdx][0]);
                    //            }
                    //            if (nADPComponents == 6)
                    //            {
                    //                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    //                {
                    //                    double* multipliers = &adpMultipliers[symOpIdx][0];
                    //                    double new_temp_factor = exp(-multipliers[0] * mAdp[atomIdx][0] - multipliers[1] * mAdp[atomIdx][1]
                    //                        - multipliers[2] * mAdp[atomIdx][2] - multipliers[3] * mAdp[atomIdx][3]
                    //                        - multipliers[4] * mAdp[atomIdx][4] - multipliers[5] * mAdp[atomIdx][5]);
                    //                    temp_factor_multiplier_iter[atomIdx][symOpIdx] = new_temp_factor / temp_factor[atomIdx][symOpIdx];
                    //                    temp_factor[atomIdx][symOpIdx] = new_temp_factor;
                    //                    double exponent = -2.0 * 
                    //                        ( step[0] * step[0] * mAdp[atomIdx][0] +
                    //                          step[1] * step[1] * mAdp[atomIdx][1] +
                    //                          step[2] * step[2] * mAdp[atomIdx][2] +
                    //                          2.0 * step[0] * step[1] * mAdp[atomIdx][3] + 
                    //                          2.0 * step[0] * step[2] * mAdp[atomIdx][4] +
                    //                          2.0 * step[1] * step[3] * mAdp[atomIdx][5] );
                    //                    temp_factor_multiplier_multiplier[atomIdx][symOpIdx] = exp(exponent);
                    //                }
                    //            }

                    //        }
                    //    }


                    //    //int lineBreak = hklLine[idxInLine][lineDirection] - hklLine[idxInLine - 1][lineDirection] - 1;
                    //    //if(lineBreak == 0)
                    //    //    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    //    //        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    //    //            line_phase_factor[atomIdx][symOpIdx] *= phase_factor_multiplier[atomIdx][symOpIdx];
                    //    //else
                    //    //    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    //    //        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    //    //        {
                    //    //            double phase_angle = rotated_h_2pi[symOpIdx] * mR_at[atomIdx] + translation_factor_2pi[symOpIdx]; 
                    //    //            line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                    //    //        }

                    //}

                    calculate_line_temperature_factors(
                        //in:
                        hkl,
                        hklIndex,
                        hklLine,
                        idxInLine,
                        lineDirection,
                        rotated_h,
                        temp_factor_multiplier_iter,
                        temp_factor_multiplier_multiplier,
                        adpMultipliers,
                        //out:
                        temp_factor);

                    calculate_line_phase_factors(
                        //in:
                        hkl,
                        hklIndex,
                        hklLine,
                        idxInLine,
                        lineDirection,
                        rotated_h,
                        phase_factor_multiplier,
                        //rotated_h_2pi,
                        //translation_factor_2pi,

                        //out:
                        line_phase_factor);


                    realFContrib = 0;
                    imagFContrib = 0;

 
                    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    {


                        if (!countAtomContribution[atomIdx])
                            continue;


                        double atomWeight = atomWeights[atomIdx];// mCrystal.atoms[atomIdx].occupancy* mAtomicMultiplicityWeight[atomIdx];

                        n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                        if (n_adp_components > 0)
                            adps = &mAdp[atomIdx][0];

                        unweightedTransformedAtomFF_Sum = 0.0;

                        for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        {

                            complex<REAL> atomic_position_phase_factor = line_phase_factor[atomIdx][symOpIdx];

                            atom_ff = atomic_ff[symOpIdx][atomIdx] + mAnomalous[atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);

                            if (n_adp_components == 6)
                            {
                                double* multipliers = &adpMultipliers[symOpIdx][0];
                                //temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                //    - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                //    - multipliers[4] * adps[4] - multipliers[5] * adps[5]);
                                temperature_factor = temp_factor[atomIdx][symOpIdx];
                                unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;
                            }
                            else
                                unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

                            unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

                        } // symmetry operations


                        if (n_adp_components == 1)
                        {
                            temperature_factor = temp_factor[atomIdx][0];
                            //temperature_factor = exp(-hVectorLength2 * (*adps));
                            unweightedTransformedAtomFF_Sum *= temperature_factor;
                        }

                        realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
                        imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

                    } // symmetrically independent atoms

                    f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);
                    //cin >> hklIndex;
                } // h vectors in line
            } // lines


        }


        void calculateStructureFactors_iso_aniso_split(
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution)
        {

        }


		void AnyScattererStructureFactorCalculator2::calculateStructureFactors(
			const std::vector<Vector3i>& _hkl,
			std::vector<std::complex<double> >& f,
			const std::vector<bool>& countAtomContribution)
		{
            calculateStructureFactors_ordered_hkl_alg(_hkl, f, countAtomContribution); 
            return;

			std::vector<Vector3d> hkl;
			makeCartesianHkl(_hkl, hkl);



			const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
			
			const REAL two_pi = 2 * REAL(M_PI);
			const REAL two_pi_sqare = 2 * REAL(M_PI * M_PI);
			REAL temperature_factor, hVectorLength, hVectorLength2;
			REAL const* adps;
			vector<Vector3<REAL> > rotated_h(mSymOps.size());

			vector<REAL> translation_factor(mSymOps.size());
			complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum;

			REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
			REAL realFContrib, imagFContrib;


			//REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];

			complex<REAL> anomalousScattering, atom_ff, aux;

			//--

			int nSymmOps = mSymOps.size();
			int n_adp_components, nAtoms, nHklVectors = hkl.size();
			vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));
			nAtoms = mR_at.size();

            vector<double> atomWeights(nAtoms);
            for (int i = 0; i < nAtoms; i++)
                atomWeights[i] = mCrystal.atoms[i].occupancy * mAtomicMultiplicityWeight[i];
            
			f.assign(nHklVectors, 0.0);


			//

			vector < vector<complex<double> > > atomic_ff(nSymmOps, vector < complex<double> >(nAtoms));
            vector < vector<complex<double> > > atomic_ff_anom(nSymmOps, vector < complex<double> >(nAtoms));
            vector<double> phaseFactor(nAtoms * nSymmOps);
            vector<Vector3d> rotated_h_2pi(nSymmOps);
            vector<double> translation_factor_2pi(nSymmOps);
			for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
			{

                //cout << hklIndex << "/" << nHklVectors << endl;
				hVectorLength2 = hkl[hklIndex] * hkl[hklIndex];
				//hVectorLength = sqrt(hVectorLength2);

				pre_atom_loop_sf_calc( /*  in: */  mSymOps, hkl[hklIndex],
					/* out: */ rotated_h, translation_factor, adpMultipliers);

                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                {
                    mManager->calculateCart(rotated_h[symOpIdx], atomic_ff[symOpIdx], countAtomContribution);
                    rotated_h_2pi[symOpIdx] = rotated_h[symOpIdx] * two_pi;
                    translation_factor_2pi[symOpIdx] = translation_factor[symOpIdx] * two_pi;
                }

                int idx = 0;
                
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        phaseFactor[idx] = rotated_h_2pi[symOpIdx] * mR_at[atomIdx] + translation_factor_2pi[symOpIdx];
                        idx++;
                    }
                //for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                //{
                //    Vector3d rotated_h_2pi = rotated_h[symOpIdx] * two_pi;
                //    double tf_2pi = translation_factor[symOpIdx] * two_pi;
                //    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                //    {
                //        phaseFactor[idx] = rotated_h_2pi * mR_at[atomIdx] + tf_2pi;
                //        ++idx;
                //    }
                //}


                /*
                atomic_ff_anom = atomic_ff;
                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                        atomic_ff_anom[symOpIdx][atomIdx] += mAnomalous[atomIdx];
                        */
				realFContrib = 0;
				imagFContrib = 0;

                idx = 0;
				for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
				{

					//cout << endl << endl << mCrystal.atoms[atomIdx].label << endl;

					if (!countAtomContribution[atomIdx])
						continue;


                    double atomWeight = atomWeights[atomIdx];// mCrystal.atoms[atomIdx].occupancy* mAtomicMultiplicityWeight[atomIdx];
                    
					n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

					if (n_adp_components > 0)
						adps = &mAdp[atomIdx][0];

					unweightedTransformedAtomFF_Sum = 0.0;
                    
					for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
					{
                        /*
						const Vector3d& rotated_h_ref = rotated_h[symOpIdx];

						atomic_phase_factor_phase = two_pi * (rotated_h_ref * mR_at[atomIdx] +
							translation_factor[symOpIdx]);
                            */
                        atomic_phase_factor_phase = phaseFactor[idx];
                        ++idx;
                        atomic_phase_factor_real = cos(atomic_phase_factor_phase);
						atomic_phase_factor_im = sin(atomic_phase_factor_phase);
                        //atomic_phase_factor_im = sqrt(max(0.0,1.0 - atomic_phase_factor_real * atomic_phase_factor_real));
                        //if (atomic_phase_factor_phase < 0)
                          //  atomic_phase_factor_im = -atomic_phase_factor_im;
                        //atomic_phase_factor_im = atomic_phase_factor_real + atomic_phase_factor_phase;
                        

						complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);
                        
                        //complex<REAL> atomic_position_phase_factor(1.0, 1.0);// omic_phase_factor_real, atomic_phase_factor_im);
                        //atom_ff = atomic_ff_anom[symOpIdx][atomIdx];// atomic_ff[symOpIdx][atomIdx] + mAnomalous[atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);
                        atom_ff = atomic_ff[symOpIdx][atomIdx] + mAnomalous[atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);

						if (n_adp_components == 6)
						{
							double* multipliers = &adpMultipliers[symOpIdx][0];
                            temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
								- multipliers[2] * adps[2] - multipliers[3] * adps[3]
								- multipliers[4] * adps[4] - multipliers[5] * adps[5]);

							unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;
						}
						else
							unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

						unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

					} // symmetry operations


					if (n_adp_components == 1)
					{
						temperature_factor = exp(-hVectorLength2 * (*adps));
						unweightedTransformedAtomFF_Sum *= temperature_factor;
					}

					realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
					imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

				} // symmetrically independent atoms

				f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);
				//cin >> hklIndex;
			} // h vectors

		}

// parallel
        void AnyScattererStructureFactorCalculator2::calculateStructureFactors2(
            const std::vector<Vector3i>& _hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution)
        {
            std::vector<Vector3d> hkl;
            makeCartesianHkl(_hkl, hkl);

            const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
            const REAL two_pi = 2 * REAL(M_PI);
            const REAL two_pi_sqare = 2 * REAL(M_PI * M_PI);
            const int nSymmOps = mSymOps.size();
            const int nHklVectors = int(hkl.size());
            const int nAtoms = mR_at.size();

            f.assign(nHklVectors, 0.0);
            

#pragma omp parallel for
            for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
            {

                //
                vector < vector<complex<double> > > atomic_ff(nSymmOps, vector < complex<double> >(nAtoms));
                REAL temperature_factor, hVectorLength, hVectorLength2;
                REAL const* adps;
                vector<Vector3<REAL> > rotated_h(nSymmOps);

                vector<REAL> translation_factor(nSymmOps);
                complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum, anomalousScattering, atom_ff, aux;

                REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
                REAL realFContrib, imagFContrib, atomWeight;

                //--


                int n_adp_components;
                vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));

                 
                //


                hVectorLength2 = hkl[hklIndex] * hkl[hklIndex];
                hVectorLength = sqrt(hVectorLength2);

                pre_atom_loop_sf_calc( /*  in: */  mSymOps, hkl[hklIndex],
                    /* out: */ rotated_h, translation_factor, adpMultipliers);

                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    mManager->calculateCart(rotated_h[symOpIdx], atomic_ff[symOpIdx], countAtomContribution);


                realFContrib = 0;
                imagFContrib = 0;


                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {

                    //cout << endl << endl << mCrystal.atoms[atomIdx].label << endl;

                    if (!countAtomContribution[atomIdx])
                        continue;


                    atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];

                    n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                    if (!n_adp_components == 0)
                        adps = &mAdp[atomIdx][0];

                    unweightedTransformedAtomFF_Sum = 0.0;

                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {

                        const Vector3d& rotated_h_ref = rotated_h[symOpIdx];

                        atomic_phase_factor_phase = two_pi * (rotated_h_ref * mR_at[atomIdx] +
                            translation_factor[symOpIdx]);

                        atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                        atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                        complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


                        atom_ff = atomic_ff[symOpIdx][atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);


                        if (n_adp_components == 6)
                        {
                            double* multipliers = &adpMultipliers[symOpIdx][0];
                            temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                - multipliers[4] * adps[4] - multipliers[5] * adps[5]);

                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;
                        }
                        else
                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

                        unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

                    } // symmetry operations


                    if (n_adp_components == 1)
                    {
                        temperature_factor = exp(-hVectorLength * hVectorLength * (*adps));
                        unweightedTransformedAtomFF_Sum *= temperature_factor;
                    }

                    realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
                    imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

                } // symmetrically independent atoms

                f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);
                //cin >> hklIndex;
            } // h vectors

        }

//


        /** \brief Calculates structure factors.

        \param atoms carry atomic structural parameters
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
        \param f the resulting structure factors
        */
        
        //void AnyScattereStructureFactorCalculator::calculateStructureFactors(
        //    const std::vector<AtomInCrystal> &atoms,
        //    const std::vector<Matrix3d> &localCoordinateSystems,
        //    const std::vector<Vector3i> &hkl,
        //    std::vector<std::complex<double> > &f) {}


        void AnyScattererStructureFactorCalculator2::calculateStructureFactorsAndDerivatives_ordered_hkl_alg(
            const std::vector<Vector3i>& _hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& _dTarget_df,
            const std::vector<bool>& countAtomContribution)
        {
            std::vector<Vector3d> hkl;
            makeCartesianHkl(_hkl, hkl);
            // local copies

            //vector<Vector3d> atomic_positions = _atomic_positions;
            //vector<vector<double> > atomic_displacement_parameters = _atomic_displacement_parameters;
            //vector<REAL> atomic_occupancy = _atomic_occupancy;
            //vector<REAL> atomic_multiplicity_weight = _atomic_multiplicity_weight;

            vector<vector<Vector3i> > orderedHklLines;
            vector<vector<int> > mapToOriginalSetIndices;
            int lineDirection = scattering_utilities::findPreferredHklOrderingDirection(_hkl, orderedHklLines, mapToOriginalSetIndices);
            int nLines = orderedHklLines.size();



            //--


            const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
            //complex<REAL> fAtomSphericalAndAnomalous;
            const REAL two_pi = 2 * REAL(M_PI);
            const REAL two_pi_sqare = 2 * REAL(M_PI * M_PI);
            REAL temperature_factor, hVectorLength, hVectorLength2, multiplier;
            //REAL temperature_factor, multiplier, hVectorLength;
            REAL const* adps;
            vector<Vector3<REAL> > rotated_h(mSymOps.size());

            vector<REAL> translation_factor(mSymOps.size());
            complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum, dTargetDf;

            complex<REAL> adp_derivatives[6];
            complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
            REAL realFContrib, imagFContrib;


            REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];

            complex<REAL> anomalousScattering, atom_ff, aux;

            //--

            int nSymmOps = mSymOps.size();
            int n_adp_components, nAtoms, nHklVectors = hkl.size();
            vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));
            nAtoms = mR_at.size();

            vector<complex<double> > dTarget_df = _dTarget_df;
            for (int i = 0; i < nHklVectors; i++)
                dTarget_df[i] = conj(dTarget_df[i]);

            // set dTarget_dparam to zero ..

            f.assign(nHklVectors, 0.0);

            dTarget_dparam.resize(mR_at.size());

            for (int atom_index = 0; atom_index < nAtoms; atom_index++)
            {
                dTarget_dparam[atom_index].adp_derivatives.assign(mAdp[atom_index].size(), 0.0);
                dTarget_dparam[atom_index].atomic_position_derivatives = Vector3d(0, 0, 0);
                dTarget_dparam[atom_index].occupancy_derivatives = 0.0;
            }

            //

            vector < vector<complex<double> > > atomic_ff(nSymmOps, vector < complex<double> >(nAtoms));


            for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
            {

                //cout << _hkl[hklIndex][0] << " " << _hkl[hklIndex][1] << " " << _hkl[hklIndex][2] << "\n";

                dTargetDf = dTarget_df[hklIndex];

                hVectorLength2 = hkl[hklIndex] * hkl[hklIndex];
                hVectorLength = sqrt(hVectorLength2);

                pre_atom_loop_sf_calc( /*  in: */  mSymOps, hkl[hklIndex],
                    /* out: */ rotated_h, translation_factor, adpMultipliers);

                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    mManager->calculateCart(rotated_h[symOpIdx], atomic_ff[symOpIdx], countAtomContribution);

                realFContrib = 0;
                imagFContrib = 0;


                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {

                    //cout << endl << endl << mCrystal.atoms[atomIdx].label << endl;

                    if (!countAtomContribution[atomIdx])
                        continue;


                    atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];
                    //anomalousScattering = wfnParams[atomWfnIdx].anomalous_scattering;
                    //fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

                    n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                    if (!n_adp_components == 0)
                        adps = &mAdp[atomIdx][0];

                    if (n_adp_components == 6)
                        for (int i = 0; i < 6; i++)
                            adp_derivatives[i] = 0.0;

                    xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

                    unweightedTransformedAtomFF_Sum = 0.0;

                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {

                        const Vector3d& rotated_h_ref = rotated_h[symOpIdx];

                        atomic_phase_factor_phase = two_pi * (rotated_h_ref * mR_at[atomIdx] +
                            translation_factor[symOpIdx]);

                        atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                        atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                        complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


                        atom_ff = atomic_ff[symOpIdx][atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);


                        if (n_adp_components == 6)
                        {
                            double* multipliers = &adpMultipliers[symOpIdx][0];
                            temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                - multipliers[4] * adps[4] - multipliers[5] * adps[5]);

                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;

                            adp_derivatives[0] -= multipliers[0] * unweightedTransformedAtomFF;
                            adp_derivatives[1] -= multipliers[1] * unweightedTransformedAtomFF;
                            adp_derivatives[2] -= multipliers[2] * unweightedTransformedAtomFF;
                            adp_derivatives[3] -= multipliers[3] * unweightedTransformedAtomFF;
                            adp_derivatives[4] -= multipliers[4] * unweightedTransformedAtomFF;
                            adp_derivatives[5] -= multipliers[5] * unweightedTransformedAtomFF;

                        }
                        else
                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

                        unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

                        xyz_derivatives[0] += rotated_h_ref[0] * unweightedTransformedAtomFF;
                        xyz_derivatives[1] += rotated_h_ref[1] * unweightedTransformedAtomFF;
                        xyz_derivatives[2] += rotated_h_ref[2] * unweightedTransformedAtomFF;

                        //cout << atom_ff << " " << unweightedTransformedAtomFF << endl;

                    } // symmetry operations


                    if (n_adp_components == 1)
                    {
                        temperature_factor = exp(-hVectorLength * hVectorLength * (*adps));
                        unweightedTransformedAtomFF_Sum *= temperature_factor;
                        dTarget_dparam[atomIdx].adp_derivatives[0] -= hVectorLength2 * (dTargetDf.real() * unweightedTransformedAtomFF_Sum.real() -
                            dTargetDf.imag() * unweightedTransformedAtomFF_Sum.imag());
                    }
                    else
                        if (n_adp_components == 6)
                            for (int i = 0; i < 6; ++i)
                                dTarget_dparam[atomIdx].adp_derivatives[i] += dTargetDf.real() * adp_derivatives[i].real() -
                                dTargetDf.imag() * adp_derivatives[i].imag();


                    dTarget_dparam[atomIdx].occupancy_derivatives += unweightedTransformedAtomFF_Sum.real() * dTargetDf.real() -
                        unweightedTransformedAtomFF_Sum.imag() * dTargetDf.imag();
                    //cout << mCrystal.atoms[atomIdx].label << " " << unweightedTransformedAtomFF_Sum * atomWeight << endl;
                    realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
                    imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

                    //cout << atomIdx << " " << unweightedTransformedAtomFF_Sum.real()* atomWeight << " " << unweightedTransformedAtomFF_Sum.imag()* atomWeight << "\n";

                    n_adp_components == 1 ? aux = temperature_factor * dTargetDf : aux = dTargetDf;

                    for (int i = 0; i < 3; ++i)
                        dTarget_dparam[atomIdx].atomic_position_derivatives[i] -= aux.real() * xyz_derivatives[i].imag() +
                        aux.imag() * xyz_derivatives[i].real();

                } // symmetrically independent atoms

                f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);
                //cin >> hklIndex;
            } // h vectors

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];
                dTarget_dparam[atomIdx].occupancy_derivatives *= mAtomicMultiplicityWeight[atomIdx];
                multiplier = two_pi_sqare * atomWeight;
                for (int i = 0; i < dTarget_dparam[atomIdx].adp_derivatives.size(); i++)
                    dTarget_dparam[atomIdx].adp_derivatives[i] *= multiplier;
                multiplier = two_pi * atomWeight;
                dTarget_dparam[atomIdx].atomic_position_derivatives.x *= multiplier;
                dTarget_dparam[atomIdx].atomic_position_derivatives.y *= multiplier;
                dTarget_dparam[atomIdx].atomic_position_derivatives.z *= multiplier;
            }

            convertDerivatives(dTarget_dparam);
        }

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
        void AnyScattererStructureFactorCalculator2::calculateStructureFactorsAndDerivatives(
            const std::vector<Vector3i> &_hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &_dTarget_df,
            const std::vector<bool> &countAtomContribution) 
        {
            std::vector<Vector3d> hkl;
            makeCartesianHkl(_hkl, hkl);
            // local copies
            
            //vector<Vector3d> atomic_positions = _atomic_positions;
            //vector<vector<double> > atomic_displacement_parameters = _atomic_displacement_parameters;
            //vector<REAL> atomic_occupancy = _atomic_occupancy;
            //vector<REAL> atomic_multiplicity_weight = _atomic_multiplicity_weight;

            


            //--


            const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
            //complex<REAL> fAtomSphericalAndAnomalous;
            const REAL two_pi = 2 * REAL(M_PI);
            const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
            REAL temperature_factor, hVectorLength, hVectorLength2, multiplier;
            //REAL temperature_factor, multiplier, hVectorLength;
            REAL const *adps;
            vector<Vector3<REAL> > rotated_h(mSymOps.size());

            vector<REAL> translation_factor(mSymOps.size());
            complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum, dTargetDf;

            complex<REAL> adp_derivatives[6];
            complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
            REAL realFContrib, imagFContrib;

            
            REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];
            
            complex<REAL> anomalousScattering, atom_ff, aux;

            //--

            int nSymmOps = mSymOps.size();
            int n_adp_components, nAtoms, nHklVectors = hkl.size();
            vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));
            nAtoms = mR_at.size();
            
            vector<complex<double> > dTarget_df = _dTarget_df;
            for (int i = 0; i < nHklVectors; i++)
                dTarget_df[i] = conj(dTarget_df[i]);

            // set dTarget_dparam to zero ..

            f.assign(nHklVectors, 0.0);

            dTarget_dparam.resize(mR_at.size());

            for (int atom_index = 0; atom_index < nAtoms ; atom_index++)
            {
                dTarget_dparam[atom_index].adp_derivatives.assign(mAdp[atom_index].size(), 0.0);
                dTarget_dparam[atom_index].atomic_position_derivatives = Vector3d(0, 0, 0);
                dTarget_dparam[atom_index].occupancy_derivatives = 0.0;
            }

            //

			vector < vector<complex<double> > > atomic_ff(nSymmOps, vector < complex<double> >(nAtoms));
			

            for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
            {

				//cout << _hkl[hklIndex][0] << " " << _hkl[hklIndex][1] << " " << _hkl[hklIndex][2] << "\n";

                dTargetDf = dTarget_df[hklIndex];

                hVectorLength2 = hkl[hklIndex] * hkl[hklIndex];
                hVectorLength = sqrt(hVectorLength2);
                            
                pre_atom_loop_sf_calc( /*  in: */  mSymOps, hkl[hklIndex],
                                       /* out: */ rotated_h, translation_factor, adpMultipliers);

				for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
					mManager->calculateCart(rotated_h[symOpIdx], atomic_ff[symOpIdx], countAtomContribution);

                realFContrib = 0;
                imagFContrib = 0;


                for (int atomIdx = 0; atomIdx <nAtoms; atomIdx++)
                {

                    //cout << endl << endl << mCrystal.atoms[atomIdx].label << endl;

                    if (!countAtomContribution[atomIdx])
                        continue;


                    atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];
                    //anomalousScattering = wfnParams[atomWfnIdx].anomalous_scattering;
                    //fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

                    n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                    if (!n_adp_components==0)
                        adps = &mAdp[atomIdx][0];

                    if (n_adp_components == 6)
                        for (int i = 0; i< 6; i++)
                            adp_derivatives[i] = 0.0;

                    xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

                    unweightedTransformedAtomFF_Sum = 0.0;

                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {

                        const Vector3d & rotated_h_ref = rotated_h[symOpIdx];

                        atomic_phase_factor_phase = two_pi*(rotated_h_ref * mR_at[atomIdx] +
                            translation_factor[symOpIdx]);

                        atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                        atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                        complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


						atom_ff = atomic_ff[symOpIdx][atomIdx];// mManager->calculateCart(atomIdx, rotated_h_ref);


                        if (n_adp_components == 6)
                        {
                            double *multipliers = &adpMultipliers[symOpIdx][0];
                            temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                - multipliers[4] * adps[4] - multipliers[5] * adps[5]);

                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;

                            adp_derivatives[0] -= multipliers[0] * unweightedTransformedAtomFF;
                            adp_derivatives[1] -= multipliers[1] * unweightedTransformedAtomFF;
                            adp_derivatives[2] -= multipliers[2] * unweightedTransformedAtomFF;
                            adp_derivatives[3] -= multipliers[3] * unweightedTransformedAtomFF;
                            adp_derivatives[4] -= multipliers[4] * unweightedTransformedAtomFF;
                            adp_derivatives[5] -= multipliers[5] * unweightedTransformedAtomFF;

                        }
                        else
                            unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

                        unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

                        xyz_derivatives[0] += rotated_h_ref[0] * unweightedTransformedAtomFF;
                        xyz_derivatives[1] += rotated_h_ref[1] * unweightedTransformedAtomFF;
                        xyz_derivatives[2] += rotated_h_ref[2] * unweightedTransformedAtomFF;

                        //cout << atom_ff << " " << unweightedTransformedAtomFF << endl;

                    } // symmetry operations


                    if (n_adp_components == 1)
                    {
                        temperature_factor = exp(-hVectorLength * hVectorLength * (*adps));
                        unweightedTransformedAtomFF_Sum *= temperature_factor;
                        dTarget_dparam[atomIdx].adp_derivatives[0] -= hVectorLength2 * (dTargetDf.real() * unweightedTransformedAtomFF_Sum.real() -
                            dTargetDf.imag() * unweightedTransformedAtomFF_Sum.imag());
                    }
                    else
                        if (n_adp_components == 6)
                            for (int i = 0; i<6; ++i)
                                dTarget_dparam[atomIdx].adp_derivatives[i] += dTargetDf.real()*adp_derivatives[i].real() -
                                dTargetDf.imag()*adp_derivatives[i].imag();


                    dTarget_dparam[atomIdx].occupancy_derivatives += unweightedTransformedAtomFF_Sum.real()*dTargetDf.real() -
                        unweightedTransformedAtomFF_Sum.imag()*dTargetDf.imag();
                    //cout << mCrystal.atoms[atomIdx].label << " " << unweightedTransformedAtomFF_Sum * atomWeight << endl;
                    realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
                    imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

					//cout << atomIdx << " " << unweightedTransformedAtomFF_Sum.real()* atomWeight << " " << unweightedTransformedAtomFF_Sum.imag()* atomWeight << "\n";
                    
                    n_adp_components == 1 ? aux = temperature_factor * dTargetDf : aux = dTargetDf;

                    for (int i = 0; i<3; ++i)
                        dTarget_dparam[atomIdx].atomic_position_derivatives[i] -= aux.real()*xyz_derivatives[i].imag() +
                        aux.imag()*xyz_derivatives[i].real();

                } // symmetrically independent atoms

                f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);
                //cin >> hklIndex;
            } // h vectors

            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++)
            {
                atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];
                dTarget_dparam[atomIdx].occupancy_derivatives *= mAtomicMultiplicityWeight[atomIdx];
                multiplier = two_pi_sqare * atomWeight;
                for (int i = 0; i<dTarget_dparam[atomIdx].adp_derivatives.size(); i++)
                    dTarget_dparam[atomIdx].adp_derivatives[i] *= multiplier;
                multiplier = two_pi * atomWeight;
                dTarget_dparam[atomIdx].atomic_position_derivatives.x *= multiplier;
                dTarget_dparam[atomIdx].atomic_position_derivatives.y *= multiplier;
                dTarget_dparam[atomIdx].atomic_position_derivatives.z *= multiplier;
            }

            convertDerivatives(dTarget_dparam);
        }




        void AnyScattererStructureFactorCalculator2::pre_atom_loop_sf_calc(
            //in:
            const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
            const Vector3<double> &hVector,
            //out:
            vector<Vector3<double> > &rotated_h,
            std::vector<double> &translation_factor,
            std::vector<std::vector<double> > &adp_multipliers)
        {
            for (int symmOpIdx = 0; symmOpIdx< symOps.size(); symmOpIdx++)
            {
                translation_factor[symmOpIdx] = hVector*symOps[symmOpIdx].translation;
                rotated_h[symmOpIdx] = hVector*symOps[symmOpIdx].rotation;
                

                // sets mAdpMultipliers
                Vector3<double> &h = rotated_h[symmOpIdx];
                double *adpMultipliers = &adp_multipliers[symmOpIdx][0];

                adpMultipliers[0] = h.x*h.x;
                adpMultipliers[1] = h.y*h.y;
                adpMultipliers[2] = h.z*h.z;
                adpMultipliers[3] = 2.0*h.x*h.y;
                adpMultipliers[4] = 2.0*h.x*h.z;
                adpMultipliers[5] = 2.0*h.y*h.z;
            }

        }


        void AnyScattererStructureFactorCalculator2::convertDerivatives(
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam)
            const
        {
            int atomIdx, nAtoms = dTarget_dparam.size();
            vector<complex<double> > adpIn(6), adpOut(6);
            Vector3<complex<double> > xyzIn, xyzOut;
            int i;

            if ( mAdpConvention != spc::AdpConvention::U_cart)
            {
                for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
                    if (dTarget_dparam[atomIdx].adp_derivatives.size() == 6)
                    {
                        for (i = 0; i<6; i++)
                            adpIn[i] = dTarget_dparam[atomIdx].adp_derivatives[i];
                        mConverter.convertDerivativesADP(adpIn, adpOut, spc::AdpConvention::U_cart, mAdpConvention);
                        for (i = 0; i<6; i++)
                            dTarget_dparam[atomIdx].adp_derivatives[i] = adpOut[i].real();
                    }
            }
            if (mXyzConvention != spc::XyzCoordinateSystem::cartesian)
                for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
                {
                    for (i = 0; i<3; i++)
                        xyzIn[i] = dTarget_dparam[atomIdx].atomic_position_derivatives[i];
                    mConverter.convertDerivativesXyz(xyzIn, xyzOut, spc::XyzCoordinateSystem::cartesian, mXyzConvention);
                    for (i = 0; i<3; i++)
                        dTarget_dparam[atomIdx].atomic_position_derivatives[i] = xyzOut[i].real();
                }
        }

        void AnyScattererStructureFactorCalculator2::convertDerivatives(
            discamb::SfDerivativesAtHkl &derivatives)
            const
        {
            int atomIdx, nAtoms = derivatives.occupancyDerivatives.size();
            vector<complex<double> > adpIn(6), adpOut(6);
            Vector3<complex<double> > xyzIn, xyzOut;
            int i;

            if (mAdpConvention != spc::AdpConvention::U_cart)
            {
                for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
                    if (derivatives.adpDerivatives[atomIdx].size() == 6)
                    {
                        for (i = 0; i<6; i++)
                            adpIn[i] = derivatives.adpDerivatives[atomIdx][i];
                        mConverter.convertDerivativesADP(adpIn, adpOut, spc::AdpConvention::U_cart, mAdpConvention);
                        derivatives.adpDerivatives[atomIdx] = adpOut;
                        //for (i = 0; i<6; i++)
                          //  dTarget_dparam[atomIdx].adp_derivatives[i] = adpOut[i].real();
                    }
            }
            if (mXyzConvention != spc::XyzCoordinateSystem::cartesian)
                for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
                {
                    xyzIn = derivatives.atomicPostionDerivatives[atomIdx];

                    //mConverter.convertDerivativesXyz(xyzIn, xyzOut, spc::cartesian, mXyzConvention);
                    mConverter.convertDerivativesXyz(xyzIn, derivatives.atomicPostionDerivatives[atomIdx], spc::XyzCoordinateSystem::cartesian, mXyzConvention);
                    //for (i = 0; i<3; i++)
                    //    dTarget_dparam[atomIdx].atomic_position_derivatives[i] = xyzOut[i].real();
                }

        }

        void AnyScattererStructureFactorCalculator2::calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &structureFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution)
        {
            Vector3d hklCart;

            //mCrystal
            mReciprocalLatticeUnitCell.fractionalToCartesian(hkl, hklCart);
            
            //--

            const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
            const REAL two_pi = 2 * REAL(M_PI);
            const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
            REAL temperature_factor, hVectorLength, hVectorLength2, multiplier;
            REAL const *adps;
            vector<Vector3<REAL> > rotated_h(mSymOps.size());

            vector<REAL> translation_factor(mSymOps.size());
            complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum, weightedTransformedAtomFF_Sum;

            complex<REAL> adp_derivatives[6];
            complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
            REAL realFContrib, imagFContrib;


            REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];

            complex<REAL> anomalousScattering, atom_ff, aux;


            //--

            int nSymmOps = mSymOps.size();
            int n_adp_components, nAtoms;
            vector<vector<REAL> > adpMultipliers(nSymmOps, vector<double>(6));
            nAtoms = mR_at.size();

            // set dTarget_dparam to zero ..

            
            derivatives.adpDerivatives.resize(nAtoms);
            derivatives.atomicPostionDerivatives.resize(nAtoms);
            derivatives.occupancyDerivatives.resize(nAtoms);
			
			//cout << " ------- \n";
            
			for (int atom_index = 0; atom_index < nAtoms; atom_index++)
            {

                derivatives.adpDerivatives[atom_index].assign(mAdp[atom_index].size(), 0.0);
                derivatives.atomicPostionDerivatives[atom_index] = Vector3d(0, 0, 0);
                derivatives.occupancyDerivatives[atom_index] = 0.0;
            }

            //

            hVectorLength2 = hklCart * hklCart;
            hVectorLength = sqrt(hVectorLength2);

            pre_atom_loop_sf_calc( /*  in: */  mSymOps, hklCart,
                    /* out: */ rotated_h, translation_factor, adpMultipliers);


			// -- precalculate form factors

			for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
			{
				mManager->calculateCart(rotated_h[symOpIdx], mAtomsFormFactor, countAtomContribution);
				
				for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
					mAtomsFormFactors[atomIdx][symOpIdx] = mAtomsFormFactor[atomIdx];
			}


			// --

            realFContrib = 0;
            imagFContrib = 0;


            for (int atomIdx = 0; atomIdx <nAtoms; atomIdx++)
            {

                if (!countAtomContribution[atomIdx])
                    continue;


                atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];

                n_adp_components = mAdp[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                if (!n_adp_components == 0)
                    adps = &mAdp[atomIdx][0];

                if (n_adp_components == 6)
                    for (int i = 0; i< 6; i++)
                        adp_derivatives[i] = 0.0;

                xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

                unweightedTransformedAtomFF_Sum = 0.0;

				const vector<complex<double> >& atom_ffs = mAtomsFormFactors[atomIdx];

                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                {

                    const Vector3d & rotated_h_ref = rotated_h[symOpIdx];

                    atomic_phase_factor_phase = two_pi*(rotated_h_ref * mR_at[atomIdx] +
                        translation_factor[symOpIdx]);

                    atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                    atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                    complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


                    //atom_ff = mManager->calculateCart(atomIdx, rotated_h_ref) + mAnomalous[atomIdx];
					atom_ff = atom_ffs[symOpIdx] + mAnomalous[atomIdx];


                    if (n_adp_components == 6)
                    {
                        double *multipliers = &adpMultipliers[symOpIdx][0];
                        temperature_factor = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                            - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                            - multipliers[4] * adps[4] - multipliers[5] * adps[5]);

                        unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor * temperature_factor;

                        adp_derivatives[0] -= multipliers[0] * unweightedTransformedAtomFF;
                        adp_derivatives[1] -= multipliers[1] * unweightedTransformedAtomFF;
                        adp_derivatives[2] -= multipliers[2] * unweightedTransformedAtomFF;
                        adp_derivatives[3] -= multipliers[3] * unweightedTransformedAtomFF;
                        adp_derivatives[4] -= multipliers[4] * unweightedTransformedAtomFF;
                        adp_derivatives[5] -= multipliers[5] * unweightedTransformedAtomFF;

                    }
                    else
                        unweightedTransformedAtomFF = atom_ff * atomic_position_phase_factor;

                    unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;

                    xyz_derivatives[0] += rotated_h_ref[0] * unweightedTransformedAtomFF;
                    xyz_derivatives[1] += rotated_h_ref[1] * unweightedTransformedAtomFF;
                    xyz_derivatives[2] += rotated_h_ref[2] * unweightedTransformedAtomFF;


                } // symmetry operations


                if (n_adp_components == 1)
                {
                    temperature_factor = exp(-hVectorLength * hVectorLength * (*adps));
                    unweightedTransformedAtomFF_Sum *= temperature_factor;
                    //dTarget_dparam[atomIdx].adp_derivatives[0] -= hVectorLength2 * (dTargetDf.real() * unweightedTransformedAtomFF_Sum.real() -
                      //  dTargetDf.imag() * unweightedTransformedAtomFF_Sum.imag());
                    derivatives.adpDerivatives[atomIdx][0] = -hVectorLength2 * unweightedTransformedAtomFF_Sum;
                }
                else
                    if (n_adp_components == 6)
                        for (int i = 0; i < 6; ++i)
                            derivatives.adpDerivatives[atomIdx][i] = adp_derivatives[i];
                            //dTarget_dparam[atomIdx].adp_derivatives[i] += dTargetDf.real()*adp_derivatives[i].real() -
                            //dTargetDf.imag()*adp_derivatives[i].imag();

                derivatives.occupancyDerivatives[atomIdx] = unweightedTransformedAtomFF_Sum;
                //dTarget_dparam[atomIdx].occupancy_derivatives += unweightedTransformedAtomFF_Sum.real()*dTargetDf.real() -
                  //  unweightedTransformedAtomFF_Sum.imag()*dTargetDf.imag();

                realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
                imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;
                //cout << atomIdx << " " << unweightedTransformedAtomFF_Sum.real() * atomWeight
                //     << " " << unweightedTransformedAtomFF_Sum.imag() * atomWeight << endl;
                n_adp_components == 1 ? aux = temperature_factor : aux = 1;

                for (int i = 0; i < 3; ++i)
                    derivatives.atomicPostionDerivatives[atomIdx][i] += xyz_derivatives[i] * aux;
                    //dTarget_dparam[atomIdx].atomic_position_derivatives[i] -= aux.real()*xyz_derivatives[i].imag() +
                    //aux.imag()*xyz_derivatives[i].real();

            } // symmetrically independent atoms
			//cout << " ------- \n";
            structureFactor = complex<REAL>(realFContrib, imagFContrib);

            

            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++)
            {
                atomWeight = mCrystal.atoms[atomIdx].occupancy * mAtomicMultiplicityWeight[atomIdx];
                //dTarget_dparam[atomIdx].occupancy_derivatives *= mAtomicMultiplicityWeight[atomIdx];
                derivatives.occupancyDerivatives[atomIdx] *= mAtomicMultiplicityWeight[atomIdx];
                multiplier = two_pi_sqare * atomWeight;
                
                //for (int i = 0; i<dTarget_dparam[atomIdx].adp_derivatives.size(); i++)
                for (int i = 0; i<derivatives.adpDerivatives[atomIdx].size(); i++)
                    derivatives.adpDerivatives[atomIdx][i] *= multiplier;
                complex<double> xyz_multiplier = two_pi_i * atomWeight;
                derivatives.atomicPostionDerivatives[atomIdx].x *= xyz_multiplier;
                derivatives.atomicPostionDerivatives[atomIdx].y *= xyz_multiplier;
                derivatives.atomicPostionDerivatives[atomIdx].z *= xyz_multiplier;

                //dTarget_dparam[atomIdx].atomic_position_derivatives.x *= multiplier;
                //dTarget_dparam[atomIdx].atomic_position_derivatives.y *= multiplier;
                //dTarget_dparam[atomIdx].atomic_position_derivatives.z *= multiplier;
            }

            convertDerivatives(derivatives);
            
        }



} // namespace discamb

