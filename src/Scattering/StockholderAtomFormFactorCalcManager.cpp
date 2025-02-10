//#include "discamb/BasicChemistry/chemical_element_data.h"
#include "discamb/MathUtilities/SphericalHarmonics.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include "discamb/BasicUtilities/constants.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/config.h"
//#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/IO/mol2_io.h"
//#include "discamb/IO/vtk_io.h"
#include "discamb/IO/wfn_io.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/MathUtilities/lebedev_laikov.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/MathUtilities/radial_grid.h"
#include "discamb/MathUtilities/real_spherical_harmonics.h"
//#include "discamb/QuantumChemistry/electric_multipoles.h"
#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/QuantumChemistry/HirshfeldPartition.h"
#include "discamb/QuantumChemistry/MbisPartition.h"
#include "discamb/QuantumChemistry/PromoleculeElectronDensity.h"
#include "discamb/Scattering/StockholderAtomFormFactorCalcManager.h"
#include "discamb/Scattering/taam_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"



#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <cmath>
#include <filesystem>

#include <omp.h>

using namespace std;
using namespace discamb::string_utilities;



namespace {

	double delta(int i, int j)
	{
		if (i == j)
			return 1.0;
		return 0.0;
	}


    void removeFileIfExists(
        const std::string &fileName)
    {
		if(filesystem::exists(fileName))
			filesystem::remove(fileName);
    }

    bool splitAtomAndSymmOp(
        const string& s,
        string& atomLabel,
        string& symmOp)
    {
        vector<string> words;
        discamb::string_utilities::split(s,words,',');
        if (words.size() != 4)
            return false;
        
        atomLabel = words[0];
        symmOp = s.substr(words[0].size() + 1);
        return true;
    }

	bool isIdentity(
		const discamb::Matrix3d& m)
	{
		double d = fabs(1.0 - m(0, 0)) + fabs(0.0 - m(0, 1)) + fabs(0.0 - m(0, 2)) +
			       fabs(0.0 - m(1, 0)) + fabs(1.0 - m(1, 1)) + fabs(0.0 - m(1, 2)) +
			       fabs(0.0 - m(2, 0)) + fabs(0.0 - m(2, 1)) + fabs(1.0 - m(2, 2));
		return (d < 1e-8);
	}

	discamb::Vector3d dipoleMoment(
		const std::vector<double>& p1,
		const discamb::Matrix3d& lcs,
		int n_l,
		double kappa,
		double zeta)
	{
		discamb::Vector3d result, p_global = lcs * discamb::Vector3d(p1);
		result = (-1.0) * p_global * double(4.0 * (n_l + 3.0) / (kappa * zeta * 3.0));
		return result;
	}



}

namespace discamb {

	StockholderAtomFormFactorCalcManager::StockholderAtomFormFactorCalcManager(
		const Crystal& crystal,
        const HirshfeldAtomModelSettings& settings)
	{
        mSettings = settings;

#ifndef HAS_SPH_BESSEL
        if (mSettings.sphericalHarmonicsExpansionLevel >= 0)
            on_error::throwException("this compilation does not support multipole expansion for HAR", __FILE__, __LINE__);
#endif


        mUseSpehricalDensitySubstraction = false;
		mCrystal = crystal;
		mReciprocalSpaceUnitCell.set(mCrystal.unitCell);
        mJobName = "noname";
        setElementalIntegrationGrids();
        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            mAtomLabel2IdxInAsymmUnit[crystal.atoms[atomIdx].label] = atomIdx;

		mWfnGeneratorRunner = shared_ptr< WaveFunctionDataGeneratorRunner>(WaveFunctionDataGeneratorRunner::create(settings.wfnCalculation.qmProgram));
        mWfnGeneratorRunner->set(settings.wfnCalculation.qmProgramSpecificData);
		crystal_structure_utilities::atomicNumbers(crystal, mAtomicNumbers);
        
        setSubsystems(settings.crystalFragments, settings.representatives);
        setDistributedMultipoleCenters();
        printDiagnosticInfo();
	}

    void StockholderAtomFormFactorCalcManager::printDiagnosticInfo()
    {
        if (mSettings.diagnostics.printSubsystemsToXyz)
            for (int i = 0; i < mSubsystemAtomicNumbers.size(); i++)
                xyz_io::writeXyz(mSettings.crystalFragments[i].label + ".xyz", mSubsystemAtomicNumbers[i], mSubsystemAtomicPositions[i]);


        if (mSettings.diagnostics.printSubsystemsToMol2)
            for (int i = 0; i < mSubsystemAtomicNumbers.size(); i++)
                mol2_io::write(mSettings.crystalFragments[i].label + ".mol2", mSubsystemAtomicNumbers[i], mSubsystemAtomicPositions[i], true);

        if (mSettings.diagnostics.printMultipoleClustersToXyz || mSettings.diagnostics.printMultipoleClusterToMol2)
            for (int i = 0; i < mSubsystemAtomicNumbers.size(); i++)
            {
                vector<Vector3d> positions;
                vector<string> symbols;
                vector<int> atomicNumbers;
                int n = mMultipoleBearingAtoms[i].size();
                for (auto& center : mMultipoleBearingAtoms[i])
                {
                    positions.push_back(crystal_structure_utilities::atomPosition(center.first, center.second, mCrystal));
                    symbols.push_back(mCrystal.atoms[center.first].type);
                    atomicNumbers.push_back(periodic_table::atomicNumber(symbols.back()));
                }
                if (mSettings.diagnostics.printMultipoleClustersToXyz || mSettings.diagnostics.printMultipoleClusterToMol2)
                    xyz_io::writeXyz(string("multipole_cluster_") + mSettings.crystalFragments[i].label + ".xyz", symbols, positions);
                if (mSettings.diagnostics.printMultipoleClustersToXyz || mSettings.diagnostics.printMultipoleClusterToMol2)
                    mol2_io::write(string("multipole_cluster_")+ mSettings.crystalFragments[i].label + ".mol2", atomicNumbers, positions, false);
            }

    }

    void StockholderAtomFormFactorCalcManager::setDistributedMultipoleCenters()
    {
        mMultipoleBearingAtoms.clear();
        mCustomWeightMultipoles.clear();

        if (!mSettings.multipoleExpansion.calculate)
            return;

        UnitCellContent ucContent(mCrystal);
        int subsystemIdx, nSubsystems = mSettings.crystalFragments.size();
        vector<optional<DistributedMultipoleCentersSettings> > subsystemSpecificDmSettings(nSubsystems);
        if (!mSettings.fragmentWfnCalculation.empty())
            for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
                subsystemSpecificDmSettings[subsystemIdx] = mSettings.fragmentWfnCalculation[subsystemIdx].distributedMultipoleCluster;
        vector<FragmentAtoms> fragmentAtoms;
        for (auto const subsystem : mSettings.crystalFragments)
            fragmentAtoms.push_back(subsystem.atoms);

        vector< vector<pair<string, string > > > multipoleClusterAtoms(nSubsystems);

        distributed_multipoles::get_fragments_multipole_centers(
            subsystemSpecificDmSettings,
            mSettings.multipoleExpansion.clusterThreshold,
            mCrystal,
            fragmentAtoms,
            multipoleClusterAtoms, mCustomWeightMultipoles);

        mMultipoleBearingAtoms.resize(nSubsystems);
        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {
            vector < pair < int, string > > atomList;
            crystal_structure_utilities::convertAtomList(
                mCrystal, multipoleClusterAtoms[subsystemIdx], atomList);
            for (auto atom : atomList)
                mMultipoleBearingAtoms[subsystemIdx].push_back({ atom.first, SpaceGroupOperation(atom.second) });
        }
        

    }

	void StockholderAtomFormFactorCalcManager::setTransforms()
	{
		int nAtoms = mCrystal.atoms.size();
        
		mTransforms.clear(); 
		mTransforms.resize(nAtoms);
		pair<string, string> clusterAtom;
		string symmOpAsString, clusterAtomLabel;
		for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
			for (auto& rep : mSettings.representatives[atomIdx])
			{
                clusterAtom = mSettings.crystalFragments[rep.fragmentIdx].atoms.atomList[rep.idxInSubsystem];
				mTransforms[atomIdx].push_back(shared_ptr<AtomicDensityTransform>(
					AtomicDensityTransform::create(rep.transformType, mCrystal, rep.transformSetting, clusterAtom, rep.atomLabel )));
			}
	}


    void StockholderAtomFormFactorCalcManager::setSubsystems(
        const std::vector<QmFragmentInCrystal>& subsystems,
        const std::vector<std::vector<AtomRepresentativeInfo> >& representatives)
    {
        
    	map<string, int> labelToIndex;
        vector<int> atomicNumbers;

        crystal_structure_utilities::atomicNumbers(mCrystal, atomicNumbers);
        
    	int idx = 0;
    	for (auto& atom : mCrystal.atoms)
    		labelToIndex[atom.label] = idx++;

        int subsystemIdx, nSubsystems = subsystems.size();
        
        mSubsystemAtomicNumbers.clear();
        mSubsystemAtomicNumbers.resize(nSubsystems);
        mSubsystemAtomicPositions.clear(); // angstroms
        mSubsystemAtomicPositions.resize(nSubsystems);

        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {

            for (auto& atom : mSettings.crystalFragments[subsystemIdx].atoms.atomList)
            {
                mSubsystemAtomicNumbers[subsystemIdx].push_back(atomicNumbers[labelToIndex[atom.first]]);
                mSubsystemAtomicPositions[subsystemIdx].push_back(
                    crystal_structure_utilities::atomPosition(atom.first, atom.second, mCrystal));
            }

            for (auto& capH : mSettings.crystalFragments[subsystemIdx].atoms.cappingHydrogens)
            {
                mSubsystemAtomicNumbers[subsystemIdx].push_back(1);
                mSubsystemAtomicPositions[subsystemIdx].push_back(fragmentation::capping_h_position(mCrystal, capH));
            }
        }

    	setElementalIntegrationGrids();

        setAndPrintRepresentatives(representatives);
        setTransforms();
        allocateElectronDensityContainers();
    }

	void StockholderAtomFormFactorCalcManager::setAndPrintRepresentatives(
		const std::vector < std::vector<AtomRepresentativeInfo> >& representatives) 
    {
        int nAtoms, subsystemIdx, nSubsystems =  mSettings.crystalFragments.size();
        
		nAtoms = mCrystal.atoms.size();
		
		mListOfRepresentativesInSubsystem.resize(nSubsystems);
		mAtomsInSubsystemWhichAreRepresentatives.resize(nSubsystems);

		for (int atomIdx = 0; atomIdx < mSettings.representatives.size(); atomIdx++)
			for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
                mListOfRepresentativesInSubsystem[mSettings.representatives[atomIdx][i].fragmentIdx].push_back({ atomIdx,i });
		
		mSubsystemAtomsWhichAreRepresentatives.resize(nSubsystems);
		for (auto &representatives: mSettings.representatives)
			for (auto &representative: representatives)
                mSubsystemAtomsWhichAreRepresentatives[representative.fragmentIdx].insert(representative.idxInSubsystem);


        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {
            string fileName = string("rep_atoms_") + to_string(subsystemIdx + 1);
            ofstream out(fileName);


            out << mSettings.representatives.size() << endl;
				

            for (int atomIdx = 0; atomIdx < mSettings.representatives.size(); atomIdx++)
            {
                int nRepresentatives = 0;
                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
					if (mSettings.representatives[atomIdx][i].fragmentIdx == subsystemIdx)
						nRepresentatives++;
					

                out << nRepresentatives << endl;
                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
					if (mSettings.representatives[atomIdx][i].fragmentIdx == subsystemIdx)
					{
						out << mSettings.representatives[atomIdx][i].idxInSubsystem << endl;
						mAtomsInSubsystemWhichAreRepresentatives[subsystemIdx].push_back(mSettings.representatives[atomIdx][i].idxInSubsystem);
					}
            }
            out.close();
        }
	}





	void StockholderAtomFormFactorCalcManager::setPositions(
		const Crystal& crystal)
	{
		int atomIdx, nAtoms, capH_Idx, nCapH, nSubsystems, subsystemIdx;

        nSubsystems = mSettings.crystalFragments.size();

        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {
            auto& atoms = mSettings.crystalFragments[subsystemIdx].atoms.atomList;
            nAtoms = atoms.size(); 
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                mSubsystemAtomicPositions[subsystemIdx][atomIdx] = crystal_structure_utilities::atomPosition( mAtomLabel2IdxInAsymmUnit[atoms[atomIdx].first], atoms[atomIdx].second, crystal);
            auto& capH = mSettings.crystalFragments[subsystemIdx].atoms.cappingHydrogens;
            nCapH = capH.size();
            for(capH_Idx=0; capH_Idx<nCapH; capH_Idx++)
                mSubsystemAtomicPositions[subsystemIdx][nAtoms+capH_Idx] = fragmentation::capping_h_position(crystal, capH[capH_Idx]);

        }

	}

	StockholderAtomFormFactorCalcManager::~StockholderAtomFormFactorCalcManager()
	{

	}

	void StockholderAtomFormFactorCalcManager::update(
		const std::vector<AtomInCrystal>& atoms)
	{
		mCrystal.atoms = atoms;
		setPositions(mCrystal);
	}


	void StockholderAtomFormFactorCalcManager::calculateMultipoles()
	{
		vector<int> atomicNumbers;
		crystal_structure_utilities::atomicNumbers(mCrystal, atomicNumbers);
		mMultipoles.assign(mCrystal.atoms.size(), ElectricMultipoles());
		for (int atomIdx = 0; atomIdx < mCrystal.atoms.size(); atomIdx++)
			calculateMultipoles(atomIdx, atomicNumbers[atomIdx], mMultipoles[atomIdx]);

	}

	void StockholderAtomFormFactorCalcManager::calculateMultipoles(
		int atomIdx, 
		int z,
		ElectricMultipoles& multipolesExpansion)
	{

		Matrix3d densityTransform, hTransform;
		double weight, weightsSum = 0;
		Vector3d hklCart, hklRotated, r_rotated;
		int subsystemIdx, atomInSubsystemIdx;
		vector<Vector3d> points;
		ElectricMultipoles representativeMultipoles;// , multipoles;
		vector<double> weightedCharges;

		multipolesExpansion.charge = 0;
		multipolesExpansion.dipole.set(0, 0, 0);
		multipolesExpansion.quadrupole.set(0, 0, 0, 0, 0, 0, 0, 0, 0);
		
		

		for (int representativeIdx = 0; representativeIdx < mSettings.representatives[atomIdx].size(); representativeIdx++)
		{
			densityTransform = mTransforms[atomIdx][representativeIdx]->getTransformMatrix(mCrystal);

            subsystemIdx = mSettings.representatives[atomIdx][representativeIdx].fragmentIdx;
			atomInSubsystemIdx = mSettings.representatives[atomIdx][representativeIdx].idxInSubsystem;
			


			vector<Vector3d> points = mElementGridPoints[z];
			const vector<double>& density = mAtomInClusterElectronDensity[subsystemIdx][atomInSubsystemIdx];
			const vector<double>& weights = mElementIntegrationWeights[z];

			for (auto& r : points)
			{
				r_rotated = densityTransform * r;
				r = r_rotated/constants::Angstrom;
			}

			//#######################


			int pointIdx, nIntegrationPoints = density.size();
			weightedCharges = density;

			for (pointIdx = 0; pointIdx < nIntegrationPoints; pointIdx++)
				weightedCharges[pointIdx] *= -weights[pointIdx];

            distributed_multipoles::calculate_multipoles(
                points, 
                weightedCharges, 
                representativeMultipoles, 
                mSettings.multipoleExpansion.multipoleExpansionLevel);

            if (!mUseSpehricalDensitySubstraction)
                representativeMultipoles.charge += z;

			if (mSettings.representatives[atomIdx][representativeIdx].isWeightFixed)
				weight = mSettings.representatives[atomIdx][representativeIdx].fixedWeightValue;
			else
			{
				weight = 0;
				on_error::not_implemented(__FILE__, __LINE__);
			}

			multipolesExpansion.charge += weight * representativeMultipoles.charge;
			multipolesExpansion.dipole += weight * representativeMultipoles.dipole;
			multipolesExpansion.quadrupole += weight * representativeMultipoles.quadrupole;

			weightsSum += weight;
		}

		multipolesExpansion.charge /= weightsSum;
		multipolesExpansion.dipole /= weightsSum;
		multipolesExpansion.quadrupole /= weightsSum;

	}


    void StockholderAtomFormFactorCalcManager::calculateFracWithSphericalHarmonics(
        const std::vector<Vector3d>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom) const
    {
#ifndef HAS_SPH_BESSEL
        if (mSettings.sphericalHarmonicsExpansionLevel >= 0)
            on_error::throwException("this compilation does not support multipole expansion for HAR", __FILE__, __LINE__);
#endif

        
        cout << "calculateFracWithSphericalHarmonics\n";
        //constexpr int maxL = 7;
        int maxL = mSettings.sphericalHarmonicsExpansionLevel;

        int nAtoms = includeAtom.size();
        int hklIdx, nHkl = hkl.size();
        formFactors.assign(nHkl, vector<complex<double> >(nAtoms, 0.0));
        // atom,l,l+m,r
        vector<vector<vector<vector<double> > > > radialFunctions;
        //int nAng = 590;
        //int nRad = 48;

        int nAng = mSettings.formFactorsCalculation.integrationGrid.angularGridSize;
        int nRad = mSettings.formFactorsCalculation.integrationGrid.radialGridSize;


        //-------------###########
        // z, point
        vector<vector<double> > radialGrid(113);
        vector<vector<double> > radialWeights(113);
        vector<Vector3d> angularGrid;
        vector<double> angularWeights;

        std::set<int> uniqueZ;

        for (auto& clusterZ : mSubsystemAtomicNumbers)
            uniqueZ.insert(clusterZ.begin(), clusterZ.end());

        vector<Vector3d> angularGridPoints, gridPoints;
        Vector3d r, atomPosition;

        lebedev_laikov::get_grid(nAng, angularGrid, angularWeights);


        for (auto z : uniqueZ)
            radial_grid::mura_knowles(z, nRad, radialGrid[z], radialWeights[z]);

        
        vector< vector<vector<double> > > sph(nAng);

        for (int i=0;i<nAng; i++)
        {
            sph[i].resize(maxL + 1, vector<double>(2*maxL+1));
            MChLib::SphericalHarmonics::calculate(maxL, angularGrid[i].x, angularGrid[i].y, angularGrid[i].z, sph[i]);
            //real_spherical_harmonics::getWfnNormalized<maxL>(angularGrid[i], sph[i]);
            
        }



        
        radialFunctions.resize(nAtoms);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (includeAtom[atomIdx])
            {
                radialFunctions[atomIdx].resize(maxL + 1, vector<vector<double> >(2 * maxL + 1, vector<double>(nRad,0.0)));
                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
                {
                    if (i > 0)
                        on_error::throwException("it works only for 1 representative - coordination systems adjustment is not introduced", __FILE__, __LINE__);

                    int subsystemIdx = mSettings.representatives[atomIdx][i].fragmentIdx;
                    int atomInSubsystemIdx = mSettings.representatives[atomIdx][i].idxInSubsystem;

                    const vector<double>& density = mAtomInClusterElectronDensity[subsystemIdx][atomInSubsystemIdx];

                    for (int l = 0; l <= maxL; l++)
                        for (int m = -l; m <= l; m++)
                        {
                            int idx = l + m;
                            for (int radGridIdx = 0; radGridIdx < nRad; radGridIdx++)
                            {
                                double d = 0;
                                for (int angGridIdx = 0; angGridIdx < nAng; angGridIdx++)
                                {
                                    double density_value = density[angGridIdx * nRad + radGridIdx];
                                    d += 4 * M_PI * angularWeights[angGridIdx] * sph[angGridIdx][l][idx] * density_value;
                                    //d += 4 * M_PI * angularWeights[angGridIdx] * sph[angGridIdx][l][idx] * density[angGridIdx * nRad + radGridIdx];
                                }
                                radialFunctions[atomIdx][l][idx][radGridIdx] += mSettings.representatives[atomIdx][i].fixedWeightValue * d;
                            }
                        }

                }
            }
        
        //++
        
        ofstream out("rad");
        for (int radGridIdx = 0; radGridIdx < nRad; radGridIdx++)
        {
            out << radGridIdx << " ";
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                out << radialFunctions[atomIdx][0][0][radGridIdx] << " ";
            out << "\n";
        }

        out.close();

        // 1, i, -1, -1*i,1, i, -1, -1*i,1
        const complex<double> im(0.0, 1.0);
        const vector<complex<double> > i_power_l = { 1.0, im, -1.0, -im , 1.0, im, -1.0, -im , 1.0, im };
        vector<vector<double> > realSph_of_hkl(maxL+1);
        for (int i = 0; i <= maxL; i++)
            realSph_of_hkl[i].resize(2 * i + 1);

        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {



            Vector3d hklCart, twoPiHklCartAu; 

            
            mReciprocalSpaceUnitCell.fractionalToCartesian(hkl[hklIdx], hklCart);

            twoPiHklCartAu = 2.0 * M_PI * hklCart / constants::Angstrom;
            
            double twoPiS = sqrt(twoPiHklCartAu * twoPiHklCartAu);

            {
                Vector3d hklNormalized;
                if (hkl[hklIdx] == Vector3d(0,0,0))
                    hklNormalized.set(1, 0, 0);
                else
                    hklNormalized = hklCart / sqrt(hklCart * hklCart);
                //real_spherical_harmonics::getWfnNormalized<maxL>(hklNormalized, realSph_of_hkl);
                MChLib::SphericalHarmonics::calculate(maxL, hklNormalized.x, hklNormalized.y, hklNormalized.z, realSph_of_hkl);
            }



            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                double weightsSum = 0;
                Vector3d twoPiHklRotated;
                Vector3d hklNormalized;
                
                int z = mAtomicNumbers[atomIdx];

                formFactors[hklIdx][atomIdx] = 0.0;
                
                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
                {
                    twoPiHklRotated = twoPiHklCartAu;
                    hklNormalized = twoPiHklRotated / sqrt(twoPiHklRotated * twoPiHklRotated);

                    
                    complex<double> f=0.0;

                    for (int l = 0; l <= maxL; l++)
                        for (int m = -l; m <= l; m++)
                        {
                            int l_plus_m = l + m;
                            double fourierBesselTransform = 0;
                            double x = 0;

                            for (int radGridIdx = 0; radGridIdx < nRad; radGridIdx++)
                            {
                                double twoPiSr = twoPiS * radialGrid[z][radGridIdx];
#ifdef HAS_SPH_BESSEL
                                fourierBesselTransform += sph_bessel(l, twoPiSr) * radialGrid[z][radGridIdx] * radialGrid[z][radGridIdx] *
                                    radialFunctions[atomIdx][l][l_plus_m][radGridIdx] * radialWeights[z][radGridIdx];
#endif
                            }
                            f += 4.0 * M_PI * fourierBesselTransform * i_power_l[l] * realSph_of_hkl[l][l_plus_m];
                        }
                    formFactors[hklIdx][atomIdx] += mSettings.representatives[atomIdx][i].fixedWeightValue * f;

                }

                if (mUseSpehricalDensitySubstraction)
                    formFactors[hklIdx][atomIdx] += mSphericalFormFactors[mAtomToSphericalTypeMap[atomIdx]].calculate(sqrt(hklCart * hklCart));
            }
        }
        out.open("frad");
        for (int i = 1; i <= 50; i++)
        {
            double step = 0.05; 
            double twoPiS = 2.0 * M_PI * step * i / constants::Angstrom;
            out << setw(14) << setprecision(6) << fixed << step * i;
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                int z = mAtomicNumbers[atomIdx];

                

                for (int l = 0; l <= maxL; l++)
                    for (int m = -l; m <= l; m++)
                    {
                        int l_plus_m = l + m;
                        double fourierBesselTransform = 0;

                        for (int radGridIdx = 0; radGridIdx < nRad; radGridIdx++)
                        {
                            double twoPiSr = twoPiS * radialGrid[z][radGridIdx];
#ifdef HAS_SPH_BESSEL
                            fourierBesselTransform += sph_bessel(l, twoPiSr) * radialGrid[z][radGridIdx] * radialGrid[z][radGridIdx] *
                                radialFunctions[atomIdx][l][l_plus_m][radGridIdx] * radialWeights[z][radGridIdx];
#endif
                        }
                        
                        out << setw(14) << setprecision(6) << fourierBesselTransform;

                    }
            }
            out << "\n";
        }
        out.close();
    }


    void StockholderAtomFormFactorCalcManager::calculateFrac(
        const std::vector<Vector3d>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
     
        if (mSettings.sphericalHarmonicsExpansionLevel >= 0)
        {
            calculateFracWithSphericalHarmonics(hkl, formFactors, includeAtom);
            return;
        }

        int nAtoms = includeAtom.size();
        int nHkl = hkl.size();
        int nThreads = omp_get_max_threads();

        //formFactors[hkl idx][atom idx]
        formFactors.assign(nHkl, vector<complex<double> >(nAtoms, 0.0));

        // set transforms

        vector<vector<Matrix3d> > hTransforms(nAtoms);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (includeAtom[atomIdx])
            {
                hTransforms[atomIdx].resize(mSettings.representatives[atomIdx].size());
                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
                {
                    auto densityTransform = mTransforms[atomIdx][i]->getTransformMatrix(mCrystal);
                    hTransforms[atomIdx][i] = algebra3d::transpose3d(densityTransform);
                }
            }

        // allocate weightedCosValues/weightedSinValues - cos(r)*integrationWeight(r)
        // [thread][z][point]
        vector< vector<vector<double> > > weightedCosValues(nThreads, vector<vector<double> >(113)), weightedSinValues(nThreads, vector<vector<double> >(113));
        std::set<int> uniqueZ(mAtomicNumbers.begin(), mAtomicNumbers.end());

        for (auto z : uniqueZ)
            for (int thread = 0; thread < nThreads; thread++)
            {
                weightedCosValues[thread][z].resize(mElementGridPoints[z].size() / 2);
                weightedSinValues[thread][z].resize(mElementGridPoints[z].size() / 2);
            }

        int nHklPerThread = nHkl / mSettings.hardware.nCores;// mN_Threads;
        int thread0_hkl_processed = 0;
#pragma omp parallel for
        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            Vector3d hklCart, twoPiHklCartAu;
            double realPart, imaginaryPart, phase;
            int threadIdx = omp_get_thread_num();
            int halfN_Points;

            mReciprocalSpaceUnitCell.fractionalToCartesian(hkl[hklIdx], hklCart);

            twoPiHklCartAu = 2.0 * M_PI * hklCart / constants::Angstrom;

            int pointIdx;

            bool useInversionSymmetry = true;

            // fill cosValues, sinValues

            double w;

            for (int z : uniqueZ)
            {
                const vector<Vector3d>& points = mElementGridPoints[z];
                const vector<double>& weights = mElementIntegrationWeights[z];
                vector<double>& cosTable = weightedCosValues[threadIdx][z];
                vector<double>& sinTable = weightedSinValues[threadIdx][z];
                halfN_Points = points.size() / 2;

                for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
                {
                    phase = twoPiHklCartAu * points[pointIdx];
                    w = weights[pointIdx];
                    cosTable[pointIdx] = cos(phase) * w;
                    sinTable[pointIdx] = sin(phase) * w;
                }
            }

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                double weight, weightsSum = 0;
                Vector3d twoPiHklRotated;
                int clusterIdx, atomInClusterIdx, z = mAtomicNumbers[atomIdx];

                for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
                {
                    twoPiHklRotated = hTransforms[atomIdx][i] * twoPiHklCartAu;

                    clusterIdx = mSettings.representatives[atomIdx][i].fragmentIdx;
                    atomInClusterIdx = mSettings.representatives[atomIdx][i].idxInSubsystem;

                    const vector<double>& density = mAtomInClusterElectronDensity[clusterIdx][atomInClusterIdx];
                    halfN_Points = density.size() / 2;

                    realPart = 0;
                    imaginaryPart = 0;


                    if (!isIdentity(hTransforms[atomIdx][i]))
                    {

                        const vector<Vector3d>& points = mElementGridPoints[z];
                        const vector<double>& weights = mElementIntegrationWeights[z];

                        for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
                        {
                            phase = twoPiHklRotated * points[pointIdx];
                            realPart += cos(phase) * weights[pointIdx] * (density[pointIdx] + density[pointIdx + halfN_Points]);
                            imaginaryPart += sin(phase) * weights[pointIdx] * (density[pointIdx] - density[pointIdx + halfN_Points]);
                        }
                    }
                    else
                    {
                        vector<double>& c = weightedCosValues[threadIdx][z];
                        vector<double>& s = weightedSinValues[threadIdx][z];
                        for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
                        {
                            realPart += c[pointIdx] * (density[pointIdx] + density[pointIdx + halfN_Points]);
                            imaginaryPart += s[pointIdx] * (density[pointIdx] - density[pointIdx + halfN_Points]);
                        }

                    }

                    if (mSettings.representatives[atomIdx][i].isWeightFixed)
                        weight = mSettings.representatives[atomIdx][i].fixedWeightValue;
                    else
                    {
                        on_error::not_implemented(__FILE__, __LINE__);
                        weight = 0;
                    }

                    formFactors[hklIdx][atomIdx] += weight * complex<double>(realPart, imaginaryPart);
                    weightsSum += weight;

                }

                formFactors[hklIdx][atomIdx] /= weightsSum;

                if (mUseSpehricalDensitySubstraction)
                    formFactors[hklIdx][atomIdx] += mSphericalFormFactors[mAtomToSphericalTypeMap[atomIdx]].calculate(sqrt(hklCart * hklCart));
            }
        }

    }


	void StockholderAtomFormFactorCalcManager::calculateFrac(
		const std::vector<Vector3i>& hkl,
		std::vector < std::vector<std::complex<double> > >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        

        if (mSettings.sphericalHarmonicsExpansionLevel >= 0)
        {
            vector<Vector3d> hkl_d;
            for (const auto& h : hkl)
                hkl_d.push_back(h);
            calculateFracWithSphericalHarmonics(hkl_d, formFactors, includeAtom);
            return;
        }

		int nAtoms = includeAtom.size();
		int nHkl = hkl.size();
		int nThreads = omp_get_max_threads();
		
		//formFactors[hkl idx][atom idx]
		formFactors.assign(nHkl, vector<complex<double> >(nAtoms, 0.0));

		// set transforms
		
		vector<vector<Matrix3d> > hTransforms(nAtoms);

		for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
			if (includeAtom[atomIdx])
			{
				hTransforms[atomIdx].resize(mSettings.representatives[atomIdx].size());
				for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
				{
					auto densityTransform = mTransforms[atomIdx][i]->getTransformMatrix(mCrystal);
					hTransforms[atomIdx][i] = algebra3d::transpose3d(densityTransform);
				}
			}

		// allocate weightedCosValues/weightedSinValues - cos(r)*integrationWeight(r)
		// [thread][z][point]
		vector< vector<vector<double> > > weightedCosValues(nThreads, vector<vector<double> >(113)), weightedSinValues(nThreads, vector<vector<double> >(113));
		std::set<int> uniqueZ(mAtomicNumbers.begin(), mAtomicNumbers.end());

		for (auto z : uniqueZ)
			for (int thread = 0; thread < nThreads; thread++)
			{
				weightedCosValues[thread][z].resize(mElementGridPoints[z].size()/2);
				weightedSinValues[thread][z].resize(mElementGridPoints[z].size()/2);
			}
        
        int nHklPerThread = nHkl / mSettings.hardware.nCores;// mN_Threads;
        int thread0_hkl_processed = 0;
#pragma omp parallel for
		for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
		{
			Vector3d hklCart, twoPiHklCartAu;
			double realPart, imaginaryPart, phase;
			int threadIdx = omp_get_thread_num();
			int halfN_Points;

			mReciprocalSpaceUnitCell.fractionalToCartesian(hkl[hklIdx], hklCart);

			twoPiHklCartAu = 2.0 * M_PI * hklCart / constants::Angstrom;

			int pointIdx;

			bool useInversionSymmetry = true;

			// fill cosValues, sinValues

			double w;

			for (int z : uniqueZ)
			{
				const vector<Vector3d>& points = mElementGridPoints[z];
				const vector<double>& weights = mElementIntegrationWeights[z];
				vector<double>& cosTable = weightedCosValues[threadIdx][z];
				vector<double>& sinTable = weightedSinValues[threadIdx][z];
				halfN_Points = points.size() / 2;

				for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
				{
					phase = twoPiHklCartAu * points[pointIdx];
					w = weights[pointIdx];
					cosTable[pointIdx] = cos(phase) * w;
					sinTable[pointIdx] = sin(phase) * w;
				}
			}

			for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
			{
				double weight, weightsSum = 0;
				Vector3d twoPiHklRotated;
				int clusterIdx, atomInClusterIdx, z = mAtomicNumbers[atomIdx];

				for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
				{
					//hklRotated = (hTransforms[atomIdx][i] * hklCart) / constants::Angstrom;
					twoPiHklRotated = hTransforms[atomIdx][i] * twoPiHklCartAu;

					clusterIdx = mSettings.representatives[atomIdx][i].fragmentIdx;
					atomInClusterIdx = mSettings.representatives[atomIdx][i].idxInSubsystem;

					const vector<double>& density = mAtomInClusterElectronDensity[clusterIdx][atomInClusterIdx];
					halfN_Points = density.size()/2;

					realPart = 0;
					imaginaryPart = 0;


					if (!isIdentity(hTransforms[atomIdx][i]))
					{

						const vector<Vector3d>& points = mElementGridPoints[z];
						const vector<double>& weights = mElementIntegrationWeights[z];

						for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
						{
							phase = twoPiHklRotated * points[pointIdx];
							realPart += cos(phase) * weights[pointIdx] * (density[pointIdx] + density[pointIdx+ halfN_Points]);
							imaginaryPart += sin(phase) * weights[pointIdx] * (density[pointIdx] - density[pointIdx + halfN_Points]);
						}
					}
					else
					{
						vector<double>& c = weightedCosValues[threadIdx][z];
						vector<double>& s = weightedSinValues[threadIdx][z];
						for (pointIdx = 0; pointIdx < halfN_Points; pointIdx++)
						{
							realPart += c[pointIdx] * (density[pointIdx] + density[pointIdx+ halfN_Points]);
							imaginaryPart += s[pointIdx] * (density[pointIdx] - density[pointIdx + halfN_Points]);
						}

					}

					if (mSettings.representatives[atomIdx][i].isWeightFixed)
						weight = mSettings.representatives[atomIdx][i].fixedWeightValue;
					else
					{
						on_error::not_implemented(__FILE__, __LINE__);
						weight = 0;
					}

					formFactors[hklIdx][atomIdx] += weight * complex<double>(realPart, imaginaryPart);
					weightsSum += weight;

				}

				formFactors[hklIdx][atomIdx] /= weightsSum;

				if (mUseSpehricalDensitySubstraction)
					formFactors[hklIdx][atomIdx] += mSphericalFormFactors[mAtomToSphericalTypeMap[atomIdx]].calculate(sqrt(hklCart * hklCart));
			}
		}       

	}





	void StockholderAtomFormFactorCalcManager::calculateCart(
		const std::vector <Vector3d>& hkl,
		std::vector < std::vector<std::complex<double> > >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
		Vector3d frac;
		vector<Vector3i> hklFracVec;
		Vector3i hklFrac;

		for (auto& h : hkl)
		{
			mReciprocalSpaceUnitCell.cartesianToFractional(h, frac);
			hklFrac[0] = math_utilities::roundInt(frac[0]);
			hklFrac[1] = math_utilities::roundInt(frac[1]);
			hklFrac[2] = math_utilities::roundInt(frac[2]);
			hklFracVec.push_back(hklFrac);
		}


		calculateFrac(hklFracVec, formFactors, includeAtom);

	}


	std::complex<double> StockholderAtomFormFactorCalcManager::calculateFrac(
		int atomIdx,
		const Vector3i& hkl)
		const
	{
		complex<double> result = 0.0;
		double realPart, imaginaryPart;
		Matrix3d densityTransform, hTransform;
		double phase, weight, weightsSum = 0;
		Vector3d hklCart, hklRotated;
		int clusterIdx, atomInClusterIdx, pointIdx, nPoints, z = mAtomicNumbers[atomIdx];
        
        mReciprocalSpaceUnitCell.fractionalToCartesian(hkl, hklCart);
		
        for (int i = 0; i < mSettings.representatives[atomIdx].size(); i++)
		{
			densityTransform = mTransforms[atomIdx][i]->getTransformMatrix(mCrystal);
			hTransform = algebra3d::transpose3d(densityTransform);
			hklRotated = (hTransform * hklCart) / constants::Angstrom;
			
			clusterIdx = mSettings.representatives[atomIdx][i].fragmentIdx;
			atomInClusterIdx = mSettings.representatives[atomIdx][i].fragmentIdx;

			const vector<Vector3d>& points = mElementGridPoints[z];
			const vector<double>& density = mAtomInClusterElectronDensity[clusterIdx][atomInClusterIdx];
			const vector<double>& weights = mElementIntegrationWeights[z];

			nPoints = points.size();

			realPart = 0;
			imaginaryPart = 0;

			for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
			{
				phase = 2.0 * M_PI * hklRotated * points[pointIdx];
				realPart += cos(phase) * density[pointIdx] * weights[pointIdx];
				imaginaryPart += sin(phase) * density[pointIdx] * weights[pointIdx];
			}

			if (mSettings.representatives[atomIdx][i].isWeightFixed)
				weight = mSettings.representatives[atomIdx][i].fixedWeightValue;
			else
				on_error::not_implemented(__FILE__, __LINE__);

			result += weight * complex<double>(realPart, imaginaryPart);
			weightsSum += weight;
		}

		result /= weightsSum;

        if (mUseSpehricalDensitySubstraction)
            result += mSphericalFormFactors[mAtomToSphericalTypeMap[atomIdx]].calculate(sqrt(hklCart * hklCart));

		return result;

	}

    void StockholderAtomFormFactorCalcManager::setSphericalAtomicDensities(
        const ProatomDB& db)
    {
        mUseSpehricalDensitySubstraction = true;

        map<int, int> atomicNumber2TypeIndex;
        std::set<int> uniqueZ;
        vector<int> atomicNumbers;

        crystal_structure_utilities::atomicNumbers(mCrystal, atomicNumbers);
        uniqueZ.insert(atomicNumbers.begin(), atomicNumbers.end());

        int index = 0;
        mSphericalTypeAtomicNumber.clear();
        mSphericalFormFactors.clear();
        mSphericalDensities.resize(uniqueZ.size());

        for (auto& z : uniqueZ)
        {
            atomicNumber2TypeIndex[z] = index;
            mSphericalTypeAtomicNumber.push_back(z);
            db.getSphericalAtom(z, 0, mSphericalDensities[index]);
            mSphericalFormFactors.push_back(NumericalSphericalAtomFormFactor(mSphericalDensities[index]));
            index++;
        }

        mAtomToSphericalTypeMap.clear();
        for (auto z : atomicNumbers)
            mAtomToSphericalTypeMap.push_back(atomicNumber2TypeIndex[z]);
    }

	std::complex<double> StockholderAtomFormFactorCalcManager::calculateCart(
		int atomIdx,
		const Vector3d& hkl)
		const
	{
		Vector3d frac;
		Vector3i hklFrac;
		mReciprocalSpaceUnitCell.cartesianToFractional(hkl, frac);

		hklFrac[0] = math_utilities::roundInt(frac[0]);
		hklFrac[1] = math_utilities::roundInt(frac[1]);
		hklFrac[2] = math_utilities::roundInt(frac[2]);

		return calculateFrac(atomIdx, hklFrac);
	}




	// for each atom q_charge, q1_dipol, q2_dipol, q1_quadrupole, ..
	// and the same with positions
	// positions with respect to atom position
	void StockholderAtomFormFactorCalcManager::calculateAsymmetricUnitPointChargeMultipoleRepresentation(
		std::vector< std::vector<double> > &chargeAsymm,
		std::vector< std::vector<Vector3d> > &positionAsymm)
	{
		chargeAsymm.clear();
		positionAsymm.clear();
				

		// charges around atoms in asymmetric unit, positions in atom's local coordinate system
		chargeAsymm.resize(mCrystal.atoms.size());
		positionAsymm.resize(mCrystal.atoms.size());

		vector<double> chargesDipol, chargesQuadrupol;
		vector<Vector3d> positionsDipol, positionsQuadrupol;

        double distance = mSettings.multipoleExpansion.multipoleChargeDistance;// mMultipolePointChargeToAtomDistance;// 0.01;

		for (int atomIdx = 0; atomIdx < mCrystal.atoms.size(); atomIdx++)
		{
			chargeAsymm[atomIdx].push_back(mMultipoles[atomIdx].charge);
			positionAsymm[atomIdx].push_back(Vector3d());

            distributed_multipoles::dipole_as_point_charges(mMultipoles[atomIdx].dipole, distance, positionsDipol, chargesDipol);
            distributed_multipoles::quadrupole_as_point_charges(mMultipoles[atomIdx].quadrupole, distance, positionsQuadrupol, chargesQuadrupol);

			chargeAsymm[atomIdx].insert(chargeAsymm[atomIdx].end(), chargesDipol.begin(), chargesDipol.end());
			chargeAsymm[atomIdx].insert(chargeAsymm[atomIdx].end(), chargesQuadrupol.begin(), chargesQuadrupol.end());

			positionAsymm[atomIdx].insert(positionAsymm[atomIdx].end(), positionsDipol.begin(), positionsDipol.end());
			positionAsymm[atomIdx].insert(positionAsymm[atomIdx].end(), positionsQuadrupol.begin(), positionsQuadrupol.end());
		}


	}


	void StockholderAtomFormFactorCalcManager::calculatePointChargeMultipoleRepresentation(
		int clusterIdx,
		const std::vector< std::vector<double> >& chargeAsymm,
		const std::vector< std::vector<Vector3d> >& positionAsymm,// with respect to atom
		const std::vector<Vector3d> &asymmUnitAtomCartesian,
		std::vector<double>& charge,
		std::vector<Vector3d>& position)
	{
		charge.clear();
		position.clear();
				
		Matrix3d rotation;
		Vector3d translation;


		//int max_n_charges[] = { 1, 3, 10 };

		double weight;
		int idx = 0;
  
		for (auto& atom : mMultipoleBearingAtoms[clusterIdx])
		{
			int atomIdx = atom.first;
			crystal_structure_utilities::symmetryOperationCartesian(atom.second, mCrystal.unitCell, rotation, translation);
			//mWeightMultipolesWithOccupancy[clusterIdx][idx] ? weight = mCrystal.atoms[atomIdx].occupancy : weight = 1.0;

            weight = mCustomWeightMultipoles[clusterIdx][idx].value_or(mCrystal.atoms[atomIdx].occupancy);

            //if (mCustomWeightMultipoles[clusterIdx][idx].first)
            //    weight = mCustomWeightMultipoles[clusterIdx][idx].second;
            //else
            //    weight = mCrystal.atoms[atomIdx].occupancy;


            //if (mCustomWeightMultipoles[clusterIdx][idx].first)
            //    weight = mCustomWeightMultipoles[clusterIdx][idx].second;
            //else
            //    weight = mCrystal.atoms[atomIdx].occupancy;
			

            //packPointChargeRepresentationOfMultipoles(r0, q0, r, q);

            const vector<discamb::Vector3d> &r = positionAsymm[atomIdx];
            const vector<double> & q = chargeAsymm[atomIdx];


            for (int i = 0; i < q.size(); i++)
            {
            	charge.push_back(weight * q[i]);
            	Vector3d rRot = rotation * (r[i] + asymmUnitAtomCartesian[atomIdx]) + translation;

            	position.push_back(rRot);
            }

            //-------- end new

			idx++;
		}

	}

	void StockholderAtomFormFactorCalcManager::calculatePointChargeMultipoleRepresentation(
		std::vector< std::vector<double> > & charge,
		std::vector< std::vector<Vector3d> > & position)
	{

        int atomIdx, nAtoms = mCrystal.atoms.size();

        vector< vector<double> > chargeAsymm(nAtoms);
		vector< vector<Vector3d> > positionAsymm(nAtoms);

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            distributed_multipoles::multipoles_as_point_charges(mMultipoles[atomIdx], mSettings.multipoleExpansion.multipoleChargeDistance,
                positionAsymm[atomIdx], chargeAsymm[atomIdx]);
		

		int subsystemIdx, nSubsystems = mSettings.crystalFragments.size();
		
		charge.resize(nSubsystems);
		position.resize(nSubsystems);

        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {
            vector<double> weigths;
            vector<pair<string, string> > clusterAtoms;
            for (auto const& atom : mMultipoleBearingAtoms[subsystemIdx])
                clusterAtoms.push_back({ mCrystal.atoms[atom.first].label, atom.second.string() });

            for (int i = 0; i < mMultipoleBearingAtoms[subsystemIdx].size(); i++)
                weigths.push_back(mCustomWeightMultipoles[subsystemIdx][i].value_or(
                    mCrystal.atoms[mMultipoleBearingAtoms[subsystemIdx][i].first].occupancy));



            distributed_multipoles::point_charges_cluster(mCrystal, clusterAtoms, weigths, positionAsymm, chargeAsymm,
                position[subsystemIdx], charge[subsystemIdx] );


        }
		

	}



	bool StockholderAtomFormFactorCalcManager::multipoleCalculationConverged()
		const
	{
		if (mMultipoles.size() != mMultipolesPrevious.size())
			return false;

		int i, j, k, n = mMultipoles.size();
		
		double diff, maxDiff = 0;
		string what_diffs_max;
		char xyz[] = { 'x','y','z' };
		for (i = 0; i < n; i++)
		{
			diff = fabs(mMultipoles[i].charge - mMultipolesPrevious[i].charge);
			if (diff > maxDiff)
			{
				maxDiff = diff;
				what_diffs_max = string("q") + to_string(i + 1);
			}
			for (j = 0; j < 3; j++)
			{
				diff = fabs(mMultipoles[i].dipole[j] - mMultipolesPrevious[i].dipole[j]);
				if (diff > maxDiff)
				{
					maxDiff = diff;
					what_diffs_max = string("q") + to_string(i + 1) + xyz[j];
				}
				for (k = 0; k < 3; k++)
				{
					diff = fabs(mMultipoles[i].quadrupole(j, k) - mMultipolesPrevious[i].quadrupole(j, k));
					if (diff > maxDiff)
					{
						maxDiff = diff;
						what_diffs_max = string("q") + to_string(i + 1) + xyz[j];
						what_diffs_max += xyz[k];
					}

				}
			}

		}
		cout << "max diff " << what_diffs_max << " " << maxDiff << endl;
        return maxDiff <= mSettings.multipoleExpansion.multipoleConvergenceThreshold;
	}





	void StockholderAtomFormFactorCalcManager::setMultipolesTo0(
		std::vector<ElectricMultipoles>& multipoles)
	{
        for (auto& m : multipoles)
            m = ElectricMultipoles();
	}

	std::string StockholderAtomFormFactorCalcManager::numberTo000string(
		int n)
	{
		string s = to_string(n);
		if (s.size() == 1)
			s = string("00") + s;
		if (s.size() == 2)
			s = string("0") + s;
		return s;

	}


    void StockholderAtomFormFactorCalcManager::setSubsystemWfnCalcData(
        int subsystemIdx, 
        WaveFunctionCalculationData& data,
        const std::vector<double> pointCharges,
        const std::vector<Vector3d> pointChargePosition,
        const std::string& outputWfnFileName)
    {
        data = WaveFunctionCalculationData();
        
        data.wfnFileName = outputWfnFileName;
        data.qmSystem.pointChargeValue = pointCharges;
        data.qmSystem.pointChargePosition = pointChargePosition;

        if(!mSettings.fragmentWfnCalculation.empty())
            data.atomIdx2BasisSetMap = mSettings.fragmentWfnCalculation[subsystemIdx].atomicIdx2BasisSetMap;

        auto wfnCalculation = mSettings.wfnCalculation;
        
        if (!mSettings.fragmentWfnCalculation.empty())
            if (mSettings.fragmentWfnCalculation[subsystemIdx].wfnCalculation)
                wfnCalculation = mSettings.fragmentWfnCalculation[subsystemIdx].wfnCalculation.value();

        auto& subsystem = mSettings.crystalFragments[subsystemIdx];

        data.hardware = wfnCalculation.hardwareSettings.value_or(mSettings.hardware);


        data.jobName = subsystem.label;

        data.qmProgramSpecificData = wfnCalculation.qmProgramSpecificData;
        data.qmSettings = wfnCalculation.qmSettings;
        vector<ChemicalElement> elements;
        crystal_structure_utilities::convertToXyzAndElementList(mCrystal, subsystem.atoms.atomList, elements, data.qmSystem.positions);
        for (auto& element : elements)
            data.qmSystem.atomicNumbers.push_back(element.atomicNumber());
        data.qmSystem.charge = subsystem.charge;
        data.qmSystem.spin_multilicity = subsystem.spin_multiplicity;
    }

	void StockholderAtomFormFactorCalcManager::calculateWfns(
		int indexAddedToFileNames)
	{       
        int subsystemIdx, nSubsystems = mSettings.crystalFragments.size();
		string wfnFileName;
        string format = WaveFunctionDataGeneratorRunner::wfnFormatAsString(getFormat());
		vector< vector<double> > charge(nSubsystems);
		vector< vector<Vector3d> > position(nSubsystems);
        vector<DistributedMultipolesSet> dm(nSubsystems);
		
        if (!mMultipoles.empty() && !mMultipoleBearingAtoms.empty())
            calculatePointChargeMultipoleRepresentation(charge, position);



        // GAMESS -> wfn
        // ORCA -> molden -> wfx
        // NWChemm -> molden -> wfx
        // Q-Chem -> molden -> wfx
        // GAUSSIAN -> wfx
        // PSI4 -> fchk -> wfx ( ? or ? ) PSI4 -> molden -> wfx

        vector< WaveFunctionCalculationData> wfnCalcData(nSubsystems);
		for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
		{

			wfnFileName = mJobName + string("_") + mSettings.crystalFragments[subsystemIdx].label +
				          string("_wfn_") + numberTo000string(indexAddedToFileNames) + "." + format;

            setSubsystemWfnCalcData(subsystemIdx, wfnCalcData[subsystemIdx], charge[subsystemIdx],
                position[subsystemIdx], wfnFileName);
                
		}
        HardwareResources hardware = mSettings.wfnCalculation.hardwareSettings.value_or(mSettings.hardware);
        mWfnGeneratorRunner->runMultipleJobs(wfnCalcData, hardware.nCores, hardware.totalMemoryMB / hardware.nCores);

	}





	void StockholderAtomFormFactorCalcManager::setElectronDensityPartition(
		int subsystemIdx,
		const std::string& wfnFileName,
		const std::string& wfnFormat,
		std::shared_ptr<ElectronDensityPartition>& partition)
	{

        if (mSettings.electronDensityPartition.type == ElectronDensityPartitionType::MBIS)
        {
            partition = std::shared_ptr<ElectronDensityPartition>(new MbisPartition(wfnFileName, mSettings.electronDensityPartition.partitionSpecificData));
            return;
        }


        // Hirshfeld partition specific


        if (mSettings.electronDensityPartition.type == ElectronDensityPartitionType::Hirshfeld)
        {
            string sphericalAtomsFile;
            if(!mSettings.electronDensityPartition.partitionSpecificData.is_null())
                sphericalAtomsFile = mSettings.electronDensityPartition.partitionSpecificData.value("atoms file", string());
            if (!sphericalAtomsFile.empty())
                partition = shared_ptr<HirshfeldPartition>(new HirshfeldPartition(wfnFileName, sphericalAtomsFile));
            else
            {
                map<int, int> atomIdx2BasisSetMap, atomicNumber2BasisSetMap;
                vector<string> basis_sets;
                std::set<string> uniqueBasisSets;
                map<string, int> basisSet2Idx;
                auto wfnSettings = mSettings.wfnCalculation;
                if (!mSettings.fragmentWfnCalculation.empty())
                    if (mSettings.fragmentWfnCalculation[subsystemIdx].wfnCalculation.has_value())
                        wfnSettings = mSettings.fragmentWfnCalculation[subsystemIdx].wfnCalculation.value();
                for (auto item : wfnSettings.qmSettings.atomicNumber2BasisSetMap)
                    uniqueBasisSets.insert(item.second);

                if(!mSettings.fragmentWfnCalculation.empty())
                    for (auto item : mSettings.fragmentWfnCalculation[subsystemIdx].atomicIdx2BasisSetMap)
                        uniqueBasisSets.insert(item.second);


                basis_sets.insert(basis_sets.end(), uniqueBasisSets.begin(), uniqueBasisSets.end());
                for (int i = 0; i < basis_sets.size(); i++)
                    basisSet2Idx[basis_sets[i]] = i;

                for (auto item : wfnSettings.qmSettings.atomicNumber2BasisSetMap)
                    atomicNumber2BasisSetMap[item.first] = basisSet2Idx[item.second];


                if(!mSettings.fragmentWfnCalculation.empty())
                    for (auto item : mSettings.fragmentWfnCalculation[subsystemIdx].atomicIdx2BasisSetMap)
                        atomIdx2BasisSetMap[item.first] = basisSet2Idx[item.second];


                partition = shared_ptr<HirshfeldPartition>(new HirshfeldPartition(
                    wfnFileName, 
                    wfnSettings.qmSettings.qmMethod, 
                    wfnSettings.qmSettings.basisSet, 
                    basis_sets, 
                    atomicNumber2BasisSetMap, 
                    atomIdx2BasisSetMap, 
                    true,
                    mWfnGeneratorRunner->name(),
                    mWfnGeneratorRunner->getExecFolder(),
                    wfnSettings.qmSettings.relativisticMethod, 
                    mSettings.hardware.nCores, 
                    mSettings.hardware.totalMemoryMB));
            }

            static_pointer_cast<HirshfeldPartition>(partition)->setPower(mSettings.electronDensityPartition.power);
            partition->applySettings(mSettings.electronDensityPartition.partitionSpecificData);

            return;
        }

        on_error::throwException("unsupported electron density partition in specified for Hirshfeld Atom Model", __FILE__, __LINE__);
	}
	
	void StockholderAtomFormFactorCalcManager::electronDensityFromWfn(
		int indexAddedToFileNames,
		int subsystemIdx,
		int molecularOrbitalIdx)
	{
        WallClockTimer timer;
        timer.start();
        string format = WaveFunctionDataGeneratorRunner::wfnFormatAsString(getFormat());

		string wfnFileName = mJobName + string("_") + mSettings.crystalFragments[subsystemIdx].label +
			string("_wfn_") + numberTo000string(indexAddedToFileNames) + "." + format;

		
        vector<shared_ptr<ElectronDensityPartition> > edPartitions(1);
        setElectronDensityPartition(subsystemIdx, wfnFileName, format, edPartitions[0]);

		int nThreads = omp_get_max_threads();

		for (int i = 1; i < nThreads; i++)
			edPartitions.push_back(shared_ptr<ElectronDensityPartition>(edPartitions[0]->clone()));



		if (nThreads == 1)
		{
			for (auto atomInCluster : mSubsystemAtomsWhichAreRepresentatives[subsystemIdx])
			{
                if (mSettings.electronDensityPartition.edCalcAtomIncludeRange)
                {
                    edPartitions[0]->setIncludeRange(mSettings.electronDensityPartition.edCalcAtomIncludeRange.value());
                    edPartitions[0]->setAtom(atomInCluster);
                }

				int z = mSubsystemAtomicNumbers[subsystemIdx][atomInCluster];

				vector<double>& densities = mAtomInClusterElectronDensity[subsystemIdx][atomInCluster];
				vector<Vector3d>& points = mElementGridPoints[z];

				int pointIdx, nPoints = points.size();

				densities.resize(nPoints);
				Vector3d atomPosition = mSubsystemAtomicPositions[subsystemIdx][atomInCluster] * constants::Angstrom;
				for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
					densities[pointIdx] = edPartitions[0]->calculate(atomInCluster, atomPosition + points[pointIdx], molecularOrbitalIdx);
			}
		}
		else
		{

			vector<int> atomIndices(mSubsystemAtomsWhichAreRepresentatives[subsystemIdx].begin(), mSubsystemAtomsWhichAreRepresentatives[subsystemIdx].end());
			int nAtomsToBeConsidered = int(atomIndices.size());
#pragma omp parallel for
			for (int i = 0; i < nAtomsToBeConsidered; i++)
			{
				int threadNumber = omp_get_thread_num();
				int atomInCluster = atomIndices[i];
				int z = mSubsystemAtomicNumbers[subsystemIdx][atomInCluster];

				vector<double>& densities = mAtomInClusterElectronDensity[subsystemIdx][atomInCluster];
				vector<Vector3d>& points = mElementGridPoints[z];

				int pointIdx, nPoints = points.size();

				densities.resize(nPoints);
				Vector3d atomPosition = mSubsystemAtomicPositions[subsystemIdx][atomInCluster] * constants::Angstrom;

                if (mSettings.electronDensityPartition.edCalcAtomIncludeRange)
                {
                    edPartitions[threadNumber]->setIncludeRange(mSettings.electronDensityPartition.edCalcAtomIncludeRange.value());
                    edPartitions[threadNumber]->setAtom(atomInCluster);
                }


				for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
					densities[pointIdx] = edPartitions[threadNumber]->calculate(atomInCluster, atomPosition + points[pointIdx], molecularOrbitalIdx);
			}

		}

        cout << "electron density calculation, wfn " << subsystemIdx + 1 << " " << timer.stop() / 1000.0 << " s\n";
	}

    WaveFunctionFileFormat StockholderAtomFormFactorCalcManager::getFormat()
        const
    {
        vector<WaveFunctionFileFormat> supportedFormats;
        mWfnGeneratorRunner->supportedFormats(supportedFormats);


        if (find(supportedFormats.begin(), supportedFormats.end(), WaveFunctionFileFormat::wfx) != supportedFormats.end())
            return WaveFunctionFileFormat::wfx;
        if (find(supportedFormats.begin(), supportedFormats.end(), WaveFunctionFileFormat::wfn) != supportedFormats.end())
            return WaveFunctionFileFormat::wfn;
        if (find(supportedFormats.begin(), supportedFormats.end(), WaveFunctionFileFormat::molden) != supportedFormats.end())
            return WaveFunctionFileFormat::molden;
        if (find(supportedFormats.begin(), supportedFormats.end(), WaveFunctionFileFormat::fchk) != supportedFormats.end())
            return WaveFunctionFileFormat::fchk;

        return WaveFunctionFileFormat::mkl;
    }


	void StockholderAtomFormFactorCalcManager::electronDensityFromWfn(
		int indexAddedToFileNames,
		int subsystemIdx)
	{
        WallClockTimer timer;
        timer.start();

		// Hirshfeld density
        string format = WaveFunctionDataGeneratorRunner::wfnFormatAsString(getFormat());// mWfnGeneratorRunner->defaultFileFormat());

		string wfnFileName = mJobName + string("_") + mSettings.crystalFragments[subsystemIdx].label +
			string("_wfn_") + numberTo000string(indexAddedToFileNames) + "." + format;

		vector<shared_ptr<ElectronDensityPartition> > edPartitions(1);
        
        setElectronDensityPartition(subsystemIdx, wfnFileName, format, edPartitions[0]);
		//setHirshfeldPartitition(clusterIdx, wfnFileName, format, hPartitions[0]);
		
		int nThreads = omp_get_max_threads();

		for (int i = 1; i < nThreads; i++)
            edPartitions.push_back(shared_ptr<ElectronDensityPartition>(edPartitions[0]->clone()));
			//hPartitions.push_back(shared_ptr<HirshfeldPartition>(hPartitions[0]->clone()));

        if (mSettings.electronDensityPartition.edCalcAtomIncludeRange)
            for (int i = 0; i < nThreads; i++)
                edPartitions[i]->setIncludeRange(mSettings.electronDensityPartition.edCalcAtomIncludeRange.value());


		if (nThreads == 1)
		{

			for (auto atomInCluster : mSubsystemAtomsWhichAreRepresentatives[subsystemIdx])
			{
				int z = mSubsystemAtomicNumbers[subsystemIdx][atomInCluster];

                if (mSettings.electronDensityPartition.edCalcAtomIncludeRange)
                    edPartitions[0]->setAtom(atomInCluster);

				vector<double>& densities = mAtomInClusterElectronDensity[subsystemIdx][atomInCluster];
				vector<Vector3d>& points = mElementGridPoints[z];

				int pointIdx, nPoints = points.size();

				densities.resize(nPoints);
				Vector3d atomPosition = mSubsystemAtomicPositions[subsystemIdx][atomInCluster] * constants::Angstrom;
				for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
					densities[pointIdx] = edPartitions[0]->calculate(atomInCluster, atomPosition + points[pointIdx]);
			}
		}
		else
		{
			
			vector<int> atomIndices(mSubsystemAtomsWhichAreRepresentatives[subsystemIdx].begin(), mSubsystemAtomsWhichAreRepresentatives[subsystemIdx].end());
			int nAtomsToBeConsidered = int(atomIndices.size());
#pragma omp parallel for
			for (int i=0; i< nAtomsToBeConsidered; i++)
			{
				int threadNumber = omp_get_thread_num();
				int atomInCluster = atomIndices[i];
				int z = mSubsystemAtomicNumbers[subsystemIdx][atomInCluster];

				vector<double>& densities = mAtomInClusterElectronDensity[subsystemIdx][atomInCluster];
				vector<Vector3d>& points = mElementGridPoints[z];

				int pointIdx, nPoints = points.size();

				densities.resize(nPoints);
				Vector3d atomPosition = mSubsystemAtomicPositions[subsystemIdx][atomInCluster] * constants::Angstrom;

                if (mSettings.electronDensityPartition.edCalcAtomIncludeRange)
                    edPartitions[threadNumber]->setAtom(atomInCluster);



				for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
					densities[pointIdx] = edPartitions[threadNumber]->calculate(atomInCluster, atomPosition + points[pointIdx]);
			}

		}
		
        cout << "setting partition and electron density calculation, wfn " << subsystemIdx + 1 << " " << timer.stop() / 1000.0 << " s\n";
	
	}

	void StockholderAtomFormFactorCalcManager::setOrbital(
		int idx)
	{
		mOrbitalIdx=idx;
	}

	void StockholderAtomFormFactorCalcManager::electronDensityFromWfn(
		int indexAddedToFileNames)
	{

		int subsystemIdx, nSubsystems = mSettings.crystalFragments.size();
		string wfnFileName, densityFileName, representativesFileName;
		string format = WaveFunctionDataGeneratorRunner::wfnFormatAsString(mWfnGeneratorRunner->defaultFileFormat());

		for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
		{
			if(mOrbitalIdx<0)
				electronDensityFromWfn(indexAddedToFileNames, subsystemIdx);
			else
				electronDensityFromWfn(indexAddedToFileNames, subsystemIdx, mOrbitalIdx);
		}

	}


	//void StockholderAtomFormFactorCalcManager::setMultipolesfromTaam()
	//{
 //       MATTS_BankReader bankReader;
	//	vector<AtomType> atomTypes;
	//	vector<AtomTypeHC_Parameters> hcParameters;
	//	BankSettings bankSettings;
	//	string bankString;
	//	stringstream bankStream;
	//	default_ubdb_bank_string(bankString);
	//	bankStream << bankString;
	//	bankReader.read(bankStream, atomTypes, hcParameters, bankSettings, true);


	//	CrystalAtomTypeAssigner assigner;
	//	assigner.setAtomTypes(atomTypes);
	//	assigner.setDescriptorsSettings(DescriptorsSettings());
	//	vector < LocalCoordinateSystem<AtomInCrystalID> > lcs;
	//	vector<int> types;
	//	assigner.assign(mCrystal, types, lcs);

	//	ofstream out("TAAM_assignment");
	//	assigner.printAssignment(out, mCrystal, types, lcs);
	//	

	//	HC_ModelParameters multipoleModelPalameters;


	//	vector<int> atomicNumbers;
	//	vector<int> nonMultipolarAtoms;
	//	crystal_structure_utilities::atomicNumbers(mCrystal, atomicNumbers);
 //       taam_utilities::type_assignment_to_unscaled_HC_parameters(
	//		hcParameters, types, atomicNumbers,
	//		multipoleModelPalameters, true, nonMultipolarAtoms);

 //       vector<bool> nonMultipolar(atomicNumbers.size(), false);
 //       for (auto atomIdx : nonMultipolarAtoms)
 //           nonMultipolar[atomIdx] = true;

	//	vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcaCalculators;
	//	for (auto coordinateSystem : lcs)
	//		lcaCalculators.push_back(
	//			shared_ptr<LocalCoordinateSystemInCrystal>(
	//				new LocalCoordinateSystemCalculator(coordinateSystem, mCrystal)));

	//	
	//	vector<double> coreElectronCharges;
	//	for (auto& wfn : multipoleModelPalameters.wfn_parameters)
	//	{
	//		double q_core_el = 0;
	//		for (auto idx : wfn.core_orbitals_indices)
	//			q_core_el -= wfn.orbital_occupancy[idx];
	//		coreElectronCharges.push_back(q_core_el);
	//	}

	//	int atomIdx, nAtoms = mCrystal.atoms.size();
	//	int typeIdx, wfnTypeIdx;
	//	Matrix3d localCoordinatesMatrix;
 //       mMultipoles.assign(nAtoms, ElectricMultipoles());

	//	out << "\n\n Multipoles from TAAM (charge, dipole):\n\n";

	//	for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
	//	{
	//		out << setw(8) << mCrystal.atoms[atomIdx].label << " ";

 //           if (nonMultipolar[atomIdx])
 //           {
 //               mMultipoles[atomIdx].charge = 0;
 //               mMultipoles[atomIdx].dipole = Vector3d();
 //           }
 //           else
 //           {
 //               typeIdx = multipoleModelPalameters.atom_to_type_map[atomIdx];
 //               wfnTypeIdx = multipoleModelPalameters.atom_to_wfn_map[atomIdx];

 //               mMultipoles[atomIdx].charge = atomicNumbers[atomIdx] -
 //                   multipoleModelPalameters.type_parameters[typeIdx].p_val +
 //                   coreElectronCharges[wfnTypeIdx];
 //               lcaCalculators[atomIdx]->calculate(localCoordinatesMatrix, mCrystal);
 //               if (multipoleModelPalameters.type_parameters[typeIdx].p_lm.size() > 1)
 //               {
 //                   vector<double>& p = multipoleModelPalameters.type_parameters[typeIdx].p_lm[1];
 //                   mMultipoles[atomIdx].dipole =
 //                       dipoleMoment(
 //                           { p[2],p[0],p[1] },
 //                           localCoordinatesMatrix,
 //                           multipoleModelPalameters.wfn_parameters[wfnTypeIdx].deformation_valence_power[1],
 //                           multipoleModelPalameters.type_parameters[typeIdx].kappa_deformation_valence,
 //                           multipoleModelPalameters.wfn_parameters[wfnTypeIdx].deformation_valence_exponent);
 //               }
 //               else
 //                   mMultipoles[atomIdx].dipole.set(0.0, 0.0, 0.0);
 //           }

	//		out << setw(12) << setprecision(4) << fixed << mMultipoles[atomIdx].charge
	//			<< setw(12) << mMultipoles[atomIdx].dipole[0]
	//			<< setw(12) << mMultipoles[atomIdx].dipole[1]
	//			<< setw(12) << mMultipoles[atomIdx].dipole[2] << "\n";
	//	}


	//	out.close();
	//}

//	bool StockholderAtomFormFactorCalcManager::tryToFindAndReadMultipolesFile(
//		int &previousStep, string &fName)

    bool StockholderAtomFormFactorCalcManager::tryToFindAndReadMultipolesFile(
        int& previousStep,
        std::string& fName,
        std::vector< ElectricMultipoles>& multipoles)
        const
    {
        multipoles.clear();
        string fileCore = mJobName + string("_multipoles_");
        vector<pair<filesystem::file_time_type, string> > multipole_files;
        int index;
        for (auto& entry : filesystem::directory_iterator(filesystem::current_path()))
            if (filesystem::is_regular_file(entry.path()))
            {
                string fileName = entry.path().filename().string();
                if (fileName.find(fileCore) == 0)
                {
                    multipole_files.push_back({ filesystem::last_write_time(fileName), fileName });
                }
            }

        if (multipole_files.empty())
            return false;

        sort(multipole_files.begin(), multipole_files.end());

        readMultipolesFromFile(multipole_files.back().second, multipoles);
        fName = multipole_files.back().second;
        // get step number

        string numberStr = fName.substr(mJobName.size() + 12);
        index = 0;
        if (numberStr.size() == 3)
        {
            if (numberStr[0] == '0')
            {
                if (numberStr[1] == '0')
                    index = stoi(numberStr.substr(2));
                else
                    index = stoi(numberStr.substr(1));
            }
            else
                index = stoi(numberStr);

        }

        previousStep = index;

        return true;

    }


    bool StockholderAtomFormFactorCalcManager::tryToSetMultipolesFromFile(
        int& previousStep, 
        std::string& fName)
	{
        return tryToFindAndReadMultipolesFile(
                previousStep,
                fName,
                mMultipoles);
		//string fileCore = mJobName + string("_multipoles_");
		////vector<pair<int, string> > multipole_files;
  //      //filesystem::file_time_type
  //      vector<pair<fs::file_time_type, string> > multipole_files;
		//int index;
		//for(auto &entry: fs::directory_iterator(fs::current_path()))
		//	if (fs::is_regular_file(entry.path()))
		//	{
		//		string fileName = entry.path().filename().string();
		//		if (fileName.find(fileCore) == 0)
		//		{
		//			//string numberStr = fileName.substr(mJobName.size() + 12);
		//			//if (numberStr.size() == 3)
		//			//{
		//			//	if (numberStr[0] == '0')
		//			//	{
		//			//		if (numberStr[1] == '0')
		//			//			index = stoi(numberStr.substr(2));
		//			//		else
		//			//			index = stoi(numberStr.substr(1));
		//			//	}
		//			//	else
		//			//		index = stoi(numberStr);

		//			//	multipole_files.push_back({ index, fileName });
		//			//}
  //                  multipole_files.push_back({fs::last_write_time(fileName), fileName });
		//		}
		//	}
		// 
		//if (multipole_files.empty())
		//	return false;
  //      
		//sort(multipole_files.begin(), multipole_files.end());

		//readMultipolesFromFile(multipole_files.back().second, mMultipoles);
		//fName = multipole_files.back().second;
  //      // get step number

  //      string numberStr = fName.substr(mJobName.size() + 12);
  //      index = 0;
  //      if (numberStr.size() == 3)
  //      {
  //      	if (numberStr[0] == '0')
  //      	{
  //      		if (numberStr[1] == '0')
  //      			index = stoi(numberStr.substr(2));
  //      		else
  //      			index = stoi(numberStr.substr(1));
  //      	}
  //      	else
  //      		index = stoi(numberStr);

  //      }
  //              
		//previousStep = index;

		//return true;
	}

	void StockholderAtomFormFactorCalcManager::calculateElectronDensityWithPredefinedMultipoles(
		bool useTaamMultipolesIfNoOtherProvided,
		bool lookForMultipolesFile,
		const std::vector<ElectricMultipoles>& multipoles)
	{
		//cout << "calculate electron density with predefined multipoles" << endl;
		bool hasMultipolesFromFile = false;
        //double d;
        //cin >> d;
        //cout << d * d << endl;
		int previousStep = 0;
		if (multipoles.empty())
		{
			if (lookForMultipolesFile)
			{
				string multipolesFile;
				//hasMultipolesFromFile = tryToFindAndReadMultipolesFile(previousStep, multipolesFile);
                hasMultipolesFromFile = tryToSetMultipolesFromFile(previousStep, multipolesFile);
				if(hasMultipolesFromFile)
					cout << "read multipoles from file " << multipolesFile << endl;
			}

			if (useTaamMultipolesIfNoOtherProvided && !hasMultipolesFromFile)
			{
				//setMultipolesfromTaam();
			}
		}
		else
			mMultipoles = multipoles;

		int step = previousStep + 1; 
		string indexString = numberTo000string(step);
		calculateWfns(step);
		electronDensityFromWfn(step);
		calculateMultipoles();
		saveMultipolesToFile(mJobName + string("_multipoles_") + indexString, mMultipoles);


	}


	void StockholderAtomFormFactorCalcManager::calculateElectronDensityInSelfConsistentMultipoleEmbedding(
		const vector<ElectricMultipoles> &multipoles,
        bool multipolesAreConverged,
		int previousStep)
	{

		mMultipoles = multipoles;

		int step = previousStep;
		bool noIterationsNeeded = mMultipoleBearingAtoms.empty() || multipolesAreConverged;
		int iteration = 0;

		do
		{
			step++;
			string indexString = numberTo000string(step);
			calculateWfns(step);
			electronDensityFromWfn(step);
			mMultipolesPrevious = mMultipoles;
			calculateMultipoles();
			
			saveMultipolesToFile(mJobName + string("_multipoles_") + indexString, mMultipoles);

			if (noIterationsNeeded)
				break;
			
			iteration++;
			if (iteration >= mSettings.multipoleExpansion.maxN_Steps)// mMaxN_EE_SCF_Cycles)
				break;
		} 
		while (!multipoleCalculationConverged());

		//readElectronDensityIntegrationGrids(step);


	}

	void StockholderAtomFormFactorCalcManager::allocateElectronDensityContainers()
	{
		//int nAtoms = mRepresentatives.size();

		//mIntegrationGridPoints.resize(nAtoms);
		//mAtomicDensities.resize(nAtoms);
		//mIntegrationWeights.resize(nAtoms);


		//for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
		//{
		//	nAtomRepresentatives = mRepresentatives[atomIdx].size();

		//	mIntegrationGridPoints[atomIdx].resize(nAtomRepresentatives);
		//	mAtomicDensities[atomIdx].resize(nAtomRepresentatives);
		//	mIntegrationWeights[atomIdx].resize(nAtomRepresentatives);
		//}

		// new things
        int subsystemIdx, nSubsystems = mSettings.crystalFragments.size();
		mAtomInClusterElectronDensity.resize(nSubsystems);
        for (subsystemIdx = 0; subsystemIdx < nSubsystems; subsystemIdx++)
        {
            int nAtomsInSubsystem = mSettings.crystalFragments[subsystemIdx].atoms.atomList.size() +
                mSettings.crystalFragments[subsystemIdx].atoms.cappingHydrogens.size();
            mAtomInClusterElectronDensity[subsystemIdx].resize(nAtomsInSubsystem);
        }
		
		setElementalIntegrationGrids();
	}

	// [z][point idx]
//	std::vector < std::vector<Vector3d> > mElementGridPoints;
	// [z][point idx]
//	std::vector < std::vector<double> > mElementIntegrationWeights;
	void StockholderAtomFormFactorCalcManager::setElementalIntegrationGrids()
	{
        
        //const pair<int, int> defaultGrid({ 590, 75 });
        const int maxZ = 113;
        //vector<string> gridSpecs, words;
        //vector<pair<int, int> > gridInfo(113, defaultGrid);
        vector<ham_settings::IntegrationGrid> gridInfo(maxZ, mSettings.formFactorsCalculation.integrationGrid);
        for (auto& element_specif_grid : mSettings.formFactorsCalculation.elementSpecificIntegrationGrid)
            gridInfo[element_specif_grid.first] = element_specif_grid.second;


        //string_utilities::split(mGrid, gridSpecs);

        //if (!gridSpecs.empty())
        //{
        //    string_utilities::split(gridSpecs[0] , words, ',');
        //    if (words.size() == 2)
        //        gridInfo.assign(113, { stoi(words[0]), stoi(words[1]) });

        //    for (auto item : gridSpecs)
        //    {
        //        string_utilities::split(item, words, ',');
        //        
        //        if (words.size() != 3 && words.size() != 2)
        //            on_error::throwException(string("invalid specification of integration grid"), __FILE__, __LINE__);
        //        
        //        if (words.size() == 3)
        //        {
        //            vector<int> atomicNumbers;
        //            basic_chemistry_utilities::getElementsList(words[0], atomicNumbers);
        //            int nAng = stoi(words[1]);
        //            int nRad = stoi(words[2]);
        //            for(int z: atomicNumbers)
        //                gridInfo[z] = { nAng, nRad };
        //        }
        //    }
        //}

        mElementGridPoints.clear();
        mElementIntegrationWeights.clear();


        mElementGridPoints.resize(113);
        mElementIntegrationWeights.resize(113);

        std::set<int> uniqueZ;

        for (auto& subsetZ : mSubsystemAtomicNumbers)
            uniqueZ.insert(subsetZ.begin(), subsetZ.end());

        vector<Vector3d> angularGridPoints, gridPoints;
        vector<double> radialPoints, radialWeights, angularWeights, weights;
        Vector3d r, atomPosition;


        for (auto z : uniqueZ)
        {
            int nAng = gridInfo[z].angularGridSize;
            int nRad = gridInfo[z].radialGridSize;
            lebedev_laikov::get_grid(nAng , angularGridPoints, angularWeights);

            gridPoints.resize(nAng * nRad);
            weights.resize(nAng * nRad);


            //radial_grid::treutler_ahlrichs(z, nRad, radialPoints, radialWeights);
            radial_grid::mura_knowles(z, nRad, radialPoints, radialWeights);
            int pointIdx = 0;
            for (int angularIdx = 0; angularIdx < nAng; angularIdx++)
                for (int radialIdx = 0; radialIdx < nRad; radialIdx++)
                {
                    weights[pointIdx] = radialWeights[radialIdx] * angularWeights[angularIdx] *
                        radialPoints[radialIdx] * radialPoints[radialIdx] * M_PI * 4.0;
                    gridPoints[pointIdx] = radialPoints[radialIdx] * angularGridPoints[angularIdx];
                    pointIdx++;
                }
            mElementGridPoints[z] = gridPoints;
            mElementIntegrationWeights[z] = weights;
        }


        return;

	}



	void StockholderAtomFormFactorCalcManager::calculateAtomicDensities()
	{
        ham_settings::WfnCalculation::Restart restartOption = mSettings.wfnCalculation.restart; 
        int restartStep = mSettings.wfnCalculation.restartStep;
        string previousMultipolesFile = mSettings.multipoleExpansion.previousMultipolesFile;


        // set multipoles to 0
        for (auto& m : mMultipoles)
            m = ElectricMultipoles();
        for (auto& m : mMultipolesPrevious)
            m = ElectricMultipoles();
        //
        bool noMultipoeBearingAtoms = true;
        for (auto& atoms : mMultipoleBearingAtoms)
            if (!atoms.empty())
                noMultipoeBearingAtoms = false;
		//mElectrostaticEmbeddingSCF_CycleIndex
		// option (1) start electrostatic embedding scf from scratch
		switch (restartOption)
		{
			case  ham_settings::WfnCalculation::Restart::from_scratch:
				if (mSettings.multipoleExpansion.maxN_Steps == 0 || noMultipoeBearingAtoms)
					calculateElectronDensityWithPredefinedMultipoles(mSettings.multipoleExpansion.startWithTaam, mTryToReadMultipolesFromFile);
                else
                {
                    vector<ElectricMultipoles>  multipoles;
                    int previousStep = 0;
                    string fName; 

                    tryToFindAndReadMultipolesFile(previousStep, fName, multipoles);
                    calculateElectronDensityInSelfConsistentMultipoleEmbedding(multipoles, false, previousStep);
                }
				return;
			case  ham_settings::WfnCalculation::Restart::from_converged_wfn:
				electronDensityFromWfn(restartStep);
				return;
		}


		

		if (restartOption == ham_settings::WfnCalculation::Restart::from_non_converged_wfn)
		{
			electronDensityFromWfn(restartStep);
			calculateMultipoles();
		}
		
		if (!previousMultipolesFile.empty())
		{
			readMultipolesFromFile(previousMultipolesFile, mMultipolesPrevious);
			if (multipoleCalculationConverged())
				return;
		}

        bool multipolesConverged = false;

		if (restartOption == ham_settings::WfnCalculation::Restart::from_not_converged_multipoles)
			readMultipolesFromFile(mJobName + string("_multipoles_") + numberTo000string(restartStep), mMultipoles);

        if (restartOption == ham_settings::WfnCalculation::Restart::from_converged_multipoles)
        {
            readMultipolesFromFile(mJobName + string("_multipoles_") + numberTo000string(restartStep), mMultipoles);
            multipolesConverged = true;
        }

		calculateElectronDensityInSelfConsistentMultipoleEmbedding(mMultipoles, multipolesConverged, restartStep);
		

	}


	void StockholderAtomFormFactorCalcManager::saveMultipolesToFile(
		const std::string fileName,
		const std::vector<ElectricMultipoles>& multipoles)
	{
		ofstream out(fileName);
		if (!out.good())
			on_error::throwException(string("cannot open file '") + fileName +
				string("' for saving multipole moments for electrostatic embedding SCF calculations"),
				__FILE__, __LINE__);

		out << "n_centers " << multipoles.size() << "\n";
		for (auto const& m : multipoles)
			out << fixed << setprecision(8) << setw(15) << m.charge
				<< setw(15) << m.dipole[0]
				<< setw(15) << m.dipole[1]
				<< setw(15) << m.dipole[2]
				<< setw(15) << m.quadrupole(0, 0)
				<< setw(15) << m.quadrupole(1, 1)
				<< setw(15) << m.quadrupole(2, 2)
				<< setw(15) << m.quadrupole(0, 1)
				<< setw(15) << m.quadrupole(0, 2)
				<< setw(15) << m.quadrupole(1, 2) << "\n";
		
		out.close();
	}

	void StockholderAtomFormFactorCalcManager::readMultipolesFromFile(
		const std::string fileName,
		std::vector<ElectricMultipoles>& multipoles)
	{
		multipoles.clear();
		ifstream in(fileName);
		if (!in.good())
			on_error::throwException(string("cannot open file '") + fileName +
				string("' for reading multipole moments for electrostatic embedding SCF calculations"),
				__FILE__, __LINE__);
		
		string s;
		int i, n;
		in >> s >> n;
		multipoles.resize(n);
		for (i = 0; i < n; i++)
		{
			in >> multipoles[i].charge >> multipoles[i].dipole[0] >> multipoles[i].dipole[1] >> multipoles[i].dipole[2]
				>> multipoles[i].quadrupole(0, 0) >> multipoles[i].quadrupole(1, 1) >> multipoles[i].quadrupole(2, 2)
				>> multipoles[i].quadrupole(0, 1) >> multipoles[i].quadrupole(0, 2) >> multipoles[i].quadrupole(1, 2);
			multipoles[i].quadrupole(1, 0) = multipoles[i].quadrupole(0, 1);
			multipoles[i].quadrupole(2, 0) = multipoles[i].quadrupole(0, 2);
			multipoles[i].quadrupole(2, 1) = multipoles[i].quadrupole(1, 2);
		}
		in.close();
	}


}
