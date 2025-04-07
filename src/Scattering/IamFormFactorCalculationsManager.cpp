#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/Scattering/NGaussianFormFactorsTable.h"
#include "discamb/BasicUtilities/on_error.h"

#include <map>

using namespace std;

namespace discamb {

     
    
    IamFormFactorCalculationsManager::IamFormFactorCalculationsManager()
    {
    }

    IamFormFactorCalculationsManager::~IamFormFactorCalculationsManager()
    {
    }

    //void IamFormaFactorCalculationsManager::setStructure(
    //    const Crystal &crystal)


    IamFormFactorCalculationsManager::IamFormFactorCalculationsManager(
        const Crystal &crystal,
        const std::string &table_type)
    {
        mReciprocalSpaceUnitCell.set(crystal.unitCell);
        int atomIdx, nAtoms, typeIdx;
        map<string, int> atomTypeToIndexMap;
        map<string, int>::const_iterator it;
        string type;

        nAtoms = crystal.atoms.size();
        mAtomToFormFactorMap.resize(nAtoms);

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            type = crystal.atoms[atomIdx].type;
            if (type.empty())
            {
                int atomicNumber=0;
                if (!crystal.atoms[atomIdx].label.empty())
                    atomicNumber = basic_chemistry_utilities::atomicNumberFromLabel(crystal.atoms[atomIdx].label);
                if (atomicNumber==0)
                    on_error::throwException("unknown atom type when initializing independent atom model calculations engine", __FILE__, __LINE__);
                type = periodic_table::symbol(atomicNumber);
            }

            it = atomTypeToIndexMap.find(type);
            if (it == atomTypeToIndexMap.end())
            {
                typeIdx = mFormFactors.size();
                
                if (!n_gaussian_form_factors_table::hasFormFactor(type, table_type))
                    on_error::throwException(string("undefined scatterer type: '") + type + string("'"), __FILE__, __LINE__);

                mFormFactors.push_back(n_gaussian_form_factors_table::getFormFactor(type, table_type));
                atomTypeToIndexMap[type] = typeIdx;
            }
            

            mAtomToFormFactorMap[atomIdx] = atomTypeToIndexMap[type];

        }
        mFormFactorValues.resize(mFormFactors.size());
    }

    IamFormFactorCalculationsManager::IamFormFactorCalculationsManager(
        const Crystal &crystal,
        const std::map<std::string, NGaussianFormFactor> &formFactors)
    {
        mReciprocalSpaceUnitCell.set(crystal.unitCell);
        int atomIdx, nAtoms, typeIdx;
        map<string, int> atomTypeToIndexMap;
        map<string, int>::const_iterator it;
        string type;

        nAtoms = crystal.atoms.size();
        mAtomToFormFactorMap.resize(nAtoms);

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            type = crystal.atoms[atomIdx].type;
            if (type.empty())
            {
                int atomicNumber = 0;
                if (!crystal.atoms[atomIdx].label.empty())
                    atomicNumber = basic_chemistry_utilities::atomicNumberFromLabel(crystal.atoms[atomIdx].label);
                if (atomicNumber == 0)
                    on_error::throwException("unknown atom type when initializing independent atom model calculations engine", __FILE__, __LINE__);
                type = periodic_table::symbol(atomicNumber);
            }

            it = atomTypeToIndexMap.find(type);
            if (it == atomTypeToIndexMap.end())
            {
                typeIdx = mFormFactors.size();
                if(formFactors.find(type) == formFactors.end())
                    on_error::throwException(string("undefined scatterer type: '") + type + string("'"), __FILE__, __LINE__);

                mFormFactors.push_back(formFactors.find(type)->second);
                atomTypeToIndexMap[type] = typeIdx;
            }


            mAtomToFormFactorMap[atomIdx] = atomTypeToIndexMap[type];

        }
        mFormFactorValues.resize(mFormFactors.size());
    }

    void IamFormFactorCalculationsManager::update(
        const std::vector<AtomInCrystal> &atoms)
    {
    }

    std::complex<double> IamFormFactorCalculationsManager::calculateFrac(int atomIdx, const Vector3i &hkl)
        const
    {
        return calculateCart(atomIdx, convertMillerIndexToCartesian(mReciprocalSpaceUnitCell, hkl));
    }

    std::complex<double> IamFormFactorCalculationsManager::calculateCart(
        int atomIdx, 
        const Vector3d &hkl)
        const
    {
        return mFormFactors[mAtomToFormFactorMap[atomIdx]].calculate_h(sqrt(hkl*hkl));
    }
    

    void IamFormFactorCalculationsManager::calculateFrac(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom) 
        const
    {
        Vector3d hklCart;
        mReciprocalSpaceUnitCell.fractionalToCartesian(hkl, hklCart);
        calculateCart(hklCart, formFactors, includeAtom);
    }

    void IamFormFactorCalculationsManager::calculateCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int nTypes = mFormFactors.size();
        int nAtoms = mAtomToFormFactorMap.size();
        formFactors.resize(nAtoms);
        double h = sqrt(hkl * hkl);
        for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
            mFormFactorValues[typeIdx] = mFormFactors[typeIdx].calculate_h(h);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (includeAtom[atomIdx])
                formFactors[atomIdx] = mFormFactorValues[mAtomToFormFactorMap[atomIdx]];//calculateCart(atomIdx, q);
    }


    void IamFormFactorCalculationsManager::calculateFrac(
        const std::vector<Vector3d>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom) const
    {
        int hklIdx, nHkl, atomIdx, nAtoms;
        nHkl = hkl.size();
        nAtoms = mAtomToFormFactorMap.size();
        formFactors.resize(nHkl);
        int nTypes = mFormFactors.size();
        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {

            formFactors[hklIdx].resize(nAtoms);
            calculateFrac(hkl[hklIdx], formFactors[hklIdx], includeAtom);
            //Vector3d q;
            //mReciprocalSpaceUnitCell.fractionalToCartesian(hkl[hklIdx], q);
            //double h = sqrt(q*q);
            //for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
            //    mFormFactorValues[typeIdx] = mFormFactors[typeIdx].calculate_h(h);

            //for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            //    if (includeAtom[atomIdx])
            //        formFactors[hklIdx][atomIdx] = mFormFactorValues[mAtomToFormFactorMap[atomIdx]];//calculateCart(atomIdx, q);
                
        }
    }

}
