#include "discamb/StructuralProperties/CovalentRadiousBondDetector.h"
#include "discamb/BasicChemistry/chemical_element_data.h"

#include "discamb/BasicUtilities/constants.h"

#include <cmath>

namespace discamb {

    CovalentRadiousBondDetector::CovalentRadiousBondDetector()
    {
        setThreshold(0.1);
    }

    CovalentRadiousBondDetector::~CovalentRadiousBondDetector()
    {
    }


    void CovalentRadiousBondDetector::setThreshold(
        double threshold,
        bool inAngstroms)
    {
        mThreshold = threshold;
        if (!inAngstroms)
            mThreshold /= constants::Angstrom;
        //constants::Angstrom
        //inAngstroms ? mThreshold = structural_utilities::angstromsToBohrs(threshold) : mThreshold = threshold;
    }


    double CovalentRadiousBondDetector::getThreshold(
        bool inAngstroms)
        const
    {
        if (inAngstroms)
            return mThreshold;// structural_utilities::bohrsToAngstroms(mThreshold);
        else return mThreshold * constants::Angstrom;
    }

    void CovalentRadiousBondDetector::set(
        const std::string &thresholdInAngstroms)
    {
        setThreshold(stod(thresholdInAngstroms));
    }

    bool CovalentRadiousBondDetector::areBonded(
        int atomic_number_1,
        const Vector3d atom1_position,
        int atomic_number_2,
        const Vector3d atom2_position)
        const
    {
        Vector3d diff = atom1_position - atom2_position;
        double interatomicDistance = sqrt(diff*diff);
        if (interatomicDistance > chemical_element_data::covalentRadius(atomic_number_1) + chemical_element_data::covalentRadius(atomic_number_2) + mThreshold)
            return false;
        return true;
    }
    
    // distance in Angstroms
    bool CovalentRadiousBondDetector::areBonded(
        int atomic_number_1,
        int atomic_number_2,
        double distance)
        const
    {
        return (distance < chemical_element_data::covalentRadius(atomic_number_1) + chemical_element_data::covalentRadius(atomic_number_2) + mThreshold);
    }

}// namespace discamb

