#pragma once

#ifndef _DISCAMB_STRUCTURALPROPERTIES_COVALENTRADIOUSBONDDETECTOR_H_
#define _DISCAMB_STRUCTURALPROPERTIES_COVALENTRADIOUSBONDDETECTOR_H_

#include "discamb/MathUtilities/Vector3.h"

#include <string>

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    /**

    detects bond between atoms (A and B) if the interatomic R distance fulfills the following condition:

    R <= covalent radious (atom A) + covalent radious (atom B) + threshold

    */

    class CovalentRadiousBondDetector {
    public:

        CovalentRadiousBondDetector();
        ~CovalentRadiousBondDetector();

        // threshold in Angstroms
        void setThreshold(double threshold, bool inAngstroms = true);
        double getThreshold(bool inAngstroms = true) const;
        void set(const std::string &thresholdInAngstroms);

        // positions in Angstroms
        bool areBonded(int atomic_number_1, const Vector3d atom1_position, int atomic_number_2, const Vector3d atom2_position)	const;
        // distance in Angstroms
        bool areBonded(int atomic_number_1, int atomic_number_2, double distance) const;
    private:
        double mThreshold; // in Angstroms
    };
    /** @}*/
} // namespace discamb


#endif
