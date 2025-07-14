#ifndef _DISCAMB_HC_MODEL_SLATERORBITALWFNDATA_H_
#define _DISCAMB_HC_MODEL_SLATERORBITALWFNDATA_H_
#include "discamb/HC_Model/HC_WfnData.h"

#include <utility>
#include <vector>
#include <string>
#include <map>

struct LocalDef_SlaterTerm
{
	double coefficient;
	int powerR;
	double exponent;
};

struct LocalDef_Orbital {
	int principal_number;
	int azimutal_number;
	int nTerms;
	LocalDef_SlaterTerm term[20];
};

struct LocalDef_Wfn {
	std::string label;
	int atomicNumber;
	int charge;
	int nOrbs;
	LocalDef_Orbital orbs[20];
	
};

namespace discamb {

    /** \ingroup HC_Model 

        \brief Selected atomic wavefunctions data parent class.

        Parent class providing atomic wavefunction data (discamb::HC_WfnBankEntry).
        The data include also info on assignment of the orbital to valence and core 
        (for detail on the alassignment algorithm see \ref algorithms_for_parameterization "algorithms section" )
    */

    class SlaterOrbitalWfnData
    {
    public:
        SlaterOrbitalWfnData(LocalDef_Wfn[], std::string[], int);
        virtual ~SlaterOrbitalWfnData();
        /** \brief Get all wavefunction data entries.*/
        void getEntries(std::vector<discamb::HC_WfnBankEntry> &entries);
        /** \brief Get wavefunction data entry.
        
        The entry is specified with label \p atomType, \p hasType is set to true or false depending if 
        the entry is present, if it is not then empty entry is returned.*/
        const discamb::HC_WfnBankEntry &getEntry(const std::string &atomType, bool &hasType) const;

        /** \brief Get wavefunction data entry. 
        
            The specified with label \p atomType, 
            if it is not available then discamb::Exception is thrown.*/
        const discamb::HC_WfnBankEntry &getEntry(const std::string &atomType) const;

    private:
        std::map<std::string, discamb::HC_WfnBankEntry> mEntries;
        discamb::HC_WfnBankEntry mEmptyEntry;
    };

}

#endif /*_DISCAMB_HC_MODEL_SLATERORBITALWFNDATA_H_*/