#ifndef _DISCAMB_HC_MODEL_HC_MODELPARAMETERS_H_
#define _DISCAMB_HC_MODEL_HC_MODELPARAMETERS_H_

#include "HC_WfnData.h"
#include "HC_AtomTypeParameters.h"


namespace discamb {

/** \ingroup HC_Model 

\brief Container for data parameterizing Hansen-Coppens multipole model of electron density 
       not including local coordinate systems.




*/



struct HC_ModelParameters
{
    std::vector<HC_WfnType> wfn_parameters;
    std::vector<HC_AtomTypeParameters> type_parameters;
    std::vector<int> atom_to_wfn_map;
    std::vector<int> atom_to_type_map;
};

} //namespace discamb

#endif /*_DISCAMB_HC_MODEL_HC_MODELPARAMETERS_H_*/

