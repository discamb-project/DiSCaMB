#ifndef _DISCAMB_HC_MODEL_SUCOPPENSMACCHIDATA_H_
#define _DISCAMB_HC_MODEL_SUCOPPENSMACCHIDATA_H_


#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/HC_Model/SlaterOrbitalWfnData.h"

#include <utility>
#include <vector>
#include <string>
#include <map>


namespace discamb {

    /** \ingroup HC_Model 

        \brief Selected atomic wavefunctions data from Su & Coppens (1998) and Macchi & Coppens (2001).

        Class providing atomic wavefunction data (discamb::HC_WfnBankEntry) based on 
        Su & Coppens \cite Su_Coppens_1998 and Macchi & Coppens publications \cite Macchi_Coppens_2001 .
        The data include also info on assignment of the orbital to valence and core 
        (for detail on the alassignment algorithm see \ref algorithms_for_parameterization "algorithms section" )
    */

    class SuCoppensMacchiData : public SlaterOrbitalWfnData
    {
    public:
        SuCoppensMacchiData();
    };

}

#endif /*_DISCAMB_HC_MODEL_SUCOPPENSMACCHIDATA_H_*/


