#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/BasicUtilities/OnError.h"

namespace discamb {
    
    int Crystal::atomIdx(
        const std::string& atomLabel)
        const
    {
        int idx=0;
        if (!atomIdx(atomLabel, idx))
            on_error::throwException("no atom with label '" + atomLabel + "found in the crystal structure", __FILE__, __LINE__);
        return idx;
    }

    bool Crystal::atomIdx(
        const std::string& atomLabel,
        int& idx)
        const
    {
        idx = 0;
        for (int i = 0; i < this->atoms.size(); i++)
            if (atoms[i].label == atomLabel)
            {
                idx = i;
                return true;
            }
        return false;
    }


}

