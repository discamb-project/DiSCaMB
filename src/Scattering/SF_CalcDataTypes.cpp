#include "discamb/Scattering/SF_CalcDataTypes.h"

#include <map>
#include <algorithm>

using namespace std;

namespace discamb {

void ScatteringParameters::groupAtomTypes(
    const std::vector<std::string> &atomType,
    std::vector<std::string> &types, 
    std::vector<int> &atomToTypeMap)
{
    int i, nAtoms;
    map<string, int> atomTypeIndices;
    types.clear();
    

    nAtoms = atomType.size();
    atomToTypeMap.resize(nAtoms);

    for (i = 0; i<nAtoms; i++)
    {

        if (find( types.begin(), types.end(), atomType[i]) == types.end())
        {
            atomTypeIndices[atomType[i]] = types.size();
            types.push_back(atomType[i]);
        }

        atomToTypeMap[i] = atomTypeIndices[atomType[i]];
    }
}

} // namespace discamb
