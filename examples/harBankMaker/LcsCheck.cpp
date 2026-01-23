#include "LcsCheck.h"
#include "AssignAtomTypes.h"

#include "discamb/BasicUtilities/on_error.h"

#include <fstream>

using namespace discamb;
using namespace std;

LcsCheck::LcsCheck()
{
}

LcsCheck::~LcsCheck()
{
}
void LcsCheck::set()
{
    nlohmann::json data;
    readSettings(data);

    // reads atom types
    string bankFile = data.value("bank file", string(""));
    if (bankFile.empty())
        on_error::throwException("expected \"bank file\" field in setting.json", __FILE__, __LINE__);

    string selectedAtomIdsFile = data.value("selected atom type ids file", string(""));

    AssignAtomTypes::readBank(bankFile, selectedAtomIdsFile, mAtomTypes, mDescriptorsSettings);

}
void LcsCheck::run()
{
    ofstream out(mOutputFileName);

    out << "\nAtoms types with 1 neighbour and potentially colinear lcs directions\n\n";

    for (auto& type : mAtomTypes)
        if (!type.localCoordinateSystem.automaticallyDefined)
            continue;
        else
        {
            if (type.connectivity[0].size() == 1)
            {
                int neighbourIdx = type.connectivity[0][0];

                if (type.atoms[neighbourIdx].nNeighbours<3)
                    out << type.id << "\n";
            }
        }

    out << "\nAtoms types with 2 neighbours\n\n";

    for (auto& type : mAtomTypes)
        if (!type.localCoordinateSystem.automaticallyDefined)
            continue;
        else
        {
            if (type.connectivity[0].size() == 2)
            {
                int neighbourIdx = type.connectivity[0][0];

                if (type.atoms[neighbourIdx].nNeighbours < 3)
                    out << type.id << "\n";
            }
        }

    out << "\nAtoms types with more than 4 neighbours\n\n";
    for (auto& type : mAtomTypes)
        if (!type.localCoordinateSystem.automaticallyDefined)
            continue;
        else
        {
            if (type.connectivity[0].size() > 4)

                out << type.id << "\n";
        }
    
    out.close();
}
