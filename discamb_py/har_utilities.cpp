#include "har_utilities.h"
#include "discamb/Scattering/gar_utilities.h"

#include <iostream>

using namespace std;

namespace har_utilities {

    std::vector<Representative>  find_default_representatives(
        const CrystalStructure& crystal_structure,
        const std::vector<std::vector<std::vector<int> > >& fragments)
        {
            std::vector<Representative> rerpresenataives;

            int nFragments = fragments.size();
            vector<vector<pair<string, string> > > subsystemAtoms(nFragments);
            for (int fragIdx = 0; fragIdx < nFragments; fragIdx++)
            {
                for (auto const& atom : fragments[fragIdx])
                {
                    discamb::UnitCellContent::AtomID atomId(atom[0], discamb::Vector3i(atom[1], atom[2], atom[3]));
                    string label, symmOpStr;
                    crystal_structure.getUnitCellContent().interpreteAtomID(atomId, label, symmOpStr);
                    subsystemAtoms[fragIdx].push_back({ label, symmOpStr });
                }
            }

            vector<vector<discamb::AtomRepresentativeInfo> > reps;
            discamb::gar_utilities::findDefaultRepresentatives(
                crystal_structure.getCrystal(), 
                subsystemAtoms,
                reps,
                false);
            for (int atomIdx=0; atomIdx< reps.size(); atomIdx++)
            {
                for (auto& rep : reps[atomIdx])
                {
                    Representative representative;
                    representative.atomIdx = rep.idxInSubsystem;
                    representative.substructureIdx = rep.fragmentIdx;
                    representative.weight = rep.fixedWeightValue;
                    rerpresenataives.push_back(representative);
                }
            }
            return rerpresenataives;
        }

}
