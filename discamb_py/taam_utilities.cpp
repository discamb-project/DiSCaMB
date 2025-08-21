#include "taam_utilities.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"

using namespace discamb;
using namespace std;

namespace {
    void readBank(
        const string& bank_path, 
        vector<AtomType> &types,
        BankSettings &bankSettings)
    {
        MATTS_BankReader bankReader;
        vector<AtomTypeHC_Parameters> hcParameters;
        bankReader.read(bank_path, types, hcParameters, bankSettings);
    }
}

namespace taam_utilities {

    std::vector<int> find_atom_types(
        const CrystalStructure& crystal_structure,
        const std::string& bank_path)
    {
        vector<AtomType> types;
        BankSettings bankSettings;

        readBank(bank_path, types, bankSettings);

        CrystalAtomTypeAssigner assigner;
        assigner.setAtomTypes(types);
        assigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
        vector<int> typeId;
        vector<LocalCoordinateSystem<AtomInCrystalID> > lcs;
        assigner.assign(crystal_structure.getCrystal(), typeId, lcs);
        return typeId;
    }

    std::vector<std::vector<int> > find_atom_types_for_fragments(
        const CrystalStructure& crystal_structure,
        const std::string& bank_path,
        const std::vector< std::vector<std::vector<int> > >& fragments)
    {
        vector<AtomType> types;
        BankSettings bankSettings;

        readBank(bank_path, types, bankSettings);

        int fragmentIdx, nFragments = fragments.size();
        vector < vector <pair<int, string> > > fragmentAtoms(nFragments);
        Crystal const& crystal = crystal_structure.getCrystal();
        map<string, int> label2idx;
        for (int i = 0; i < crystal.atoms.size(); i++)
            label2idx[crystal.atoms[i].label] = i;

        vector<vector<int> > atomsToAssign(nFragments);

        for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
        {
            int atomIdx, nAtoms = fragments[fragmentIdx].size();
            //atomsToAssign[fragmentIdx].resize(nAtoms, true);
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                atomsToAssign[fragmentIdx].push_back(atomIdx);
                string label = crystal_structure.getAtomLabel(fragments[fragmentIdx][atomIdx][0]);
                
                pair<int, string> atom{ label2idx[label], crystal_structure.getSymmetryOperationStr(fragments[fragmentIdx][atomIdx]) };
                fragmentAtoms[fragmentIdx].push_back(atom);
            }
        }

        CrystalAtomTypeAssigner assigner;
        assigner.setAtomTypes(types);
        assigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
        vector< vector<int> > typeIdx;
        vector< vector<LocalCoordinateSystem<AtomInCrystalID> > > lcs;
        assigner.assign(crystal_structure.getCrystal(), fragmentAtoms, atomsToAssign, typeIdx, lcs);
        return typeIdx;
    }


    std::vector < std::string > get_atom_type_names(const std::string& bank_path)
    {
        vector<AtomType> types;
        BankSettings bankSettings;

        readBank(bank_path, types, bankSettings);
        vector<string> ids;
        for (auto& type : types)
            ids.push_back(type.id);
        return ids;
    }

}
