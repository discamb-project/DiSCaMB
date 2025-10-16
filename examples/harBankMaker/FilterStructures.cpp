#include "FilterStructures.h"

#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/shelx_io.h"

#include "json.hpp"


#include <fstream>

using namespace std;
using namespace discamb;

FilterStructures::FilterStructures()
{

}

FilterStructures::~FilterStructures()
{

}

void FilterStructures::set()
{
    nlohmann::json data;
    readSettings(data);

    mChosenResFolder = data.value("chosen res folder", "chosen");
    
    mHighAdpsStructurePercent = data.value("% structures with too high Ueq", mHighAdpsStructurePercent);
       
    mMaxUeqivalent = data.value("max Ueq", mMaxUeqivalent);

    string folder;
    if (mResFolder.empty())
        folder = filesystem::current_path().string();
    else
        folder = mResFolder.string();

    mResFolder = data.value("res folder", folder);

    cout << "Filtering out structures, settings:\n"
        << "  max Ueq - " << mMaxUeqivalent << endl
        << "  % structures with highest Ueq to be removed - " << mHighAdpsStructurePercent << "\n"
        << "  res folder - " << mResFolder.string() << "\n";
}

void FilterStructures::run()
{
    vector<string> files;
    vector<string> final_structures;
    file_system_utilities::find_files("res", mResFolder.string(), files, false);

    cout << __LINE__ << "\n";

    ofstream log("filter.log");

    vector<double> u_eqivalent, v;
    vector<pair<double, int> > max_u_eq;
    // structure_replicas["GLYALA"] vector of indices of max_u_eq
    // refering to "GLYALA" entries
    map<string, vector<pair<double, int> > > structure_replicas;
    int nStructuresUani = 0;

    

    log << "removed structures:\n";

    for (int idx=0; idx< files.size(); idx++)
    {
        cout << '\r' << files[idx] << "\n";
        auto& file = files[idx];

        Crystal crystal;
        shelx_io::read((mResFolder  / file).string(), crystal);
        
        // skip if u_iso for non-H atom
        for (auto& atom : crystal.atoms)
            if (atom.type != "H")
                if (atom.adp.size() != 6)
                {
                    log << file << " non-H atom with Uiso\n";
                    continue;
                }

        nStructuresUani++;

        // find U_eq of non-H atoms and find the largest one

        crystal_structure_utilities::u_eq(crystal, v);
        u_eqivalent.clear();
        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            if (crystal.atoms[atomIdx].type != "H")
                u_eqivalent.push_back(v[atomIdx]);

        double u_max = *max_element(u_eqivalent.begin(), u_eqivalent.end());

        if (u_max > mMaxUeqivalent)
        {
            log << file << " max U_eq exceeded\n";
            continue;
        }

        max_u_eq.push_back({ u_max, idx });
        structure_replicas[files[idx].substr(0, 6)].push_back(max_u_eq.back());
    }

    // remove multiple occurences of the same structure

    vector<pair<double, int> > max_u_eq_tmp;
    for (auto& group : structure_replicas)
    {
        auto it = min_element(group.second.begin(), group.second.end());
        max_u_eq_tmp.push_back(*it);
        for (auto iter = group.second.begin(); iter != group.second.end(); iter++)
            if (it != iter)
                log << files[iter->second] << " removed repeating structure with higher max U_eq\n";
    }

    int nRepeatingRemoved = max_u_eq.size() - max_u_eq_tmp.size();

    max_u_eq_tmp.swap(max_u_eq);

    // n0 - size of initial set of structures when repeats are removed
    // and structures with U-iso for non-H atoms are removed
    int n0 = nStructuresUani - nRepeatingRemoved;
    // n removed due to U_eq larger than threshold
    int nRemoved = nStructuresUani - max_u_eq.size();
    int totalNtoRemove = int(n0 * mHighAdpsStructurePercent/100.0);
    int nStillToRemove = totalNtoRemove - nRemoved;
    if (nStillToRemove > 0)
    {
        sort(max_u_eq.begin(), max_u_eq.end());
        int nToKeep = max_u_eq.size() - nStillToRemove;
        for (int i = nToKeep; i < max_u_eq.size(); i++)
            log << files[max_u_eq[i].second] << " large U_eq group\n";
        max_u_eq.resize(nToKeep);
    }

    log << "structures to use:\n";

    for (auto& item : max_u_eq)
    {
        final_structures.push_back(files[item.second]);
        log << files[item.second] << "\n";
    }
    log.close();

    ofstream out("copy.bat");
    if (filesystem::exists(mChosenResFolder))
        filesystem::remove_all(mChosenResFolder);
    filesystem::create_directory(mChosenResFolder);
    
    for (auto& file : final_structures)
        out << "copy " << (mResFolder / file).string() << " " << mChosenResFolder << "\n";
    out.close();
}


