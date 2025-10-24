#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/AtomTyping/TypeMatchAlgorithm.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"
#include "discamb/IO/atom_type_io.h"
#include "discamb/IO/MATTS_BankReader.h"



#include "json.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>

using namespace discamb;
using namespace std;

void type_subtype(const string& bankFile)
{
    vector<AtomType> types;
    vector<nlohmann::json> typeData;
    DescriptorsSettings descriptorsSettings;

    atom_type_io::readAtomTypesJson(bankFile, types, typeData, descriptorsSettings);

    int nTypes = types.size();
    vector<TypeMatchAlgorithm> typeMatchAlgorithms(nTypes);
    for (int i = 0; i < nTypes; i++)
        typeMatchAlgorithms[i].setType(types[i]);
    ofstream out("generalize.log");

    vector<vector<int> > generalizedTypeIdx(nTypes);
    vector<vector<int> > generalizedBy(nTypes);

    for (int i = 0; i < nTypes; i++)
        for (int j = 0; j < nTypes; j++)
            if (i != j)
                if (typeMatchAlgorithms[i].generalize(typeMatchAlgorithms[j]))
                {
                    generalizedTypeIdx[i].push_back(j);
                    generalizedBy[j].push_back(i);
                }
    
    vector<pair<int, int> > equivalentTypes;

    for (int i = 0; i < nTypes; i++)
        for (int otherType: generalizedTypeIdx[i])
        {
            if (find(generalizedTypeIdx[otherType].begin(),
                generalizedTypeIdx[otherType].end(), i) !=
                generalizedTypeIdx[otherType].end())
                    equivalentTypes.push_back({ i, otherType });

        }

    if (equivalentTypes.empty())
        out << "no two types are euqivalent\n\n";
    else
    {
        out << "equivalent type pairs:\n";
        for (auto& p : equivalentTypes)
            out << types[p.first].id << " " << types[p.second].id << "\n";
        out << "\n";
    }

    for (int i = 0; i < nTypes; i++)
    {
        if (generalizedTypeIdx[i].empty())
            out << "type " << types[i].id << " does not generalize any other type\n";
        else
        {
            out << "type " << types[i].id << " generalize types:\n";
            for (int idx: generalizedTypeIdx[i])
                out << "   " << types[idx].id << "\n";
        }
    }
    out.close();
    // type_hierarchy[level] - the higher level the more general types are
    vector<vector<int> > type_hierarchy;
    vector<bool> type_at_lower_level(nTypes, false);
    int nTypesInHierarchy = 0;
    bool all_levels_found = false;
    bool errors_in_hierarchy = false;
    while (!all_levels_found)
    {
        vector<int> level;
        for (int i = 0; i < nTypes; i++)
            if (!type_at_lower_level[i])
            {
                if (generalizedTypeIdx[i].empty())
                {
                    level.push_back(i);
                    type_at_lower_level[i] = true;
                }
                else
                {
                    bool all_generalized_types_at_lower_level = true;
                    for (int generalized_type_idx : generalizedTypeIdx[i])
                        if (!type_at_lower_level[generalized_type_idx])
                            all_generalized_types_at_lower_level = false;
                    if (all_generalized_types_at_lower_level)
                    {
                        level.push_back(i);
                        type_at_lower_level[i] = true;
                    }
                }
            }
        if(!level.empty())
            type_hierarchy.push_back(level);
        else{
            all_levels_found = true;
            if (nTypesInHierarchy == nTypes)
                errors_in_hierarchy = true;
        }
    }
    out.open("type_hierarchy.log");
    if (errors_in_hierarchy)
        out << "errors found in type hierarchy\n\n";
    out << "type hierarchy (the higher level the more general types are):\n";
    for (int l = 0; l < type_hierarchy.size(); l++)
    {
        out << "level " << l << ":\n";
        for (int type_idx : type_hierarchy[l])
            out << "   " << types[type_idx].id << "\n";
    }
    out.close();

    out.open("type_hierarchy_with_defs.log");
    if (errors_in_hierarchy)
        out << "errors found in type hierarchy\n\n";
    out << "type hierarchy (the higher level the more general types are):\n";
    for (int l = 0; l < type_hierarchy.size(); l++)
    {
        out << "#################################################################\n\n";
        out << "level " << l << ":\n\n";
        out << "#################################################################\n\n";

        for (int type_idx : type_hierarchy[l])
        {
            out << "   type " << types[type_idx].id << "\n";

            out << typeData[type_idx].dump(4) << "\n\n";
            if (l > 0)
            {
                out << "GENERALIZE TYPES:\n\n";
                for (int generalized_type_idx : generalizedTypeIdx[type_idx])
                {
                    if (find(type_hierarchy[l - 1].begin(), type_hierarchy[l - 1].end(), generalized_type_idx)!= type_hierarchy[l - 1].end())
                    {
                        out << "      type " << types[generalized_type_idx].id << "\n";
                        out << typeData[generalized_type_idx].dump(4) << "\n\n";
                    }
                }
            }
            out << "\n-------------------------------------------------\n\n";
        }
    }
    out.close();


    //out.open("type_generalization_with_defs.log");
    //out << "only types generaation with definitions:\n\n";
    out.open("types_and_named_rings");
    map<int, vector<int> > named_ring_count_to_type_idx;
    for (int i = 0; i < nTypes; i++)
        named_ring_count_to_type_idx[types[i].ringLabels.size()].push_back(i);
    
    for (auto& type_set : named_ring_count_to_type_idx)
    {
        out << "types with " << type_set.first << " named rings:\n";
        for (int type_idx : type_set.second)
        {
            out << "   " << types[type_idx].id;
            for(auto &label: types[type_idx].ringLabels)
                out << " " << label;
            out << "\n";
        }
        out << "\n\n";
    }

    out.close();


    // generalization via different branches
    vector<int> typeLevel(nTypes);
    for (int l = 0; l < type_hierarchy.size(); l++)
        for (int i = 0; i < type_hierarchy[l].size(); i++)
            typeLevel[type_hierarchy[l][i]] = l;
    out.open("multigeneralization");
    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        //generalizedTypeIdx[i].push_back(j);
        if (generalizedBy[typeIdx].size() > 1)
        {
            vector<pair<int,int> > levels_and_types;
            for (int i = 0; i < generalizedBy[typeIdx].size(); i++)
                levels_and_types.push_back({ typeLevel[generalizedBy[typeIdx][i]],generalizedBy[typeIdx][i] });
            
            sort(levels_and_types.begin(), levels_and_types.end());
            bool type_name_printed = false;

            for(int j=1;j< levels_and_types.size();j++)
                if (levels_and_types[j].first == levels_and_types[j - 1].first)
                {
                    if (!type_name_printed)
                    {
                        out << types[typeIdx].id << " generalized with same level types:" << "\n";
                        out << "    " << types[levels_and_types[j-1].second].id;
                        type_name_printed = true;
                    }
                    out << "    " << types[levels_and_types[j].second].id;
                }
            if (type_name_printed)
                out << "\n";
        }

    }
    out.close();


    // cross-testing
//sortTypesByGenarality
    vector<vector<string> > type_hierarchy_str(type_hierarchy.size());

    for (int l = 0; l < type_hierarchy.size(); l++)
        for (int i = 0; i < type_hierarchy[l].size(); i++)
            type_hierarchy_str[l].push_back(types[type_hierarchy[l][i]].id);


    vector<vector<string> > type_hierarchy_str2(type_hierarchy.size());
    atom_typing_utilities::sortTypesByGenarality(types);
    int type_count = 0;
    for (int l = 0; l < type_hierarchy.size(); l++)
        for (int i = 0; i < type_hierarchy[l].size(); i++)
            type_hierarchy_str2[l].push_back(types[type_count++].id);

    bool hierarchiesTheSame = true;
    for (int l = 0; l < type_hierarchy.size(); l++)
    {
        sort(type_hierarchy_str[l].begin(), type_hierarchy_str[l].end());
        sort(type_hierarchy_str2[l].begin(), type_hierarchy_str2[l].end());
        if (type_hierarchy_str[l] != type_hierarchy_str2[l])
            hierarchiesTheSame = false;

    }
    cout << "hierarchies generated with two algorithms are the same: " << boolalpha << hierarchiesTheSame << "\n";

}



int main(int argc, char *argv[])
{


    try {
        try {

            if (argc != 2)
                on_error::throwException("expected dabank file as input", __FILE__, __LINE__);
            type_subtype(argv[1]);
        }
        catch (nlohmann::json::parse_error& e)
        {
            std::cout << "message: " << e.what() << '\n'
                << "exception id: " << e.id << '\n'
                << "byte position of error: " << e.byte << std::endl;

            stringstream ss;

            ss << "message: " << e.what() << '\n'
               << "exception id: " << e.id << '\n'
               << "byte position of error: " << e.byte << std::endl;

            string message = "Error when parsing JSON file\n" + ss.str();
            on_error::throwException(message, __FILE__, __LINE__);
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
        ofstream out("discamb_error.log" , ostream::app);
        out << e.what();
        out.close();
    }
}


