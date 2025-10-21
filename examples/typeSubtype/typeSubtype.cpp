#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/AtomTyping/TypeMatchAlgorithm.h"
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
    DescriptorsSettings descriptorsSettings;

    atom_type_io::readAtomTypes(bankFile, types, descriptorsSettings);

    int nTypes = types.size();
    vector<TypeMatchAlgorithm> typeMatchAlgorithms(nTypes);
    for (int i = 0; i < nTypes; i++)
        typeMatchAlgorithms[i].setType(types[i]);
    ofstream out("generalize.log");

    vector<vector<int> > generalizedTypeIdx(nTypes);

    for (int i = 0; i < nTypes; i++)
        for (int j = 0; j < nTypes; j++)
            if (i != j)
                if (typeMatchAlgorithms[i].generalize(typeMatchAlgorithms[j]))
                    generalizedTypeIdx[i].push_back(j);
    
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
    //out.open("type_generalization_with_defs.log");
    //out << "only types generaation with definitions:\n\n";
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


