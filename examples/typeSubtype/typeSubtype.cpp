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



    //for (int i = 0; i < nTypes; i++)
    //{
    //    vector<string> generalizedTypes;
    //    for (int j = 0; j < nTypes; j++)
    //        if (i != j)
    //            if (typeMatchAlgorithms[i].generalize(typeMatchAlgorithms[j]))
    //            {
    //                generalizedTypes.push_back(types[j].id);
    //                generalizedTypeIdx[i].push_back(j);
    //            }
    //    if (generalizedTypes.empty())
    //        out << "type " << types[i].id << " does not generalize any other type\n";
    //    else
    //    {
    //        out << "type " << types[i].id << " generalize types:\n";
    //        for (string& generalizedType : generalizedTypes)
    //            out << "   " << generalizedType << "\n";
    //    }
    //}
    out.close();
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


