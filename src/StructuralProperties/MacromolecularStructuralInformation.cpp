#include "discamb/StructuralProperties/MacromolecularStructuralInformation.h"

#include "discamb/BasicUtilities/string_utilities.h"

using namespace std;

namespace discamb {

    void MacromolecularStructuralInformation::set(
        const nlohmann::json& data) 
    {
        connectivity.clear();
        planes.clear();
        altlocs.clear();

        if (data.find("atom data") != data.end() && data["atom data"].is_array()) {
            for (const auto& atom_data : data["atom data"]) {
                if (atom_data.find("altloc") != atom_data.end() && atom_data["altloc"].is_string()) {
                    std::string altloc_str = atom_data["altloc"].get<std::string>();
                    if (!altloc_str.empty())
                        altlocs.push_back(altloc_str[0]);
                }
                if (atom_data.find("plane") != atom_data.end() && atom_data["altloc"].is_array()) {
                    vector<pair<int, string> > plane_list;
                    for (const auto& atom_in_plane : atom_data["plane"]) {
                        if (atom_in_plane.is_number_integer())
                            plane_list.push_back(make_pair<int, string>(atom_in_plane.get<int>(), string("x,y,z")));
                        if (atom_in_plane.is_string()) {
                            vector<string> words;
                            string_utilities::split(atom_in_plane.get<string>(), words, ',');
                            if (words.size() == 2) {
                                plane_list.push_back({ stoi(words[0]), words[1] });
                            }
                            else
                                on_error::throwException("Error parsing atom plane information from JSON data.", __FILE__, __LINE__);
                        }

                    }
                    planes.push_back(plane_list);
                }
            }
        }

    }
}