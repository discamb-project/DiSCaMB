#include "discamb/AtomTyping/TypeTree.h"

#include "discamb/AtomTyping/TypeMatchAlgorithm.h"
#include "discamb/BasicUtilities/OnError.h"

#include <algorithm>

using namespace std;

namespace discamb{

    void TypeTree::set(
        const std::vector<AtomType>& types)
    {
        nodes.clear();
        roots.clear();

        int nTypes = types.size();
        vector<TypeMatchAlgorithm> typeMatchAlgorithms(nTypes);
        for (int i = 0; i < nTypes; i++)
            typeMatchAlgorithms[i].setType(types[i]);
        
        map<string, vector<string> > generalizedTypes;

        for (int i = 0; i < nTypes; i++)
        {
            generalizedTypes[types[i].id] = vector<string>();
            for (int j = 0; j < nTypes; j++)
                if (i != j)
                    if (typeMatchAlgorithms[i].generalize(typeMatchAlgorithms[j]))
                        generalizedTypes[types[i].id].push_back(types[j].id);
        }
        

        // each node is generalized by number of nodes which 
        // corresponds to its depth in tree
        map<string, vector<string> > generalizingTypes;

        for (int i = 0; i < nTypes; i++)
            generalizingTypes[types[i].id] = vector<string>();
        for (auto& item : generalizedTypes)
        {
            string typeLabel = item.first;
            vector<string>& generalizedTypes = item.second;
            for (string& type : generalizedTypes)
                generalizingTypes[type].push_back(typeLabel);
        }

        for (auto& item : generalizingTypes)
        {
            string typeName = item.first;
            vector<pair<int, string> > depth_label;
            for (auto generalizingType : item.second)
                depth_label.push_back({ generalizingTypes[generalizingType].size(), generalizingType });
            sort(depth_label.begin(), depth_label.end());
            for (int i = 1; i < depth_label.size(); i++)
                if (depth_label[i].first == depth_label[i - 1].first)
                {
                    string message = typeName + " generalized by two types: " + depth_label[i].second +
                        " and " + depth_label[i - 1].second;
                    on_error::throwException(message, __FILE__, __LINE__);
                }

            nodes[typeName].depth = item.second.size();
            nodes[typeName].id = typeName;
            if (!depth_label.empty())
                nodes[typeName].parent = depth_label.back().second;
        }

        for (auto& node : nodes)
        {
            if (node.second.parent.empty())
                roots.push_back(node.second.id);
            else
                nodes[node.second.parent].children.insert(node.first);
        }


    }

}