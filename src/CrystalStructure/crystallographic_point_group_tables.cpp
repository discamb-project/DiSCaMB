#include "discamb/CrystalStructure/crystallographic_point_group_tables.h"
#include "discamb/CrystalStructure/SpaceGroupOperation.h"
#include "discamb/BasicUtilities/string_utilities.h"

#include <map>

using namespace std;

namespace discamb {
    namespace {
    
        map<string, vector<string> > point_groups = {
            {"1", {"x,y,z"}},
            {"-1", {"x,y,z", "-x,-y,-z"}},
            {"2", {"x,y,z", "-x,y,-z"}},
            {"m", {"x,y,z", "x,-y,z"}},
            {"2/m", {"x,y,z", "-x,y,-z", "x,-y,z", "-x,-y,-z"}},
            {"222", {"x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z"}}
            
        };
    
    }
    namespace crystallographic_point_group_tables {
        std::vector<Matrix3i> getPointGroupOperations(const std::string& pointGroupSymbol) {
            // Implementation to retrieve the point group operations based on the symbol

            return {};
        }
        std::vector<std::string> getPointGroupOperationsStr(const std::string& pointGroupSymbol) {
            // Implementation to retrieve the point group operations as strings based on the symbol
            // This is a placeholder for the actual implementation
            return {};
        }
        std::vector<std::string> getAvailablePointGroups() {
            // Implementation to retrieve all available point groups
            // This is a placeholder for the actual implementation
            return {};
        }

        std::string findPointGroup(const std::vector<Matrix3i>& symmOps)
        {
            // Implementation to find the point group based on symmetry operations

            vector<Matrix3i> symmOpsCopy = symmOps;
            sort(symmOpsCopy.begin(), symmOpsCopy.end());
            int nSymmOps = symmOpsCopy.size();
            for (const auto& point_group : point_groups) {
                if (nSymmOps == point_group.second.size())
                {
                    vector<Matrix3i> pgOps = getPointGroupOperations(point_group.first);
                    sort(pgOps.begin(), pgOps.end());
                    if (pgOps == symmOpsCopy)
                        return point_group.first;
                }
            }
            return "";
        }

        std::string findPointGroup(
            const vector<Matrix3i>& symmOps,
            vector<int>& canonicalOrder)
        {
            int nSymmOps = symmOps.size();
            canonicalOrder.assign(nSymmOps, -1);

            for (const auto& point_group : point_groups) {
                if (point_group.second.size() == nSymmOps)
                {
                    vector<Matrix3i> pg_operations(nSymmOps);
                    for(int i=0;i<nSymmOps;i++)
                    {
                        SpaceGroupOperation spaceGroupOp(point_group.second[i]);
                        Matrix3i rotationMatrix;
                        Vector3<CrystallographicRational> translationVector;
                        spaceGroupOp.get(rotationMatrix, translationVector);
                        pg_operations[i] = rotationMatrix;
                    }

                    vector<bool> symmOpFound(nSymmOps, false);
                    bool foundPointGroup = true;
                    for (int idxInTable = 0; idxInTable < nSymmOps; idxInTable++)
                    {
                        SpaceGroupOperation spaceGroupOp(point_group.second[idxInTable]);
                        Matrix3i rotationMatrix;
                        Vector3<CrystallographicRational> translationVector;
                        spaceGroupOp.get(rotationMatrix, translationVector);
                        Matrix3i pg_operation_table = rotationMatrix;

                        bool foundOperation = false;
                        for (int idxInArg = 0; idxInArg < nSymmOps; idxInArg++)
                            if (symmOps[idxInArg] == pg_operation_table)
                            {
                                foundOperation = true;
                                canonicalOrder[idxInArg] = idxInTable;
                            }
                        
                        if (!foundOperation)
                            foundPointGroup = false;
                    }
                    if (foundPointGroup)
                        return point_group.first;
                }
            }
            return string();
        }

        std::string findPointGroup(
            const std::vector<std::string>& symmetryOperations,
            std::vector<int>& canonicalOrder)
        {
            vector<Matrix3i> symmOps;
            for (auto const& op : symmetryOperations) {
                SpaceGroupOperation spaceGroupOp(op);
                Matrix3i rotationMatrix;
                Vector3<CrystallographicRational> translationVector;
                spaceGroupOp.get(rotationMatrix, translationVector);
                symmOps.push_back(rotationMatrix);
            }
            return findPointGroup(symmOps, canonicalOrder);

        }

        std::string findPointGroup(const std::vector <std::string> &symmetryOperation)
        {
            // Implementation to find the point group based on symmetry operations
            
            vector<Matrix3i> symmOps;
            for (auto const& op : symmetryOperation) {
                SpaceGroupOperation spaceGroupOp(op);
                Matrix3i rotationMatrix;
                Vector3<CrystallographicRational> translationVector;
                spaceGroupOp.get(rotationMatrix, translationVector);
                symmOps.push_back(rotationMatrix);
            }
            return findPointGroup(symmOps);
        }

        int pointGroupIdx(const std::string& pointGroupSymbol) {
            auto pointGroups = getAvailablePointGroups();
            for (int i = 0; i < pointGroups.size(); ++i) {
                if (pointGroups[i] == pointGroupSymbol) {
                    return i;
                }
            }
            return -1; // Not found
        }

    }
}
