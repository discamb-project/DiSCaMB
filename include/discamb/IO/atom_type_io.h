#pragma once

#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"


#include "json.hpp"


namespace discamb {
    namespace atom_type_io {

        void readAtomTypes(const std::string& fileName, std::vector<AtomType>& atomTypes, DescriptorsSettings& descriptorsSettings);

        void readAtomTypesJson(const std::string& jsonFileName, std::vector<AtomType> &atomTypes, DescriptorsSettings &descriptorsSettings);
        void readAtomTypesJson(const std::string& jsonFileName, std::vector<AtomType> &atomType, std::vector<nlohmann::json> &typeData, DescriptorsSettings& descriptorsSettings);
        
        void json2atomType(
            const nlohmann::json &data, 
            const std::map<std::string, std::set<int> >& element_groups, 
            const std::string atomTypeId,
            AtomType& atomType,
            const AtomDescriptors &defaultAtomDescriptors,
            const AtomDescriptors &defaultCetralAtomDescriptors);

        //void parseStructuralFormulaTxt(const std::vector<std::string>& lines, std::vector<AtomDescriptors> &atoms, std::vector<std::vector<int> > &connectivity);
        void parseStructuralFormulaTxt(
            const std::vector<std::string>& lines,
            const std::map<std::string, std::set<int> >& element_groups,
            const AtomDescriptors& defaultDescriptors,
            const AtomDescriptors& centralAtomDescriptors,
            std::vector<AtomDescriptors>& atoms,
            std::vector<std::vector<int> >& connectivityFinal);

        void json2atomDescriptors(const nlohmann::json& data, AtomDescriptors& atomDescriptors);
        void json2elementGroups(const nlohmann::json& data, std::map<std::string, std::set<int> > &elementGroups);

        void writeAtomTypesTxt(
            const std::string &fileName, 
            const std::vector<AtomType>& atomTypes,
            const std::map<std::string, std::set<int> >& element_groups = std::map<std::string, std::set<int> >());
    }
}

