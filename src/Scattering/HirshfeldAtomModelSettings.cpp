#include "discamb/Scattering/HirshfeldAtomModelSettings.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/IO/fragmentation_io.h"
#include "discamb/IO/discamb_io.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/Scattering/gar_utilities.h"

#include <optional>
#include <set>

using namespace std;

namespace discamb {

    namespace ham_settings {

        namespace {
            optional<int> indexOfContainingSet(const vector<vector<pair<string, string> > >& sets, const string& s)
            {
                optional<int> result;
                for (int i = 0; i < int(sets.size()); i++)
                    for (auto const& item : sets[i])
                        if (item.first == s)
                            return optional<int>(i);
                
                return optional<int>();
            }

        }


        void PartitionData::set(
            const nlohmann::json& data)
        {
            //-----------
            // set type
            //-----------

            string harTypeString = data.value("HAR type", "Hirshfeld");
            //alternatively
            harTypeString = data.value("partition", harTypeString);
            type = stockholderPartitionTypeFromString(harTypeString);

            if (type == ElectronDensityPartitionType::none)
                on_error::throwException("invalid name of electron density partition: " + harTypeString, __FILE__, __LINE__);

            //-----------
            // set partitionSpecificData
            //-----------

            partitionSpecificData.clear();

            if (data.find("partition settings") != data.end())
                partitionSpecificData = data.find("partition settings").value();

            //-----------
            // set power
            //-----------

            power = data.value("power", 1.0); 

            // set edCalcAtomIncludeRange

            if (data.contains("density calculation range"))
                edCalcAtomIncludeRange = data.find("density calculation range")->get<double>();
            else
                edCalcAtomIncludeRange.reset();

        }

        
        void WfnCalculation::set(
            const nlohmann::json& data)
        {
            *this = WfnCalculation();

            //-- qmProgram

            if (data.find("qm program") != data.end())
                qmProgram = WaveFunctionDataGeneratorRunner::string2program(data["qm program"].get<string>());
            else
                on_error::throwException("\"qm program\" field missing when defining wave function calculation settings", __FILE__, __LINE__);

            //-- qmProgramSpecificData

            

            if (data.find("qm program settings") != data.end())
                qmProgramSpecificData = data["qm program settings"];

            if (data.find("qm folder") != data.end())
                qmProgramSpecificData["folder"] = data["qm folder"].get<string>();
            if (data.find("molden2aim folder") != data.end())
                qmProgramSpecificData["molden2aim folder"] = data["molden2aim folder"].get<string>();
            //if (data.find("qm template file") != data.end())
            //{
            //    string templateString;
            //    file_system_utilities::file2string(data["qm template"].get<string>(), templateString);
            //    data["qm template"] = templateString;
            //}
                        
            //-- method, relativisticMethod, basisSet
            
            qmSettings.set(data);

            // multipoles



            //-- restart & restartStep

            string restartOptionStr = data.value("restart", "from scratch");

            map<string, Restart> str2restartOption = { {"from scratch", Restart::from_scratch}, 
                {"from converged wfn", Restart::from_converged_wfn}, {"from converged density", Restart::from_converged_density},
                {"from not converged wfn", Restart::from_converged_wfn}, {"from not converged density", Restart::from_not_converged_density},
                {"from converged multipoles", Restart::from_converged_multipoles}, 
                {"from not converged multipoles", Restart::from_not_converged_multipoles} };

            string validNamesStr;
            for (auto item : str2restartOption)
                validNamesStr += item.first + "\n";

            if (str2restartOption.find(restartOptionStr) == str2restartOption.end())
                on_error::throwException(string("invalid restart option '") + restartOptionStr +
                    string("' for electron density calculation, valid options:\n") + validNamesStr,
                    __FILE__, __LINE__);
            else
                restart = str2restartOption[restartOptionStr];

            if (restart != Restart::from_scratch)
            {
                if (data.find("restart step") == data.end())
                    on_error::throwException("missing entry 'restart step' in wave function calculation specification", __FILE__, __LINE__);

                restartStep = data.find("restart step").value().get<int>();
            }

        };

        void CrystalFragmentWfnCalculation::set(
            const nlohmann::json& _data,
            double clusterThreshold)
        {
            nlohmann::json data = _data;
            
            *this = CrystalFragmentWfnCalculation();
            if (data.find("multipole sites") != data.end())
            {
                //auto const& multipoleSitesData = *data.find("multipole sites");
                auto multipoleSitesData = *data.find("multipole sites");
                if (multipoleSitesData.find("threshold") == multipoleSitesData.end())
                    multipoleSitesData["threshold"] = clusterThreshold;
                distributedMultipoleCluster = DistributedMultipoleCentersSettings();
                distributedMultipoleCluster->set(multipoleSitesData);
            }
        }

        void MultipoleExpansionCalculation::set(
            const nlohmann::json& data)
        {
            *this = MultipoleExpansionCalculation();

            if (data.find("multipoles") != data.end())
            {
                if (!data["multipoles"].is_boolean())
                {
                    auto multipolesData = data.find("multipoles");
                    //multipoleClusterThreshold = multipolesData->value("threshold", multipoleClusterThreshold);
                    clusterThreshold = multipolesData->value("threshold", clusterThreshold);
                    multipoleExpansionLevel = multipolesData->value("l max", multipoleExpansionLevel);
                    multipoleChargeDistance = multipolesData->value("charge distance", multipoleChargeDistance);
                    multipoleConvergenceThreshold = multipolesData->value("convergence", multipoleConvergenceThreshold);
                    lookForMultipoleFile = multipolesData->value("find file", lookForMultipoleFile);
                    maxN_Steps = multipolesData->value("max n steps", maxN_Steps);
                    startWithTaam = multipolesData->value("start with taam", startWithTaam);
                    previousMultipolesFile = multipolesData->value("previous multipoles file", previousMultipolesFile);
                    calculate = true;
                }
                else
                    calculate = data["multipoles"].get<bool>();

            }
            else
                calculate = false;
        }

        void Diagnostics::set(
            const nlohmann::json& data)
        {
            *this = Diagnostics();
            prerun = data.value("prerun", prerun);
            printSubsystemsToXyz = data.value("print xyz", printSubsystemsToXyz);
            printMultipoleClustersToXyz = data.value("print multipole cluster to xyz", printMultipoleClustersToXyz);
            printSubsystemsToMol2 = data.value("print mol2", printSubsystemsToMol2);
            printMultipoleClusterToMol2 = data.value("print multipole cluster to mol2", printMultipoleClusterToMol2);
        }

        

        void FormFactorsCalculation::set(
            const nlohmann::json& data)
        {
            *this = FormFactorsCalculation();          

            integrationGrid.angularGridSize = data.value("angular grid", integrationGrid.angularGridSize);
            integrationGrid.radialGridSize = data.value("radial grid", integrationGrid.radialGridSize);
            string gridStr;
            if (data.find("grid") != data.end())
                gridStr = data.find("grid").value().get<string>();

            vector<string> gridSpecs, words;
            string_utilities::split(gridStr, gridSpecs);

            for (auto item : gridSpecs)
            {
                string_utilities::split(item, words, ',');

                if (words.size() != 3 && words.size() != 2)
                    on_error::throwException(string("invalid specification of integration grid"), __FILE__, __LINE__);

                if (words.size() == 3)
                {
                    std::set<int> atomicNumbers;
                    basic_chemistry_utilities::getElementsList(words[0], atomicNumbers);
                    int nAng = stoi(words[1]);
                    int nRad = stoi(words[2]);
                    for (int z : atomicNumbers)
                    {
                        elementSpecificIntegrationGrid[z].angularGridSize = nAng;
                        elementSpecificIntegrationGrid[z].radialGridSize = nRad;
                    }
                }
                else
                {
                    integrationGrid.angularGridSize = stoi(words[0]);
                    integrationGrid.radialGridSize = stoi(words[1]);
                }
            }

        };



        void setCrystalFragments(
            const nlohmann::json& data,
            const Crystal& crystal,
            std::vector<QmFragmentInCrystal>& crystalFragments)
        {
            crystalFragments.clear();
            string errorMessage = "invalid format of subsystems definition for Hirshfeld Atom Model calculations";
            

            //if (subsystemsJson.empty())
            //{
            //    vector<vector<pair<string, string> > > subsystemAtoms;
            //    vector<string> subsystemLabels;
            //    vector<int> subsystemCharges, subsystemSpinMultiplicity;

            //    gar_utilities::findDeaultSubsystems(
            //        crystal, subsystemAtoms, subsystemLabels,
            //        subsystemCharges, subsystemSpinMultiplicity);
            //    for()
            //}

            /*
            (1) no specification - the structure is automatically divided into subunits 
                               with charge 0 and spin multiplicity 0
            
            (2) subsystems:
            ["C1,1","C2,-1,1,picrate"] - list spcifying chare, spin multiplicity (optional, default = 1)
                                        and name (optional, default = subsystem_n, n enumerates subsystems)
                                        for susbsystems identified by the name one of the atoms they include.
                                        the structure is automatically divided into subunits
            (3) subsystems: "filename" - subsystems are defined in the specified file

            (4) subsystems: 
            [ 
                {
                    "name": "water",
                    "charge": 0,
                    "spin multiplicity": 1,
                    "atoms": "O3 H2 H3"
                },
                {
                    "name": "oxalic_acid",
                    "charge": 0,
                    "spin multiplicity": 1,
                    "atoms": "O1 O2 C1 H1 O1,1-X,1-Y,1-Z O2,1-X,1-Y,1-Z C1,1-X,1-Y,1-Z H1,1-X,1-Y,1-Z"
                 }
            ]
            name, charge and spin multiplicity are optional
            "atoms" can be also specified as an array:
            [
                {"add": "C1 C2"},
                {
                    "connect": "C3",
                    "remove": "C4 C9",
                    "remove and cap": "C5 C6"
                },
                {"connect": "C12 C17"}
            ]
            or object
                {
                    "connect": "C3",
                    "remove": "C4 C9",
                    "remove and cap": "C5 C6"
                }
            */

            nlohmann::json subsystemsData;
            if (data.find("qm structure") != data.end())
                subsystemsData = data["qm structure"];
            if (data.find("subsystems") != data.end())
                subsystemsData = data["subsystems"];

            enum class DefinitionType {AUTOMATIC, FROM_FILE, STRING_ARRAY, OBJECT_ARRAY};

            // detect definition type

            DefinitionType definitionType = DefinitionType::AUTOMATIC;
            //nlohmann::json subsystemsData;
            //if (data.find("subsystems") != data.end())
            if(!subsystemsData.is_null())
            {
                //subsystemsData = data["subsystems"];
                if (subsystemsData.is_string())
                    definitionType = DefinitionType::FROM_FILE;
                else
                {
                    if (!subsystemsData.is_array())
                        on_error::throwException(errorMessage, __FILE__, __LINE__);
                    if (subsystemsData.size() == 0)
                        return;
                    if (subsystemsData.front().is_string())
                        definitionType = DefinitionType::STRING_ARRAY;
                    else
                    {
                        if(subsystemsData.front().is_object())
                            definitionType = DefinitionType::OBJECT_ARRAY;
                        else
                            on_error::throwException(errorMessage, __FILE__, __LINE__);
                    }
                }
            }
            
            // end of detect definition type

            //AUTOMATIC, STRING, STRING_ARRAY, OBJECT_ARRAY
            if (definitionType == DefinitionType::AUTOMATIC || definitionType == DefinitionType::STRING_ARRAY)
            {
                vector<vector<pair<string, string> > > subsystemAtoms;
                vector<string> subsystemLabels;
                vector<int> subsystemCharges;
                vector<int> subsystemSpinMultiplicity;
                fragmentation::intermolecular(crystal, subsystemAtoms, 
                    subsystemLabels, subsystemCharges, subsystemSpinMultiplicity);

                if (definitionType == DefinitionType::STRING_ARRAY)
                {
                    vector<string> words;
                    for (auto it = subsystemsData.begin(); it != subsystemsData.end(); ++it)
                    {
                        string_utilities::split(it->get<string>(), words, ',');
                        optional<int> idx = indexOfContainingSet(subsystemAtoms, words[0]);
                        if (idx)
                        {
                            if (words.size() > 1)
                                subsystemCharges[idx.value()] = stoi(words[1]);
                            if (words.size() > 2)
                                subsystemSpinMultiplicity[idx.value()] = stoi(words[2]);
                            if (words.size() > 3)
                                subsystemLabels[idx.value()] = words[3];
                        }
                    }
                }
                size_t fragIdx, nFrag = subsystemAtoms.size();
                crystalFragments.resize(nFrag);
                for (fragIdx = 0; fragIdx < nFrag; fragIdx++)
                {
                    crystalFragments[fragIdx].atoms.atomList = subsystemAtoms[fragIdx];
                    crystalFragments[fragIdx].charge = subsystemCharges[fragIdx];
                    crystalFragments[fragIdx].label = subsystemLabels[fragIdx];
                    crystalFragments[fragIdx].spin_multiplicity = subsystemSpinMultiplicity[fragIdx];
                }
                return;
            }


            std::vector<FragmentConstructionData> fragmentsConstructionData;

            if (definitionType == DefinitionType::FROM_FILE)
                fragmentation_io::readPlainText(subsystemsData.get<string>(), fragmentsConstructionData);

            if (definitionType == DefinitionType::OBJECT_ARRAY)
            {
                fragmentsConstructionData.resize(subsystemsData.size());
                int subsystemIdx = 0;
                for (auto& sub : subsystemsData)
                {
                    fragmentsConstructionData[subsystemIdx].set(sub);
                    if(fragmentsConstructionData[subsystemIdx].label.empty())
                        fragmentsConstructionData[subsystemIdx].label = "subsystem_" + to_string(subsystemIdx + 1);
                    subsystemIdx++;
                }
            }
            fragmentation::make_qm_fragments(crystal, fragmentsConstructionData, crystalFragments);
            // fragmentPartConstructionData to atom list

        //    UnitCellContent unitCellContent;
        //    unitCellContent.set(crystal);
        //    set<UnitCellContent::AtomID> fragmentAtoms;
        //    set<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > cappingHydrogens;
        //    vector<UnitCellContent::AtomID> atomList;
        //    vector<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > cappingHydrogenList;
        //    vector<vector<UnitCellContent::AtomID> > connectivity;
        //    
        //    structural_properties::calcUnitCellConnectivity(unitCellContent, connectivity, 0.4);
        //    string label, symmetryOperationStr;

        //    for (auto& fragment : fragments)
        //    {
        //        fragmentAtoms.clear();
        //        atomList.clear();
        //        QmFragmentInCrystal subsystem;
        //        subsystem.charge = fragment.charge;
        //        subsystem.label = fragment.label;
        //        subsystem.spin_multiplicity = fragment.spin_multiplicity;

        //        for (auto& fragmentConstructionData : fragment.fragmentPartConstructionData)
        //        {
        //            fragmentConstructionData.getAtomList(unitCellContent, connectivity, atomList, cappingHydrogenList);
        //            fragmentAtoms.insert(atomList.begin(), atomList.end());
        //            cappingHydrogens.insert(cappingHydrogenList.begin(), cappingHydrogenList.end());
        //        }
        //        for (auto& atom : fragmentAtoms)
        //        {
        //            unitCellContent.interpreteAtomID(atom, label, symmetryOperationStr);
        //            subsystem.atoms.atomList.push_back({ label, symmetryOperationStr });
        //        }
        //        for (auto& cappingH : cappingHydrogens)
        //        {
        //            subsystem.atoms.cappingHydrogens.resize(subsystem.atoms.cappingHydrogens.size() + 1);
        //            auto& capH = subsystem.atoms.cappingHydrogens.back();
        //            unitCellContent.interpreteAtomID(cappingH.first, label, symmetryOperationStr);
        //            capH.bondedAtom = label;
        //            capH.bondedAtomSymmOp = symmetryOperationStr;
        //            unitCellContent.interpreteAtomID(cappingH.second, label, symmetryOperationStr);
        //            capH.directingAtom = label;
        //            capH.directingAtomSymmOp = symmetryOperationStr;
        //        }
        //        crystalFragments.push_back(subsystem);
        //    }
        }

        /*
        sets wfnData to default value and then sets 
        its members 
        */

        void wfnCalcData(
            const WfnCalculation& wfnCalc,
            WaveFunctionCalculationData& wfnData)
        {
            wfnData = WaveFunctionCalculationData();
            if (wfnCalc.hardwareSettings)
                wfnData.hardware = wfnCalc.hardwareSettings.value();
            wfnData.qmProgramSpecificData = wfnCalc.qmProgramSpecificData;
            wfnData.qmSettings = wfnCalc.qmSettings;
        }


        //void subsystemWfnCalcData(
        //    const Crystal& crystal,
        //    const Subsystem& subsystem,
        //    const WfnCalculation& wfnCalc,
        //    const SubsystemWfnCalculation& subsystemWfnCalc,
        //    WaveFunctionCalculationData& wfnData)
        //{
        //    wfnData = WaveFunctionCalculationData();
        //    bool hasSubsystemWfnSettings = subsystemWfnCalc.wfnCalculation.has_value();
        //    const auto & wfnCalcSettings = (hasSubsystemWfnSettings ? subsystemWfnCalc.wfnCalculation.value() : wfnCalc);

        //    wfnCalcData(wfnCalcSettings, wfnData);
        //    wfnData.atomIdx2BasisSetMap = subsystemWfnCalc.atomicIdx2BasisSetMap;

        //    if(wfnCalc.hardwareSettings)
        //        wfnData.hardware = wfnCalc.hardwareSettings.value();
        //    if (hasSubsystemWfnSettings)
        //        if (subsystemWfnCalc.wfnCalculation.value().hardwareSettings)
        //            wfnData.hardware = subsystemWfnCalc.wfnCalculation.value().hardwareSettings.value();

        //    wfnData.jobName = subsystem.label;
        //    vector<ChemicalElement> chemicalElements;
        //    vector<Vector3d> positions;
        //    vector<int> atomicNumbers;
        //    subsystem.toXyzMol(crystal, chemicalElements, positions);
        //    for (auto& element : chemicalElements)
        //        atomicNumbers.push_back(element.atomicNumber());
        //    wfnData.qmSystem.positions = positions;
        //    wfnData.qmSystem.charge = subsystem.charge;
        //    
        //    //wfnData.qmSystem.pointChargePosition
        //}

        void setRepresentativesFromJsonArray(
            const nlohmann::json& json_array,
            const Crystal& crystal,
            const std::vector<QmFragmentInCrystal>& subsystems,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives)
        {
            map<string, int> subsystemLabel2idx, atomLabel2idx;
            
            representatives.clear();
            representatives.resize(crystal.atoms.size());
            for (int i = 0; subsystems.size(); i++)
                subsystemLabel2idx[subsystems[i].label] = i;
            for (int i = 0; i < crystal.atoms.size(); i++)
                atomLabel2idx[crystal.atoms[i].label] = i;

            for (auto const& item : json_array)
            {
                int susbsystem_idx;
                if(item.find("label")==item.end())
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, missing 'label' field",
                        __FILE__, __LINE__);
                if (subsystemLabel2idx.find(item["label"].get<string>()) == subsystemLabel2idx.end())
                    {
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, unknown subsystem label: " +
                        item["label"].get<string>(), __FILE__, __LINE__);
                }
                
                susbsystem_idx = subsystemLabel2idx[item["label"].get<string>()];
                
                vector<double> weights;

                if (item.find("atom weights") == item.end())
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, missing 'atom weights' field",
                        __FILE__, __LINE__);
                if (!item["atom weights"].is_array())
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, expected array for 'atom weights' field",
                        __FILE__, __LINE__);

                
                for (auto const& weight_item : item["atom weights"])
                    if(!weight_item.is_number_float())
                        on_error::throwException(
                            "invalid format of atom representatives definition in json data, expected float for 'atom weights' field",
                            __FILE__, __LINE__);
                    else
                        weights.push_back(weight_item.get<double>());
                if(weights.size() != subsystems[susbsystem_idx].atoms.atomList.size())
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, number of weights does not match number of atoms in subsystem",
                        __FILE__, __LINE__);
                for (int idxInSubsystem = 0; idxInSubsystem < weights.size(); idxInSubsystem++)
                {
                    AtomRepresentativeInfo atomRepresentativeInfo;
                    atomRepresentativeInfo.atomLabel = subsystems[susbsystem_idx].atoms.atomList[0].first;
                    atomRepresentativeInfo.fixedWeightValue = weights[susbsystem_idx];
                    atomRepresentativeInfo.fragmentIdx = susbsystem_idx;
                    atomRepresentativeInfo.isWeightFixed = true;
                    atomRepresentativeInfo.symmetryCode = subsystems[susbsystem_idx].atoms.atomList[0].second;
                    atomRepresentativeInfo.symmetryCode = "invert";
                    int atom_idx = atomLabel2idx[atomRepresentativeInfo.atomLabel];
                    representatives[atom_idx].push_back(atomRepresentativeInfo);
                }
            }
        }


        void setRepresentatives(
            const nlohmann::json& data,
            const Crystal& crystal,
            const std::vector<QmFragmentInCrystal>& subsystems,
            std::vector<std::vector<AtomRepresentativeInfo> > &representatives)
        {
            representatives.clear();

            if(data.find("atom representatives") != data.end())
            {
                if (data["atom representatives"].is_array())
                    setRepresentativesFromJsonArray(
                        data["atom representatives"],
                        crystal,
                        subsystems,
                        representatives);
                else
                    on_error::throwException(
                        "invalid format of atom representatives definition in json data, expected array",
                        __FILE__, __LINE__);
                return;
            }

            string atomFile = data.value("qm atoms use", string());

            vector<vector<pair<string, string> > > subsystemAtoms;
            vector<string> subsystemLabels;
            for (auto& subsystem : subsystems)
            {
                subsystemAtoms.push_back(subsystem.atoms.atomList);
                subsystemLabels.push_back(subsystem.label);
            }


            if (atomFile.empty())
                gar_utilities::findDefaultRepresentatives(crystal, subsystemAtoms, representatives);
            else
                discamb_io::read_representatives(atomFile, crystal, subsystemLabels, subsystemAtoms, representatives);

        }

        void setSubsystemsWfnCalculation(
            const nlohmann::json& data,
            const std::vector<QmFragmentInCrystal>& subsystems,
            double clusterThreshold,
            std::vector<CrystalFragmentWfnCalculation> & subsystemWfnCalculation)
        {
            subsystemWfnCalculation.clear();
            subsystemWfnCalculation.resize(subsystems.size());
            map<string, int> label2idx;
            for (int i = 0; i < subsystems.size(); i++)
                label2idx[subsystems[i].label] = i;
            if (data.find("qm settings per subsystem") != data.end())
            {
                for (auto const& item : data.find("qm settings per subsystem").value().items())
                {
                    //subsystemWfnCalculation[label2idx[item.key()]]. = defaultMultipoleExpansionSettings;
                    subsystemWfnCalculation[label2idx[item.key()]].set(item.value(), clusterThreshold);
                }
            }
        }


    }

    void HirshfeldAtomModelSettings::set(
        const nlohmann::json& data,
        const Crystal& crystal)
    {
        ham_settings::setCrystalFragments(data, crystal, crystalFragments);
        electronDensityPartition.set(data);
        wfnCalculation.set(data);
        multipoleExpansion.set(data);
        ham_settings::setSubsystemsWfnCalculation(data, crystalFragments, multipoleExpansion.clusterThreshold, fragmentWfnCalculation);
        formFactorsCalculation.set(data);
        hardware.set(data);
        diagnostics.set(data);
        ham_settings::setRepresentatives(data, crystal, crystalFragments, representatives);
        sphericalHarmonicsExpansionLevel = data.value("spherical harmonics expansion level", sphericalHarmonicsExpansionLevel);
    }
    
}