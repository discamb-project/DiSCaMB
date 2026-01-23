#pragma once
#include "Program.h"

#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"
#include "discamb/AtomTyping/MolecularAtomTypeAssigner.h"
#include "discamb/IO/mol2_io.h"

#include <vector>
#include <filesystem>
#include <map>


class AssignAtomTypes: public Program
{
public:
    AssignAtomTypes();
    virtual ~AssignAtomTypes();
    virtual void set();
    virtual void run();
    static void readBank(
        const std::string &bankFile, 
        const std::string &selectedTypesFile, 
        std::vector< discamb::AtomType>& atomTypes,
        discamb::DescriptorsSettings & descriptorsSettings);
private:
    std::filesystem::path mChosenResFolder = std::filesystem::current_path() / std::string("chosen");
    std::filesystem::path mMolFolder = std::filesystem::current_path() / std::string("mol");
    std::string mOutputFileName = std::string("assignment.log");
    std::string mLcsCheckFileName = std::string("lcs_angle_check.log");
    double mLcsAbsCosAngleThreshold = 0.86602540378; // cos(30 degrees)
    std::vector< discamb::AtomType> mAtomTypes;
    discamb::DescriptorsSettings mDescriptorsSettings;
    //int mMinNumberOfInstances = 1;
    discamb::MolecularAtomTypeAssigner mAssigner;
    
    static void findFormulas(
        const discamb::mol2_io::Mol2Data& mol2Data,
        const std::vector<int>& atomicNumbers,
        std::vector<std::map<int, int> >& formulas);

    struct StructuralFormula{
        std::vector<int> z;
        std::vector < std::vector <int> > connectivity;
    };

    static void findSubstructures(
        const discamb::mol2_io::Mol2Data& mol2Data,
        const std::vector<int>& atomicNumbers,
        std::vector<std::vector<int> >& atomIndices,
        std::vector<StructuralFormula>& structuralFormulas);
        


    static int findStructuralFormula(const StructuralFormula& sf, const std::vector<StructuralFormula>& sfs);
};

