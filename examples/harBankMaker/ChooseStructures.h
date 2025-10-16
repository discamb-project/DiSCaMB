#pragma once
#include "Program.h"

#include <filesystem>
#include <string>
#include <vector>
#include <set>


class ChooseStructures : public Program
{
public:
    ChooseStructures();
    virtual ~ChooseStructures();
    virtual void set();
    virtual void run();
    // reads chosen_mol file, mol_idx starts from 1
    static void readChosenMol(std::vector<std::string>& structureName, std::vector<int> &mol_idx);
private:
    std::filesystem::path mMolFolder = std::filesystem::current_path() / std::string("mol");
    std::string mAtomTypesLog = std::string("assignment.log");
    void readAtomTypeAssignement();
    std::vector<std::string> mParentStructureNames;
    std::vector<std::string> mStructuresToExclude;
    struct Substructure {
        int parentStructure;
        std::vector<std::string> atomLabels;
        std::vector<int> z;
        std::vector<std::string> types;
        int idxInMol2;
        std::map<int, int> formula;
        int formulaIdx;
        std::string formulaString;
    };
    std::vector<Substructure> mSubstructures;
    struct TypeStatistics {
        std::string typeLabel="----";
        int nOccurences = 0;
        int nContainingStructures = 0;
        int nOccurencesInMoleculesWithDiffrentFormula = 0;
    };
    std::vector<TypeStatistics> mTypeStatistics;
    int mTargetNumberOfContainingMolecules=1;
    bool mAllowLessThanTarget = false;
    bool search_further(const std::vector<int>& nContainingStructures) const;
    
    double calcScore(
        int structureIdx, 
        const std::vector<std::vector<int> >& containingMolecules, 
        const std::set<int> &typesInStructure,
        const std::vector<bool>& typeInstancesFound,
        const std::vector<bool>& formulaUsed) const;

    static double formulaSimilarity(const std::map<int, int>& formula1, const std::map<int, int>& formula2);
    
};

