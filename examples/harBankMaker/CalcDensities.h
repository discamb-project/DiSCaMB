#pragma once
#include "Program.h"
#include "discamb/QuantumChemistry/ElectronDensityPartition.h"
#include "discamb/MathUtilities/Matrix3.h"

#include "json.hpp"

#include <memory>

class CalcDensities : public Program
{
public:
    CalcDensities();
    virtual ~CalcDensities();
    virtual void set();
    virtual void run();
private:
    
    bool mSaveAtomicDensityFiles = false;
    //static std::shared_ptr<ElectronDensityPartition> createElectronDensityPartition(
    //    const std::string &partitionName, 
    //    const std::string &wfxFileName, 
    //    const nlohmann::json& data);

    nlohmann::json mJsonSettings;
    std::string mBankFile, mSelectedAtomTypesFile;

    struct WfnSystem {
        std::vector<int> atomTypes;
        std::vector<int> atomicNumbers;
        std::vector<std::string> atomLabels;
        std::vector<discamb::Vector3d> positions;
        std::vector<discamb::Matrix3d> lcs;
        std::string structureName;
        int idxInStructure = 1;
        bool wfxExists = false;
    };

    std::vector<WfnSystem> mWfnSystems;

    //std::vector<std::vector<int> > mAtomIndices;
    //std::vector<std::vector<int> > mAtomTypes;
    //std::vector<std::vector<int> > mAtomicNumbers;
    //std::vector<std::vector<discamb::Matrix3d> > mLcs;
    int mRadialGridSize = 0;
    int mAngularGridSize = 0;
    bool mSkipDensityCalculationIfWfxAbsent = false;
    bool mCalcDipole = true;
    std::vector<std::string> mTypeNames;
    // atomic number
    std::map<int, std::vector<double> > mRadialGrids;
    std::vector<discamb::Vector3d> mAngularGrid;
    std::vector<std::vector<double> > mRadialWeights;
    std::vector<double> mAngularWeights;
    // [type][instance].first - WfnSystem index, .second - atom index
    std::vector<std::vector<std::pair<int, int> > > mTypeOccurences;
    
    void readStructures();
    void setGrids();

    static void printAtomicDensity(const std::string& fileName, const std::vector<double>& density, int nAngular, int nRadial, int atomicNumber);
};
