#include "CalcDensities.h"
#include "ChooseStructures.h"
#include "AssignAtomTypes.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/IO/vtk_io.h"
#include "discamb/IO/wfn_io.h"
#include "discamb/BasicUtilities/Constants.h"
#include "discamb/AtomTyping/MolecularLcsCalculator.h"
#include "discamb/MathUtilities/radial_grid.h"
#include "discamb/MathUtilities/lebedev_laikov.h"
#include "discamb/QuantumChemistry/HirshfeldPartition.h"

#include <fstream>

using namespace discamb;
using namespace std;


CalcDensities::CalcDensities()
{
}

CalcDensities::~CalcDensities() 
{
}

void CalcDensities::set()
{
    readSettings(mJsonSettings);
    mRadialGridSize = mJsonSettings.value("radial grid", mRadialGridSize);
    mAngularGridSize = mJsonSettings.value("angular grid", mAngularGridSize);

    if (mRadialGridSize == 0 || mAngularGridSize == 0)
        on_error::throwException("define 'radial grid' and 'angular grid', at least one of them is missing", __FILE__, __LINE__);

    mBankFile = mJsonSettings.value("bank file", string(""));
    if (mBankFile.empty())
        on_error::throwException("expected \"bank file\" field in setting.json", __FILE__, __LINE__);

    mSelectedAtomTypesFile = mJsonSettings.value("selected atom type ids file", string(""));

    mSkipDensityCalculationIfWfxAbsent = mJsonSettings.value("skip density calculation if wfx absent", mSkipDensityCalculationIfWfxAbsent);

    mSaveAtomicDensityFiles = mJsonSettings.value("save atomic densities", mSaveAtomicDensityFiles);
}

void CalcDensities::readStructures()
{
    vector<string> structureName;
    vector<int> idxInStructure;

    ChooseStructures::readChosenMol(structureName, idxInStructure);

    int substructureIdx, nSubstructures = idxInStructure.size();

    mWfnSystems.resize(nSubstructures);

    vector<AtomType> atomTypes;
    DescriptorsSettings descriptorsSettings;
    AssignAtomTypes::readBank(mBankFile, mSelectedAtomTypesFile, atomTypes, descriptorsSettings);
    int nTypes = atomTypes.size();
    mTypeNames.resize(nTypes);
    for (int i = 0; i < nTypes; i++)
        mTypeNames[i] = atomTypes[i].id;
    mTypeOccurences.clear();
    mTypeOccurences.resize(nTypes);
    MolecularAtomTypeAssigner assigner;
    assigner.setAtomTypes(atomTypes);
    assigner.setDescriptorsSettings(descriptorsSettings);

    
    for (substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        
        filesystem::path p = "./chosen_mol_dir";
        p /= structureName[substructureIdx] + ".mol2";
        mol2_io::Mol2Data data;
        mol2_io::read(p.string(), data);
        vector<mol2_io::Mol2Data> substructures;
        data.split(substructures);

        mol2_io::Mol2Data& substructure = substructures[idxInStructure[substructureIdx]-1];
        
        StructureWithDescriptors structure;
        vector<int> atomicNumbers;
        mol2_io::atomicNumbers(substructure, atomicNumbers);

        int nAtoms = atomicNumbers.size();

        mWfnSystems[substructureIdx].lcs.resize(nAtoms);
        mWfnSystems[substructureIdx].atomTypes.resize(nAtoms);
        mWfnSystems[substructureIdx].structureName = structureName[substructureIdx];
        mWfnSystems[substructureIdx].atomLabels = substructure.atomName;
        mWfnSystems[substructureIdx].idxInStructure = idxInStructure[substructureIdx];
        mWfnSystems[substructureIdx].atomicNumbers = atomicNumbers;
        mWfnSystems[substructureIdx].positions = substructure.atomPosition;
        //mWfnSystems[substructureIdx].a

        string wfxFileName = structureName[substructureIdx] + "_" + to_string(idxInStructure[substructureIdx]) + ".wfx";
        string filePath = (filesystem::current_path() / string("wfnCalc") / wfxFileName).string();
        mWfnSystems[substructureIdx].wfxExists = filesystem::exists(filePath);
        if (mWfnSystems[substructureIdx].wfxExists)
        {
            wfn_io::WfnFileData wfnData;
            wfn_io::read_wfx(filePath, wfnData);
            mWfnSystems[substructureIdx].positions = wfnData.center_position;
            for (auto& r : mWfnSystems[substructureIdx].positions)
                r /= constants::Angstrom;
        }
        else
            continue;

        structure.set(atomicNumbers, substructure.atomPosition, substructure.atomName);
        
        vector<LocalCoordinateSystem<int> > lcs;
        assigner.assign(structure, mWfnSystems[substructureIdx].atomTypes, lcs);
        const vector<int> &typeId = mWfnSystems[substructureIdx].atomTypes;
        MolecularLcsCalculator lcsCalculator;
        
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            if (typeId[atomIdx] >= 0)
            {
                bool sameChirality = true;
                lcsCalculator.set(lcs[atomIdx]);
                lcsCalculator.calculate(mWfnSystems[substructureIdx].lcs[atomIdx], substructure.atomPosition, sameChirality);
                mTypeOccurences[typeId[atomIdx]].push_back({ substructureIdx , atomIdx});
            }
        }

    }
}

void CalcDensities::printAtomicDensity(
    const std::string& fileName,
    const std::vector<double>& density,
    int nAngular,
    int nRadial,
    int atomicNumber)
{
    ofstream out(fileName);
    if (!out.good())
        on_error::throwException("cannot write to file " + fileName, __FILE__, __LINE__);

    out << "atomic number " << atomicNumber << "\n"
        << "angular grid " << nAngular << "\n"
        << "radial grid " << nRadial << "\n";
    out << setprecision(12);
    for (double v : density)
        out << v << "\n";
    out.close();
}

void CalcDensities::setGrids()
{
    std::set<int> unique_atomic_numbers;
    for (auto& system : mWfnSystems)
        unique_atomic_numbers.insert(system.atomicNumbers.begin(), system.atomicNumbers.end());

    vector<double> radialGrid, weights;
//    mRadialGridSize,  mAngularGridSize 
    mRadialWeights.resize(120);
    for (int z : unique_atomic_numbers)
    {
        radial_grid::mura_knowles(z, mRadialGridSize, radialGrid, weights);
        mRadialGrids[z] = radialGrid;
        mRadialWeights[z] = weights;
    }
    lebedev_laikov::get_grid(mAngularGridSize, mAngularGrid, mAngularWeights);
}

void CalcDensities::run()
{
    readStructures();

    setGrids();

    if (mRadialGrids.empty())
        return;

    // reads atom types

    vector<string> wfxFiles;
    
    int substructureIdx, nSubstructures = mWfnSystems.size();

    string partitionName = mJsonSettings.value("partition", string("Hirshfeld"));

    int nTypes = mTypeNames.size();
    int radialIdx, angularIdx;
    int radialGridSize = mRadialGrids.begin()->second.size();
    int angularGridSize = mAngularGrid.size();
    int gridSize = radialGridSize * angularGridSize;
    vector<int> nTypeOccurences(nTypes, 0);
    vector<vector<double> > typeElectronDensitites(nTypes, vector<double>(gridSize, 0.0));
    vector<double> atomicElectronDensity(gridSize), diffDensity(gridSize);
    
    bool calcDiffDensity = (string_utilities::toLower(partitionName) == string("hirshfeld"));

    ofstream out;
    if (mCalcDipole)
        out.open("atomic_multipoles");

    vector<vector<double> > atomicCharges(nSubstructures);
    vector<vector<Vector3d> > atomicDipoles(nSubstructures);

    for (substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        WfnSystem& system = mWfnSystems[substructureIdx];

        if (mCalcDipole)
            out << system.structureName << " " << system.idxInStructure << "\n";

        shared_ptr<ElectronDensityPartition> partition;
        shared_ptr<HirshfeldPartition> hPartition;
        shared_ptr<ElectronDensityCalculator> edCalculator = shared_ptr<ElectronDensityCalculator>(new ElectronDensityCalculator);
        string wfxFileName = system.structureName + "_" + to_string(system.idxInStructure) + ".wfx";
        string filePath = (filesystem::current_path() / string("wfnCalc") / wfxFileName).string();

        if (!filesystem::exists(filePath))
        {
            if (mSkipDensityCalculationIfWfxAbsent)
                continue;
            else
                on_error::throwException("missing file '" + filePath + "'", __FILE__, __LINE__);
        }

        

        cout << "processing " << wfxFileName << "\n";
        edCalculator->setFromWavefunctionFile(filePath);

        partition = shared_ptr<ElectronDensityPartition>(ElectronDensityPartition::create(partitionName));
        if(mJsonSettings.contains("density calculation range"))
            partition->setIncludeRange(mJsonSettings["density calculation range"].get<double>());
        vector<Vector3d> positions = system.positions;
        for (auto& r : positions)
            r *= constants::Angstrom;
        //partition->set(mJsonSettings, system.atomicNumbers, system.positions, edCalculator);
        partition->set(mJsonSettings, system.atomicNumbers, positions, edCalculator);
        hPartition = static_pointer_cast<HirshfeldPartition>(partition);

        int atomIdx, nAtoms = system.atomicNumbers.size();
        
        atomicCharges[substructureIdx].resize(nAtoms);
        atomicDipoles[substructureIdx].resize(nAtoms);
        

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if(system.atomTypes[atomIdx]>=0)
            {
                nTypeOccurences[system.atomTypes[atomIdx]]++;
                Matrix3d& lcs = system.lcs[atomIdx];
                int z = system.atomicNumbers[atomIdx];
                int typeIdx = system.atomTypes[atomIdx];
                
                vector<double>& radialGrid = mRadialGrids[z];
                vector<Vector3d> gridPointPositions;
                int gridIdx = 0;
                for (angularIdx = 0; angularIdx < angularGridSize; angularIdx++)
                {
                    partition->setAtom(atomIdx);
                    Vector3d r = lcs * mAngularGrid[angularIdx];
                    for (radialIdx = 0; radialIdx < radialGridSize; radialIdx++)
                    {
                        Vector3d gridPoint = r * radialGrid[radialIdx] + system.positions[atomIdx] * constants::Angstrom;
                        double v = partition->calculate(gridPoint);
                        typeElectronDensitites[typeIdx][gridIdx] += v;
                        atomicElectronDensity[gridIdx] = v;
                        if (calcDiffDensity)
                            diffDensity[gridIdx] = hPartition->calculateDifferenceDensity(atomIdx, gridPoint);
                        gridPointPositions.push_back(gridPoint/constants::Angstrom);
                        gridIdx++;
                    }
                }
                if (mSaveAtomicDensityFiles)
                {
                    string fileNameCore = system.structureName + "_" + to_string(system.idxInStructure) + "_" + system.atomLabels[atomIdx];
                    string fileName = fileNameCore + ".den";
                    printAtomicDensity(fileName, atomicElectronDensity, angularGridSize, radialGridSize, z);
                    for (double& d : atomicElectronDensity)
                        d *= constants::Angstrom * constants::Angstrom * constants::Angstrom;
                    vtk_io::save_point_data(fileNameCore + ".vtk", gridPointPositions, atomicElectronDensity, fileNameCore);
                    vtk_io::save_point_data(fileNameCore + "_diff.vtk", gridPointPositions, diffDensity, fileNameCore);
                    for (double& d : atomicElectronDensity)
                        d /= constants::Angstrom * constants::Angstrom * constants::Angstrom;
                }
                if (mCalcDipole)
                {
                    int idx = 0;
                    Vector3d dipole;
                    double charge = z;
                    
                    for (angularIdx = 0; angularIdx < angularGridSize; angularIdx++)
                    {
                        Vector3d r = lcs * mAngularGrid[angularIdx];
                        for (radialIdx = 0; radialIdx < radialGridSize; radialIdx++)
                        {
                            Vector3d r = radialGrid[radialIdx] * mAngularGrid[angularIdx];
                            double w = 4.0 * M_PI * r * r * mRadialWeights[z][radialIdx] * mAngularWeights[angularIdx];
                            charge -= atomicElectronDensity[idx] * w;
                            dipole -= r * atomicElectronDensity[idx] * w;
                            idx++;
                        }
                    }
                    if (mCalcDipole)
                    {
                        out << system.atomLabels[atomIdx] << " " << charge << "\n"
                            << dipole.x << " " << dipole.y << " " << dipole.z << "\n";
                        atomicCharges[substructureIdx][atomIdx] = charge;
                        atomicDipoles[substructureIdx][atomIdx] = dipole;
                        dipole *= 1 / sqrt(dipole * dipole);
                        out<< dipole.x << " " << dipole.y << " " << dipole.z << "\n\n";
                    }

                }
            }

    }



    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        
        if (nTypeOccurences[typeIdx] > 0)
        {
            for (double& v : typeElectronDensitites[typeIdx])
                v /= nTypeOccurences[typeIdx];
            auto typeInstance = mTypeOccurences[typeIdx][0];
            int z = mWfnSystems[typeInstance.first].atomicNumbers[typeInstance.second];
            printAtomicDensity(string("type_") + mTypeNames[typeIdx], typeElectronDensitites[typeIdx], angularGridSize, radialGridSize, z);
        }
    }


    if (mCalcDipole)
    {
        out << setprecision(6) << fixed;

        for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {
            
            if (nTypeOccurences[typeIdx] > 0)
            {
                out << mTypeNames[typeIdx] << "\n";
                // charge
                out << "\ncharges\n\n";
                for (auto item : mTypeOccurences[typeIdx])
                    if (!atomicCharges[item.first].empty())
                    {
                        out << setw(12) << atomicCharges[item.first][item.second]
                            << "   " << mWfnSystems[item.first].structureName << "   " << mWfnSystems[item.first].idxInStructure
                            << "   " << mWfnSystems[item.first].atomLabels[item.second] << "\n";
                    }
                // dipole
                out << "\ndipoles\n\n";
                for (auto item : mTypeOccurences[typeIdx])
                    if (!atomicCharges[item.first].empty())
                    {
                        out << setw(12) << atomicDipoles[item.first][item.second].x << " "
                            << setw(12) << atomicDipoles[item.first][item.second].y << " "
                            << setw(12) << atomicDipoles[item.first][item.second].z << " "
                            << "   " << mWfnSystems[item.first].structureName << "   " << mWfnSystems[item.first].idxInStructure
                            << "   " << mWfnSystems[item.first].atomLabels[item.second] << "\n";
                    }
                // normalized dipole
                out << "\nnormalized dipoles\n\n";
                for (auto item : mTypeOccurences[typeIdx])
                    if (!atomicCharges[item.first].empty())
                    {
                        Vector3d d = atomicDipoles[item.first][item.second];
                        d /= sqrt(d * d);
                        out << setw(12) << d.x << " "
                            << setw(12) << d.y << " "
                            << setw(12) << d.z << " "
                            << "   " << mWfnSystems[item.first].structureName << "   " << mWfnSystems[item.first].idxInStructure
                            << "   " << mWfnSystems[item.first].atomLabels[item.second] << "\n";
                    }
                out << "\n";
            }
        }
    }

    if (mCalcDipole)
        out.close();
}

