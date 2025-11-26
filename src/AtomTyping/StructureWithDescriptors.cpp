#include "discamb/AtomTyping/StructureWithDescriptors.h"

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/StructuralProperties/CovalentRadiousBondDetector.h"
#include "discamb/StructuralProperties/RingCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include <algorithm>
#include <cmath>


using namespace std;

namespace {

	int maxRingSize = 8;
}


namespace discamb{

/*
    double covalentBondThreshold = 0.1;
    double atomPlanarityThreshold = 0.1;
    double ringPlanarityThreshold = 0.1;
    double atomInRingPlanarityThreshold = 0.1;
    int atomInRingMaxNeighbourCount = 3;
*/

void DescriptorsSettings::readFromJson(
    const nlohmann::json& data)
{
    covalentBondThreshold = data.value("covalent bond threshold", covalentBondThreshold);
    atomPlanarityThreshold = data.value("atom planarity threshold", atomPlanarityThreshold);
    ringPlanarityThreshold = data.value("ring planarity threshold", ringPlanarityThreshold);
    atomInRingPlanarityThreshold = data.value("atom in ring planarity threshold", atomInRingPlanarityThreshold);
    atomInRingMaxNeighbourCount = data.value("atom in ring max neighbour count", atomInRingMaxNeighbourCount);
    if(data.find("max bond length in aromatic ring")!= data.end())
    {
        maxBondLengthAromaticRing.clear();
        for (auto& [key, value] : data["max bond length in aromatic ring"].items())
        {
            vector<string> elements;
            string_utilities::split(key, elements, ',');
            if (elements.size() != 2)
                throw std::runtime_error("Invalid key in max bond length aromatic ring settings: " + key);
            int z1 = periodic_table::atomicNumber(elements[0]);
            int z2 = periodic_table::atomicNumber(elements[1]);
            double length = value.get<double>();
            maxBondLengthAromaticRing[{z1, z2}] = length;
        }
    }
}

void StructureWithDescriptors::set(
    const MoleculeData& m)
{
    set(m.atomicNumbers, m.atomPositions, m.atomLabels);
}


void StructureWithDescriptors::set(
    const std::vector<int> &atomicNumbers, 
    const std::vector<Vector3d> &positions,
    const std::vector<std::string> &labels)
{
    CovalentRadiousBondDetector bondDetector;
    int nAtoms, i, j;
    vector<vector<int> > connectivity;

    bondDetector.setThreshold(settings.covalentBondThreshold);//angstroms

    nAtoms = atomicNumbers.size();

    connectivity.resize(nAtoms);

    for (i = 0; i < nAtoms; i++)
        for (j = 0; j < i; j++)
        if (bondDetector.areBonded(atomicNumbers[i], positions[i], atomicNumbers[j], positions[j]))
            {
                connectivity[i].push_back(j);
                connectivity[j].push_back(i);
            }
    
    set(atomicNumbers, positions, connectivity, labels);
}


void StructureWithDescriptors::set(
    const std::vector<int> &atomicNumbers,
    const std::vector<Vector3d> &positions,
    const std::vector<std::vector<int> > &predefinedConnectivity,
    const std::vector<std::string> &_labels)
{
    int nAtoms, atomIndex, ringIndex,nRings,ringSize,j,nNeighbors;
    RingCalculator ringCalculator;
    //PlanarityCalculator planarityCalculator;

    vector<string> labels = _labels;
    if (labels.empty())
        for (int z : atomicNumbers)
            labels.push_back(periodic_table::symbol(z));


    ringCalculator.set(
        settings.ringPlanarityThreshold,
        settings.atomInRingPlanarityThreshold, 
        settings.atomInRingMaxNeighbourCount,
        settings.maxBondLengthAromaticRing);



    connectivity = predefinedConnectivity;

    nAtoms = connectivity.size();

    

    
    atomDescriptors.clear();
    atomDescriptors.resize(nAtoms);

    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {

		atomDescriptors[atomIndex].position = positions[atomIndex];
        atomDescriptors[atomIndex].label = labels[atomIndex];
        atomDescriptors[atomIndex].atomicNumber = atomicNumbers[atomIndex];
        atomDescriptors[atomIndex].nNeighbors = connectivity[atomIndex].size();

        // get neighbors formula
        nNeighbors = connectivity[atomIndex].size();
        //atomDescriptors[atomIndex].neighborsFormula.resize(nNeighbors);
        //for (j = 0; j < nNeighbors; j++)
        //    atomDescriptors[atomIndex].neighborsFormula[j] = atomicNumbers[connectivity[atomIndex][j]];

        atomDescriptors[atomIndex].neighborsFormula.clear();
        for (j = 0; j < nNeighbors; j++)
            atomDescriptors[atomIndex].neighborsFormula.insert(atomicNumbers[connectivity[atomIndex][j]]);


        //sort(atomDescriptors[atomIndex].neighborsFormula.begin(),atomDescriptors[atomIndex].neighborsFormula.end());

        atomDescriptors[atomIndex].planar = structural_properties::atomPlanarity(atomIndex, atomicNumbers, positions, connectivity[atomIndex], atomDescriptors[atomIndex].planarityDisplacementEsd);
        atomDescriptors[atomIndex].planarRingsIndices.clear();
        atomDescriptors[atomIndex].planarRingsSizes.clear();
        //atomDescriptors[i].
    }

    vector<vector<int> > rings;
	vector<double> atomPlanarityEsd(nAtoms);

	//----
	

	for (int i = 0; i < nAtoms; i++)
		atomPlanarityEsd[i] = atomDescriptors[i].planarityDisplacementEsd;
	//----
  //  cout << __LINE__ << endl;
    ringCalculator.calculateRings(connectivity, positions, atomPlanarityEsd, maxRingSize, planarRings, planarRingsPlanarityEsd, atomicNumbers);
    //ringCalculator.calculateRings(connectivity, positions, atomPlanarityEsd, maxRingSize, planarRings, planarRingsPlanarityEsd);
    nRings = planarRings.size();
    vector<vector<pair<int, int> > > atomsRings(nAtoms);

//    cout << __LINE__ << endl;

    for (ringIndex = 0; ringIndex < nRings; ringIndex++)
    {
        ringSize = planarRings[ringIndex].size();
        for (int i = 0; i < ringSize; i++)
        {
            //atomDescriptors[planarRings[ringIndex][i]].ringsIndices.push_back(ringIndex);
            //atomDescriptors[planarRings[ringIndex][i]].ringsSizes.push_back(ringSize);
            atomIndex = planarRings[ringIndex][i];
            atomsRings[atomIndex].push_back({ ringSize , ringIndex });
        }
    }

    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {
        nRings = atomsRings[atomIndex].size();
        sort(atomsRings[atomIndex].begin(), atomsRings[atomIndex].end());
        for (ringIndex = 0; ringIndex < nRings; ringIndex++)
        {
            atomDescriptors[atomIndex].planarRingsIndices.push_back(atomsRings[atomIndex][ringIndex].second);
            atomDescriptors[atomIndex].planarRingsSizes.push_back(atomsRings[atomIndex][ringIndex].first);
        }
    }
    
    // setting threeMemberRings & fourMemberRings;
    RingCalculator smallRingCalculator;
    //smallRingCalculator.set(100, 100, 100);
    vector<vector<int> > smallRings;
    //vector<double> smallRingsPlanarityEsd;

    //smallRingCalculator.calculateRings(connectivity, positions, atomPlanarityEsd, 4, smallRings, smallRingsPlanarityEsd);

    
    std::set<std::set<int> > rings3, rings4;
    RingCalculator::calculate34Rings(connectivity, rings3, rings4);
    for (auto &ring : rings3)
        smallRings.push_back(vector<int>(ring.begin(), ring.end()));
    for (auto &ring : rings4)
        smallRings.push_back(vector<int>(ring.begin(), ring.end()));
    


    for (auto &atom : atomDescriptors)
    {
        atom.n3memberRings = 0;
        atom.n4memberRings = 0;
    }


    for (auto const &ring : smallRings)
        if (ring.size() == 3)
            for (auto atomIdx : ring)
                atomDescriptors[atomIdx].n3memberRings++;
        else
            for (auto atomIdx : ring)
                atomDescriptors[atomIdx].n4memberRings++;

    //-----------

    Vector3d interatomicVector;
    nAtoms = atomicNumbers.size();
    bondLengths.resize(nAtoms);


    for (int i = 0; i < nAtoms; i++)
    {
        nNeighbors = connectivity[i].size();
        bondLengths[i].resize(nNeighbors);

        for (j = 0; j < nNeighbors; j++)
        {
            interatomicVector = positions[i] - positions[connectivity[i][j]];
            bondLengths[i][j] = sqrt(interatomicVector*interatomicVector);
        }
    }

}

void StructureWithDescriptors::setFromCrystalFragment(
    const std::vector<int> &atomicNumbers,
    const std::vector<Vector3d> &positions,
    const DataForCrystalFragment &data,
    const std::vector<std::string> &labels)
{
    CovalentRadiousBondDetector bondDetector;
    int nAtoms, i, j;
    vector<vector<int> > connectivity;

    bondDetector.setThreshold(settings.covalentBondThreshold);//angstroms

    nAtoms = atomicNumbers.size();

    connectivity.resize(nAtoms);

    for (i = 0; i < nAtoms; i++)
        for (j = 0; j < i; j++)
            if (bondDetector.areBonded(atomicNumbers[i], positions[i], atomicNumbers[j], positions[j]))
            {
                connectivity[i].push_back(j);
                connectivity[j].push_back(i);
            }

    setFromCrystalFragment(atomicNumbers, positions, connectivity, data, labels);


}

void StructureWithDescriptors::setFromCrystalFragment(
    const std::vector<int> &atomicNumbers,
    const std::vector<Vector3d> &positions,
    const std::vector<std::vector<int> > &predefinedConnectivity,
    const DataForCrystalFragment &data,
    const  std::vector<std::string> &_labels)
{
    int nAtoms, atomIndex, ringIndex, nRings, ringSize, j, nNeighbors;
    RingCalculator ringCalculator;
    //PlanarityCalculator planarityCalculator;

    vector<string> labels = _labels;
    if (labels.empty())
        for (int z : atomicNumbers)
            labels.push_back(periodic_table::symbol(z));


    ringCalculator.set(settings.ringPlanarityThreshold, settings.atomInRingPlanarityThreshold, settings.atomInRingMaxNeighbourCount, settings.maxBondLengthAromaticRing);

    connectivity = predefinedConnectivity;

    nAtoms = connectivity.size();

    atomDescriptors.clear();
    atomDescriptors.resize(nAtoms);

    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {

        atomDescriptors[atomIndex].position = positions[atomIndex];
        atomDescriptors[atomIndex].label = labels[atomIndex];
        atomDescriptors[atomIndex].atomicNumber = atomicNumbers[atomIndex];
        atomDescriptors[atomIndex].nNeighbors = connectivity[atomIndex].size();

        // get neighbors formula
        nNeighbors = connectivity[atomIndex].size();
        //atomDescriptors[atomIndex].neighborsFormula.resize(nNeighbors);
        atomDescriptors[atomIndex].neighborsFormula.clear();// resize(nNeighbors);
        for (j = 0; j < nNeighbors; j++)
            atomDescriptors[atomIndex].neighborsFormula.insert(atomicNumbers[connectivity[atomIndex][j]]);
            //atomDescriptors[atomIndex].neighborsFormula[j] = atomicNumbers[connectivity[atomIndex][j]];
       // sort(atomDescriptors[atomIndex].neighborsFormula.begin(), atomDescriptors[atomIndex].neighborsFormula.end());

        atomDescriptors[atomIndex].planar = structural_properties::atomPlanarity(atomIndex, atomicNumbers, positions, connectivity[atomIndex], atomDescriptors[atomIndex].planarityDisplacementEsd);
        atomDescriptors[atomIndex].planarRingsIndices.clear();
        atomDescriptors[atomIndex].planarRingsSizes.clear();
        //atomDescriptors[i].
    }

    //---- planar rings

    vector<vector<int> > rings;
    vector<double> atomPlanarityEsd(nAtoms);
    vector<int> atomsToTakeIntoAccount;
    int nShellsToTakeIntoAccount = min(int(data.shells.size()), data.namedNeighboursRange+1);
    for (int shellIdx = 0; shellIdx < nShellsToTakeIntoAccount ; shellIdx++)
        atomsToTakeIntoAccount.insert(atomsToTakeIntoAccount.end(), data.shells[shellIdx].begin(), data.shells[shellIdx].end());

    for (int i = 0; i < nAtoms; i++)
        atomPlanarityEsd[i] = atomDescriptors[i].planarityDisplacementEsd;
    // ringCalculator.calculateRings(connectivity, positions, atomPlanarityEsd, maxRingSize, planarRings, planarRingsPlanarityEsd);
    ringCalculator.calculateRings(atomsToTakeIntoAccount, connectivity, positions, atomPlanarityEsd, maxRingSize, planarRings, planarRingsPlanarityEsd, atomicNumbers);
    nRings = planarRings.size();
    vector<vector<pair<int, int> > > atomsRings(nAtoms);



    for (ringIndex = 0; ringIndex < nRings; ringIndex++)
    {
        ringSize = planarRings[ringIndex].size();
        for (int i = 0; i < ringSize; i++)
        {
            atomIndex = planarRings[ringIndex][i];
            atomsRings[atomIndex].push_back({ ringSize , ringIndex });
        }
    }

    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {
        nRings = atomsRings[atomIndex].size();
        sort(atomsRings[atomIndex].begin(), atomsRings[atomIndex].end());
        for (ringIndex = 0; ringIndex < nRings; ringIndex++)
        {
            atomDescriptors[atomIndex].planarRingsIndices.push_back(atomsRings[atomIndex][ringIndex].second);
            atomDescriptors[atomIndex].planarRingsSizes.push_back(atomsRings[atomIndex][ringIndex].first);
        }
    }

    // setting threeMemberRings & fourMemberRings;

    //SimpleRingCalculator smallRingCalculator;
    //smallRingCalculator.set(100, 100, 100);
    vector<vector<int> > smallRings;
    //vector<double> smallRingsPlanarityEsd;

    atomsToTakeIntoAccount = data.shells[0];
    
    std::set<std::set<int> > rings3, rings4;
    //smallRingCalculator.calculate34Rings(connectivity, rings3, rings4);
    RingCalculator::calculate34Rings(atomsToTakeIntoAccount, connectivity, rings3, rings4);
    //smallRingCalculator.calculateRings(connectivity, positions, atomPlanarityEsd, 4, planarRings, planarRingsPlanarityEsd);
    
    for (auto &ring : rings3)
        smallRings.push_back(vector<int>(ring.begin(), ring.end()));
    for (auto &ring : rings4)
        smallRings.push_back(vector<int>(ring.begin(), ring.end()));
        
    //    cout << __LINE__ << endl;

    for (auto &atom : atomDescriptors)
    {
        atom.n3memberRings = 0;
        atom.n4memberRings = 0;
        //atom.threeMemberRingsIndices.clear();
        //atom.fourMemberRingsIndices.clear();
    }

    //fourMemberRings.clear();
    //threeMemberRings.clear();

    //int threeMemberRingIndex = 0;
    //int fourMemberRingIndex = 0;

    for (auto const &ring : smallRings)
        if (ring.size() == 3)
        {
            //threeMemberRings.push_back(ring);
            for (auto atomIdx : ring)
                atomDescriptors[atomIdx].n3memberRings++;
                //atomDescriptors[atomIdx].threeMemberRingsIndices.push_back(threeMemberRingIndex);
            //threeMemberRingIndex++;
        }
        else
        {
            //fourMemberRings.push_back(ring);
            for (auto atomIdx : ring)
                atomDescriptors[atomIdx].n4memberRings++;
                //atomDescriptors[atomIdx].fourMemberRingsIndices.push_back(fourMemberRingIndex);
            //fourMemberRingIndex++;
        }

    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    {
        int representative = data.mapToCoreShellAtoms[atomIndex];

        atomDescriptors[atomIndex].n3memberRings = atomDescriptors[representative].n3memberRings;
        atomDescriptors[atomIndex].n4memberRings = atomDescriptors[representative].n4memberRings;
    }


    // bond lengths

    Vector3d interatomicVector;
    nAtoms = atomicNumbers.size();
    bondLengths.resize(nAtoms);


    for (int i = 0; i < nAtoms; i++)
    {
        nNeighbors = connectivity[i].size();
        bondLengths[i].resize(nNeighbors);

        for (j = 0; j < nNeighbors; j++)
        {
            interatomicVector = positions[i] - positions[connectivity[i][j]];
            bondLengths[i][j] = sqrt(interatomicVector*interatomicVector);
        }
    }

}


} //namespace discamb

