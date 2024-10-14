
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"

#include "discamb/IO/xyz_io.h"

#include <iomanip>
#include <map>
#include <set>

using namespace std;

namespace discamb {
    namespace {
        bool isSphericalType(const string& s)
        {
            int nZeros=0;
            int nDigits = 0;
            for (auto c : s)
            {
                if (c == '0')
                    nZeros++;
                if (isdigit(c))
                    nDigits++;
            }
            if (nDigits == nZeros && nDigits == 3)
                return true;
            return false;
        }
    }

    CrystalAtomTypeAssigner::CrystalAtomTypeAssigner()
    {
        mDescriptorsSettings = DescriptorsSettings();
        mPlanarRingsRange = mRings34range = mTotalRange = mLabeledNeighbourRange = 0;
    }

    CrystalAtomTypeAssigner::~CrystalAtomTypeAssigner()
    {
    }

    std::string CrystalAtomTypeAssigner::typeLabel(
        int typeId) 
        const
    {
        if (typeId >= 0)
            return mAtomTypes[typeId].id;

        return "none";
    }

    void CrystalAtomTypeAssigner::assign(
        const Crystal& crystal,
        std::vector<int>& typeID,
        std::vector<std::vector<Vector3d> >& lcs)
        const
    {
        vector<LocalCoordinateSystem<AtomInCrystalID> > _lcs;
        assign(crystal, typeID, _lcs);
        int atomIdx, nAtoms = typeID.size();
        vector<Vector3d> lcs0 = { Vector3d(1,0,0), Vector3d(0,1,0) ,Vector3d(0,0,1) };
        lcs.assign(nAtoms, lcs0);

        LocalCoordinateSystemCalculator lcsCalculator;

        for(atomIdx=0; atomIdx<nAtoms; atomIdx++)
            if (typeID[atomIdx] >= 0)
            {
                lcsCalculator.set(_lcs[atomIdx], crystal);
                lcsCalculator.calculate(lcs[atomIdx][0], lcs[atomIdx][1], lcs[atomIdx][2], crystal);
            }
    }

    void CrystalAtomTypeAssigner::printAssignment(
        std::ostream &out,
        const Crystal &crystal,
        const std::vector<int> &typeID,
        const std::vector< LocalCoordinateSystem<AtomInCrystalID> > &lcs)
        const
    {
        int atomIdx, nAtoms = crystal.atoms.size();
        int maxAtomLabelWidth = 0;
        int maxTypeIdWidth = 0;
        vector<string> atomLabels;
        int nAssigned = 0;
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            atomLabels.push_back(crystal.atoms[atomIdx].label);
            if (atomLabels[atomIdx].size() > maxAtomLabelWidth)
                maxAtomLabelWidth = atomLabels[atomIdx].size();
            string typeLabels;
            if (typeID[atomIdx] >= 0)
            {
                typeLabels = mTypeLabels[typeID[atomIdx]];
                if(!isSphericalType(typeLabels))
                    nAssigned++;
            }
            else
                typeLabels = "undefined";
            if (maxTypeIdWidth < typeLabels.size())
                maxTypeIdWidth = typeLabels.size();
        }

        out << "Atom type assigned to " << nAssigned << " of " << nAtoms << ".\n";
        if (nAssigned != nAtoms)
        {
            out << "Atoms with unassigned atom types :\n";
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (typeID[atomIdx] < 0)
                    out << setw(maxAtomLabelWidth + 3) << crystal.atoms[atomIdx].label << "\n";
                else
                    if (isSphericalType(mTypeLabels[typeID[atomIdx]]))
                        out << setw(maxAtomLabelWidth + 3) << crystal.atoms[atomIdx].label << "\n";
        }
        if (nAssigned > 0)
        {
            out << "Atoms with assigned atom types and local coordinate systems:\n";
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                
                if (typeID[atomIdx] >= 0)
                    if (!isSphericalType(mTypeLabels[typeID[atomIdx]]))
                        out << setw(maxAtomLabelWidth + 3) << crystal.atoms[atomIdx].label
                            << setw(maxTypeIdWidth + 3) << mTypeLabels[typeID[atomIdx]]
                            << "    " << ubdbLcsAsString(lcs[atomIdx], atomLabels) << "\n";
            }
        }
        out << "\n";
    }

    void CrystalAtomTypeAssigner::setRanges()
    {
        mTotalRange = atom_typing_utilities::atomTypesRange(mAtomTypes, mDescriptorsSettings.maxPlanarRing, mLabeledNeighbourRange, mPlanarRingsRange, mRings34range);
    }

    void CrystalAtomTypeAssigner::setAtomTypes(
        const std::vector<AtomType> &atomTypes)
    {
        mAtomTypes = atomTypes;
        mTypeLabels.clear();
        for (auto const &type : atomTypes)
            mTypeLabels.push_back(type.id);
        mAssigner.setAtomTypes(atomTypes);
        setRanges();
    }

    void CrystalAtomTypeAssigner::setDescriptorsSettings(
        const DescriptorsSettings &settings)
    {
        mAssigner.setDescriptorsSettings(settings);
        mDescriptorsSettings = settings;
        setRanges();
    }

    void CrystalAtomTypeAssigner::assign(
        const Crystal& crystal,
        std::vector<int> &typeID,
        std::vector<LocalCoordinateSystem<AtomInCrystalID> > &lcs)
        const
    {
        UnitCellContent unitCellContent;
        vector<UnitCellContent::AtomID> asymmetricUnit, graph;
        unitCellContent.set(crystal);
        SpaceGroupOperation spaceGroupOperation;
        int atomIdx, nAtoms = crystal.atoms.size();
        Vector3<CrystallographicRational> translation;
        Matrix3i rotation;
        Vector3i latticeTranslation;

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            /*
            UnitCellContent atoms fractional coordinates are adjusted to be between 0 and 1,
            the atom pushed back into asymmetricUnit corresponds to original atom from asymmetric unit
            (with no cordinates adjusted)
            */
            spaceGroupOperation = unitCellContent.getGeneratingOperation(atomIdx, 0);
            spaceGroupOperation.getTranslation(translation);
            latticeTranslation.set(-translation[0].numerator(), -translation[1].numerator(), -translation[2].numerator());
            asymmetricUnit.push_back(UnitCellContent::AtomID(atomIdx, latticeTranslation));
        }

        //structural_properties::graphToNthNeighbour(unitCellContent, asymmetricUnit, graph, 8, mDescriptorsSettings.covalentBondThreshold);
        structural_properties::graphToNthNeighbour(unitCellContent, asymmetricUnit, graph, mTotalRange, mDescriptorsSettings.covalentBondThreshold);

        // graph as list of atoms and generating symmetry operations
        // atom order is preserved

        //vector<pair<int, SpaceGroupOperation> > graphAtoms;
        vector<AtomInCrystalID> graphAtoms;
        vector<int> atomicNumbers, atomicNumbersASU;
        vector<Vector3d> positions;
        Vector3d fractional, cartesian;
        int atomIndexInAsymmetricUnit;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbersASU);

        for (auto &atom : graph)
        {
            // atom symmetry operation
            unitCellContent.getGeneratingOperation(atom.atomIndex, 0).get(rotation, translation);
            translation += atom.unitCellPosition;
            spaceGroupOperation.set(rotation, translation);
            // --
            atomIndexInAsymmetricUnit = unitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex);
            graphAtoms.push_back({atomIndexInAsymmetricUnit, spaceGroupOperation});
            // generate coordinates
            spaceGroupOperation.apply(crystal.atoms[atomIndexInAsymmetricUnit].coordinates, fractional);
            crystal.unitCell.fractionalToCartesian(fractional, cartesian);
            positions.push_back(cartesian);
            // --
            atomicNumbers.push_back(atomicNumbersASU[atomIndexInAsymmetricUnit]);
        }

        vector<LocalCoordinateSystem<int> > lcsMolecule;
        StructureWithDescriptors structureWithDescriptors;
        //vector<int> typeId;
        vector<int> atomToAssign(crystal.atoms.size());

        for (int i = 0; i < crystal.atoms.size(); i++)
            atomToAssign[i] = i;

        structureWithDescriptors.set(atomicNumbers, positions);
        mAssigner.assign(structureWithDescriptors, atomToAssign, typeID, lcsMolecule);

        // convert lcs molecule to lcs crystal
        int lcsIdx, nLcs = lcsMolecule.size();
        lcs.clear();
        lcs.resize(nLcs);
        for(lcsIdx=0;lcsIdx<nLcs;lcsIdx++)
            convertUbdbLcs(lcsMolecule[lcsIdx], graphAtoms, lcs[lcsIdx]);
    }

    void CrystalAtomTypeAssigner::assign(
        const Crystal& crystal,
        std::vector<int> &typeID,
        std::vector<LocalCoordinateSystem<AtomInCrystalID> > &lcs,
        StructureWithDescriptors &structureWithDescriptors)
        const
    {
        vector<int> shellSizes;
        vector<pair<int, string> >  asuWithNeighbours;
        vector<int> atomicNumbers;
        vector<string> labels;
        vector<Vector3d> positions;

        structural_properties::assymetricUnitWithNeighbours(crystal, asuWithNeighbours, atomicNumbers, positions, labels, mTotalRange, mDescriptorsSettings.covalentBondThreshold, shellSizes);

        StructureWithDescriptors::DataForCrystalFragment data;
        data.namedNeighboursRange = mLabeledNeighbourRange;
        data.shells.resize(shellSizes.size());
        int counter = 0;
        for (int shellIdx = 0; shellIdx < shellSizes.size(); shellIdx++)
            for (int i = 0; i < shellSizes[shellIdx]; i++)
                data.shells[shellIdx].push_back(counter++);

        for (auto &atom : asuWithNeighbours)
            data.mapToCoreShellAtoms.push_back(atom.first);

        structureWithDescriptors.setFromCrystalFragment(atomicNumbers, positions, data, labels);

        vector<LocalCoordinateSystem<int> > lcsMolecule;
        vector<int> atomToAssign(crystal.atoms.size());

        for (int i = 0; i < crystal.atoms.size(); i++)
            atomToAssign[i] = i;
                
        mAssigner.assign(structureWithDescriptors, atomToAssign, typeID, lcsMolecule);

        vector<AtomInCrystalID> graphAtoms;

        for (auto &atom : asuWithNeighbours)
            graphAtoms.push_back(AtomInCrystalID(atom.first, SpaceGroupOperation(atom.second)));

        // convert lcs molecule to lcs crystal
        int lcsIdx, nLcs = lcsMolecule.size();
        lcs.clear();
        lcs.resize(nLcs);
        for (lcsIdx = 0; lcsIdx < nLcs; lcsIdx++)
            convertUbdbLcs(lcsMolecule[lcsIdx], graphAtoms, lcs[lcsIdx]);


    }


    void CrystalAtomTypeAssigner::assign(
        const Crystal& crystal,
        const std::vector<std::vector<AtomInCrystalID> >& fragments,
        const std::vector< std::vector<AtomRepresentativeInfo> >& atomRepresentatives,
        std::vector< std::vector<int> >& typeId,
        std::vector< std::vector<LocalCoordinateSystem<AtomInCrystalID> > >& lcs)
        const
    {
        int atomIdx, fragmentIdx, nFragments = fragments.size();

        typeId.clear();
        lcs.clear();
        typeId.resize(nFragments);
        lcs.resize(nFragments);

        map<string, int> label2atomIdx;
        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            label2atomIdx[crystal.atoms[atomIdx].label] = atomIdx;
        
        vector<set<int> > atom2Assign(nFragments);
        for (auto const& representatives : atomRepresentatives)
            for (auto const& representative : representatives)
                atom2Assign[representative.fragmentIdx].insert(representative.idxInSubsystem);


        vector<int> atomicNumbersAsymmetricUnit, atomicNumbersFragment;
        vector<Vector3d> positions;


        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbersAsymmetricUnit);

        for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
        {
            
            int nAtoms = fragments[fragmentIdx].size();
         
            positions.clear();
            atomicNumbersFragment.clear();

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                const SpaceGroupOperation &spaceGroupOperation = fragments[fragmentIdx][atomIdx].getSymmetryOperation();
                int atomIdxInCrystal = fragments[fragmentIdx][atomIdx].index();
                Vector3d r_cart, r_frac, r0_frac = crystal.atoms[atomIdxInCrystal].coordinates;
                spaceGroupOperation.apply(r0_frac, r_frac);
                crystal.unitCell.fractionalToCartesian(r_frac, r_cart);
                positions.push_back(r_cart);
                atomicNumbersFragment.push_back(atomicNumbersAsymmetricUnit[atomIdxInCrystal]);
            }

            //xyz_io::writeXyz(string("fragment_") + to_string(fragmentIdx + 1), atomicNumbersFragment, positions);

            vector<LocalCoordinateSystem<int> > lcsMolecule;
            StructureWithDescriptors structureWithDescriptors;
            vector<int> atomToAssign(atom2Assign[fragmentIdx].begin(), atom2Assign[fragmentIdx].end());
            
            structureWithDescriptors.set(atomicNumbersFragment, positions);
            mAssigner.assign(structureWithDescriptors, atomToAssign, typeId[fragmentIdx], lcsMolecule);

            // convert lcs molecule to lcs crystal
            int lcsIdx, nLcs = lcsMolecule.size();
            lcs[fragmentIdx].resize(nLcs);
            for (lcsIdx = 0; lcsIdx < nLcs; lcsIdx++)
                convertUbdbLcs(lcsMolecule[lcsIdx], fragments[fragmentIdx], lcs[fragmentIdx][lcsIdx]);

        }
    }

}


