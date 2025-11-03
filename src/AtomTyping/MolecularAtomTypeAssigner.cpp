#include "discamb/AtomTyping/MolecularAtomTypeAssigner.h"
#include "discamb/AtomTyping/TypeMatchAlgorithm.h"

using namespace std;

namespace discamb {

    MolecularAtomTypeAssigner::MolecularAtomTypeAssigner()
    {
    }

    MolecularAtomTypeAssigner::~MolecularAtomTypeAssigner()
    {
    }

    void MolecularAtomTypeAssigner::setDescriptorsSettings(
        const DescriptorsSettings &settings)
    {
        mDescriptorsSettings = settings;
    }

    void MolecularAtomTypeAssigner::setAtomTypes(
        const std::vector<AtomType> &atomTypes)
    {
        int nTypes = atomTypes.size();
        mAtomTypeMatchalgorithms.resize(nTypes);
        for (int i = 0; i < nTypes; i++)
            mAtomTypeMatchalgorithms[i].setType(atomTypes[i]);
    }

    void MolecularAtomTypeAssigner::assign(
        const std::vector<int> &atomicNumbers,
        const std::vector<Vector3d> &positions,
        const std::vector<std::string> &labels,
        std::vector<int> &atomsToAssign,
        std::vector<int> &typeID,
        std::vector<LocalCoordinateSystem<int> > &lcs)
        const
    {
        StructureWithDescriptors descriptors;
        descriptors.settings = mDescriptorsSettings;
        descriptors.set(atomicNumbers, positions, labels);
        assign(descriptors, atomsToAssign, typeID, lcs);
    }


    void MolecularAtomTypeAssigner::assign(
        const std::vector<int> &atomicNumbers, 
        const std::vector<Vector3d> &positions,
        const std::vector<std::string> &labels,
        std::vector<int> &typesID,
        std::vector<LocalCoordinateSystem<int> > &lcs)
        const
    {
        StructureWithDescriptors descriptors;
        descriptors.settings = mDescriptorsSettings;
        descriptors.set(atomicNumbers, positions, labels);
        assign(descriptors, typesID, lcs);
    }

    void MolecularAtomTypeAssigner::assign(
        const StructureWithDescriptors &descriptors,
        std::vector<int> &typeID,
        std::vector<LocalCoordinateSystem<int> > &lcs)
        const
    {
        int i, nAtoms = descriptors.atomDescriptors.size();
        vector<int> atomsToAssign;
        atomsToAssign.resize(nAtoms);
        for (i = 0; i < nAtoms; i++)
            atomsToAssign[i] = i;

        assign(descriptors, atomsToAssign, typeID, lcs);

    }

    void MolecularAtomTypeAssigner::assign_all_possible(
        const StructureWithDescriptors& descriptors,
        std::vector< std::vector<int> >& typeIDs)
    {
        int atomIndex, nAtoms;

        //SharedPointer<DescriptionDetailTreeNode>::type searchTreeNode, childNode;
        vector<int> atomMatch;


        nAtoms = descriptors.atomDescriptors.size();
        
        typeIDs.clear();
        typeIDs.resize(nAtoms);

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            // --- new code ---
            for (int i = 0; i < mAtomTypeMatchalgorithms.size(); i++)
            {
                LocalCoordinateSystem<int> lcs;
                if (mAtomTypeMatchalgorithms[i].match(i, descriptors, lcs))
                    typeIDs[atomIndex].push_back(i);
            }

        }

    }

    void MolecularAtomTypeAssigner::assign(
        const StructureWithDescriptors &descriptors,
        const std::vector<int> &atomsToAssign,
        std::vector<int> &typesID,
        std::vector<LocalCoordinateSystem<int> > &coordinateSystems)
        const
    {
        int atomIndex, nAtoms;

        //SharedPointer<DescriptionDetailTreeNode>::type searchTreeNode, childNode;
        vector<int> atomMatch;


        nAtoms = atomsToAssign.size();

        typesID.assign(nAtoms, -1);
        coordinateSystems.clear();
        coordinateSystems.resize(nAtoms);
           
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            // --- new code ---
            for (int i = 0; i < mAtomTypeMatchalgorithms.size(); i++)
            {
                if (mAtomTypeMatchalgorithms[i].match(atomsToAssign[atomIndex], descriptors, coordinateSystems[atomIndex]))
                {
                    typesID[atomIndex] = i;
                    break;
                }
            }

        }


            // --- old code ---
   //         // find atom type
			////cout << atomIndex + 1 << endl;
   //         searchTreeNode = mAtomTypeDescriptionTree.getRoot();
   //         nChildNodes = searchTreeNode->childNodes.size();

   //         foundType = false;

   //         while (nChildNodes != 0) // iterate over 'levels' of match tree
   //         {
   //             foundType = false;

   //             for (childNodeIndex = 0; childNodeIndex < nChildNodes; childNodeIndex++)
   //             {
   //                 childNode = searchTreeNode->childNodes[childNodeIndex];

   //                 isLeafNode = (childNode->typeID >= 0); // if true then the child node corresponds to indexed type from data bank

   //                 if (isLeafNode)
   //                     foundType = childNode->matchAlgorithm->match(atomsToAssign[atomIndex], descriptors, coordinateSystems[atomIndex]);
   //                 else
   //                     foundType = childNode->matchAlgorithm->match(atomsToAssign[atomIndex], descriptors);

   //                 if (foundType)
   //                 {
   //                     searchTreeNode = childNode;
   //                     break;
   //                 }
   //             }

   //             if (!foundType)
   //                 onAtomTypeNotFound(atomIndex, descriptors.connectivity, descriptors.atomDescriptors);

   //             //nChildNodes = searchTreeNode->childNodes.size();
   //             nChildNodes = 0;
   //         }

   //         if (foundType)
   //         {
   //             typesID[atomIndex] = childNode->typeID;
   //             onAtomTypeFound(atomIndex, typesID[atomIndex]);
   //         }
   //     }

    }

    
    
    //void MolecularAtomTypeAssigner::onAtomTypeNotFound(
    //    int atomIndex,
    //    const std::vector<std::vector<int> > &connectivity,
    //    const std::vector<AtomInStructureDescriptors> &atomsDescriptors)
    //    const
    //{
    //    //cout << "ATOM TYPE NOT FOUND " << atomIndex << endl;
    //}

    //void MolecularAtomTypeAssigner::onAtomTypeFound(
    //    int atomIndex,
    //    int typeID)
    //    const
    //{
    //    //cout << "ATOM TYPE FOUND " << atomIndex << endl;
    //    //cout << "     type index - " << typeID << endl;
    //}


}


