#include "discamb/Scattering/taam_utilities.h"
#include "discamb/MathUtilities/GraphAlgorithms.h"
#include "discamb/StructuralProperties/RingCalculator.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"
#include "discamb/BasicChemistry/PeriodicTable.h"

#include <set>
#include <map>
#include <algorithm>

using namespace std;

namespace {

    void checkIfAtomicNumbersRepresentedInWfnBank(
        const vector<int> &atomicNumbers,
        set<int> &presentZ,
        set<int> &missingZ)
    {
        discamb::ClementiRoettiData clementiRoettiData;
        discamb::HC_WfnBankEntry wfnEntry;
        bool hasType;
        for (auto z : atomicNumbers)
        {
            wfnEntry = clementiRoettiData.getEntry(discamb::periodic_table::symbol(z), hasType);
            if (!hasType)
                missingZ.insert(z);
            else
                presentZ.insert(z);
        }
    }

}

namespace discamb {
    namespace taam_utilities {

        void type_assignment_to_unscaled_HC_parameters(
            const vector<AtomTypeHC_Parameters> &bankMultipoleParameters,
            const vector<int> &atomTypeAssignment,
            const vector<int> &atomicNumbers,
            HC_ModelParameters &parameters,
            bool notAssignedRepresentedWithSlaters,
            vector<int> &nonMultipolarAtoms)
        {
            int atomIdx, nAtoms = atomTypeAssignment.size();
            

            // set the first type to void type if there is a need

            set<int> atomicNumbersInWfnBank, atomicNumbersNotInWfnBank;
            map<int, int> zToWfnIdx;
            checkIfAtomicNumbersRepresentedInWfnBank(atomicNumbers, atomicNumbersInWfnBank, atomicNumbersNotInWfnBank);
            // chceck if there are some atoms for which there are no slater wfn
            // make a void wfn type for them
            bool hasVoidType = !(atomicNumbersNotInWfnBank.empty());

            if (hasVoidType)
            {
                parameters.wfn_parameters.resize(1);
                parameters.type_parameters.resize(1);
                void_atom_type(parameters.wfn_parameters[0], parameters.type_parameters[0]);
            }
            
            // if some atoms have not assigned type
            // and they should be represented with Hansen-Coppens model
            // then types for them needs to be created
            // (or alternatively they are put in bank of type on the end of it)
            // the latter is used at the moment


            // set wfn types
            int nVoidTypes = hasVoidType ? 1 : 0;
            parameters.wfn_parameters.resize(atomicNumbersInWfnBank.size() + nVoidTypes);
            int counter = nVoidTypes;
            discamb::ClementiRoettiData clementiRoettiData;
            DeformationValenceParameters def_val(DeformationValenceParameters::ParametersType::UBDB);

            for (auto z : atomicNumbersInWfnBank)
            {
                zToWfnIdx[z] = counter;
                parameters.wfn_parameters[counter] = clementiRoettiData.getEntry(periodic_table::symbol(z));
                def_val.getParameters(periodic_table::symbol(z), parameters.wfn_parameters[counter].deformation_valence_exponent,
                    parameters.wfn_parameters[counter].deformation_valence_power);
                counter++;
            }

            // set bank types

            set<int> bankTypes(atomTypeAssignment.begin(), atomTypeAssignment.end());
            map<int, int> type2idx;
            bankTypes.erase(-1);
            counter = nVoidTypes;
            parameters.type_parameters.resize(bankTypes.size() + nVoidTypes);

            for (int idx : bankTypes)
            {
                type2idx[idx] = counter;
                parameters.type_parameters[counter].kappa_deformation_valence = bankMultipoleParameters[idx].kappa_prime;
                parameters.type_parameters[counter].kappa_spherical_valence = bankMultipoleParameters[idx].kappa;
                parameters.type_parameters[counter].p_val = bankMultipoleParameters[idx].p_val;
                int maxL = 0;
                bool hasPlm = !bankMultipoleParameters[idx].p_lm_indices.empty();
                for (auto p : bankMultipoleParameters[idx].p_lm_indices)
                    maxL = std::max(maxL, p.first);

                vector<vector<double> > &plm = parameters.type_parameters[counter].p_lm;
                if (hasPlm)
                {
                    plm.resize(maxL + 1);
                    for (int i = 0; i < int(maxL + 1); i++)
                        plm[i].assign(2 * i + 1, 0);
                }
                const vector<pair<int, int> > &plm_indices = bankMultipoleParameters[idx].p_lm_indices;

                for (int i = 0; i < bankMultipoleParameters[idx].p_lm_indices.size(); i++)
                {
                    int l_plus_m = plm_indices[i].second + int(plm_indices[i].first);
                    plm[plm_indices[i].first][l_plus_m] = bankMultipoleParameters[idx].p_lms[i];
                }
                counter++;
            }

            // atom to types map 
            parameters.atom_to_type_map.resize(nAtoms);
            parameters.atom_to_wfn_map.resize(nAtoms);
            nonMultipolarAtoms.clear();

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if (atomTypeAssignment[atomIdx] >= 0)
                {
                    parameters.atom_to_type_map[atomIdx] = type2idx[atomTypeAssignment[atomIdx]];
                    parameters.atom_to_wfn_map[atomIdx] = zToWfnIdx[atomicNumbers[atomIdx]];

                }
                else
                {  
                    parameters.atom_to_type_map[atomIdx] = 0;
                    
                    // atomicNumbersInWfnBank, atomicNumbersNotInWfnBank
                    if (notAssignedRepresentedWithSlaters && 
                        atomicNumbersInWfnBank.find(atomicNumbers[atomIdx]) != atomicNumbersInWfnBank.end())
                    {
                        parameters.atom_to_wfn_map[atomIdx] = zToWfnIdx[atomicNumbers[atomIdx]];
                    }
                    else
                    {
                        parameters.atom_to_wfn_map[atomIdx] = 0;
                        nonMultipolarAtoms.push_back(atomIdx);
                    }
                }
            }

            int nTypes = parameters.type_parameters.size();
            vector<double> typePval(nTypes);
            vector<double> typePvalSigma(nTypes);
            if (hasVoidType)
            {
                typePval[0] = 0.0;
                typePvalSigma[0] = 0.0;
            }

            counter = (hasVoidType ? 1 : 0);
            for (int idx : bankTypes)
            {
                typePval[counter] = bankMultipoleParameters[idx].p_val;
                typePvalSigma[counter] = bankMultipoleParameters[idx].p_val_sigma;
                counter++;
            }

            counter = (hasVoidType ? 1 : 0);
            for (int idx : bankTypes)
            {
                parameters.type_parameters[counter].p_val = typePval[counter];
                counter++;
            }


        }


        void type_assignment_to_HC_parameters(
            const vector<AtomTypeHC_Parameters>& bankMultipoleParameters,
            const vector<int>& atomTypeAssignment,
            const std::vector<double>& multiplicityTimesOccupancy,
            const vector<int>& atomicNumbers,
            double targetCharge,
            HC_ModelParameters& parameters,
            bool notAssignedRepresentedWithSlaters,
            vector<int>& nonMultipolarAtoms)
        {
            int atomIdx, nAtoms = atomTypeAssignment.size();


            // set the first type to void type if there is a need

            set<int> atomicNumbersInWfnBank, atomicNumbersNotInWfnBank;
            map<int, int> zToWfnIdx;
            checkIfAtomicNumbersRepresentedInWfnBank(atomicNumbers, atomicNumbersInWfnBank, atomicNumbersNotInWfnBank);
            // chceck if there are some atoms for which there are no slater wfn
            // make a void wfn type for them
            bool hasVoidType = !(atomicNumbersNotInWfnBank.empty());

            if (hasVoidType)
            {
                parameters.wfn_parameters.resize(1);
                parameters.type_parameters.resize(1);
                void_atom_type(parameters.wfn_parameters[0], parameters.type_parameters[0]);
            }

            // if some atoms have not assigned type
            // and they should be represented with Hansen-Coppens model
            // then types for them needs to be created
            // (or alternatively they are put in bank of type on the end of it)
            // the latter is used at the moment


            // set wfn types
            int nVoidTypes = hasVoidType ? 1 : 0;
            parameters.wfn_parameters.resize(atomicNumbersInWfnBank.size() + nVoidTypes);
            int counter = nVoidTypes;
            discamb::ClementiRoettiData clementiRoettiData;
            DeformationValenceParameters def_val(DeformationValenceParameters::ParametersType::UBDB);

            for (auto z : atomicNumbersInWfnBank)
            {
                zToWfnIdx[z] = counter;
                parameters.wfn_parameters[counter] = clementiRoettiData.getEntry(periodic_table::symbol(z));
                def_val.getParameters(periodic_table::symbol(z), parameters.wfn_parameters[counter].deformation_valence_exponent,
                    parameters.wfn_parameters[counter].deformation_valence_power);
                counter++;
            }

            // set bank types

            set<int> bankTypes(atomTypeAssignment.begin(), atomTypeAssignment.end());
            map<int, int> type2idx;
            bankTypes.erase(-1);
            counter = nVoidTypes;
            parameters.type_parameters.resize(bankTypes.size() + nVoidTypes);

            for (int idx : bankTypes)
            {
                type2idx[idx] = counter;
                parameters.type_parameters[counter].kappa_deformation_valence = bankMultipoleParameters[idx].kappa_prime;
                parameters.type_parameters[counter].kappa_spherical_valence = bankMultipoleParameters[idx].kappa;
                parameters.type_parameters[counter].p_val = bankMultipoleParameters[idx].p_val;
                int maxL = -1;
                for (auto p : bankMultipoleParameters[idx].p_lm_indices)
                    maxL = std::max(maxL, int(p.first));

                vector<vector<double> >& plm = parameters.type_parameters[counter].p_lm;
                if (maxL > -1)
                {
                    int maxL_plus_one = maxL + 1;
                    plm.resize(maxL_plus_one);
                    for (int i = 0; i < int(maxL_plus_one); i++)
                        plm[i].assign(2 * i + 1, 0);
                }
                const vector<pair<int, int> >& plm_indices = bankMultipoleParameters[idx].p_lm_indices;

                for (int i = 0; i < bankMultipoleParameters[idx].p_lm_indices.size(); i++)
                {
                    int l_plus_m = plm_indices[i].second + int(plm_indices[i].first);
                    plm[plm_indices[i].first][l_plus_m] = bankMultipoleParameters[idx].p_lms[i];
                }
                counter++;
            }

            // atom to types map 
            parameters.atom_to_type_map.resize(nAtoms);
            parameters.atom_to_wfn_map.resize(nAtoms);
            nonMultipolarAtoms.clear();

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if (atomTypeAssignment[atomIdx] >= 0)
                {
                    parameters.atom_to_type_map[atomIdx] = type2idx[atomTypeAssignment[atomIdx]];
                    parameters.atom_to_wfn_map[atomIdx] = zToWfnIdx[atomicNumbers[atomIdx]];

                }
                else
                {
                    parameters.atom_to_type_map[atomIdx] = 0;

                    // atomicNumbersInWfnBank, atomicNumbersNotInWfnBank
                    if (notAssignedRepresentedWithSlaters &&
                        atomicNumbersInWfnBank.find(atomicNumbers[atomIdx]) != atomicNumbersInWfnBank.end())
                    {
                        parameters.atom_to_wfn_map[atomIdx] = zToWfnIdx[atomicNumbers[atomIdx]];
                    }
                    else
                    {
                        parameters.atom_to_wfn_map[atomIdx] = 0;
                        nonMultipolarAtoms.push_back(atomIdx);
                    }
                }
            }

            int nTypes = parameters.type_parameters.size();
            vector<double> typePval(nTypes);
            vector<double> typePvalSigma(nTypes);
            if (hasVoidType)
            {
                typePval[0] = 0.0;
                typePvalSigma[0] = 0.0;
            }

            counter = (hasVoidType ? 1 : 0);
            for (int idx : bankTypes)
            {
                typePval[counter] = bankMultipoleParameters[idx].p_val;
                typePvalSigma[counter] = bankMultipoleParameters[idx].p_val_sigma;
                counter++;
            }

            //vector<int> atomToTypeMap(parameters.atom_to_type_map.begin(), parameters.atom_to_type_map.end());
            vector<int> atomToTypeMap;
            for (auto a2t : parameters.atom_to_type_map)
                atomToTypeMap.push_back(static_cast<int>(a2t));
            
            //electroneutrality_Faerman_Price(typePval, typePvalSigma, atomToTypeMap, atomicNumbers);
            electroneutrality_Faerman_Price(typePval, typePvalSigma, atomToTypeMap, multiplicityTimesOccupancy, atomicNumbers, targetCharge);


            counter = (hasVoidType ? 1 : 0);
            for (int idx : bankTypes)
            {
                parameters.type_parameters[counter].p_val = typePval[counter];
                counter++;
            }


        }



        void void_atom_type(
            HC_WfnType &wfnType,
            HC_AtomTypeParameters &type_parameters)
        {
            discamb::ClementiRoettiData clementiRoettiData;

            wfnType = clementiRoettiData.getEntry("H");
            wfnType.deformation_valence_exponent = 1.0;

            type_parameters.kappa_deformation_valence = 1.0;
            type_parameters.kappa_spherical_valence = 1.0;
            type_parameters.p_val = 0.0;
            type_parameters.p_lm.clear();
        }

        //void electroneutrality_Faerman_Price(
        //    std::vector<double> &typePval,
        //    const std::vector<double> &_typePvalSigma,
        //    const std::vector<int> &atomToTypeMap,
        //    const std::vector <int> atomicNumbers)
        //{
        //    int atomIdx, nAtoms = atomToTypeMap.size();
        //    double dQ, sigmaSum;
        //    vector<double> valenceOccupancy(37, 0);
        //    std::vector<double> typePvalSigma = _typePvalSigma;
        //    vector<bool> typeUsed(typePval.size(), false);
        //    discamb::ClementiRoettiData clementiRoettiData;
        //    vector<HC_WfnBankEntry> entries;
        //    clementiRoettiData.getEntries(entries);
        //    set<int> allowedZ;
        //    for (auto &entry : entries)
        //        if (entry.charge == 0)
        //        {
        //            for (auto orbIdx : entry.valence_orbitals_indices)
        //                valenceOccupancy[entry.atomic_number] += entry.orbital_occupancy[orbIdx];
        //            allowedZ.insert(entry.atomic_number);
        //        }
        //    
        //    dQ = 0;
        //    sigmaSum = 0;
        //    vector<double> atomDq(nAtoms);
        //    vector<double> atomPval(nAtoms);
        //    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        //        if (atomToTypeMap[atomIdx] >= 0 && (allowedZ.find(atomicNumbers[atomIdx])!=allowedZ.end()) )
        //        {
        //            atomPval[atomIdx] = typePval[atomToTypeMap[atomIdx]];
        //            atomDq[atomIdx] = valenceOccupancy[atomicNumbers[atomIdx]] - typePval[atomToTypeMap[atomIdx]];
        //            dQ += valenceOccupancy[atomicNumbers[atomIdx]] - typePval[atomToTypeMap[atomIdx]];
        //            sigmaSum += typePvalSigma[atomToTypeMap[atomIdx]];
        //            typeUsed[atomToTypeMap[atomIdx]] = true;
        //        }

        //    if (dQ == 0)
        //        return;

        //    if (sigmaSum == 0)
        //    {
        //        sigmaSum = nAtoms;
        //        typePvalSigma.assign(typePvalSigma.size(), 1.0);
        //    }
        //    
        //    for (int i = 0; i < typePval.size(); i++)
        //        if (typeUsed[i])
        //            typePval[i] += dQ * typePvalSigma[i] / sigmaSum;


        //}

        void electroneutrality_Faerman_Price(
            std::vector<double>& typePval,
            const std::vector<double>& _typePvalSigma,
            const std::vector<int>& atomToTypeMap,
            const std::vector<double>& multiplicity,
            const std::vector <int> atomicNumbers,
            double targetCharge)
        {
            int atomIdx, nAtoms = atomToTypeMap.size();
            double dQ, sigmaSum;
            vector<double> valenceOccupancy(37, 0);
            std::vector<double> typePvalSigma = _typePvalSigma;
            vector<bool> typeUsed(typePval.size(), false);
            discamb::ClementiRoettiData clementiRoettiData;
            vector<HC_WfnBankEntry> entries;
            clementiRoettiData.getEntries(entries);
            set<int> allowedZ;
            for (auto& entry : entries)
                if (entry.charge == 0)
                {
                    for (auto orbIdx : entry.valence_orbitals_indices)
                        valenceOccupancy[entry.atomic_number] += entry.orbital_occupancy[orbIdx];
                    allowedZ.insert(entry.atomic_number);
                }

            dQ = 0;
            sigmaSum = 0;
            double currentQ = 0;
            vector<double> atomDq(nAtoms);
            vector<double> atomPval(nAtoms);
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (atomToTypeMap[atomIdx] >= 0 && (allowedZ.find(atomicNumbers[atomIdx]) != allowedZ.end()))
                {
                    atomPval[atomIdx] = typePval[atomToTypeMap[atomIdx]];
                    atomDq[atomIdx] = valenceOccupancy[atomicNumbers[atomIdx]] - typePval[atomToTypeMap[atomIdx]];
                    //dQ += (valenceOccupancy[atomicNumbers[atomIdx]] - typePval[atomToTypeMap[atomIdx]]) * multiplicity[atomIdx];
                    currentQ += (valenceOccupancy[atomicNumbers[atomIdx]] - typePval[atomToTypeMap[atomIdx]]) * multiplicity[atomIdx];
                    sigmaSum += typePvalSigma[atomToTypeMap[atomIdx]] * multiplicity[atomIdx];
                    typeUsed[atomToTypeMap[atomIdx]] = true;
                }
            dQ = targetCharge - currentQ;
            if (dQ == 0)
                return;

            if (sigmaSum == 0)
            {
                sigmaSum = double(nAtoms);
                typePvalSigma.assign(typePvalSigma.size(), 1.0);
            }

            for (int i = 0; i < typePval.size(); i++)
                if (typeUsed[i])
                    //typePval[i] += dQ * typePvalSigma[i] / sigmaSum;
                    typePval[i] -= dQ * typePvalSigma[i] / sigmaSum;


        }


        int atomTypeRange(
            const AtomType &type,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range)
        {
            namedNeighboursRange = 0;
            planarRingsRange = 0;
            ring34range = 0;
            


            int shellIdx, nShells, nAtoms = type.atoms.size();
            vector<vector<int> > graphShells;
            vector<int> n3rings(nAtoms,0), n4rings(nAtoms,0);
            set<set<int> > rings3, rings4;
            graph_algorithms::breadth_first_search(type.connectivity, 0, graphShells);
            int range = 0;
            nShells = graphShells.size();
            RingCalculator::calculate34Rings(type.connectivity, rings3, rings4);

            for (auto &ring : rings3)
                for (int atomIdx : ring)
                    n3rings[atomIdx]++;

            for (auto &ring : rings4)
                for (int atomIdx : ring)
                    n4rings[atomIdx]++;

            range = nShells - 1;

            
            // rings, neighbours, shell

            /*
            1. ringInfo.in3Ring and ringInfo.n3rings
            2. ringInfo.in4Ring and ringInfo.n4rings
            3. planar rings
            4. shell
            5. neighbours
            */
            int atom3ringRange;
            int atom4ringRange;
            int atomPlanarRingRange;

            for(shellIdx=0; shellIdx<nShells; shellIdx++)
                for (int atomIdx : graphShells[shellIdx])
                {
                    AtomRingInfo const &ringInfo = type.atoms[atomIdx].ringInfo;
                    // 3 ring
                    atom3ringRange = 0;
                    // if not in ring we have to check if there is a ring
                    // if undefined we do not check anything
                    // if in ring we check for ring unless nRings = number of rings from connectivity
                    if (ringInfo.in3Ring != Tribool::Undefined &&
                        !(ringInfo.n3rings >= 0 && n3rings[atomIdx] == ringInfo.n3rings))
                        atom3ringRange = shellIdx + 1;
                    // 4 ring
                    atom4ringRange = 0;
                    if (ringInfo.in4Ring != Tribool::Undefined &&
                        !(ringInfo.n4rings >= 0 && n4rings[atomIdx] == ringInfo.n4rings))
                        atom3ringRange = shellIdx + 2;
                    // planar ring
                    atomPlanarRingRange = 0;
                    if (ringInfo.inRing != Tribool::Undefined)
                    {
                        int maxRing;
                        if (ringInfo.inRing == Tribool::False)
                            maxRing = maxPlanarRing;
                        else
                        {
                            if (ringInfo.inAnyAdditionalRing)
                                maxRing = maxPlanarRing;
                            else
                            {
                                int n1, n2;
                                n1 = n2 = 0;
                                if (!ringInfo.labeledContainingRings.empty())
                                    n1 = std::max_element(ringInfo.labeledContainingRings.begin(), ringInfo.labeledContainingRings.end())->first;
                                if (!ringInfo.nonLabeledContainingRings.empty())
                                    n2 = *std::max_element(ringInfo.nonLabeledContainingRings.begin(), ringInfo.nonLabeledContainingRings.end());
                                maxRing = std::max(n1, n2);
                            }
                        }
                        atomPlanarRingRange += (maxRing + 1) / 2 + 1 + shellIdx;
                    }

                    planarRingsRange = max(planarRingsRange, atomPlanarRingRange);
                    ring34range = max( {ring34range, atom3ringRange, atom4ringRange });
                }

            namedNeighboursRange = nShells - 1;
            return max({ planarRingsRange, ring34range, namedNeighboursRange+1 });
        }

        int atomTypesRange(
            const std::vector<AtomType> &types,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range)
        {

            int atomRange34ring, atomPlanarRingRange, atomNeighbourRange, atomTotalRange;

            ring34range = 0;
            planarRingsRange = 0;
            namedNeighboursRange = 0;
            int totalRange = 0;

            for (int i = 0; i < types.size(); i++)
            {
                atomTotalRange = taam_utilities::atomTypeRange(types[i], maxPlanarRing,
                    atomNeighbourRange, atomPlanarRingRange, atomRange34ring);
                
                namedNeighboursRange = max(atomNeighbourRange, namedNeighboursRange);
                planarRingsRange = max(atomPlanarRingRange, planarRingsRange);
                ring34range = max(atomRange34ring, ring34range);
                totalRange = max(atomTotalRange, totalRange);
            }
            return totalRange;
        }

    }

}