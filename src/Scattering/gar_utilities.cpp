#include "discamb/Scattering/gar_utilities.h"

#include "discamb/BasicUtilities/utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/StructuralProperties/atom_selection.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/StructuralProperties/SimpleAIF_ScoreCalculator.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/MathUtilities/lebedev_laikov.h"
#include "discamb/MathUtilities/radial_grid.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {
    namespace gar_utilities {

        void findDeaultSubsystems(
            const Crystal& crystal,
            std::vector<std::vector<std::pair<std::string, std::string> > >& clusterAtoms,
            std::vector<std::string>& clustersLabels,
            std::vector<int>& clustersCharges,
            std::vector<int>& clustersSpinMultiplicity,
            int spinMultiplicityHint)
        {
            UnitCellContent uc;
            vector<vector<UnitCellContent::AtomID> > molecules;
            vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
            uc.set(crystal);

            structural_properties::splitAsymmetricUnitIntoMolecules(uc, molecules, networkBonds, 0.4);
            int nMolecules = molecules.size();
            int nAtomsInAsymmetricUnit = crystal.atoms.size();
            clusterAtoms.resize(nMolecules);

            string atomLabel, symmOp;


            for (int moleculeIdx = 0; moleculeIdx < nMolecules; moleculeIdx++)
            {
                for (auto& atom : molecules[moleculeIdx])
                {
                    uc.interpreteAtomID(atom, atomLabel, symmOp);
                    clusterAtoms[moleculeIdx].push_back({ atomLabel, symmOp });
                }
                clustersLabels.push_back("mol_" + to_string(moleculeIdx + 1));
                clustersCharges.push_back(0);

                if (spinMultiplicityHint == 0)
                    clustersSpinMultiplicity.push_back(1);
                else
                    clustersSpinMultiplicity.push_back(spinMultiplicityHint);

            }
        }

        void findDefaultRepresentatives(
            const Crystal& crystal,
            const std::vector<std::vector<std::pair<std::string, std::string> > >& subsystemAtoms,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives)
        {
            vector<vector<pair<string, string> > > representativesPerCluster;
            findDefaultRepresentatives(crystal, subsystemAtoms, representatives, representativesPerCluster);
        }


        void findDefaultRepresentatives(
            const Crystal& crystal,
            const vector<vector<pair<string, string> > >& clusterAtoms,
            vector<vector<AtomRepresentativeInfo> >& representatives,
            vector<vector<pair<string, string> > >& representativesPerCluster)
        {
            representatives.clear();
            representativesPerCluster.clear();
            vector<vector<double> > scores;
            SimpleAIF_ScoreCalculator scoreCalculator;
            scoreCalculator.assignScore(crystal, clusterAtoms, scores);

            int atomInClusterIdx, clusterIdx, nClusters = clusterAtoms.size();
            int nAtomsInAsymmetricUnit = crystal.atoms.size();


            // [atom idx][i] - {score, {cluster, atomIdx} }
            vector<vector<pair<double, pair<int, int> > > > representativesAndScores(nAtomsInAsymmetricUnit);

            representatives.resize(nAtomsInAsymmetricUnit);
            representativesPerCluster.resize(nClusters);
            //vector<bool> hasXyzRepresentativeInClusters(nAtomsInAsymmetricUnit, false);

            map<string, int> label2idx;
            for (int i = 0; i < crystal.atoms.size(); i++)
                label2idx[crystal.atoms[i].label] = i;

            for (int clusterIdx = 0; clusterIdx < nClusters; clusterIdx++)
                for (int atomInClusterIdx = 0; atomInClusterIdx < clusterAtoms[clusterIdx].size(); atomInClusterIdx++)
                    if (clusterAtoms[clusterIdx][atomInClusterIdx].first.find('@') == string::npos)
                        representativesAndScores[label2idx[clusterAtoms[clusterIdx][atomInClusterIdx].first]].push_back(
                            { scores[clusterIdx][atomInClusterIdx], { clusterIdx, atomInClusterIdx } });


            for (int atomIdx = 0; atomIdx < nAtomsInAsymmetricUnit; atomIdx++)
            {
                if (representativesAndScores[atomIdx].empty())
                    on_error::throwException("no represrentative for atom " + crystal.atoms[atomIdx].label, __FILE__, __LINE__);

                auto representative = *max_element(representativesAndScores[atomIdx].begin(), representativesAndScores[atomIdx].end());

                double bestScore = representative.first;

                vector<pair<int, int> > bestReps;

                for (auto rep : representativesAndScores[atomIdx])
                    if (rep.first == bestScore)
                        bestReps.push_back(rep.second);

                for (auto rep : bestReps)
                {
                    atomInClusterIdx = rep.second;
                    clusterIdx = rep.first;

                    auto& atom = clusterAtoms[clusterIdx][atomInClusterIdx];

                    representativesPerCluster[clusterIdx].push_back(atom);

                    int idx = representatives[atomIdx].size();
                    representatives[atomIdx].resize(idx + 1);
                    representatives[atomIdx][idx].subsystemIdx = clusterIdx;
                    representatives[atomIdx][idx].idxInSubsystem = atomInClusterIdx;
                    representatives[atomIdx][idx].atomLabel = atom.first;
                    representatives[atomIdx][idx].fixedWeightValue = 1.0;
                    representatives[atomIdx][idx].isWeightFixed = true;
                    representatives[atomIdx][idx].symmetryCode = atom.second;
                    representatives[atomIdx][idx].transformType = string("symmetry");
                    representatives[atomIdx][idx].transformSetting = string("invert ") + atom.second;
                }
            }



        }

        void readRepresentatives(
            const string& fileName,
            const Crystal& crystal,
            const std::vector<std::string>& subsystemLabels,
            const std::vector< std::vector<std::pair<std::string, std::string> > >& subsystemAtoms,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives)
        {

        }


        //[z][point idx]
        void elementalIntegrationGrids(
            const std::set<int> &atomicNumbers,
            std::vector<std::vector<Vector3d> > &elementslGridPoints,
            std::vector<std::vector<double> >& elementalIntegrationWeights)
        {

            int nAngular = 590;
            int nRadial = 99;

            elementslGridPoints.clear();
            elementalIntegrationWeights.clear();


            elementslGridPoints.resize(113);
            elementalIntegrationWeights.resize(113);

            vector<Vector3d> angularGridPoints, gridPoints;
            vector<double> radialPoints, radialWeights, angularWeights, weights;
            Vector3d r, atomPosition;

            lebedev_laikov::get_grid(nAngular, angularGridPoints, angularWeights);

            gridPoints.resize(nAngular * nRadial);
            weights.resize(nAngular * nRadial);

            for (auto z : atomicNumbers)
            {
                radial_grid::treutler_ahlrichs(z, nRadial, radialPoints, radialWeights);
                int pointIdx = 0;
                for (int angularIdx = 0; angularIdx < nAngular; angularIdx++)
                    for (int radialIdx = 0; radialIdx < nRadial; radialIdx++)
                    {
                        weights[pointIdx] = radialWeights[radialIdx] * angularWeights[angularIdx] *
                            radialPoints[radialIdx] * radialPoints[radialIdx] * M_PI * 4.0;
                        gridPoints[pointIdx] = radialPoints[radialIdx] * angularGridPoints[angularIdx];
                        pointIdx++;
                    }
                elementslGridPoints[z] = gridPoints;
                elementalIntegrationWeights[z] = weights;
            }


        }

        /*for one fragment*/
        //void findMultipoleClasterAtoms(
        //    const std::optional<DistributedMultipoleCentersSettings>& settings,
        //    double threshold,
        //    const UnitCellContent& ucContent,
        //    const std::vector<std::vector<UnitCellContent::AtomID> >& unitCellMolecules,
        //    const FragmentAtoms& fragmentAtoms,
        //    std::vector<std::pair<std::string, std::string > >& multipoleClusterAtoms,
        //    std::vector<std::optional<double> >& multipoleCustomWeights)
        //{
        //    multipoleCustomWeights.clear();
        //    multipoleClusterAtoms.clear();
        //    
        //    pair<string, string> atomLabelAndSymmOp;
        //    vector<UnitCellContent::AtomID> fragment, expandedFragment, multipoleCluster;
        //    UnitCellContent::AtomID atomId;

        //    //  filling expandedFragment with qm system atoms 
        //    // (capping H atoms replaced by directing atoms)
        //    //  and similarily fragment (but no capping H related atoms included)

        //    crystal_structure_utilities::convertAtomList(ucContent, fragmentAtoms.atomList, fragment);
        //    expandedFragment = fragment;

        //    for (auto& capH : fragmentAtoms.cappingHydrogens)
        //    {
        //        ucContent.findAtom(capH.directingAtom, capH.directingAtomSymmOp, atomId);
        //        expandedFragment.push_back(atomId);
        //    }          

        //    // multipole bearing atoms and cluster atoms will be saved in multipoleCluster
        //    
        //    bool hasPredefinedAtoms;
        //    if(settings)
        //        hasPredefinedAtoms = !(settings.value().predefinedMultipoleClusterAtoms.empty());

        //    if (!hasPredefinedAtoms)
        //        structural_properties::makeCluster(ucContent, fragment, unitCellMolecules, multipoleCluster, threshold, false);
        //    else // multipole bearing atoms only taken from predefined list will be saved in multipoleCluster
        //        crystal_structure_utilities::convertAtomList(ucContent, settings.value().predefinedMultipoleClusterAtoms, multipoleCluster);

        //    // remove atoms from multipoleCluster which belong to qm system

        //    AtomSubsetSelectionData qmSystemSelection;

        //    for (auto const& atom : expandedFragment)
        //    {
        //        ucContent.interpreteAtomID(atomId, atomLabelAndSymmOp.first, atomLabelAndSymmOp.second);
        //        qmSystemSelection.atomList.push_back(atomLabelAndSymmOp);
        //    }
        //    atom_selection::remove_subset(ucContent, multipoleCluster, qmSystemSelection);

        //    if (settings)
        //        for (auto& subsetSelectionArgs : settings.value().atomsToOmit)
        //            atom_selection::remove_subset(ucContent, multipoleCluster, subsetSelectionArgs);

        //    // it is final set of atoms in cluster, translating the list to other format
        //    crystal_structure_utilities::convertAtomList(ucContent, multipoleCluster, multipoleClusterAtoms);
        //    multipoleCustomWeights.resize(multipoleCluster.size());

        //    // assign custom weights 
        //    if(settings)
        //        for (auto& wght : settings.value().atomsWithCustomWeights)
        //        {
        //            vector<UnitCellContent::AtomID> selection;
        //            atom_selection::select_subset(ucContent, multipoleCluster, wght.first, selection);
        //            for (auto item : selection)
        //                multipoleCustomWeights[utilities::index(multipoleCluster, item)] = wght.second;
        //        }
        //}

        /*for multiple fragments*/
        //void findMultipoleClasterAtoms(
        //    const vector< optional<DistributedMultipoleCentersSettings> >& settings,
        //    double threshold,
        //    const Crystal& crystal,
        //    const vector<FragmentAtoms>& fragmentAtoms,
        //    vector< vector< pair< string, string > > >& multipoleClusterAtoms,
        //    vector< vector< optional<double> > >& multipoleCustomWeight)
        //{
        //    multipoleClusterAtoms.clear();
        //    multipoleCustomWeight.clear();

        //    UnitCellContent ucContent(crystal);
        //    vector< vector< UnitCellContent::AtomID > > unitCellMolecules;
        //    vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
        //    structural_properties::splitUnitCellIntoMolecules(ucContent, unitCellMolecules, networkBonds);
        //    int fragmentIdx, nFragments = fragmentAtoms.size();
        //    multipoleClusterAtoms.resize(nFragments);
        //    multipoleCustomWeight.resize(nFragments);

        //    for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
        //        findMultipoleClasterAtoms(
        //            settings[fragmentIdx], 
        //            threshold,
        //            ucContent,
        //            unitCellMolecules,
        //            fragmentAtoms[fragmentIdx],
        //            multipoleClusterAtoms[fragmentIdx],
        //            multipoleCustomWeight[fragmentIdx]);
        //}


    } 
}

