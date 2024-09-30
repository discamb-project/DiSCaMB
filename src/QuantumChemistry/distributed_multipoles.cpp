#include "discamb/QuantumChemistry/distributed_multipoles.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/utilities.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include <fstream>

using namespace std;

namespace {
    double delta(int i, int j)
    {
        if (i == j)
            return 1.0;
        return 0.0;
    }

    bool similarVectors(
        const discamb::Vector3d& v1,
        const discamb::Vector3d& v2,
        double threshold)
    {
        if (fabs(v1[0] - v2[0]) > threshold)
            return false;
        if (fabs(v1[1] - v2[1]) > threshold)
            return false;
        if (fabs(v1[2] - v2[2]) > threshold)
            return false;
        return true;
    }


    void packPointCharges(
        const vector<discamb::Vector3d>& r0,
        const vector<double>& q0,
        vector<discamb::Vector3d>& r,
        vector<double>& q,
        double threshold)
    {
        r.clear();
        q.clear();

        for (int i = 0; i < r0.size(); i++)
        {
            bool merged = false;
            for (int j = 0; j < r.size(); j++)
            {
                if (similarVectors(r0[i], r[j], threshold))
                {
                    q[j] += q0[i];
                    merged = true;
                    break;
                }
            }

            if (!merged)
            {
                q.push_back(q0[i]);
                r.push_back(r0[i]);
            }

        }

    }

}

namespace discamb {

    void DistributedMultipoleCentersSettings::set(
        nlohmann::json const& jsonData)
    {
        *this = DistributedMultipoleCentersSettings();
        //
        multipoleClusterThreshold = jsonData.value("threshold", multipoleClusterThreshold);
        if (jsonData.find("atoms file") != jsonData.end())
        {
            string fileName = jsonData["atoms file"].get<string>();
            ifstream in(fileName);
            string line;
            vector<string> words;
            if (!in.good())
                on_error::throwException("cannot read multipole centers \"atoms file\": \"" + fileName + "\"", __FILE__, __LINE__);
            while(getline(in, line))
            {
                
                string_utilities::split(line, words);
                if (words.size() == 2)
                    predefinedMultipoleClusterAtoms.push_back({ words[0], words[1] });
                if (words.size() == 1)
                    predefinedMultipoleClusterAtoms.push_back({ words[0], string("X,Y,Z") });

            }
            in.close();
        }
        //
        if (jsonData.find("remove") != jsonData.end())
        {
            if(jsonData.find("remove")->is_array())
                for (auto const& item : jsonData.find("remove").value())
                {
                    AtomSubsetSelectionData subsetSelectionData;
                    subsetSelectionData.set(item);
                    atomsToOmit.push_back(subsetSelectionData);
                }
            else
            {
                AtomSubsetSelectionData subsetSelectionData;
                subsetSelectionData.set(jsonData.find("remove").value());
                atomsToOmit.push_back(subsetSelectionData);
            }
        }


            
            //this->atomsToOmit
    }

    namespace distributed_multipoles {

        void get_fragment_multipole_centers(
            const DistributedMultipoleCentersSettings& settings,
            const UnitCellContent& ucContent,
            const std::vector<std::vector<UnitCellContent::AtomID> >& unitCellMolecules,
            const FragmentAtoms& fragmentAtoms,
            std::vector<std::pair<std::string, std::string> >& multipoleClusterAtoms,
            std::vector<std::optional<double> >& multipoleCustomWeights)
        {
            multipoleCustomWeights.clear();
            multipoleClusterAtoms.clear();

            pair<string, string> atomLabelAndSymmOp;
            vector<UnitCellContent::AtomID> fragment, expandedFragment, multipoleCluster;
            UnitCellContent::AtomID atomId;

            //  filling expandedFragment with qm system atoms 
            // (capping H atoms replaced by directing atoms)
            //  and similarily fragment (but no capping H related atoms included)

            crystal_structure_utilities::convertAtomList(ucContent, fragmentAtoms.atomList, fragment);
            expandedFragment = fragment;

            for (auto& capH : fragmentAtoms.cappingHydrogens)
            {
                ucContent.findAtom(capH.directingAtom, capH.directingAtomSymmOp, atomId);
                expandedFragment.push_back(atomId);
            }

            // multipole bearing atoms and cluster atoms will be saved in multipoleCluster

            //bool hasPredefinedAtoms = !(settings.predefinedMultipoleClusterAtoms.empty());

            if (settings.predefinedMultipoleClusterAtoms.empty())
                structural_properties::makeCluster(ucContent, fragment, unitCellMolecules, multipoleCluster, settings.multipoleClusterThreshold, false);
            else // multipole bearing atoms only taken from predefined list will be saved in multipoleCluster
                crystal_structure_utilities::convertAtomList(ucContent, settings.predefinedMultipoleClusterAtoms, multipoleCluster);

            // remove atoms from multipoleCluster which belong to qm system

            AtomSubsetSelectionData qmSystemSelection;

            for (auto const& atomId : expandedFragment)
            {
                ucContent.interpreteAtomID(atomId, atomLabelAndSymmOp.first, atomLabelAndSymmOp.second);
                qmSystemSelection.atomList.push_back(atomLabelAndSymmOp);
            }
            atom_selection::remove_subset(ucContent, multipoleCluster, qmSystemSelection);

            for (auto& subsetSelectionArgs : settings.atomsToOmit)
                atom_selection::remove_subset(ucContent, multipoleCluster, subsetSelectionArgs);

            // it is final set of atoms in cluster, translating the list to other format
            crystal_structure_utilities::convertAtomList(ucContent, multipoleCluster, multipoleClusterAtoms);
            multipoleCustomWeights.resize(multipoleCluster.size());

            // assign custom weights 

            for (auto& wght : settings.atomsWithCustomWeights)
            {
                vector<UnitCellContent::AtomID> selection;
                atom_selection::select_subset(ucContent, multipoleCluster, wght.first, selection);
                for (auto item : selection)
                    multipoleCustomWeights[utilities::index(multipoleCluster, item)] = wght.second;
            }
        }

        void get_fragments_multipole_centers(
            const std::vector<std::optional<DistributedMultipoleCentersSettings> >& settings,
            double threshold,
            const Crystal& crystal,
            const std::vector<FragmentAtoms>& fragmentAtoms,
            std::vector< std::vector<std::pair<std::string, std::string > > >& multipoleClusterAtoms,
            std::vector < std::vector<std::optional<double> > >& multipoleCustomWeight)
        {
            multipoleClusterAtoms.clear();
            multipoleCustomWeight.clear();

            int fragmentIdx, nFragments = fragmentAtoms.size();
            
            UnitCellContent ucContent(crystal);
            vector< vector< UnitCellContent::AtomID > > unitCellMolecules;
            vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
            structural_properties::splitUnitCellIntoMolecules(ucContent, unitCellMolecules, networkBonds);
            multipoleClusterAtoms.resize(nFragments);
            multipoleCustomWeight.resize(nFragments);
            DistributedMultipoleCentersSettings generalSettings;
            generalSettings.multipoleClusterThreshold = threshold;
            for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
            {
                bool hasCustomSettings = false;
                if (!settings.empty())
                    hasCustomSettings = settings[fragmentIdx].has_value();
                auto const& dmSettings = (hasCustomSettings ? settings[fragmentIdx].value() : generalSettings);
                get_fragment_multipole_centers(
                    dmSettings,
                    ucContent,
                    unitCellMolecules,
                    fragmentAtoms[fragmentIdx],
                    multipoleClusterAtoms[fragmentIdx],
                    multipoleCustomWeight[fragmentIdx]);
            }
        }



        void calculate_dipole(
            const std::vector<Vector3d>& r,
            const std::vector<double>& charges,
            Vector3d& dipole)
        {
            int i, n = r.size();
            dipole.set(0, 0, 0);
            for (i = 0; i < n; i++)
                dipole += charges[i] * r[i];
        }

        void calculate_quadrupole(
            const std::vector<Vector3d>& r,
            const std::vector<double>& charges,
            Matrix3d& quadrupole)
        {
            int i, u, w, n = r.size();
            double r2;
            quadrupole.set(0, 0, 0, 0, 0, 0, 0, 0, 0);

            for (i = 0; i < n; i++)
            {
                r2 = r[i] * r[i];
                for (u = 0; u < 3; u++)
                    for (w = 0; w < 3; w++)
                        quadrupole(u, w) += charges[i] *
                        (3 * r[i][u] * r[i][w] - delta(u, w) * r2);
            }

        }

        void calculate_multipoles(
            const std::vector<Vector3d>& r,
            const std::vector<double>& charges,
            ElectricMultipoles& multipoles,
            int max_l)
        {
            int i, u, w, n = int(r.size());
            double r2;
            multipoles = ElectricMultipoles();

            for (i = 0; i < n; i++)
            {
                multipoles.charge += charges[i];
                if (max_l > 0)
                    multipoles.dipole += charges[i] * r[i];
                if (max_l > 1)
                {
                    r2 = r[i] * r[i];
                    for (u = 0; u < 3; u++)
                        for (w = 0; w < 3; w++)
                            multipoles.quadrupole(u, w) += charges[i] *
                            (3 * r[i][u] * r[i][w] - delta(u, w) * r2);
                }
            }

        }

        // distance - distance from center to charge
        void dipole_as_point_charges(
            const Vector3d& dipole,
            double distance,
            std::vector<Vector3d>& r,
            std::vector<double>& charges)
        {


            double abs_value = sqrt(dipole * dipole);
            Vector3d v;
            if (abs_value == 0.0)
                v = distance * Vector3d(1, 0, 0);
            else
                v = distance * dipole / abs_value;
            double q = abs_value / (2.0 * distance);
            r.clear();
            r.push_back(v);
            r.push_back(-v);
            charges.clear();
            charges.push_back(q);
            charges.push_back(-q);
        }

        // distance from center to charge, 7 point model
        void quadrupole_as_point_charges(
            const Matrix3d& quadrupole,
            double distance,
            std::vector<Vector3d>& r,
            std::vector<double>& charges)
        {
            Vector3d eigVec1, eigVec2, eigVec3, q_diag;
            double d = distance;

            algebra3d::eigensystemRealSymm(quadrupole, eigVec1, eigVec2, eigVec3, q_diag.x, q_diag.y, q_diag.z);
            Vector3d v1, v2;
            r.clear();
            charges.clear();

            v1.set(2, -1, -1);
            v1 /= sqrt(6.0);
            v2.set(0, 1, -1);
            v2 /= sqrt(2.0);

            double c1, c2, q1, q2;

            c1 = v1 * q_diag;
            c2 = v2 * q_diag;

            Vector3d q_diff = q_diag - c1 * v1 - c2 * v2;

            q1 = c1 / (d * d * 2 * sqrt(6.0));
            q2 = c2 / (d * d * 6 * sqrt(2.0));

            r.push_back(d * eigVec1);
            r.push_back(-d * eigVec1);

            charges.push_back(q1);
            charges.push_back(q1);


            r.push_back(d * eigVec2);
            r.push_back(-d * eigVec2);
            charges.push_back(q2);
            charges.push_back(q2);
            r.push_back(d * eigVec3);
            r.push_back(-d * eigVec3);
            charges.push_back(-q2);
            charges.push_back(-q2);

            r.push_back({ 0,0,0 });
            charges.push_back(-2 * q1);

        }


        void multipoles_as_point_charges(
            const ElectricMultipoles& multipoles,
            double distance,
            std::vector<Vector3d>& r,
            std::vector<double>& charges,
            bool pack,
            double packThreshold)
        {
            r.clear();
            charges.clear();

            std::vector<double> q;
            std::vector<Vector3d> position;

            if (multipoles.dipole != Vector3d(0, 0, 0))
                dipole_as_point_charges(multipoles.dipole, distance, r, charges);

            if (multipoles.charge != 0)
            {
                r.push_back(Vector3d(0, 0, 0));
                charges.push_back(multipoles.charge);
            }

            if (!multipoles.quadrupole.isZero())
            {
                quadrupole_as_point_charges(multipoles.quadrupole, distance, position, q);
                r.insert(r.end(), position.begin(), position.end());
                charges.insert(charges.end(), q.begin(), q.end());
            }

            if (pack)
            {
                packPointCharges(r, charges, position, q, packThreshold);
                charges.swap(q);
                r.swap(position);
            }

        }

        void point_charges_cluster(
            const Crystal& crystal,
            const std::vector<std::pair<std::string, std::string> > _clusterAtoms,
            const std::vector<double>& clusterAtomWeigth,
            const std::vector< std::vector<Vector3d> >& asymmR,
            const std::vector < std::vector<double> >& asymmQ,
            std::vector<Vector3d>& r,
            std::vector<double>& q)
        {
            r.clear();
            q.clear();

            vector<pair<int, string> > clusterAtoms;
            crystal_structure_utilities::convertAtomList(crystal, _clusterAtoms, clusterAtoms);

            vector<Vector3d> asymmUnitAtomCartesian(crystal.atoms.size());
            for (int i = 0; i < crystal.atoms.size(); i++)
                crystal.unitCell.fractionalToCartesian(crystal.atoms[i].coordinates, asymmUnitAtomCartesian[i]);

            Matrix3d rotation;
            Vector3d translation;

            int idx = 0;

            for (auto& atom : clusterAtoms)
            {
                int atomIdx = atom.first;
                crystal_structure_utilities::symmetryOperationCartesian(atom.second, crystal.unitCell, rotation, translation);

                const vector<discamb::Vector3d>& qPositions = asymmR[atomIdx];
                const vector<double>& qValues = asymmQ[atomIdx];

                for (int i = 0; i < qValues.size(); i++)
                {
                    q.push_back(clusterAtomWeigth[idx] * qValues[i]);
                    Vector3d rRot = rotation * (qPositions[i] + asymmUnitAtomCartesian[atomIdx]) + translation;
                    r.push_back(rRot);
                }

                idx++;
            }

        }


    }
}
