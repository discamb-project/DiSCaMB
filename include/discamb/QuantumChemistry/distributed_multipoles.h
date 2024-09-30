#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/StructuralProperties/atom_selection.h"

#include "json.hpp"

namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    struct DistributedMultipoleCentersSettings {
        double multipoleClusterThreshold = 8.0;
        // if predefinedMultipoleClusterAtoms is not empty no atoms will be generated 
        // using threshold, instead 
        std::vector<std::pair<std::string, std::string> > predefinedMultipoleClusterAtoms;
        std::vector<AtomSubsetSelectionData> atomsToOmit;
        std::vector<std::pair<AtomSubsetSelectionData, double> > atomsWithCustomWeights;
        void set(nlohmann::json const& jsonData);
    };


    struct ElectricMultipoles {
        double charge = 0.0;
        Vector3d dipole;
        Matrix3d quadrupole;
    };


    struct DistributedMultipolesSet {
        std::vector<Vector3d> positions;// Angstroms
        std::vector<ElectricMultipoles> multipoles;
    };


    namespace distributed_multipoles {

        /** for one fragment*/
        void get_fragment_multipole_centers(
            const DistributedMultipoleCentersSettings& settings,
            const UnitCellContent& ucContent,
            const std::vector<std::vector<UnitCellContent::AtomID> >& unitCellMolecules,
            const FragmentAtoms& fragmentAtoms,
            std::vector<std::pair<std::string, std::string> >& multipoleCenters,
            std::vector<std::optional<double> >& weights);

        /** for multiple fragments*/
        void get_fragments_multipole_centers(
            const std::vector<std::optional<DistributedMultipoleCentersSettings> >& settings,
            double threshold,
            const Crystal& crystal,
            const std::vector<FragmentAtoms>& fragmentAtoms,
            std::vector< std::vector<std::pair<std::string, std::string > > >& multipoleClusterAtoms,
            std::vector < std::vector<std::optional<double> > >& multipoleCustomWeight);


        void calculate_dipole(const std::vector<Vector3d>& r, const std::vector<double>& charges, Vector3d& dipole);

        void calculate_quadrupole(const std::vector<Vector3d>& r, const std::vector<double>& charges, Matrix3d& quadrupole);

        void calculate_multipoles(const std::vector<Vector3d>& r, const std::vector<double>& charges, ElectricMultipoles& multipoles, int max_l);
        // distance from center to charge
        void dipole_as_point_charges(const Vector3d& dipole, double distance, std::vector<Vector3d>& r, std::vector<double>& charges);
        // distance from center to charge, 7 point model

        void quadrupole_as_point_charges(const Matrix3d& quadrupole, double distance, std::vector<Vector3d>& r, std::vector<double>& charges);

        void multipoles_as_point_charges(
            const ElectricMultipoles& multipoles,
            double distance,
            std::vector<Vector3d>& r,
            std::vector<double>& charges,
            bool pack = true,
            double = 1e-6);

        /*
        calculates cluster of point charges (points in r, charges in q)
        representing point multipoles centered at atoms given in clusterAtoms
        charges and positions of the point charges for each atom in asymmetric
        unit have to given in asymmQ and asymmR
        */
        void point_charges_cluster(
            const Crystal& crystal,
            const std::vector<std::pair<std::string, std::string> > clusterAtoms,
            const std::vector<double>& clusterAtomWeigth,
            const std::vector< std::vector<Vector3d> >& asymmR,
            const std::vector < std::vector<double> >& asymmQ,
            std::vector<Vector3d>& r,
            std::vector<double>& q);
    }

    /**@}*/

}

