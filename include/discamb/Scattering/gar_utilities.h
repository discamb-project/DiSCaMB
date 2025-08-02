#pragma once

#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/Scattering/AtomRepresentativeInfo.h"

#include <optional>
#include <set>


namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    namespace gar_utilities {

        void findDeaultSubsystems(
            const Crystal& crystal,
            std::vector<std::vector<std::pair<std::string, std::string> > >& clusterAtoms,
            std::vector<std::string>& clustersLabels,
            std::vector<int>& clustersCharges,
            std::vector<int>& clustersSpinMultiplicity,
            int spinMultiplicityHint = 0);

        // false if not all atoms in crysta have representatives
        bool findDefaultRepresentatives(
            const Crystal& crystal,
            const std::vector<std::vector<std::pair<std::string, std::string> > >& subsystemAtoms,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives,
            std::vector<std::vector<std::pair<std::string, std::string> > >& representativesPerSubsystem,
            bool throw_if_false = true);

        // false if not all atoms in crysta have representatives
        bool findDefaultRepresentatives(
            const Crystal& crystal,
            const std::vector<std::vector<std::pair<std::string, std::string> > >& subsystemAtoms,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives,
            bool throw_if_false = true);


        void elementalIntegrationGrids(
            const std::set<int>& atomicNumbers,
            std::vector<std::vector<Vector3d> >& gridPoints,
            std::vector<std::vector<double> >& weights);

        /*for one fragment*/
        //void findMultipoleClasterAtoms(
        //    const std::optional<DistributedMultipoleCentersSettings>& settings,
        //    double threshold,
        //    const UnitCellContent& ucContent,
        //    const std::vector<std::vector<UnitCellContent::AtomID> >& unitCellMolecules,
        //    const FragmentAtoms& fragmentAtoms,
        //    std::vector<std::pair<std::string, std::string > >& multipoleClusterAtoms,
        //    std::vector<std::optional<double> >& multipoleCustomWeight);

        /*for multiple fragments*/
        //void findMultipoleClasterAtoms(
        //    const std::vector<std::optional<DistributedMultipoleCentersSettings> >& settings,
        //    double threshold,
        //    const Crystal &crystal,
        //    const std::vector<FragmentAtoms>& fragmentAtoms,
        //    std::vector<std::pair<std::string, std::string > >& multipoleClusterAtoms,
        //    std::vector<std::optional<double> >& multipoleCustomWeight);

    }
    /** @}*/
}
