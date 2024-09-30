#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/Scattering/AtomRepresentativeInfo.h"

#include <iosfwd>

namespace discamb {
    /**
    * \addtogroup IO IO
    * @{
    */

    namespace discamb_io {

        void read_qm_systems(
            std::istream& in,
            const Crystal& crystal,
            std::vector<int>& charge,
            std::vector<int>& spinMultiplicity,
            std::vector<std::string>& systemLabels,
            std::vector< std::vector<std::pair<std::string, std::string> > >& clusters,
            std::vector<std::map<int, std::string> >& atomIdx2BasisSetMap);

        void read_qm_systems_HirshFrag(
            const std::string &file_name,
            std::vector<int>& charge,
            std::vector<int>& spinMultiplicity,
            std::vector<std::string>& systemLabels,
            std::vector< std::vector<std::pair<std::string, std::string> > >& clusters);

        void read_representatives(
            const std::string& fileName,
            const Crystal& crystal,
            const std::vector<std::string>& subsystemLabels,
            const std::vector< std::vector<std::pair<std::string, std::string> > >& subsystemAtoms,
            std::vector<std::vector<AtomRepresentativeInfo> >& representatives);


    }
    /**@}*/
}

