#pragma once

#include <vector>

namespace discamb{

    /**
    
    the data imported from Molden2AIM https://github.com/zorkzou/Molden2AIM 
    (Molden2AIM/src/edflib-pbe0.f90)

    ! EDFLIB: EDF data library by PBE0 (Ver.1, 2019.05.26)
    !
    ! Authors: W. Zou and C. Gao
    
    see also  https://doi.org/10.1002/jcc.25214

    */
    namespace ecp_electron_density_tables{
        void ecp_electron_density(
            int atomic_number,
            int n_core_electrons,
            std::vector<double> &exponents,
            std::vector<double> &coefficients);
    }

}

