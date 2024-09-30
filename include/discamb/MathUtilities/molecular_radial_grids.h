#pragma once

#include <vector>
#include <stddef.h>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    namespace molecular_radial_grids {
        /**
        Oliver Treutler and Reinhart Ahlrichs
        Efficient molecular numerical integration schemes
        J. Chem. Phys. 102 (1), 1 January 1995 p. 346
        */

        void getRadialGrid(int atomic_number, int grid_size, std::vector<double> &points, std::vector<double> &weights);

        /**
        as described in:
        J.M. Perez-Jorda, E. San-Fabian, F. Moscardo, A simple, reliable and efficient scheme for automatic numerical
        integration, Comput. Phys. Commun. 70(2) (1992) 271-284.
        (also mentioned in Journal of Computational and Applied Mathematics 112 (1999) 243-251)
        */

        void getGaussChebyshevGrid(int n, std::vector<double> &r, std::vector<double> &weights);
    }
    /** @}*/
}