#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "discamb/MathUtilities/molecular_radial_grids.h"

#include <cmath>

using namespace std;

namespace {
    //                                 0     1     2     3     4     5     6     7     8     9
    double xi_Treutler_Ahlrichs[] = { 0.0 , 0.8 , 0.9 , 1.8 , 1.4 , 1.3 , 1.1 , 0.9 , 0.9 , 0.9 , 
                                      0.9 , 1.4 , 1.3 , 1.3 , 1.2 , 1.1 , 1.0 , 1.0 , 1.0 , 1.5 ,
                                      1.4 , 1.3 , 1.2 , 1.2 , 1.2 , 1.2 , 1.2 , 1.2 , 1.1 , 1.1 ,
                                      1.1 , 1.1 , 1.0 , 0.9 , 0.9 , 0.9 , 0.9 };
}

namespace discamb {
    namespace molecular_radial_grids {


        void getRadialGrid(
            int atomic_number,
            int grid_size,
            vector<double> &points,
            vector<double> &weights)
        {
            vector<double> x_gc;
            vector<double> weights_gc;
            const double  bohr_over_angstrom = 0.52917721092;
            int i;
            double xi;
            double x, dr_dx;


            points.resize(grid_size);
            weights.resize(grid_size);
            getGaussChebyshevGrid(grid_size, x_gc, weights_gc);

            xi = xi_Treutler_Ahlrichs[atomic_number];//*bohr_over_angstrom;

            for (i = 0; i < grid_size; i++)
            {

                x = x_gc[i];
                // formula M4 from Treutler & Ahlrichs
                points[i] = (xi / log(2.0)) * pow(1.0 + x, 0.6) * log(2.0 / (1.0 - x));
                dr_dx = points[i] * (0.6 / (1 + x) + 1.0 / ((1 - x)*log(2.0 / (1.0 - x))));
                weights[i] = weights_gc[i] * dr_dx;
            }

        }



        void getGaussChebyshevGrid(
            int n,
            std::vector<double> &r,
            std::vector<double> &weights)
        {
            int i;
            double x, p;
            r.resize(n);
            weights.resize(n);

            for (i = 1; i <= n; i++)
            {
                p = M_PI / (n + 1);
                x = cos(i*p);
                weights[i - 1] = 16.0 / 3.0*pow(sin(i*p), 4.0) / (n + 1.0);
                r[i - 1] = 1.0 + 2.0 / M_PI * ((1.0 + 2.0 / 3.0*(1 - x * x))*x*sqrt(1 - x * x) - acos(x));
            }
        }


    }
}