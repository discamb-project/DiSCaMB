#include "discamb/MathUtilities/radial_grid.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cstddef>
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
    namespace radial_grid{
        void treutler_ahlrichs(
            int atomic_number,
            int grid_size,
            vector<double>& points,
            vector<double>& weights)
        {
            vector<double> x_gc;
            vector<double> weights_gc;
            //const double  bohr_over_angstrom = 0.52917721092;
            int i;
            double xi;
            double x, dr_dx;

            points.resize(grid_size);
            weights.resize(grid_size);
            perez_jorda(grid_size, x_gc, weights_gc);
            if (atomic_number <= 36)
                xi = xi_Treutler_Ahlrichs[atomic_number];//*bohr_over_angstrom;
            else
                xi = 1.0;

            for (i = 1; i <= grid_size; i++)
            {

                x = x_gc[i - 1];
                // formula M4 from Treutler & Ahlrichs
                points[i - 1] = (xi / log(2.0)) * pow(1.0 + x, 0.6) * log(2.0 / (1.0 - x));
                dr_dx = points[i - 1] * (0.6 / (1 + x) + 1.0 / ((1 - x) * log(2.0 / (1.0 - x))));
                weights[i - 1] = weights_gc[i - 1] * dr_dx;
            }

        }

        void mura_knowles(
            int atomic_number,
            int grid_size,
            vector<double>& points,
            vector<double>& weights)
        {
            //on_error::throwException("Mura-Knowles integratin scheme does not work\n", __FILE__, __LINE__);
            int i;
            double x;

            points.resize(grid_size);
            weights.resize(grid_size);

            double R3, R=5;
            int group = periodic_table::group(atomic_number);
            if (group == 1 || group == 2)
                R = 7.0;

            R3 = R * R * R;

            double x3, v;

            for (i = 1; i <= grid_size; i++)
            {
                
                x = double(i) / (grid_size + 1.0);
                
                x3 = x * x * x;
                
                points[i - 1] = -R * log(1 - x3);

                v = x * log(1 - x3);
                //weights[i - 1] = 3 * v * v * R3 / ((grid_size + 1) * (1 - x3));
                //weights[i - 1] = 3 * x*x*log(1-x3)*log(1-x3) * R3 / ((grid_size + 1) * (1 - x3));
                weights[i - 1] = 3 * x * x * R / ((grid_size + 1) * (1 - x3));
                
            }

        }


        void becke(
            int atomic_number,
            int grid_size,
            vector<double>& points,
            vector<double>& weights)
        {
            int i;
            double x;

            points.resize(grid_size);
            weights.resize(grid_size);

            double R = 3.0;

            for (i = 1; i <= grid_size; i++)
            {

                x = cos(i*M_PI/(grid_size+1.0));
                points[i - 1] = (1.0 + x) / (1.0 - x) * R;

                //weights[i - 1] = 2*M_PI/(grid_size+1.0) * (1.0+x)* (1.0 + x) * sqrt(1.0 + x) / 
                //                 ((1.0 - x) * (1.0 - x) * (1.0 - x) * sqrt(1.0 - x)) * R*R*R;

                weights[i - 1] = 2*M_PI/(grid_size+1.0) * (1.0+x)* sqrt(1.0 + x) / 
                                 ((1.0 - x) * (1.0 - x) * sqrt(1.0 - x)) * R*R*R;

            }

        }


        void perez_jorda(
            int n,
            std::vector<double>& r,
            std::vector<double>& weights)
        {
            int i;
            double x, p;
            r.resize(n);
            weights.resize(n);

            for (i = 1; i <= n; i++)
            {
                p = M_PI / (n + 1);
                x = cos(i * p);
                weights[i - 1] = 16.0 / 3.0 * pow(sin(i * p), 4.0) / (n + 1.0);
                r[i - 1] = 1.0 + 2.0 / M_PI * ((1.0 + 2.0 / 3.0 * (1 - x * x)) * x * sqrt(1 - x * x) - acos(x));
            }
        }

}
}
