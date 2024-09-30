#pragma once

#include <vector>

namespace discamb
{

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    namespace statistics {
        /**
        returns 0 i v is empty
        */
        double average(const std::vector<double> &v);

        double sum_of_squares(const std::vector<double> &v);

        /**
        residual sum of squares, returns 0 if v is empty
        */
        double rss(const std::vector<double> &v, double average);
        /**
        estimated standard deviation, returns 0 if size of v is less than 2
        */
        double esd(const std::vector<double> &v, double average);
        void average_and_esd(const std::vector<double> &v, double &average, double &esd);
        void average_and_rss(const std::vector<double> &v, double &average, double &rss);
        
    }
    /** @}*/
}
