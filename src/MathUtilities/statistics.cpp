#include "discamb/MathUtilities/statistics.h"


#include <cmath>

using namespace std;

namespace discamb
{
    namespace statistics {

        double average(
            const vector<double> &v)
        {
            if (v.empty())
                return 0;

            double result = 0;
            for (auto x : v)
                result += x;
            return result / v.size();
        }
        
        double rss(
            const std::vector<double> &v,
            double average)
        {
            double diff, rss = 0;
            for (double x : v)
            {
                diff = x - average;
                rss += diff * diff;
            }
            return rss;
        }

        double esd(
            const vector<double> &v,
            double average)
        {
            int n = v.size();
            if (n < 2)
                return 0;

            double rss = statistics::rss(v, average);
            return sqrt(rss / (n - 1));
        }

        void average_and_esd(
            const std::vector<double> &v,
            double &average,
            double &esd)
        {
            average = statistics::average(v);
            esd = statistics::esd(v, average);
        }

        double sum_of_squares(
            const std::vector<double> &v)
        {
            double d = 0;
            for (double x : v)
                d += x * x;
            return d;
        }

        void average_and_rss(
            const std::vector<double> &v,
            double &average,
            double &rss)
        {
            average = statistics::average(v);
            rss = statistics::rss(v, average);
        }
    }

}
