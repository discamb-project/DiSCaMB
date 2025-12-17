#pragma once

#include "json.hpp"

#include <algorithm>
#include <complex>

namespace discamb {

    /**
* \addtogroup BasicUtilities BasicUtilities
* @{
*/


    struct HardwareResources {
        int nCores = 1;
        int totalMemoryMB = 1000;
        int diskSpaceGB = 10;
        void set(const nlohmann::json& data, bool setMissingToDefault = false);
    };


    namespace utilities {

        template<typename T>
        bool is_complex_type(const T&)
        {
            return false;
        }

        template<typename T>
        bool is_complex_type(const std::complex<T> &t)
        {
            return true;
        }

        template<typename T, typename U>
        bool hasElement(
            const T& container,
            const U& value)
        {
            auto it = std::find(container.begin(), container.end(), value);
            return (it != container.end());
        }


        template<typename T, typename U>
        bool hasElement(
            const T& container,
            const U& value,
            int& idx)
        {
            auto it = std::find(container.begin(), container.end(), value);
            if (it != container.end())
            {
                idx = std::distance(container.begin(), it);
                return true;
            }
            return false;
        }

        // -1 if the value is not present in the container
        template<typename T, typename U>
        int index(const T& container, const U& value)
        {
            int idx;
            if (hasElement(container, value, idx))
                return int(idx);
            return -1;
        }



    }
    /**@}*/
}

