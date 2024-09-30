#include "discamb/BasicUtilities/utilities.h"

using namespace std;

namespace discamb {

    //struct HardwareResources {
    //    int nCores;
    //    int memoryMB;
    void HardwareResources::set(const nlohmann::json& data)
    {
        *this = HardwareResources();
        nCores = data.value("n cores", nCores);
        nCores = data.value("n threads", nCores);

        if (data.find("memory") != data.end())
        {
            if (data.find("memory")->is_string())
            {
                string memoryStr(data["memory"].get<string>());
                totalMemoryMB = stoul(memoryStr);
                if (memoryStr.find('G') != string::npos)
                    totalMemoryMB *= 1000;
                if (memoryStr.find('T') != string::npos)
                    totalMemoryMB *= 1000000;
            }
            else
                totalMemoryMB = data["memory"].get<int>();
        }
        else
            totalMemoryMB = 1000;


        if (data.find("disk") != data.end())
        {
            if (data.find("disk")->is_string())
            {
                string diskStr(data["disk"].get<string>());
                diskSpaceGB = stoul(diskStr);
                if (diskStr.find('M') != string::npos)
                    totalMemoryMB /= 1000;
                if (diskStr.find('T') != string::npos)
                    totalMemoryMB *= 1000;
            }
            else
                diskSpaceGB = data["disk"].get<int>();
        }
        else
            diskSpaceGB = 0;

    }
    
}

