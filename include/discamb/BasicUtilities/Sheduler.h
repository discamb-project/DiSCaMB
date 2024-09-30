#pragma once

#include "Task.h"

#include <vector>
#include <tuple>
#include <memory>
#include <atomic>

namespace discamb
{
    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


    struct Scheduler {
        void onTaskFinished(int nReleasedCPU, std::string);
        //n_cpu, name, to be run
        std::vector<std::tuple<int, std::string, bool, std::shared_ptr<Task> > > tasks;
        
        std::atomic_int nToFinish;
        std::atomic_int nFreeCPU;
        void run();
        bool runNextTask();
    };
    /**@}*/
}