#pragma once

#include <string>

namespace discamb {

    struct Scheduler;

    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


    struct Task {
        std::string name;
        int n;
        Scheduler* scheduler;
        //virtual void operator()() = 0;
        void operator()();
        virtual void run() = 0;
        Task(const Task& t);
        Task();
    };
    /**@}*/
}