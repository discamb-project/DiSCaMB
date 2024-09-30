#include "discamb/BasicUtilities/Sheduler.h"

#include <thread>
#include <iostream>

using namespace std;

namespace discamb {

    void Scheduler::run()
    {
        //while (runNextTask())
        //{
        //}
        int i = 0;
        while (true)
        {
            runNextTask();
            this_thread::sleep_for(chrono::milliseconds(250));
            i++;
            //cout << i * 0.25 << "\n";
            if (nToFinish == 0)
                break;
        }

    }

    bool Scheduler::runNextTask()
    {
        
        bool run = false;
        for (auto& item : tasks)
            if (get<2>(item))
                if (get<0>(item) <= nFreeCPU)
                {

                    nFreeCPU -= get<0>(item);
                    get<2>(item) = false;
                    get<3>(item)->n = get<0>(item);
                    get<3>(item)->name = get<1>(item);
                    get<3>(item)->scheduler = this;
                    cout << get<1>(item) << endl;
                    thread t(&Task::operator(), get<3>(item));
                    t.detach();
                    run = true;
                    break;
                }
        if (!run)
            return false;
        return true;
    }


    void Scheduler::onTaskFinished(int nReleasedCPU, string taskId)
    {
        for (int i = 0; i < tasks.size(); i++)
            if (get<1>(tasks[i]) == taskId)
                get<2>(tasks[i]) = false;
        nFreeCPU += nReleasedCPU;
        nToFinish--;
        
    }


    //void Scheduler::onTaskFinished(int nReleasedCPU, string taskId)
    //{
    //    for (int i = 0; i < tasks.size(); i++)
    //        if (get<1>(tasks[i]) == taskId)
    //            get<2>(tasks[i]) = false;
    //    nFreeCPU += nReleasedCPU;
    //    nToFinish--;
    //    run();
    //}


}
