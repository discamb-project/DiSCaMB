#include "discamb/BasicUtilities/Task.h"
#include "discamb/BasicUtilities/Sheduler.h"

namespace discamb {

Task::Task()
{
    n = 1;
    scheduler = nullptr;
}

Task::Task(const Task& t)
{
    name = t.name;
    n = t.n;
    scheduler = t.scheduler;
}

void Task::operator()()
{
    this->run();
    scheduler->onTaskFinished(n, name);
}

}

