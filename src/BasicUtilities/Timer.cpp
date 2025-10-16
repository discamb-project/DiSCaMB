#include "discamb/BasicUtilities/Timer.h"

using namespace std;

namespace discamb {

WallClockTimer::WallClockTimer()
{
}

WallClockTimer::~WallClockTimer()
{
}

void WallClockTimer::start()
{
    mStart = chrono::steady_clock::now();
}

// returns time in ms
double WallClockTimer::stop()
{
    mStop = chrono::steady_clock::now();
    return double(chrono::duration_cast<chrono::milliseconds>(mStop-mStart).count());
}

double WallClockTimer::elapsedTime()
const
{
    auto time_now = chrono::steady_clock::now();
    return double(chrono::duration_cast<chrono::milliseconds>(time_now - mStart).count());
}

std::string WallClockTimer::type()
const
{
    return string("std::chrono::steady_clock");
}




CpuTimer::CpuTimer()
{
}

CpuTimer::~CpuTimer()
{
}

void CpuTimer::start()
{
#ifdef _WIN32
    GetProcessTimes(GetCurrentProcess(), &mAux1, &mAux2, &mAux3, &mStart);
#else
    mStart = clock();
#endif

}

double CpuTimer::stop()
{
#ifdef _WIN32
    GetProcessTimes(GetCurrentProcess(), &mAux1, &mAux2, &mAux3, &mStop);
    LARGE_INTEGER stop, start;
    stop.HighPart = mStop.dwHighDateTime;
    stop.LowPart = mStop.dwLowDateTime;
    start.HighPart = mStart.dwHighDateTime;
    start.LowPart = mStart.dwLowDateTime;
    return (stop.QuadPart - start.QuadPart)*0.0001;
#else
    return double(clock() - mStart)*1000.0 / CLOCKS_PER_SEC;
#endif
}

std::string CpuTimer::type()
const
{
#ifdef _WIN32
    return string("Windows GetProcessTimes");
#else
    return string("std::clock");
#endif

}

}//namespace discamb
