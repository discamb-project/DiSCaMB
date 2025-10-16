
#ifndef _DISCAMB_BASICUTILITIES_TIMER_H_
#define _DISCAMB_BASICUTILITIES_TIMER_H_

#include "discamb/config.h"
#include <string>

//unistd.h


#include <chrono>
#include <time.h>

#ifdef _WIN32
    #define NOMINMAX
    #include <Windows.h>
#endif


namespace discamb {

    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


/**
\brief Measure wall clock time.

 Implementation depends on the environment:
If USE_CPP11 is defined then uses std::chrono::steady_clock, else if HAS_SYS_TIME and HAS_UNISTD are defined the uses
POSIX gettimeofday, else if _WIN32 is defined Windows QueryPerformanceFrequency is used, else C time_t time(time_t *) is used


*/

class WallClockTimer
{
public:
    WallClockTimer();
    ~WallClockTimer();
    /**\brief Starts time measurement.*/
    void start();
    /**\brief Stops time measurement and returns time in ms.*/
    double stop();
    double elapsedTime() const;
    /** \brief Returns type of the clock.
        Result is one of the following: std::chrono:: 'steady_clock', 'POSIX <sys/time.h> gettimeofday',
        'Windows QueryPerformanceFrequency', 'C time_t time(time_t *)'
 */
    std::string type() const;

private:

    std::chrono::steady_clock::time_point mStart, mStop;
};


/** \brief Measures CPU time. 
Implementation depends on environment: if _WIN32 is defined GetProcessTimes is used otherwise \<ctime\> clock() is used.
*/

class CpuTimer
{
public:
    CpuTimer();
    ~CpuTimer();
    /** \brief Starts execution time measurement.*/
    void start();
    /** \brief Returns execution time measurement.*/
    double stop(); 

    /** \brief Returns type of the clock.
    Result is one of the following: std::chrono:: 'high_resolution_clock', 'POSIX <sys/time.h> gettimeofday',
    'Windows QueryPerformanceFrequency', 'C time_t time(time_t *)'
    */

    std::string type() const;

private:
#ifdef _WIN32
    FILETIME mStart, mStop, mAux1, mAux2, mAux3;
#else
    clock_t mStart;
#endif
};

/**@}*/

}// namespace discamb

#endif /*_DISCAMB_BASICUTILITIES_TIMER_H_*/

