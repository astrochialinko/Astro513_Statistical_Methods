/*
  Time-stamp: <STimer.cc on Monday, 11 October, 2021 at 21:42:41 MST (pinto)>

  Simple class to support timing code (serial version)

 */

#ifndef INCLUDE_STIMER
#define INCLUDE_STIMER

#include <sys/time.h>
#include <cassert>

class STimer {
public:
    STimer();
    ~STimer();
    void Start();
    void Stop();
    double Elapsed();
    void Clear();
    void CStart();
    void Report(const char * msg);

    int timeron;

    struct timeval tuse, tstart, timer;
};

STimer::STimer(void) {
    timeron = false;
    timerclear(&timer);
    timeron = 0;
}

STimer::~STimer() {
}

void STimer::Start() {
    assert(!timeron);
    assert( gettimeofday( &(tstart), (struct timezone *)NULL ) == 0 );
    timeron = 1;
}

void STimer::Stop() {
    assert( timeron );
    struct timeval dt;
    assert( gettimeofday( &(tuse), (struct timezone *)NULL ) == 0 );
    timersub(&(tuse), &(tstart), &dt);
    timeradd(&dt, &(timer), &(timer));
    timeron = 0;
}

void STimer::Clear() {
    assert(!timeron);
    timerclear(&(timer));
}

void STimer::CStart() {
    Clear(); Start();
}

double STimer::Elapsed() {
    return  timer.tv_sec + 1e-6*timer.tv_usec;
}

void STimer::Report(const char * msg) {
    printf("%s %e\n", msg, Elapsed() );
}

#endif // INCLUDE_STIMER
