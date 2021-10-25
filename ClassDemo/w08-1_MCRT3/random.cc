#ifndef INCLUDE_MYRNG
#define INCLUDE_MYRNG

#include <vector>
#include <random>
#include <omp.h>

class MyRNG {
public:

    MyRNG()  {
        // Create max_threads instances of RNG
        // seed them identically, and then discard i*sep numbers to get max_threads independent sequences
        nthreads = omp_get_max_threads();
    }

    ~MyRNG() {
        if(init) cleanUp();
    }

    void cleanUp() {
        for(int i=0; i<nthreads; i++) {
            delete rng[i];
            delete gen[i];
        }
        delete[] rng;
        delete[] gen;
    }

    void seed(int seed) {

        if(init) cleanUp();

        init = true;
        rng = new std::mt19937_64*[nthreads];
        gen = new std::uniform_real_distribution<double>*[nthreads];
        for(int i=0; i<nthreads; i++) {
            rng[i] = new std::mt19937_64();
            gen[i] = new std::uniform_real_distribution<double>(0.0,1.0);
            rng[i]->seed(seed*i);  // Bad! Bad! Bad! Bad! Bad!
        }
    }

    double operator() () {
        int thread = omp_get_thread_num();
        return gen[thread]->operator()(*rng[thread]);
    }

    std::vector<double> vec(int n) {
        int thread = omp_get_thread_num();
        std::vector<double> result(n);
        for( auto &r : result ) r = gen[thread]->operator()(*rng[thread]);
        return result;
    }

    bool init=false;
    int nthreads;
    std::mt19937_64 **rng;
    std::uniform_real_distribution<double> **gen;
};

#if 0

#include <cstdio>
#include "misc.cc"
#include "tallys.cc"

int main() {
    int nthreads = omp_get_max_threads();

    int sep = 1;
    int N = 100000;

    MyRNG xi;
    xi.seed(98123927);

    int ne = 121;
    TallyGrid foo(-0.1, 1.1, ne);

#pragma omp parallel for
    for (int i=0; i<N; i++) {
        foo.tally(1, xi());
    }

    foo.finish();

    for(int i=0; i<ne; i++) {
        printf("%e %e %e\n", foo.xgrid[i], foo.mean[i], foo.sig[i]);
    }

}
#endif

#endif  // INCLUDE_STIMER
