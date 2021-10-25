#ifndef INCLUDE_TALLYS
#define INCLUDE_TALLYS

#include <cmath>
#include <cstdio>
#include <vector>
#include <omp.h>

class TallyBase {
public:

    typedef std::vector<double> dVec;
    typedef std::vector<uint64_t> uVec;
    typedef std::vector<std::vector<double>> ddVec;
    typedef std::vector<std::vector<uint64_t>> uuVec;


    TallyBase(int nx) : nx(nx) {
        nthreads = omp_get_max_threads();
        sum.resize(nthreads, dVec(nx, 0.0));
        sum2.resize(nthreads, dVec(nx, 0.0));
        cnt.resize(nthreads, uVec(nx, 0));
        btotal.resize(nthreads, 0.0);
    }

    void add(double q, double x, int i) {
        int thread = omp_get_thread_num();
        int ii = std::min(nx-1, std::max(0, i));
        sum[thread][ii] += q;
        sum2[thread][ii] += q*q;
        cnt[thread][ii] += 1;
        btotal[thread] += q*x;
    }

    template <class T>
    void sumThreads(std::vector<std::vector<T>> &vec) {
        for(int t=1; t<nthreads; t++) {
            for(int i=0; i<nx; i++) vec[0][i] += vec[t][i];
        }
    }

    template <class T>
    void sumThreads(std::vector<T> &vec) {
        for(int t=1; t<nthreads; t++) {
            vec[0] += vec[t];
        }
    }

    void done(double scale) {

        sumThreads(sum);
        sumThreads(sum2);
        sumThreads(cnt);
        sumThreads(btotal);

        nsamp = 0;
        for(int i=0; i<nx; i++) nsamp += cnt[0][i];

        mean.resize(nx);
        sig.resize(nx);

        for(int j=0; j<nx; j++) {
            if( cnt[0][j] > 0) {
                mean[j] = sum[0][j] / cnt[0][j];
                sig[j] = std::sqrt( (sum2[0][j]/cnt[0][j] - mean[j]*mean[j]) / cnt[0][j] );
            }
            else {
                mean[j] = 0;
                sig[j] = 0;
            }
            mean[j] *= nsamp*scale;
            sig[j] *=nsamp*scale;
        }

        total = btotal[0] * scale;
    }

    int nthreads, nx;
    double total;
    uint64_t nsamp;
    dVec xgrid, mean, sig, btotal;
    ddVec sum, sum2;
    uuVec cnt;
};

class TallyGrid : public TallyBase {
public:

    TallyGrid(double xMin, double xMax, int nx) : TallyBase(nx) {
        xgrid.resize(nx);
        double dx = (xMax-xMin)/(nx-1);
        for(int i=0; i<nx; i++) xgrid[i] = xMin + i*dx;
    }

    void tally(double q, double x) {
        double dx = xgrid[1]-xgrid[0];
        int i = bound(0, std::floor((x-xgrid[0])/dx), nx-1);
        add(q, x, i);
    }

    void finish(double scale) {
        done(scale);
    }

    dVec xgrid;
};


class TallyZone : public TallyBase {
public:

    TallyZone(int nx) : TallyBase(nx) {}

    void tally(double q, int i) {
        add(q, 1.0, i);
    }

    void finish(double scale) {
        done(scale);
    }

};

class TallyScalar : public TallyBase {
public:

    TallyScalar() : TallyBase(1) {}

    void tally(double q) {
        add(q, 1.0, 0);
    }

    void finish(double scale) {
        done(scale);
        mean = TallyBase::mean[0];
        sig = TallyBase::sig[0];
        sum = TallyBase::sum[0][0];
    }

    double mean, sig, sum;
};


#endif // INCLUDE_TALLYS
