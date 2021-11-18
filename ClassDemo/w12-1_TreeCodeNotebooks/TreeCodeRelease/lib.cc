/*
Time-stamp: <lib.cc on Tuesday, 16 November, 2021 at 09:19:55 MST (pinto)>

extern "C" inteface functions for the Python library.

*/

//#include <gperftools/profiler.h>

#include <cassert>
#include <omp.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "pprint.cc"
#include "smallvec.cc"
#include "pprint.cc"
#include "util.cc"
#include "STimer.cc"
#include "permute.cc"

#include "ptype.cc"
#include "particleset.cc"
#include "octree.cc"
#include "avxdirect.cc"
#include "bhtree.cc"

#include <cstring>

extern "C" {

    ParticleSet<double> *newParticleSetDouble() {
        return new ParticleSet<double>(); // create without allocating arrays
    }

    const char *getTypeString(ParticleSet<double> *me) {
        std::string s = me->getTypeString();
        return strdup(s.c_str());
    }

    void setPSDStorage(ParticleSet<double> *me, double *r, double *v, double *a,
                       double *pot, double *mass, int *id, int N) {

        std::get<0>(me->base) = reinterpret_cast<SmallVec<double, 3> *>(r);
        std::get<1>(me->base) = reinterpret_cast<SmallVec<double, 3> *>(v);
        std::get<2>(me->base) = reinterpret_cast<SmallVec<double, 3> *>(a);
        std::get<3>(me->base) = pot;
        std::get<4>(me->base) = mass;
        std::get<5>(me->base) = id;
        me->N = N;
    }

    void freeParticleSet(ParticleSet<double> *me) {
        me->N = 0; // prevent freeing Python's memory
        delete me;
    }

    void computeBoundingBox(ParticleSet<double> *me) {
        me->computeBoundingBox();
    }

    void getBoundingBox(ParticleSet<double> *me, double *centre, double *halfWidth) {
        for(int d=0; d<3; d++) centre[d] = me->centre[d];
        *halfWidth = me->halfWidth;
    }

    void read(ParticleSet<double> *me, char *fname) {
        me->read(std::string(fname), false);
    }

    void write(ParticleSet<double> *me, char *fname) {
        me->write(std::string(fname));
    }

    BHtree<double> *newBHTreeDouble(ParticleSet<double> *pset, int maxLeafSize, double epsSmooth, int maxSources) {
        return new BHtree<double>( *pset, maxLeafSize, epsSmooth, maxSources);
    }

    void BHfree(BHtree<double> *ptr) {
        delete ptr;
    }

    void makeTree(BHtree<double> *me, bool verbose, bool check, bool getStats) {
        me->verbose = verbose;
        me->check = check;
        me->getStats = getStats;
        me->makeTree();
    }

    void rspAcc(BHtree<double> *me, double theta, int maxILlength) {
        SmallVec<double, 3> delta(0.0);
        me->resetStats();
        me->recursiveSubsetPartition(theta, maxILlength, delta);
        me->reportStats();
    }

    void BHsubsets(BHtree<double> *me, double theta, int maxILlength) {
        SmallVec<double, 3> delta(0.0);
        me->resetStats();
        me->BHsubsets(theta, maxILlength, delta);
        me->reportStats();
    }

    void accAll(BHtree<double> *me, double theta) {
        SmallVec<double, 3> delta(0.0);
        me->resetStats();
        me->accAll(theta, delta);
        me->reportStats();
    }

    void putParticlesInOriginalOrder(BHtree<double> *me) {
        me->putParticlesInOriginalOrder();
    }

    /*
      Allows one to register a print function which will print in a Jupyter Notebook.
     */
    void registerPrint(BHtree<double> *me, void (*print_callback)(const char *s) ) {
        me->registerPrint(print_callback);
    }

}
