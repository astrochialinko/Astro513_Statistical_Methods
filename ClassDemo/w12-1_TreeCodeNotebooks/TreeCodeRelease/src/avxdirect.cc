#ifndef __AVXDIRECT_CC__
#define __AVXDIRECT_CC__

/*

  NBody direct interactions using the AVX registers.

  Source points can be added to an interaction set using various means.
  The accelerations on a sink set can them be computed using the "compute" calls.

  See bhtree.cc for examples of use (in the traversals near the end of the file).

*/

#define NUM_PIPE 4
#define NEWTON_ITER

#include "smallvec.cc"
#include "pprint.cc"
#include "particleset.cc"

template <typename T>
class AVXDirect {
private:
    template <typename U>
    struct srcpoint {
        U x, y, z, m;
    };

    template <typename U, int N>
    struct sinkstruct {
        U x[N], y[N], z[N], eps2[N];
    };

    template <typename U, int N>
    struct accstruct {
        U x[N], y[N], z[N], pot[N];
    };

public:
    AVXDirect();
    ~AVXDirect();
    void initAVX(int maxsrc);
    void compute( SmallVec<T,3> *sinkr, T *sinkmass, int nsink, SmallVec<T, 3> &delta, T eps2, SmallVec<T,3> *acc, T *pot );
    void compute( ParticleSet<T> &ps, int sinkBegin, int sinkEnd, SmallVec<T, 3> &delta, T eps2);

    void resetSrcCount() {
        nsrc = 0;
    }

    template <typename U>
    inline void addSrcPoint(SmallVec<U, 3> &point, U &mass) {
        srcdata[nsrc++] = {point[0], point[1], point[2], mass};
    }

    template <typename U>
    inline void addSrcPoint(U &x, U &y, U &z, U &mass) {
        srcdata[nsrc++] = {x, y, z, mass};
    }

    void addSources( SmallVec<T, 3> *srcr, T *srcmass, int beg, int end ) {
        for(int i=beg; i<=end; i++) {
            srcdata[nsrc++] = {srcr[i][0], srcr[i][1], srcr[i][2], srcmass[i]};
        }
    }

    void addSources( ParticleSet<T> &ps, int beg, int end ) {
        for(int i=beg; i<=end; i++) {
            srcdata[nsrc++] = {ps.arrays().r[i][0], ps.arrays().r[i][1], ps.arrays().r[i][2], ps.arrays().mass[i]};
        }
    }

    int maxsrc;
    int nsrc;
    srcpoint<T> *srcdata;

private:
    void setconpot();

    void KernelAccPot(sinkstruct<T, NUM_PIPE> *sinkdata, srcpoint<T> *srcdata, int nsrc,
                      sinkstruct<T, NUM_PIPE> *deltas, accstruct<T, NUM_PIPE> *accdata);

    sinkstruct<T, NUM_PIPE> *sinkdata;
    sinkstruct<T, NUM_PIPE> *deltas;
    accstruct<T, NUM_PIPE> *accdata;
    T conpot;
};

template <typename T>
AVXDirect<T>::AVXDirect() {}

template <typename T>
void AVXDirect<T>::initAVX(int _maxsrc) {
    maxsrc = _maxsrc;
    srcdata  = (srcpoint<T> *)             aligned_alloc(128, sizeof(srcpoint<T>)*(maxsrc+4));  assert( srcdata != NULL );
    sinkdata = (sinkstruct<T, NUM_PIPE> *) aligned_alloc(128, sizeof(sinkstruct<T, NUM_PIPE>)); assert( srcdata != NULL );
    deltas   = (sinkstruct<T, NUM_PIPE> *) aligned_alloc(128, sizeof(sinkstruct<T, NUM_PIPE>)); assert( deltas  != NULL );
    accdata  = (accstruct<T, NUM_PIPE> *)  aligned_alloc(128, sizeof(accstruct<T, NUM_PIPE>));  assert( accdata != NULL );
    setconpot();
    resetSrcCount();
}

template <typename T>
AVXDirect<T>::~AVXDirect() {
    if(srcdata != NULL) {
        free(srcdata); free(sinkdata); free(deltas); free(accdata);
    }
}

template <typename T>
void AVXDirect<T>::setconpot() {
    pprint("AVXDirect::setconpot: type error\n");
    exit(1);
}

template <>
void AVXDirect<float>::setconpot() {
    conpot = -0.5;
}

template <>
void AVXDirect<double>::setconpot() {
    conpot = sqrt(2.0);
}

/*
  Apply a particleset src of length nsrc to a particleset sink of length nsink with coordinate systems separated
  by delta, to produce a particleset of acceleration and potential.
*/

template <typename T>
void AVXDirect<T>::compute( ParticleSet<T> &ps, int sinkBegin, int sinkEnd, SmallVec<T, 3> &delta, T eps2 ) {

#pragma omp critical
    if(nsrc > maxsrc) pprint("AVXDirect::compute: nsrc: %d > maxsrc: %d\n", nsrc, maxsrc);
    assert( nsrc <= maxsrc );

#ifndef NEWTON_ITER
    conpot = 1.0;
#endif
    T conacc = conpot*conpot*conpot;

    for(int i=0; i<NUM_PIPE; i++) {
        deltas->x[i] = (float)delta[0];
        deltas->y[i] = (float)delta[1];
        deltas->z[i] = (float)delta[2];
    }

    for(int p=sinkBegin; p<=sinkEnd; p+=NUM_PIPE) {
        //pprint("p: %d\n", p);
        int e = NUM_PIPE; if(p+e>sinkEnd) e = sinkEnd-p+1;
        //pprint("e: %d\n", e);
        for(int j=0; j<e; j++) {
            //pprint("j: %e\n", j);
            sinkdata->x[j] = ps.arrays().r[p+j][0];
            sinkdata->y[j] = ps.arrays().r[p+j][1];
            sinkdata->z[j] = ps.arrays().r[p+j][2];
            sinkdata->eps2[j] = eps2;
        }

        // apply all of the sources to these sinks
        KernelAccPot(sinkdata, srcdata, nsrc, deltas, accdata);

        for(int j=0; j<e; j++) {
            ps.arrays().a[p+j][0] = accdata->x[j]*conacc;
            ps.arrays().a[p+j][1] = accdata->y[j]*conacc;
            ps.arrays().a[p+j][2] = accdata->z[j]*conacc;
            ps.arrays().pot[p+j]  = accdata->pot[j]*conpot;
        }
    }
}

/*
  Apply a particleset src of length nsrc to a particleset sink of length nsink with coordinate systems separated
  by delta, to produce a particleset of acceleration and potential.
 */
template <typename T>
void AVXDirect<T>::compute(SmallVec<T,3> *sinkr, T *sinkmass, int nsink,
                           SmallVec<T, 3> &delta, T eps2,
                           SmallVec<T,3> *acc, T *pot) {

#pragma omp critical
    if(nsrc > maxsrc) pprint("AVXDirect::compute: nsrc: %d > maxsrc: %d\n", nsrc, maxsrc);
    assert( nsrc <= maxsrc );

#ifndef NEWTON_ITER
    conpot = 1.0;
#endif
    T conacc = conpot*conpot*conpot;

    for(int i=0; i<NUM_PIPE; i++) {
        deltas->x[i] = (float)delta[0];
        deltas->y[i] = (float)delta[1];
        deltas->z[i] = (float)delta[2];
    }

    for(int p=0; p<nsink; p+=NUM_PIPE) {

        int e = NUM_PIPE; if(p+e>nsink) e = nsink-p;
        for(int j=0; j<e; j++) {
            sinkdata->x[j] = sinkr[p+j][0];
            sinkdata->y[j] = sinkr[p+j][1];
            sinkdata->z[j] = sinkr[p+j][2];
            sinkdata->eps2[j] = eps2;
        }

        // applied all of the sources to these sinks
        KernelAccPot(sinkdata, srcdata, nsrc, deltas, accdata);

        for(int j=0; j<e; j++) {
            acc[p+j][0] = accdata->x[j]*conacc;
            acc[p+j][1] = accdata->y[j]*conacc;
            acc[p+j][2] = accdata->z[j]*conacc;
            pot[p+j]    = accdata->pot[j]*conpot;
        }
    }
}


#define AX    YMM08
#define AY    YMM09
#define AZ    YMM10
#define PHI   YMM11

#define DX    YMM12
#define DY    YMM13
#define DZ    YMM14
#define MJ    YMM07

#define J1    YMM00
#define J2    YMM01
#define X2    YMM02
#define xX2   XMM02
#define Y2    YMM03

#define XI    YMM04
#define YI    YMM05
#define ZI    YMM06
#define EPS2  YMM15

#include "avxsseabrev.h"
#define ALIGN64 __attribute__ ((aligned(64)))

// do one xj on four xi's per loop
template <>
void AVXDirect<double>::KernelAccPot(sinkstruct<double, NUM_PIPE> *sinkdata, srcpoint<double> *srcdata, int nsrc,
                                     sinkstruct<double, NUM_PIPE> *deltas,
                                     accstruct<double, NUM_PIPE> *accdata) {
    int j;

#ifdef NEWTON_ITER
    double threehalf ALIGN64;
    threehalf = 1.5;
#endif

    PREFETCH(*srcdata);

    VZEROALL; // to zero out acc registers for accumulation

#ifdef NEWTON_ITER
    VBROADCASTSD(threehalf, J1);
#endif

    VLOADPD(*sinkdata->x, XI);    // load 4 xi's
    VLOADPD(*sinkdata->y, YI);
    VLOADPD(*sinkdata->z, ZI);
    VLOADPD(*sinkdata->eps2, EPS2);

    VADDPD_M(*deltas->x, XI, XI); // add deltas to sink position
    VADDPD_M(*deltas->y, YI, YI);
    VADDPD_M(*deltas->z, ZI, ZI);

    VBROADCASTSD(srcdata->x, X2); // load 1 xj into four copies
    VBROADCASTSD(srcdata->y, Y2);
    VBROADCASTSD(srcdata->z, J2);
    VBROADCASTSD(srcdata->m, MJ);
    srcdata++;

    for(j=0; j<nsrc; j++) {

        PREFETCH(*srcdata);  // does this do anything -- unlikely!

        VSUBPD(XI, X2, DX);
        VSUBPD(ZI, J2, DZ);
        VSUBPD(YI, Y2, DY);

        VMULPD(DX, DX, X2);      // X2 = DX^2
        VMULPD(DZ, DZ, J2);
        VMULPD(DY, DY, Y2);

        VADDPD(X2, J2, J2);   // J2 = X2 + J2 = DX^2 + DZ^2
        VADDPD(EPS2, Y2, Y2); // Y2 = Y2 + EPS2 = DY^2 + eps^2
        VADDPD(J2, Y2, Y2);   // Y2 = Y2 + J2 = DX^2 + DY^2 + DZ^2 + eps^2 = R^2

#ifdef NEWTON_ITER
        VADDPD(Y2, Y2, X2);
        VCVTPD2PS(X2, xX2);   // convert to float to use rsqrt
#else
        VCVTPD2PS(Y2, xX2);
#endif

        VRSQRTPS(xX2, xX2);   // 1/sqrt(2*R^2)
        VCVTPS2PD(xX2, X2);

#ifdef NEWTON_ITER
        VMULPD(X2, X2, J2);   // J2 = x0^2
        VMULPD(J2, Y2, J2);   // J2 = x0^2 R^2
        VSUBPD(J2, J1, J2);   // J2 = 1.5 - x0^2 R^2
        VMULPD(J2, X2, X2);   // X2 = 1.5 x0 - x0^3 R^2

        VMULPD(X2, X2, J2);   // J2 = x0^2
        VMULPD(J2, Y2, J2);   // J2 = x0^2 R^2
        VSUBPD(J2, J1, J2);   // J2 = 1.5 - x0^2 R^2
        VMULPD(J2, X2, X2);   // X2 = 1.5 x0 - x0^3 R^2
#endif

        VMULPD(X2, MJ, MJ);   // MJ = MJ/R = m/R
        VSUBPD(MJ, PHI, PHI); // PHI = m/R - PHI

        VMULPD(X2, X2, Y2);   // Y2 = 1/R^2
        VMULPD(MJ, Y2, Y2);   // Y2 = m/R * 1/R^2 = m/R^3

        VMULPD(Y2, DX, DX);   // DX = DX m / R^3 = dx m / R^3
        VMULPD(Y2, DY, DY);
        VMULPD(Y2, DZ, DZ);

        VBROADCASTSD(srcdata->x, X2);
        VBROADCASTSD(srcdata->y, Y2);
        VBROADCASTSD(srcdata->z, J2);
        VBROADCASTSD(srcdata->m, MJ);
        srcdata++;

        VADDPD(DX, AX, AX);
        VADDPD(DY, AY, AY);
        VADDPD(DZ, AZ, AZ);
    }

    VSTORPD(AX, *accdata->x);
    VSTORPD(AY, *accdata->y);
    VSTORPD(AZ, *accdata->z);
    VSTORPD(PHI, *accdata->pot);

}

//#define NUNROLL 2 // do two xj's on four xi's per loop
#define NUNROLL 4 // do four xj's on four xi's per loop

template <>
void AVXDirect<float>::KernelAccPot(sinkstruct<float, NUM_PIPE> *sinkdata, srcpoint<float> *srcdata, int nsrc,
                                    sinkstruct<float, NUM_PIPE> *deltas,
                                    accstruct<float, NUM_PIPE> *accdata) {
    int j;

#ifdef NEWTON_ITER
    float three ALIGN64;
    three = -3.0;
#endif

    PREFETCH(srcdata[0]);

    VZEROALL; // to zero out acc registers for accumulation

    LOADPS(*sinkdata->x, XMM04);    // load 4 xi's into XMM04 (aliased into XI)
    LOADPS(*sinkdata->y, XMM05);
    LOADPS(*sinkdata->z, XMM06);
    LOADPS(*sinkdata->eps2, XMM15);

    VADDPS_M(*deltas->x, XMM04, XMM04); // add deltas to sink position
    VADDPS_M(*deltas->y, XMM05, XMM05);
    VADDPS_M(*deltas->z, XMM06, XMM06);

    VPERM2F128(XI, XI, XI, 0x00); // copy low 4 sp into high 4 sp for xi, yi, zi, and eps
    VPERM2F128(YI, YI, YI, 0x00);
    VPERM2F128(ZI, ZI, ZI, 0x00);
    VPERM2F128(EPS2, EPS2, EPS2, 0x00);
    // so that, e.g., XI = (xi0, xi1, xi2, xi3, xi0, xi1, xi2, xi3) and so on...

#if (2 == NUNROLL)  // 2 j's on 4 i's

    VLOADPS(*(srcdata), J1);       // load 2 (x,y,z,m)'s into J1
    srcdata += 2;

    VSHUFPS(J1, J1, X2, 0x00);    // X2 = (xj0,xj0,xj0,xj0, xj1,xj1,xj1,xj1)
    VSHUFPS(J1, J1, J2, 0xaa);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, Y2, 0x55);

    for(j = 0; j < nsrc; j += 2){

        VSUBPS(XI, X2, DX);       // DX = (xi0-xj0, xi1-xj0, xi2-xj0, xi3-xj0, xi0-xj1, xi1-xj1, xi2-xj1, xi3-xj1)
        VSUBPS(ZI, J2, DZ);
        VSUBPS(YI, Y2, DY);

        VLOADPS(*(srcdata), J1);
        srcdata += 2;

        VMULPS(DX, DX, X2);      // X2 = DX^2
        VMULPS(DZ, DZ, J2);
        VMULPS(DY, DY, Y2);

        VADDPS(X2, J2, J2);      // J2 = X2 + J2 = DX^2 + DZ^2
        VADDPS(EPS2, Y2, Y2);    // Y2 = Y2 + EPS2 = DY^2 + eps^2
        VADDPS(J2, Y2, Y2);      // Y2 = Y2 + J2 = DX^2 + DY^2 + DZ^2 + eps^2 = R^2

#ifdef NEWTON_ITER
        VBROADCASTSS(three, J2);
#endif

        VRSQRTPS(Y2, X2);

#ifdef NEWTON_ITER
        VMULPS(X2, Y2, Y2);
        VMULPS(X2, Y2, Y2);
        VADDPS(J2, Y2, Y2);
        VMULPS(Y2, X2, X2);
#endif

        VMULPS(X2, MJ, MJ);      // MJ = MJ/R = m/R
        VMULPS(X2, X2, Y2);      // Y2 = 1/R^2

        VMULPS(MJ, Y2, Y2);      // Y2 = m/R * 1/R^2 = m/R^3
        VSUBPS(MJ, PHI, PHI);    // PHI = m/R - PHI

        VMULPS(Y2, DX, DX);      // DX = DX m / R^3 = dx m / R^3
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VSHUFPS(J1, J1, X2, 0x00);
        VSHUFPS(J1, J1, J2, 0xaa);
        VSHUFPS(J1, J1, MJ, 0xff);
        VSHUFPS(J1, J1, Y2, 0x55);

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);
    }

#elif (4 == NUNROLL)  // 4 j's on 4 i's

    VLOADPS(*(srcdata), J1);
    VLOADPS(*(srcdata+2), J2);

    srcdata += 4;

    VSHUFPS(J1, J1, X2, 0x00);
    VSHUFPS(J1, J1, Y2, 0x55);
    VSHUFPS(J1, J1, MJ, 0xff);
    VSHUFPS(J1, J1, J1, 0xaa);

    for(j = 0 ; j < nsrc; j += 4) {

        VSUBPS(XI, X2, DX);
        VSUBPS(YI, Y2, DY);
        VSUBPS(ZI, J1, DZ);

        VMULPS(DX, DX, X2);
        VMULPS(DZ, DZ, J1);
        VMULPS(DY, DY, Y2);

        VADDPS(J1, X2, X2);
        VADDPS(EPS2, Y2, Y2);
        VADDPS(Y2, X2, Y2);

#ifdef NEWTON_ITER
        VBROADCASTSS(three, J1);
#endif
        VRSQRTPS(Y2, X2);

#ifdef NEWTON_ITER
        VMULPS(X2, Y2, Y2);
        VMULPS(X2, Y2, Y2);
        VADDPS(J1, Y2, Y2);
        VMULPS(Y2, X2, X2);
#endif
        VLOADPS(*(srcdata), J1);

        VMULPS(X2, MJ, MJ);
        VMULPS(X2, X2, Y2);

        VMULPS(MJ, Y2, Y2);
        VSUBPS(MJ, PHI, PHI);

        VMULPS(Y2, DX, DX);
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VSHUFPS(J2, J2, X2, 0x00);
        VSHUFPS(J2, J2, MJ, 0xff);
        VSHUFPS(J2, J2, Y2, 0x55);
        VSHUFPS(J2, J2, J2, 0xaa);

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);

        VSUBPS(XI, X2, DX);
        VSUBPS(YI, Y2, DY);
        VSUBPS(ZI, J2, DZ);

        VMULPS(DX, DX, X2);
        VMULPS(DZ, DZ, J2);
        VMULPS(DY, DY, Y2);

        VADDPS(J2, X2, X2);
        VADDPS(EPS2, Y2, Y2);
        VADDPS(Y2, X2, Y2);

#ifdef NEWTON_ITER
        VBROADCASTSS(three, J2);
#endif
        VRSQRTPS(Y2, X2);

#ifdef NEWTON_ITER
        VMULPS(X2, Y2, Y2);
        VMULPS(X2, Y2, Y2);
        VADDPS(J2, Y2, Y2);
        VMULPS(Y2, X2, X2);
#endif

        VLOADPS(*(srcdata+2), J2);

        VMULPS(X2, MJ, MJ);
        VMULPS(X2, X2, Y2);

        srcdata += 4;
        PREFETCH(*(srcdata));

        VMULPS(MJ, Y2, Y2);
        VSUBPS(MJ, PHI, PHI);

        VMULPS(Y2, DX, DX);
        VMULPS(Y2, DY, DY);
        VMULPS(Y2, DZ, DZ);

        VSHUFPS(J1, J1, X2, 0x00);
        VSHUFPS(J1, J1, MJ, 0xff);
        VSHUFPS(J1, J1, Y2, 0x55);
        VSHUFPS(J1, J1, J1, 0xaa);

        VADDPS(DX, AX, AX);
        VADDPS(DY, AY, AY);
        VADDPS(DZ, AZ, AZ);
    }
#else
#error "must define NUNROLL in AVXDirect<float>::KernelAccPot"
#endif

    VEXTRACTF128(AX, XMM00, 0x01);   // extract upper 128 from AX into XMM00
    VADDPS(AX, YMM00, AX);           // add to lower 128
    VEXTRACTF128(AY, XMM01, 0x01);
    VADDPS(AY, YMM01, AY);
    VEXTRACTF128(AZ, XMM02, 0x01);
    VADDPS(AZ, YMM02, AZ);
    VEXTRACTF128(PHI, XMM03, 0x01);
    VADDPS(PHI, YMM03, PHI);

    STORPS(XMM08, *accdata->x);
    STORPS(XMM09, *accdata->y);
    STORPS(XMM10, *accdata->z);
    STORPS(XMM11, *accdata->pot);

}
#undef AX
#undef AY
#undef AZ
#undef PHI
#undef DX
#undef DY
#undef DZ
#undef MJ
#undef J1
#undef J2
#undef X2
#undef Y2
#undef XI
#undef YI
#undef ZI
#undef EPS2
#undef NUNROLL

#endif // __AVXDIRECT_CC__
