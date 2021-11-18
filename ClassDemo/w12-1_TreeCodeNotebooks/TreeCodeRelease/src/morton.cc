#ifndef __MORTON_INCLUDE__
#define __MORTON_INCLUDE__

#include "inttype.h"
#include "smallvec.cc"
#include "pprint.cc"

//--------------------------------------------------------------------------------
// Some useful ways to print 126-bit Morton keys
//--------------------------------------------------------------------------------
#define OCTL(x, oct) (uint64)(( (x)>>(oct*3) )&7UL)
void out126Octal(uint128 r) {
    for(int l=41; l>=0; l--)printf("%lo", OCTL(r, l));
}
#undef OCTL

#define BITL(x, bit) (uint64)(( (x)>>(bit) )&1UL)
void out64Binary(uint64 r) {
    for(int l=63; l>=0; l--) printf("%lu", BITL(r,l));
}
#undef BITL

void out128Binary(uint128 r) {
    uint64 l,h;
    h = (r>>64);
    l = r;
    out64Binary(h);
    out64Binary(l);
}

void out128Hex(uint128 r) {
    uint64 l,h;
    h = (r>>64);
    l = r;
    printf("%016lx%016lx",h,l);
}

//--------------------------------------------------------------------------------
// Simple but slower Morton key class for testing
//--------------------------------------------------------------------------------
#if 1
template <int D, int keyLen>
class MortonBinary {
public:

    typedef SmallVec<double, D> dvec;

    MortonBinary(dvec &boxMin, dvec &boxMax) : boxMin(boxMin), boxMax(boxMax) {
        scale = 1.0/(boxMax-boxMin).maxcomponent();
    }

    uint128 makeKey(dvec v){
        v = scale*(v-boxMin);

        dvec c = {0.5,0.5,0.5};
        double r = 0.25;

        uint128 ONE = 0x1UL;

        uint32 offset = D*(keyLen-1);
        uint128 key = 0;
        for (int level=0; level<keyLen; level++) {
            uint128 oct = 0;
            for(int d=0; d<D; d++) {
                if( v[d] >= c[d] ) {
                    oct |= (ONE<<d);
                    c[d] += r;
                }
                else {
                    c[d] -= r;
                }
            }
            key |= (oct<<(offset-3*level));

            r *= 0.5;
        }
        return key;
    }

    dvec boxMin, boxMax;
    double scale;
};
#endif

/*
While gcc and clang have support for 128-bit integers, they do not support
128 bit constants! We have therefore to load up the constants from two 64 bit
values.
 */
uint128 load128(uint128 h, uint128 l) {
    uint128 q = (h<<64);
    q |= l;
    return q;
}

class Morton {
public:

    using dvec = SmallVec<double, 3>;
    using mortonKey = uint128;

    Morton() {}

    Morton(dvec &boxMin, dvec &boxMax) {
        initMorton( boxMin, boxMax );
    }

    template <class U>
    void initMorton(U _boxmin, U _boxmax) {
        boxMin = (SmallVec<double, 3>) _boxmin;
        boxMax = (SmallVec<double, 3>) _boxmax;
        for(int d=0; d<3; d++) scale[d] = 1.0/(boxMax[d]-boxMin[d]);

        m1  = 0x3ffffffffffUL;
        c64 = load128(    0x3ff00000000UL, 0x00000000ffffffffUL);
        c32 = load128(    0x3ff00000000UL, 0xffff00000000ffffUL);
        c16 = load128(0x30000ff0000ff00UL, 0x00ff0000ff0000ffUL);
        c8  = load128(0x300f00f00f00f00UL, 0xf00f00f00f00f00fUL);
        c4  = load128(0x30c30c30c30c30cUL, 0x30c30c30c30c30c3UL);
        c2  = load128(0x924924924924924UL, 0x9249249249249249UL);
        upper42mask0 = load128(0x0UL, 0x3ffffffffffUL);
        upper42mask0 = ~(upper42mask0<<124);
    }

    mortonKey makeKey(dvec v) {
        mortonKey key = 0;
        for(int d=0; d<3; d++) {
            double q = scale[d]*(v[d]-boxMin[d]) * (((uint128)1<<42));
            uint128 r = (uint128)q;
            r &= m1;
            r = (r | (r << 64)) & c64;
            r = (r | (r << 32)) & c32;
            r = (r | (r << 16)) & c16;
            r = (r | (r << 8))  & c8;
            r = (r | (r << 4))  & c4;
            r = (r | (r << 2))  & c2;

            //r = r & upper42mask0;
            key |= (r<<d);
        }

        return key;
    }

    int levelKey(mortonKey m, int l) {
        int shr = 123-3*l;
        return (m>>shr) & 7UL;
    }

    void printKey(mortonKey m) {
        out126Octal(m);
    }

    bool AM(mortonKey x, mortonKey y) {
        return ( x < y );
    }

    SmallVec<int, 3> mask[8] = { {-1,-1,-1},
                                 { 1,-1,-1},
                                 {-1, 1,-1},
                                 { 1, 1,-1},
                                 {-1,-1, 1},
                                 { 1,-1, 1},
                                 {-1, 1, 1},
                                 { 1, 1, 1} };

    uint128 m1, c64, c32, c16, c8, c4, c2, upper42mask0;
    dvec boxMin, boxMax, scale;
};


#endif // __MORTON_INCLUDE__

#ifdef TEST_MORTON

#include <cstdlib>

int main() {

    using dvec =  SmallVec<double, 3>;

    dvec zeros = {0.0,0.0,0.0};
    double eps =  pow(2.0, -42);
    dvec ones = {1.0+eps,1.0+eps,1.0+eps};

    MortonBinary<3,42> mortonA = MortonBinary<3,42>(zeros, ones);
    Morton morton = Morton(zeros, ones);

    for (int i=1; i<44; i++) {
        double q = pow(2.0, -i);
        q = 1.0 - q;
        dvec pos = {q, q, q};

        uint128 keyA = mortonA.makeKey(pos);
        uint128 key = morton.makeKey(pos);

        out126Octal(keyA); printf("\n");
        out126Octal(key); printf("\n\n");
        assert (key == keyA);
    }


    dvec pos = {0.125, 0.25, 0.5};
    uint128 key = morton.makeKey(pos);
    out126Octal(key); printf("\n\n");
    for(int l=0; l<42; l++) {
        printf("%d: %d\n",l, morton.levelKey(key, l));
    }

    for(int i=0; i<8; i++) {
        pprint("% d\n", morton.mask[i]);
    }

    srand(30495803);

    for(int i=0; i<10; i++) {
        pos = {drand48(), drand48(), drand48()};
        uint128 key = morton.makeKey(pos);
        out126Octal(key); printf("\n\n");
    }

    pos = {1.0, 1.0, 1.0};
    key = morton.makeKey(pos);
    out126Octal(key); printf("\n");
    out128Binary(key); printf("\n");
    for(int i=0; i<42; i++) {
        printf("%d %d\n", i, morton.levelKey(key, i));
    }

}

#endif // TEST_MORTON
