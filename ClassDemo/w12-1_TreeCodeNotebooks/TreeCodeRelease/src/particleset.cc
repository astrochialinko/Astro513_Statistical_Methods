#ifndef __PARTICLE_SET_CC__
#define __PARTICLE_SET_CC__

/*

  Implements a set of particles defined by the Ptype<T> using AOS
  for the particle data.

 */


#include <cstring>
#include <numeric>

#include "util.cc"
#include "pprint.cc"
#include "AOS.cc"

template <class T>
class ParticleSet : public AOS<Ptype<T>> {
public:

    using ptype = Ptype<T>;

    using AOS<ptype>::operator[];

    ParticleSet() : AOS<ptype>(0), bboxComputed(false) {}
    ParticleSet(int N) : AOS<ptype>(N), bboxComputed(false) {}

    void computeBoundingBox() {
        assert( this->N > 0 );

        double rmin[3], rmax[3];

        for(int d=0; d<3; d++) {
            rmin[d] = std::numeric_limits<double>::max();
            rmax[d] = -rmin[d];
        }

        for(int i=0; i<this->N; i++)
            for(int d=0; d<3; d++) {
                rmin[d] = std::min(rmin[d], AOS<ptype>::operator[](i).r[d]);
                rmax[d] = std::max(rmax[d], AOS<ptype>::operator[](i).r[d]);
            }

        for(int d=0; d<3; d++) centre[d] = 0.5 * (rmin[d] + rmax[d]);

        halfWidth = 0;
        for(int d=0; d<3; d++)
            halfWidth = maxof( halfWidth, (double)(fabs(rmax[d]-centre[d])), (double)(fabs(rmin[d]-centre[d])) );

        bboxComputed = true;
    }

    inline bool hasBoundingBox() {
        return bboxComputed;
    }

    template <typename U>
    void permute(U *origPerm, U len) {
        U i,j,k;
        Ptype<T> q;

        U *perm = new U[len];
        for(i=0; i<len; i++) perm[i] = origPerm[i];

        for(i=0; i<len; i++ ) {
            q = AOS<ptype>::operator[](i);
            for(k=i; perm[k]!=i; k=perm[j], perm[j]=j) {
                j = k;
                AOS<ptype>::operator[](k) = AOS<ptype>::operator[](perm[k]);
            }
            AOS<ptype>::operator[](k) = q;
            perm[k] = k;
        }
        delete[] perm;
    }

    // Provide a string describing the Ptype<T> structure.
    std::string getTypeString() {
        std::string line;
        int len = Ptype<T>::defs.size();
        for(int i=0; i<len-2; i+=2) {
            line += std::string(Ptype<T>::defs[i]) + ":" + std::string(Ptype<T>::defs[i+1]) + ";";
        }
        line += std::string(Ptype<T>::defs[len-2]) + ":" + std::string(Ptype<T>::defs[len-1]);
        return line;
    }


    void write(const std::string &basename) {
        std::string fname = basename + ".pset";
        std::ofstream stream(fname.c_str(), std::ios::out | std::ios::binary);
        if( !stream )
            FATAL("particleSet: write: cannot open file '%s' for writing: %s\n", fname, strerror(errno));

       // write descriptor of data types in file
        std::string line = getTypeString();
        int len = line.size();
        const char *str = line.c_str();
        int bytes = 0;
        stream.write((char*)&len,        sizeof(len));
        bytes += sizeof(len);
        stream.write((char*)str,                 len+1);
        bytes += len+1;
        stream.write((char*)&(this->N),  sizeof(this->N));
        bytes += sizeof(this->N);
        stream.write((char*)&centre,     sizeof(centre));
        bytes += sizeof(centre);
        stream.write((char*)&halfWidth,  sizeof(halfWidth));
        bytes += sizeof(halfWidth);
        pprint("header length: %d\n", bytes);
        AOS<Ptype<T>>::write(stream);
        bytes+= this->N*sizeof(Ptype<T>);
        stream.close();

        pprint("ParticleSet: wrote %d particles, %d bytes to '%s'\n", this->N, bytes, fname);
    }

    void read(const std::string &basename, bool resize) {
        std::string fname = basename + ".pset";

        std::ifstream stream(fname.c_str(), std::ios::in | std::ios::binary);
        if( !stream )
            FATAL("particleSet: read: cannot open file '%s' for reading: %s\n", fname, strerror(errno));

        int len;
        stream.read((char*)&len,        sizeof(len));
        char str[len];
        stream.read(str,                len+1);
        // be sure we are reading the same kind of data which was written
        std::string line = getTypeString();
        if( line != std::string(str) )
            FATAL("Pset: read: current Ptype does not have same format as in file: %s\n", fname);

        decltype(this->N) np;
        stream.read((char*)&np,        sizeof(np));
        stream.read((char*)&centre,    sizeof(SmallVec<double, 3>));
        stream.read((char*)&halfWidth, sizeof(double));
        if( resize )
            this->reserve(np);
        else if( np != this->N ) {
            FATAL("ParticleSet: read: won't resize due to flag: is: %d  was: %d\n", np, this->N);
        }
        AOS<Ptype<T>>::read(stream);
        stream.close();

        pprint("ParticleSet: read %d particles from '%s'\n", this->N, fname);
    }

    bool bboxComputed;
    SmallVec<T, 3> centre;
    T halfWidth;
};

#endif // __PARTICLE_SET_CC__


#ifdef TEST_PARTICLESET
/*
  Compile with:
     g++ -I../lib -std=c++2a -Ofast -DTEST_PARTICLESET particleset.cc 2>&1 | less
 */


template <typename T>
void printP(const PtypeRef<T> &p) {
    pprint("%4d  % .4e % .4e % .4e % .4e %.2f\n", p.id,  p.r, p.v, p.a, p.pot, p.mass);
}

template <class T>
Ptype<T> func(int i, T x) {

    using vec = SmallVec<T, 3>;
    double fi = (double)i;

    vec r = vec(fi, 2*fi, 3*fi);
    vec v = vec(sqrt(fi), sqrt(2*fi), sqrt(3*fi));
    vec a = vec(fi*fi, fi*fi, 9*fi*fi);
    T pot = v.norm();
    T mass = r.norm();
    uint id = i;

    Ptype<T> Q = {r, v, a, pot, mass, i};
    return Q;
}

int main() {

    using Real = double;
    using ptype = Ptype<Real>;
    using vec = SmallVec<Real, 3>;

    int N = 1000;
    ParticleSet<Real> s(N);

    pprint("N: %d\n", s.N);

    // fill with data
    int sign = -1;
    for(int i=0; i<N; i++) {
        ptype Q = func(i, (Real)(1.0));
        s.set(i, Q);
    }

    s.computeBoundingBox();
    pprint("Bounding Box: % e % e\n", s.centre, s.halfWidth);

    // save set
    {
        s.write("outfoo");
    }

    // fill set with zeros
    for(int i=0; i<N; i++) {
        s[i] = {vec(0,0,0), vec(0,0,0), vec(0,0,0), 0.0, 0.0, 0};
    }

    // read
    {
        s.read("outfoo", true);
    }

    s.computeBoundingBox();
    pprint("Bounding Box: % e % e\n", s.centre, s.halfWidth);

    // permute
    int *findex = new int[s.N];
    for(int i=0; i<s.N; i++) findex[i] = s.N-i-1;
    int *rindex = new int[s.N];
    for(int i=0; i<s.N; i++) rindex[findex[i]] = i;

    s.permute(findex, (int)s.N);

    for(int i=0; i<s.N; i++) {
        ptype Q = func(findex[i], (Real)(1.0) );
        assert (Q.r == s[i].r);
    }
    pprint("OK\n");
    s.computeBoundingBox();
    pprint("Bounding Box: % e % e\n", s.centre, s.halfWidth);

    pprint("%s\n", s.getTypeString() );

}

#endif // TEST_PARTICLE_SET
