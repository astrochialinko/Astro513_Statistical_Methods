/*
  Simple implementation of an N-body code entirely in C++.
 */



#include "STimer.cc"
#include "pprint.cc"
#include "util.cc"
#include "particleset.cc"
#include "bhtree.cc"

/*
  Hamiltonian class provides the update for the timestepper
 */
template <typename T>
class Hamiltonian {
public:

    using vec = SmallVec<T, 3>;

    Hamiltonian(double theta, int maxLeafSize, int maxSrcLen, int maxInteractionListLength, double epsSmooth)
        : theta(theta), maxLeafSize(maxLeafSize), maxSrcLen(maxSrcLen),
          maxInteractionListLength(maxInteractionListLength),
          epsSmooth(epsSmooth) {  }

    void getAcc(ParticleSet<T> &ps) {

        ps.computeBoundingBox();

         BHtree<T> BH = BHtree<T>(ps, maxLeafSize, epsSmooth, maxSrcLen);
        //BH->getStats = true;
        //BH->check = true;
        //BH->verbose = true;
        BH.makeTree();

        BH.BHsubsets(theta,  maxInteractionListLength, vec(0.0));
        BH.reportStats();
    }

    void update_v(ParticleSet<T> &ps) {
        // nothing to do
    }

    void update_a(ParticleSet<T> &ps) {
        // update acceleration from positions
        getAcc(ps);
    }

    //BHtree<T> *BH;

    int maxLeafSize, maxSrcLen, maxInteractionListLength;
    double theta, epsSmooth;

};

/*
  State holds the state of the system and implements the kick and drift operators
  using an instance of Hamiltonian.
 */
template <typename T>
class State {
public:
    State(ParticleSet<T> &ps, Hamiltonian<T> &H) : ps(ps), H(H), time(0) {}

    State *kick(double dt, bool recalc) {
        if(recalc) H.update_a(ps);
        for(int i=0; i<ps.N; i++) ps[i].v += dt * ps[i].a;
        return this;
    }

    State *drift(double dt, bool recalc) {
        if(recalc) H.update_v(ps);
        for(int i=0; i<ps.N; i++) ps[i].r += dt * ps[i].v;
        return this;
    }

    double time;
    int step=0;
    ParticleSet<T> &ps;
    Hamiltonian<T> &H;
};

/*
  step implements KDK or DKD leapfrog integrators.
*/
enum Method { KDK, DKD };
template <Method M, typename T>
void step(double dt, State<T> *s) {
    if constexpr (M == Method::KDK) {
            s = s->kick(0.5*dt,false)->drift(dt,false)->kick(0.5*dt,true);
            s->time += dt;
            s->step++;
        }
    if constexpr (M == Method::DKD) {
            s = s->drift(0.5*dt,false)->kick(dt,true)->drift(0.5*dt,false);
            s->time += dt;
            s->step++;
   }
}

/*
  getE computes the system totak kinetic and potential energy.
 */
template <typename T>
class Conserve {
public:
    Conserve(ParticleSet<T> &ps) : ps(ps) {
        getE(K0, P0, E0);
    }

    void getE(double &K, double &P, double &E) {
        K = P = 0;
        for(int i=0; i<ps.N; i++) {
            K += 0.5*ps[i].v.norm2()*ps[i].mass;
            P += ps[i].pot*ps[i].mass;
        }
        P *= 0.5;
        E = K+P;
    }

    void printCons() {
        double K, P, E;
        getE(K, P, E);
        pprint("E error: %.2e\n", (K+P-E0)/E0);
    }

    ParticleSet<T> &ps;
    double K0, P0, E0;
};

template <typename T>
void readDubinski(ParticleSet<T> &ps) {

    using vec = SmallVec<T,3>;

    std::string fname("dubinski.tab");

    // count number of lines in file
    std::ifstream f(fname,std::ios::in);
    int N = 0;
    std::string line;
    while(!f.eof()) { std::getline(f,line); N++; }
    f.close();
    N--;
    pprint("read %d lines from %s\n", N, fname);

    // read in data
    ps.reserve(N);
    f.open(fname,std::ios::in);
    for(int i=0; i<N; i++) {
        double m, x, y, z, vx, vy, vz;
        f >> m >> x >> y >> z >> vx >> vy >> vz;
        ps[i].mass = m;
        ps[i].r = vec(x,y,z);
        ps[i].v = vec(vx,vy,vz);
    }
    f.close();
}

int main() {

    using vec = SmallVec<double, 3>;

    int N = pow(2,20);
    ParticleSet<double> ps(0);
    readDubinski(ps);
    ps.computeBoundingBox();
    pprint("center: %e  halfWidth: %e\n\n", ps.centre, ps.halfWidth);

    double theta = 0.6;
    int maxLeafSize = 32;
    int maxSrcLen = ps.N;
    int maxInteractionListLength =5000;
    double epsSmooth = 0.98*pow(ps.N, -0.26) * ps.halfWidth;
    Hamiltonian<double> H(theta, maxLeafSize, maxSrcLen, maxInteractionListLength, epsSmooth);

    State s(ps, H);

    // zero-dt step to get potential
    step<KDK>(0.0, &s);
    Conserve C(ps);

    double dt = 0.1;

    STimer steptime;
    STimer fulltime;
    fulltime.Clear();
    while(s.time<20.0) {
        steptime.CStart();
        fulltime.Start();
        step<KDK>(dt, &s);
        steptime.Stop();
        fulltime.Stop();
        double rate = ps.N/steptime.Elapsed();
        if(s.step%10==0) { pprint("%5d rate: %e  time: %e  ", s.step, rate, fulltime.Elapsed()); C.printCons(); }
    }

}
