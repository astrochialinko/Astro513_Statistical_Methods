#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <omp.h>
#include <cassert>
#include "vformat.cc"

// include utility functions
#include "misc.cc"
#include "tallys.cc"
#include "random.cc"
#include "STimer.cc"

const double  clight = 2.99792458e8;    // m/s
const double hplanck = 6.62607004e-34;  // m^2 kg / s
const double  kboltz = 1.38064852e-23;  // m^2 kg / s^2 / K

const double      pi = M_PI;
const double   twopi = 2*pi;

#define PRINT(args...) printn(vformat(args).c_str())
void printit(const char *s) {
    fprintf(stderr,"%s", s);
}

// Container for a photon "packet"
class Packet {
public:
    Packet(double r,    // radius
           double mu,   // direction cosine
           double E,    // energy
           double W,    // weight
           int z) :     // current zone
           r(r), mu(mu), E(E), W(W), z(z) {}

    Packet(const Packet &p) : r(p.r), mu(p.mu), E(p.E), W(p.W), z(p.z) {} // copy constructor

    double r, mu, E, W;
    int z;
};

//--------------------------------------------------------------------------------
//
// Monte Carlo photon transport class
//
//  TODO:
//     scale tallies to correct fluxes
//     correct error bars
//
//--------------------------------------------------------------------------------

class RT {
public:

    typedef std::vector<double> Dvec;

    RT(Dvec &rr, Dvec &vv, Dvec &rhorho, int &nz) : r(rr), v(vv), rho(rhorho), nz(nz) {
        // set up tallying for energy deposition, etc.
        Dep = new TallyZone(nz);
        PdV = new TallyZone(nz);
        Inner = new TallyScalar();
        Emit = new TallyScalar();
    }

    ~RT() {
        delete Dep;
        delete PdV;
        delete Inner;
        delete Emit;
        delete Out;
        delete EmitSpec;
    }

    // Set up tallying for emission and emergent spectra
    void setSpectrumGrid(double eMin, double eMax, int ne) {
        Out = new TallyGrid(eMin, eMax, ne);
        EmitSpec = new TallyGrid(eMin, eMax, ne);
    }

    // Run the simulation
    void doSimulation(int Npackets, uint64_t seed) {

        N = Npackets;

        totalDepth(1);
        totalDepth(0.8);
        PRINT("\n");

        xi.seed(seed); // seed rng

        STimer all;  // execution timing
        all.CStart();

        // Parallelism is achieved by following packet histories independently,
        //   each with its own thread.
        // This for-loop is the only parallel construct.
        //   All global variables (r, v, rho) are read-only.
        //   The tally classes work on a per-thread basis.
        // The optimum chunk-size will depend upon the work per history
#pragma omp parallel for schedule(dynamic, 10)
        for(int i=0; i<N; i++) {                 // follow N photon histories
            if(i%10000==0) {
                PRINT("\r%7d K histories",i/1000);
                //printn(vformat("\r%7d K histories", i/1000).c_str());
            }
            Packet p = newPhoton(i);
            doPhoton(p);
        }
        PRINT("\n\n");

        all.Stop();
        PRINT("total time following histories: %e seconds\n", all.Elapsed());
        PRINT("rate: %e packets/second\n\n", N/all.Elapsed());

        // Sum over thread-sums and form statistics
        double scale = 1.0/N;
        Out->finish(scale);
        Dep->finish(scale);
        PdV->finish(scale);
        Inner->finish(scale);
        Emit->finish(scale);
        EmitSpec->finish(scale);

        report();
    }

    void totalDepth(double E) {
        Packet p = Packet(r[0], 1.0, E, 1.0, 0);
        double old_tauBoundaryCut = tauBoundaryCut;
        tauBoundaryCut = 1e30;
        double tau = tauToEscape(p);
        tauBoundaryCut = old_tauBoundaryCut;
        PRINT("total optical depth at E: %.4f  : %.4e\n", E, tau);
    }

    // Generate a new photon packet
    Packet newPhoton(int i) {
        double E = (1.2-0.8)*( (double)i/N) + 0.8;
        double mu = xi();
        double rr = r[0];
        double weight = 1;
        int zone = 0;

        // tally total emitted energy and emitted spectrum
        Emit->tally(E*weight);
        EmitSpec->tally(weight, E);

        return Packet(rr, mu, E, weight, zone);
    }

    //--------------------------------------------------------------------------------
    //
    // Main physics section: insert code to handle local interactions here
    //
    //--------------------------------------------------------------------------------

    // return a vector with extinction coefficient to various processes <<in lab frame>>
    // NB: last entry must be total of all previous
    // FIXME: devise a more-general way to keep track of this; a "process table"?
    Dvec extinctionCoefficient(Packet p) {

        double line = 1.0 * gaussian(p.E, 1.0, 0.025) * rho[p.z];
        double scatter = 1 * rho[p.z];
        double absorb = 0.0 * rho[p.z];
        double total = line + scatter + absorb;

        return std::vector<double>({line, scatter, absorb, total});
    }

    // perform a photon interaction with the material
    // kill is whether to terminate the photon
    void interactChoice(Packet &p, Dvec &chi, bool &kill) {
        // get total extinction coefficient
        double chi_total = chi.back();

        // Random choice of interaction
        double choose = xi();

        if(choose < chi[0]/chi_total) {   // line scatter
            double beta = getBeta(p);     // get velocity at photon's position
            double Ebefore = p.E;
            completeRedistribution(p, beta, 1.0, 0.025);
            //coherentScatter(p, beta);
            double dE = Ebefore - p.E;    // get PdV work
            double PdVdep = p.W * dE;     // PdV work to tally
            PdV->tally(PdVdep, p.z);       // Tally PdV work

            kill = false;                 // Don't kill photon
        }
        else if(choose < (chi[0]+chi[1])/chi_total) { // continuous scatter
            double beta = getBeta(p);     // get velocity at photon's position
            double Ebefore = p.E;
            coherentScatter(p, beta);
            double dE = Ebefore - p.E;    // get PdV work
            double PdVdep = p.W * dE;     // PdV work to tally
            PdV->tally(PdVdep, p.z);       // Tally PdV work

            kill = false;                 // Don't kill photon
        }
        else {                            // absorb
            double WEdep = p.W * p.E;     // Energy to deposit
            Dep->tally(WEdep, p.z);        // tally into zone p.z

            kill = true;                  // Packet will be terminated
        }
    }

    // perform a photon interaction with the material
    // kill is whether to terminate the photon
    void interactSurvival(Packet &p, Dvec &chi, bool &kill) {
        // get total extinction coefficient
        double chi_total = chi.back();

        // Survival biasing (aka implicit absorption)

        double rat = chi[1]/chi_total;    // Fraction to absorb
        double WEdep = p.E * p.W * rat;   // Absorbed energy
        Dep->tally(WEdep, p.z);           // tally into zone p.z

        p.W *= (1-rat);                   // Weight remaining after absorption

        double beta = getBeta(p);         // get velocity at photon's position
        double Ebefore = p.E;
        coherentScatter(p, beta);         // do the scatter
        double dE = Ebefore - p.E;        // change in lab frame photon energy
        double PdVdep = p.W * dE;         // add to absorbed energy in lab frame
                                          // don't reduce weight; energy comes from dE
        PdV->tally(PdVdep, p.z);          // Tally PdV work

        kill = false;                     // never kill packet for survival biasing
    }

    //--------------------------------------------------------------------------------
    //
    // Particle moving section
    //
    //--------------------------------------------------------------------------------

    // Find optical depth to problem boundary along straight line
    //   Make a copy of current photon & follow it in a straight line to the boundary,
    //   integrating as we go.
    //   If optical depth becomes larger than tauBoundaryCut, don't bother going further.
    double tauToEscape(Packet p) {

        double path;
        Packet pfake(p);                                               // make copy of packet

        double tau = 0;
        while (pfake.r < r.back() ) {               // while we don't arrive at boundary
            distToNextBoundary(pfake, path);                           // get remaining distance in this zone
            std::vector<double> chi = extinctionCoefficient(pfake);
            double chi_total = chi.back();
            tau += chi_total * path;                                   // integrate tau
            if( tau > tauBoundaryCut )                                 // if tau too large, quit
                return 1e38;
            move(pfake, path);                                         // otherwise move to next zone boundary
        }
        if (pfake.r <= r[0]) tau = 1e38;                               // run into inner boundary, don't force

        return tau;
    }

    // Move packet to next event
    //   zone boundary crossing or interaction
    void travel(Packet &p, bool &kill) {
        double path, pathToBoundary;

        // do forced scattering ("peel off")
        //   find fraction of weight which makes it to surface without further interaction,
        //   and tally that weight
        double tauEscape = tauToEscape(p);    // find optical depth to surface
        double escape = 0;
        if( tauEscape < tauBoundaryCut ) {    // if small enough...
            escape = std::exp(-tauEscape);    // determine escape fraction
            double wtescape = p.W * escape;   // and tally outgoing energy
            Out->tally(wtescape, p.E);
            p.W *= (1.0-escape);              // reduce weight accordingly
        }

        double dtau = -std::log( 1.0 - xi()*(1-escape) );    // draw from P(tau | doesn't get to surface)
        Dvec chi = extinctionCoefficient(p);
        double chi_total = chi.back();
        path = dtau/chi_total;                               // path length to next interaction, assuming chi
                                                             // from this zone; packet won't get any further...
        distToNextBoundary(p, pathToBoundary);               // path length to next zone boundary

        if(path < pathToBoundary) {                          // move whichever is less
            move(p, path);
#ifdef SURVIVAL
            interactSurvival(p, chi, kill);                  // do interaction via survival
#else
            interactChoice(p, chi, kill);                    // interact using random choice
#endif
        }
        else {
            move(p, pathToBoundary);                         // nothing to do at zone crossing
            kill = false;
        }

        return;
    }

    // Follow packet until it is terminated
    void doPhoton(Packet &p) {

        bool kill;
        double WE;

        while(true) {

            travel(p, kill);                   // Move to next event and do it

            if(kill) {
                break;                    // If the interaction terminated the packet, quit loop.
            }

            if(p.W < lowWeight) {              // If packet weight falls too low, use Russian Roulette.
                if (xi() < p.W/avgWeight)      // Terminate packet if xi > p.W/avgWeight, else give avgWeight.
                    p.W = avgWeight;
                else {
                    break;
                }
            }

            if (p.r <= r[0]) {                 // If packet intersects inner boudary, absorb it there.
                Inner->tally(p.W*p.E);
                break;
            }

            if( p.r >= r.back() ) {            // A packet should never arrive at the surface!
                //assert(1==0);
                PRINT("bad escape!\n");
                break;
            }
        }
    }

    // Find distance to next zone boundary along packet's direction cosine
    // Probably a cleaner way exists:
    //    special case where photon is epsilon away from inner zone radius.
    //    problems with doing if-tests on floating point!
    void distToNextBoundary(Packet &p, double &path) {
        double disc, rb, sign;
        double mu = p.mu;

        int z = p.z;

        if(p.r == 0)
            disc = -2;                       // special case for r=0 inner boundary;
                                             // i.e. no inner boundary
        else {
            double rat = r[z]/p.r;
            disc = -sqrt( 1 - rat*rat );
        }

        if (mu > disc) {                     // next zone is r[p.z+1] & we are outgoing
            rb = r[p.z+1];
            sign = 1;
        }
        else {
            if(p.r > r[p.z]) {
                rb = r[p.z];                 // ingoing
                sign = -1;
            }
            else {                           // special case of p.r==r[p.z]
                double rat = r[p.z-1]/p.r;
                disc = -sqrt( 1 - rat*rat );
                if(mu>disc) {                // will intersect r[p.z] again through zone p.z-1
                    rb = r[p.z];             // technically, outgoing, but back to same zone boundary
                    sign = 1;
                }
                else {
                    rb = r[p.z-1];           // will intersect r[p.z-1] through zone p.z-1
                    sign = -1;               // ingoing
                }
            }
        }

        path = -mu*p.r + sign*sqrt(mu*mu*p.r*p.r - p.r*p.r + rb*rb);
        path = std::max(path, 1.0e-8*p.r);  // force a move of epsilon
    }

    // update packet r, mu, and zone by path along current mu
    void move(Packet &p, double path) {
        double r0 = p.r;
        double mu = p.mu;

        p.r = sqrt(p.r*p.r + path*path + 2*p.r*path*p.mu);
        if(p.r>0)
            p.mu = (path+r0*p.mu)/p.r;
        else
            p.mu = 1.0;

        p.z = idx(r, p.r);  // get new zone number
    }

    //--------------------------------------------------------------------------------
    //
    // Interaction processes:
    //
    //--------------------------------------------------------------------------------

    // perform isotropic coherent scattering in co-moving frame moving at beta w/r to lab frame
    void coherentScatter(Packet &p, double beta) {
        invlorentz(p, beta);
        p.mu = 2*xi()-1;
        lorentz(p, beta);
    }

    // Next three functions sample zero-temperature Compton scattering in co-moving
    //    frame moving at beta w/r to lab frame
    // NB: these expect photon energy in units of electron rest mass (0.511 MeV)
    double r0sq = 10*(2.82*2.82);
    double r0sqpi = r0sq/2.0 * twopi;

    // cumulative distribution for scattering through a = cos(theta)
    double aklein(double ei, double a) {
        double bax = 1.+ei+ei*a;
        double fa=2.*(a+1.)/(ei*ei)+1./ei*std::log(bax)*(1.-2./(ei*ei)*(1.+ei));
        fa=fa-(2.*(a*a-1.)/bax+1./(bax*bax)-1.)/2./ei;
        return fa*r0sqpi;
    }

    // Solve
    //    xi = aklein(energy in, cosine of scattering angle)
    // using bisection (pretty crude!)
    void sampleCompton(double &mu, double &E) {
        double anew = 0;
        double lo = -1, hi = 1;

        double sigma = aklein(E, 1.0);
        double area = xi();

        double diff = 1;
        while(diff > 1e-5) {
            double a = anew;
            double fa = aklein(E,a)/sigma;
            if(fa>area)
                hi = a;
            else
                lo = a;
            anew = 0.5*(hi+lo);
            diff = fabs(anew-a);
        }

        double sin2 = 1 - mu*mu;
        double cosphi = 2*xi()-1;
        mu = -mu*anew - sqrt(sin2*(1.0-anew*anew))*cosphi;
        E = E/(1+E+E*anew);
    }

    // do the scattering and return energy given to the scattering medium
    void comptonScatter(Packet &p, double beta, double &dE) {
        double E0 = p.E;
        invlorentz(p, beta);
        sampleCompton(p.mu, p.E);
        lorentz(p, beta);
        dE = E0 - p.E;
    }

    // perform complete redistribution scattering from a Gaussian (Doppler) line profile
    //   with line center energy e0 and width w
    // NB: extinction coefficient needs to use same e0 & w
    void completeRedistribution(Packet &p, double beta, double e0, double w) {
        invlorentz(p, beta);
        p.mu = 2*xi()-1;
        p.E = sampleNormal(e0, w);
        lorentz(p, beta);
    }

    // Sample from a normal distribution by Box-Muller
    double sampleNormal(double mean, double sigma) {
        double x = std::sqrt( -2*std::log( xi() ) );
        double y = x * std::cos(twopi * xi());
        return y*sigma + mean;
    }

    // interpolate in the velocity profile to get beta at site of interaction
    // currently used in the Lorentz transform functions below
    double getBeta(Packet p) {
        // interpolate veocity
        double eps = std::max(0.0, std::min(1.0, ( p.r - r[p.z] ) / ( r[p.z+1] - r[p.z] )));
        double beta = ( eps*v[p.z+1] + (1-eps)*v[p.z] ) / clight;
        assert( (eps>=0) && (eps<1) );
        return beta;
    }

    // Lorentz transformation from lab frame at rest into co-moving frame moving at beta
    void lorentz(Packet &p, double beta) {
        p.mu = (p.mu+beta)/(1+beta*p.mu);
        p.E = p.E*sqrt(1-beta*beta)/(1-beta*p.mu);
    }

    // Lorentz transformation from co-moving frame moving at beta into lab frame at rest
    void invlorentz(Packet &p, double beta) {
        p.mu = (p.mu - beta)/(1-beta*p.mu);
        p.E = p.E*(1-beta*p.mu)/sqrt(1-beta*beta);
    }

    //--------------------------------------------------------------------------------
    // reporting function
    //--------------------------------------------------------------------------------
    void report() {
        double Escape = Out->total;
        double Deposit = Dep->total;
        double InnerDep = Inner->total;
        double PdVdep = PdV->total;
        double EmitTot = Emit->total;
        PRINT("Deposit: %e %e  Escape: %e  PdV: %e   Emit: %e\n",
                       Deposit, InnerDep,  Escape, PdVdep, EmitTot);
        PRINT("1-sum: %e\n", 1 - (Deposit+InnerDep+Escape+PdVdep)/EmitTot);
    }

    int nz, N;
    Dvec r, v, rho;
    MyRNG xi;

    TallyGrid *Out = 0;
    TallyZone *Dep = 0;
    TallyZone *PdV = 0;
    TallyScalar *Inner = 0;
    TallyScalar *Emit = 0;
    TallyGrid *EmitSpec = 0;

    // parameters:
    // Should make getter/setter for these

    // stop integration to boundary (for forced scattering)
    //   if tau > tauBoundaryCut
    double tauBoundaryCut = 14; // e^{-14} < 1e-6

    // Russian Roulette for low packet weight:
    // if p.W < lowWeight, packet killed if xi < lowWeight/avgWeight, else p.W = avgWeight
    double lowWeight = 0.001;
    double avgWeight = 0.01;

    void (*printn)(const char *s) = printit;

};



//--------------------------------------------------------------------------------
// Getters & Setters for Python interface.
//--------------------------------------------------------------------------------
// function to take a c-style array double *a and fill a std::vector<double> v
std::vector<double> toVector(double *a, int na) {
    std::vector<double> v;
    v.assign(a, a+na);
    return v;
}

// function to take a std::vector<double> v and fill a c-style array double *a
void fromVector(double *a, std::vector<double> vec) {
    std::copy( vec.begin(), vec.end(), a);
}

// These must use "C" (not C++) naming to avoid C++'s name-mangling
extern "C" {
    // Call the RT constructor and return a pointer to the instance
    RT* newRT(double *r, double *v, double *rho, int nz) {

        // FIXME: vector's don't stick around after function exits!
        std::vector<double> rvec = toVector(r, nz);
        std::vector<double> vvec = toVector(v, nz);
        std::vector<double> rhovec = toVector(rho, nz);

        RT* me = new RT(rvec, vvec, rhovec, nz);
        return me;
    }

    void setSpectrumGrid(RT* me, double eMin, double eMax, int ne) {
        me->setSpectrumGrid(eMin, eMax, ne);
    }

    void doSimulation(RT* me, int N, uint64_t seed) {
        me->doSimulation(N, seed);
    }

    void getSpectrum(RT* me, double *egrid, double *mean, double *sigma) {
        fromVector(egrid, me->Out->xgrid);
        fromVector(mean, me->Out->mean);
        fromVector(sigma, me->Out->sig);
    }

    void getDeposition(RT* me, double *mean, double *sigma) {
        fromVector(mean, me->Dep->mean);
        fromVector(sigma, me->Dep->sig);
    }

    void getPdV(RT* me, double *mean, double *sigma) {
        fromVector(mean, me->PdV->mean);
        fromVector(sigma, me->PdV->sig);
    }

    void getInner(RT* me, double *mean, double *sigma) {
        *mean = me->Inner->mean;
        *sigma = me->Inner->sig;
    }

    void getEmitSpec(RT* me, double *egrid, double *mean, double *sigma) {
        fromVector(egrid, me->EmitSpec->xgrid);
        fromVector(mean, me->EmitSpec->mean);
        fromVector(sigma, me->EmitSpec->sig);
    }

    void registerPrint(RT *me, void (*print_callback)(const char *s)) {
        me->printn = print_callback;
    }

}

// compile w/ main unless PYTHONLIB is set
#ifndef PYTHONLIB
int main() {

    // make model
    int nz = 11;
    double rmin = 0, rmax = 10;
    double vmin = 0, vmax = 0.0;//1*clight;
    double rhomin = 1, rhomax = 1;

    std::vector<double> r = linspace(rmin, rmax, nz);
    std::vector<double> v = linspace(vmin, vmax, nz);
    std::vector<double> rho = linspace(rhomin, rhomax, nz);

    // get instance of RT
    RT sim(r, v, rho, nz);

    // define emergent spectrum grid
    double eMin = 0.8, eMax = 1.2;
    int ne = 101;
    sim.setSpectrumGrid(eMin, eMax, ne);

    // run simulation
    int N = 1<<10;
    uint64_t seed = 2873642343;
    sim.doSimulation(N, seed);

    sim.report();

    FILE *sout = fopen("spectrum", "w");
    for(int i=0; i<ne; i++) {
        fprintf(sout, "%4d   %6.4f   %.4e  %.4e  %0.4e\n", i, sim.Out->xgrid[i], sim.Out->mean[i],
                sim.Out->sig[i],
                sim.EmitSpec->mean[i]);
    }
    fclose(sout);

    FILE *sdep = fopen("deposition", "w");
    for(int i=0; i<nz; i++) {
        fprintf(sdep, "%4d  %6.4f  %.4e  %.4e\n", i, r[i], sim.Dep->mean[i], sim.Dep->sig[i]);
    }
    fclose(sdep);

}
#endif
