#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <cstdio>
#include <cmath>
#include <vector>

// Return evenly-spaced vector on [a,b] with n elements. NB: type is that of lo and hi
template <class T>
std::vector<T> linspace(T lo, T hi, int n) {
    std::vector<T> result(n);

    T dv = (hi-lo)/(n-1);
    for(int i=0; i<n; i++) result[i] = lo + i*dv;
    return result;
}

// bound i on [lo,hi]
inline int bound(int i, int lo, int hi) {
    return std::max(std::min(i,hi), lo);
}

/*
  Given a value val and a vector vec, returns index i such that
  vec[i] <= val < vec[i+1]
  returns -1 if val<vec[0]
  return length of vec if val>=vec[last]

  WARNING: does floating-point comparisons!
 */
int idx(std::vector<double> vec, double val) {
    if(val < vec[0]) return -1;
    if(val > vec.back()) return vec.size();
    double dx = vec[1]-vec[0];
    int i = std::floor((val-vec[0])/dx);
    return i;
}

double gaussian(double x, double mu, double s) {
    return 1/(std::sqrt(M_PI)*s) * std::exp(-0.5*((x-mu)*(x-mu)/(s*s)));
}


#endif // INCLUDE_MISC

#if 0

int main() {

    int nx = 101;
    std::vector<double> xgrid = linspace(0.5, 2.5, nx);
    for(int i=0; i<nx; i++) printf("%4d %e\n", i, xgrid[i]);

    int i = idx(xgrid, 1.0);
    printf("%d\n", i);
}




#endif
