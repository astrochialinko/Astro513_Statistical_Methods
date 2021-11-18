/*
  Time-stamp: <util.cc on Friday, 12 November, 2021 at 11:18:50 MST (pinto)>

  Various useful functions

 */

#ifndef __UTIL_INCLUDE__
#define __UTIL_INCLUDE__

#include <cstdio>    // for fprintf
#include <algorithm> // for std::min, std::max
#include "smallvec.cc" // for SmallVec
#include "pprint.cc"

// Allocate "num" elements of the appropriate type to the pointer "name"
template <typename T>
void alloc(T *&name, int num) {
    using pointertype = typename std::remove_reference<decltype(name)>::type;
    using basetype = typename std::remove_pointer<pointertype>::type;
    name = new basetype[num];
}

// Error reporting;  the arguments to FATAL are identical to those for printf
// The function name, file, and line number are reported, then the arguments are passed to fprintf(stderr,...)
// and exit(1) is called.
#define FATAL(...) fatalError( __PRETTY_FUNCTION__, __FILE__, __LINE__, __VA_ARGS__)
template <typename... Args>
void fatalError(const char *func, const char *file, int line, const char *format, Args... args) {
    fpprint(std::cerr, format, args...);
    fpprint(std::cerr, "Error: %s in %s at line %d: ",func, file, line);
    fpprint(std::cerr, "\n");
    exit(1);
}

#define WARN(...) warnError( __PRETTY_FUNCTION__, __FILE__, __LINE__, __VA_ARGS__)
template <typename... Args>
void warnError(const char *func, const char *file, int line, const char *format, Args... args) {
    fpprint(std::cerr, "Warning: %s in %s at line %d: ",func, file, line);
    fpprint(std::cerr, format, args...);
    fpprint(std::cerr, "\n");
}

// functions to take max of an arbitrary number of arguments > 1
template <typename T>
T maxof( T a, T b) {
    return std::max(a, b);
}
template <typename T, typename... Args>
T maxof(T a, T b, Args... args) {
    return maxof(std::max(a, b), args...);
}

// function to take min of an arbitrary number of arguments > 1
template <typename T>
T minof( T a, T b) {
    return std::min(a, b);
}
template <typename T, typename... Args>
T minof(T a, T b, Args... args) {
    return minof(std::min(a, b), args...);
}

// Standard sign transfer function (note: if argument >=0 , return 1; <0 return -1)
template <class T>
T sign(T x) {
    return x<0 ? (T)(-1) : (T)1;
}

// square of POD argument
template <class T>
inline T sqr(T x) {return x*x;}

// dot product of two POD arrays
template <class T>
T dot(T *x, T *y, int len) {
    T sum = 0;
    for(int i=0; i<len; i++) sum += x[i] * y[i];
    return sum;
}

// CLOSE function for SmallVec's
template <typename T, typename U, int D>
inline int CLOSE(const SmallVec<T, D> a, const SmallVec<T, D> b, const U epsilon) {
    int yes = 1;
    for(int d=0; d<D; d++) {
        double foo = 0.5*( fabs((a[d] - b[d])) / (fabs(b[d]) + fabs(a[d])) );
        yes = yes && (foo < epsilon || fabs((a[d] - b[d])) < epsilon);
    }
    if(!yes) {
        fpprint(std::cerr, "CLOSE failure: % e % e : % e\n", a, b, a-b);
    }

    return yes;
}

// CLOSE function for POD
template <typename T, typename U>
inline int CLOSE(const T a, const T b, const U epsilon) {
    int yes = 1;
    T relerr = 0.5*( fabs((a - b)) / (fabs(b) + fabs(a)) );
    yes = yes && (relerr < epsilon || fabs((a - b)) < epsilon);
    if(!yes) {
        fpprint(std::cerr, "CLOSE failure: % e % e : % e\n", a, b, a-b);
    }
    return yes;
}

#endif // __UTIL_INCLUDE__
