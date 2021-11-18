/*
Time-stamp: <smallvec.cc on Tuesday, 30 October, 2018 at 00:38:55 MST (pinto)>

SmallVec class: a class for small vectors of static length

NOTE:
   This code relies on various template magic; you will likely need to use the
    -std=c++0x or -std=c++11 or -std=c++17 etc. flag when compiling

 */

#ifndef __SMALLVECTOR_CC__
#define __SMALLVECTOR_CC__
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <initializer_list>
#include <typeinfo>

// Smallvec class

template <class T, int D>
class SmallVec {
public:
    T data[D];

    // constructor to set data to zero:
    //   SmallVec<double, 3> x; -> x = (0,0,0)
    SmallVec() : data{0} {}

    // constructor to broadcast a scalar to all components:
    //    SmallVec<double, 3> x(3); -> x = (3,3,3)
    // (how to do this at compile time?)
    template <class U>
    explicit
    SmallVec( const U& scalar) {
        for(int i=0; i<D; i++) data[i] = T(scalar);
    }

    // initialzation constructor (e.g., SmallVec<double, 3>(3, 4.5, 5): )
    template <typename... Tail>
    SmallVec(typename std::enable_if<sizeof...(Tail)+1==D, T>::type head, Tail... tail)
        : data{ head, T(tail)... } {}

    // converting copy constructor:
    //    SmallVec<U, D> q(1,2,3);
    //    SmallVec<T, D> = q;
    template <typename U>
    SmallVec(const SmallVec<U,D>& other) {
        for(int i=0; i<D; i++) data[i] = T(other.data[i]);
    }

    // default copy constructor:
    //    SmallVec<double, 3> q(1,2,3)
    //    SmallVec<double, 3> z = q;
    SmallVec(const SmallVec<T,D>& other) = default;

    // allow access to components with usual vector notation
    T operator[] ( const size_t i ) const {
        assert( i < D );
        return *(data+i);
    }
    T& operator [] ( const size_t i ) {
        assert( i < D );
        return *(data+i);
    }

    // Copy assignment, possibly with type promotion
    template <class U>
    SmallVec& operator=(const SmallVec<U,D>& rhs) {
        for(int i=0; i<D; i++) data[i] = rhs.data[i];
        return *this;
    }

    // zero the components
    SmallVec<T, D> zero() {
        for(int i=0; i<D; i++) data[i] = 0;
        return *this;
    }

    // --------------------------------------------------------------------------------

    //unary minus
    SmallVec operator -() {
        SmallVec<T,D> result;
        for(int i=0; i<D; i++) result[i] = -data[i];
        return result;
    }

    // --------------------------------------------------------------------------------

    // combinations of arithmetic functions and =
    template <class U>
    SmallVec<T,D>& operator += ( const SmallVec<U,D>& rhs ) {
        for(int i=0; i<D; i++) data[i] += rhs.data[i];
        return *this;
    }

    template <class U>
    SmallVec<T,D>& operator -= ( const SmallVec<U,D>& rhs ) {
        for(int i=0; i<D; i++) data[i] -= rhs.data[i];
        return *this;
    }

    template <class U>
    SmallVec<T,D>& operator *= ( const U rhs ) {
        for(int i=0; i<D; i++) data[i] *= rhs;
        return *this;
    }

    template <class U>
    SmallVec<T, D>& operator /= ( const U rhs ) {
        for(int i=0; i<D; i++) data[i] /= rhs;
        return *this;
    }

    // --------------------------------------------------------------------------------

    // functions of components:

    // 2-norm:
    typename std::common_type<T, double>::type norm() const {
        typename std::common_type<T, double>::type sum = 0;
        for(int i=0; i<D; i++) sum += data[i]*data[i];
        return sqrt(sum);
    }

    // square of 2-norm:
    typename std::common_type<T, double>::type norm2() const {
        typename std::common_type<T, double>::type sum = 0;
        for(int i=0; i<D; i++) sum += data[i]*data[i];
        return sum;
    }

    T maxcomponent() const {
        T maxc = data[0];
        for(int i=1; i<D; i++) maxc = (data[i]>maxc)?data[i]:maxc;
        return maxc;
    }
    T mincomponent() const {
        T minc = data[0];
        for(int i=1; i<D; i++) minc = (data[i]<minc)?data[i]:minc;
        return minc;
    }

    // return true if all components are on [a,b]
    template <class U, class V>
    bool inrange(U low, V hi) {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] >= low && data[i] <= hi);
        return val;
    }

    // return true if all components are on [a,b], component-wise
    template <class U, class V>
    bool inrange(SmallVec<U, D> low, SmallVec<V, D> hi) {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] >= low[i] && data[i] <= hi[i]);
        return val;
    }

    // --------------------------------------------------------------------------------

    // dot and cross products:

    // dot product
    template<class U>
    typename std::common_type<T, U>::type dot(const SmallVec<U, D>& rhs) const {
        typename std::common_type<T, U>::type tmp = data[0]*rhs[0];
        for(int i=1; i<D; i++) tmp += data[i]*rhs[i];
        return tmp;
    }

    // cross product only defined if D==3:
    template <class U>
    inline SmallVec<typename std::common_type<T, U>::type, 3>
    cross( const SmallVec<U, 3>& rhs ) const {
        SmallVec<typename std::common_type<T, U>::type, 3> tmp;
        tmp[0] = data[1] * rhs[2] - data[2] * rhs[1];
        tmp[1] = data[2] * rhs[0] - data[0] * rhs[2];
        tmp[2] = data[0] * rhs[1] - data[1] * rhs[0];
        return tmp;
    }

    // --------------------------------------------------------------------------------

    // absolute difference between *this and rhs is less than epsilon in each component
    template <class U, class V>
    bool absclose( const SmallVec<U, D>& rhs, const V epsilon ) const {
        SmallVec<typename std::common_type<T, U>::type, D> diff;
        diff = *this - rhs;
        bool val = true;
        for(int i=0; i<D; i++) val = val && (fabs(diff[i]) < epsilon);
        return val;
    }

    // relative difference between *this and rhs is less than epsilon in each component
    template <class U, class V>
    int relclose( const SmallVec<U, D>& rhs, const V epsilon ) const {
        SmallVec<typename std::common_type<T, U>::type, D> sum, diff;
        for(int i=0; i<D; i++) sum[i] = fabs(data[i]) + fabs(rhs[i]);
        diff = *this - rhs;
        bool val = true;
        for(int i=0; i<D; i++) val = val && (( 2*fabs(diff[i]) / sum[i] ) < epsilon);
        return val;
    }

    // --------------------------------------------------------------------------------

    bool is_finite() const {
        bool val = true;
        for(int i=0; i<D; i++) val = val && std::isfinite(data[i]);
        return val;
    }

    // --------------------------------------------------------------------------------

    // relational operators:

    template <class U>
    bool operator == ( const SmallVec<U, D>& rhs ) const {
        bool val = true;
        for(int i=0; i<D; i++) val = val && (data[i] == rhs[i]);
        return val;
    }

    template <class U>
    bool operator != ( const SmallVec<U, D>& rhs ) const {
        return !(*this == rhs);
    }

    // define < as in sorted order by dimension
    template <class U>
    bool operator < ( const SmallVec<U, D>& rhs ) const {
        if(data[0] < rhs[0]) return true;
        for(int i=1; i<D; i++)
            if( data[i-1]==rhs[i-1] && data[i] < rhs[i] ) return true;
        return false;
    }

    template <class U>
    bool operator <= ( const SmallVec<U, D>& rhs ) const {
        if(  (*this < rhs) || (*this==rhs) ) return true;
        return false;
    }

    // define > as in reverse sorted order by dimension
    template <class U>
    bool operator > ( const SmallVec<U, D>& rhs ) const {
        if(data[0] > rhs[0]) return true;
        for(int i=1; i<D; i++)
            if( data[i-1]==rhs[i-1] && data[i] > rhs[i] ) return true;
        return false;
    }

    template <class U>
    bool operator >= ( const SmallVec<U, D>& rhs ) const {
        if( (*this>rhs) || (*this==rhs) ) return true;
        return false;
    }

    // --------------------------------------------------------------------------------

    // stream operator: keep the format settings from being destroyed by the
    // non-numeric characters output
    friend std::ostream& operator <<( std::ostream& o, const SmallVec<T, D>& v ) {
        std::streamsize tmpw = o.width();
        std::streamsize tmpp = o.precision();
        char tmps = o.fill();
        std::ios::fmtflags tmpf = o.flags();  // format flags like "scientific" and "left" and "showpoint"
        o << std::setw(1);
        o << "(";
        for(int i=0; i<D-1; i++) {
            o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
            o << v.data[i];
            o << ",";
        }
        o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
        o << v.data[D-1];
        o << ")";
        return o;
    }

};

// componentwise addition and subtraction:
template <class T, class U, int D>
inline SmallVec<typename std::common_type<T, U>::type, D>
operator + (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename std::common_type<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i] + rhs.data[i];
    return tmp;
}

template <class T, class U, int D>
inline SmallVec<typename std::common_type<T, U>::type, D>
operator - (const SmallVec<T, D>& lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename std::common_type<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i] - rhs.data[i];
    return tmp;
}

// left and right multiplication by a scalar:
template <class T, class U, int D>
inline SmallVec<typename std::common_type<T, U>::type, D>
operator * ( const SmallVec<T, D>& lhs, const U rhs ) {
    SmallVec<typename std::common_type<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i]*rhs;
    return tmp;
}

template <class T, class U, int D>
inline SmallVec<typename std::common_type<T, U>::type, D>
operator * ( const T lhs, const SmallVec<U, D>& rhs ) {
    SmallVec<typename std::common_type<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs*rhs.data[i];
    return tmp;
}

// right division by a scalar:
template <class T, class U, int D>
inline SmallVec<typename std::common_type<T, U>::type, D>
operator / ( const SmallVec<T, D>& lhs, const U rhs ) {
    SmallVec<typename std::common_type<T, U>::type, D> tmp;
    for(int i=0; i<D; i++) tmp.data[i] = lhs.data[i]/rhs;
    return tmp;
}

#endif // __SMALLVECTOR_CC__


#ifdef TESTSMALLVEC

// Not really a test -- just some examples of how to use smallvec

#include <iostream>
#include <complex>
using namespace std;
typedef complex<double> Complex;

// abbreviations for the typographically challenged
typedef SmallVec<int, 3> ivec;
typedef SmallVec<double, 3> fvec;
typedef SmallVec<double, 3> dvec;
typedef SmallVec<Complex, 3> cvec;

struct body {
    dvec pos, vel, acc;
    double mass;
};


int main() {

    fvec a;
    dvec b, c;
    a = b + c;

    // some simple vector algebra (only works in 3D)
    dvec xhat(1,0,0);
    dvec yhat = dvec(0,1,0);
    dvec zhat = xhat.cross(yhat);
    cout << "zhat = " << zhat << endl;

    double projx = dvec(1.2, 3.1, -2.7).dot(xhat);
    cout << "projx = " << projx << endl;

    // complex arithmetic
    cvec cfoo;
    cfoo = cvec(Complex(1,0), Complex(0,1), Complex(1,1));
    cout << "Complex vector = " <<cfoo << endl;
    cout << "norm of that vector = " << cfoo.norm() << endl;

    // make some particles in the box (0,0,0) to (1,1,1)
    int np = 1000;
    body *p = (body *)malloc(sizeof(body)*np); assert(p != (body *)NULL);
    for(int i=0; i<np; i++) {
        p[i].pos = dvec(drand48(), drand48(), drand48());
        p[i].mass = 1;
        p[i].acc.zero();
    }

    // now compute the gravitational acceleration with Plummer softening
    double eps = 0.0001;
    for(int i=0; i<np; i++) {
        for(int j=0; j<np; j++) {
            dvec dr = p[i].pos - p[j].pos;
            double denom = 1.0/sqrt(dr.norm2() + eps*eps);
            denom = denom*denom*denom;
            p[i].acc += dr*denom*p[j].mass;
        }
    }

    double dt = 0.001;
    // kick the velocities
    for(int i=0; i<np; i++) p[i].vel += dt*p[i].acc;

    // and drift the positions
    for(int i=0; i<np; i++) p[i].pos += dt*p[i].vel;

    cout << "done pseudo-nbody\n";
}

#endif
