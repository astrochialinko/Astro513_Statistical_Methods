#ifndef __AOS_CC__
#define __AOS_CC__

#include <tuple>
#include <iostream>
#include <fstream>

#include "ptype.cc"

//--------------------------------------------------------------------------------
//
//  Array-of-structs implemented as struct-of-arrays
//
//  To use, define a struct, e.g. in ptype.in, and then run mkAOStype on that file.
//  This will create the type P used to define this class.
//
//--------------------------------------------------------------------------------
template <class P>
class AOS {
public:

    AOS( ) : N(0) { }
    AOS( int N ) : N(0) { reserve(N); }
    ~AOS() {  free();  }

    // release memory
    void free() {
        if( N>0 ) {
            freeAll(base);
            N = 0;
        }
    }

    // allocate memory
    void reserve(int nNew) {
        if(N>0) free();
        N = nNew;
        if(N>0) reserveAll(base, N);
    }

    int size() {
        return N;
    }

    // allow referring to elements with AOS notation:
    // if a is a member of class P, then it can be referred to as
    //      pset[i].a
    inline typename P::reftype operator[](int index) {
        return { elementsAll(base, index) };
    }

    // set all of the data at index to a class P instance Q:
    // Since pset[i] is P::reftype, it can be used as is:
    //      pset[i] = pset[j]
    //      P A = pset[10];
    //      pset[12] = A
    inline void set(int index, P& Q) {
        (*this)[index].set(Q);
    }

    // swap data at positions i and j
    inline void swap(int i, int j) {
        typename P::tuptype pJ = (*this)(j);
        (*this)(j) = (*this)(i);
        (*this)(i) = pJ;
    }

    // allow referring to entire arrays in SOA notation:
    // If a is a member of class P, then pset.arrays().a is
    // a pointer to the underlying array.
    // E.g., pset.arrays().a[i] is the same as pset[i].a
    inline typename P::ptrtype arrays() {
        return arraysAll(base);
    }

#if 0
    // allow axccessing entire P instances:
    // for example,  pset(i) = pset(j)
    // returns a reftuple at index
    auto operator()(int index) {
        return tupAll(base, index);
    }
#endif

    void print(int index) {
        print(tupAll(base, index));
    }

    // write the entire AOS instance to the stream
    void write(std::ofstream &stream) {
        writeAll(base, N, stream);
    }

    // read the entire AOS instance to the stream
    // NB: one must have reserved a length of N
    void read(std::ifstream &stream) {
        readAll(base, N, stream);
    }

    int N = 0;
    typename P::basetype base;


//private:

    //--------------------------------------------------------------------------------
    // allocate np elements of type T to the array given from tuple of arrays
    template <class T>
    void reserveOne(T** base1, int np ) {
        *base1 = new T[np];
    }

    template<typename T, std::size_t... Is>
    void reserve_impl( T& t, int np, std::index_sequence<Is...> ) {
        auto x = { (reserveOne(&(std::get<Is>(t)), np),0)... };
    }

    template<class T>
    void reserveAll( T &t, int index ) {
        constexpr auto size = std::tuple_size<T>{};
        reserve_impl(t, index, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // free storage from arrays in tuple of arrays
    template<typename T, std::size_t... Is>
    void free_impl( T& t, std::index_sequence<Is...> ) {
        auto x = { (delete[] std::get<Is>(t),0)... };
    }

    template<class T>
    void freeAll( T &t ) {
        constexpr auto size = std::tuple_size<T>{};
        free_impl(t, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // get std::initializer_list of arrays from tuple of arrays
    template<typename T, std::size_t... Is>
    inline typename P::ptrtype arraysImpl(const T& t, std::index_sequence<Is...>) {
        return { (std::get<Is>(t))... };
    }

    template<class T>
    inline typename P::ptrtype arraysAll(const T &t) {
        constexpr auto size = std::tuple_size<T>{};
    return arraysImpl(t, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // get std::initializer_list of array elements at index from tuple of arrays
    template<typename T, std::size_t... Is>
    inline typename P::reftype elementsImpl(const T& t, int index, std::index_sequence<Is...>) {
        return { (std::get<Is>(t))[index]... };
    }

    template<class T>
    inline typename P::reftype elementsAll(const T &t, int index) {
        constexpr auto size = std::tuple_size<T>{};
        return elementsImpl(t, index, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // get reftuple of array elements at index from tuple of arrays
    template<typename T, std::size_t... Is>
    inline typename P::tupref tupImpl(const T& t, int index, std::index_sequence<Is...>) {
        return { (std::get<Is>(t))[index]... };
    }

    template<class T>
    inline typename P::tupref tupAll(const T &t, int index) {
        constexpr auto size = std::tuple_size<T>{};
        return tupImpl(t, index, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // print the contents of a tuple
    template<class TupType, size_t... I>
    void print(const TupType& _tup, std::index_sequence<I...>) {
        std::cout << "[ ";
        (..., (std::cout << (I == 0? "" : ", ") << std::get<I>(_tup)));
        std::cout << " ]\n";
    }

    template<class... T>
    void print (const std::tuple<T...>& _tup) {
        print(_tup, std::make_index_sequence<sizeof...(T)>());
    }

    //--------------------------------------------------------------------------------
    // read np elements into one array in a tuple of arrays
    template <class T>
    void readOne(T a, int np, std::ifstream& istream) {
        using type = typename std::remove_pointer<T>::type;
        istream.read((char*)(a), np*sizeof(type));
    }

    template<typename T, std::size_t... Is>
    void read_impl( T& t, int np, std::ifstream& istream, std::index_sequence<Is...>) {
        auto x = { (readOne((std::get<Is>(t)), np, istream),0)... };
    }

    template<class T>
    void readAll( T &t, int np, std::ifstream& istream) {
        constexpr auto size = std::tuple_size<T>{};
        read_impl(t, np, istream, std::make_index_sequence<size>{});
    }

    //--------------------------------------------------------------------------------
    // write np elements from one array in a tuple of arrays
    template <class T>
    void writeOne(T a, int np, std::ofstream& ostream) {
        using type = typename std::remove_pointer<T>::type;
        ostream.write((char*)(a), np*sizeof(type));
    }

    template<typename T, std::size_t... Is>
    void write_impl( T& t, int np, std::ofstream& ostream, std::index_sequence<Is...>) {
        auto x = { (writeOne((std::get<Is>(t)), np, ostream),0)... };
    }

    template<class T>
    void writeAll( T &t, int np, std::ofstream& ostream) {
        constexpr auto size = std::tuple_size<T>{};
        write_impl(t, np, ostream, std::make_index_sequence<size>{});
    }

};

#endif // __AOS_CC__


#ifdef TEST_AOS

/*
  Compile with:
     g++ -Ilib -std=c++2a -Ofast -DTEST_AOS AOS.cc 2>&1 | less
 */

#include "STimer.cc"
#include "util.cc"

template <typename T>
void printP(const PtypeRef<T> &p) {
    pprint("%4d  % .4e % .4e % .4e %6.2f\n", p.id,  p.r, p.v, p.a, p.mass);
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

template <class T>
bool equal( const PtypeRef<T> &a, const PtypeRef<T> &b ) {
    bool result = true;
    result &= ( a.r == b.r );
    result &= ( a.v == b.v );
    result &= ( a.a == b.a );
    result &= ( a.pot == b.pot );
    result &= ( a.mass == b.mass );
    result &= ( a.id == b.id );

    if(!result) {
        pprint("ee % e % e % e % e % e %4d\n", a.r, a.v, a.a, a.pot, a.mass, a.id);
        pprint("ee % e % e % e % e % e %4d\n", b.r, b.v, b.a, b.pot, b.mass, b.id);
    }

    return result;
}

template <class T>
bool equal( const PtypeRef<T> &a, const Ptype<T> &b ) {
    bool result = true;
    result &= ( a.r == b.r );
    result &= ( a.v == b.v );
    result &= ( a.a == b.a );
    result &= ( a.pot == b.pot );
    result &= ( a.mass == b.mass );
    result &= ( a.id == b.id );

    if(!result) {
        pprint("eq % e % e % e % e % e %4d\n", a.r, a.v, a.a, a.pot, a.mass, a.id);
        pprint("eq % e % e % e % e % e %4d\n", b.r, b.v, b.a, b.pot, b.mass, b.id);
    }

    return result;
}

int main() {

    using Real = double;
    using ptype = Ptype<Real>;
    using vec = SmallVec<Real, 3>;

    STimer timeit;

    uint N = 10000;

#if 1
    // will need to do this for the Python wrapper
    AOS<ptype> s = AOS<ptype>(0);
    vec *rr = new vec[N];
    vec *vv = new vec[N];
    vec *aa = new vec[N];
    double *ppot = new double[N];
    double *mmass = new double[N];
    int *iid = new int[N];

    std::get<0>(s.base) = rr;
    std::get<1>(s.base) = vv;
    std::get<2>(s.base) = aa;
    std::get<3>(s.base) = ppot;
    std::get<4>(s.base) = mmass;
    std::get<5>(s.base) = iid;
    s.N = N;
#else

    AOS<ptype> s = AOS<ptype>(N);

#endif


    // access as AOS
    timeit.CStart();
    for(int i=0; i<N/2; i++) {
        Ptype<Real> Q = func( i, (Real)(1.0) );
        s[i].r    = Q.r;
        s[i].v    = Q.v;
        s[i].a    = Q.a;
        s[i].pot  = Q.pot;
        s[i].mass = Q.mass;
        s[i].id   = Q.id;
    }
    timeit.Stop();
    timeit.Report("AOS style:");

    // Access as SOA:
    vec *r = s.arrays().r;
    vec *v = s.arrays().v;
    vec *a = s.arrays().a;

    timeit.CStart();
    for(int i=N/2; i<N; i++) {
        Ptype<Real> Q = func(i, (Real)(1.0) );
        r[i] = Q.r;
        v[i] = Q.v;
        a[i] = Q.a;
        s.arrays().pot[i] = Q.pot;
        s.arrays().mass[i] = Q.mass;
        s.arrays().id[i] = Q.id;
    }
    timeit.Stop();
    timeit.Report("SOA style:");

    // Write out data
    {
        std::ofstream file("outfoo");
        s.write(file);
        file.close();
    }

    // make a copy
    AOS<ptype> sc = AOS<ptype>();
    sc.reserve(N);

    for(int i=0; i<s.N; i++) sc[i] = s[i];

    // compare
    pprint("copy compare: ");
    for(int i=0; i<s.N; i++) {
        assert( equal(s[i], sc[i]) );
        Ptype<Real> Q = func( i, (Real)(1.0) );
        assert( equal(s[i], Q) );
    }
    printf("OK\n");

    // put garbage in original
    for(int i=0; i<s.N; i++) {
        s[i].r.zero();
        s[i].v    = {0.0, 0.0, 0.0};
        s[i].a    = vec(0.0);
        s[i].pot  = 0;
        s[i].mass = 0;
        s[i].id   = 0;
    }

    // read back in
    {
        std::ifstream file("outfoo");
        s.read(file);
        file.close();
    }

    // compare
    pprint("read compare: ");
    for(int i=0; i<s.N; i++) {
        assert( equal(s[i], sc[i]) );
        Ptype<Real> Q = func( i, (Real)(1.0) );
        assert( equal(s[i], Q) );
    }
    printf("OK\n");

    // exercise set
    s.free();
    s.reserve(sc.N);
    for(int i=0; i<N; i++) {
        Ptype<Real> Q = func( i, (Real)(1.0) );
        s.set(i, Q);
    }

    // compare
    pprint("set compare: ");
    for(int i=0; i<s.N; i++) {
        assert( equal(s[i], sc[i]) );
        Ptype<Real> Q = func( i, (Real)(1.0) );
        assert( equal(s[i], Q) );
    }
    printf("OK\n");

    for(int i=0; i<10; i++) {
        s.print(i); pprint("\n");
    }

}

#endif // TEST_AOS
