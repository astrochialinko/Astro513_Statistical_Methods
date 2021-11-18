/*
  Time-stamp: <range.cc on Tuesday, 2 October, 2018 at 13:15:55 MST (pinto)>

  Range provides an iterator similar to Python's range.

  Typical usage:
     for( auto var : Range(begin, end, step) ) ...
     for( auto var : Range(begin, end) ) ...  step defaults to one
     for( auto var : Range(end) ) ...  begin defaults to zero, step defaults to one

 */

#include <cassert>
#include <type_traits>

#ifndef __RANGE_INCLUDE__
#define  __RANGE_INCLUDE__

template <class T>
class Range {
public:
  Range(const T first, const T last, const T step) : first_(first), last_(last), step_(step) {
      static_assert(std::is_integral<T>::value, "Range: only integral types allowed in range");
      assert( step != 0 );
      assert( ( (step<0) && ( last <= first ) ) || ( (step>0) && (first <= last) ) );
  }

    class iterator {
    public:
        iterator(const Range &r, const T val) : owner_(r), val_(val) {}
        T operator*() const { return val_; }
        bool operator!=(const iterator &other) const { return val_ != other.val_; }
        iterator operator++() {
            val_ += owner_.step_;
            return *this;
        }

    private:
        const Range &owner_;
        T val_;
    };

    iterator begin() const { return iterator(*this, first_); }
    iterator end() const { return iterator(*this, last_); }

private:
    T first_, last_, step_;
};

template <class T>
Range<T> range(T first, T last, T step) {
  return Range<T>(first, last, step);
}

template <class T>
Range<T> range(T first, T last) {
  return Range<T>(first, last, 1);
}

template <class T>
Range<T> range(T last) {
  return Range<T>(0, last, 1);
}

//----------------------------------------------------------------------------------------------------
template <typename T>
class mrange {
public:

    class iterator {
    public:
        T operator *() const { return i_; }

        bool operator !=(const iterator &other) const { return i_[0] != other.i_[0]; }

        const iterator &operator ++() {
            for(int d=D-1; d>=0; d--) {
                i_[d]+=step_[d];
                if( d>0 && i_[d] == hi_[d] ) {
                    i_[d] = lo_[d];
                }
                else
                    break;
            }
            return *this;
        }

        iterator(T start, T end, T step) : i_(start), lo_(start), hi_(end), step_(step), D(end.size()) {}

    private:
        T i_;
        T lo_, hi_, step_;
        int D;
    };

    mrange(T begin, T end, T step) : begin_(begin, end, step), end_(end, end, step) { }
    iterator begin() const { return begin_; }
    iterator end() const { return end_; }

private:
    iterator begin_;
    iterator end_;
};
#endif //  __RANGE_INCLUDE__

//----------------------------------------------------------------------------------------------------

#if 0
#include <vector>
#include "smallvec.cc"
#include "pprint.cc"

int main() {

    for( auto i : range(10) ) pprint("%d ", i); pprint("\n");
    pprint("\n");
    for( auto i : range(0, 10) ) pprint("%d ", i); pprint("\n");
    for( auto i : range(0, 10, 2) ) pprint("%d ", i); pprint("\n");
    pprint("\n");
    for( auto i : range(10, 0, -1) ) pprint("%d ", i); pprint("\n");
    for( auto i : range(10, 0, -2) ) pprint("%d ", i); pprint("\n");
    pprint("\n");

    //for( auto i : range(0.0, 1.0, 0.5) ) pprint("%e ", i); pprint("\n");

    {
        std::vector<int> low = {0,0,0}, hi = {4,4,4}, step = {2,2,2}, nstep{-2, -2, -2};
        for( auto xyz : mrange( low, hi, step ) ) pprint("(%d %d %d)", xyz[0], xyz[1], xyz[2]); pprint("\n");
        for( auto xyz : mrange( hi, low, nstep ) ) pprint("(%d %d %d)", xyz[0], xyz[1], xyz[2]); pprint("\n");
        pprint("\n");
        //        for( auto xyz : mrange(low, hi) ) pprint("(%d %d %d)", xyz[0], xyz[1], xyz[2]); pprint("\n");

    }
    using ivec = SmallVec<int, 3>;
    for( auto xyz : mrange( ivec(0), ivec(4,4,4), ivec(2,1,2) ) ) pprint("(%d %d %d)", xyz[0], xyz[1], xyz[2]); pprint("\n");
    for( auto xyz : mrange( ivec(4,0,4), ivec(0,4,0), ivec(-2,1,-2) ) ) pprint("(%d %d %d)", xyz[0], xyz[1], xyz[2]); pprint("\n");
    pprint("\n");

}
#endif
