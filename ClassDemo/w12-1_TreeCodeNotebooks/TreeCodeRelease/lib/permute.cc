/*
Time-stamp: <permute.cc on Monday, 15 November, 2021 at 23:12:30 MST (pinto)>

Sorting and permuting functions

compile with -D_GLIBCXX_PARALLEL to enable parallel std::sort
 */

#ifndef __PERMUTE_CC__
#define  __PERMUTE_CC__

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <algorithm>
#include <omp.h>

/*
  Return the permutation which sorts the given data, which is left unaltered.
*/
template <typename T, typename I, typename Func>
void sortWithPermutation(T *data, I *perm, I len, Func ascending) {
    std::iota( perm, perm + len, 0);
    std::sort( perm, perm + len, ascending );
}

template <typename T, typename I>
void sortWithPermutation(T *data, I *perm, I len) {
    std::iota( perm, perm + len, 0);
    std::sort( perm, perm + len,
               [&](const I &a, const I &b) {
                   return data[a] < data[b];
               }
        );
}

/*
  In-situ permutation: slower but O(1) extra space.
  The data is permuted according to index; the index is returned to iota(0,n) (i.e. 0,1,2,...,n-1)
*/
template <typename T, typename I>
void inSituPermuteDestroyPermutation(T *data, I *perm, I len) {
    I i,j,k;
    T q;

    for( i=0; i<len; i++ ) {
        q = data[i];
        for(k=i; perm[k]!=i; k=perm[j], perm[j]=j) {
            j = k;
            data[k] = data[perm[k]];
        }
        data[k] = q;
        perm[k] = k;
    }
}

/*
  In-situ permutation: slower but O(1) extra space.
  The data is permuted according to origPerm, which is preserved intact
*/
template <typename T, typename I>
void inSituPermute(T *data, I *origPerm, I len) {
    I *perm = new I[len];
    for(int i=0; i<len; i++) perm[i] = origPerm[i];
    inSituPermuteDestroyPermutation(data, perm, len);
    delete[] perm;
}

/*
  direct permutation; requires O(N) extra space.
*/
template <typename T, typename I>
void permute(T *data, I *perm, I len) {
    T *tmp = new T[len];
    for(int i=0; i<len; i++) tmp[i] = data[i];
    for(int i=0; i<len; i++) data[i] = tmp[perm[i]];
    delete[] tmp;
}

#endif  // __PERMUTE_CC__
