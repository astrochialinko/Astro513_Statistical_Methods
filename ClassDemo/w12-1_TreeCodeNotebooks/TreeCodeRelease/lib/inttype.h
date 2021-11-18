#ifndef __INTTYPE_CC__
#define __INTTYPE_CC__

// portable integer classes:
// This defines integer types in terms of clang and g++
// and possibly defines a 128 bit type for use on 32-bit machines

#include <cstdint>

#ifdef ARCH32
   #include "uint128_t.cc"
   using uint128 = uint128_t;
#else
   using uint128 = __uint128_t;
#endif

using uint64 = uint64_t;
using uint32 = uint32_t;
using uint64 = uint64_t;
using uint32 = uint32_t;
using uint16 = uint16_t;
using uint8  = uint8_t;

#endif // __INTTYPE_CC__
