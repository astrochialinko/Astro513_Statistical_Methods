/*
  To test: g++ -D__TEST_PPRINT__ pprint.cc
 */

/*
  C++ replacement for printf

  Time-stamp: <pprint.cc on Thursday, 11 November, 2021 at 16:39:35 MST (pinto)>

  The format specification is just like normal printf, but for a few
  peculiarities arising from the fact that the format does not <<need>> to tell
  the function about the type of the data in the argument list, since the actual
  printing mechanism is the << operator of c++.  Any object which overloads <<
  can be given as an argument. If that argument does not override the previous
  format state (i.e., std::iomanip), then the format will apply to the object.

  As in printf, % followed by any single character (except another %) signals
  printing the next item in the argument list.

  % followed by an integer followed by a single character signifies a width (in
  characters) which is the <<minimum>> width in which to formate the next item
  in the argument list (if the item takes more space to format naturally, then
  that length is used instead).

  Which character follows % only matters in four instances:
    %f -> fixed-point format
    %e -> scientific notation
    %x -> hex
    %o -> octal

  %e and %f can be used alone (%e or %f), with a width (e.g. %10f, %15e), or
  with a width and a precision (e.g. %10.4f, %15.7e). Any format can take a
  precision, but it is ignored if the format is not %f or %e.

  A "-" sign following % means "left-justify the next item."

  A "+" following % (or following a "-" sign) means "show the plus sign
  for positive numbers (for numeric data)."

  Of course, if the width isn't enough for a given format and precision, then the
  intrinsic with is used; thus %0.7e will format a value in scientific notation
  with seven digits after the decimal point, and whatever width is required.

  Thus, a format is one of:

  %[-][+][ ][efox...]
  %[-][+][ ][0-9]+[efox...]
  %[-][+][ ][0-9]+.[0-9]+[efox...]

  or any character other than % instead of e, f, o, or x if one doesn't care about numeric format

  See the test code at the end of this file for examples; run it compiled with -D__TESTPPRINT__

  NOTE:
     This code relies on variadic templates which are defined in the new C++ standard.
     You will likely need to use the -std=c++0x flag to g++ when compiling with pprint.


 */

#ifndef __PPRINT_CC_INCLUDED__
#define __PPRINT_CC_INCLUDED__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <cstdlib>

#define RESET                                           \
    out.unsetf(std::ios_base::floatfield);              \
    out.unsetf(std::ios_base::showpos);                 \
    out.unsetf(std::ios_base::showbase);                \
    out.unsetf(std::ios_base::left);                    \
    out << std::setbase(10) << std::setprecision(6);

//---------------------------------
#include <type_traits>

namespace detail
{

template <typename T, typename NameGetter>
struct has_member_impl
{
    typedef char matched_return_type;
    typedef long unmatched_return_type;

    template <typename C>
    static matched_return_type f(typename NameGetter::template get<C>*);

    template <typename C>
    static unmatched_return_type f(...);

public:
    static const bool result = (sizeof(f<T>(0)) == sizeof(matched_return_type));
};

}

template <typename T, typename NameGetter>
struct has_member :
        std::integral_constant<bool, detail::has_member_impl<T, NameGetter>::result>
{ };

struct check_has_print
{
    template <typename T,
              void (T::*)(std::ostream&, std::string&) const = &T::print
             >
    struct get
    { };
};

template <typename T>
struct has_print :
        has_member<T, check_has_print>
{ };

template <bool> struct printness {};
typedef printness<true> print_tag;
typedef printness<false> non_print_tag;

template <typename T>
void f2(std::ostream &o, T x, std::string &fmt, print_tag) { /* POD */
    x.print(o, fmt);
}

template <typename T>
void f2(std::ostream &o, T x, std::string &fmt, non_print_tag) { /* non-POD */
    o << x;
}

#if 0
template <typename T>
void f(T x) {
    // Dispatch to f2 based on tag.
    f2(x, printness<has_print<T>::value>());
}
#endif

//---------------------------------




// set width, precision, and format of fields using std::iomanip
// only %e and %f specify a particular format -- anything else will print anything!
// argument is a copy of the rest of the format line following the current %
// Returns length of format specifier to be gobbled up in this line.
template <typename IO>
size_t pprint_parse_flags(const char *s, IO &out, std::string &fmt) {
    bool ok = true;
    int width = 0;
    int precision = 6;
    long flags = 0;
    char fill = ' ';
    bool alternate = false;
    int extra = 0;

    int i = 0;
    bool more = true;
    while (more)
    {  switch (s[i])
        {
        case '+':
            flags |= std::ios::showpos;
            break;
        case '-':
            flags |= std::ios::left;
            break;
        case '0':
            flags |= std::ios::internal;
            fill = '0';
            break;
        case '#':
            alternate = true;
            break;
        case ' ':
            flags |= std::ios::right;
            extra = 1;
            break;
        default:
            more = false;
            break;
        }
        if (more) i++;
    }

    if (isdigit(s[i])) {
        width += atoi(s+i);
        do i++; while (isdigit(s[i]));
    }
    if (s[i] == '.') {
        i++;
        precision = atoi(s+i);
        while (isdigit(s[i])) i++;
    }

    switch (s[i]) {
    case 'd':
        flags |= std::ios::dec;
                     break;
    case 'x':
        flags |= std::ios::hex;
        if (alternate) flags |= std::ios::showbase;
        break;
    case 'X':
        flags |= std::ios::hex | std::ios::uppercase;
        if (alternate) flags |= std::ios::showbase;
        break;
    case 'o':
        flags |= std::ios::hex;
        if (alternate) flags |= std::ios::showbase;
        break;
    case 'f':
        flags |= std::ios::fixed;
        if (alternate) flags |= std::ios::showpoint;
        width += extra;
        break;
    case 'e':
        flags |= std::ios::scientific;
        if (alternate) flags |= std::ios::showpoint;
        //        if(width==0) width = std::max(12, precision+6);
        if(width==0) width = precision+6;
        width += extra;
        break;
    case 'E':
        flags |= std::ios::scientific | std::ios::uppercase;
        if (alternate) flags |= std::ios::showpoint;
        break;
    case 'g':
        if (alternate) flags |= std::ios::showpoint;
        break;
    case 'G':
        flags |= std::ios::uppercase;
        if (alternate) flags |= std::ios::showpoint;
        break;
    case 's':
        break;
    default:
        std::cout << "unknown: " << s[i] << std::endl;
        ok = false;
        break;
    }
    if(!ok)
        throw std::runtime_error("invalid format: unknown descriptor\n");

    out.unsetf(std::ios::adjustfield | std::ios::basefield | std::ios::floatfield);
    out.setf((std::ios_base::fmtflags) flags);
    out.width(width);
    out.precision(precision);
    out.fill(fill);

    int len = i;

    fmt = "%" + std::string(s).substr(0,len+1);

    return len;
}

// final recursion of pprint: print everything after the last format specifier
template <typename IO>
void fpprint(IO &out, const char* s) {
    while (s && *s) {
        if (*s=='%' && *++s!='%') throw std::runtime_error("invalid format: missing arguments\n");
        out << *s++;
    }
}

// recursive pprint: print up to the first format specifier, interpret the specifier,
// print the first value on the argument list, then make recursive call with remaining
// format string and remaining args
template<typename T, typename... Args>
void fpprint(std::ostream &out, const char* s, T value, Args... args) {
    std::string fmt;
    while (s && *s) {
        if (*s=='%' && *++s!='%') {          // a format specifier
            s += pprint_parse_flags(s, out, fmt);
            f2(out, value, fmt, printness<has_print<T>::value>());
            RESET;
            return fpprint(out, ++s, args...);
        }
        out << *s++;
    }
    throw std::runtime_error("extra arguments provided to printf\n");
}
#undef RESET

template<typename... Args>
void pprint(const char *s, Args... args) {
    fpprint(std::cout, s, args...);
}
void pprint(const char* s) {
    fpprint(std::cout, s);
}

template<typename... Args>
std::string spprint(const char *s, Args... args) {
    std::ostringstream os;
    fpprint(os, s, args...);
    return os.str();
}

#endif // __PPRINT_CC_INCLUDED__


#if 0
#include "Threevector.cc"
#include <sstream>

template <class T>
void xx(char *s, T x) {
    char c[1024];
    std::stringstream cpp;

    sprintf(c, s, x);
    fpprint(cpp, s, x);

    std::cout << c << std::endl;
    std::cout << cpp.str() << std::endl;

    if (std::string(c) != cpp.str() ) {
        std::cout << "FAIL\n";
        exit(1);
    }
}

int main() {

    int i = 42;
    double x = 0.1234567890e-02;
    ThreeVector<double> V(0.1, 0.2, 0.3);

#if 1

    xx("%.12e", x);
    xx("% .12e", x);
    xx("% .12e", -x);

    std::cout << "---------------------------------\n";

    xx("%e", x);
    xx("%12.3e", x);
    xx("%12e", x);
    xx("% e", x);
    xx("% e", -x);
    xx("% 9.3e", x);
    xx("% 9.3e", -x);

    std::cout << "---------------------------------\n";

    xx("%g", x);
    xx("%#14.11g", x);

    xx("%12.6g", 9.999999e+05);
    xx("%12.7g", 9.999999e+05);

    std::cout << "---------------------------------\n";

    xx("%f", x);
    xx("%5.2f", 19.05);
    xx("% 5.2f", 19.05);
    xx("% 5.2f", -19.05);

    std::cout << "---------------------------------\n";

    xx("%d", i);
    xx("%5d", i);
    xx("%05d", i);
    xx("%5d", -i);
    xx("% 5d", -i);
    //    xx("% d", i); // this fails -- how to fix w/o knowing the magnitude of the value?
    xx("% d", -i);

    std::cout << "---------------------------------\n";




    //    pprint("the values are: %d and %e and %10.3e so there!\n", i, x, V);
        pprint("%d %e %10.3g\n", i, x, V);
        pprint("%d %e %#10.3g\n", i, x, V);

    pprint("% e\n", x);
    pprint("% e\n", -x);
    printf("% e\n", x);
    printf("% e\n", -x);


    pprint("%04d %e %20.10e\n", i, x, x);
    pprint("%04d %e %20.10e\n", i, -x, -x);

    //    std::cout << std::scientific << x << std::endl;

    //    std::cout << std::is_fundamental<ThreeVector<double>>::value << std::endl;
    //    std::cout << std::is_fundamental<int>::value << std::endl;
#endif


        pprint("%d %e %10.3g\n", i, x, V);
        pprint("%d %e %#10.3g\n", i, x, V);

        pprint("%d %e\n", i, V);
        pprint("%d % e\n", i, V);
        pprint("%d % e %e\n", i, -V, V);

        pprint("%e\n", V);

        for(int i=0; i<5; i++) {
            double3 Q(drand48()-0.5, drand48()-0.5, drand48()-0.5);
            pprint("% .2e\n", Q);
        }

}
#endif
