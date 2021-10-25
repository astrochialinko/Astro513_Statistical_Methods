#include <string>
#include <cstdarg>
#include <vector>

// Produce a std::string from C-style printf call
// requires at least C++11
const std::string vformat(const char * const zcFormat, ...) {

    // initialize use of the variable argument array
    va_list vaArgs;
    va_start(vaArgs, zcFormat);

    // reliably acquire the size
    // from a copy of the variable argument array
    // and a functionally reliable call to mock the formatting
    va_list vaArgsCopy;
    va_copy(vaArgsCopy, vaArgs);
    const int iLen = std::vsnprintf(NULL, 0, zcFormat, vaArgsCopy);
    va_end(vaArgsCopy);

    // return a formatted string without risking memory mismanagement
    // and without assuming any compiler or platform specific behavior
    std::vector<char> zc(iLen + 1);
    std::vsnprintf(zc.data(), zc.size(), zcFormat, vaArgs);
    va_end(vaArgs);
    return std::string(zc.data(), iLen); }

#if 0
#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdio>

// demonstration of use
int main() {

    std::cout << vformat("Int 1 is %d, Int 2 is %d, Int 3 is %d\n", 11, 22, 33);

    printf("%s", vformat("Int 1 is %d, Int 2 is %d, Int 3 is %d\n", 11, 22, 33).c_str() );

    return 0;
}
#endif
