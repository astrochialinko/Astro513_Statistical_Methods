#ifndef __TYPE_NAME_CC__
#define __TYPE_NAME_CC__

#include <string_view>

template <typename T>
constexpr auto type_name() {
  std::string_view name, prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void)";
#endif

  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

#endif // __TYPE_NAME_CC__

#if 0
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <array>

template <typename T>
struct MyStruct {

    T a, b, c;
    int k;

    static constexpr auto foo = type_name<T>();

    static constexpr std::array<const std::basic_string_view<char>, 2> goo = { type_name<decltype(a)>(), type_name<decltype(b)>() };
};



int main() {

    using U = double;

    MyStruct<U> x = MyStruct<U>();

    std::cout << x.foo << "\n";

    std::string l(x.goo[0]);
    for(int i=1; i<x.goo.size(); i++ ) l += ", " + std::string(x.goo[i]);
    std::cout << l << "\n";

    auto goo = type_name<decltype(x.goo)>();

    std::cout << goo << "\n";


}
#endif
