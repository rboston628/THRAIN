#ifndef MATHLIBH
#define MATHLIBH

#include <functional>

namespace math {

// good ol' pi!
double constexpr PI = 3.1415926535897932384626433832795;

// 8-point Gauss-Legendre on [a,b]
double gauss_legendre_8_point(double const a, double const b, std::function<double(double)> const &func);

} // namespace math


#endif