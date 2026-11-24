#ifndef MATHLIBCPP
#define MATHLIBCPP

#include "mathlib.h"

namespace { 
// abscissae (x) and weights (w) for classic 8-point Gauss-Legendre integration
static double constexpr xGL[] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
static double constexpr wGL[] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
} // namespace


// 8-point Gauss-Legendre on [a,b]
double math::gauss_legendre_8_point(double const a, double const b, std::function<double(double)> const &func){
	double const mid = 0.5*(a+b), half = 0.5*(b-a);
	double sum = 0.0;
	for(int i=0; i<4; i++){
		sum += wGL[i]*(func(mid + half * xGL[i]) + func(mid - half*xGL[i]));
	}
	return half*sum;
}

#endif