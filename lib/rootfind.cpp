#ifndef ROOTFINDINGMETHODS_C
#define ROOTFINDINGMETHODS_C

#include "rootfind.h"
#include "logger.h"
#include <complex>
#include <stdio.h>

double rootfind::pseudo_unif(){
	static int a=53, b=122, r=17737, dart = 39;//for pseudorandom positioning
	dart = (a*dart+b)%r; //generates a psuedo-random integer in (0,r)
	return (double(dart)/double(r));
}

double rootfind::pseudo_unif(double xmin, double xmax){
	double dart = xmin;
	if (xmin < xmax){
		dart = xmin + rootfind::pseudo_unif() * (xmax - xmin);
	}
	else if (xmin > xmax){
		dart = xmax + rootfind::pseudo_unif() * (xmin - xmax);
	}
	return dart;
}

namespace {
double left(double x0, double dx) {
	if (x0 > 0.0) return (1. - dx) * x0;
	if (x0 < 0.0) return (1. + dx) * x0;
	return -dx;
}
double right(double x0, double dx) {
	if (x0 > 0.0) return (1. + dx) * x0;
	if (x0 < 0.0) return (1. - dx) * x0;
	return dx;
}
int fail(double &xmin, double &xmax) {
	xmin = nan("");
	xmax = nan("");
	ThrainLogger::error("Unable to find a root near this guess\n");
	return 1;
}
double fail(double &x) {
	x = nan("");
	ThrainLogger::error("Unable to find a root near this guess\n");
	return nan("");
}

enum class Slope {NONE, UP, DOWN, FLAT, CONCAVE_UP, CONCAVE_DOWN};
Slope determine_behavior(double yleft, double ymid, double yright) {
	// in ideal case, 
	// either y1 < y < y2 (UP)
	// or     y1 > y > y2 (DOWN)
	// in non-ideal cases 
	// could be y1 = y = y2 (FLAT)
	// could be y>y1 and y>y2, (CONCAVE)
	// or       y<y1 and y<y2
	if( (yleft==ymid) && (yright==ymid) ) return Slope::FLAT;
	else if( (yleft<=ymid && ymid<yright) || (yleft<ymid && ymid<=yright) ) return Slope::UP;
	else if( (yleft>=ymid && ymid>yright) || (yleft>ymid && ymid>=yright) ) return Slope::DOWN;
	else if( (yleft< ymid && yright< ymid) ) return Slope::CONCAVE_DOWN;
	else if( (yleft>ymid && yright> ymid) ) return Slope::CONCAVE_UP;
	return Slope::NONE;
}

enum class MoveDirection : char {LEFT = 0, RIGHT = 1, INDETERMINATE = -1};
MoveDirection move_direction(double yleft, double ymid, double yright) {
	// determine the slope of the function around ymid
	Slope behavior = determine_behavior(yleft, ymid, yright);
	
	// now decide which direction to move the brackets based on slope of function
	bool move_left = false;
	if (ymid > 0.0) {
		switch(behavior) {
			case Slope::UP: move_left = true; break;
			case Slope::DOWN: move_left = false; break;
			case Slope::CONCAVE_DOWN: move_left = false; break;
			default: return MoveDirection::INDETERMINATE;
		}
	} else if (ymid < 0.0) {
		switch(behavior) {
			case Slope::UP: move_left = false; break;
			case Slope::DOWN: move_left = true; break;
			case Slope::CONCAVE_UP: move_left = false; break;
			default: return MoveDirection::INDETERMINATE;
		}
	}
	return move_left ? MoveDirection::LEFT : MoveDirection::RIGHT;
}
} // namespace

// *************************************************************************************
//					BISECTION METHODS
//  These functions implement bisection searches in one parameter
// *************************************************************************************

// this will find brackets using a simple expansion/movement search
// one bracket is moved to one side until the function changes sign, to find brackets
int rootfind::bisection_find_brackets_move(
	std::function<double(double)> const& func,    //the function to find zero of
	double const x0,                        //an initial guess for the zero
	double &xmin,                           //the lower bracket -- will be returned
	double &xmax                            //the upper bracket -- will be returned
){
	static std::size_t constexpr MAX_ITER = 100; 
	// test two points near first guess
	double ymin, ymax;

	// happy path 1
	// the initial guess is the zero
	double y = func(x0);
	if(y==0.0){
		xmin = x0;
		xmax = x0;
		return 0;
	}

	// test two points near first guess
	double DX = 0.01;
	double x1=left(x0, DX), x2=right(x0, DX);
	double y1=func(x1), y2=func(x2);

	// unhappy path: any evaluations are non-finite
	if (!std::isfinite(y) || !std::isfinite(y1) || !std::isfinite(y2)) return fail(xmin, xmax);

	// happy path 2
	// either y1 or y2 are brackets
	if(y1*y < 0.0) { // if y1 and y are on opposite side
		xmin = x1;
		xmax = x0;
		return 0;
	}
	if(y2*y < 0.0) { // if y2 and y are on opposite sides
		xmax = x2;
		xmin = x0;
		return 0;
	}

	MoveDirection move = move_direction(y1, y, y2);
	if (move == MoveDirection::INDETERMINATE) {
		return fail(xmin, xmax);
	}
	bool move_left = (move == MoveDirection::LEFT);

	// move one bracket until the function changes sign
	double xm = (move_left ? x1 : x2), ym = (move_left ? y1 : y2);
	for (std::size_t i = 0; ym * y > 0.0; ++i) {
		if (i >= MAX_ITER) return fail(xmin, xmax);
		DX *= 2.0;
		xm = (move_left ? left(x0, DX) : right(x0, DX));
		ym = func(xm);
		if (!std::isfinite(ym)) return fail(xmin, xmax);
	}

	// make sure in order
	if (move_left) {
		xmin = xm;
		xmax = x0;
	} else {
		xmin = x0;
		xmax = xm;
	}
	return 0;
}

// this will find brackets using newton's method to look for the next zero, then a bit beyond it
// why use one zero-finding method to prepare another zero-finding method? 
//    because Newton's method can fail, but bisection searches can go to nearly arbitrary accuracy
int rootfind::bisection_find_brackets_newton(
	std::function<double(double)> const& func, //the function to find zero of
	double const x0,                     //an initial guess for the zero
	double& xmin,                        //the lower bracket -- will be returned
	double& xmax                         //the upper bracket -- will be returned
){
	static std::size_t constexpr MAX_ITER = 100; 
	double y = func(x0);

	// unhappy path: initial guess is non-finite
	if (!std::isfinite(y)) return fail(xmin, xmax);
	
	// happy path: x0 is already a zero
	if (y == 0.0 ) {
		xmin = left(x0, 0.01);
		xmax = right(x0, 0.01);
		return 0;
	}

	// sample on either side of the initial guess
	double x1 = left(x0, 0.01), x2 = right(x0, 0.01);
	double y1 = func(x1), y2 = func(x2);

	// happy path: the left and right form a bracket
	if (y1 * y2 < 0.0) {
		xmin = x1;
		xmax = x2;
		return 0;
	}

	// determine whether to move left or right
	MoveDirection move = move_direction(y1, y, y2);
	if (move == MoveDirection::INDETERMINATE) {
		return fail(xmin, xmax);
	}
	bool move_left = (move == MoveDirection::LEFT);

	//while the two ys are on same side of axis, keep reposition until zero is bound
	//we use Newton's method to the nearest zero
	double xm = (move_left ? x1 : x2);
	double ym = (move_left ? y1 : y2);
	if (!std::isfinite(ym)) return fail(xmin, xmax);
	int iter = 0;
	while(ym * y >= 0.0){
		iter++;
		if(iter > MAX_ITER) return fail(xmin, xmax);
		//compute numerical derivative
		double xdx = left(xm, 0.01);
		double ydx = func(xdx);
		xdx = xm - ym * (xdx - xm) / (ydx - ym);
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( xdx == xm ) {
			//if we are moving left, we want to move further left
			//if we are moving right, we want to move further right
			if (move_left) xdx = left(xdx, 0.1);
			else xdx = right(xdx, 0.1);
		}
		//limit amount of change allowed in single step
		double maxstep = 0.2 * std::max(fabs(xm), fabs(x0));
		if (move_left) xdx = std::max(xdx, xm - maxstep);
		else xdx = std::min(xdx, xm + maxstep);
		ym = func(xdx);
		xm = xdx;
	}
	// now bracket found, order correctly
	if(move_left) {
		xmax = x2;
		xmin = xm; 
	}
	else {
		xmax = xm;
		xmin = x1;
	}
	return 0;
}

//given brackets bounding a single zero, find the zero
// the brackets xmin, xmax MUST bound a single zero
double rootfind::bisection_search(
	std::function<double(double)> const& func, //the function to find zero of
	double &x,                           //the location of zero -- will be returned
	double &xmin,                        //the lower bracket
	double &xmax                         //the upper bracket
){
	static std::size_t constexpr MAX_ITER = 100;
	static std::size_t constexpr MAX_BRACKET_TRIES = 3;

	// if the brackets do not enclose a root, find better brackets
	double ymin = func(xmin), ymax=func(xmax);
	if (ymin * ymax > 0){
		double xlo = xmin, xhi = xmax;
		for(std::size_t t = 0; t < MAX_BRACKET_TRIES; t++){
			int status = bisection_find_brackets_newton(func, x, xlo, xhi);
			if (status) status = bisection_find_brackets_move(func, x, xlo, xhi);
			if (status) break;
			ymin = func(xlo);
			ymax = func(xhi);
			// for later tries, perturb the guess
			x = pseudo_unif(left(x, 0.5), right(x, 0.5));
		}
		if (ymin * ymax > 0) return fail(x);
		xmin = xlo;
		xmax = xhi;
	}

	// happy path: either side is a root
	if (ymin == 0.0) {
		x = xmin;
		return 0.0;
	}
	if (ymax == 0.0) {
		x = xmax;
		return 0.0;
	}

	// now begin bisection search
	double y = ymin;
	for (std::size_t iter = 0; iter < MAX_ITER; iter++) {
		x = 0.5 * (xmin + xmax);
		y = func(x);

		// handle possible non-finite evaluation of function at x
		if (!std::isfinite(y)) {
			for (int retry = 0; retry < MAX_BRACKET_TRIES && !std::isfinite(y); retry++) {
				x = pseudo_unif(xmin, xmax);
				y = func(x);
			}
			if (!std::isfinite(y)) return fail(x);
		}

		// we found it!
		if (y == 0.0) return 0.0;

		// otherwise, keep going
		if (y * ymin > 0.0) {
			xmin = x;
			ymin = y;
		} else {
			xmax = x;
			ymax = y;
		}
	}
	// return the midpoint of the final brackets
	x = 0.5 * (xmin + xmax);
	return func(x);
}


#endif