#ifndef ROOTFINDINGMETHODS_C
#define ROOTFINDINGMETHODS_C

#include "rootfind.h"
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

// *************************************************************************************
//					BISECTION METHODS
//  These functions implement bisection searches in one parameter
// *************************************************************************************

// this will find brackets using a simple expansion/movement search
// one bracket is moved to one side until the function changes sign, to find brackets
int rootfind::bisection_find_brackets_move(
	std::function<double(double)>& func,    //the function to find zero of
	double const x0,                        //an initial guess for the zero
	double &xmin,                           //the lower bracket -- will be returned
	double &xmax                            //the upper bracket -- will be returned
){
	// test two points near first guess
	double DX = 0.01;
	double x=x0, x1=(1.-DX)*x0, x2=(1.+DX)*x0;
	double y=func(x), y1=func(x1), y2=func(x2);
	double ymin, ymax;

	// happy path 1
	// y is the zero
	if(y==0.0){
		xmin = x0;
		xmax = x0;
		return 0;
	}

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

	// all sampled points on same side

	// in ideal case, 
	// either y1 < y < y2 (UP)
	// or     y1 > y > y2 (DOWN)
	// in non-ideal cases 
	// could be y1 = y = y (FLAT)
	// could be y>y1 and y>y2, (CONCAVE)
	// or       y<y1 and y<y2
	enum class Slope {NONE, UP, DOWN, FLAT, CONCAVE_UP, CONCAVE_DOWN};
	Slope behavior = Slope::NONE;

	if( (y1==y) && (y2==y) ) behavior = Slope::FLAT;
	else if( (y1<=y && y<y2) || (y1<y && y<=y2) ) behavior = Slope::UP;
	else if( (y1>=y && y>y2) || (y1>y && y>=y2) ) behavior = Slope::DOWN;
	else if( (y1< y && y2<y) ) behavior = Slope::CONCAVE_DOWN;
	else if( (y1>y && y2> y) ) behavior = Slope::CONCAVE_UP;
	
	// above zero
	if(y>0){
		switch(behavior){
		case Slope::UP: {
			while(y1 > 0.0){
				DX *= 2.0;
				x1 = (1.-DX)*x0;
				y1 = func(x1);
			}
			xmin = x1;
			xmax = x0;
			return 0;
		} break;
		case Slope::CONCAVE_DOWN:
		case Slope::DOWN: {
			while(y2 > 0.0){
				DX *= 2.0;
				x2 = (1.+DX)*x0;
				y2 = func(x2);
			}
			xmax = x2;
			xmin = x0;
			return 0;
		} break;
		default:{
			xmin = nan("");
			xmax = nan("");
			perror("Unable to find a root near this guess");
			return 1;
		}
		}
	}
	// if below zero
	else if (y<0.0){
		switch(behavior){
		case Slope::UP: 
		case Slope::CONCAVE_UP: {
			while(y2 < 0.0){
				DX *= 2.0;
				x2 = (1.+DX)*x0;
				y2 = func(x2);
			}
			xmax = x2;
			xmin = x0;
			return 0;
		} break;
		case Slope::DOWN: {
			while(y1 < 0.0){
				DX *= 2.0;
				x1 = (1.-DX)*x0;
				y1 = func(x1);
			}
			xmin = x1;
			xmax = x0;
			return 0;
		} break;
		default:{
			xmin = nan("");
			xmax = nan("");
			perror("Unable to find a root near this guess");
			return 1;
		}
		}
	}

	//swap if they are backwards
	// this conditionn should be unreachable
	if(xmin > xmax){
		double temp = xmax;
		xmax = xmin; xmin = temp;
	}
	return 0;
}

// this will find brackets using newton's method to look for the next zero, then a bit beyond it
// why use one zero-finding method to prepare another zero-finding method? 
//    because Newton's method can fail, but bisection searches can go to nearly arbitrary accuracy
int rootfind::bisection_find_brackets_newton(
	std::function<double(double)>& func, //the function to find zero of
	double const x0,                     //an initial guess for the zero
	double& xmin,                        //the lower bracket -- will be returned
	double& xmax                         //the upper bracket -- will be returned
){
	double y1=func(x0), y2=y1;	
	double x1=x0, x2=x0;
	double xdx, ydx;
	double ymax, ymin;
	
	// happy path: x0 is already a zero
	if (y1 == 0.0 ) {
		xmin = 0.99*x0;
		xmax = 1.01*x0;
		y1 = func(xmin);
		y2 = func(xmax);
	}

	//while the two ys are on same side of axis, keep reposition until zero is bound
	//we use Newton's method to the nearest zero
	while(y1*y2 >= 0.0){
		//compute numerical derivative
		xdx = 1.01*x2;
		ydx = func(xdx);
		xdx = x2 - y2*(xdx-x2)/(ydx-y2);
		//Newton's method can fail at this if it only approaches the zero from one direction
		//In such a case, slightly broaden the bracket to get to other side of zero
		if( xdx == x2 ) {
			if(xdx>x1) xdx *= 1.01;
			else if(xdx<x1) xdx *= 0.99;
			else xdx = pseudo_unif(0.0, 2.0) * x1;
		}
		//limit amount of change allowed in single step
		if(xdx > 1.1*x2) xdx = 1.1*x2;	//if we increased, don't increase too much
		if(xdx < 0.9*x2) xdx = 0.9*x2;	//if we decreased, don't decrease too much
		ydx = y2;
		y2 = func(xdx);
		x2 = xdx;
	}
	//sometimes the brackets will be on either end of hyperbolic divergence
	//check that solution changes continuously between the two brackets
	double x3 = 0.5*(x1+x2), y3 = func(x3);
	double scale = 1.01;
	while( (y3*y1>=0.0 & fabs(y3)>fabs(y1)) | (y3*y2>=0.0 & fabs(y3)>fabs(y2))){
		//this can sometimes be solved by taking broader steps in secant method
		scale *= 1.1;
		x2=x1;
		y2=y1;
		while(y1*y2 >= 0.0){
			//compute numerical derivative
			xdx = scale*x2;
			ydx = func(xdx);
			xdx = x2 - y2*(xdx-x2)/(ydx-y2);
			//Newton's method can fail at this if it only approaches the zero from one direction
			//In such a case, slightly broaden the bracket to get to other side of zero
			if( xdx == x2 ) {
				if(xdx>x1) xdx *= 1.01;
				if(xdx<x1) xdx *= 0.99;
			}
			ydx = y2;
			y2 = func(xdx);
			x2 = xdx;
		}
		x3 = 0.5*(x1+x2);
		y3 = func(x3);
	}
	
	if(x2 > x1){
		xmax = x2; ymax = y2;
		xmin = x1; ymin = y1; 
	}
	else {
		xmax = x1; ymax = y1;
		xmin = x2; ymin = y2;
	}
	//swap if they are backwards
	if(xmin > xmax){
		double temp = xmax;
		xmax = xmin; xmin = temp;
		temp = ymax;
		ymax = ymin; ymin = temp;
	}
	return 0;
}

//given brackets bounding a single zero, find the zero
// the brackets xmin, xmax MUST bound a single zero
double rootfind::bisection_search(
	std::function<double(double)>& func, //the function to find zero of
	double &x,                           //the location of zero -- will be returned
	double &xmin,                        //the lower bracket
	double &xmax                         //the upper bracket
){
	// if the brackets do not enclose a root, find better brackets
	double ymin = func(xmin), ymax=func(xmax);
	while (ymin * ymax > 0){
		bisection_find_brackets_newton(func, x, xmin, xmax);
		ymin = func(xmin);
		ymax = func(xmax);
	}

	x = 0.5*(xmin+xmax);
	double xold = fabs(xmax-xmin), xnew = xold;
	double y = 1.0, y2 = 0.0;

	// now begin bisection search
	while( fabs(y)>0.0 || std::isnan(y) ){	
		// continue trying new brackts until they move
		while ( (xnew == xold) || std::isnan(y) ) {
			y = func(x);
			if ( y*ymax > 0.0 ) {
				xmax = x;
				ymax = y;
			}
			else if ( y*ymin > 0.0 ) {
				xmin = x;
				ymin = y;
			}
			xnew = fabs(xmax - xmin);
			x = rootfind::pseudo_unif(xmin, xmax);
			// if it is truly stuck, y will not change either; then quit
			if ( y2==y ) break;
			y2 = y;
		}

		x = 0.5*(xmin+xmax);
		y = func(x);
			
		//if the brackets are not moving, stop the search
		if(y2==y) break;
		xold = xnew;
	}
	return y;
}


#endif