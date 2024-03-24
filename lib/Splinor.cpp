#ifndef SPLINORCLASS
#define SPLINORCLASS

#include "Splinor.h"
#include "matrix.h"
#include <stdio.h>
#include <cassert>

Splinor::Splinor(
	double const *const x, 
	double const *const y, 
	int const L,
	bc bc_type, // default natural
	double const yprime0, // default zero
	double const yprimeN  // default zero
) :
	len(L), xlast(0), xa(x[0]), xb(x[L-1]), bc_type(bc_type), yprimea(yprime0), yprimeb(yprimeN)
{	
	//the coefficient and interpolation methods assume that x is strictly increasing
	//if x is decreasing, then read it in backwards.
	xarr = new double[len];
	yarr = new double[len];
	if(xa<xb){
		for(int i=0; i<len; i++){
			xarr[i] = x[i];
			yarr[i] = y[i];
		}
	}
	else if (xa>xb){
		double xt = xb;
		xb = xa;
		xa = xt;
		for(int i=0; i<len; i++){
			xarr[i] = x[len-1-i];
			yarr[i] = y[len-1-i];
		}	
	}
	// solve for the coefficients
	S = new double[len];
	calculateSplineCoefficients(x, y, bc_type);
}

Splinor::~Splinor(){
	delete[] S;
	delete[] xarr;
	delete[] yarr;
}

std::size_t Splinor::findPosition(double const xpos){
	std::size_t xind = xlast;	//begin at position we last checked
	// two happy paths -- position is at beginning or end
	if      (xpos >= xb) xind = len-2;	//if at the end, give the end
	else if (xpos <= xa) xind = 0;		//if at beginning, give beginning
	// two happy paths -- index is the last one, or one past the last one
	else if (xarr[xlast  ] <= xpos && xpos < xarr[xlast+1]) xind = xlast;
	else if (xarr[xlast+1] <= xpos && xpos < xarr[xlast+2]) xind = xlast+1; 
	// otherwise perform a search for it
	else {
		std::size_t imin = 0, imax = len-2;
		// need to find xind such that
		// xarr[xind] <= xpos < xarr[xind+1]
		while (xpos < xarr[xind] || xpos >= xarr[xind+1]){
			// the position is greater than this range
			if (xarr[xind+1] <= xpos) {
				imin = xind+1;
			}
			// the position is smaller than this range
			else if (xpos < xarr[xind]) {
				imax = xind;
			}
			xind = std::size_t((imax+imin)/2);
		}
		assert(xarr[xind] <= xpos);
		assert(xpos < xarr[xind+1]);
	}
	xlast = xind; //remember position for next time we look - normally fit in order
	return xind;
}

std::tuple<double,double,double,double> Splinor::getCoefficients(std::size_t const xind){
	double a, b, c, d, h1;
	h1 = xarr[xind+1]-xarr[xind  ];
	a = (S[xind+1]-S[xind])/(6.*h1);
	b = S[xind]/2;
	c = (yarr[xind+1]-yarr[xind])/h1 - h1*(2.*S[xind]+S[xind+1])/6.;
	d = yarr[xind];
	return std::make_tuple(a, b, c, d);
}

double Splinor::operator()(double const xpos){
	return interp(xpos);
}

double Splinor::interp(double const xpos){
	std::size_t xind = findPosition(xpos);
	double a, b, c, d;
	std::tie(a, b, c, d) = getCoefficients(xind);
	double t = xpos - xarr[xind];
	return t * ( t * (a * t + b) + c ) + d;
}

double Splinor::deriv(double const xpos){
	std::size_t xind = findPosition(xpos);
	double a, b, c, d;
	std::tie(a, b, c, d) = getCoefficients(xind);
	double t = xpos - xarr[xind];
	return t * ( 3. * a * t + 2. * b ) + c;//formula for derivative of splinor
}

void Splinor::makeNaturalSpline(
	double const *const x,
	double const *const y
){
	// natural splines form an (N-2)X(N-2) system
	// with y'' = 0 at both ends
	S[0] = S[len-1] = 0.0;
	// the tridiagonal matrix has same form at each row
	std::size_t const N = len-2;
	// the matrix coefficients
	double *left = new double[N];
	double *diag = new double[N];
	double *rite = new double[N];
	// now fill in the rest
	for(std::size_t i=1; i<=N; i++){
		left[i-1] = x[i] - x[i-1];
		diag[i-1] = 2. * (x[i+1] - x[i-1]);
		rite[i-1] = x[i+1] - x[i];
	}
	// there is no left-diagonal on first row, nor right-diagonal on last row
	left[0] = rite[N-1] = 0.0;
	// create the RHS of the equation
	double *RHS = new double[N];
	for(std::size_t i=0; i<N; i++){
		RHS[i] = 6.0 * ( (y[i+2]-y[i+1])/(x[i+2]-x[i+1]) - (y[i+1]-y[i])/(x[i+1]-x[i]) );
	}
	// solve  the equation -- note this finds S_1, ... S_{len-2}
	matrix::invertTridiagonal(left, diag, rite, RHS, S+1, N);
	// clear the matrix coefficients
	delete[] left;
	delete[] rite;
	delete[] diag;
	delete[] RHS;
}

void Splinor::makeClampedSpline(
	double const *const x,
	double const *const y
){
	// the matrix coefficients
	double *left = new double[len];
	double *diag = new double[len];
	double *rite = new double[len];
	left[0] = rite[len-1] = 0.0;
	// setup the first row {b_0 c_0 0 ... 0}
	diag[0] = (x[1] - x[0]) * 2.0;
	rite[0] = (x[1] - x[0]);
	// setup the last row {0 ... 0 a_{N-1} b_{N-1}}
	diag[len-1] = (x[len-1] - x[len-2]) * 2.0;
	left[len-1] = (x[len-1] - x[len-2]);
	// now fill in the rest
	for(std::size_t i=1; i<len-1; i++){
		left[i] = x[i] - x[i-1];
		diag[i] = 2. * (x[i+1] - x[i-1]);
		rite[i] = x[i+1] - x[i];
	}
	// create the RHS of the equation
	double *RHS = new double[len];
	RHS[0] = 6.0 * ( (y[1]-y[0])/(x[1]-x[0]) - yprimea);
	for(std::size_t i=1; i<len-1; i++){
		RHS[i] = 6.0 * ( (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]) );
	}
	RHS[len-1] = 6.0 * (yprimeb - (y[len-1]-y[len-2])/(x[len-1]-x[len-2]));
	// solve the equation
	matrix::invertTridiagonal(left, diag, rite, RHS, S, len);
	// clear the matrix coefficients
	delete[] left;
	delete[] rite;
	delete[] diag;
	delete[] RHS;
}

void Splinor::makeQuadraticSpline(
	double const *const x,
	double const *const y
){
	// the matrix coefficients
	double *left = new double[len];
	double *diag = new double[len];
	double *rite = new double[len];
	left[0] = rite[len-1] = 0.0;
	// setup the first row {b_0 c_0 0 ... 0}
	diag[0] = 1.0;
	rite[0] = -1.0;
	// setup the last row {0 ... 0 a_{N-1} b_{N-1}}
	diag[len-1] = 1.0;
	left[len-1] = -1.0;
	// now fill in the rest
	for(std::size_t i=1; i<len-1; i++){
		left[i] = x[i] - x[i-1];
		diag[i] = 2. * (x[i+1] - x[i-1]);
		rite[i] = x[i+1] - x[i];
	}
	// create the RHS of the equation
	double *RHS = new double[len];
	RHS[0] = 0.0;
	for(std::size_t i=1; i<len-1; i++){
		RHS[i] = 6.0 * ( (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]) );
	}
	RHS[len-1] = 0.0;
	// solve  the equation
	matrix::invertTridiagonal(left, diag, rite, RHS, S, len);
	// clear the matrix coefficients
	delete[] left;
	delete[] rite;
	delete[] diag;
	delete[] RHS;
}

void Splinor::calculateSplineCoefficients(
	double const *const x,
	double const *const y,
	bc const bc_type
){
	switch(bc_type){
		case bc::NATURAL:
			makeNaturalSpline(x,y);
			break;
		case bc::CLAMPED:
			makeClampedSpline(x,y);
			break;
		case bc::QUADRATIC:
			makeQuadraticSpline(x,y);
			break;
		default:
			// if for some reason this is unset, make a natural spline
			makeNaturalSpline(x,y);
	}
}

#endif