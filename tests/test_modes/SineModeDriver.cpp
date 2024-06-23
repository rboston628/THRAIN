//**************************************************************************************
//						SINE MODE DRIVER
// SineModeDriver.cpp
//     Solves simple sinusoidal equation for sine waves.
//     The equation solved is
//       x dy1/dx = + x * w * y2
//       x dy2/dx = - x * w * y1
//     which should exactly yield y1 = sin(w*x), y2 = cos(w*x)
//     note: for consistency with other modes, y1 = x^(l-2) ytilde1
//	THIS IS NOT PART OF THRAIN -- ONLY FOR TESTS OF MODE.cpp AS SOLVER
//  Reece Boston, Mar 14, 2024
//***************************************************************************************

#ifndef SINEMODEDRIVERCLASS
#define SINEMODEDRIVERCLASS

#include "SineModeDriver.h"

//constructor
SineModeDriver::SineModeDriver(Star *st) : ModeDriver(num_var, st)
{	
	// set lengths
	len_star = star->length();
	len = (len_star%2==1? (len_star+1)/2 : len_star/2+1);
}

//destructor
SineModeDriver::~SineModeDriver(){}

void SineModeDriver::varnames(std::string* names) {
    names[0] = "sine"; names[1] = "cosine";
}

std::size_t SineModeDriver::length(){return len;}

double SineModeDriver::Gamma1(){return 0.0;}

double SineModeDriver::rad(std::size_t X){return double(X)/double(len-1);}

//these functions come from coupled wave equations -- see Unno et al. Ch 18
void SineModeDriver::getCoeff(
	double *CCI,
	std::size_t X,
	int b,
	double sig2,
	int L
)
{	
	double (*CC)[num_var] =(double (*)[num_var]) CCI;
	double a11 = double(2-L) * x(X,b);
	double a12 = sqrt(sig2) * x(X,b);
	CC[sin][sin] =  a11;
	CC[sin][cos] =  a12;
	CC[cos][sin] = -a12;
	CC[cos][cos] =  a11;
}

double SineModeDriver::x(std::size_t const X, int const b){
	return double(2*X + b)/double(len_star-1);
}

void SineModeDriver::getBoundaryMatrix(int nv, double **y, int* indexOrder){
	const double yy[num_var][num_var] = {
		{ 1.0, 1.0 },
		{ 1.0, 1.0 }
	};
	const int indices[num_var] = {0,0};
	for(int i=0; i<num_var; i++)
		for(int j=0; j<num_var; j++)
			y[i][j] = yy[i][j];
	for(int i=0; i<num_var; i++) indexOrder[i] = indices[i];
}

void SineModeDriver::setupBoundaries(){}

std::size_t SineModeDriver::CentralBC(double **ymode, double *y0, double omeg2, int l, int m){
	ymode[0][sin] = 0.0;
	ymode[0][cos] = 1.0;
	double X = sqrt(omeg2) * rad(1);
	ymode[1][sin] = X - pow(X,3)/6.0 + pow(X,5)/120.0 - pow(X,7)/5040.0;
	ymode[1][cos] = 1.- pow(X,2)/2.0 + pow(X,4)/24.0  - pow(X,6)/720.0;
	for(std::size_t x=0; x<=1; x++){
		ymode[x][sin] *= pow(rad(x), 2-l);
		ymode[x][cos] *= pow(rad(x), 2-l);
	}
	return 1;
}

std::size_t SineModeDriver::SurfaceBC(double **ymode, double *ys, double omeg2, int l, int m){
	ymode[len-1][sin] = 0.0;
	ymode[len-1][cos] = 1.0;
	return len-1;
}

double SineModeDriver::SSR(double w2, int L, ModeBase* mode){
        double sum = 0.0;
        double y1 , y2;
        for(std::size_t X=1; X<len; X++){
            y1 = mode->getY(sin, X);
            y2 = mode->getY(cos, X);
            sum += std::abs(y1*y1 + y2*y2 - 1.0);
        }
        return (sum/len);
}

double SineModeDriver::tidal_overlap(ModeBase* mode){return 1.0;}

double SineModeDriver::innerproduct(ModeBase* mode1, ModeBase* mode2){
	double dx = 1.0/double(len-1);
	double sum1 = 1.0, sum2 = 0.0, integral = 0.0;
	double sin1, sin2, cos1, cos2;
	for(std::size_t i=1; i<len; i++){
		sin1 = mode1->getY(sin, i);
		sin2 = mode2->getY(sin, i);
		cos1 = mode1->getY(cos, i);
		cos2 = mode2->getY(cos, i);
		sum2 = sum1;
		sum1 = (sin1*sin2 + cos1*cos2);
		integral += 0.5*dx*(sum1+sum2);
	}
	return integral;
}

#endif