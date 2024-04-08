//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
//  ChandrasekharWD++.h
//  		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//	 		based on a cold degenerate electron gas equation of state
//	 		Chandrasekhar defines a variable y in equation (12)
//	 	This model assumes T=0, and ignores Coulombic and other effects
//	 	The surface is not treated in any special way
//	 	Updated to include composition gradient, indicated by mu_e
//  Reece Boston, Mar 24, 2022
// **************************************************************************************/

#ifndef ChandrasekharWDCLASS
#define ChandrasekharWDCLASS

#include "ChandrasekharWD++.h"

void ChandrasekharWD::chemical_gradient(const double xi, const double dx, double& mu, double& dmu){
	if(acore==1.){
		mu = mu0;
		dmu = 0.0;
		return;
	}
	double xcore = dx*(acore*len);
	double xswap = dx*(aswap*len);
	double EXP = exp(k*(xi-xswap));
	if(xi<xcore){
		mu = mu0;
		dmu= 0.0;
	}
	else {
		mu = (mu0-1.) + 1.0/(1.+EXP);
		dmu= -k*EXP*pow(1.+EXP,-2);
	}
}

void ChandrasekharWD::basic_setup(){
 	//exclude unphysical values of Y0
 	if(Y0 < 1.) Y0 = 1.;
 	Y02 = Y0*Y0;
	X0 = sqrt(Y02-1.);
 	X02 = Y02-1.;
 	//name this model for files
 	name = strmakef("ChandrasekharWD.%1.3lf", Y0);

 	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;

 	// prepare stellar arrays
	xi= new double[len]; //normalized radius
	x = new double[len]; //the relativity factor x = pF/mc
	y = new double[len]; //Chandrasekhar's y, y^2=1+x^2
	z = new double[len]; //derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
 }

 void ChandrasekharWD::init_arrays(){
 	mass = new double[len];
 	f    = new double[len];
 	mue  = new double[len];
 	dmue = new double[len];
	mass[0] = 0.0; mue[0] = mu0; dmue[0] = 0.0;
	x[0] = X0;
	f[0] = Chandrasekhar::factor_f(x[0]);
 	for(std::size_t X=1; X<len; X++){
		x[X] = sqrt(y[X]*y[X]-1.);
		if(y[X]<1.0) x[X] = 0.0;
		f[X] = Chandrasekhar::factor_f(x[X]);
 		chemical_gradient(xi[X], dx, mue[X], dmue[X]);
 		mass[X] = -4.*m_pi*xi[X]*xi[X]*z[X]/mue[X];
 	}
 }

//initalize white dwarf from central value of y and length
ChandrasekharWD::ChandrasekharWD(
	double Y0, std::size_t L, 
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(L), mu0(mu0), k(k), acore(acore), aswap(aswap),
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0)
{

	basic_setup();

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star
	dx = sqrt(6.0)/(len-1);
	double yS=1.0, ddx = dx;
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	//this is used to control bisection search to guarantee it terminates
	double dxold = 1.0;
		
	yS = RK4integrate(len, dx);
	
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !std::isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(std::isnan(yS)){
				dx = 1.0/len; ddx *= 0.1; yS = 1.0;
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx *= 0.5;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	if(dx<0) {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}
	
	//now use bisection to find dx so that yS=0.0
	while( fabs(yS)>0.0 || std::isnan(yS) ){
		dx = 0.5*(dxmin+dxmax);
		yS = RK4integrate(len, dx);
		
		if(std::isnan(yS)){
			yS = -1.0;
		}
		if( (yS*ySmax>0.0) ){
			dxmax = dx;
			ySmax = yS;
		}
		else if( (yS*ySmin>0.0) ){
			dxmin = dx;
			ySmin = yS;
		}
		//if the brackets are not moving, stop the search
		if(dxold == fabs(dxmin-dxmax)) break;
		dxold = fabs(dxmin-dxmax);
	}
	RK4integrate(len, dx, 1);

	//now set physical properties of the white dwarf
	init_arrays();

	// find the fitting point for modes
	indexFit = len/2;
	for(std::size_t X=1; X<len; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(x[X-1]>0.5 & x[X+1]<=0.5) indexFit = X;
	}
	indexFit /= 2;

	// prepare boundary expansions
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}

//initalize white dwarf from central value of y and length, with specific step size
//this should only be used for testing scaling issues
ChandrasekharWD::ChandrasekharWD(
	double Y0, std::size_t L, const double dx, 
	double mu0, double k, double acore, double aswap
)
	: Y0(Y0), len(L), dx(dx), mu0(mu0), k(k), acore(acore), aswap(aswap),
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0)
{	
	basic_setup();
		
	//we know dx, so no need for bisection search
	RK4integrate(len, dx, 1);

	//now set physical properties of the white dwarf
	
	//set initial density, pressure
	init_arrays();
	
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %lu\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}

//initialize a white dwarf with specific values of A0, B0, constant chemical profile
//this can allow for matching to tabulated values of white dwarfs
ChandrasekharWD::ChandrasekharWD(
	double Y0, std::size_t L, 
	const double A0, const double B0
)
	: Y0(Y0), len(L), A0(A0), B0(B0), 
	acore(1.0),  aswap(0.0), mu0(1.0), k(0.0)
{
	basic_setup();

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	//an initial guess -- based on constant volume star
	dx = sqrt(6.0)/(len-1);
	double yS=1.0, ddx = dx;
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	//this is used to control bisection search to guarantee it terminates
	double dxold = 1.0;
		
	yS = RK4integrate(len, dx);
	
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !std::isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(std::isnan(yS)){
				dx = 1.0/len; ddx *= 0.1; yS = 1.0;
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx *= 0.5;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	if(dx<0) {printf("somehow dx is negative...\n"); dx=-dx;}
	if(ySmin*ySmax > 0.0) {printf("big problem, chief\n"); exit(EXIT_FAILURE);}
	
	//now use bisection to find dx so that yS=0.0
	while( fabs(yS)>0.0 || std::isnan(yS) ){
		dx = 0.5*(dxmin+dxmax);
		yS = RK4integrate(len, dx);
		
		if(std::isnan(yS)){
			yS = -1.0;
		}
		if( (yS*ySmax>0.0) ){
			dxmax = dx;
			ySmax = yS;
		}
		else if( (yS*ySmin>0.0) ){
			dxmin = dx;
			ySmin = yS;
		}
		//if the brackets are not moving, stop the search
		if(dxold == fabs(dxmin-dxmax)) break;
		dxold = fabs(dxmin-dxmax);
	}
	RK4integrate(len, dx, 1);

	//now set physical properties of the white dwarf
	init_arrays();

	indexFit = len/2;
	for(std::size_t X=1; X<len; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(x[X-1]>0.5 & x[X+1]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}

ChandrasekharWD::~ChandrasekharWD(){
	delete[] xi;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] mass;
	delete[] f;
	delete[] mue;
	delete[] dmue;
}

//integrate the polytrope up to Len using RK4
double ChandrasekharWD::RK4integrate(const std::size_t Len, double dx){
	//the YC, ZC, XC are corrected values of y, z, x to be used in equations
	double YC,ZC,XC,   YC12, YHe4, YH1;
	double MUC, DMUC;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	double K[4], L[4], K_He4[4], K_H1[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	//near the surface, y can be <1, then (y^1-1)^3/2 can lead to NaNs
	//to avoid this, make YCN complex and only take the real part
	std::complex<double> YCN(0.0,0.0);//helps avoid NaNs near edge
	double MC14 = 0.0, MHe4 = 0.0, MH1 = 0.0;
	
	//set our initial conditions
	xi[0] = 0.0; y[0] = Y0; z[0] = 0.0;
	for(std::size_t X = 0; X<Len-1; X++){
		XC = xi[X]; YC = y[X]; ZC = z[X];
		chemical_gradient(XC,dx, MUC, DMUC);
		for(std::size_t a = 0; a<4; a++){
			//find current calues
			YCN = YC*YC - 1.;
			YCN = pow(YCN, 1.5);
			// K = dy = dxi*(dy/dxi) = dxi*z
			K[a] =  dx*ZC;
			// L = dz = dxi*(dy^2/dxi^2) = dxi*[ -(y^2-1)^3/2 - 2(dy/dxi)/xi ]
			L[a] = -dx*( YCN.real()*MUC*MUC + 2.0*ZC/XC -  ZC*DMUC/MUC);
			if(XC==0) L[0] = -dx/3.0*YCN.real()*MUC*MUC; //must use analytic expression at very center
			//calculate "corrected" positions using previous shift vectors
			XC = (double(X)+B[a])*dx;//x[X] + B[a]*dx;
			YC = y[X] + B[a]*K[a];
			ZC = z[X] + B[a]*L[a];
			chemical_gradient(XC,dx, MUC, DMUC);
		}		
		xi[X+1] = double(X+1)*dx;
		y[X+1] = y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		z[X+1] = z[X] + L[0]/6.0 + L[1]/3.0 + L[2]/3.0 + L[3]/6.0;
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(y[X+1]<1.0) {return y[X+1]-1.;}
	}
	return y[Len-1]-1.;
}

std::size_t ChandrasekharWD::RK4integrate(const std::size_t Len, double dx, int grid){
	grid=1;
	//the YC, ZC, XC are mid-point values of y, z, x to be used in equations
	double YC,ZC,XC;
	double MUC, DMUC;
	//the arrays K and L contain the corrections to y and z
	double K[4], L[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	std::complex<double> YCN(0.0,0.0);// = (0.0i);
	
	//set our initial conditions
	x[0] = 0.0; y[0] = Y0; z[0] = 0.0;
	for(std::size_t X = 0; X<Len-1; X++){
		XC = double(X)*dx; YC = y[X]; ZC = z[X];
		chemical_gradient(XC,dx, MUC, DMUC);
		for(std::size_t a = 0; a<4; a++){
			//to handle NaNs at surface, use real part of complex root
			YCN = YC*YC - 1.;
			YCN = pow(YCN, 1.5);
			K[a] = dx*ZC;
			L[a] = -dx*( YCN.real()*MUC*MUC + 2.0*ZC/XC - ZC*DMUC/MUC);
			if(X==0) L[0] = -dx/3.0*YCN.real()*MUC*MUC;
			//calculate "corrected" positions using previous shift vectors
			XC = (double(X)+B[a])*dx;
			YC = y[X] + B[a]*K[a];
			ZC = z[X] + B[a]*L[a];
			chemical_gradient(XC,dx, MUC, DMUC);
		}
		xi[X+1] = double(X+1)*dx;
		y[X+1] = y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		z[X+1] = z[X] + L[0]/6.0 + L[1]/3.0 + L[2]/3.0 + L[3]/6.0;
	}
	return Len;
}

//Here we define functions to access radius, pressure, etc.
double ChandrasekharWD::rad(std::size_t X){
	return Rn*xi[X];
}
double ChandrasekharWD::rho(std::size_t X){
	using namespace Chandrasekhar;
	return B0*mue[X]*pow(x[X],3);
}
double ChandrasekharWD::drhodr(std::size_t X){
	using namespace Chandrasekhar;
	return 3.*B0*mue[X]*x[X]*y[X]*z[X]/Rn + B0*pow(x[X],3)*dmue[X]/Rn; //note: dx/dxi = z*y/x
}
double ChandrasekharWD::P(std::size_t X){
	using namespace Chandrasekhar;
	return A0*f[X];
}
double ChandrasekharWD::dPdr(std::size_t X){
	using namespace Chandrasekhar;
	return 8.*A0*pow(x[X],3)*z[X]/Rn; //see eqn (9) in Chandrasekhar 1935
}
double ChandrasekharWD::Phi(std::size_t X){
	using namespace Chandrasekhar;
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return -8.*A0/B0 * (y[X] - 1. + z[len-1]*xi[len-1])/mue[X];
}
double ChandrasekharWD::dPhidr(std::size_t X){
	using namespace Chandrasekhar;
	return -8.*A0/B0 * z[X]/Rn/mue[X];
}
double ChandrasekharWD::mr(std::size_t X){
	using namespace Chandrasekhar;
	return B0*pow(Rn,3)*mass[X];
}

double ChandrasekharWD::Schwarzschild_A(std::size_t X, double GamPert){
	if(GamPert==0.0) return dmue[X]/mue[X]/Rn;
	else        	 return dmue[X]/mue[X]/Rn-z[X]*(8.*pow(x[X],3)/f[X]/GamPert-3.*y[X]/pow(x[X],2))/Rn;
}

double ChandrasekharWD::getAstar(std::size_t X, double GamPert){
	if(GamPert==0.0) return -dmue[X]/mue[X]*xi[X];
	else        	 return -dmue[X]/mue[X]*xi[X] + z[X]*xi[X]*(8.*pow(x[X],3)/f[X]/GamPert-3.*y[X]/pow(x[X],2));
}

double ChandrasekharWD::Ledoux(std::size_t X, double GamPert){
	return dmue[X]/mue[X]/Rn;
}

double ChandrasekharWD::getU(std::size_t X){
	if(X==0) return 3.0;
	return -mue[X]*mue[X]*xi[X]*pow(x[X],3)/z[X];
}

double ChandrasekharWD::getVg(std::size_t X, double GamPert){
	if(GamPert==0.0) return -8.*pow(x[X],3)*xi[X]*z[X]/f[X]/Gamma1(X);
	else			 return -8.*pow(x[X],3)*xi[X]*z[X]/f[X]/GamPert;
}

double ChandrasekharWD::getC(std::size_t X){
	if(X==0) return z[len-1]*xi[len-1]/2./yc[1]*mue[0]/mue[len-1];
	return (z[len-1]/z[X])*(xi[X]/xi[len-1])*(mue[X]/mue[len-1]);
}

double ChandrasekharWD::Gamma1(std::size_t X){
	return 8.*pow(x[X],5)/(3.*y[X]*f[X]);
}

double ChandrasekharWD::sound_speed2(std::size_t X, double GamPert){
	using namespace Chandrasekhar;
	if(GamPert == 0.0) return Gamma1(X)*A0/B0*f[X]/pow(x[X],3)/mue[X];
	else               return GamPert  *A0/B0*f[X]/pow(x[X],3)/mue[X];
}

double ChandrasekharWD::Radius(){return Rn*xi[len-1];}	//total radius
double ChandrasekharWD::Mass(){return mr(len-1);}//total mass
double ChandrasekharWD::Gee(){return G_CGS;}


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void ChandrasekharWD::setupCenter(){	
	//the central coefficients -- expanded in terms of r/R
	double z1 =  z[len-1];
	double x1 = xi[len-1];
	double mu2 = mue[0]*mue[0];
	//
	yc[0]=	Y0;
	yc[1]=	-X02*X0/6.*pow(x1,2)*mu2;
	yc[2]=	X02*X02*Y0/40.*pow(x1,4)*mu2*mu2;
	yc[3]=	X0*(-14.*X02*X02-19.*X02*X02*X02)/5040.*pow(x1,6)*pow(mu2,3);
	//
	xc[0]=	X0;
	xc[1]=	-X02*Y0/6.*pow(x1,2)*mu2;
	xc[2]=  X02*X0*(4.+9.*X02)/360.*pow(x1,4);
	//
	fc[0]=	Chandrasekhar::factor_f(X0);
	fc[1]=	xc[1]*X0*(2. - 3.*X0 + 3.*X02)/3.*pow(x1,2);
}

void ChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double xi2 = xi[len-1]*xi[len-1];
 	double ac1 = (8.*X0*X02/Gam1/fc[0] - 3.*Y0/X02);
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) Ac[0] = 0.0;
	if(maxPow>=2) Ac[1] = 2.*yc[1]*ac1;
	if(maxPow>=4) Ac[2] = 4.*yc[2]*ac1
 			- yc[1]*16.*X02*X0*fc[1]/Gam1/fc[0]/fc[0]
 			+ yc[1]*48.*X02*xc[1]/Gam1/fc[0]
 			+ yc[1]*12.*xc[1]*Y0/pow(X0,3)
 			- 6.*yc[1]*yc[1]/X02;	
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
 	double x1 = xi[len-1];	
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -16.*X02*X0*yc[1]/Gam1/fc[0];
 	if(maxPow>=4) Vc[2] =-8.*   pow(x1,2)*pow(X0,6)*(5.*fc[1]+4.*pow(x1,2)*fc[0]*X0*Y0)/(15.*Gam1*fc[0]*fc[0]);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void ChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = xi[len-1];
	double mu2 = mue[0]*mue[0];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = -0.6*mu2*x1*x1*X0*Y0;
	if(maxPow>=4) Uc[2] = 9.*xc[2]/X0 - pow(x1*x1*mu2,2)*X02*(0.08+0.134*X02);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getC1Center(double *cc, int& maxPow){
	double x1 = xi[len-1];
	double z1 = z[len-1];
	double murat = mue[0]/mue[len-1];
	if(maxPow>=0) cc[0] = z1*x1/2./yc[1]*murat;
	if(maxPow>=2) cc[1] = -z1*x1*yc[2]/yc[1]/yc[1]*murat;
	if(maxPow>=4) cc[2] = -z1*x1*(1.5*yc[1]*yc[3]-2.*yc[2]*yc[2])/pow(yc[1],3)*murat;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}


// **************************  SURFACE BOUNDARY  ***************************************
// the following provide coefficients for surface expansions of A*, Vg, U, c1 in terms of t=1-r/R
//	Note that A*, Vg require a power -1
//	Note that up to order N, we only require:
//		A*, Vg, c1 up to order N-1
//		U          up to order N
//	Here maxPow represents the maixmum order of expansions of y1,..,y4 in LAWE
//	If maxPow = 0, we need terms -1
//	If maxPow = 1, we need terms -1, 0
//	If maxPow = 2, we need terms -1, 0, 1
void ChandrasekharWD::setupSurface(){}

void ChandrasekharWD::getAstarSurface(double *As, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 = xi[len-1];
	double a1 =-x1*z[len-1];
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) As[O-1] = 0.0;
	if(maxPow>= 0) As[O  ] = 0.75*a1;
	if(maxPow>= 1) As[O+1] = 0.75*a1 - (8./Gam1 + 0.375)*a1*a1;
	if(maxPow>= 2) As[O+2] = 0.0;
	if(maxPow>= 3) As[O+3] = 0.0;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void ChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(len-2) : g);
	double x1 = xi[len-1];
	double a1 =-x1*z[len-1];
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) Vs[O-1] = 1.5;
	if(maxPow>= 0) Vs[O  ] = 0.75*a1;
	if(maxPow>= 1) Vs[O+1] = 0.75*a1 + (8./Gam1 - 0.375)*a1*a1;
	if(maxPow>= 2) Vs[O+2] = a1*3./16.*pow(a1-2.,2)+4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) Vs[O+3] = -3.*a1*pow(a1-2,3)/32.+(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void ChandrasekharWD::getUSurface(double *Us, int& maxPow){
	//these happen to all be zero at surface
	for(int j=0; j<=maxPow; j++){
		Us[j] = 0.0;
	}
}

void ChandrasekharWD::getC1Surface(double *cs, int& maxPow){
	if(maxPow >= 0) cs[0] =  1.;
	if(maxPow >= 1) cs[1] = -3.;
	if(maxPow >= 2) cs[2] =  3.;
	if(maxPow >= 3) cs[3] = -1.;
	if(maxPow >= 4) cs[4] =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow  > 4) maxPow = 4;
}

void ChandrasekharWD::writeStar(char const *const c){
	//create names for files to be opened
	std::string outputdir;
 	if(c==NULL)	outputdir = strmakef("./out/%c", name.c_str());
 	else        outputdir = strmakef("./%s/star",c);

 	system( ("mkdir -p "+outputdir).c_str() );

 	printStar(outputdir.c_str());
 	printDeg(outputdir.c_str());
 	printBV(outputdir.c_str());
 	printChem(outputdir.c_str());
 	printCoefficients(outputdir.c_str());
}

void ChandrasekharWD::printDeg(char const *const outputdir){
 	std::string txtname, outname;
 	txtname = addstring(outputdir, "/degeneracy.txt");
 	outname = addstring(outputdir, "/degeneracy.png");

 	FILE *fp = fopen(txtname.c_str(), "w");
 	for(std::size_t X=0; X < length(); X++){
 		fprintf(fp, "%0.16le\t", xi[X]);
 		fprintf(fp, "%le\t%le\t%le", x[X], y[X], f[X]);
 		fprintf(fp, "\n");
 	}
 	fclose(fp);	

 	FILE *gnuplot = popen("gnuplot -persist", "w");
 	fprintf(gnuplot, "reset\n");
 	fprintf(gnuplot, "set term png size 1600,800\n");
 	fprintf(gnuplot, "set samples %lu\n", length());
 	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
 	fprintf(gnuplot, "set title 'Degenerate functions for %s'\n", graph_title().c_str());
 	fprintf(gnuplot, "set xlabel 'xi'\n");
 	fprintf(gnuplot, "set ylabel 'x, y, f'\n");
  	fprintf(gnuplot, "set xrange [-0.01:1.01]\n");
  	fprintf(gnuplot, "set logscale y 10\n");
  	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'x'", txtname.c_str());
  	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'y'", txtname.c_str());
  	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'f'", txtname.c_str());
 	fprintf(gnuplot, "\n");
 	fprintf(gnuplot, "exit\n");
 	pclose(gnuplot);
 }

void ChandrasekharWD::printChem(char const *const outputdir){

 	std::string txtname, outname;
 	txtname = addstring(outputdir, "/chemical.txt");
 	outname = addstring(outputdir, "/chemical.png");

 	//print the chemical gradient
 	FILE *fp  = fopen(txtname.c_str(), "w");
 	double H,He;
 	for(std::size_t X=1; X< length()-1; X++){
 		He = 2.*(mue[X]-1.)/mue[X];
 		H = 1.-He;
 		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
 			xi[X],
 			(1.-mr(X)/Mass()),
 			H,
 			He
 		);
 	}
 	fclose(fp);	
 	//plot file in png in gnuplot
 	FILE *gnuplot = popen("gnuplot -persist", "w");
 	fprintf(gnuplot, "reset\n");
 	fprintf(gnuplot, "set term png size 1000,800\n");
 	fprintf(gnuplot, "set samples %lu\n", length());
 	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
 	fprintf(gnuplot, "set title 'Chemical composition for %s'\n", graph_title().c_str());
 	fprintf(gnuplot, "set logscale x 10\n");
 	fprintf(gnuplot, "set format x '%%L'\n");
 	fprintf(gnuplot, "set xlabel 'log_{10} (1-m/M)\n");
 	fprintf(gnuplot, "set ylabel 'abundance\n");
 	fprintf(gnuplot, "set yrange [-0.01: 1.01]\n");
 	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
 	fprintf(gnuplot, "plot '%s' u 2:3 w l t 'X (H1)'",  txtname.c_str());
 	fprintf(gnuplot, ",    '%s' u 2:4 w l t 'Y (He4)'", txtname.c_str());
 	fprintf(gnuplot, "\n");
 	fprintf(gnuplot, "exit\n");
 	pclose(gnuplot);	
 }

#endif