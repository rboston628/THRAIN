//**************************************************************************************
//							POLYTROPIC STAR
// Polytrope.cpp
//		A simple stellar object in Newtonian physics obeying everywhere a polytropic
//			equation of state of the form 
//				P~rho^Gamma
//		Solves Lane-Emden equation, y=theta, x=xi, z=dy/dxi in more usual notation
//		See Hansen & Kawaler Chapter 7 for further information
// Reece Boston, Mar 24, 2022
//**************************************************************************************

#ifndef POLYTROPECLASS
#define POLYTROPECLASS

#include "ThrainMain.h"
#include "Polytrope.h"

//initalize polytrope from index and length, fit to mass M, radius R
Polytrope::Polytrope(double BigM, double BigR, double n, std::size_t L)
	: n(n), len(L), Gamma(1.0+1.0/n)
{
	std::size_t Nedge = 0; //deprecated functionality for finer surface resolution
	//if index is out of range, fail
	if( n >= 5.0 || n < 0.0 ){
		n = nan(""); len = 0;
		printf("\nInvalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	//name this polytrope for files
	name = strmakef("polytrope.%1.1f", n);
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess
	double dx = sqrt(6.0)/(len-1), yS=1.0, ddx = dx;
	Y = new double*[len];
	for(std::size_t i=0;i<len;i++) 
		Y[i] = new double[numvar];
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	std::size_t stop=0;
	double dxold = 1.0;
		
	yS = RK4integrate(len, dx);
	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !std::isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(std::isnan(yS)){
				if(n!= int(n)) yS = -1.0;
				else {
					dx = 1.0/len; ddx *= 0.1; yS = 1.0;
				}
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx = 0.99*dx;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	//some possible errors that might occur
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
		stop++;
	}	
	RK4integrate(len, dx, 1);
		
	//now set physical properties of the polytrope
	
	//set initial density, pressure
	// the physical units can be included through homology
	GG = G_CGS;
	rho0 = BigM/(4.*m_pi)*pow(BigR,-3)*(-Y[len-1][x]/Y[len-1][z]);
	P0   = G_CGS*pow(BigM,2)*pow(BigR,-4)/(4.*m_pi)*pow(Y[len-1][z],-2)/(n+1.);
	Rn   = BigR/Y[len-1][x];
	
	mass = new double[len];
	base = new double[len];
	base[0] = 1.0;
	mass[0] = 0.0;
	indexFit = len/2;
	for(std::size_t X=1; X<len; X++){
		//as we scan through x,y,z, set matching point as where y[X] = 0.5
		if(Y[X-1][y]>0.5 & Y[X][y]<=0.5) indexFit = X;
		base[X] = pow(Y[X][y], n-1.0);
		if(Y[X][y]<0.0){
			std::complex<double> YN = Y[X][y];
			YN = pow(YN,n-1);
			base[X] = YN.real();
		}
		//there is a formula (see H&K pg 268)
		//multiply by Rn^3*rho0 when mr(X) is called
		mass[X] = -4.*m_pi*Y[X][x]*Y[X][x]*Y[X][z];
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
}

Polytrope::Polytrope(double n, std::size_t L)
	: n(n), len(L), Gamma(1.0+1.0/n)
{
	std::size_t Nedge = 0; //deprecated functionality for finer surface resolution
	if(n<5.0) Nedge = 0;
	//if index is out of range, fail
	if( n > 5.0 || n < 0.0 ){
		n = nan(""); len = 0;
		printf("\nInvalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}
	//name this polytrope for files
	name = strmakef("polytrope.%1.1f", n);
	
	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with y[len-1]=0.0
	//we will find the proper dx with a bisection search
	
	//an initial guess
	double dx = sqrt(6.0)/(len-1), yS=1.0, ddx = dx;
	Y = new double*[len];
	for(std::size_t i=0;i<len;i++) 
		Y[i] = new double[numvar];
	double dxmax=1.0, dxmin=0, ySmax=-1.0, ySmin=1.0;
	std::size_t stop=0;
	double dxold = 1.0;
	
	//for n=5, there is no edge, so skip the search for dx
	//this is really the easiest way to do this.
	//yes, I do feel dirty for it
	//the n=5 edge case should only come up in testing
	if(n==5.0) goto skipfindsurface;
	
	yS = RK4integrate(len, dx);

	//find brackets on dx that bound a zero in yS
	if(yS > 0){
		dxmin = dx; ySmin = yS;
		while(yS > 0 && !std::isnan(yS)){
			dx += ddx;
			yS = RK4integrate(len, dx);
			if(std::isnan(yS)){
				if(n!= int(n)) yS = -1.0;
				else {
					dx = 1.0/len; ddx *= 0.1; yS = 1.0;
				}
			}
		}
		dxmax = dx; ySmax = yS;
	}
	else if (yS < 0){
		dxmax = dx; ySmax = yS;
		while(yS < 0){
			dx = 0.99*dx;
			yS = RK4integrate(len,dx);
		}
		dxmin = dx; ySmin = yS;
	}
	dx = 0.5*(dxmin+dxmax);
	yS = RK4integrate(len, dx);
	
	//some possible errors that might occur
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
		stop++;
	}
	//we will skip to here if n=5
	skipfindsurface:
	RK4integrate(len, dx, 1);
		
	//now set physical properties of the polytrope
	
	//set initial density, pressure
	// specifying K, rho0 not necessary, just rescales the solutions
	// the physical units can be included through homology
	GG = 1.0;
	rho0 = 1.0;
	P0 = 1.0;
	Rn = sqrt( (n+1.)/(4.*m_pi) );
	
	mass = new double[len];
	base = new double[len];
	base[0] = 1.0;
	mass[0] = 0.0;
	indexFit = len/2;
	for(std::size_t X=1; X<len; X++){
		//as we scan through x,y,z, set matching point as where y[X] = 0.5
		if(Y[X-1][y]>0.5 & Y[X][y]<=0.5) indexFit = X;
		base[X] = pow(Y[X][y], n-1.0);
		if(Y[X][y]<0.0){
			std::complex<double> YN = Y[X][y];
			YN = pow(YN,n-1);
			base[X] = YN.real();
		}
		//there is a formula (see H&K pg 268)
		//multiply by Rn^3*rho0 when mr(X) is called
		mass[X] = -4.*m_pi*Y[X][x]*Y[X][x]*Y[X][z];
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
}


//initalize polytrope from index, length, and dx
//this method is intended only for use in testing scaling relations
Polytrope::Polytrope(double n, std::size_t L, const double dx)
	: n(n), len(L), Gamma(1.0+1.0/n)
{
	std::size_t Nedge = 0;
	if(n<5.0) Nedge = 20;
	//if index is out of range, fail
	if( n > 5.0 || n < 0.0 ){
		n = nan(""); len = 0;
		printf("\nInvalid polytropic index.  Polytrope not initialized.\n");
		exit(EXIT_FAILURE);
	}

	//name this polytrope for files
	name = strmakef("polytrope.%1.1f", n);
	//reserve room for arrays
	Y = new double*[len];
	for(std::size_t i=0;i<len;i++) 
		Y[i] = new double[numvar];
	
	//we know dx, so no need for bisection search
	RK4integrate(len, dx, 1);
	
	//now set physical properties of the polytrope
	mass = new double[len];
	base = new double[len];
	//set initial density, pressure
	// specifying K, rho0 not necessary, just rescales the solutions
	rho0 = 1.0;
	P0 = 1.0;
	Rn = sqrt( (n+1.)/(4.*m_pi) );
	
	//need to be replaced with a smoother averaging technique
	//for polytrope, take advantage of having z = dtheta/dx to avoid differencing
	base[0] = 1.0;
	mass[0] = 0.0;
	indexFit = std::size_t(len/2);
	for(std::size_t X=1; X<len; X++){
		base[X] = pow(Y[X][y], n-1.);
		if(Y[X][y]<0.0){
			std::complex<double> YN = Y[X][y];
			YN = pow(YN,n-1);
			base[X] = YN.real();
		}
		//there is a formula (see H&K pg 268)
		//multiply by Rn^3*rho0 when mr(X) is called
		mass[X] = -4.*m_pi*Y[X][x]*Y[X][x]*Y[X][z];
	}
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %d\n", indexFit);
	printf("r[indexFit] = %0.32le\t", rad(indexFit));
	printf("y[indexFit] = %0.32le\n", Y[indexFit][y]);
	setupCenter();
	setupSurface();
}

Polytrope::~Polytrope(){
	for(std::size_t i=0; i<len; i++)
		delete[] Y[i];
	delete[] mass;
	delete[] base;
}


void Polytrope::centerInit(double ycenter[numvar]){
	ycenter[x] = 0.0;
	ycenter[y] = 1.0;
	ycenter[z] = 0.0;
}

void Polytrope::RK4step(double dx, double yin[numvar], double yout[numvar]){
	double YC[numvar] = {yin[x], yin[y], yin[z]};
	double K[numvar][4];
	static const double B[4] = {0.5,0.5,1.0,0.0};
	std::complex<double> YCN(0.0,0.0);
	
	for(std::size_t a = 0; a<4; a++){
		//use complex y^n to avoid NaNs when y<0
		YCN = YC[y];
		YCN = pow(YCN,n);
		//now from these, calculate next shift
		//K = dy = dx*(dy/dx) = dx*z
		K[y][a] =  dx*YC[z];
		//L = dz = d(dy/dx) = dx*(d^2y/dx^2) = dx*[-y^n - 2(dy/dx)/x]
		K[z][a] = -dx*( YCN.real() + 2.0*YC[z]/YC[x] );	
		if(YC[x]==0) K[z][a] = -dx/3.0; //at very center, need analytic expression for dz/dx
		//calculate "corrected" positions using previous shift vectors
		YC[x] = yin[x] + B[a]*dx;
		for(std::size_t b=1; b<numvar; b++)
			YC[b] = yin[b] + B[a]*K[b][a];
	}		
	yout[x] = yin[x] + dx;
	for(std::size_t b=1; b<numvar; b++)
		yout[b] = yin[b] + K[b][0]/6.0 + K[b][1]/3.0 + K[b][2]/3.0 + K[b][3]/6.0;	
}

//integrate the polytrope up to Len using RK4
//the return is used to refine choice of dx to ensure the model ends at surface
double Polytrope::RK4integrate(const std::size_t Len, double dx){

	//set initial conditions
	centerInit(Y[0]);

	for(std::size_t X = 0; X<Len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
		//sometimes, for noninteger n, at surface values y<0 lead to nans
		//if these occur, it is safe to terminate the integration, as we found surface
		if(Y[X+1][y]<0.0) {return Y[X+1][y];}
	}
	return Y[Len-1][y];
}

//integrate the polytrope up to Len using RK4
//to be used after the correct value of dx has been found 
std::size_t Polytrope::RK4integrate(const std::size_t Len, double dx, int grid){
	grid=1;

	//set initial conditions
	centerInit(Y[0]);

	//integrastr with RK4
	for(std::size_t X = 0; X<Len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
	}
	return Len;
}

//Uses bisection method to find edge
void Polytrope::findEdge(std::size_t Nedge){
	std::size_t R = len-1;
	double dxmax, dxmin;
	double fmax, fmin, f1, f2;
	
	double dx1 = (Y[R-Nedge-1][x]-Y[R-Nedge-2][x])/Nedge;
	double dx2 = 1.1*dx1;
	
	f1 = integrateEdge(Nedge, dx1);
	f2 = integrateEdge(Nedge, dx2);
	
	//find brackets
	double incr = (fabs(f2)>fabs(f1)? -1.1*dx1 : 1.1*dx1);
	//while the two Ws are on same side of axis, keep increasing s and trying again
	std::size_t stop = 0, bnd = 20;
	while(f1*f2 >= 0.0){
		dx2 += incr;
		f2 = integrateEdge(Nedge, dx2);
		if( (++stop) > bnd | dx2<0.0 ) {
			dx2 = 1.03*dx1;
			incr = -incr; stop = 0; bnd = 2*bnd;
		}
	}
	//if we were decreasing dx, then dx2 is min bracket and dx1 is max bracket
	if(incr < 0){
		dxmin = dx2; fmin = f2;
		dxmax = dx1; fmax = f1;
	}
	//if we were increasing dx, then dx2 is max bracket and dx2 is min bracket
	else{
		dxmin = dx1; fmin = f1;
		dxmax = dx2; fmax = f2;	
	}//
	if(dxmin > dxmax){
		double temp = dxmax;
		dxmax = dxmin; dxmin = temp;
		temp = fmax;
		fmax = fmin; fmin = temp;
	}
	
	//now find solution
	dx1 =  0.5*(dxmin + dxmax);
	f2 = integrateEdge(Nedge, dx1);
	stop=0;
	while( fabs(f2-f1) >0.0 ){
		f1 = f2;
		if( f1*fmax > 0.0 ){
			if( dx1 < dxmax ){
				dxmax = dx1;
				fmax = f1;
			}
		}
		else if( f1*fmin > 0.0 ){
			if( dx1 > dxmin ){
				dxmin = dx1;
				fmin = f1;
			}
		}
		dx1 = 0.5*(dxmin+dxmax);
		f2 = integrateEdge(Nedge, dx1);
	}	
}

double Polytrope::integrateEdge(std::size_t Nedge, double dx){
	for(std::size_t X=len-2-Nedge; X<len-1; X++){
		RK4step(dx, Y[X],Y[X+1]);
	}
	return Y[len-1][y];
}

//Here we define functions to access radius, pressure, etc.
double Polytrope::rad(std::size_t X){
	return (Rn*Y[X][x]);
}
double Polytrope::rho(std::size_t X){
	return rho0*Y[X][y]*base[X];
}
double Polytrope::drhodr(std::size_t X){
	return n*rho0*base[X]*(Y[X][z]/Rn);
}
double Polytrope::P(std::size_t X){
	return P0*base[X]*Y[X][y]*Y[X][y];
}
double Polytrope::dPdr(std::size_t X){
	return (n+1.)*P0*base[X]*Y[X][y]*(Y[X][z]/Rn);
}
double Polytrope::Phi(std::size_t X){
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return (n+1.)*P0/rho0 *(Y[len-1][x]*Y[len-1][z]-Y[X][y]);
}
double Polytrope::dPhidr(std::size_t X){
	return -(n+1.)*P0/rho0*Y[X][z]/Rn;
}
double Polytrope::mr(std::size_t X){
	//if Rn, rh0 are changed, this will have correct units
	return pow(Rn,3)*rho0*mass[X];
}

double Polytrope::Schwarzschild_A(std::size_t X, double GamPert){
	if(GamPert==0.0) return Y[X][z]/Y[X][y]/Rn* (n -(n+1.)/Gamma);
	else        	 return Y[X][z]/Y[X][y]/Rn* (n -(n+1.)/GamPert);
}

double Polytrope::getAstar(std::size_t X, double GamPert){
	if(GamPert==0.0) return Y[X][x]*Y[X][z]/Y[X][y]* ((n+1.)/Gamma   - n);
	else        	 return Y[X][x]*Y[X][z]/Y[X][y]* ((n+1.)/GamPert - n);
}

double Polytrope::getU(std::size_t X){
	if(X==0) return 3.0;
	return - Y[X][x]*base[X]*Y[X][y]/Y[X][z];
}

double Polytrope::getVg(std::size_t X, double GamPert){
	if(GamPert==0.0) return -Y[X][x]*Y[X][z]/Y[X][y]* (n+1.)/Gamma;
	else			 return -Y[X][x]*Y[X][z]/Y[X][y]* (n+1.)/GamPert;
}

double Polytrope::getC(std::size_t X){
	if(X==0) return -3.*Y[len-1][z]/Y[len-1][x];
	return (Y[len-1][z]/Y[X][z]) * (Y[X][x]/Y[len-1][x]);
}

double Polytrope::Gamma1(std::size_t X){
	return Gamma;
}

double Polytrope::sound_speed2(std::size_t X, double GamPert){
	if(GamPert == 0.0) return Gamma  *P0*Y[X][y]/rho0;
	else               return GamPert*P0*Y[X][y]/rho0;
}


double Polytrope::Radius(){return rad(len-1);}	//total radius
double Polytrope::Mass(){return mr(len-1);}//total mass


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void Polytrope::setupCenter(){
	//these are coefficients y(x) = 1 + th[1]xi^2 + th[2]xi^4 + th[3] xi^6 + ...
	ac[0] = 1.0; ac[1] = -1.0/6.0; ac[2] = n/120.; ac[3] = n*(5.-8.*n)/15120.;
}

void Polytrope::getAstarCenter(double *AC, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][x];
	double AV = (n*Gam1)/(n+1.) - 1.;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=0) AC[0] = 0.0;
	if(maxPow>=2) AC[1] = AV*(n+1.)/(3.*Gam1)*pow(x1,2);
	if(maxPow>=4) AC[2] = AV*2.*(n+1.)*(1./36.-n/60.)/Gam1*pow(x1,4);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void Polytrope::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][x];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = (n+1.)/(3.*Gam1)*pow(x1,2);
	if(maxPow>=4) Vc[2] = 2.*(n+1.)*(1./36.-n/60.)/Gam1*pow(x1,4);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void Polytrope::getUCenter(double *Uc, int& maxPow){
	double x1 = Y[len-1][x];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = -n/5.*pow(x1,2);
	if(maxPow>=4) Uc[2] = 54.*n*(7.*n/8100. - 1./1296. + (5.-8.*n)/15120.)*pow(x1,4);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void Polytrope::getC1Center(double *cc, int& maxPow){
	double x1 = Y[len-1][x];
	double z1 = Y[len-1][z];
	if(maxPow>=0) cc[0] = -3.*z1/x1;
	if(maxPow>=2) cc[1] = 3.*n*z1/(10.*x1)*pow(x1,2);
	if(maxPow>=4) cc[2] = n*(2.*n+25.)*z1/(1400.0*x1)*pow(x1,4);
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
void Polytrope::setupSurface(){
	double x1 = Y[len-1][x];
	double a1 =-Y[len-1][z]*x1;
	//first, we must find the coefficients in expansion theta = a1*t+a2*t^2 + ...
	as[0] = 0.0;
	for(int i=1; i<6;i++) as[i] = a1;
	if(n==0){
		as[1] = a1 = 2.0;
		as[2] = -1.0;
		as[3] = as[4] = as[5] = 0.0;
	}
	else if(n==1){
		as[3] -= as[1]*x1*x1/6.;
		as[4] -= as[1]*x1*x1/6.;
		as[5] -= as[1]*x1*x1/6.*(1.-x1*x1/20.);
	}
	else if(n==2){
		as[4] -= as[1]*as[1]*x1*x1/12.;
		as[5] -= 2.*as[1]*as[1]*x1*x1/15.;
	}
	else if(n==3){
		as[5] -= as[1]*as[1]*as[1]*x1*x1/20.;
	}
}

void Polytrope::getAstarSurface(double *As, int& maxPow, double g){
	//we make use of the fact  that A* and Vg are simply related in polytropes
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double AV = (n*Gam1)/(n+1.) - 1.;
	int O=1;
	getVgSurface(As, maxPow, g);
	for(int k=-1; k<=maxPow; k++){
		As[O+k] = AV*As[O+k];
	}
}

void Polytrope::getVgSurface(double *Vs, int& maxPow, double g){
// coefficients of Vg must extend up to maxPow-1
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][x];
	double a1 =-Y[len-1][z]*x1;
	int O=1;
	if(maxPow>= -1) Vs[O-1] = (n+1.)/Gam1;
	if(maxPow>=  0) Vs[O  ] = (n+1.)/Gam1*(as[2]/a1 - 1.);
	if(maxPow>=  1) Vs[O+1] = (n+1.)/Gam1*(2.*as[1]*as[3]-as[2]*as[1]-as[2]*as[2])/(a1*a1);
	if(maxPow>=  2) Vs[O+2] = (n+1.)/Gam1*( pow(as[2]/a1,3) + pow(as[2]/a1,2) -2.*as[3]/a1 -3.*as[2]*as[3]/(a1*a1)+3.*as[4]/a1);
	if(maxPow>=  3) Vs[O+3] = (n+1.)/Gam1*(4.*as[5]/a1 - 4.*as[2]*as[4]/(a1*a1) - 3.*as[4]/a1 - 2.*pow(as[3]/a1,2) + 4.*pow(as[2]/a1,2)*as[3]/a1 +3.*as[2]*as[3]/(a1*a1) - pow(as[2]/a1,4) - pow(as[2]/a1,3));
	//if more  terms than this requested, cap number of terms
	if(maxPow > 3) maxPow = O+3;
}

void Polytrope::getUSurface(double *Us, int& maxPow){
// coefficients of U must extend up to order maxPow
	double x1 =  Y[len-1][x];
	double a1 = -Y[len-1][z]*x1;	
	for(int a=0; a<=maxPow; a++) Us[a] = 0.0;
	
	if(n==0){
		if(maxPow>=0) Us[0] =  3.0;
	}
	else if(n==1){
		if(maxPow>=1) Us[1] =  x1*x1;
		if(maxPow>=2) Us[2] = -2.0*x1*x1;
		if(maxPow>=3) Us[3] =  x1*x1*(1.0 + x1*x1/3.0);
		if(maxPow>=4) Us[4] = -x1*x1*x1*x1;
	}
	else if(n==2){
		if(maxPow>=2) Us[2] =  a1*x1*x1;
		if(maxPow>=3) Us[3] = -a1*x1*x1;
		if(maxPow>=4) Us[4] =  0.0;
	}
	else if(n==3){
		if(maxPow>=3) Us[3] = a1*a1*x1*x1;
		if(maxPow>=4) Us[4] = 0.0;
	}
	else if(n==4){
		if(maxPow>=4) Us[4] = a1*a1*a1*x1*x1;
	}
	//if more  terms than this requested, cap number of terms
	if(maxPow > 4) maxPow = 4;
}

void Polytrope::getC1Surface(double *cs, int& maxPow){
// coefficients of c1 are only needed up to order maxPow-1
	double x1 = Y[len-1][x];
	double z1 = Y[len-1][z];
	if(maxPow>=0) cs[0]  =  1.;
	if(maxPow>=1) cs[1]  = -1. - 2.*as[2]/as[1];
	if(maxPow>=2) cs[2]  =  2.*as[2]/as[1] + 4.*pow(as[2]/as[1],2) - 3.*as[3]/as[1];
	if(maxPow>=3) cs[3]  = -4.*pow(as[2]/as[1],2)-8.*pow(as[2]/as[1],3) + 3.*as[3]/as[1] + 12.*as[2]*as[3]/(as[1]*as[1]) - 4.*as[4]/as[1];
	if(maxPow>=4) cs[4]  =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}


void Polytrope::writeStar(const char *const c){
	//create names for files to be opened
	std::string pathname;
	if(c==NULL)	pathname = addstring("./out/", name);
	else pathname = addstring("./",c)+"/star/";

	system( addstring("mkdir -p ", pathname).c_str() );
	
	printStar(pathname.c_str());
	printBV(pathname.c_str(), 5./3.);
	printCoefficients(pathname.c_str(), 5./3.);
}


int Polytrope::read_star_input(FILE* input_file, Calculation::InputData& calcdata){
 	char input_buffer[512];
	std::string instring;
	
	calcdata.input_params.reserve(2);
	double index=0.0, grid_size=0.0;
	fscanf(input_file, " %lf", &index);	//read in the index
	fscanf(input_file, "%lf\n", &grid_size);	//read in the grid size
	calcdata.input_params.push_back(index);
	calcdata.input_params.push_back(grid_size);
	printf("n=%lf\n", calcdata.input_params[0]);
	calcdata.Ngrid = int(calcdata.input_params[1]);
	
	//now read in desired physical properties of star -- specify two at a time
	//for a polytrope, can specify two of mass, radius, logg, or surface z
	double temp;
	//read in first physical parameter -- name as string, value as double
	fscanf(input_file, "Params: %s %lf\t", input_buffer, &temp);
	instring=std::string(input_buffer);
	//save value in appropriate slot use bit-masking to keep track of variables
	if(     !instring.compare("mass"))   {calcdata.mass = temp; calcdata.params|=units::ParamType::pmass;}
	else if(!instring.compare("radius")) {calcdata.radius=temp; calcdata.params|=units::ParamType::pradius;}
	else if(!instring.compare("zsurf"))  {calcdata.zsurf =temp; calcdata.params|=units::ParamType::pzsurf;}
	else if(!instring.compare("logg"))   {calcdata.logg  =temp; calcdata.params|=units::ParamType::plogg;}
	char tempparam = calcdata.params;	//value after first read-in, for comparison later
	//read in second physical parameters -- name as string, value as double
	fscanf(input_file, "%s %lf\n", input_buffer, &temp);
	instring=std::string(input_buffer);
	//save valeu in appropriate slot, again use bit-masking so that binary value of params indicates which were used
	if(     !instring.compare("mass"))   {calcdata.mass = temp; calcdata.params|=units::ParamType::pmass;}
	else if(!instring.compare("radius")) {calcdata.radius=temp; calcdata.params|=units::ParamType::pradius;}
	else if(!instring.compare("zsurf"))  {calcdata.zsurf =temp; calcdata.params|=units::ParamType::pzsurf;}
	else if(!instring.compare("logg"))   {calcdata.logg  =temp; calcdata.params|=units::ParamType::plogg;}
	else {
		printf("ERROR IN INPUT: invalid parameter to polytope:\n allowed mass, radius, zsurf, or logg\n");
		printf("This error is fatal.  Quitting.\n");
		return 1;
	}
	//make sure we did not specify the same parameter twice
	if(calcdata.params == tempparam){
		printf("ERROR IN INPUT: double specification of stellar parameters, star underdefined\n");
		printf("This error is fatal.  Quitting.\n");
		switch(units::ParamType(calcdata.params)){
			case units::ParamType::pmass:   printf("\tmass given twice\n"); break;
			case units::ParamType::pradius: printf("\tradius given twice\n"); break;
			case units::ParamType::pzsurf:  printf("\tzsurf given twice\n"); break;
			case units::ParamType::plogg:   printf("\tlog g given twice\n"); break;
			case units::ParamType::pteff:   printf("\tteff given twice\n"); break;
		}
		return 1;
	}	
	if((calcdata.params&units::ParamType::pteff)){
		printf("ERROR IN INPUT: polytropes cannot be assigned temperature\n");
		printf("This error is fatal.  Quitting.\n");
		return 1;
	}

	return 0;
}


#endif