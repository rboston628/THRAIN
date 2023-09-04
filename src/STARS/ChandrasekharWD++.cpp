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
#include "../lib/rootfind.h"
#include "../src/ThrainMain.h"


void ChandrasekharWD::basic_setup(){
	//exclude unphysical values of Y0
	if(Y0 < 1.) Y0 = 1.;
	Y02 = Y0*Y0;
	X02 = Y02-1.;
	X0  = sqrt(X02);
	//name this model for files
	name = strmakef("ChandrasekharWD.%1.3lf", Y0);

	Rn = sqrt(2.*A0/(m_pi*Gee()))/B0;

	// prepare stellar arrays
	Y = new double*[len];
	for(std::size_t i=0; i<len; i++){
		Y[i] = new double[numvar];
	}
}

void ChandrasekharWD::init_arrays(){
	mass = new double[len];
	x3   = new double[len];
	mue  = new double[len];
	dmue = new double[len];
	for(int X=0; X<len; X++){
		Y[X][x] = sqrt(std::complex<double>(Y[X][y]*Y[X][y]-1.0)).real();
		Y[X][f] = Chandrasekhar::factor_f(Y[X][x]);
		chemical_gradient(Y[X][x], Y[X][z], mue[X], dmue[X]);
		mass[X] = -4.*m_pi*Y[X][xi]*Y[X][xi]*Y[X][z]/mue[X];
		x3[X] = pow(Y[X][x],3);
	}
}

//initalize white dwarf from central value of y and length
ChandrasekharWD::ChandrasekharWD( double Y0, std::size_t L, ChemicalGrad chem)
	: Y0(Y0), len(L), 
	chemical_gradient(chem),
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0)
{
	basic_setup();

	//we find an appropriate grid spacing for array holding star data
	//we need for dx to be such that the integration ends with x[len-1]=0.0
	//we will find the proper dx with a bisection search

	double dxi = 1./double(len); //2.*sqrt(6.0)/(len-1); //an initial guess from constant volume star
	double dximax=1.0, dximin=0, yS;
	std::function<double(double)> surface = [this](double step) -> double {
		return RK4integrate(step, SurfaceBehavior::STOP_AT_ZERO);
	};
	rootfind::bisection_find_brackets_newton(surface, dxi, dximin, dximax);
	yS = rootfind::bisection_search(surface, dxi, dximin, dximax);
	RK4integrate(dxi, SurfaceBehavior::CONTINUE_FULL_LENGTH);
	printf("edge  =\t%le %le\n", yS, Y[len-1][y]-1.0);

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = len/2;
	for(std::size_t X=1; X<len-1; X++){
		//scan through x, set matching point where x[X] = 0.5
		if(Y[X-1][x]>0.5 & Y[X+1][x]<=0.5) indexFit = X;
	}
	indexFit /= 2;
	setupCenter();
	setupSurface();
	printf("%0.8lf\t%le\t%le\n", 1./Y02, Mass()/MSOLAR, Radius());
}

//initalize white dwarf from central value of y and length, with specific step size
//this should only be used for testing scaling issues
ChandrasekharWD::ChandrasekharWD( double Y0, std::size_t L, const double dxi, ChemicalGrad chem)
	: Y0(Y0), len(L), 
	A0(Chandrasekhar::A0), B0(Chandrasekhar::B0),
	chemical_gradient(chem)
	// mu0(mu0), k(k), acore(acore), aswap(aswap)
{	
	basic_setup();
	//we know dx, so no need for bisection search
	RK4integrate(dxi, SurfaceBehavior::CONTINUE_FULL_LENGTH);

	//now set physical properties of the white dwarf
	init_arrays();
	indexFit = 512*round(double(len)/1024.0);
	indexFit /= 2;
	printf("  indexFit  = %lu\n", indexFit);
	printf("r[indexFit] = %0.32le\n", rad(indexFit));
	setupCenter();
	setupSurface();
	printf("%0.2lf\t%le\t%le\n", 1./Y02, mr(len-1)/MSOLAR, rad(len-1));
}

ChandrasekharWD::~ChandrasekharWD(){
	for(std::size_t i=0; i<len; i++){
		delete[] Y[i];
	}
	delete[] mass;
	delete[] mue;
	delete[] dmue;
	delete[] x3;
}

void ChandrasekharWD::centerInit(double ycenter[numvar]){
	ycenter[xi] = 0.0;
	ycenter[y]  = Y0;
	ycenter[x]  = X0;
	ycenter[z]  = 0.0;
	ycenter[f]  = Chandrasekhar::factor_f(X0);
}

void ChandrasekharWD::RK4step(double dxi, double yin[numvar], double yout[numvar]){
	double YC[numvar] = {0.};
	double K[4][numvar] = {0.};
	static const double B[4] = {0.5,0.5,1.0,0.0};
	static const double d[4] = {6.0,3.0,3.0,6.0};

	std::complex<double> YCN(0.0,0.0);
	double MUC = 1.0, DMUC = 0.0;

	for(VarName b = xi; b<numvar; b = VarName(b+1)){
		YC[b] = yin[b];
	}
	chemical_gradient(YC[x], YC[z], MUC, DMUC);
	for(int a = 0; a<4; a++){
		//find current calues
		YCN = YC[y]*YC[y] - 1.;
		YCN = pow(YCN, 1.5);
		K[a][xi] = dxi;
		K[a][y] =  dxi*YC[z];
		K[a][z] = -dxi*( YCN.real()*MUC*MUC + 2.0*YC[z]/YC[xi] -  YC[z]*DMUC/MUC );
		if(YC[xi]==0) K[a][z] = -dxi/3.0*YCN.real()*MUC*MUC; //must use analytic expression at very center
		//calculate "corrected" positions using previous shift vectors
		for(VarName b=xi; b<=z; b=VarName(b+1)){
			YC[b] = yin[b] + B[a]*K[a][b];
		}
		if(YC[y]>=1.0) YC[x] = sqrt(YC[y]*YC[y]-1.0);
		chemical_gradient(YC[x],YC[z], MUC, DMUC);
	}
	yout[xi] = yin[xi] + dxi;
	for(VarName b=y; b<=z; b=VarName(b+1)){
		yout[b] = yin[b];
		for(int a=0; a<4; a++) yout[b] += K[a][b]/d[a];
	}
	yout[x] = sqrt(yout[y]*yout[y]-1.);
	if(yout[y]<1.) yout[x] = 0.0;
}

//integrate the polytrope up to Len using RK4
double ChandrasekharWD::RK4integrate(double dx, SurfaceBehavior sb){
	//set our initial conditions
	centerInit(Y[0]);

	double ysurface = Y0;
	std::size_t X;
	for( X = 0; X<len-1; X++){
		RK4step(dx, Y[X], Y[X+1]);
		//sometimes, for noninteger n, if array is too large values y<0 lead to nans
		//if these occur, it is safe to terminate the integration
		if(Y[X+1][y]<=1.0 && bool(sb)) break;
	}
	if(X<len-1)	ysurface = Y[X+1][y]   - 1.;
	else        ysurface = Y[len-1][y] - 1.;
	return ysurface;
}

//Here we define functions to access radius, pressure, etc.
double ChandrasekharWD::rad(std::size_t X){
	return Rn*Y[X][xi];
}
double ChandrasekharWD::rho(std::size_t X){
	using namespace Chandrasekhar;
	return B0*mue[X]*x3[X];
}
double ChandrasekharWD::drhodr(std::size_t X){
	using namespace Chandrasekhar;
	return 3.*B0*mue[X]*Y[X][x]*Y[X][y]*Y[X][z]/Rn + B0*x3[X]*dmue[X]/Rn; //note: dx/dxi = z*y/x
}
double ChandrasekharWD::P(std::size_t X){
	using namespace Chandrasekhar;
	return A0*Y[X][f];
}
double ChandrasekharWD::dPdr(std::size_t X){
	using namespace Chandrasekhar;
	return 8.*A0*x3[X]*Y[X][z]/Rn; //see eqn (9) in Chandrasekhar 1935
}
double ChandrasekharWD::Phi(std::size_t X){
	using namespace Chandrasekhar;
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return -8.*A0/B0 * (Y[X][y] - 1. + Y[len-1][z]*Y[len-1][xi])/mue[X];
}
double ChandrasekharWD::dPhidr(std::size_t X){
	using namespace Chandrasekhar;
	return -8.*A0/B0 * Y[X][z]/Rn/mue[X];
}
double ChandrasekharWD::mr(std::size_t X){
	using namespace Chandrasekhar;
	return B0*pow(Rn,3)*mass[X];
}

double ChandrasekharWD::Schwarzschild_A(std::size_t X, double GamPert){
	if(GamPert==0.0) return dmue[X]/mue[X]/Rn;
	else        	 return dmue[X]/mue[X]/Rn
							- Y[X][z]*(8.*x3[X]/Y[X][f]/GamPert-3.*Y[X][y]/pow(Y[X][x],2))/Rn;
}

double ChandrasekharWD::getAstar(std::size_t X, double GamPert){
	if(GamPert==0.0) return -dmue[X]/mue[X]*Y[X][xi];
	else        	 return -dmue[X]/mue[X]*Y[X][xi] 
							+ Y[X][z]*Y[X][xi]*(8.*x3[X]/Y[X][f]/GamPert-3.*Y[X][y]/pow(Y[X][x],2));
}

double ChandrasekharWD::Ledoux(std::size_t X, double GamPert){
	return dmue[X]/mue[X]/Rn;
}

double ChandrasekharWD::getU(std::size_t X){
	if(X==0) return 3.0;
	return -mue[X]*mue[X]*Y[X][xi]*x3[X]/Y[X][z];
}

double ChandrasekharWD::getVg(std::size_t X, double GamPert){
	if(GamPert==0.0) return -8.*x3[X]*Y[X][xi]*Y[X][z]/Y[X][f]/Gamma1(X);
	else			 return -8.*x3[X]*Y[X][xi]*Y[X][z]/Y[X][f]/GamPert;
}

double ChandrasekharWD::getC(std::size_t X){
	if(X==0) return Y[len-1][z]*Y[len-1][xi]/2./yc[1]*mue[0]/mue[len-1];
	return (Y[len-1][z]/Y[X][z])*(Y[X][xi]/Y[len-1][xi])*(mue[X]/mue[len-1]);
}

double ChandrasekharWD::Gamma1(std::size_t X){
	if(Y[X][x]==0.) return 0.0;
	return 8.*pow(Y[X][x],5)/(3.*Y[X][y]*Y[X][f]);
}

double ChandrasekharWD::sound_speed2(std::size_t X, double GamPert){
	using namespace Chandrasekhar;
	if(GamPert == 0.0) return Gamma1(X)*A0/B0*Y[X][f]/x3[X]/mue[X];
	else               return GamPert  *A0/B0*Y[X][f]/x3[X]/mue[X];
}


double ChandrasekharWD::Radius(){return Rn*Y[len-1][xi];}	//total radius
double ChandrasekharWD::Mass(){return mr(len-1);}//total mass
double ChandrasekharWD::Gee(){return G_CGS;}

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
void ChandrasekharWD::setupCenter(){	
	//the central coefficients -- expanded in terms of r/R
	double z1 = Y[len-1][z];
	double x1 = Y[len-1][xi];
	double mu2 = mue[0]*mue[0];
	//
	yc[0]=  Y0;
	yc[1]= -X02*X0/6.*pow(x1,2)*mu2;
	yc[2]=  X02*X02*Y0/40.*pow(x1,4)*mu2*mu2;
	yc[3]=  X0*(-14.*X02*X02-19.*X02*X02*X02)/5040.*pow(x1,6)*pow(mu2,3);
	//
	xc[0]=  X0;
	xc[1]= -X02*Y0/6.*pow(x1,2)*mu2;
	xc[2]=  X02*X0*(4.+9.*X02)/360.*pow(x1,4);
	//
	fc[0]=  Chandrasekhar::factor_f(X0);
	fc[1]= -X02*X0*(2. - 3.*X0 + 3.*X02)/3.*pow(x1,2);
}

void ChandrasekharWD::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double xi2 = Y[len-1][xi]*Y[len-1][xi];
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
	double x1 = Y[len-1][xi];
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -16.*X02*X0*yc[1]/Gam1/fc[0];
	if(maxPow>=4) Vc[2] =-8.*   pow(x1,2)*pow(X0,6)*(5.*fc[1]+4.*pow(x1,2)*fc[0]*X0*Y0)/(15.*Gam1*fc[0]*fc[0]);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4; 
}

void ChandrasekharWD::getUCenter(double *Uc, int& maxPow){
	double x1 = Y[len-1][xi];
	double mu2 = mue[0]*mue[0];
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] =-0.6*mu2*x1*x1*X0*Y0;
	if(maxPow>=4) Uc[2] = 
		9.*xc[2]/X0 - pow(x1*x1*mu2,2)*X02*(0.08+0.134*X02);
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}

void ChandrasekharWD::getC1Center(double *cc, int& maxPow){
	double x1 = Y[len-1][xi];
	double z1 = Y[len-1][z];
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
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][xi];
	double a1 =-Y[len-1][z]*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) As[O-1] = 1.5;
	if(maxPow>= 0) As[O  ] = 0.75*a1;
	if(maxPow>= 1) As[O+1] = 0.75*a1 + 8.*a1*a1/Gam1 - 3.*a1*a1/8.;
	if(maxPow>= 2) As[O+2] = a1*3./16.*pow(a1-2.,2)+4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) As[O+3] = -3.*a1*pow(a1-2,3)/32.+(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow> 3) maxPow = O+3;
}

void ChandrasekharWD::getVgSurface(double *Vs, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	double x1 = Y[len-1][xi];
	double a1 =-Y[len-1][z]*x1;
	int O=1;
	//depending on power requested, return appropriate number of terms
	if(maxPow>=-1) Vs[O-1] = 0.0;
	if(maxPow>= 0) Vs[O  ] = 0.0;
	if(maxPow>= 1) Vs[O+1] = -8.*a1*a1/Gam1;
	if(maxPow>= 2) Vs[O+2] = -4./3.*a1*a1*(12.+5.*a1)/Gam1;
	if(maxPow>= 3) Vs[O+3] = -(724.*a1*a1/45. + 20.*a1 + 24.)*a1*a1/Gam1;
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
		fprintf(fp, "%0.16le\t", Y[X][xi]);
		fprintf(fp, "%le\t%le\t%le", Y[X][x], Y[X][y], Y[X][f]);
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
			Y[X][xi],
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


int ChandrasekharWD::read_star_input(Calculation::InputData& calcdata, FILE* input_file){
	char input_buffer[512];
	std::string instring;
	
	calcdata.input_params.resize(3);
	fscanf(input_file, " %lf", &calcdata.input_params[0]);	//read in the central y value
	printf("y0=%lf\n", calcdata.input_params[0]);
	//read in an integer designating the chemical composition to use for mu
	fscanf(input_file, " %lf", &calcdata.input_params[1]);
	switch((int) calcdata.input_params[1]){
		case 0:
			printf("standard Chandrasekhar WD\n");
			break;
		case 1:
			printf("using logistic composition\n");
			break;
	}
	fscanf(input_file, "%lf\n", &calcdata.input_params[2]);	//read in the grid size
	calcdata.Ngrid = std::size_t(calcdata.input_params[2]);

	return 0;
}

// TODO echo, write methods


#endif