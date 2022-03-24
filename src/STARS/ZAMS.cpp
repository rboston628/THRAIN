//**************************************************************************************
//							ZERO AGE MAIN SEQUENCE STAR
// ZAMS.cpp
//	A Newtonian star at the beginning of the main sequence lifetime
//	Chemically homogeneous (X,Y,Z)
//  The interior equation of state includes pressure contributions from
//		P = P_ions + P_rad
//	Integrates with Schwarzschild's dimensionless logarithmic variables
//  This code draws on inspiration from:
//		ZAMS.f by C.J. Hansen & S.D. Kawaler -- http://astro.if.ufrgs.br/evol/evolve/hansen/
//  	StatStar by Bradley W. Carroll & Dale A. Ostlie,  2007.
//  Reece Boston, Mar 24, 2022
// **************************************************************************************/

#ifndef ZAMSCLASS
#define ZAMSCLASS

#include "ZAMS.h"
#include "Polytrope.cpp"

double radiative_opacity(const StellarVar& ly, const Abundance& X){
//the below is from Hansen & Kawaler (also found in Shapiro & Teukolsky 1983 and Schwarzschild 1958)
	double meanZA = X.H1 + X.He4 + 3.*X.C12 + 4.*X.O16; // average <Z^2/A>
	//electron (Thomson) scattering
	double ke = 0.195*(1.+X.H1);
	//see Schwarzschild 1958 (pg 237 or eq16.2), Schwarzschild 1946 (eq9), and Shapiro & Teuskolsky (eq4.1.8)
	//note assumption of guillotine (Gaunt) factor g ~ 10
	double kbf = 4.34e24*(X.C12+X.O16)*(1.+X.H1)*exp(ly[dens]-3.5*ly[temp]);
	//see Schearzschild 1858, or Hansen and Kawaler eq4.37
	double kff = 3.68e22*(X.H1 +X.He4 + (X.C12+X.O16)*meanZA)*(1.+X.H1)*exp(ly[dens]-3.5*ly[temp]);
	//from Iben 1975, Appendix B
	double t6 = ly[temp]-6.*log(10.0);
	double C  = pow(2.019 + 1e-4*exp(ly[dens]-1.7*t6),2.425);
	double Ap = 1. + C*(1.+C/24.55);
	double mu = meanZA - 1.;
	double B  =  3.86 + 0.252*sqrt(mu) + 0.018*mu;
	double A  = (1.25 + 0.488*sqrt(mu) + 0.092*mu)/0.67;
	double lr = -A + B*t6;
	double k = Ap*exp(0.67*(ly[dens]-lr));
	return k + ke + kbf + kff;
}

double conductive_opacity(const StellarVar& ly, const Abundance& X){
// below is based on Cassis 2007, which relies on Yakovlev & Urpin 1980 and Potekhin et al 1999
	double rho = exp(ly[dens]);
	double x = pow(rho/Chandrasekhar::B0/X.mu_e(), 1./3.);
	//the Coulomb factor, which also controls crystallization
	double GC = 2.275e5*(X.mean_coulomb())*exp(ly[dens]/3.-ly[temp]);
	//Coloumb logarithm from Yakovlev&Urpin 1980, eq 11, vastly simplified from Potekhin 1999
	double Coulomb_Logarithm = log(2.*m_pi*X.mean_Z()/3.)/3. + log(1.5+3./GC)/2. - 0.5*x*x/(1.+x*x); 
	// now get opacity, by combining Cassisi 2007 eq 1, 2, 8
	double opaccoeff = 512.*boltzmann_sigma*pow(electron.mass_CGS,2)*pow(electron.charge_CGS,4)
				*proton.mass_CGS/(3.*pow(plank_h_CGS,3)*pow(boltzmann_k,2));
	double opacity = opaccoeff*X.mean_A()*sqrt(1.+x*x)*exp(2.*(ly[temp]-ly[dens]))*Coulomb_Logarithm;
	return opacity;//*/
}

double ZAMS::opacity(const StellarVar& ly, const Abundance& X){
	double k_rad = radiative_opacity(ly,X);
	double k_cond= conductive_opacity(ly,X);
	return 1./(1./k_rad+1./k_cond);
}

//initalize ZAMS model
ZAMS::ZAMS(
	double M,   //mass, in solar masses
	double R,   //radius, in solar radii
	double Teff,//effective temperature, in kelvin
	int Len,
	double X,
	double Y
)	: Msolar(M), Rsolar(R), Teff(Teff), Ntot(Len), Xtot(X,Y,(1.-X-Y)/2.,(1.-X-Y)/2.)
{
	dXelem = Abundance(0.,0.,0.,0.);
	//begin by assigning values to EOS and chemical abundance based on file
	setup();
	//prepare a starting model from a ChandrasekharWD, to guess R, P0
	initFromPolytrope();
	
	//set indices for different regions
	//the boundary of the core
	Ncore = (int) (0.7*double(Ntot));  //the core uses a fifth of the grid points
	if( Ncore%2==0) Ncore = Ncore+1;
	//the boundary of atmosphre
	Natm = (int) (0.3*Ntot);   //the atmosphere and envelope use the remaining
	if( Natm%2==1) Natm = Natm+1;
	Ntot = Ncore + Natm;
	double Mcore = 0.99;
	double Rcore = 0.90;
	//initilize the log-radial grid
	logQ = new double[Ntot];	//logQ could be either r or m
	setupGrid(Mcore, Ncore);
	//initialize equilibrium structure arrays
	logY = new StellarVar[Ntot];
	dlogY= new StellarVar[Ntot];
	Xmass  = massFraction();
	
	printf("star  :\t M=%le\tR=%le\tL=%le\n", Msolar, Rsolar, Lsolar);
	sprintf(name, "WD_M%0.2f_L%0.2f_X%0.2f_Y%0.2f", Msolar, Rsolar, Xtot[0], Xtot[1]);
	// find relevant scales for re-dimensionalizing the quantities
	Dscale = Mstar*pow(Rstar,-3)/(4.*m_pi);
	Pscale = G_CGS*pow(Mstar,2)*pow(Rstar,-4)/(4.*m_pi);
	Tscale = Xtot.mean_A()*proton.mass_CGS/boltzmann_k*G_CGS*Mstar/Rstar;
	printf("scales:\t M=%le\tR=%le\tL=%le\n", Mstar, Rstar, Lstar);
	printf("scales:\t D=%le\tP=%le\tT=%le\n", Dscale, Pscale, Tscale);
	printf("Xtot  :\t %0.8lf %0.8lf %0.8lf %0.8lf\n", Xtot[0], Xtot[1], Xtot[2], Xtot[3]);
	printf("Xmass :\t %0.8lf %0.8lf %0.8lf %0.8lf\n", Xmass[0],Xmass[1],Xmass[2],Xmass[3]);
	
	Yscale = StellarVar(Dscale,Rstar,Pscale,Mstar,Tscale,Lstar);
	logYscale = log(Yscale);
	Ysolar = StellarVar(  0.0 , Rsolar, 0.0 , Msolar, 0.0 , Lsolar);
	Ystar  = StellarVar(  0.0 , Rstar,  0.0 , Mstar,  0.0 , Lstar );
			
	//begin the calculation with constant step size throughout entire star to tune P0,T0
	printf("Converging model to P0, T0, R\n");
	double P0 = Ystart0[pres]/Pscale;//1.e24/Pscale;// //dimensionless
	double T0 = 1.e3*Teff/Tscale; //dimensionless
	double qs = 1.0e-2;
	//variables of the Newtonian gradient-descent method
	double x[numv] = {P0, T0, qs};
	double z[numv] = {0., 0.,-1.e0};
	double f[numv], dx[numv], F=1.0;
	double M1, M2=0.0;
	double dFdx[numv], dfdx[numv][numv];
	double x2[numv], f2[numv], varied[numv], used[numv], F2;
	//first iteration
	int count=1;
	double maxDF = 1.0;
	while( maxDF > 1.0e-10 ){
		printf("ITERATION: %d\n", count++);
		//recompute difference, and Jacobian matrix
		joinAtCenter(x, f, F);
		printf("\tX1:\t"); for(int i=0; i<numv; i++) printf("%le ", x[i]); printf("\n");
		printf("\tf1:\t"); for(int i=0; i<numv; i++) printf("%le ", f[i]); printf("\n");
		printf("\tF1:\t%le\n", F);
		
		//calculate a Jacobian matrix by varying each variable
		for(int i=0; i<numv; i++) varied[i] = 1.01*x[i];
		for(int i=0; i<numv; i++){
			for(int j=0; j<numv; j++) used[j] = x[j];
			used[i] = varied[i];
			joinAtCenter(used, f2, F2);
			for(int j=0; j<numv; j++) dfdx[j][i] = (f2[j]-f[j])/(varied[i]-x[i]);
			dFdx[i] = (F2-F)/(varied[i]-x[i]);
		}
		
		//Newton's algorithm calls for -f
		for(int i=0; i<numv; i++) f2[i] = -f[i];
		
		//save the past gradient in case matrix inversion fails
		double dxsave[numv]; for(int i=0;i<numv;i++) dxsave[i] = dx[i];
		//invert the matrix -- check for errors
		if(invertMatrix(dfdx,f2, dx)){
			//if the matrix is singular or otherwise fails
			// then just do something to try to recover
			printf("ERROR: Matrix inversion failed!\n");
			//use the last gradient
			for(int i=0; i<numv; i++) dx[i] = dxsave[i];
		}

		bool failed = false;
		for(int i=0; i<numv; i++) if(isnan(dx[i])) failed=true;
		//if(failed) for(int i=0; i<numv; i++) dx[i] = dxsave[i];
		bool allzero = true;
		for(int i=0; i<numv; i++) allzero &= (dx[i] == 0.0);
		if(allzero) for(int i=0; i<numv; i++) {
			printf("all zero:\t%le %le\n", dx[0], dx[1]);
			dx[i] = 0.1*x[i];
		}
	
		double L = 1.0;
		for(int i=0; i<numv; i++) x2[i] = x[i] + dx[i];
		bool negative = false;
		for(int i=0; i<numv; i++) if(x2[i]<z[i]) negative=true;
		while(negative){
			for(int i=0; i<numv; i++) dx[i] *= 0.1;
			for(int i=0; i<numv; i++) x2[i] = x[i]+dx[i];
			printf("\tX2:\t"); for(int i=0; i<numv; i++) printf("%le ", x2[i]); printf("\n");
			negative = false;
			for(int i=0; i<numv; i++) if(x2[i]<z[i]) negative=true;	
		}
		
		joinAtCenter(x2, f2, F2);
		int freeze = 0;
		while(F2 > F){
			if(L<1.e-3) {L = 1.e-3; break;}
			L *= 0.1;
			for(int i=0; i<numv; i++) x2[i] = x[i] + L*dx[i];
			joinAtCenter(x2, f2, F2);
		}
		printf("\tdx:\t"); for(int i=0; i<numv; i++) printf("%le ", L*dx[i]); printf("\n");
		F = F2;
		for(int i=0; i<numv; i++) f[i] = f2[i];
		for(int i=0; i<numv; i++) x[i] = x2[i];
		
		printf("\tX2:\t"); for(int i=0; i<numv; i++) printf("%le ", x[i]); printf("\n");
		printf("\tf2:\t"); for(int i=0; i<numv; i++) printf("%le ", f[i]); printf("\n");
		printf("\tF2:\t%le\n", F);
				
		// determine the max size of the differences
		maxDF = -1.0;
		for(int i=0; i<numv; i++) if(fabs(f[i]) > maxDF) maxDF = fabs(f[i]);		
		printf("MAX DIF = %le\n", maxDF);
	}
	
	printf("MODEL CONVERGED!\n");
	joinAtCenter(x,f,F);
	
	//all done!
	populateBruntVaisala();
	indexFit = Ncore/2;
	
	printf("Teff=%le\tTsurf=%le\n", Teff, Tscale*exp(logY[Ntot-1][temp]));
	printf("Xtot :\t%0.8lf %0.8lf %0.8lf %0.8lf\n", Xtot.H1 , Xtot.He4 , Xtot.C12 , Xtot.O16);
}

ZAMS::~ZAMS(){
	delete[] logQ;
	delete[] logY;
	delete[] dlogY;
	delete[] adiabatic_1;
	delete[] nabla;
	delete[] nabla_ad;
	delete[] brunt_vaisala;
	delete[] ledoux;
	delete[] kappa;
}


void ZAMS::setup(){
	FILE *input_file;
	if(!(input_file=fopen("zams.txt", "r"))){
		printf("ERROR: EOS file not found... using a default\n");
		//do something
		std::vector<PartialPressure> corePres{ rad_gas, ideal ,coul, deg_zero };
		std::vector<PartialPressure> atmPres{  rad_gas, ideal };
		for(PartialPressure p : corePres) core_pressure.push_back(p);
		for(PartialPressure p : atmPres)   atm_pressure.push_back(p);
		return;
	}
	size_t buffer_size = 256;
	ssize_t line_size;
	char *input_buffer = NULL, *pressure;
	std::string instring;
	EOS *pres = NULL;
	line_size = getline(&input_buffer, &buffer_size, input_file);
	printf("%s", input_buffer);fflush(stdout);
	while(line_size > 0){
		//fscanf(input_file, "%s\n", input_buffer);
		line_size = getline(&input_buffer, &buffer_size, input_file);
		if(line_size >1) printf("%s", input_buffer);
		if(     !strcmp(input_buffer, "core:\n")) pres = &core_pressure;
		else if(!strcmp(input_buffer, "atm:\n"))  pres = &atm_pressure;
		if(pres!=NULL){
			getline(&input_buffer, &buffer_size, input_file);
			pressure = strtok(input_buffer, " \t\n");
			while(pressure != NULL){
				printf("\t%s", pressure);
				if(     !strcmp(pressure, "rad"))
					pres->push_back(rad_gas);
				else if(!strcmp(pressure, "ideal"))
					pres->push_back(ideal);
				else if(!strcmp(pressure, "coul"))
					 pres->push_back(coul);
				else if(!strcmp(pressure, "deg_zero"))
					pres->push_back(deg_zero);
				else if(!strcmp(pressure, "deg_partial"))
					pres->push_back(deg_partial);
				else if(!strcmp(pressure, "deg_finite"))
					pres->push_back(deg_finite);
				else if(!strcmp(pressure, "deg_trap"))
					pres->push_back(deg_trap);
				//else if(!strcmp(pressure, "deg_large"))
				//	pres->push_back(pressure_deg_large);
				else printf("ERROR: partial pressure term unrecognized!\n");
				pressure = strtok(NULL, " \t\n");
			}
			printf("\n");
		}
		pres = NULL;
	}
}

void ZAMS::initFromPolytrope(){
	printf("Preparing starting values from polytrope model\n");
	Mstar = Msolar*MSOLAR;
	Rstar = Rsolar*RSOLAR;
	
	Xtot.mu_e();
	
	Polytrope *testStar = new Polytrope(Mstar, Rstar, 1.5, Ntot);
	
	Rstar = testStar->Radius();
	Rsolar = Rstar/RSOLAR;
	Lstar = 4.*m_pi*boltzmann_sigma*pow(Rstar,2)*pow(Teff,4);
	Lsolar = Lstar/LSOLAR;
	
	//make initial guesses for core pressure based on this model
	// we must OVER-estimate, to avoid radiation pressure dominating
	Ystart0[dens] = 1.5*testStar->rho(0);
	Ystart0[pres] = 1.5*testStar->P(0);
	//guesses of surface -- not very good, don't use
	YstartS[dens] = testStar->rho(Ntot-2);
	YstartS[pres] = testStar->P(  Ntot-2);
	
	printf("Radius in CM: %le\n", Rstar);
	printf("Radius in RE: %le\n", Rsolar);
	
	printf("Core guess: %le \t %le\n", Ystart0[pres], Ystart0[dens]);
	printf("Surf guess: %le \t %le\n", YstartS[pres], YstartS[dens]);
	
	return;
}


//**************************************************************************************/
//  The Simple WD Grid
//    we will use a simple core+envelope grid, similar to that used by Hansen & Kawaler
//        in their ZAMS.f program
//    the core is considered to comprise the inner 85% of the star's mass
//        within the core, we assume a constant step size
//    the atmosphere uses a decreasing step size
//        the step size is modelled after a geometric series, and slowly approaches 0.9999
//**************************************************************************************/
void ZAMS::setupGrid(double Qcore, int Ncenter){
	//in the core, just use constant step size
	double dq = Qcore/(Ncenter-1);
	double logQcore = log(Qcore);
		
	//in the core, use constant step size
	double q = 0.0; //the center
	int n=0;
	for( ; n<Ncenter; n++){
		logQ[n] = log(q);
		q = q + dq;
	}
	
	//in the atmosphere, use a decreasing step size based on geometric series
	double dq1 = 1. - dq/(1.-Qcore);
	// fixed-point iteration to ensure self-consistency of dx1
	for(int a=0;a<20;a++) dq1 = 1.0 - dq/(1.-Qcore)*(1.-pow(dq1, (Ntot-Ncenter-1)/2));
	//now form the grid, starting at end of core
	for( ; n<=Ntot-1; n++){
		logQ[n] = log(q);
		if(n%2==0) dq = dq*dq1;
		q = q + dq;
	}
	
	// FOR SOME REASON, THIS HAS SLIGHTLY WRONG LIMIT
	for(n=0; n<Ntot; n++){
		logQ[n] = logQ[n]-logQ[Ntot-1];
	}
	indexFit = Ncenter;
}

//**************************************************************************************/
//  The ZAMS Core
//    the opacity is taken to be some constant small value, entirely due to conduction.
//**************************************************************************************/
double ZAMS::calculateCore(const double x[numv], int Nmax){
	//prepare variables for RK4
	double dlogQ;
	StellarVar logYC;
	double nabla, eps;
	double rholast = 10.0;
	
	//the array K contains the corrections to logY
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	//the first step must be handled specially
	int X = firstCoreStep(x, rholast, Nmax);
	//now begin integrations
	for(; X<Nmax; X++){
		dlogQ = logQ[X+1]-logQ[X];
		logYC = logY[X];
		for(int a = 0; a<4; a++){
			eps   = energyProduction(logYC, Xtot);
			nabla = energyTransport( logYC, Xtot);
			//now from these, calculate next shift
			//use radius as independent variable
			K[a] = dlogYdlogM(logYC, nabla, eps)*dlogQ;
			//calculate "corrected" positions using previous shift vectors
			logYC = logY[X] + K[a]*B[a];
			logYC[dens] = equationOfState(logYC, Xtot, rholast);
		}
		//update the structure variables
		logY[X+1] = logY[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;		
		//update the density from EOS
		logY[X+1][dens] = equationOfState(logY[X+1], Xtot,rholast);
		//save the derivatives
		eps   = energyProduction(logY[X+1], Xtot);
		nabla = energyTransport( logY[X+1], Xtot);
		dlogY[X+1] = dlogYdlogR(logY[X+1], nabla, eps);
		dlogY[X][dens] = (logY[X+1][dens]-logY[X][dens])/(logY[X+1][radi]-logY[X][radi]);
	}
	return rholast;
}

//**************************************************************************************/
//  First Step of Core
//    this method performs the first step in the integration from the core
//    this must be done using non-log variables, beause log variables are singular
//    this step must be integrated in radius, because all equations in mass are singular
//**************************************************************************************/
int ZAMS::firstCoreStep(const double x[numv], double& rholast, int Nmax){
	//prepare variables
	double P0 = (x[0]), T0 = (x[1]);

	//variables to use in RK4
	double dQ;
	StellarVar YC, logYC;
	double nabla, eps;
	//the arrays K and L contain the corrections to y and z to be applied
	// when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	//set initial conditions
	rholast = Ystart0[dens]; //sqrt(P0*Pscale); //an initial guess
	StellarVar Ycenter(rholast/Dscale, 0.0, P0, 0.0, T0, 0.0);
	StellarVar logYcenter = log(Ycenter);
	logYcenter[dens] = equationOfState(logYcenter, Xtot, rholast);
	Ycenter[dens] = exp(logYcenter[dens]);
	Ystart0[dens] = rholast;
	
	//begin RK4 at center
	int start = 1;	//number of points to compute in this way
	StellarVar Y[start+1];		//the array of points computed this way
	//begin at the center
	Y[0]     = Ycenter;			//the central values
	logY[0]  = logYcenter;		//the central log-values
	dlogY[0] = StellarVar();	//the central log-derivative is undefined
	//setting the initial step size
	dQ = exp((log(3.) + logQ[1] - logY[0][dens])/3.);//
	//now begin integrations
	for(int X = 0; X<start; X++){
		YC = Y[X];
		for(int a = 0; a<4; a++){
			eps   = energyProduction(log(YC), Xtot);
			nabla = energyTransport( log(YC), Xtot);
			K[a] = dYdR(YC, nabla, eps)*dQ;
			//at very center, account for r=0
			if(YC[radi]==0.0) K[a][pres] = K[a][temp] = 0.0;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[X] + K[a]*B[a];
			YC[dens] = exp(equationOfState(log(YC), Xtot, rholast));
		}
		//update the structure variables
		Y[X+1] = Y[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//calculate the log versions
		logY[X+1] = log(Y[X+1]);
		logY[X+1][mass] = logQ[X+1];
		logY[X+1][dens] = equationOfState(logY[X+1], Xtot, rholast);
		Y[X+1][dens] = exp(logY[X+1][dens]);
		//calculte derivatives at this location
		eps   = energyProduction(logY[X+1], Xtot);
		nabla = energyTransport( logY[X+1], Xtot);
		dlogY[X+1] = dlogYdlogR(logY[X+1], nabla, eps);
		dlogY[X+1][dens] = (logY[X+1][dens]-logY[X][dens])/(logY[X+1][radi]-logY[X][radi]);
	}
	return start;
}

//**************************************************************************************/
//  The ZAMS Helium Atmosphere
//    the EOS incudes electrons, ions, Coulomb corrections, partial ionization, and radiation
//    either radiative or convective heat transport allowed to take place
//    the composition is a mixture of Helium, Carbon, and Hydrogen
//	Must integrate from surface (photosphere) to the core
//**************************************************************************************/
void ZAMS::calculateAtmosphere(const double x[numv]){	
	//prepare variables
	double dlogQ;
	StellarVar logYC;
	double nabla, eps;
	double rholast = 1.0;
	
	//the array K contains the corrections to logY to be applied
	//    when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	
	//set our initial conditions
	int X = firstAtmosphereStep(x, rholast);
	for(; X > Ncore ; X--){
		dlogQ = logQ[X-1] - logQ[X];
		logYC = logY[X];		
		for(int a = 0; a<4; a++){
			//find energy production, heat transport
			eps   = energyProduction(logYC, Xtot);
			nabla = energyTransport( logYC, Xtot);
			//now from these, calculate next shift
			//use radius as independent variable
			K[a] = dlogYdlogM(logYC, nabla, eps)*dlogQ;
			//calculate "corrected" positions using previous shift vectors
			logYC = logY[X] + K[a]*B[a];
			logYC[dens] = equationOfState(logYC, Xtot, rholast);
		}
		//update the chemical composition and the density
		logY[X-1] = logY[X] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//update the chemical abundances and calculate the density from EOS
		logY[X-1][dens] = equationOfState(logY[X-1], Xtot, rholast);
		//calculate the derivatives at this location
		eps   = energyProduction(logY[X-1], Xtot);
		nabla = energyTransport( logY[X-1], Xtot);
		dlogY[X-1] = dlogYdlogR(logY[X-1], nabla, eps);
		dlogY[X-1][dens] = (logY[X-1][dens]-logY[X][dens])/(logY[X-1][radi]-logY[X][radi]);
	}
	return;
}

//**************************************************************************************/
//  First Step of Atmosphere
//    this method performs the first step in the integration from the surface
//    this must be done using non-log variables, beause log variables are singular
//**************************************************************************************/
int ZAMS::firstAtmosphereStep(const double x[numv], double& rholast){
	//prepare variables
	double rs = (1.0 + x[2]);
	double ms = (1.0);
	double LS = pow(1.0 + x[2],2);
	
	//set our initial conditions
	//set the surface values of logY
	//
	rholast = 1.e-9; //an initial guess (in CGS)
	StellarVar Ysurf(rholast/Dscale, rs,  1.0, ms, Teff/Tscale, LS);
	StellarVar ly = log(Ysurf)+logYscale;
	double kp = opacity(ly, Xtot);
	double gs = G_CGS*Mstar*pow(Rstar*rs,-2);
	double A = Xtot.mean_A()/(Teff*N_Avogadro*boltzmann_k);
	double a =2.*gs/3.;
	double b = radiation_a*pow(Teff,4)/6.;
	double f1 = atm_pressure(rholast, Teff, Xtot)-a/kp-b, f2=1.0e2;
	double x1=rholast, x2 = A*(a/kp - b);
//	printf("rho guess: %le %le %le\n", x1, x2, f1);
	while(fabs(x1-x2)/x1 > 1e-10){
		x1 = x2;
		x2 = 1.01*x1;
		ly[dens] = log(x2);
		kp = opacity(ly,Xtot);
		f2 = atm_pressure(x2,Teff,Xtot) - a/kp-b;
		x2 = x1 - f1*(x2-x1)/(f2-f1);
		ly[dens] = log(x2);
		kp = opacity(ly,Xtot);
		f1 = atm_pressure(x2,Teff,Xtot) - a/kp-b;
//		printf("rho guess: %le %le %le\n", x1, x2, f1);
	}	
	rholast = x2;
	ly[dens] = log(rholast);
	//kp = opacity(ly, Xtot);//*/
	double PS = atm_pressure(rholast, Teff, Xtot);
	
	//oldie but goodie
	PS = (1.e10);			//an initial guess
	rholast = sqrt(PS);     //an initial guess
	
	/*rholast = A*(a/kp - b);
	//for(int c=0; c<50; c++){
	while(fabs(exp(ly[dens])-rholast)/rholast >1e-4){
		printf("rho guess: %le %le %le\n", exp(ly[dens]), rholast, fabs(exp(ly[dens])-rholast)/rholast);
		ly[dens] = log(rholast);
		kp = opacity(ly, Xtot);
		rholast = A*(a/kp - b);
	}
	double PS = a/kp + b;//*/
	
	//set surface variables
	Ysurf = StellarVar( rholast/Dscale, rs, PS/Pscale, ms, Teff/Tscale, LS );
	logY[Ntot-1] = log(Ysurf);
	logY[Ntot-1][dens] = equationOfState(logY[Ntot-1], Xtot, rholast);
	Ysurf[dens] = exp(logY[Ntot-1][dens]);
	//calculate surface derivative
	double eps   = energyProduction(log(Ysurf), Xtot);
	double nabla = energyTransport( log(Ysurf), Xtot);
	dlogY[Ntot-1] = dlogYdlogR(logY[Ntot-1], nabla, eps);

	//set up RK4 integration
	double dQ;
	StellarVar YC, logYC;
	//the array K contains the corrections to logY to be applied
	//    when generating the next step of the correction
	StellarVar K[4];
	//Butcher tableau for RK4
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};

	int start = Ntot-2;
	StellarVar Y[(Ntot-start)];
	Y[0] = Ysurf;
	for(int X = Ntot-1; X > start; X--){
		YC = Y[Ntot-1-X];
		dQ = exp(logQ[X-1]) - exp(logQ[X]);
		for(int a = 0; a<4; a++){
			eps   = energyProduction(log(YC), Xtot);
			nabla = energyTransport( log(YC), Xtot);
			//now from these, calculate next shift
			K[a] = dYdM(YC, nabla, eps)*dQ;
			//calculate "corrected" positions using previous shift vectors
			YC = Y[Ntot-1-X] + K[a]*B[a];
			YC[dens] = exp(equationOfState(log(YC), Xtot, rholast));
			if(isnan(YC[pres])) printf("%d %d NAN pres\n", X, a);
			if(isnan(YC[dens])) printf("%d %d NAN dens\n", X, a);
			if(isnan(YC[mass])) printf("%d %d NAN mass\n", X, a);
			if(isnan(YC[temp])) printf("%d %d NAN temp\n", X, a);
			if(isnan(YC[lumi])) printf("%d %d NAN lumi\n", X, a);
		}
		//update the chemical composition and the density
		Y[Ntot-X] = Y[Ntot-X-1] + K[0]/6.0 + K[1]/3.0 + K[2]/3.0 + K[3]/6.0;
		//calculate the log versions
		logY[X-1] = log(Y[Ntot-X]);
		logY[X-1][dens] = equationOfState(logY[X-1], Xtot, rholast);
		Y[Ntot-X][dens] = exp(logY[X-1][dens]);
		//calculate the derivatives at this location
		eps   = energyProduction(logY[X-1], Xtot);
		nabla = energyTransport( logY[X-1], Xtot);
		dlogY[X-1] = dlogYdlogR(logY[X-1], nabla, eps);
		dlogY[X-1][dens] = (logY[X-1][dens]-logY[X][dens])/(logY[X-1][radi]-logY[X][radi]);
	}
	return start;
}


//**************************************************************************************/
//  JOIN AT CENTER -- Double-Shooting Algorithm 
//    Here we begin integration from core and from surface, joining at central location
//    Will return the difference in two integrations at the core
//    Will also return Jacobian matrix, how differences change with dependent parameters
//**************************************************************************************/
void ZAMS::joinAtCenter(double x[numv], double f[numv], double& F){
	StellarVar Y1, Y2;
	//calculate first using the current best guesses
	calculateCore(x,Ncore);
	Y1 = exp(logY[Ncore]);
	calculateAtmosphere(x);
	Y2 = exp(logY[Ncore]);
	Abundance Xm1 = massFraction();
	
	f[0] = (Y1[pres]-Y2[pres]);
	f[1] = (Y1[temp]-Y2[temp]);
	f[2] = (Y1[radi]-Y2[radi]);
//	f[3] = (Y1[lumi]-Y2[lumi]);
	
	F = 0.0;
	for(int i=0; i<numv; i++) F += 0.5*( f[i]*f[i] );

	return;
}

// differential equations
StellarVar ZAMS::dYdR(const StellarVar& y, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] = 0.0;
	deriv[pres] =-y[mass]*y[dens]*pow(y[radi],-2);
	deriv[temp] =-y[mass]*y[dens]*pow(y[radi],-2)*y[temp]/y[pres]*nabla;
	deriv[mass] = y[dens]*pow(y[radi],2);
	deriv[lumi] = y[dens]*pow(y[radi],2)*epsilon;
	deriv[radi] = 1.0;
	return deriv;
}

// logarithmic derivatives of variables
StellarVar ZAMS::dlogYdlogR(  const StellarVar& logy, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] =  0.0;
	deriv[pres] = -exp(logy[mass] + logy[dens] - logy[pres] - logy[radi]);
	deriv[temp] = -exp(logy[mass] + logy[dens] - logy[pres] - logy[radi])*nabla;
	deriv[mass] =  exp(logy[dens] - logy[mass]          + 3.0*logy[radi]);
	deriv[lumi] =  exp(logy[dens] - logy[lumi]          + 3.0*logy[radi])*epsilon;
	deriv[radi] =  exp(0.0);
	return deriv;
}

// differential equations using mass
StellarVar ZAMS::dYdM(const StellarVar& y, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] = 0.0;
	deriv[pres] =-y[mass]*pow(y[radi],-4);
	deriv[temp] =-y[mass]*pow(y[radi],-4)*y[temp]/y[pres]*nabla;
	deriv[radi] = pow(y[radi],-2)/y[dens];
	deriv[lumi] = epsilon;
	deriv[mass] = 1.0;
	return deriv;
}

// logarithmic derivatives of variables
StellarVar ZAMS::dlogYdlogM(  const StellarVar& logy, const double& nabla, const double epsilon){
	StellarVar deriv;
	deriv[dens] =  0.0;
	deriv[pres] = -exp(2.*logy[mass] - logy[pres] - 4.*logy[radi]);
	deriv[temp] = -exp(2.*logy[mass] - logy[pres] - 4.*logy[radi])*nabla;
	deriv[radi] =  exp(   logy[mass] - logy[dens] - 3.*logy[radi]);
	deriv[lumi] =  exp(   logy[mass] - logy[lumi])*epsilon;
	deriv[mass] =  exp(0.0);
	return deriv;
}


//**************************************************************************************/
//  Constitutive Physics
//**************************************************************************************/

// EPSILON  the source of luminosity
// Borrows from 
//		ZAMS.f by C.J. Hansen & S.D. Kawaler -- http://astro.if.ufrgs.br/evol/evolve/hansen/
//  	StatStar by Bradley W. Carroll & Dale A. Ostlie,  2007.
double ZAMS::energyProduction(const StellarVar& logy, const Abundance& X){
	//in our setting, dm/dr = rho*r^2,     with BCs 0 and 1
	//but also,       dL/dr = rho*r^2*eps, with BCs 0 and 1
	const double fpp = 1, f3a = 1;
	double T   = Tscale*exp(logy[temp]);
	double rho = Dscale*exp(logy[dens]);
	const double A_pp = 0.241, A_CNO = 8.67e20, A_He = 50.9;
    double lT6 = logy[temp] - 6.*log(10.) + logYscale[temp];
    double lT8 = logy[temp] - 8.*log(10.) + logYscale[temp];
    double lrho= logy[dens] + logYscale[dens];
	//PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
    double psipp = 1 + 1.412e8*(1/X.H1 - 1)*exp(-49.98*exp(-lT6/3.));
    double Cpp = 1 + 0.0123*exp(lT6/3.) + 0.0109*exp(2.*lT6/3.) + 0.000938*exp(lT6);
    double eps_pp = A_pp*X.H1*X.H1*fpp*psipp*Cpp*exp(-2.*lT6/3.+lrho)*exp(-33.80*exp(-lT6/3.));
    //CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
    double XCNO = (X.C12+X.O16);
    double CCNO = 1 + 0.0027*exp(lT6/3.) - 0.00778*exp(2.*lT6/3.) - 0.000149*exp(lT6);
    double eps_CNO = A_CNO*X.H1*XCNO*CCNO*exp(-2.*lT6/3.+lrho)*exp(-152.28*exp(-lT6/3.));
    //Helium burning (Kippenhahn and Weigert, Eq. 18.67)
    double eps_He = A_He*pow(X.He4, 3)*exp(-3.*lT8 + 2.*lrho)*f3a*exp(-44.027/exp(lT8));
    //return epsilon; combined energy generation rate
    return (eps_pp + eps_CNO + eps_He)*Mstar/Lstar;
    
    //this part directly from ZAMS.FOR
	//if(exp(ly[temp]) < 14.0*log(10)){
	//	return 1.0e-20;
	//}
	//T9=DEXP(YGO(2))/1.0D9
	//T913=T9**(1.0D0/3.0D0)
	//T923=T913*T913
	//EPSPP=2.4D4*RHO*X*X*DEXP(-3.380/T913)/T923
	//EPSCNO=4.4D25*RHO*X*Z*DEXP(-15.228D0/T913)/T923
	//EPS=EPSPP+EPSCNO
}

// NABLA    the heat transportation, calculting from radiative/conductive transport
double ZAMS::energyTransport( const StellarVar& logy, const Abundance& X)
{
	//rescale variables to CGS
	StellarVar ly = logy + logYscale;
	double kappa = opacity(ly, X);
	double D_coeff = 3./(16.*m_pi*radiation_a*C_CGS*G_CGS);
	double delta = D_coeff*kappa*exp(ly[lumi]+ly[pres]-ly[mass]-4.*ly[temp]);
	if(isnan(delta)) delta = D_coeff*kappa*exp(ly[pres]-4.*ly[temp])*Lstar/Mstar;
	
	double delta_ad = getEOS(exp(ly), X)->nabla_ad(exp(ly[dens]), exp(ly[temp]), X);
		
	if(delta_ad < delta) {
		return delta_ad;
	}
	else if(delta < delta_ad) {
		return delta;
	}
	else return delta;
}

// EQUATION OF STATE  relates pressure and density
EOS* ZAMS::getEOS(const StellarVar &y, const Abundance& X){
//	return &core_pressure;
	double Patm  =  atm_pressure(y[dens], y[temp], X);
	double Pcore = core_pressure(y[dens], y[temp], X);
	if      (Pcore <  Patm) return &atm_pressure;
	else if (Pcore >= Patm) return &core_pressure;
	else return &core_pressure; //*/
}

double ZAMS::equationOfState(const StellarVar& logy, const Abundance& chem, double& rho_last){
	StellarVar y = exp(logy+logYscale);
	//find the radiation pressure
	double Prad, rho=0.0;
	Prad = radiation_a/3.*pow(y[temp],4);
	if(Prad > y[pres]) {
		printf("\nERROR: RADIATION PRESSURE TOO HIGH!\n");
	//	y[temp] = 0.0;
	}
	
	rho = getEOS(y, chem)->invert(rho_last, y[pres], y[temp], chem);
	
	if(rho<0.0) rho = -rho;
	if(isnan(rho)) printf("RHO IS NAN!\n");
	if(isnan(log(rho))){
		printf("%le %le %le %le %le %le\n", exp(logy[radi]), y[pres], y[temp], rho, rho/Dscale, log(rho/Dscale));
	}
	
	//save the density	
	rho_last = rho;
	//return the log density
	return log(rho/Dscale);
}

//**************************************************************************************/
//  Chemical Composition
//**************************************************************************************/
Abundance ZAMS::findAbundance(
	const double logm,
	const double logr,
	Abundance& dX
)
{
	dX = Abundance(0.,0.,0.,0.);
	Xtot.e = chemical::h;		
	return Xtot;
}

Abundance ZAMS::massFraction(){
	return Xtot;
}

//**************************************************************************************/
//  ACCESSORS
//    Here we define functions to access radius, pressure, etc.
//**************************************************************************************/
double ZAMS::rad(int X){
	return Rstar*exp(logY[X][radi]);
}
double ZAMS::rho(int X){
	return Dscale*exp(logY[X][dens]);
}
double ZAMS::drhodr(int X){
	 // dlogD/dlogr = dD/dr*r/D --> dD/dr = D/r * dlogD/dlogr
	return Dscale/Rstar*exp(logY[X][dens]-logY[X][radi])*dlogY[X][dens];
}
double ZAMS::P(int X){
	return Pscale*exp(logY[X][pres]);
}
double ZAMS::dPdr(int X){
	return Pscale/Rstar*exp(logY[X][pres]-logY[X][radi])*dlogY[X][pres];
}
double ZAMS::Phi(int X){
	//zeroed to join exterior solution at surface, where Phi->0 at infty
	return 0.0;
}
double ZAMS::dPhidr(int X){
	return G_CGS*Mstar*pow(Rstar,-2)*exp(logY[X][mass] - 2.*logY[X][radi]);
}
double ZAMS::mr(int X){
	return Mstar*exp(logY[X][mass]);
}

double ZAMS::Schwarzschild_A(int X, double GamPert){
	if(GamPert==0.0) return -brunt_vaisala[X]*exp(2.*logY[X][radi]-logY[X][mass])/Rstar;
	else        	 return (dlogY[X][dens] - dlogY[X][pres]/GamPert)*exp(-logY[X][radi])/Rstar;
}

double ZAMS::getAstar(int X, double GamPert){
	if(GamPert==0.0) return brunt_vaisala[X]*exp(3.*logY[X][radi]-logY[X][mass]);
	else        	 return (dlogY[X][pres]/GamPert - dlogY[X][dens]);
}

//the Ledoux part of the Schwarzschild discriminant
double ZAMS::Ledoux(int X, double GamPert){
	if(GamPert==0.0) return ledoux[X];
	else             return 0.0;
}

double ZAMS::getU(int X){
	if(X==0) return 3.0;
	else     return dlogY[X][mass];	//this is already dlog(m)/dlog(r)
}

double ZAMS::getVg(int X, double GamPert){
	if(GamPert==0.0) return -dlogY[X][pres]/Gamma1(X);//we already have dlog(P)/dlog(r)
	else			 return -dlogY[X][pres]/GamPert;
}

double ZAMS::getC(int X){
	if(X==0) return 3.*exp(-logY[0][dens]);
	else     return exp(3.*logY[X][radi] - logY[X][mass]);
}

double ZAMS::Gamma1(int X){
	return adiabatic_1[X];
}

double ZAMS::sound_speed2(int X, double GamPert){
	if(GamPert == 0.0) return Gamma1(X) * exp(logY[X][pres]-logY[X][dens])*Pscale/Dscale;
	else               return GamPert   * exp(logY[X][pres]-logY[X][dens])*Pscale/Dscale;
}

void ZAMS::populateBruntVaisala(){
	adiabatic_1 = new double[Ntot];
	nabla = new double[Ntot];
	nabla_ad = new double[Ntot];
	brunt_vaisala = new double[Ntot];
	ledoux = new double[Ntot];
	kappa = new double[Ntot];
	
	StellarVar YY = exp(logY[0]+logYscale);
	StellarVar ly;
	double g, chiT, chir, nbad;
	double P, rho, Pcore, Patm;
	double kappa_rad, kappa_cond, kappa0, D_coeff;
	EOS* myEOS;

	for(int X=0; X<Ntot;X++){
		ly = logY[X] + logYscale;
		YY = exp(ly);
			
		kappa[X] = opacity(ly, Xtot);
		D_coeff = 3./(16.*m_pi*radiation_a*C_CGS*G_CGS);
		nabla[X] = D_coeff*kappa[X]*exp(ly[pres]+ly[lumi]-ly[mass]-4.*ly[temp]);
		myEOS = getEOS(YY, Xtot);

		adiabatic_1[X] = myEOS->Gamma1(  YY[dens],YY[temp],Xtot);
		nabla_ad[X]    = myEOS->nabla_ad(YY[dens],YY[temp],Xtot);
		ledoux[X]      = myEOS->Ledoux(  YY[dens],YY[temp],Xtot,dXelem);
		chiT           = myEOS->chiT(    YY[dens],YY[temp],Xtot);
		chir           = myEOS->chiRho(  YY[dens],YY[temp],Xtot);
		if(nabla_ad[X] < nabla[X]) nabla[X] = nabla_ad[X];
		brunt_vaisala[X] = exp(2.*logY[X][mass]-4.*logY[X][radi]+logY[X][dens]-logY[X][pres])
				*chiT/chir*(nabla_ad[X]-nabla[X]+ledoux[X]);
	}
}

double ZAMS::Radius(){return Rstar*exp(logY[Ntot-1][radi]);}	//total radius
double ZAMS::Mass(){  return Mstar*exp(logY[Ntot-1][mass]);}	//total mass
double ZAMS::Gee(){return G_CGS;}
//in Newtonian, light speed is infinity...
double ZAMS::light_speed2(){return C_CGS*C_CGS;}

// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R
void ZAMS::setupCenter(){
	nc = 1./(Gamma1(1)-1.);
	double rho0 = Dscale*exp(logY[0][dens]);
	double P0   = Pscale*exp(logY[0][pres]);
	double term = m_pi*rho0*rho0/P0/(nc+1.);
	ac[0] = 1.0; ac[1] = 0.0;
	ac[2] = -2./3.*term; ac[3] = 0.0;
	ac[4] = 2.*nc/15.*pow(term,2); ac[5] = 0.0;
	ac[6] = -4.*nc*(8.*nc-5.)/(945.)*pow(term,3); ac[7] = 0.0;
}

void ZAMS::getAstarCenter(double *AC, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(1) : g);
	double AV = (nc*gam1)/(nc+1.)-1.;
	getVgSurface(AC, maxPow, g);
	for(int k=0; k<=maxPow; k++) AC[k] *= AV;
}

void ZAMS::getVgCenter(double *Vc, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(1) : g);
	if(maxPow>=0) Vc[0] = 0.0;
	if(maxPow>=2) Vc[1] = -2.*(1.+nc)*ac[2]/gam1;
	if(maxPow>=4) Vc[2] = (2.*ac[2]*ac[2]-4.*ac[4])*(1.+nc)/gam1;
	if(maxPow> 4) maxPow = 4;
}

void ZAMS::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = 3.0;
	if(maxPow>=2) Uc[1] = 6.*nc/5.*ac[2];
	if(maxPow>=4) Uc[2] = 12./7.*nc*ac[4]+(24./175.*nc-6./7.)*nc*ac[2]*ac[2];
	if(maxPow> 4) maxPow = 4;	
}

void ZAMS::getC1Center(double *cc, int& maxPow){
	double term = 1./(4.*m_pi*exp(logY[0][dens])*Dscale);
	if(maxPow>=0) cc[0] = getC(0.0);
	if(maxPow>=2) cc[1] = -9.*nc/5.*ac[2]*term;
	if(maxPow>=4) cc[2] = ((9./14.+153.*nc/350.)*nc*ac[2]*ac[2] - 9.*nc/7.*ac[4])*term;
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
void ZAMS::setupSurface(){}

void ZAMS::getAstarSurface(double *As, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(Ntot-2) : g);
	double AV = (1.5*gam1)/(1.5+1.) - 1.;
	int O=1;
	getVgSurface(As, maxPow, g);
	for(int k=-1; k<=maxPow; k++){
		As[O+k] = AV*As[O+k];
	}
}

void ZAMS::getVgSurface(double *Vs, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(Ntot-2) : g);
	int O=1;
	if(maxPow>= -1) Vs[O-1] = (1.5+1.)/gam1;
	if(maxPow>=  0) Vs[O  ] = 0.0;
	if(maxPow>=  1) Vs[O+1] = 0.0;
	if(maxPow>=  2) Vs[O+2] = -315./16./gam1;
	if(maxPow>=  3) Vs[O+3] = 1575./32./gam1;
	//if more  terms than this requested, cap number of terms
	if(maxPow > 3) maxPow = O+3;
}

void ZAMS::getUSurface(double *Us, int& maxPow){
// coefficients of U must extend up to order maxPow	
	for(int a=0; a<=maxPow; a++) Us[a] = 0.0;
}

void ZAMS::getC1Surface(double *cs, int& maxPow){
// coefficients of c1 are only needed up to order maxPow-1
	if(maxPow>=0) cs[0]  =  1.;
	if(maxPow>=1) cs[1]  = -3.;
	if(maxPow>=2) cs[2]  =  3.;
	if(maxPow>=3) cs[3]  = -1.;
	if(maxPow>=4) cs[4]  =  0.; //does not actually appear in equations
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;	
}

// **************************  WRITE STAR  ***************************************
// will print the following graphs of the star, using gnuplot
// 1) graph of T,m,rho vs. radius
// 2) graph of chemical composition vs. inverse mass fraction
// 3) graph of chemical composition vs. radius
// 4) graph of Brunt-Vaisala and Lamb frequency vs. logP (to compare to Goldreich & Wu)
// optional argument c[] is calculation directory where files should be written
// if no argument, files written to ./out/[name], where name is shown in constructor
void ZAMS::writeStar(char *c){
	printf("PRINTING STAR\n");
	//create names for files to be opened
	char filename[256];
	char rootname[256];
	char txtname[256];
	char outname[256];
	if(c==NULL)	sprintf(filename, "./out/%s", name);
	else{
		sprintf(filename, "./%s/star", c);
	}
	sprintf(txtname, "%s/star.txt", filename);
	sprintf(outname, "%s/star.png", filename);

	char title[256]; graph_title(title);

	//prepare the output directory, making sure it exists
	FILE *fp;
	if(!(fp = fopen(txtname, "w")) ){
		if(c==NULL){
			system("mkdir ./out");
			char command[256]; sprintf(command, "mkdir ./out/%s", name);
			system(command);
			fp = fopen(txtname, "w");
		}
		else {
			char command[256]; sprintf(command, "mkdir ./%s", c);
			system(command);
			sprintf(command, "mkdir %s", filename);
			system(command);
			if(!(fp = fopen(filename, "w"))) printf("big trouble, boss\n");		
		}
	}
	
	//print results to text file
	// radius rho pressure gravity
	double logT = logY[0][temp];
	double logL = logY[Ntot-1][lumi];
	double logM = logY[Ntot-1][mass];
	double logR = logY[Ntot-1][radi];
	for(int X=0; X< Ntot; X++){
		fprintf(fp, "%d\t", X);
		fprintf(fp, "%0.16le\t", (1.-exp(logY[X][mass])));	//col 2
		fprintf(fp, "%0.16le\t", exp(logY[X][dens])*Dscale);	//col 3
		fprintf(fp, "%0.16le\t", exp(logY[X][radi]-logR));		//col 4
		fprintf(fp, "%0.16le\t", exp(logY[X][pres])*Pscale);	//col 5
		fprintf(fp, "%0.16le\t", exp(logY[X][mass]));			//col 6
		fprintf(fp, "%0.16le\t", exp(logY[X][temp]-logT));		//col 7
		fprintf(fp, "%0.16le\t", exp(logY[X][lumi]));			//col 8
		fprintf(fp, "%0.16le\n", exp(logQ[X]));					//col 9
		//fflush(fp);
	}
	fclose(fp);	
	//plot file in png in gnuplot, and open png
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Profile for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'rho, P'\n");
	fprintf(gnuplot, "set logscale y 10\n set ytics 10 nomirror\nset format y '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale x 10\n set format x '%%L'\n");
	fprintf(gnuplot, "set yrange [1e0:1e24]\n");
//	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	fprintf(gnuplot, "set y2range [-0.01:1.01]\nset y2tics 0.2 nomirror\n");
	fprintf(gnuplot, "set y2label 'normalized temperature, radius'\n");
	fprintf(gnuplot, "plot ");
	fprintf(gnuplot, "     '%s' u 2:7 w l t 'T' axes x1y2", txtname);
	fprintf(gnuplot, ",    '%s' u 2:3 w l t 'rho'", txtname);
	fprintf(gnuplot, ",    '%s' u 2:5 w l t 'P'", txtname);
	fprintf(gnuplot, ",    '%s' u 2:4 w l t 'r' axes x1y2", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
		
	//make a graph of T, nabla, opacity
	sprintf(txtname, "%s/opacity.txt", filename);
	sprintf(outname, "%s/opacity.png", filename);
	fp  = fopen(txtname, "w");
	Abundance X;
	StellarVar ly;
	for(int x=0; x < Ntot; x++){
		ly = logY[x];
		fprintf(fp, "%0.24le\t%0.16le", (1.-exp(ly[mass]-logM)), exp(ly[temp]-logT));	//col1, col2
		//get the layer, for using in energyTransport()
		fprintf(fp, "\t%0.16le\t%0.16le", nabla[x], nabla_ad[x]); //col3, col4
		//form the approximation for opacity
		ly = logY[x] + logYscale;
		fprintf(fp, "\t%0.16le", kappa[x]); //col5
		fprintf(fp, "\t%0.16le", adiabatic_1[x]); //col6
		fprintf(fp, "\n");
	}
	fclose(fp);	
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'opacity, temperature for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log(1-m/M)'\n");
	fprintf(gnuplot, "set y2label 'temperature'\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set logscale y 10\n");
//	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	fprintf(gnuplot, "set yrange [1e-5:1e6]\n");
	fprintf(gnuplot, "set y2range [0:1]\n");
	fprintf(gnuplot, "set ytics 10 nomirror\n");
	fprintf(gnuplot, "set y2tics 0.1 nomirror\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'T' axes x1y2", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'nabla'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'nabla_{ad}'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'opacity'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:6 w l t 'Gamma1' lc 10", txtname);
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	
	//print the Brunt-Vaisala frequency
	sprintf(txtname, "%s/BruntVaisala.txt", filename);
	sprintf(outname, "%s/BruntVaisala.png", filename);
	fp  = fopen(txtname, "w");
	double N2 = -1.0;
	for(int X=1; X< length(); X++){
		N2 = brunt_vaisala[X]*G_CGS*Mstar*pow(Rstar,-3);
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			(1.-exp(logY[X][mass]-logM)),
			N2,
			2.*sound_speed2(X,0.)*exp(-2.*logY[X][radi])*pow(Rstar,-2));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	//fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Brunt-Vaisala and Lamb for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 & log_{10} L_1^2 (Hz^2)\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
//	fprintf(gnuplot, "set xrange [1e-8:1e0]\n");
	fprintf(gnuplot, "set yrange [1e-8:1e2]\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'N^2'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'L_1^2'", txtname);
	//fprintf(gnuplot, ",    '%s' u 1:(-$2) w l t '-N^2'", txtname);
	fprintf(gnuplot, "\n");
	
	//print the pulsation coeffcients frequency
	sprintf(txtname, "%s/coefficients.txt", filename);
	sprintf(outname, "%s/coefficients.png", filename);
	fp  = fopen(txtname, "w");
	int maxpow=4;
	double A,U,V,C, read[4], x0 = exp(logY[1][radi]);
	getAstarCenter(read, maxpow, 0);
	A = read[0] + read[1]*x0*x0 + read[2]*pow(x0,4);
	getVgCenter(read, maxpow, 0);
	V = read[0] + read[1]*x0*x0 + read[2]*pow(x0,4);
	getC1Center(read, maxpow);
	C = read[0] + read[1]*x0*x0 + read[2]*pow(x0,4);
	getUCenter(read, maxpow);
	U = read[0] + read[1]*x0*x0 + read[2]*pow(x0,4);
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		0.0, A, U, V, C);
	for(int X=1; X< length()-1; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			exp(logQ[X]),
			-exp(logY[X][radi])*Schwarzschild_A(X, 0.)*Rstar,
			getU(X),
			getVg(X, 0.),
			getC(X)
		);
	}
	double reads[5], t1 = 1.-exp(logY[Ntot-2][radi]);
	maxpow=4;
	getAstarSurface(reads, maxpow, 0);
	A = reads[0]/t1 + reads[1] + reads[2]*t1;
	getVgSurface(reads, maxpow, 0);
	V = reads[0]/t1 + reads[1] + reads[2]*t1;
	getC1Surface(reads, maxpow);
	C = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	getUSurface(reads, maxpow);
	U = reads[0] + reads[1]*t1 + reads[2]*t1*t1;
	fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
		1.0, A, U, V, C);
	fclose(fp);	
	//plot file in png in gnuplot
	//gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname);
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title);
	//fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'A*, U, V_g, c_1'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'A*'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'U'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'V_g'", txtname);
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'c_1'", txtname);
	fprintf(gnuplot, "\n");
	
	
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}

//calculte a scaled sum square residual from backsubstitution into physical equations
double ZAMS::SSR(){	
	double checkEuler;
	double checkPoiss;
	int len = Ntot;
			
	//sum up errors in equations
	checkEuler = 0.0;
	checkPoiss = 0.0;
	double d2Phi = 0.0;
	double e1, e2, n1, n2;
	FILE *fp = fopen("SSR.txt", "w");
	for(int X=4; X<len-4; X++){
		//Euler equation
		e1 = fabs(dPdr(X) + rho(X)*dPhidr(X) );
		n1 = fabs(dPdr(X)) + fabs(rho(X)*dPhidr(X));
		//Poisson equation
		//calculate numerical derivatives
		d2Phi = Gee()*dlogY[X][mass]*pow(rad(X),-3)*mr(X) - 2.*Gee()*mr(X)*pow(rad(X),-3);
		e2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)  -      2.0*dPhidr(X)  -      d2Phi*rad(X) );
		n2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)) + fabs(2.0*dPhidr(X)) + fabs(d2Phi*rad(X) );
		//add absolute error
		e1 = e1/n1;
		e2 = e2/n2;
		fprintf(fp, "%d\t%le\t%le\n", X, e1, e2);
		checkEuler += e1*e1;
		checkPoiss += e2*e2;
	}
	fclose(fp);
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %d\n", length());
	fprintf(gnuplot, "set output '%s'\n", "SSR.png");
	char title[256]; graph_title(title);
	fprintf(gnuplot, "set title 'error for %s'\n", title);
	fprintf(gnuplot, "set xlabel 'gird point'\n");
	fprintf(gnuplot, "set ylabel 'error'\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set arrow 1 from %d, graph 0 to %d, graph 1 lc rgb 'red' nohead\n", Ncore, Ncore);
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'Euler error'", "SSR.txt");
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'Poisson error'", "SSR.txt");
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
	return sqrt((checkPoiss+checkEuler)/double(2*len));
}

#endif