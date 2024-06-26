//**************************************************************************************
//							MESA Wrapper
// MESA.cpp
//		This is a simple wrapper, intended to work with data generated using MESA
//	
// Reece Boston, Sep 02, 2023
//**************************************************************************************

//*****
 // WARNING: must have at least one interpolated step at a midpoint, 
 // or the RK4 method being used to find modes will NOT work.
 //*****


#ifndef MESACLASS
#define MESACLASS

#include "MESA.h"

MESA::MESA(const char* filename, std::size_t L) : len(L) {
	printf("Beginning read-in of MESA data.\n");
	FILE *infile;
	if( !(infile = fopen(filename, "r")) ) {
		printf("No such file %s\n", filename);
		printf("In directory: ");fflush(stdout);
		system("pwd");
		printf("Available files:\n");fflush(stdout);
		system("ls");
		exit(1);
	}  
	fscanf(infile, "  %lu    %le     %le     %le     %*d\n", &Ntot, &Mstar, &Rstar, &Lstar);
	printf("\tMass = %0.2lg\tRadis = %0.2lg\tLuminosity = %0.2lg\n", Mstar, Rstar, Lstar);
	printf("Number of grid points in MESA calculation %lu\n", Ntot);
	name = strmakef("MESA.M%1.2g.R%1.2g.L%1.2g", Mstar, Rstar, Lstar);
			
	//read through file once to get the max values of R,M
	std::size_t firstdata = ftell(infile);
	for(std::size_t k=0; k<Ntot-1; k++){
		fscanf(infile, "%*[^\n]\n");
	}
	fscanf(infile, "%*d %lg %lg %lg %*lg %lg %*lg %lg %*[^\n]", &Rtot, &Mtot, &Ltot, ts, &dels);

	//set scale quantities
	Dscale = Mtot*pow(Rtot,-3)/(4.*m_pi);
	Pscale = G_CGS*pow(Mtot,2)*pow(Rtot,-4)/(4.*m_pi);
	Gscale = G_CGS*Mtot*pow(Rtot,-2);
	
	//create temporary arrays to hold values read from file
	double *rt = new double[Ntot]; //radius
	double *dt = new double[Ntot]; //density
	double *pt = new double[Ntot]; //pressure
	double *mt = new double[Ntot]; //mass
	double *gt = new double[Ntot]; //gravitational field
	double *Nt = new double[Ntot]; //Brunt-Vaisala frequency
	double *Gt = new double[Ntot]; //Gamma1
		
	//read from the file
	fseek(infile, firstdata, SEEK_SET);
	for(std::size_t k=0; k<Ntot; k++){
		fscanf(infile, " %*d %lg  %lg  %*lg", &rt[k], &mt[k]);
		fscanf(infile, " %lg %*lg %lg  %*lg", &pt[k], &dt[k]);
		fscanf(infile, " %lg %lg           ", &Nt[k], &Gt[k]);
		fscanf(infile, " %*[^\n]\n");
		//dis-dimensionalize the variables
		rt[k] /= Rtot;
		mt[k] /= Mtot;
		pt[k] /= Pscale;
		dt[k] /= Dscale;
		if(rt[k]>0)	gt[k] = mt[k]*pow(rt[k],-2);
		else gt[k] = 0.0;
		Nt[k] /= (G_CGS*Mtot*pow(Rtot,-3));
	}
	double Tscale = 2.*proton.mass_CGS/boltzmann_k*G_CGS*Mtot/Rtot;
	ts[0] /= Tscale;
	
	printf("Done reading in data!\n");
	fclose(infile);
	
	//divide intervals up into binary fractions
	//must use binary fractions to mesh with RK4 method requiring half-points
	double n = log2(len-1) - log2(Ntot-1);
	if(n < 1) n = 1; //must at least divide each interval in half
	len = pow(2,int(n))*(Ntot-1) + 1;
	subgrid = pow(2,int(n));
	printf("number of grid points to be used %lu\n", len);
			
	radi = new double[len];
	std::size_t kk=0;
	for(std::size_t k=0; k<Ntot-1; k++){
		double dr = (rt[k+1]-rt[k])/subgrid;
		for(std::size_t ki=0; ki<subgrid; ki++,kk++) radi[kk] = rt[k]+double(ki)*dr;
	}
	radi[len-1] = rt[Ntot-1];
		
	//now we need to make the stellar variables for pulsation
	double *At = new double[Ntot]; //Astar
	double *Ut = new double[Ntot]; //U
	double *Ct = new double[Ntot]; //c1
	double *Vt = new double[Ntot]; //Vg
	
	At[0] = 0.0; Ut[0] = 3.0; Vt[0] = 0.0; Ct[0] = 3./dt[0];
	for(std::size_t k=1; k<Ntot-1; k++){
		At[k] = Nt[k]*rt[k]/gt[k];
		Ct[k] = pow(rt[k],3)/mt[k];
		Ut[k] = dt[k]*pow(rt[k],3)/mt[k];
		Vt[k] = gt[k]*rt[k]*dt[k]/pt[k]; //intentionally exclude Gamma1
	}
	At[Ntot-1] = 0.0; Ut[Ntot-1] = 0.0; Vt[Ntot-1] = 0.0; Ct[Ntot-1] = 1.0;
	
	//now splne fit all quantitites -- remove arrays to conserve space
	dens = new Splinor(rt, dt, Ntot);
	double rhof = dt[0]/2.0;
	delete[] dt;
	pres = new Splinor(rt, pt, Ntot);
	delete[] pt;
	mass = new Splinor(rt, mt, Ntot);
	delete[] mt;
	grav = new Splinor(rt, gt, Ntot);
	delete[] gt;
	BVfq = new Splinor(rt, Nt, Ntot);
	delete[] Nt;
	Gam1 = new Splinor(rt, Gt, Ntot);
	delete[] Gt;
	//pulsation variables
	aSpline = new Splinor(rt, At, Ntot);
	delete[] At;
	uSpline = new Splinor(rt, Ut, Ntot);
	delete[] Ut;
	vSpline = new Splinor(rt, Vt, Ntot);
	delete[] Vt;
	cSpline = new Splinor(rt, Ct, Ntot);
	delete[] Ct;
	//remove the temporary radial array
	delete[] rt;
	
	indexFit = len/2;
	for(std::size_t X=1; X<len-1; X++){
		//set matching point as where density is half its central value
		if(dens->interp(radi[X-1])>rhof & dens->interp(radi[X+1])<=rhof)
			indexFit = X;
	}
	indexFit /= 2;
	
	//prepare boundary condition values
	setupCenter();
	setupSurface();
}

//destructor
MESA::~MESA(){
	delete[] radi;
	delete dens;
	delete pres;
	delete mass;
	delete grav;
	delete Gam1;
	delete BVfq;
	delete aSpline;
	delete uSpline;
	delete cSpline;
	delete vSpline;
}

double MESA::Radius(){return Rstar;} //total radius in appropriate units
double MESA::Mass()  {return Mstar;} //total mass   in appropriate units
double MESA::Gee()   {return G_CGS;} //Newton's constant in appropriate units

//Here we define functions to access radius, pressure, etc.
double MESA::rad(std::size_t X){
	if(X>=0 & X<len) return Rtot*radi[X];
	else {
		printf("\nradius out of range\n");	
		return nan("");
	}
}
double MESA::rho(std::size_t X){
	if(X>=0 & X<len) return Dscale*dens->interp(radi[X]);
	else {
		printf("\nrho out of range\n");	
		return nan("");
	}
}
double MESA::drhodr(std::size_t X){
	if(X==0) return 0.0;
	else if(X>0 & X<len) return Dscale/Rtot*dens->deriv(radi[X]);
	else {
		printf("\ndrhodr out of range\n");	
		return nan("");
	}
}
double MESA::P(std::size_t X){
	if(X>=0 & X<len) return Pscale*pres->interp(radi[X]);
	else {
		printf("\nP out of range\n");	
		return nan("");
	}
}
double MESA::dPdr(std::size_t X){
	if(X==0) return 0.0;
	else if(X>0 & X<len) return Pscale/Rtot*pres->deriv(radi[X]);
	else {
		printf("\ndPdr out of range\n");	
		return nan("");
	}
}
/** NEEDS TO BE IMPLEMENTED **/
double MESA::Phi(std::size_t X){
	return 0.0;
}
double MESA::dPhidr(std::size_t X){
	if(X>=0 & X<len) return Gscale*grav->interp(radi[X]);
	else {
		perror("\ng out of range\n");	
		return nan("");
	}
}
double MESA::mr(std::size_t X){
	if(X>=0 && X<len) return Mtot*mass->interp(radi[X]);
	else {
		perror("\nm out of range\n");	
		return nan("");
	}
}

double MESA::Schwarzschild_A(std::size_t X, double GamPert){
	/*NOTE:  A* = -r*A, but rad is not dimensionless -- so must divide by R here*/
	if(GamPert==0.0) return -aSpline->interp(radi[X])/radi[X]/Rtot;
	else             return dens->deriv(radi[X])/dens->interp(radi[X])/Rtot
						  - pres->deriv(radi[X])/pres->interp(radi[X])/Rtot/GamPert;
}

double MESA::getAstar(std::size_t X, double GamPert){
	if(GamPert==0.0) return aSpline->interp(radi[X]);
	else             return radi[X]*pres->deriv(radi[X])/pres->interp(radi[X])/GamPert
						  - radi[X]*dens->deriv(radi[X])/dens->interp(radi[X]);
}

double MESA::getU(std::size_t X){
	//U = 4pi rho r/g
	if(X==0) return 3.0;
	return uSpline->interp(radi[X]);
}

double MESA::getVg(std::size_t X, double GamPert){
	//Vg = gr rho/(Gamma1*P)
	if(GamPert==0.0) return vSpline->interp(radi[X])/Gam1->interp(radi[X]);
	else			 return vSpline->interp(radi[X])/GamPert;
}

double MESA::getC(std::size_t X){
	if(X==0) return 3./dens->interp(0.0);
	return cSpline->interp(radi[X]);
}


double MESA::Gamma1(std::size_t X){return Gam1->interp(radi[X]);}

double MESA::sound_speed2(std::size_t X, double GamPert){
	if(GamPert == 0.0) return Pscale/Dscale*pres->interp(radi[X])/dens->interp(radi[X])*Gam1->interp(radi[X]);
	else               return Pscale/Dscale*pres->interp(radi[X])/dens->interp(radi[X])*GamPert;
}


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R
//	We can use general expressions for c1, A*, Vg, U in terms of coefficients of rho, P, T,
//  The coefficients of rho, P, T are found by assuming a polytrope equation of state, P=P(rho)
//  Then coefficients for rho can be found in terms of the polytrope solution
//  This method allows for more flexibility, if we choose a different approximate EOS in center
void MESA::setupCenter(){
	nc = 1./(Gamma1(0)-1.);
	double rho0 = dens->interp(0.0);
	double P0 = pres->interp(0.0);
	double term = rho0*rho0/P0/(nc+1.)/4.;
	ac[0] = 1.0; ac[1] = 0.0;
	ac[2] = -2./3.*term; ac[3] = 0.0;
	ac[4] = 2.*nc/15.*pow(term,2); ac[5] = 0.0;
	ac[6] = -4.*nc*(8.*nc-5.)/(945.)*pow(term,3); ac[7] = 0.0;
	
	//terms x^0
	dc[0] = dens->interp(0.0);
	pc[0] = pres->interp(0.0);
	//terms x^2
	pc[1] = (nc+1.) * ac[2] * pc[0];
	dc[1] = nc      * ac[2] * dc[0];
	//terms x^4
	pc[2] = 0.5 * (nc+1.) * (2. * ac[4] +  nc     * ac[2] * ac[2]) * pc[0];
	dc[2] = 0.5 * nc      * (2. * ac[4] + (nc-1.) * ac[2] * ac[2]) * dc[0];
	
	//now the variables
	// c1
	c0[0] = 3./dc[0];
	c0[1] = -1.8*dc[1]/dc[0]/dc[0];
	c0[2] = 1.08*dc[1]*dc[1]*pow(dc[0],-3) - (9./7.)*dc[2]/dc[0]/dc[0];
	// Vg -- deliberately exclude factor 1/Gamma1
	V0[0] = 0.0;
	V0[1] = -2.*pc[1]/pc[0];
	V0[2] = 2.*pc[1]*pc[1]/pc[0]/pc[0] - 4.*pc[2]/pc[0];
	// U
	U0[0] = 3.0;
	U0[1] = 1.2*dc[1]/dc[0];
	U0[2] = -0.72*dc[1]*dc[1]/pc[0]/pc[0] + (12./7.)*dc[2]/dc[0];
	// A*
	double N20=BVfq->deriv(0.0);
	A0[0] = 0.0;
	A0[1] = N20*c0[0];
	A0[2] = N20*c0[1];	
}

void MESA::getAstarCenter(double *Ac, int& maxPow, double g){
	if(maxPow>=0) Ac[0] = A0[0];
	if(maxPow>=2) Ac[1] = A0[1];
	if(maxPow>=4) Ac[2] = A0[2];
	if(maxPow> 4) maxPow = 4;
}

void MESA::getVgCenter(double *Vc, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gamma1(1) : g);
	if(maxPow>=0) Vc[0] = V0[0]/gam1;
	if(maxPow>=2) Vc[1] = V0[1]/gam1;
	if(maxPow>=4) Vc[2] = V0[2]/gam1;
	if(maxPow> 4) maxPow = 4;
}

void MESA::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = U0[0];
	if(maxPow>=2) Uc[1] = U0[1];
	if(maxPow>=4) Uc[2] = U0[2];
	if(maxPow> 4) maxPow = 4;
}

void MESA::getC1Center(double *cc, int& maxPow){
	if(maxPow>=0) cc[0] = c0[0];
	if(maxPow>=2) cc[1] = c0[1];
	if(maxPow>=4) cc[2] = c0[2];
	//if more  terms than this requested, cap number of terms
	if(maxPow> 4) maxPow = 4;
}


// **************************  SURFACE BOUNDARY  ***************************************
// the following provide coefficients for surface expansions of A*, Vg, U, c1 in terms of t=1-r/R
// These terms are borrowed from polytrope, for n=1.5
//	Note that A*, Vg require a power -1
//	Note that up to order N, we only A*, Vg, c1 up to order N-1, but require U to order N
//	Here maxPow represents the maixmum order of expansions of y1,..,y4 in LAWE
//	If maxPow = 0, we need terms -1
//	If maxPow = 1, we need terms -1, 0
//	If maxPow = 2, we need terms -1, 0, 1
//	We still do not have very good values for this boundary
void MESA::setupSurface(){	
	double const Tscale = 2.*proton.mass_CGS/boltzmann_k*G_CGS*Mtot/Rtot;
	double const Cideal = N_Avogadro*boltzmann_k*Dscale*Tscale/Pscale;
	double const Crad   = radiation_a/3.*pow(Tscale,4)/Pscale;
	//the x^0 part
	ts[0] = ts[0]; // redundant, but for consistency
	ps[0] = pres->interp(1.0);
	ds[0] = dens->interp(1.0);
	// the x^1 part
	ts[1] = 0.0;
	ps[1] = 0.0;
	ds[1] = 0.0; //(ps[1]-4.*Crad*pow(ts[0],3)*ts[1]-Cideal*ts[1]*ds[0])/(Cideal*ts[0]);
	
	// the x^2 part
	double tempprod=1.0;
	int mp = 2;
	ts[2] =(2.*ps[0]*ts[0]*ds[0]-ps[1]*ts[0]*ds[0]+ps[0]*ts[1]*ds[0]-ps[0]*ts[0]*ds[0]*ds[0]+ps[0]*ts[0]*ds[1]);
	ts[2] = dels*ts[2]/(2.*ps[0]*ps[0]);
	ps[2] = 0.5*(2.*ds[0] - ds[0]*ds[0] + ds[1]);
	ds[2] = ps[2] - Crad*(4.*ts[0]*ts[1]*ts[1]*ts[1]+12.*ts[0]*ts[0]*ts[1]*ts[2]+ts[0]*ts[0]*ts[0]*ts[3]);
	tempprod=0.0;
	for(int a=0; a<=mp; a++){
		for(int b=0; b<=mp-a; b++){
			for(int c=0; c<=mp-a-b; c++){
				tempprod += ts[a]*ts[b]*ts[c]*ts[mp-a-b-c];
			}
		}
	}
	ds[mp] -= Crad*tempprod; 
	for(int a=0; a<mp; a++){
		ds[mp] -= Cideal*ts[mp-a]*ds[a];
	}
	ds[mp] /= (Cideal*ts[0]);
	
	//the x^3 part
	mp = 3;
	ts[3] = (2.*ps[1]*ts[1]*ds[1]-ps[2]*ts[1]*ds[1]+ps[1]*ts[2]*ds[1] + ps[1]*ts[1]*ds[2]);
	ts[3] = ts[3]*dels/(3.*ps[1]*ps[1]);
	ps[3] = ds[0] - ds[0]*ds[0]/3. + 2.*ds[1]/3. - 0.5*ds[0]*ds[1] + ds[2]/3.;
	ds[3] = ps[3];
	tempprod=0.0;
	for(int a=0; a<=mp; a++){
		for(int b=0; b<=mp-a; b++){
			for(int c=0; c<=mp-a-b; c++){
				tempprod += ts[a]*ts[b]*ts[c]*ts[mp-a-b-c];
			}
		}
	}
	ds[mp] -= Crad*tempprod;
	for(int a=0; a<mp; a++){
		ds[mp] -= Cideal*ts[mp-a]*ds[a];
	}
	ds[mp] /= (Cideal*ts[0]);
	
	// the x^4 part
	mp = 4;
	ts[4] = (dels/(24.*pow(ps[0],4) ))*
 	(-6.*pow(ps[1],3)*ts[0]*ds[0] + 6.*ps[0]*ps[1]*(2.*ps[2]*ts[0]*ds[0]
 		 + ps[1]*(ts[1]*ds[0] + ts[0]*(-(-2. + ds[0])*ds[0] + ds[1]))) + 
  		3.*ps[0]*ps[0]*(-2*ps[3]*ts[0]*ds[0] + 2.*ps[2]*((-ts[1] + ts[0]*(-2 + ds[0]))*ds[0] - ts[0]*ds[1])
  		+ ps[1]*(2*(-ts[2] + ts[1]*(-2 + ds[0]))*ds[0] - 2.*ts[1]*ds[1]
  		+ ts[0]*(ds[0]*(-6. + 2.*ds[0] + 3.*ds[1]) - 2.*(2.*ds[1] + ds[2]))))
  		+ pow(ps[0],3)*(6.*ts[3]*ds[0] + 6.*ts[2]*(-(-2. + ds[0])*ds[0] + ds[1])
  		+ 3.*ts[1]*(ds[0]*(6. - 2.*ds[0] - 3.*ds[1]) + 2.*(2.*ds[1] + ds[2])) 
  		+ ts[0]*(-8.*ds[0]*ds[0] - 3.*(-6. + ds[1])*ds[1] - 8.*ds[0]*(-3. + ds[1] + ds[2]) 
  		+ 6.*(2.*ds[2] + ds[3])))); //magical formula from mathematica
	ps[4] = (24.*ds[0] - 8.*ds[0]*ds[0] + 18.*ds[1] - 8.*ds[0]*ds[1] - 3.*ds[1]*ds[1] 
			+ 12.*ds[2] - 8.*ds[0]*ds[2] + 6.*ds[3])/24.;
	tempprod=0.0;
	for(int a=0; a<=mp; a++){
		for(int b=0; b<=mp-a; b++){
			for(int c=0; c<=mp-a-b; c++){
				tempprod += ts[a]*ts[b]*ts[c]*ts[mp-a-b-c];
			}
		}
	}
	ds[mp] -= Crad*tempprod;
	for(int a=0; a<mp; a++){
		ds[mp] -= Cideal*ts[mp-a]*ds[a];
	}
	ds[mp] /= (Cideal*ts[0]);
		
	//now set the structure variables
	// c
	c1[0] =  1.;
	c1[1] = -3.;
	c1[2] =  3. + 0.5*ds[1];
	c1[3] = -1. - 13.*ds[1]/6. + ds[2]/3.;
	c1[4] = 0.; // does not actually appear in equations
	// U
	U1[0] = 0.0;
	U1[1] = ds[1];
	U1[2] = ds[2] - 3.*ds[1];
	U1[3] = ds[3] - 3.*ds[2] + ds[1]*ds[1]/2. + 3.*ds[1];
	U1[4] = ds[4] - 3.*ds[3] + 5./6.*ds[1]*ds[2] + 3.*ds[2] - 13./6.*ds[1]*ds[1] - ds[1];
	// Vg
	int O = 1; // an anchor index, allows more natural expression of coefficients
	for( std::size_t i=2; i<=4; i++) ps[i] /= ps[1];
	// for( std::size_t i=2, i<=4; i++) ds[i] /= ds[1];
	V1[O-1] = 1.0;
	V1[O+0] = ps[2] - 1.0;
	V1[O+1] = 2.*ps[3] - ps[2]*(1.+ps[2]);
	V1[O+2] = 3.*ps[4] - 3.*ps[3]*ps[2] - 2.*ps[3] + ps[2]*(1.+ps[2]);
	V1[O+3] = 0.0; //4.*ps[5] - 4.*ps[4]*ps[2] - 3.*ps[4] - 2.*ps[3]*ps[3];
	// A*
	// actually only holding 
	double N21 = BVfq->interp(1.0);
	A1[O-1] = 0.0;
	A1[O+0] = N21*c1[0];
	A1[O+1] = N21*c1[1];
	A1[O+2] = N21*c1[2];
	A1[O+3] = N21*c1[3];
	for (std::size_t i=2; i<=4; i++) ps[i] *= ps[1];
	// for (std::size_t i=2; i<=4; i++) ds[i] *= ds[1];		
}

void MESA::getAstarSurface(double *As, int& maxPow, double g){
	maxPow = std::min(4, maxPow);
	for (std::size_t n = 0; n <= maxPow; n++) {
		As[n] = A1[n];
	}
}

void MESA::getVgSurface(double *Vs, int& maxPow, double g){
	double gam1 = (g==0.0 ? Gam1->interp(radi[len-1]) : g);
	maxPow = std::min(4, maxPow);
	for (std::size_t n = 0; n <= maxPow; n++) {
		Vs[n] = V1[n]/gam1;
	}
}

void MESA::getUSurface(double *Us, int& maxPow){
	// coefficients of U must extend up to order maxPow	
	maxPow = std::min(4, maxPow);
	for (std::size_t n = 0; n <= maxPow; n++) {
		Us[n] = U1[n];
	}
}

void MESA::getC1Surface(double *cs, int& maxPow){
	// coefficients of c1 are only needed up to order maxPow-1
	maxPow = std::min(4, maxPow);
	for (std::size_t n = 0; n <= maxPow; n++) {
		cs[n] = c1[n];
	}
}


void MESA::printBV(const char *const c, double const g){
	std::string filename = c, rootname, txtname, outname;
	std::string title = graph_title();
	
	//print the Brunt-Vaisala frequency
	txtname = filename + "/BruntVaisala.txt";
	outname = filename + "/BruntVaisala.png";
	FILE* fp  = fopen(txtname.c_str(), "w");
	double N2 = -1.0;
	double scaleN2 = G_CGS*Mtot*pow(Rtot,-3);		//put N^2 back into CGS units
	for(std::size_t X=1; X< length(); X++){
		N2 = scaleN2*BVfq->interp(radi[X]);
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\n",
			(1.-mass->interp(radi[X])),
			N2,
			2.*sound_speed2(X, 0.0)*pow(Rtot*radi[X],-2));//Lambd frequency
	}
	fclose(fp);
	//plot file in png in gnuplot
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 800,800\n");
	fprintf(gnuplot, "set samples %lu\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Brunt-Vaisala for %s'\n", title.c_str());
	fprintf(gnuplot, "set logscale x\n");
	fprintf(gnuplot, "set format x '10^{%%L}'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xlabel 'log_{10} (1-m/M)\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 (Hz^2)\n");
	fprintf(gnuplot, "set hidden3d noundefined\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'N^2'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'L_1^2'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(-$2) w l t '-N^2'", txtname.c_str());
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

void MESA::printCoefficients(const char *const c, double const g){

	std::string filename, rootname, txtname, outname;
	filename = addstring(c, "/wave_coefficient");
	system( ("mkdir -p " + filename).c_str() );
	
	//print the coefficients of the center and surface, for series analysis
	txtname = filename + "/center.txt";
	FILE *fp = fopen(txtname.c_str(), "w");
	double gam1 = Gam1->interp(radi[1]);
	fprintf(fp, "dens:\t%0.16le\t%0.16le\t%0.16le\n", dc[0],dc[1],dc[2]);
	fprintf(fp, "pres:\t%0.16le\t%0.16le\t%0.16le\n", pc[0],pc[1],pc[2]);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\t%0.16le\n", A0[0],A0[1],A0[2]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\t%0.16le\n", U0[0],U0[1],U0[2]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\t%0.16le\n", V0[0]/gam1,V0[1]/gam1,V0[2]/gam1);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\t%0.16le\n", c0[0],c0[1],c0[2]);
	fclose(fp);
	txtname = filename + "/surface.txt";
	fp = fopen(txtname.c_str(), "w");
	gam1 = Gam1->interp(radi[len-1]);
	fprintf(fp, "dens:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", ds[0],ds[1],ds[2],ds[3],ds[4]);
	fprintf(fp, "pres:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", ps[0],ps[1],ps[2],ps[3],ps[4]);
	fprintf(fp, "A*  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", A1[0],A1[1],A1[2],A1[3],A1[4]);
	fprintf(fp, "U   :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", U1[0],U1[1],U1[2],U1[3],U1[4]);
	fprintf(fp, "Vg  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", V1[0]/gam1,V1[1]/gam1,V1[2]/gam1,V1[3]/gam1,V1[4]/gam1);
	fprintf(fp, "c1  :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", c1[0],c1[1],c1[2],c1[3],c1[4]);
	fclose(fp);
	
	//print fits to those coefficients at center and surface
	std::size_t NC=15, NS=15;
	txtname = filename + "/centerfit.txt";
	fp = fopen(txtname.c_str(), "w");
	double x2 = 0.0;
	gam1 = Gam1->interp(radi[1]);
	fprintf(fp, "x           \trho         \trho_fit     \tP           \tP_fit       \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(std::size_t X=0; X<NC; X++){
		x2 = radi[X]*radi[X];
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			radi[X],
			dens->interp(radi[X]),
			dc[0]+dc[1]*x2+dc[2]*x2*x2,
			pres->interp(radi[X]),
			pc[0]+pc[1]*x2+pc[2]*x2*x2,
			-radi[X]*aSpline->interp(radi[X]),
			A0[0] + A0[1]*x2 + A0[2]*x2*x2,
			uSpline->interp(radi[X]),
			U0[0] + U0[1]*x2 + U0[2]*x2*x2,
			vSpline->interp(radi[X])/gam1,
			V0[0]/gam1 + V0[1]/gam1*x2 + V0[2]/gam1*x2*x2,
			cSpline->interp(radi[X]),
			c0[0] + c0[1]*x2 + c0[2]*x2*x2
		);
	}
	fclose(fp);
	txtname = filename + "/surfacefit.txt";
	fp = fopen(txtname.c_str(), "w");
	double t = 1.-radi[len-2];
	gam1 = Gam1->interp(radi[len-1]);
	fprintf(fp, "x           \trho         \trho_fit     \tP           \tP_fit       \tA*          \tA*_fit      \tU           \tU_fit       \tVg          \tVg_fit      \tc1          \tc1_fit\n");
	for(std::size_t X=len-1; X>=len-NS-1; X--){
		t = 1. - radi[X];
		fprintf(fp, "%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\t%+12le\n",
			radi[X],
			dens->interp(radi[X]),
			ds[0]+ds[1]*t+ds[2]*t*t+ds[3]*t*t*t+ds[4]*t*t*t*t,
			pres->interp(radi[X]),
			ps[0]+ps[1]*t+ps[2]*t*t+ps[3]*t*t*t+ps[4]*t*t*t*t,
			-radi[X]*aSpline->interp(radi[X]),
			(A1[0])/t + (A1[1]) + (A1[2])*t + (A1[3])*t*t + (A1[4])*t*t*t,
			uSpline->interp(radi[X]),
			U1[0] + U1[1]*t + U1[2]*t*t + U1[3]*t*t*t + U1[4]*t*t*t*t,
			vSpline->interp(radi[X])/gam1,
			(V1[0]/gam1)/t + (V1[1]/gam1) + (V1[2]/gam1)*t + (V1[3]/gam1)*t*t + (V1[4]/gam1)*t*t*t,
			cSpline->interp(radi[X]),
			c1[0] + c1[1]*t + c1[2]*t*t + c1[3]*t*t*t + c1[4]*t*t*t*t
		);
	}
	fclose(fp);
	
	//print the pulsation coeffcients frequency
	txtname = filename + "/coefficients.txt";
	outname = filename + "/coefficients.png";
	fp  = fopen(txtname.c_str(), "w");
	for(std::size_t X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			radi[X],
			-radi[X]*aSpline->interp(radi[X]),
			uSpline->interp(radi[X]),
			vSpline->interp(radi[X])/Gam1->interp(radi[X]),
			cSpline->interp(radi[X]));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %lu\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", graph_title().c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'A*, U, V_g, c_1'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'V_g'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'c_1'", txtname.c_str());
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
	//fits
	gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	txtname = filename + "/centerfit.txt";
	outname = filename + "/centerfit.png";
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Central Fitting by Power Series for %s'\n", graph_title().c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [0:%le]\n", radi[NC]);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/$2) w lp t 'rho'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/$4) w lp t 'P'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/$6) w lp t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/$8) w lp t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($10-$11)/$10) w lp t 'Vg'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($12-$13)/$12) w lp t 'c1'", txtname.c_str());
	fprintf(gnuplot, "\n");
	txtname = filename + "/surfacefit.txt";
	txtname = filename + "/surfacefit.png";
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Surface Fitting by Power Series for %s'\n", graph_title().c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [%le:1]\n", radi[len-NS-1]);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'rho'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'P'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($10-$11)/abs($10)) w lp t 'Vg'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($12-$13)/abs($12)) w lp t 'c1'", txtname.c_str());
	fprintf(gnuplot, "\n");
	
	//now leave gnuplot
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}


//method to print pertinent values of star to .txt, and plot them in gnuplot
void MESA::writeStar(const char *const c){
	//create names for files to be opened
	std::string filename, rootname, txtname, outname;
	if(c==NULL)	filename = strmakef("./out/%s", name.c_str());
	else filename = strmakef("./%s/star/", c);
	txtname = strmakef("%s/%s.txt", filename.c_str(), name.c_str());
	outname = strmakef("%s/%s.png", filename.c_str(), name.c_str());

	FILE *fp;
	if(!(fp = fopen(txtname.c_str(), "w")) ){
		system( ("mkdir -p " + filename).c_str() );
		if(!(fp = fopen(txtname.c_str(), "w"))) 
			printf("big trouble, boss\n");		
	}
	//print results to text file
	// radius rho pressure gravity
	double 	irc=1./dens->interp(radi[0]), 
			ipc=1./pres->interp(radi[0]),
			ig =1./grav->interp(radi[len-1]);
	for(std::size_t X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			radi[X],
			dens->interp(radi[X])*irc,  -dens->deriv(radi[X])*irc*Rtot,
			pres->interp(radi[X])*ipc,  -pres->deriv(radi[X])*ipc*Rtot,
			mass->interp(radi[X]), grav->interp(radi[X])*ig);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	std::string title = graph_title();
	fprintf(gnuplot, "set title 'Profile for %s'\n", title.c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'rho/rho_c, P/P_c, m/M, g/g_S'\n");
	fprintf(gnuplot, "set xrange [0:1]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname.c_str());
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname.c_str());
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname.c_str());
	fprintf(gnuplot, ", '%s' u 1:7 w l t 'g'", txtname.c_str());
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
		
	printBV(filename.c_str());
	printCoefficients(filename.c_str());	
}

double MESA::SSR(){	
	//sum up errors in equations
	double checkEuler = 0.0;
	double checkPoiss = 0.0;
	double d2Phi = 0.0;
	double e1, e2, n1, n2;
	FILE *fp = fopen("SSR.txt", "w");
	for(std::size_t X=2; X<len-4; X++){
		//Euler equation
		e1 = fabs(pres->deriv(radi[X])  +      dens->interp(radi[X])*grav->interp(radi[X]));
		n1 = fabs(pres->deriv(radi[X])) + fabs(dens->interp(radi[X])*grav->interp(radi[X]));
		//Poisson equation
		d2Phi = grav->deriv(radi[X]); //use derivative of spline interpolator
		e2 = fabs( dens->interp(radi[X])*radi[X]  -      2.0*grav->interp(radi[X])       - d2Phi*radi[X] );
		n2 = fabs( dens->interp(radi[X])*radi[X]) + fabs(2.0*grav->interp(radi[X])) + fabs(d2Phi*radi[X] );
		//add absolute errors
		e1 = e1/n1;
		e2 = e2/n2;
		fprintf(fp, "%le\t%le\t%le\n", radi[X], e1, e2);
		checkEuler += e1*e1;
		checkPoiss += e2*e2;
	}
	fclose(fp);
	return sqrt((checkPoiss+checkEuler)/double(2*len-5));
}

#endif