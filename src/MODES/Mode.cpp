//**************************************************************************************
//							The MODE Object
//	Mode.cpp
//		Equation agnostic 
//			-- all information specific to physics supplied by the driver
// 		Capable of Nonradial, Cowling modes by supplying different drivers
//  Reece Boston Mar 24, 2022
//  NOTE: If you make ANY changes to this file, you MUST make clean and recompile
//**************************************************************************************

#ifndef MODECLASS
#define MODECLASS

#include "Mode.h"
#include "../../lib/rootfind.h"

//This is the basic setup routine common to all constructors
// preparing all arrays for integration and matching
//The individual constructors differ only in how they initialize frequency
template <std::size_t  numvar>
void Mode<numvar>::basic_setup(){
	wronskian = [this](double x) -> double {return this->calculateWronskian(x);};
	cee2 = star->light_speed2();
	Gee  = star->Gee();

	//begin by setting all values to 1 -- will be rescaled for final solution
	for(int i=0;i<numvar;i++){
		yCenter[i]  = 1.0;
		ySurface[i] = 1.0;
	}
	//we set up stellar grid as twice the mode grid to avoid interpolation at half points
	len_star = star->length();
	len      = driver->length();
	Gamma1   = driver->Gamma1();
	
	//determine the index where fit for inward and outward integrations will occur
	xfit = star->indexFit;
	//set the arrays that define the perturbation
	rad = new double[len];
	y = new double*[len];
	for(int i=0;i<len;i++) y[i] = new double[numvar];
	//define values of radius for the mode
	for(std::size_t X=0; X<len; X++){
		rad[X] = driver->rad(X);
	}
	//define the conversion factor from dimensionless frequency to angular f in rad/s
	// f^2 = (omeg2freq) * w^2
	omeg2freq = Gee * star->Mass() * pow(star->Radius(), -3);
	
	//set up boundary matrix and index array
	double **yy = new double*[numvar];
	int *ind = new int[numvar];
	for(int yi=0; yi<numvar; yi++) yy[yi] = new double[numvar];
	driver->getBoundaryMatrix(num_var, yy, ind);
	for(int yi=0; yi<numvar; yi++) for(int yj=0; yj<numvar; yj++) boundaryMatrix[yi][yj]=yy[yi][yj];
	for(int yi=0; yi<numvar; yi++) indexOrder[yi] = ind[yi];
	delete[] yy;
	delete[] ind;
}

//This constructor is given an initial value of frequency to use
//This initial guess is used as a starting point in search for frequency
template <std::size_t numvar> 
Mode<numvar>::Mode(double omg2, int l, int m, ModeDriver *drv)
	: l(l), m(m), omega2(omg2), driver(drv), star(drv->star)
{	
	basic_setup();
	//now we  converge to solution
	converge();
	this->k = verifyMode();
}

//This constructor guesses initial frequency based on a desired k,l mode
//The initial guess is taken from the analytic solutions for n=0 polytrope
template <std::size_t numvar> 
Mode<numvar>::Mode(int K, int L, int M, ModeDriver *drv)
	: l(L), m(M), k(K), driver(drv), star(drv->star)
{	
	basic_setup();
	
	double ell = double(l);
	omega2 = 0.0;
	//for p-modes
	if(k>=0){
		//for initial guess, may as well use homogeneous frequency
		//see Cox 17.7, homogeneous model
		//homogeneous model frequencies usually larger, so reduce size of n used
		// there is no magic to this guess -- it just seems to work okay
		double ks = double(k); // - (2*l-4));		
		// 2Dn = -4 + Gamma1*[ n*(2l+2n+5) + 2l + 3]
		double gam1 = (Gamma1 == 0.0 ? star->Gamma1(0) : Gamma1);
		double Dn = gam1*( ks*(ks+ell+0.5) ) - 2.0;
		//omega^2 = Dn + sqrt( Dn^2 + l(l+1) )
		omega2 = Dn + sqrt( Dn*Dn + ell*(ell+1.0) );
		//from dimensional considerations, sigma^2 = (GM/R^3) *omega^2
	}
	//for g-modes
	else {
		//find base sigma2 when n=-1 (use above formula, simplify)
		double sig1 = -2.0 + sqrt(4. + ell*(ell+1.));
		//convert to freq
		sig1 = sqrt( sig1*nug ); //fq = fg*sigma
		//empirical formula freq_n = freq_-1 * exp(-0.4 sqrt(-n-1)) for n<0
		//determined for polytropes -- works less well with realistic WD stars
		sig1 = sig1*exp(-0.4*sqrt(-k));
		//now convert freq back to sigma2
		omega2 = (sig1*sig1)/nug; //sigma^2
	}

	//now we  converge to solution
	//by matching outward and inward solutions in the center
	converge();
	this->k = verifyMode();
}

//This constructor is given a range within which the desired frequency is known to exist
//The range must bound exactly one eigenfrequency
template <std::size_t numvar> 
Mode<numvar>::Mode(double omeg2lo, double omeg2hi, int l, int m, ModeDriver *drv)
	: l(l), m(m), driver(drv), star(drv->star)
{	
	basic_setup();

	//BEGIN BISECTION HUNT
	converged = false;
	//if omegas are in wrong order, swap 'em!
	if (omeg2lo > omeg2hi){
		double tmp = omeg2hi;
		omeg2hi = omeg2lo;
		omeg2lo = tmp;
	}	
	//the brackets given are often themselves zeros -- inch them closer to avoid this
	double dw = (omeg2hi - omeg2lo) / 100.0;
	double w2min = omeg2lo + dw;
	double w2max = omeg2hi - dw;

	//now find the initial bracketing values of the Wronskian
	double w2 = 0.5*(w2min + w2max);
	double W = rootfind::bisection_search(this->wronskian, w2, w2min, w2max);
	
	//if the Wronskian still isn't good, just find values as normal
	if(W > 1e-10) {
		omega2=w2;
		converge();
	}
	else {
		//we now know our value of omega2
		omega2=w2;
		linearMatch(omega2, yCenter, ySurface);
		converged = true;
	}
	this->k = verifyMode();
}

//destructor
template <std::size_t numvar> 
Mode<numvar>::~Mode(){
		//free space used by y1,...,y4 from stack
		delete[] rad;
		for(int i=0; i<len; i++)
			delete[] y[i];
}

template <std::size_t numvar>
void Mode<numvar>::converge(){
	converged = false;
	convergeBisect(0.0);
	linearMatch(omega2, yCenter, ySurface);
	converged = true;
}

//Find frequency using a bisection search, based on Wronskian, up to tolerance tol
template <std::size_t numvar>
void Mode<numvar>::convergeBisect(double tol){
	double w = omega2, wmin = 0.0, wmax, W;
	rootfind::bisection_find_brackets_move(this->wronskian, w, wmin, wmax);
	W = rootfind::bisection_search(this->wronskian, w, wmin, wmax);
	omega2 = w;
}	

//Newton convergence on omega2, up to tolerance tol, max number of steps term
//if term = 0, then no maximum number of steps (until integer overflow, anyway)
template <std::size_t numvar> 
void Mode<numvar>::convergeNewton(double tol, int term){
	double w = omega2, W;
	W = rootfind::newton_search(this->wronskian, w, w, tol, (std::size_t)term);
	omega2 = w;
}

//integrates from both edges toward center and returns Wronskian
//this approach based on Christensen-Dalsgaard, Astrophys.SpaceSci.316:113-120,2008
template <std::size_t numvar> 
double Mode<numvar>::calculateWronskian(double w2){		
	//prepare a matrix
	double DY[numvar][numvar];
	//store values of outward solution at fitting point in rows
	for(int i=0;i<numvar/2;i++){
		RK4out(xfit, w2, boundaryMatrix[i]);
		for(int j=0;j<numvar;j++) DY[i][j] = y[xfit][j];
	}
	//store values of inward solution at fitting point in rows
	for(int i=numvar/2; i<numvar; i++){
		RK4in( xfit, w2, boundaryMatrix[i]);
		for(int j=0;j<numvar;j++) DY[i][j] = y[xfit][j];
	}
	//the Wronskian is the determinant of this matrix
	return matrix::determinant(DY);
}

//linearly match inward and outward solutions
template <std::size_t numvar> 
void Mode<numvar>::linearMatch(double w2, double y0[numvar], double ys[numvar]){
	double DY[numvar][numvar];
	for(int i=0;i<numvar/2;i++){
		RK4out(xfit, w2, boundaryMatrix[i]  );
		for(int j=0;j<numvar;j++) DY[i][j] = y[xfit][j];
	}
	for(int i=numvar/2; i<numvar; i++){
		RK4in( xfit, w2, boundaryMatrix[i]);
		for(int j=0;j<numvar;j++) DY[i][j] = y[xfit][j];
	}
		
	//ALROGITHM TO FIND COEFFICIENTS
	double A[numvar][numvar];
	for(int i=0;i<numvar;i++){
		for(int j=0; j<numvar/2; j++)      A[i][j] = DY[j][i]; //outward solutions +
		for(int j=numvar/2; j<numvar; j++) A[i][j] =-DY[j][i]; // inward solutions -
	}
	
	double aa[numvar] = {0.0};
	double bb[numvar] = {0.0};
	matrix::invertMatrix(A, bb, aa);
	
	//for the basis BCs we chose, this will be the properly scaled physical solution
	//if we change the BCs, we must change these results to match
	for(int a=0; a<numvar/2; a++){
		y0[indexOrder[a]] *= aa[a];
	}
	for(int a=numvar/2; a<numvar-1; a++){
		ys[indexOrder[a]] *= aa[a];
	}
	ys[indexOrder[0]] = 1.0;
	
	RK4out(xfit, w2, y0);
	RK4in( xfit, w2, ys);
}

//integrates outward from interior to xmax
template<std::size_t numvar>
void Mode<numvar>::RK4out( std::size_t xmax, double w2, double y0[numvar]){
	std::size_t start = driver->CentralBC(y, y0, w2, l);
	double dx;
	for(std::size_t x = start; x < xmax; x++){
		dx = rad[x+1] - rad[x];
		RK4step(x, w2, dx, y[x], y[x+1]);
	}
}

//integrates inward from surface toward xmin
template<std::size_t numvar>
void Mode<numvar>::RK4in( std::size_t xmin, double w2, double ys[numvar]){
	std::size_t start = driver->SurfaceBC(y, ys, w2, l);
	double dx;
	for(std::size_t x = start; x > xmin; x--){
		dx = rad[x-1] - rad[x];
		RK4step(x, w2, dx, y[x], y[x-1]);
	}
}

template <std::size_t numvar>
void Mode<numvar>::RK4step(std::size_t x, double w2, double dx, double yin[numvar], double yout[numvar]){
	int sign = (std::signbit(dx) ? -1 : +1);
	double XC = rad[x], YC[numvar];
	double K[4][numvar];
	static double coeff[numvar][numvar];
	static const double B[4] = {0.5, 0.5, 1.0, 0.0};
	static const double b[4] = {6.0, 3.0, 3.0, 6.0};
	static const int    d[4] = {0, 1, 1, 2};

	for(std::size_t i=0; i<numvar; i++) YC[i] = yin[i];
	for(std::size_t a=0; a<4; a++){
		driver->getCoeff( &coeff[0][0], x, sign * d[a], w2, l);
		// calcuate shift for next step
		for(std::size_t i=0; i<numvar; i++){
			K[a][i] = 0.0;
			for(std::size_t j=0; j<numvar; j++) K[a][i] += dx * coeff[i][j]*YC[j]/XC;
		}
		//evaluate next correction from shift
		XC = rad[x] + B[a] * dx;
		for(std::size_t i=0; i<numvar; i++) YC[i] = yin[i] + B[a]*K[a][i];
	}
	for(std::size_t i=0; i<numvar; i++){
		yout[i] = yin[i];
		for(std::size_t a=0; a<4; a++) yout[i] += K[a][i]/b[a];
	}
}

//this method serves to verify that the n is indeed the desired mode number
// according to Scuflaire and Osaki classification scheme, counting zeros of xi
template <std::size_t numvar> 
int Mode<numvar>::verifyMode(){
	//trace through solution in (xi, chi) plane, parameterized by index
	int quad=0, quadP, N=0;
	for(std::size_t x=1; x<len-2; x++){
		quadP = quad;
		//determine quadrant of (xi, chi)
		quad = (y[x][0]>=0 ? (y[x][1]>0? 1 : 2 ) : (y[x][1]>=0? 4 : 3));
		if(quadP==quad) continue;
		//if solution rotates clockwise, count as negative modes (g modes)
		if(quadP == 1 & quad == 2) N--;
		if(quadP == 3 & quad == 4) N--;
		//if solution rotates counter-clockwise, positive modes (p modes)
		if(quadP == 2 & quad == 1) N++;
		if(quadP == 4 & quad == 3) N++;	
	}
	return N;
}

//print out the mode information and plot it on gnuplot
template <std::size_t numvar> 
void Mode<numvar>::writeMode(const char *const c){
	//create names for files to be opened
	std::string filename, modename = strmakef("/mode_%d.%d", l, k);
	if(c==NULL)	filename = "./out/" + star->name + "/mode";
	else filename = addstring("./", c) + "/modes";
	std::string rootname = filename + modename;
	//save data to folder to avoid clutter - make sure folder exists
	std::string txtname = rootname + ".txt";
	std::string outname = rootname + ".png";
	FILE *fp;
	if(!(fp = fopen(txtname.c_str(), "w")) ){
		system( ("mkdir -p "+filename).c_str() );
		fp = fopen(txtname.c_str(), "w");
	}
	double R = rad[len-1];
	double M = star->Mass();
	for(std::size_t x=0; x<len; x++){
		fprintf(fp, "%0.16le", rad[x]/R);
		for(int a=0; a<num_var; a++) fprintf(fp, "\t%0.16le", y[x][a]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	//plot file in png in gnuplot, and open png
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	std::string title = star->graph_title();
	fprintf(gnuplot, "set title 'full mode %d,%d in %s, period=%0.5lf s'\n", l,k, title.c_str(), getPeriod());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'log|y|'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	std::string varname[numvar]; driver->varnames(varname);
	fprintf(gnuplot, "set arrow 1 from %le, graph 0 to %le, graph 1 lc rgb 'red' nohead\n", rad[xfit]/R, rad[xfit]/R);
	fprintf(gnuplot, "plot ");
	for(int a=0; a<numvar; a++){
		fprintf(gnuplot, "%c '%s' u 1:(abs($%d)) w l t '%s'", (a==0? ' ':','),txtname.c_str(), a+2, varname[a].c_str());
	}
	fprintf(gnuplot, "\n");

	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

//write the mode, and then open the png to screen for easy viewing
template <std::size_t numvar> 
void Mode<numvar>::printMode(const char *const c){
	writeMode(c);
	std::string outname, modename = "mode_"+std::to_string(l)+"."+std::to_string(k)+".png";
	if(c==NULL) outname = "./out/"+star->name+"/mode/"+modename;
	else {
		outname = "./";
		for(std::size_t i=0; c[i]!=0; i++){
			outname += c[i];
		}
		outname += "modes/"+modename;
	}
	char openmyplot[248];
	system(std::string("open "+outname).c_str());
}

//ways to access the frequency
template <std::size_t numvar> 
double Mode<numvar>::getOmega2(){
	return omega2;
}
template <std::size_t numvar> 
double Mode<numvar>::getFreq(){
	return sqrt(omeg2freq * omega2);
}
template <std::size_t numvar> 
double Mode<numvar>::getPeriod(){
	return 2.*m_pi/getFreq();
}

//we want to be able to call SSR() on each Mode object
//however, SSR() requires equations from ModeDriver
template <std::size_t numvar> 
double Mode<numvar>::SSR(){
	return driver->SSR(omega2, l, this);
}

template <std::size_t numvar>
double Mode<numvar>::tidal_overlap(){
	return driver->tidal_overlap(this);
}

#endif