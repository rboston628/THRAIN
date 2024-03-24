//**************************************************************************************
//							STAR ABSTRACT CLASS
// Star.cpp
//		Implements the output and sum-square residual functions for equilibrium stars
//		These methods are general enough to work with any sub-class of star
//**************************************************************************************

#ifndef STARCLASS
#define STARCLASS

#include <limits>
#include "Star.h"
#include "../lib/string.h"

//in Newtonian, light speed is infinity
double Star::light_speed2() {
	return std::numeric_limits<double>::infinity();
}

//method to print pertinent values of star to .txt, and plot them in gnuplot
void Star::writeStar(const char *const c){
	//create names for files to be opened
	// char pathname[256];
	std::string outputdir;
	if(c==NULL)	outputdir = strmakef("./out/%s/", name.c_str());
	else outputdir = strmakef("./%s/star/", c);
	
	system( ("mkdir -p "+outputdir).c_str() );
	
	printStar(outputdir.c_str());
	printBV(outputdir.c_str());
	printCoefficients(outputdir.c_str());
}

//method to print pertinent values of star to .txt, and plot them in gnuplot
void Star::printStar(const char *const outputdir){
	//create names for files to be opened
	std::string txtname = strmakef("%s/%s.txt", outputdir, name.c_str());
	std::string outname = strmakef("%s/%s.png", outputdir, name.c_str());

	FILE *fp;
	if(!(fp = fopen(txtname.c_str(), "w")) ){
		system( addstring("mkdir -p ",outputdir).c_str() );
		fp = fopen(txtname.c_str(), "w");
	}
	//print results to text file
	// radius rho pressure gravity
	double irc=1./this->rho(0), ipc=1./P(0), R=Radius(), ig=1./dPhidr(length()-1);
	for(std::size_t X=0; X< length(); X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R, this->rho(X)*irc, -drhodr(X)*irc*R,
			P(X)*ipc, -dPdr(X)*ipc*R,
			mr(X)/Mass(), dPhidr(X)*ig);
		fflush(fp);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1600,800\n");
	fprintf(gnuplot, "set samples %lu\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Profile for %s'\n", graph_title().c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'rho/rho_c, P/P_c, m/M, g/g_S'\n");
	fprintf(gnuplot, "set yrange [0:1]\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'rho'", txtname.c_str());
	fprintf(gnuplot, ", '%s' u 1:4 w l t 'P'", txtname.c_str());
	fprintf(gnuplot, ", '%s' u 1:6 w l t 'm'", txtname.c_str());
	fprintf(gnuplot, "\n");
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);
}

//calculte a root-mean-square scaled residual from backsubstitution into physical equations
//  note the name SSR (sum-square residual) is a misnomer
double Star::SSR(){	
	double checkEuler;
	double checkPoiss;
	std::size_t len = length();
			
	//sum up errors in equations
	checkEuler = 0.0;
	checkPoiss = 0.0;
	double d2Phi = 0.0;
	double e1, e2, n1, n2;
	for(std::size_t X=4; X<len-4; X++){
		//Euler equation
		e1 = fabs(dPdr(X) + rho(X)*dPhidr(X) );
		n1 = fabs(dPdr(X)) + fabs(rho(X)*dPhidr(X));
		//Poisson equation
		//calculate numerical derivatives
		double a3,a2,a1,b1,b2,b3, h;
		// d2Phi/dr2 = (d/dr)(dPhi/dr)
		b3=dPhidr(X-3);
		b2=dPhidr(X-2);
		b1=dPhidr(X-1);
		a1=dPhidr(X+1);
		a2=dPhidr(X+2);
		a3=dPhidr(X+3);
		h=rad(X+1)-rad(X);
		d2Phi = (45.*a1-9.*a2+a3-45.*b1+9.*b2-b3)/(60.*h);
		e2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X) - 2.0*dPhidr(X) - d2Phi*rad(X) );
		n2 = fabs( 4.0*Gee()*m_pi*rho(X)*rad(X)) + fabs(2.0*dPhidr(X)) + fabs(d2Phi*rad(X) );
		//add absolute error
		e1 = e1/n1;
		e2 = e2/n2;
		checkEuler += e1*e1;
		checkPoiss += e2*e2;
	}
	return sqrt((checkPoiss+checkEuler)/double(2*len-8));
}

void Star::printBV(const char *const outputdir, double const gam1){
	std::string txtname = addstring(outputdir, "/BruntVaisala.txt");
	std::string outname = addstring(outputdir, "/BruntVaisala.png");
	
	//print the Brunt-Vaisala frequency
	FILE* fp  = fopen(txtname.c_str(), "w");
	double N2 = -1.0, R = Radius();
	fprintf(fp, "1-m\tN2\tL1\n");
	for(std::size_t X=1; X< length()-1; X++){
		N2 =  -Schwarzschild_A(X,gam1)*dPhidr(X);
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/R,
			(1.-mr(X)/Mass()),
			N2,
			2.*sound_speed2(X,gam1)*pow(rad(X),-2));
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE* gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Brunt-Vaisala and Lamb for %s'\n", graph_title().c_str());
	fprintf(gnuplot, "set xlabel '-log_{10}(1-m/M)'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} N^2 & log_{10} L_1^2 (Hz^2)\n");
	fprintf(gnuplot, "set logscale x 10\n");
	fprintf(gnuplot, "set logscale y 10\n");
	fprintf(gnuplot, "set format x '%%L'\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set ytics 10\n");
	fprintf(gnuplot, "plot '%s' u 2:3 w l t 'N^2'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 2:4 w l t 'L_1^2'", txtname.c_str());
	//fprintf(gnuplot, ",    '%s' u 1:(-$2) w l t '-N^2'", txtname);
	fprintf(gnuplot, "\n");
	pclose(gnuplot);
}

void Star::printCoefficients(const char *const outputdir, double const gam1){

	std::string txtname, outname;
	std::string title = graph_title();
	
	// FIRST: print out coefficients near center and surface, for series analysis
	std::string wavecoeffdir = addstring(outputdir, "/wave_coefficient");
	system( ("mkdir -p " + wavecoeffdir).c_str() ); 

	std::size_t Ntot = length();
	const int num_c=3, num_s=5;
	int bc_c = (num_c-1)*2, bc_s = (num_s-1);
	double A0[num_c] = {0.}, U0[num_c] = {0.}, V0[num_c] = {0.}, c0[num_c] = {0.};
	double A1[num_s] = {0.}, U1[num_s] = {0.}, V1[num_s] = {0.}, c1[num_s] = {0.};
	//get the central coefficients
	getAstarCenter(A0, bc_c, gam1);
	getUCenter(U0,bc_c);
	getVgCenter(V0,bc_c, gam1);
	getC1Center(c0,bc_c);
	//get the surface coefficients
	getAstarSurface(A1, bc_s, gam1);
	getUSurface(U1,bc_s);
	getVgSurface(V1,bc_s, gam1);
	getC1Surface(c1,bc_s);
	
	//print the coefficients of the center and surface, for series analysis
	txtname = wavecoeffdir + "/center.txt";
	FILE *fp = fopen(txtname.c_str(), "w");
	fprintf(fp, "A*:\t%0.16le\t%0.16le\t%0.16le\n", A0[0],A0[1],A0[2]);
	fprintf(fp, "U :\t%0.16le\t%0.16le\t%0.16le\n", U0[0],U0[1],U0[2]);
	fprintf(fp, "Vg:\t%0.16le\t%0.16le\t%0.16le\n", V0[0],V0[1],V0[2]);
	fprintf(fp, "c1:\t%0.16le\t%0.16le\t%0.16le\n", c0[0],c0[1],c0[2]);
	fclose(fp);
	txtname = wavecoeffdir + "/surface.txt";
	fp = fopen(txtname.c_str(), "w");
	fprintf(fp, "A*:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", A1[0],A1[1],A1[2],A1[3],A1[4]);
	fprintf(fp, "U :\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", U1[0],U1[1],U1[2],U1[3],U1[4]);
	fprintf(fp, "Vg:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", V1[0],V1[1],V1[2],V1[3],V1[4]);
	fprintf(fp, "c1:\t%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n", c1[0],c1[1],c1[2],c1[3],c1[4]);
	fclose(fp);
	
	//print fits to those coefficients at center and surface
	std::size_t const NC=15, NS=15;
	txtname = wavecoeffdir + "/centerfit.txt";
	fp = fopen(txtname.c_str(), "w");
	double Rtot = Radius();
	double x, x2;
	fprintf(fp, "x\tA\tA*_fit\tU\tU_fit\tVg\tVg_fit\tc1\tc1_fit\n");
	for(std::size_t X=0; X<NC; X++){
		x = rad(X)/Rtot;
		x2 = pow(x,2);
		fprintf(fp, "%0.8le\t%0.8le\t%0.8le\t%0.8le\t%0.8le\t%0.8le\t%0.8le\t%0.8le\t%0.8le\n",
			x,
			getAstar(X, gam1),
			A0[0] + A0[1]*x2 + A0[2]*x2*x2,
			getU(X),
			U0[0] + U0[1]*x2 + U0[2]*x2*x2,
			getVg(X,gam1),
			V0[0] + V0[1]*x2 + V0[2]*x2*x2,
			getC(X),
			c0[0] + c0[1]*x2 + c0[2]*x2*x2
		);
	}
	fclose(fp);
	txtname = wavecoeffdir + "/surfacefit.txt";
	fp = fopen(txtname.c_str(), "w");
	double t = 0.;
	fprintf(fp, "t%*c\tA*      \tA*_fit   \tU       \tU_fit   \tVg      \tVg_fit  \tc1      \tc1_fit  \n", 7, ' ');
	for(std::size_t X=Ntot-2; X>=Ntot-NS-1; X--){
		t = 1. - rad(X)/Rtot;
		fprintf(fp, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
			t,
			getAstar(X, gam1),
			A1[0]/t + A1[1] + A1[2]*t + A1[3]*t*t + A1[4]*t*t*t,
			getU(X),
			U1[0] + U1[1]*t + U1[2]*t*t + U1[3]*t*t*t + U1[4]*t*t*t*t,
			getVg(X, gam1),
			V1[0]/t + V1[1] + V1[2]*t + V1[3]*t*t + V1[4]*t*t*t,
			getC(X),
			c1[0] + c1[1]*t + c1[2]*t*t + c1[3]*t*t*t + c1[4]*t*t*t*t
		);
	}
	fclose(fp);
	
	//print the pulsation coeffcients
	txtname = wavecoeffdir + "/coefficients.txt";
	outname = wavecoeffdir + "/coefficients.png";
	fp  = fopen(txtname.c_str(), "w");
	fprintf(fp, "m\tA*\tU\tVg\tc1\n");
	for(std::size_t X=0; X<Ntot; X++){
		fprintf(fp, "%0.16le\t%0.16le\t%0.16le\t%0.16le\t%0.16le\n",
			rad(X)/Rtot,
			getAstar(X, gam1),
			getU(X),
			getVg(X, gam1),
			getC(X)
		);
	}
	fclose(fp);	
	//plot file in png in gnuplot
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	fprintf(gnuplot, "set samples %lu\n", length());
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Pulsation Coefficients for %s'\n", title.c_str());
	//fprintf(gnuplot, "set xlabel 'log_{10} r/R'\n");
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'A*, U, V_g, c_1'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot '%s' u 1:2 w l t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:3 w l t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:4 w l t 'V_g'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:5 w l t 'c_1'", txtname.c_str());
	fprintf(gnuplot, "\n");\
	//fits
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 1000,800\n");
	txtname = wavecoeffdir + "/centerfit.txt";
	outname = wavecoeffdir + "/centerfit.png";
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Central Fitting by Power Series for %s'\n", title.c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set yrange [1e-8:10]\n");
	fprintf(gnuplot, "set ytics ('0' 1e-8, '1e-6' 1e-6, '1e-4' 1e-4, '1e-2' 1e-2, '1' 1e0, '100' 1e2)\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [0:%le]\n", 1.01*rad(NC)/Rtot);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'Vg'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'c1'", txtname.c_str());
	fprintf(gnuplot, "\n");
	txtname = wavecoeffdir + "/surfacefit.txt";
	outname = wavecoeffdir + "/surfacefit.png";
	fprintf(gnuplot, "set output '%s'\n", outname.c_str());
	fprintf(gnuplot, "set title 'Surface Fitting by Power Series for %s'\n", title.c_str());
	fprintf(gnuplot, "set xlabel 'r/R'\n");
	fprintf(gnuplot, "set ylabel 'difference'\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set ytics 100\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "set xrange [0:%le]\n", 1.-rad(Ntot-NS-1)/Rtot);
	fprintf(gnuplot, "plot '%s' u 1:(abs($2-$3)/abs($2)) w lp t 'A*'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($4-$5)/abs($4)) w lp t 'U'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($6-$7)/abs($6)) w lp t 'Vg'", txtname.c_str());
	fprintf(gnuplot, ",    '%s' u 1:(abs($8-$9)/abs($8)) w lp t 'c1'", txtname.c_str());
	fprintf(gnuplot, "\n");
	
	//now leave gnuplot
	fprintf(gnuplot, "exit\n");
	pclose(gnuplot);	
}

#endif
