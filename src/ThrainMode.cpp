//**************************************************************************************
//							THRAIN MODE
//  ThrainMode.cpp
//	 	Handles all functionality with finding modes and calculating their properies
//  Reece Boston, Mar 24 2022
//**************************************************************************************

#ifndef THRAINMODECPP
#define THRAINMODECPP

#include "ThrainMain.h"
#include "ThrainMode.h"

namespace mode {

//this will create the correct mode driver and pass the correct driver type to the mode_finder
int create_modes(Calculation::OutputData &data_out){	
	switch(data_out.modetype){
		case modetype::cowling:
			data_out.driver = new CowlingModeDriver(data_out.star, data_out.adiabatic_index);
			mode_finder<Mode<2>,CowlingModeDriver>(data_out);
			break;
		case modetype::radial:
		case modetype::nonradial:
			data_out.driver = new NonradialModeDriver(data_out.star, data_out.adiabatic_index);			
			mode_finder<Mode<4>,NonradialModeDriver>(data_out);
			break;
	}	
	return 0;
}

void get_min_from_set(
	const std::set<int>& kfilled, 
	const std::map<int,double>& w2filled, 
	const int& ktarget,
	int& kmin, 
	double& w2min
){
	int khigh=*(kfilled.rbegin());
	// lower bound -- largest element in list smaller than target
	// and is better min than current kmin
	if(ktarget > khigh && khigh > kmin || kmin > ktarget && ktarget > khigh){
		kmin = khigh;
		w2min = w2filled.at(kmin);
	}
	// lower bound -- elements in list could bound target
	auto kmint = kfilled.lower_bound(ktarget);
	if(kmint!=kfilled.end() && (kmint)!=kfilled.begin()){
		// if it found something, use the element prior
		// will be first element equal to or lower than ktarget
		kmin = *prev(kmint);
		w2min = w2filled.at(kmin);
	}
}

void get_max_from_set(
	const std::set<int>& kfilled, 
	const std::map<int,double>& w2filled, 
	const int& ktarget,
	int& kmax,
	double& w2max 
){
	// upper bound -- find first element larger than ktarget
	auto kmaxt = kfilled.upper_bound(ktarget);
	if(kmaxt!=kfilled.end()){
		kmax = *kmaxt;
		w2max = w2filled.at(kmax);
	}
}

double compare_JCD(double n, int l, int k, double w){
	if((n!=1.5) && (n!=3.0) && (n!=4.0) ){
		return nan("");
	}
	//convert dimensionless freuqueny to same scale used in JCD-DJM paper
	w = round(w*nug*10000.0)/10000.0;
	double fJCD = 0.0;
	//for polytropes n=1.5, n=3, n=4, compare to tables,l=1,2,3, with 1<=k<=35
	if((l==1 | l==2 | l==3) & (k>0 & k<36)){
		switch(int(n*10)){
			case 15:
				fJCD = JCD1_5[l-1][k-1];
				break;
			case 30:
				fJCD = JCD3_0[l-1][k-1];
				break;
			case 40:
				fJCD = JCD4_0[l-1][k-1];
				break;
		}
		return abs(w-fJCD);
	}
	else return nan("");
}

//a single function to check against the Pekeris formula
double compare_Pekeris(double w, int l, int k, double Gam1){
	double dnl = Gam1*double(k)*(double(k+l)+0.5)-2.;
	double wPek2 = dnl + sqrt(dnl*dnl + double(l*l+l));
	//there is also a formula for f-modes due to Chandrasekhar (1964), c.f. Cox (1980) eq 17.80
	if(k==0) wPek2 = double(2.*l*(l-1))/double(2*l+1);
	return fabs(w*w - wPek2)/wPek2;
}

} // namespace mode

#endif