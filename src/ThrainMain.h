//**************************************************************************************
//
//										THRAIN
//
//**************************************************************************************

//  By Reece Boston, Mar 24 2022

//basic
#include <map>
#include <set>
#include "ThrainUnits.h"
//basic template classes
#include "constants.h"
#include "STARS/Star.h"
#include "MODES/Mode.h"
//specific Newtonian stars
#include "STARS/Polytrope.h"
#include "STARS/ChandrasekharWD++.h"
#include "STARS/MESA.h"
#include "STARS/SimpleWD.h"
//mode drivers
#include "MODES/ModeDriver.h"
#include "MODES/CowlingModeDriver.h"
#include "MODES/NonradialModeDriver.h"

#ifndef HEADER
#define HEADER

//constants labeling the physical regime being used
namespace regime { enum Regime {PN0};}
//constants labeling the stellar models used 
namespace model {enum StellarModel {polytrope, CHWD, MESA, SWD};}
namespace modetype {enum ModeType {radial, nonradial, cowling};}
namespace error {enum ErrorType {isRMSR=0, isC0, isIsopycnic, isJCD, numerror};}

//an  object specifying input parameters for a calculation
namespace Calculation{
struct InputData {
	std::string calcname;            //unique user-chosen identifier for calculation
	regime::Regime regime;           //the regime of physics to use: Newtonian, 1PN, GR
	model::StellarModel model;       //the model of star to use in calculations
	units::Units units;              //the choice of units, eg: CGS, SI, astronomical, ...
	modetype::ModeType modetype;     //the kind of modes to calculation: Cowling, 4th order, ...
	std::vector<double> input_params;//parameters needed for the stellar model calculation
	std::string str_input_param;     //a possible string parameter, such as a file name
	int mode_num;                    //the number of modes to be calculated
	int Ngrid;                       //the number of grid points to use in the calculation
	std::set<int> l;                 //the set of L in calculation
	std::map<int,std::set<int>> kl;  //for each L, the set of K to search
	double mass, radius, zsurf, logg, teff;//properties of a star
	char params;                     //records which of the above properties was specified
	double adiabatic_index;          //an adiabatic index to use in oscillations; if 0 uses natural Gamma
};

//an object specifying results for output from a calculation
struct OutputData {
	std::string calcname;			 //unique user-chosen identifier for calculation
	regime::Regime regime;			 //the regime of physics to use: Newtonian, 1PN, GR
	model::StellarModel model;		 //the model of star to use in calculations
	units::Units units;				 //the choice of units, eg: CGS, SI, astronomical, ...
	units::UnitSet unitset;			 //the base values to use in this unitset to convert back to SI
	modetype::ModeType modetype;	 //the kind of modes to calculation: Cowling, 4th order, ...
	//
	std::vector<double> input_params;//parameters needed for the stellar model calculation
	std::string str_input_param;     //a possible string parameter, such as a file name
	int mode_num, mode_writ, mode_done;//records the number of modes, which modes have been written
	int Ngrid;						 //the number of grid points to use in the calculation
	std::vector<int> l, k;           //vectors for L,K for each mode
	std::vector<double> w, f, period;//vectors of frequency and period for each mode
	std::vector<double> mode_SSR;    //the backsubstitution residual for each mode
	int i_err;						 //a number of errors to calculate for each mode
	double **err;					 //an array, for each mode listing all required errors
	double mass, radius, zsurf, logg, teff, star_SSR;
	double freq0;					 //a base frequency given the mass, radiusof star
	double adiabatic_index;			 //an adiabatic index to use in oscillates; if 0 uses natural Gamma
	char params;					 //records which of M,R,z,logg were specified
	Star* star;						 //the stellar model
	ModeDriver* driver;				 //the driver for the chosen mode type
	std::vector<ModeBase*> mode;	 //an array of modes to be calculated
	// outout error flags -- these flags indicate which columns to print to estimate numerical error
	bool error[error::numerror];
	//the deconstructor
	~OutputData(){
		delete star;
		delete driver;
		for(int e=0; e<i_err; e++)
			delete[] err[e];
	}
};

} // namespace Calculaton

//create the specified star
int create_star(Calculation::OutputData&);

#endif