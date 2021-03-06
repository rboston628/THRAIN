//**************************************************************************************
//
//										THRAIN
//
//**************************************************************************************

//  By Reece Boston, Mar 24 2022

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
//constants describing units used for calculation
namespace units {
	struct UnitSet {
		double G;
		double C;
		double base_length;
		double base_time;
		double base_mass;
	};
	enum Units {astro, geo, SI, CGS};
	enum ParamType {pmass=0b00001, pradius=0b00010, pzsurf=0b00100, plogg=0b01000, pteff=0b10000};
}
namespace modetype {enum ModeType {radial, nonradial, cowling};}
namespace error {enum ErrorType {isRMSR=0, isC0, isIsopycnic, isJCD, numerror};}

//an  object specifying input parameters for a calculation
struct CalculationInputData {
	std::string calcname;			//unique user-chosen identifier for calculation
	regime::Regime regime;			//the regime of physics to use: Newtonian, 1PN, GR
	model::StellarModel model;		//the model of star to use in calculations
	units::Units units;				//the choice of units, eg: CGS, SI, astronomical, ...
	modetype::ModeType modetype;	//the kind of modes to calculation: Cowling, 4th order, ...
	//
	double *input_params;			//parameters needed for the stellar model calculation
	int num_input_params;			//the number of parameters needed by this stellar model
	std::string str_input_param;    //a possible string parameter, such as a file name
	int mode_num;					//the number of modes to be calculated
	int Ngrid;						//the number of grid points to use in the calculation
	int *l, *k;						//pointers to arrays for L,K of modes to be calculated
	double mass, radius, zsurf, logg, teff;//properties of a star
	char params;					//records which of the above properties was specified
	double adiabatic_index;			//an adiabatic index to use in oscillations; if 0 uses natural Gamma
	//the deconstructor
	~CalculationInputData(){
		delete[] l;
		delete[] k;
		delete[] input_params;
	}
};

//an object specifying results for output from a calculation
struct CalculationOutputData {
	std::string calcname;			//unique user-chosen identifier for calculation
	regime::Regime regime;			//the regime of physics to use: Newtonian, 1PN, GR
	model::StellarModel model;		//the model of star to use in calculations
	units::Units units;				//the choice of units, eg: CGS, SI, astronomical, ...
	units::UnitSet unitset;			//the base values to use in this unitset to convert back to SI
	modetype::ModeType modetype;	//the kind of modes to calculation: Cowling, 4th order, ...
	//
	double *input_params;			//parameters needed for the stellar model calculation
	int num_input_params;			//the number of parameters needed by this stellar model
	std::string str_input_param;    //a possible string parameter, such as a file name
	int mode_num, mode_writ, mode_done;//records the number of modes, which modes have been written
	int Ngrid;						//the number of grid points to use in the calculation
	int *l, *k;						//pointers to arrays for L,K of modes to be calculated
	double *w, *f, *period;			//points to arrays of frequency and period for each mode
	int i_err;						//a number of errors to calculate for each mode
	double **err;					//an array, for each mode listing all required errors
	double mass, radius, zsurf, logg, teff, star_SSR;
	double *mode_SSR;				//an array, for each mode stating the backsubstitution residual
	double freq0;					//a base frequency given the mass, radiusof star
	double adiabatic_index;			//an adiabatic index to use in oscillates; if 0 uses natural Gamma
	char params;					//records which of M,R,x,logg were specified
	Star* star;						//the stellar model
	ModeDriver* driver;				//the driver for the chosen mode type
	std::vector<ModeBase*> mode;	//an array of modes to be calculated
	// outout error flags -- these flags indicate which columns to print to estimate numerical error
	bool error[error::numerror];
	//the deconstructor
	~CalculationOutputData(){
		delete star;
		delete driver;
		mode.clear();
		delete[] l;
		delete[] k;
		delete[] w;
		delete[] f;
		delete[] period;
		for(int e=0; e<i_err; e++)
			delete[] err[e];
		delete[] input_params;
	}
};


//create the specified star
int create_star(CalculationOutputData&);

//create the specified modes
int create_modes(CalculationOutputData&);

//format units in the chosen unit system
int format_units(CalculationOutputData&);


#endif