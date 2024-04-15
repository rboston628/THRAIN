//**************************************************************************************
//  ThrainUnits.h
//	This file contains all fuctions for handling unit conversions for the program
//  Reece Boston, Aug 20, 2023
//**************************************************************************************

#include "constants.h"
#include <unordered_map>

#ifndef HEADER

namespace Calculation {
    struct OutputData;
}
#endif

//constants describing units used for calculation
#ifndef UNITSH
#define UNITSH

namespace units {

struct UnitSet {
	double G;
	double C;
	double base_length;
	double base_time;
	double base_mass;
};
enum class Units : char {astro, geo, SI, CGS};
enum class ParamType : char {pmass=0b00001, pradius=0b00010, pzsurf=0b00100, plogg=0b01000, pteff=0b10000};

constexpr char operator| (const ParamType& a, const ParamType& b) {return char(a)|char(b);}
constexpr char operator& (const ParamType& a, const ParamType& b) {return char(a)&char(b);}
constexpr char operator| (const char& a, const ParamType& b) {return a|char(b);}
constexpr char operator& (const char& a, const ParamType& b) {return a&char(b);}
constexpr void operator|= (char& a, const ParamType& b) {a = a|char(b);}
constexpr void operator&= (char& a, const ParamType& b) {a = a&char(b);}

const std::unordered_map<Units, UnitSet> unitSets = {
	// base-length 1km = 1e5 cm
	{Units::astro, {G_astro, C_astro, 1.e5, 1.0, MSOLAR}},
	// length in light-seconds, mass to be consistent
	{Units::geo, {1.0, 1.0, 1.0/C_CGS, 1.0, C_CGS/G_CGS}},
	// standard international units
	{Units::SI, {G_CGS*1.0e-3, C_CGS*1.0e-2, 100, 1.0, 1000.0}},
	{Units::CGS, {G_CGS, C_CGS, 1.0, 1.0, 1.0}},
	// the default, erroneous case
	{Units(-85), {0.0, 0.0, 0.0, 0.0, 0.0}}
};

//format units in the chosen unit system
int format_units(Calculation::OutputData&);

// initialize ungiven params from pair given
double GMR2FromLogg(const Calculation::OutputData& data);
double getLoggFromRM (const Calculation::OutputData& data);
double getZsurfFromRM(const Calculation::OutputData& data);
double getRadiusFromZM(const Calculation::OutputData& data);
double getRadiusFromLoggM(const Calculation::OutputData& data);
double getMassFromRZ(const Calculation::OutputData& data);
double getMassFromRLogg(const Calculation::OutputData& data);

}

#endif