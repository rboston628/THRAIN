//**************************************************************************************
//  ThrainUnits.cpp
//	This file contains all fuctions for handling unit conversions for the program
//  Reece Boston, Mar 24, 2022
//**************************************************************************************

#include "ThrainMain.h"

#ifndef THRAINUNITSH
#define THRAINUNITSH

namespace units{

inline double GMR2FromLogg(const Calculation::OutputData& data){
	// logg is a dimensioned quantity -- return GM/R^2 in indicated units
	return pow(10.0,data.logg) * pow(data.unitset.base_length,2)/data.unitset.base_mass*data.unitset.G/G_CGS;
}

inline double getLoggFromRM (const Calculation::OutputData& data){
	// unfortunately this is a dimensioned quantity, and must be given in CGS units
	return log10(G_CGS*data.mass*data.unitset.base_mass*pow(data.unitset.base_length*data.radius,-2));
}
inline double getZsurfFromRM(const Calculation::OutputData& data){
	return 1./sqrt( 1. - 2.*data.unitset.G*data.mass/(data.radius*pow(data.unitset.C,2)) ) - 1.;
}

inline double getRadiusFromZM(const Calculation::OutputData& data){
	return 2.*data.unitset.G*data.mass*pow(data.unitset.C,-2)/(1.-pow(1.+data.zsurf,-2));
}

inline double getRadiusFromLoggM(const Calculation::OutputData& data){
	return sqrt(data.unitset.G*data.mass / GMR2FromLogg(data));
}

inline double getMassFromRZ(const Calculation::OutputData& data){
	return 0.5*pow(data.unitset.C,2)*data.radius/data.unitset.G*(1.-pow(1.+data.zsurf,-2));
}

inline double getMassFromRLogg(const Calculation::OutputData& data){
	return pow(data.radius,2)*GMR2FromLogg(data)/data.unitset.G;
}

int format_units(Calculation::OutputData& data){	
	//fix the unitset based on which united were specified
	data.unitset = units::unitSets.at(data.units);
	
	//the user can specify certain parameters of the star for the model; this calculates the others
	switch(data.params){
		case (ParamType::pmass|ParamType::pradius):
			data.logg = getLoggFromRM(data);
			data.zsurf = getZsurfFromRM(data);
			break;
		case (ParamType::pmass|ParamType::pzsurf):
			data.radius = getRadiusFromZM(data);
			data.logg = getLoggFromRM(data);
			break;
		case (ParamType::pmass|ParamType::plogg):
			data.radius = getRadiusFromLoggM(data);
			data.zsurf = getZsurfFromRM(data);
			break;
		case (ParamType::pradius|ParamType::pzsurf):
			data.mass = getMassFromRZ(data);
			data.logg = getLoggFromRM(data);
			break;
		case (ParamType::pradius|ParamType::plogg):
			data.mass = getMassFromRLogg(data);
			data.zsurf = getZsurfFromRM(data);
			break;
		case (ParamType::pzsurf|ParamType::plogg):
			//I don't know why this isn't implemented yet -- I don't feel like doing it now
			data.mass   = 0.0;
			data.radius = 0.0;
			throw std::runtime_error("Unimplemented paramsets zsurf and logg, try another\n");
			break;
		case (ParamType::pmass|ParamType::pteff):
			data.radius = data.zsurf = data.logg = 0.0;
			break;
		case 0:
		default:
			data.mass = data.radius = data.zsurf = data.logg = 0.0;
			break;
	}
	
	printf("\tmass=%lg\n", data.mass);
	printf("\tradius=%lg\n", data.radius);
	printf("\tzsurf=%lg\n\tlogg=%lg\n", data.zsurf, data.logg);
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = data.radius*data.unitset.base_length;
	data.freq0 = sqrt(G_CGS*mass_CGS*pow(radius_CGS,-3));
	printf("\tfrequency scale = %lf\n", data.freq0);
	
	return 0;
}
}


#endif