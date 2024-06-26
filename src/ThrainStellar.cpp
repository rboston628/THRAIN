//**************************************************************************************
//							THRAIN STELLAR
//  ThrainStellar.cpp
//		Handles all functionality with finding background models
//  Reece Boston, Mar 24, 2022
//**************************************************************************************


#include "ThrainMain.h"

//functions to create the different stellar models
int create_classical_polytrope(Calculation::OutputData &);
int create_classical_CHWD(Calculation::OutputData &);
int create_classical_MESA(Calculation::OutputData &);
int create_classical_SWD(Calculation::OutputData &);

//will select the correct creation fuction from above based on user input
int create_star(Calculation::OutputData &data_out){
	printf("Creating star...\n");
	switch(data_out.model) {
		case model::polytrope:
			create_classical_polytrope(data_out);
			break;
		case model::CHWD:
			create_classical_CHWD(data_out);
			break;
		case model::MESA:
			create_classical_MESA(data_out);
			break;
		case model::SWD:
			create_classical_SWD(data_out);
			break;
	}
	return 0;
}

//create a classical polytrope
int create_classical_polytrope(Calculation::OutputData& data){
	double mass_CGS = data.mass*data.unitset.base_mass;
	double radius_CGS = data.radius*data.unitset.base_length;
	//create star	
	data.star = new Polytrope(mass_CGS, radius_CGS, data.input_params[0], data.Ngrid);

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}


//create a classical WD following Chandrasekhar
int create_classical_CHWD(Calculation::OutputData& data){
	switch((int)data.input_params[1]){
		case 0:
			data.star = new ChandrasekharWD(data.input_params[0], data.Ngrid, Chandrasekhar::constant_mu{2.0});
			break;
		case 1: {
			double Y0 = data.input_params[0];
			double X0 = sqrt(Y0*Y0-1.);
			double F0 = Chandrasekhar::factor_f(X0);
			Chandrasekhar::sigmoidal_in_logf func {2.0, F0, 2.0, 1.0};
			data.star = new ChandrasekharWD(data.input_params[0], data.Ngrid, func);
		} break;
		default:
			data.star = new ChandrasekharWD(data.input_params[0], data.Ngrid,  Chandrasekhar::constant_mu{1.0});
	}
	
	//adjust inputs to match the actual values
	data.mass = data.star->Mass()/data.unitset.base_mass;
	data.radius = data.star->Radius()/data.unitset.base_length;
	data.params = units::ParamType::pmass|units::ParamType::pradius;
	units::format_units(data);

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}


//create a classical WD found in MESA
int create_classical_MESA(Calculation::OutputData& data){
	std::string inputname = data.str_input_param + ".dat";
	data.star = new MESA(inputname.c_str(), data.Ngrid);
	
	//adjust the inputs around the fact this is a MESA object
	data.mass = data.star->Mass()/data.unitset.base_mass;		//the mass is determined by model
	data.radius = data.star->Radius()/data.unitset.base_length;	//the radius is determined by model
	data.Ngrid = data.star->length();
	data.params = units::ParamType::pmass|units::ParamType::pradius;
	units::format_units(data);
	
	//calculate error estimate for this MESA model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}

//create a classical SimpleWD
int create_classical_SWD(Calculation::OutputData& data){
	printf("mass %lf\n", data.mass);
	printf("teff %lf\n", data.teff);
	data.star = new SimpleWD(
		data.mass, 
		data.teff,
		data.Ngrid
	);

	//adjust inputs to match the actual values
	data.mass = data.star->Mass()/data.unitset.base_mass;
	data.radius = data.star->Radius()/data.unitset.base_length;
	//calculate zsurf, logg
	data.params = units::ParamType::pmass|units::ParamType::pradius;
	units::format_units(data);
	data.params = units::ParamType::pmass|units::ParamType::pteff;

	//calculate error estimate for this model
	data.star_SSR = data.star->SSR();
	printf("SSR = %le\n", data.star_SSR);
	
	return 0;
}



