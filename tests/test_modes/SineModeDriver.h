//**************************************************************************************
//						SINE MODE DRIVER
// SineModeDriver.h
//		Solves simple sinusoidal equation for sine waves.
//	THIS IS NOT PART OF THRAIN -- ONLY FOR CERTAIN TESTS
//  Reece Boston Mar 14, 2024
//***************************************************************************************

#ifndef SINEMODEDRIVERH
#define SINEMODEDRIVERH

#include "../../src/MODES/ModeDriver.h"

class SineModeDriver : public ModeDriver {
public:
	static const std::size_t num_var=2UL;
	
	//constructor
	//initialize from a background star
	SineModeDriver(Star*);
	//destructor
	virtual ~SineModeDriver();
	
	std::size_t length() override;
	double Gamma1() override;//the adiabatic index, if set; 0 if same as star
	double rad(std::size_t) override;//radial array x = r/R
	//return the coefficient matrix A from equation x*dy/dx = A*y
	void getCoeff(double *CC, const std::size_t, const int, const double, const int) override;

private:
	std::size_t len;		//number of grid points for mode
	std::size_t len_star;	//number of grid points in star
	enum VarNames {sin=0, cos};

	double x(std::size_t const, int const);
		
	void setupBoundaries() override;
	std::size_t CentralBC(double **y, double *yo, double s2, int l, int m=0) override;
	std::size_t SurfaceBC(double **y, double *yo, double s2, int l, int m=0) override;

public:	
	void getBoundaryMatrix(int, double**, int*) override;
	void varnames(std::string *names) override;
	
	//for the Mode passed, calculate sum-square residual
	double SSR(double, int l, ModeBase*) override;
	double tidal_overlap(ModeBase*) override;
	double innerproduct(ModeBase*, ModeBase*) override;
};

#endif