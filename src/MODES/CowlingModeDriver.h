//**************************************************************************************
//						COWLING PULSATION nlm MODE DRIVER
// Cowling ModeDriver.h
//		Solves the Newtonian LAWE in Dziembowski variables in the Cowling approximation
//		Does not consider perturbtations to gravitational field, making 2nd-order modes
//  Reece Boston Mar 24, 2022
//***************************************************************************************

#ifndef COWLINGMODEDRIVERH
#define COWLINGMODEDRIVERH

#include "ModeDriver.h"

class CowlingModeDriver : public ModeDriver {
public:
	static const unsigned int num_var=2U;
	
	//constructor
	//initialize from a background star and an adiabatic index
	CowlingModeDriver(Star*, double);
	//owns raw arrays; copying would double-free
	CowlingModeDriver(CowlingModeDriver const&) = delete;
	CowlingModeDriver& operator=(CowlingModeDriver const&) = delete;
	//destructor
	virtual ~CowlingModeDriver();
	
	std::size_t length() const override final;
	double Gamma1() const override final;//the adiabatic index, if set; 0 if same as star
	double rad(std::size_t) const override final;//radial array x = r/R
	//return the coefficient matrix A from equation x*dy/dx = A*y
	void getCoeff(double *CC, const std::size_t, const int, const double, const int) override final;

private:
	std::size_t len = 0;		//number of grid points for mode
	std::size_t len_star = 0;	//number of grid points in star
	double adiabatic_index;	//adiabatic index; set to 0 to use star's Gamma1
	enum VarNames {y1=0, y2};

	//perturbation quantities
	double *r = nullptr, *A = nullptr, *U = nullptr, *C = nullptr, *V = nullptr;
	void initializeArrays();
	
	static const int BC_C = 4;	//the desired order near the center
	static const int BC_S = 4;	//the desired order near the surface
	//expansion coefficients near surface;
	double As[BC_S+1], Vs[BC_S+1], cs[BC_S+1], Us[BC_S+1], cProds[BC_S+1], k_surface = 0.0;
	//expansion coefficients near center;
	double Ac[BC_C/2+1], Vc[BC_C/2+1], cc[BC_C/2+1], Uc[BC_C/2+1], cProdc[BC_C/2+1];
	void setupBoundaries() override final;
	std::size_t CentralBC(double **y, double *yo, double s2, int l, int m=0) override final;
	std::size_t SurfaceBC(double **y, double *yo, double s2, int l, int m=0) override final;

public:	
	void getBoundaryMatrix(int, double**, int*) override final;
	void varnames(std::string *names) override final {
		names[0] = "y1"; names[1]="y2";
	}
	
	//for the Mode passed, calculate sum-square residual
	double SSR(double, int l, ModeBase const*) const override final;
	double tidal_overlap(ModeBase const*) const override final;
	double innerproduct(ModeBase const*, ModeBase const*) const override final;
};

#endif