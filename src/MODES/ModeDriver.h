// Reece Boston
// Mar 24, 2022
#ifndef MODEDRIVERH
#define MODEDRIVERH

#include "../constants.h"
#include "../STARS/Star.h"


//**************************************************************************************
//							MODE ABSTRACT BASE CLASS
//		This is an abstract base class for all Mode objects, allowing polymorphism
//		Can be used in code as a generic mode, without regard to particular physics
//		The generic Mode<N> defined in Mode.h cannot always be used polymorphically
//		When a mode is needed and exact physics unknown, use ModeBase as type
//**************************************************************************************
class ModeBase {
public:
	//file output methods to write and view plots of mode
	virtual void printMode(const char *const c = NULL) = 0;
	virtual void writeMode(const char *const c = NULL) = 0;
		
	virtual ~ModeBase(){};
	
	virtual int modeOrder() = 0;
	virtual void modeNumbers(int&, int&, int&) =0;
	virtual double getOmega2() = 0;
	virtual double SSR() = 0;
	virtual double tidal_overlap() =0;
	virtual double getFreq() =0;
	virtual double getPeriod()=0;
	
	virtual double getRad(std::size_t x) =0;
	virtual double getY(int i, std::size_t x) =0;
	virtual double getYtilde(int i, std::size_t x) =0;
};

//**************************************************************************************
//							MODE DRIVER ABSTRACT CLASS
// ModeDriver.h
//		An abstract class for all sets of mode equations to inherit
//		A ModeDriver object will provide the physical equations
//			and he boundary conditions obeyed by the individual modes
//		A single ModeDriver is capable of serving many Modes
//		Serves as the basis in code for other models of pulsations
//**************************************************************************************
class ModeDriver {
public:
	const std::size_t num_var;
	//constructor
	ModeDriver(int nv, Star *s) : num_var(nv), star(s)  {};
	virtual ~ModeDriver(){};
	virtual std::size_t length() =0;
	virtual double Gamma1() =0;
	virtual double rad(std::size_t) =0;

	virtual std::size_t CentralBC(double **y, double *yo, double s2, int l, int m=0) =0;
	virtual std::size_t SurfaceBC(double **y, double *ys, double s2, int l, int m=0) =0;
	virtual void getCoeff(double *CC, const std::size_t, const int, const double, const int) =0;
	virtual void setupBoundaries() =0;
	
	virtual double SSR(double, int, ModeBase*) =0;
	virtual double tidal_overlap(ModeBase*) =0;
	virtual double innerproduct(ModeBase*, ModeBase*) =0;
	
	//virtual double test_eigenvalue(double, int) =0;
	
	//the following two methods are added to make Mode agnostic
	//provide the order of BCs to check in forming Wronskian
	virtual void getBoundaryMatrix(int, double **, int*) =0;
	virtual void varnames(std::string*) =0;	//names of variables to print out

protected:
	Star *const star;
	int central_bc_order;
	int surface_bc_order;
	
	friend class ModeBase;
	template <std::size_t N> friend class Mode;
	friend class Star;
};

#endif