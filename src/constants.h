// **************************************************************************************
// constants.h
//  	This file describes dependencies, and constants
//  	Many constant taken from https://physics.nist.gov/cuu/Constants/index.html	
//  Reece Boston, Mar 24, 2022
// **************************************************************************************

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <functional>
#include "../lib/Splinor.h"
#include "../lib/matrix.h"
#include "../lib/string.h"
#include "../lib/logger.h"

//constants and dependencies
#ifndef CONSTANTS
#define CONSTANTS

#if defined(_WIN32)
#  define popen  _popen
#  define pclose _pclose
#endif

//good ol' pi!
constexpr double m_pi = 3.1415926535897932384626433832795;
//speed of light in CGS
constexpr double C_CGS = 2.99792458e10;// in CGS
//gravitational constant in CGS
// previous value G = 6.6725985e-8 taken from: https://iopb.res.in/~somen/constants0.html
// new value taken from: https://physics.nist.gov/cuu/Constants/index.html
constexpr double G_CGS = 6.67430e-8; // in CGS
//Avogagro's number
constexpr double N_Avogadro = 6.02214076e23;

//Boltzmann constant, Stefan-Boltzmann constant, radiation constant
constexpr double boltzmann_k = 1.3806503e-16;// in CGS
constexpr double boltzmann_sigma = 5.670374419e-5;// in CGS
constexpr double radiation_a = 7.5657e-15   ;// in CGS

//this is a fundamental frequency sqrt(GM/R^3) of the Sun, in uHz
//taken from the paper by JCD-DJM and used for comparison to their tables
constexpr double nug = 99.855377; //from JCD, in uHz
//other properties of the sun
constexpr double MSOLAR = 1.989e33;//in CGS
constexpr double RSOLAR = 6.957e10 ;//in CGS
constexpr double LSOLAR = 3.828e33;//in CGS
//radius of earth, most useful length scale for WD models
constexpr double REARTH = 6.378e8 ;//in CGS

//Planck constants in CGs units 
constexpr double planck_h_CGS = 6.62607015e-27;
constexpr double planck_hbar_CGS = 1.054571817e-27;

//properties of an electron
constexpr struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} electron = {9.1095e-28, 4.80320425e-10, planck_hbar_CGS/(9.1095e-28*C_CGS)};

//useful properties of a proton
constexpr struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} proton = {1.6726e-24, 4.80320425e-10, planck_hbar_CGS/(1.6726e-24*C_CGS)};


//in astronomical units, things are in terms of
//	base mass = MSOLAR
//  base disance = km (most appropriate for WDs and NSs, definitely not parsec)
//  base speed = km/s
constexpr double C_astro = 2.99792458e5; // km/s
constexpr double G_astro = 4.302e-3*30856775812800.0;// km (km/s)^2 /MSOLAR 
//constexpr double G_astro = 4.302e-3;// pc (km/s)^2 /MSOLAR 
//constexpr double G_astro = 1.90809e5;// RSOLAR (km/s)^2/MSOLAR
#endif
