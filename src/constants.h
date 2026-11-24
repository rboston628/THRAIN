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
#include "../lib/mathlib.h"

//constants and dependencies
#ifndef CONSTANTS
#define CONSTANTS

#if defined(_WIN32)
#  define popen  _popen
#  define pclose _pclose
#endif

//good ol' pi!
constexpr double m_pi = math::PI;

/** CONSTANTS
 * Values taken from NIST CODATA 2022 
 * https://physics.nist.gov/cuu/Constants/index.html
 * All values in CGS units, unless noted otherwise
 * 
 * Some values had originally been taken from here:
 * https://iopb.res.in/~somen/constants0.html
 * apparently based off of NIST CODATA 1998.
 * This is recorded in case of a unit discrepancy
 */

//speed of light
constexpr double C_CGS = 2.99792458e10;
//gravitational constant
constexpr double G_CGS = 6.67430e-8;
//Avogagro's number
constexpr double N_Avogadro = 6.02214076e23;

//Boltzmann constant, Stefan-Boltzmann constant, radiation constant
constexpr double boltzmann_k = 1.380649e-16;
constexpr double boltzmann_sigma = 5.670374419e-5;
constexpr double radiation_a = 4.*boltzmann_sigma/C_CGS;

//Planck constants in CGs units 
constexpr double planck_h_CGS = 6.62607015e-27;
constexpr double planck_hbar_CGS = 1.054571817e-27;

//properties of an electron
constexpr double electron_mass = 9.1093837139e-28;
constexpr double elementary_charge = 4.803204713e-10;
constexpr struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} electron = {electron_mass, elementary_charge, planck_hbar_CGS/(electron_mass*C_CGS)};

//useful properties of a proton
constexpr double proton_mass = 1.67262192595e-24;
constexpr struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} proton = {proton_mass, elementary_charge, planck_hbar_CGS/(proton_mass*C_CGS)};


/** ASTRONOMICAL CONSTANTS
 * Taken from IAU 2015
 * https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf
 * Values in CGS units, unless otherwise noted
 */
constexpr double GMSOLAR = 1.3271244e26;
constexpr double MSOLAR = GMSOLAR / G_CGS;
constexpr double RSOLAR = 6.957e10;
constexpr double LSOLAR = 3.828e33;
//radius of earth, most useful length scale for WD models
constexpr double REARTH = 6.3781e8;

//in astronomical units, things are in terms of
//	base mass = MSOLAR
//  base disance = km (most appropriate for WDs and NSs, definitely not parsec)
//  base speed = km/s
constexpr double C_astro = 2.99792458e5; // km/s
constexpr double G_astro = GMSOLAR * 1.e-15;// km (km/s)^2 /MSOLAR 
//constexpr double G_astro = 4.302e-3;// pc (km/s)^2 /MSOLAR 
//constexpr double G_astro = 1.90809e5;// RSOLAR (km/s)^2/MSOLAR

//this is a fundamental frequency sqrt(GM/R^3) of the Sun, in uHz
//taken from the paper by JCD-DJM and used for comparison to their tables
constexpr double nug = 99.855377; //from JCD, in uHz

#endif
