// **************************************************************************************
// constants.h
//  	This file describes dependencies, and constants
//  	Many constant taken from https://physics.nist.gov/cuu/Constants/index.html	
//  Reece Boston, Mar 24, 2022
// **************************************************************************************

#include <cmath>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <functional>
#include "../lib/Splinor.h"
#include "../lib/matrix.h"
#include "../lib/string.h"
// #include <time.h> //might use later to put a time stamp

//constants and dependencies
#ifndef CONSTANTS
#define CONSTANTS

//good ol' pi!
const double m_pi = 3.1415926535897932384626433832795;
//speed of light in CGS
const double C_CGS = 2.99792458e10;// in CGS
//gravitational constant in CGS
// previous value G = 6.6725985e-8 taken from: https://iopb.res.in/~somen/constants0.html
// new value taken from: https://physics.nist.gov/cuu/Constants/index.html
const double G_CGS = 6.67430e-8; // in CGS

//Avogagro's number
const double N_Avogadro = 6.02214076e23;

//Boltzmann constant, Stefan-Boltzmann constant, radiation constant
const double boltzmann_k = 1.3806503e-16;// in CGS
const double boltzmann_sigma = 5.670374419e-5;// in CGS
const double radiation_a = 7.5657e-15   ;// in CGS

//this is a fundamental frequency sqrt(GM/R^3) of the Sun, in uHz
//taken from the paper by JCD-DJM and used for comparison to their tables
const double nug = 99.855377; //from JCD, in uHz

//other properties of the sun
const double MSOLAR = 1.989e33;//in CGS
const double RSOLAR = 6.957e10 ;//in CGS
const double LSOLAR = 3.828e33;//in CGS
//radius of earth, most useful length scale for WD models
const double REARTH = 6.378e8 ;//in CGS

//Planck constants in CGs units 
const double planck_h_CGS = 6.62607015e-27;
const double planck_hbar_CGS = 1.054571817e-27;

//properties of an electron
const struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} electron = {9.1095e-28, 4.80320425e-10, planck_hbar_CGS/(9.1095e-28*C_CGS)};

//useful properties of a proton
const struct {
	const double mass_CGS;
	const double charge_CGS;
	const double compton_wavelength_CGS;
} proton = {1.6726e-24, 4.80320425e-10, planck_hbar_CGS/(1.6726e-24*C_CGS)};


//in astronomical units, things are in terms of
//	base mass = MSOLAR
//  base disance = km (most appropriate for WDs and NSs, definitely not parsec)
//  base speed = km/s
const double C_astro = 2.99792458e5; // km/s
const double G_astro = 4.302e-3*30856775812800.0;// km (km/s)^2 /MSOLAR 
//const double G_astro = 4.302e-3;// pc (km/s)^2 /MSOLAR 
//const double G_astro = 1.90809e5;// RSOLAR (km/s)^2/MSOLAR

#endif
