#ifndef DUMMYMODESC
#define DUMMYMODESC

#include "DummyMode.h"

//**************************************************************************************
//							DummyMode Class
//  	For mocking behavior of modes for use in tests
//**************************************************************************************
// needed for abstract class
void DummyMode::printMode(const char *const c){}
void DummyMode::writeMode(const char *const c){}
double DummyMode::getRad(std::size_t x){return 0.0;}
double DummyMode::getY(int i, std::size_t x){return 0.0;}
double DummyMode::getYtilde(int i, std::size_t x){return 0.0;}

double DummyMode::SSR(){return 0.0;}
double DummyMode::tidal_overlap(){return 0.0;}
double DummyMode::getFreq() {return 0.0;}
double DummyMode::getPeriod() {return 0.0;}

// the ones we intend to mock
int DummyMode::modeOrder(){return K;}
void DummyMode::modeNumbers(int& x, int& y, int& z){x=modeOrder(); y=0; z=0;}
double DummyMode::getOmega2(){return (K==0? ZEROW2 : double(K));}
DummyMode::DummyMode() {}
DummyMode::DummyMode(int k, int l, int m, ModeDriver* drv) {
    K = klist[iter];
    iter = (iter+1)%klist.size();
}
DummyMode::DummyMode(double w2, int l, int m, ModeDriver* drv) {
    K = int(floor(w2));
}
DummyMode::DummyMode(double w1, double w2, int l, int m, ModeDriver* drv) {
    K = int((w2+w1)/2);
}

int DummyMode::iter(0);
std::vector<int> DummyMode::klist({});

//**************************************************************************************
//							ControlledMode Class
//  	Like a dummy mode, but we can control the frequency outputs
//**************************************************************************************
int ControlledMode::modeOrder(){return K;}
double ControlledMode::getOmega2(){return (K==0? ZEROW2 : double(K));}
ControlledMode::ControlledMode(int k, int l, int m, ModeDriver* drv) {
    K = klist[iter];
    iter = (iter+1)%klist.size();
}
ControlledMode::ControlledMode(double w2, int l, int m, ModeDriver* drv) {
    K = klist[iter];
    iter = (iter+1)%klist.size();
}
ControlledMode::ControlledMode(double w1, double w2, int l, int m, ModeDriver* drv) {
    K = klist[iter];
    iter = (iter+1)%klist.size();
}

int ControlledMode::iter(0);
std::vector<int> ControlledMode::klist({});

//**************************************************************************************
//							DummyModeDriver Class
//  	A mostly empty mode driver, only because dummy modes needs them
//**************************************************************************************
DummyModeDriver::DummyModeDriver(Star* s, double x) : ModeDriver(num_var, s) {}
std::size_t DummyModeDriver::length() { return 0; }
double DummyModeDriver::Gamma1() { return 0; }
double DummyModeDriver::rad(std::size_t x) { return 0; }
std::size_t DummyModeDriver::CentralBC(double **y, double *yo, double s2, int l, int m) {return 0;}
std::size_t DummyModeDriver::SurfaceBC(double **y, double *ys, double s2, int l, int m) {return 0;}
void DummyModeDriver::getCoeff(double *CC, const std::size_t, const int, const double, const int) {}
void DummyModeDriver::setupBoundaries(){}
double DummyModeDriver::SSR(double x, int y, ModeBase* ) {return 0.0;}
double DummyModeDriver::tidal_overlap(ModeBase*) {return 0.0;}
double DummyModeDriver::innerproduct(ModeBase*,ModeBase*) {return 0.0;}
void DummyModeDriver::getBoundaryMatrix(int, double **, int*){}
void DummyModeDriver::varnames(std::string*){}	//names of variables to print out


//**************************************************************************************
//							SineMode Class
//  	Returns (sin wx, cos wx), where w = 2pi N, for N
//      Intended for testing, comparing to results of SineModeDriver
//**************************************************************************************
SineMode::SineMode(int N, std::size_t len): N(N), len(len), freq(2. * M_PI * N) {}
SineMode::~SineMode(){}
// to match results of SindeModeDriver, define Osaki-Scuflaire mode order = -2N
int SineMode::modeOrder(){return -2 * N;}
// relate w^2, freq, and period
double SineMode::getOmega2(){return freq * freq;}
double SineMode::getFreq(){return freq;}
double SineMode::getPeriod(){return 1./double(N);}
// radius defined on [0,1]
double SineMode::getRad(std::size_t x){return double(x)/double(len-1);}
// result is either sine or cosine
double SineMode::getY(int i, std::size_t x){
    double y = 0.0;
    switch(i){
    case sin:
        y = std::sin( freq * getRad(x) );
        break;
    case cos:
        y = std::cos( freq * getRad(x) );
        break;
    }
    return y;
}
double SineMode::getYtilde(int i, std::size_t x){return getY(i,x);}

#endif