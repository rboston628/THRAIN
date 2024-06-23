#ifndef DUMMYMODESH
#define DUMMYMODESH

//**************************************************************************************
//							The Dummy Modes
//  	These dummy modes are for use in testing.  They mimic some of the features
//      of the stellar pulsation modes, at least to the extent needed for testing.
//**************************************************************************************


#include "../src/MODES/ModeDriver.h"

const double ZEROW2 = 0.16;

//**************************************************************************************
//							DummyMode Class
//  	For mocking behavior of modes for use in tests
//      Mock out all required abstract methods with minimal implementation.
//      Set the static klist.  Created modes will then return their mode order 
//      according to the list, in that order.
//      This is intended for testing the behavior of the mode finder.
//**************************************************************************************
class DummyMode : public ModeBase {
public:
    static const unsigned int num_var=0U;
    void printMode(const char *const c) override;
    void writeMode(const char *const c) override;
    double getRad(std::size_t x) override;
    double getY(int i, std::size_t x) override;
    double getYtilde(int i, std::size_t x) override;
    void modeNumbers(int& x, int& y, int& z) override;
    double SSR() override;
    double tidal_overlap() override;
    double getFreq() override;
    double getPeriod() override;

    // the ones we intend to mock
    int modeOrder() override;
    double getOmega2() override;
    DummyMode();
    DummyMode(int k, int l, int m, ModeDriver* drv);
    DummyMode(double w2, int l, int m, ModeDriver* drv);
    DummyMode(double w1, double w2, int l, int m, ModeDriver* drv);
    static int iter;
    static std::vector<int> klist;
private:
    int K;
};


//**************************************************************************************
//							ControlledMode Class
//  	Like a dummy mode, but we can control the frequency outputs
//**************************************************************************************
class ControlledMode : public DummyMode {
public:
    static const unsigned int num_var=0U;
    // the ones we intend to mock
    int modeOrder() override;
    double getOmega2() override;
    ControlledMode(int k, int l, int m, ModeDriver* drv);
    ControlledMode(double w2, int l, int m, ModeDriver* drv);
    ControlledMode(double w1, double w2, int l, int m, ModeDriver* drv);
    static int iter;
    static std::vector<int> klist;
private:
    int K;
};


//**************************************************************************************
//							DummyModeDriver Class
//  	A mostly empty mode driver, only because dummy modes needs them
//**************************************************************************************
class DummyModeDriver : public ModeDriver {
public:
    static const unsigned int num_var=0U;
    DummyModeDriver(Star* s, double x);
    std::size_t length() override;
    double Gamma1() override;
    double rad(std::size_t x) override;
    std::size_t CentralBC(double **y, double *yo, double s2, int l, int m=0) override;
    std::size_t SurfaceBC(double **y, double *ys, double s2, int l, int m=0) override;
    void getCoeff(double *CC, const std::size_t, const int, const double, const int) override;
    void setupBoundaries() override;
    double SSR(double x, int y, ModeBase* ) override;
    double tidal_overlap(ModeBase*) override;
    double innerproduct(ModeBase*,ModeBase*) override;
    void getBoundaryMatrix(int, double **, int*) override;
	void varnames(std::string*) override;
};


//**************************************************************************************
//							SineMode Class
//  	Returns (sin wx, cos wx), where w = 2pi N, for N
//      This is intended for testing against modes form the SineModeDriver
//**************************************************************************************
class SineMode : public DummyMode {
public:
    SineMode(int N, std::size_t len);
    ~SineMode();
    int modeOrder() override;
    double getOmega2() override ;
    double getFreq() override;
    double getPeriod() override;
    double getRad(std::size_t x) override;
    double getY(int i, std::size_t x) override;
    double getYtilde(int i, std::size_t x) override;

private:
    enum VarNames {sin=0, cos};
    std::size_t len;
    int N;
    double freq;
};

#endif