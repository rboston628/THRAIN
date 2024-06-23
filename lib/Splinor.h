//A spline interpolation object.  Pass support points, creates an interpolation
//Afterwards, can be passed x-values and will return interpolate y-values

#ifndef INTERP
#define INTERP

#include <cmath>
#include <tuple>

class Interpolator {
	public:
	virtual double interp(double const) = 0;
	virtual double operator()(double const) = 0;
	virtual double deriv(double const) = 0;
	virtual ~Interpolator(){};
};


class Splinor : public Interpolator {
	public:
	enum class bc : char {NATURAL, CLAMPED, QUADRATIC};
	Splinor(
		double const *const x, 
		double const *const y, 
		int const L,
		bc const bc_type = bc::NATURAL,
		double const yprime0 = 0.0, // derivative at xa -- default natural spline
		double const yprimeN = 0.0  // derivative at xb -- default natural spline
	);
	~Splinor();
	double interp(double const) override;
	double operator()(double const) override;
	double deriv(double const) override;
	std::tuple<double,double,double,double> getCoefficients(std::size_t const);
	
	private:
	std::size_t len, xlast;
	bc const bc_type;
	double xa, xb; //boundaries
	double yprimea, yprimeb; // derivatives at boundaries
	std::size_t findPosition(double const);
	void calculateSplineCoefficients(double const *const, double const *const, bc const);
	void makeNaturalSpline(double const *const, double const *const);
	void makeClampedSpline(double const *const, double const *const);
	void makeQuadraticSpline(double const *const, double const *const);
	// TODO make a not-a-knot spline
	double *S;
	double *xarr;
	double *yarr;
};

#endif