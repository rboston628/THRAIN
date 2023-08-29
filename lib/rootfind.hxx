#ifndef ROOTFINDHXX
#define ROOTFINDHXX

#include "matrix.h"

// *************************************************************************************
//					NEWTON METHODS
//  These functions implement Newton's method in one or many dimensions
// *************************************************************************************

//single-variable search, searching for ZERO of func
//intended to work with real-valued or complex-valued numbers
template <typename T>
T rootfind::newton_search(
	std::function<T(T)>& func,	//the function to find zero of
	T &x,						//an initial guess for the zero
	T dx,						//the step to use in numerical derivatives
	double const tol,				//tolerance of search, to be this accurate
	std::size_t const max_iter	//maximum number of iterations in search
){
	T f1 = func(x), f2, ddx;
	std::size_t tracker = 1;
	
	while(std::abs(f1)>tol && tracker++!=max_iter){
		f2 = func(x+dx);
		ddx = -dx*f1/(f2-f1);
		x += ddx;
		dx = ddx;
		f1 = func(x);
	}
	return f1;
}

//single-variable search, searching for func = target
//intended to work with real-valued or complex-valued numbers
template <typename T>
T rootfind::newton_search(
	std::function<T(T)>& func,	//the function to match to target
	T const target,				//the target to be matched to
	T &x,						//an initial guess for x
	T dx,						//the step to use in numerical derivatives
	double const tol,			//tolerance of search, to be this accurate
	std::size_t const max_iter	//maximum number of iterations in search
){
	T f1 = func(x), f2, ddx;
	std::size_t tracker = 1;

	while(std::abs(target-f1)>tol && tracker++!=max_iter){
		f2 = func(x+dx);
		ddx = dx*(target-f1)/(f2-f1);
		x += ddx;
		dx = ddx;
		f1 = func(x);
	}
	return f1;
}


template <size_t np>
void rootfind::newton_search(
	std::function<void(double f[np],double x[np])>& func,	//f=vector function to zero, x=input array
	double (&x1)[np],			//an initial guess for x
	double (&dx)[np],			//the step to use in numerical derivatives
	double const tol,			//tolerance of search, to be this accurate
	std::size_t const max_iter,	//maximum number of iterations in search
	std::function<bool(double[np])>& var_limit  //a function limiting values of x1
){
	double x2[np], ddx[np], f1[np], f2[np], dfdx[np][np], dxsave[np];
	func(f1,x1);
	std::size_t tracker = 1;
	double F1=0.0, F2=0.0, maxDF = -1.0;
	for(int i=0; i<np; i++) F1 += 0.5*f1[i]*f1[i];
	
	for(int i=0; i<np; i++) ddx[i] = dx[i];
	for(int i=0; i<np; i++) maxDF = fmax(maxDF, std::abs(f1[i]));
	while( maxDF > tol && tracker++!=max_iter){
		//calculate the matrix of derivatives (df/dx)
		for(int i=0; i<np; i++) x2[i] = x1[i];
		for(int i=0; i<np; i++) {
			x2[i] = x1[i] + ddx[i];
			func(f2, x2);
			for(int j=0; j<np; j++){
				dfdx[j][i] = (f2[j]-f1[j])/(ddx[i]);
			}
			x2[i] = x1[i];
		}
		//solve the equation f1 = -(dfdx)*(dx) for dx
		for(int i=0; i<np; i++) f2[i] = -f1[i];
		for(int i=0; i<np; i++) dxsave[i] = dx[i];
		if(matrix::invertMatrix(dfdx, f2, dx)){
			//if the matrix is singular or otherwise fails
			// then do something to try to recover
			printf("ERROR: Matrix inversion failed!\nTrying to recover with past value.\n");
			//use the last gradient
			double scale = rootfind::pseudo_unif();
			for(int i=0; i<np; i++) dx[i] = scale*dxsave[i];
		}
		for(int i=0; i<np; i++) x2[i] = x1[i] + dx[i];
		//if the change causes variables to violate some limit, adjust
		// an example would be if some x must be non-negative
		while(!var_limit(x2)){
			for(int i=0; i<np; i++) dx[i] *= 0.1;
			for(int i=0; i<np; i++) x2[i] = x1[i] + dx[i];
		}
		//if the increase is too large, throttle it
		// we check this by making sure the absolute value always decreases
		func(f2, x2);
		
		F2 = 0.0;
		for(int i=0; i<np; i++) F2 += 0.5*f2[i]*f2[i];
		double L = 1.0;
		while(F2>F1){
			if(L<1.e-6) {L = 1.e-6; break;}
			L *=0.1;
			for(int i=0; i<np; i++) x2[i] = x1[i] + L*dx[i];
			func(f2, x2);
			F2 = 0.0;
			for(int i=0; i<np; i++) F2 += 0.5*f2[i]*f2[i];
		}
		F1 = F2;
		for(int i=0; i<np; i++) f1[i] = f2[i];
		for(int i=0; i<np; i++) x1[i] = x2[i];
		
		//if none of the dx are zero, use these as the new differences in numerical differentiation
		bool anyZero = false;
		for(int i=0; i<np; i++) anyZero |= (dx[i]==0.0);
		if(!anyZero){
			for(int i=0; i<np; i++) ddx[i] = dx[i];
		}
		
		//determine the max difference
		maxDF = -1.0;
		for(int i=0; i<np; i++) if(fabs(f1[i])>maxDF) maxDF = fabs(f1[i]);
	}
}

template <size_t np>
void rootfind::newton_search(
	std::function<void(double f[np],double x[np])>& func,	//f=vector function, x=input array
	double (&target)[np], 		//the target to be matched to
	double (&x1)[np], 			//an initial guess for x
	double (&dx)[np],			//the step to use in numerical derivatives
	double const tol,			//tolerance of search, to be this accurate
	std::size_t const max_iter,	//maximum number of iterations in search
	std::function<bool(double[np])>& var_limit //a function limiting values of x1
){
	double x2[np], ddx[np], f1[np], f2[np], dfdx[np][np], dxsave[np];
	func(f1,x1);
	std::size_t tracker = 1;
	double F1=0.0, F2=0.0, maxDF = 1.0;
	for(int i=0; i<np; i++) F1 += 0.5*pow((target[i]-f1[i]),2);

	for(int i=0; i<np; i++) ddx[i] = dx[i];
	while( maxDF > tol && tracker++!=max_iter){
		//calculate the matrix of derivatives (df/dx)
		for(int i=0; i<np; i++) x2[i] = x1[i];
		for(int i=0; i<np; i++) {
			x2[i] = x1[i] + dx[i];
			func(f2, x2);
			for(int j=0; j<np; j++){
				dfdx[j][i] = (f2[j]-f1[j])/(dx[i]);
			}
			x2[i] = x1[i];
		}		
		//solve the equation (target-f1) = (dfdx)*(dx) for dx
		for(int i=0; i<np; i++) f2[i] = target[i]-f1[i];
		for(int i=0; i<np; i++) dxsave[i] = dx[i];
		if(matrix::invertMatrix(dfdx, f2, dx)){
			//if the matrix is singular or otherwise fails
			// then do something to try to recover
			printf("ERROR: Matrix inversion failed!\nTrying to recover with past value.\n");
			//use a random multiple of the  last gradient
			double scale = rootfind::pseudo_unif();
			for(int i=0; i<np; i++) dx[i] = scale*dxsave[i];
		}
		for(int i=0; i<np; i++) x2[i] = x1[i] + dx[i];
		//if the change causes variables to violate some limit, adjust
		// an example would be if some x must be non-negative
		int stop=0;
		while(!var_limit(x2) && stop<4){
			for(int i=0; i<np; i++) dx[i] *= 0.1;
			for(int i=0; i<np; i++) x2[i] = x1[i] + dx[i];
			stop++;
		}
		//if the increase is too large, throttle it
		func(f2, x2);
		F2 = 0.0;
		for(int i=0; i<np; i++) F2 += 0.5*pow((target[i]-f2[i]),2);
		double L=1.0;
		while(F2>F1){
			if(L<1.e-3) {L = 1.e-3; break;}
			L *=0.1;
			for(int i=0; i<np; i++) x2[i] = x1[i] + L*dx[i];
			func(f2, x2);
			F2 = 0.0;
			for(int i=0; i<np; i++) F2 += 0.5*pow((target[i]-f2[i]),2);
		}
		F1 = F2;
		for(int i=0; i<np; i++) f1[i] = f2[i];
		for(int i=0; i<np; i++) x1[i] = x2[i];
		
		//if none of the dx are zero, use these as the new differences in numerical differentiation
		bool anyZero = false;
		for(int i=0; i<np; i++) anyZero |= (dx[i]==0.0);
		if(!anyZero){
			for(int i=0; i<np; i++) ddx[i] = dx[i];
		}
		
		maxDF = -1.0;
		for(int i=0; i<np; i++) if(fabs((target[i]-f1[i])) > maxDF) maxDF = fabs((target[i]-f1[i]));
	}
}
#endif