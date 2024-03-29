// **************************************************************************************
// matrix.h
//		This file describes certain matrix operations for the GRPulse code.  
// **************************************************************************************

#include <cmath>
#include <stdio.h>
#include <assert.h>

#ifndef MATRIX
#define MATRIX

namespace matrix {

template <std::size_t N, std::size_t M, class T>
void swap_rows(T (&m)[N][M], std::size_t i, std::size_t j){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	T temp;
	for(std::size_t k=0; k<M; k++){
		temp = m[i][k];
		m[i][k] = m[j][k];
		m[j][k] = temp;
	}
}

template <std::size_t N, class T>
void swap_rows(T (&m)[N], std::size_t i, std::size_t j){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	T temp = m[i];
	m[i] = m[j];
	m[j] = temp;
}

template <std::size_t N, std::size_t M, class T>
void add_rows(T (&m)[N][M], std::size_t i, std::size_t j, T coeff=1.0){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	for(std::size_t k=0; k<M; k++){
		m[i][k] = m[i][k] + coeff*m[j][k];
	}
}

template <std::size_t N, class T>
void add_rows(T (&m)[N], std::size_t i, std::size_t j, T coeff=1.0){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
		m[i] = m[i] + coeff*m[j];
}

template <std::size_t N, std::size_t M, class T>
void print_matrix(const T (&m)[N][M], int k=0){
	printf("MATRIX %d\n", k);
	for(std::size_t i=0; i<N; i++){
		printf("\t[");
		for(std::size_t j=0; j<M; j++){
			printf("\t%lf,", m[i][j]);
		}
		printf("]\n");
	}
}

//calculate the determinant of an NxN matrix
//uses LU decomposition, based on GSL function
//will destroy the matrix m
template <std::size_t N, class T>
T determinant(T (&m)[N][N]){
	T det=1.0, ajj=0.0, coeff=0.0;
	double temp = 0.0, big = 0.0;
	std::size_t i=0,j=0,k=0,ipiv=0;
	for(j=0; j<N-1;j++){
		//find the max in this column
		big = 0.0;
		for(i=j;i<N;i++){
			temp = std::abs(m[i][j]);
			if(temp>big){
				big=temp;
				ipiv=i;
			}
		}
		//if max element is 0, then singular matrix
		if(big==0.0){
			det = 0.0;
			break;
		};
		
		//if max is not on diagonal, swap rows
		if(ipiv != j){
			matrix::swap_rows(m, ipiv, j);
			//swapping rows swaps sign of det
			det = -det;
		}
		ajj = m[j][j];
		if(ajj!=0.0){
			// for each row below row j
			for(i=j+1;i<N;i++){
				coeff = m[i][j]/ajj;	
				for(k=j+1;k<N;k++){ // we do not need to fully diagonalize
					m[i][k] -= coeff*m[j][k];
				}
			}	
		}
		//if the diagonal element is zero, this matrix is singular
		else{
			det = 0.0;
			break;
		}
		det *= m[j][j];
	}
	det *= m[N-1][N-1];
	return det;
}

//in equation Ax=b, invert NxN matrix A to solve for x
//will destroy the matrix A
//can handle equations Ax=0
template <std::size_t N, class T>
int invertMatrix(T (&m)[N][N], T (&b)[N], T (&x)[N]){
	//a flag, in case we are soving homoegenous problem
	bool HOMOGENEOUS = true;
	T dummy = 0.0;
	for(std::size_t j=0; j<N; j++) HOMOGENEOUS &= (b[j]==0.0);

	T det=1.0, ajj=0.0, coeff=0.0;
	double temp = 0.0, big = 0.0;
	std::size_t i=0,j=0,k=0,ipiv=0;
	for(j=0; j<N;j++){
		//find the max in this column
		big = 0.0;
		for(i=j;i<N;i++){
			temp = std::abs(m[i][j]);
			if(temp>big){
				big=temp;
				ipiv=i;
			}
		}		
		//if max is not on diagonal, swap rows
		if(ipiv != j){
			matrix::swap_rows(m, ipiv, j);
			matrix::swap_rows(b, ipiv, j);
		}
		ajj = m[j][j];
		if(ajj!=0.0){
			// for each row below row j
			for(i=j+1;i<N;i++){
				coeff = -m[i][j]/ajj;	
				add_rows(m, i, j, coeff);	
				add_rows(b, i, j, coeff);		
				m[i][j] = 0.0;
			}	
		}
		T a = 1.0/m[j][j];
		for(int i=0; i<N; i++){
			m[j][i] *= a;
		}
		b[j] *= a;
	}
	//resubstitute
	std::size_t L=0;
	x[N-1] = b[N-1];
	if(HOMOGENEOUS) x[N-1] = 1.0;
	for(int i=N-2;i>=0; i--){
		L=0;
		while(m[i][L]==T(0) && L<N){	
			L++;
		}
		x[i] = (HOMOGENEOUS ? 0.0 : b[i]);
		for(std::size_t j=L+1; j<N; j++){
			x[i] -= m[i][j]*x[j];
		}
	}
	bool produced_nan = false;
	for(int i=0; i<N; i++){
		if( (x[i]-x[i]) != T(0) ) produced_nan = true;
	}
	if(produced_nan){printf("ERROR in invertMatrix: NaNs produced\n"); return 1;}
	return 0;
}

// will destroy the diagonal vectors
// this is needed by classes where N cannot be known at compile time
template <class T>
int invertTridiagonal(T *left, T *diag, T *rite, T *RHS, T *x, std::size_t const N){
	double w;
	
	// sweep foward changing coefficients
	rite[0] /= diag[0];
	RHS[0] /= diag[0];
	diag[0] = 1.0;
	for(int k=1; k<N; k++){
		w = left[k]/diag[k-1];
		diag[k] -= w * rite[k-1];
		RHS[k]  -= w * RHS[k-1];
	}

	// back-substitute
	x[N-1] = RHS[N-1]/diag[N-1];
	for(int k=(N-2); k>=0; k--){
		x[k] = (RHS[k] - rite[k]*x[k+1])/diag[k];
	}
	return 0;
}

} // namespace matrix
#endif