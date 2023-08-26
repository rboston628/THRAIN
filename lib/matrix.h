// **************************************************************************************
// matrix.h
//		This file describes certain matri operations for the GRPulse code.  
// **************************************************************************************


#ifndef MATRIX
#define MATRIX

#include <stdio.h>

namespace matrix{

template <size_t N, size_t M, class T>
void swap_rows(T (&m)[N][M], size_t i, size_t j){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	T temp;
	for(size_t k=0; k<M; k++){
		temp = m[i][k];
		m[i][k] = m[j][k];
		m[j][k] = temp;
	}
}

template <size_t N, class T>
void swap_rows(T (&m)[N], size_t i, size_t j){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	T temp = m[i];
	m[i] = m[j];
	m[j] = temp;
}

template <size_t N, size_t M, class T>
void add_rows(T (&m)[N][M], size_t i, size_t j, T coeff=1.0){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
	for(size_t k=0; k<M; k++){
		m[i][k] = m[i][k] + coeff*m[j][k];
	}
}

template <size_t N, class T>
void add_rows(T (&m)[N], size_t i, size_t j, T coeff=1.0){
	assert(0<=i && i<N);
	assert(0<=j && j<N);
		m[i] = m[i] + coeff*m[j];
}

template <size_t N, size_t M, class T>
void print_matrix(const T (&m)[N][M], int k=0){
	printf("MATRIX %d\n", k);
	for(size_t i=0; i<N; i++){
		printf("\t[");
		for(size_t j=0; j<M; j++){
			printf("\t%lf,", m[i][j]);
		}
		printf("]\n");
	}
}

//calculate the determinant of an NxN matrix
//uses LU decomposition, based on GSL function
//will destroy the matrix m
template <size_t N, class T>
T determinant(T (&m)[N][N]){
	C magnitude;
	T det=1.0, ajj=0.0, coeff=0.0;
	double temp = 0.0, big = 0.0;
	size_t i=0,j=0,k=0,ipiv=0;
	for(j=0; j<N-1;j++){
		//find the max in this column
		big = 0.0;
		for(i=j;i<N;i++){
			temp = std::fabs(m[i][j]);
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
template <size_t N, class T>
int invertMatrix(T (&m)[N][N], T (&b)[N], T (&x)[N]){
	//a flag, in case we are soving homoegenous problem
	bool HOMOGENEOUS = true;
	T dummy = 0.0;
	for(size_t j=0; j<N; j++) HOMOGENEOUS &= (b[j]==0.0);

	C bigger;
	T det=1.0, ajj=0.0, coeff=0.0;
	double temp = 0.0, big = 0.0;
	size_t i=0,j=0,k=0,ipiv=0;
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
	size_t L=0;
	x[N-1] = b[N-1];
	if(HOMOGENEOUS) x[N-1] = 1.0;
	for(int i=N-1;i>=0; i--){
		L=0;
		while(m[i][L]==T(0) && L<N-1){	
			L++;
		}
		x[i] = (HOMOGENEOUS ? 0.0 : b[i]);
		for(size_t j=L+1; j<N; j++){
			x[i] -= m[i][j]*x[j];
		}
	}
	return 0;
}

}

#endif