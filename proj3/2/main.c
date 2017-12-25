#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#define MATRIX_SIZE 1024

void init_matrix(float mat[MATRIX_SIZE][MATRIX_SIZE])
{
	    for(int i=0; i<MATRIX_SIZE; ++i)
			        for(int j=0; j<MATRIX_SIZE; ++j)
						            mat[i][j] = ( (float) rand() * 2.f / 10000 ) - 1.f;
}


void omp_MatMul(float C[MATRIX_SIZE][MATRIX_SIZE], float A[MATRIX_SIZE][MATRIX_SIZE],float B[MATRIX_SIZE][MATRIX_SIZE])
{


	double start, finish;

	start = omp_get_wtime();
#pragma omp parallel for schedule(dynamic)
	for(int r=0; r<MATRIX_SIZE; ++r) {
		for(int c=0; c<MATRIX_SIZE; ++c) {
			float res = 0.0f;
			for(int k=0; k<MATRIX_SIZE; ++k) {
				res += A[r][k] * B[k][c];
			}
			C[r][c] = res;
		}
	}
	finish = omp_get_wtime();
	printf("[OpenMP] Elapsed Time = %.3f (sec)\n", finish-start);
}

float A[MATRIX_SIZE][MATRIX_SIZE], B[MATRIX_SIZE][MATRIX_SIZE];
float C_openmp[MATRIX_SIZE][MATRIX_SIZE];
int main(int argc, char *argv[])
{
	init_matrix(A);
	init_matrix(B);

	omp_MatMul(C_openmp,A,B);

}

