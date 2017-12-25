#include <cuda_runtime.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <device_launch_parameters.h>


#define CUDA_CALL(x) { const cudaError_t a = (x); if(a != cudaSuccess) { printf("\nCuda Error: %s (err_num=%d) at line:%d\n", cudaGetErrorString(a), a, __LINE__); cudaDeviceReset(); assert(0);}}

#define MATRIX_SIZE 1024

int BLOCK_WIDTH,BLOCK_HEIGHT;
void init_mat(float mat[MATRIX_SIZE][MATRIX_SIZE]);
static __global__ void cuda_1(float C[MATRIX_SIZE][MATRIX_SIZE],float A[MATRIX_SIZE][MATRIX_SIZE] , float B[MATRIX_SIZE][MATRIX_SIZE]);

extern "C" void cuda_version_1(float C[MATRIX_SIZE][MATRIX_SIZE],float A[MATRIX_SIZE][MATRIX_SIZE] , float B[MATRIX_SIZE][MATRIX_SIZE])
{

	void *A_dev, *B_dev, *C_dev;

	cudaEvent_t start,stop;
	

	CUDA_CALL(cudaMalloc((void**)&A_dev, sizeof(float) * MATRIX_SIZE * MATRIX_SIZE));
	CUDA_CALL(cudaMalloc((void**)&B_dev, sizeof(float) * MATRIX_SIZE * MATRIX_SIZE));
	CUDA_CALL(cudaMalloc((void**)&C_dev, sizeof(float) * MATRIX_SIZE * MATRIX_SIZE));

	// take data from host to device.
	CUDA_CALL(cudaMemcpy(A_dev, A, sizeof(float) * MATRIX_SIZE * MATRIX_SIZE,cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMemcpy(A_dev, B, sizeof(float) * MATRIX_SIZE * MATRIX_SIZE,cudaMemcpyHostToDevice));

	dim3 block(BLOCK_WIDTH, BLOCK_HEIGHT);
	dim3 grid(MATRIX_SIZE/BLOCK_WIDTH, MATRIX_SIZE/BLOCK_HEIGHT);
	float dev_time;



	printf("[CUDA 1] Start Launching Kernel.\n");
	//CHECK_TIME_START_GPU();

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	cuda_1<<<grid, block>>>((float(*)[MATRIX_SIZE])C_dev, (float(*)[MATRIX_SIZE])A_dev, (float(*)[MATRIX_SIZE])B_dev);

	//CHECK_TIME_END_GPU(dev_time);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&dev_time, start, stop);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("[CUDA 1] Processing Elapsed Time : %.3f (sec).\n",dev_time/1000);




	//take result from device to host.
	CUDA_CALL(cudaMemcpy( C, C_dev, sizeof(float) * MATRIX_SIZE* MATRIX_SIZE, cudaMemcpyDeviceToHost));
	CUDA_CALL( cudaDeviceSynchronize() );


}
static __global__ void cuda_1(float C[MATRIX_SIZE][MATRIX_SIZE],float A[MATRIX_SIZE][MATRIX_SIZE] , float B[MATRIX_SIZE][MATRIX_SIZE])
{
	int col = blockIdx.x * blockDim.x + threadIdx.x;
	int row = blockIdx.y * blockDim.y + threadIdx.y;
	float result = 0.0f;

	for(int i = 0 ; i < MATRIX_SIZE; ++i)
	{
		result += A[row][i] * B[i][col];
	}
	C[row][col] = result;


}

void init_mat(float mat[MATRIX_SIZE][MATRIX_SIZE])
{
	for (int i = 0 ; i < MATRIX_SIZE; ++i)
		for(int j = 0 ; j < MATRIX_SIZE; ++j)
			mat[i][j] = ( (float) rand()*2.f/RAND_MAX ) - 1.f;

}

float A[MATRIX_SIZE][MATRIX_SIZE],B[MATRIX_SIZE][MATRIX_SIZE],C_cuda_1[MATRIX_SIZE][MATRIX_SIZE];

int main(int argc, char *argv[])
{

	printf("input BLOCK WIDTH , BLOCK HEIGHT : ");
	scanf("%d %d",&BLOCK_WIDTH,&BLOCK_HEIGHT);
	init_mat(A);
	init_mat(B);

	cuda_version_1(C_cuda_1, A, B);


}
