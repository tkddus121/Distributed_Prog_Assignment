#include <stdio.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define ARRAY_SIZE 10000
#define CUDA_CALL(x) { const cudaError_t a = (x); if(a != cudaSuccess) { printf("\nCuda Error: %s (err_num=%d) at line:%d\n", cudaGetErrorString(a), a, __LINE__); cudaDeviceReset(); assert(0);}}
void init_array(int arr[ARRAY_SIZE])
{
	for(int i  = 0 ; i < ARRAY_SIZE ; ++i)
		arr[i] = rand()/10000;

}
static __global__ void cuda_noPath(int *arr, int *output, int size)
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;

	for(int i = 1 ; i < blockDim.x ; i *= 2)
	{
		if (threadIdx.x % (i*2) == 0 && tid + i < size)
			arr[tid] = max(arr[tid], arr[tid + i]);
		__syncthreads();
	}

	if( threadIdx.x == 0)
		output[blockIdx.x] = arr[tid];
}

static __global__ void cuda_Path(int *arr, int *output, int size)
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	int base = blockIdx.x * blockDim.x;

	for(int i = 1 ; i < blockDim.x ; i *= 2)
	{
		int off = threadIdx.x * i*2;
		int idx = base + off;

		if ( off < blockDim.x  &&  idx + i < size)
			arr[tid] = max(arr[idx], arr[idx + i]);
		__syncthreads();
	}

	if( threadIdx.x == 0)
		output[blockIdx.x] = arr[tid];
}


extern "C" void cuda_2(int *res, int arr[ARRAY_SIZE] )
{
	void *output_dev, *arr_dev;

	int thread_num = 256;
	int size = ARRAY_SIZE;
	int block_num = (size + thread_num - 1) / thread_num;


	cudaEvent_t start,stop;
	
    CUDA_CALL(cudaMalloc((void**)&arr_dev, sizeof(int) * ARRAY_SIZE));
    CUDA_CALL(cudaMalloc((void**)&output_dev, sizeof(int) * block_num));

    // transfer data from host to device.
    CUDA_CALL(cudaMemcpy(arr_dev, arr, sizeof(int) * ARRAY_SIZE, cudaMemcpyHostToDevice));

    float dev_time = 0.f;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
    do {

        //CHECK_TIME_START_GPU();
        cuda_Path<<<block_num, thread_num>>>((int*)arr_dev, (int*)output_dev, size);
        //CHECK_TIME_END_GPU(device_time);

        {
            void *tmp = arr_dev;
            arr_dev = output_dev;
            output_dev = tmp;
        }
        size = block_num;
        block_num = (size + thread_num - 1) / thread_num;

    } while (size > 1);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&dev_time, start, stop);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

    printf("[CUDA Path] Elapsed Time : %.5f (sec).\n", dev_time/1000);

    // transfer result from device to host.
    CUDA_CALL(cudaMemcpy(res, arr_dev, sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL( cudaDeviceSynchronize() );

    CUDA_CALL(cudaFree(arr_dev));
    CUDA_CALL(cudaFree(output_dev));


}

extern "C" void cuda_1(int *res, int arr[ARRAY_SIZE] )
{
	void *output_dev, *arr_dev;

	int thread_num = 256;
	int size = ARRAY_SIZE;
	int block_num = (size + thread_num - 1) / thread_num;


	cudaEvent_t start,stop;
	
    CUDA_CALL(cudaMalloc((void**)&arr_dev, sizeof(int) * ARRAY_SIZE));
    CUDA_CALL(cudaMalloc((void**)&output_dev, sizeof(int) * block_num));

    // transfer data from host to device.
    CUDA_CALL(cudaMemcpy(arr_dev, arr, sizeof(int) * ARRAY_SIZE, cudaMemcpyHostToDevice));

    float dev_time = 0.f;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
    do {

        //CHECK_TIME_START_GPU();
        cuda_noPath<<<block_num, thread_num>>>((int*)arr_dev, (int*)output_dev, size);
        //CHECK_TIME_END_GPU(device_time);

        {
            void *tmp = arr_dev;
            arr_dev = output_dev;
            output_dev = tmp;
        }
        size = block_num;
        block_num = (size + thread_num - 1) / thread_num;

    } while (size > 1);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&dev_time, start, stop);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

    printf("[CUDA noPath] Elapsed Time : %.5f (sec).\n", dev_time/1000);

    // transfer result from device to host.
    CUDA_CALL(cudaMemcpy(res, arr_dev, sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CALL( cudaDeviceSynchronize() );

    CUDA_CALL(cudaFree(arr_dev));
    CUDA_CALL(cudaFree(output_dev));


}
int arr[ARRAY_SIZE],
	res_cuda1,
	res_cuda2;

int main(int argc, char *argv[])
{

	init_array(arr);
	cuda_1(&res_cuda1,arr);
	cuda_2(&res_cuda2,arr);
	return 0;

}
