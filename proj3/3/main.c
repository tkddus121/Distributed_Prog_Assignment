#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

#include <unistd.h>
#define ARRAY_SIZE 10000

void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
	Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
	Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
	Tend.tv_usec += (Tend.tv_sec*1000000);

	printf("Elapsed Time: %.5lf (sec)\n", Tend.tv_usec / 1000000.0);
}



void init_array(int arr[ARRAY_SIZE])
{
	for(int i  = 0 ; i < ARRAY_SIZE ; ++i)
		arr[i] = rand()/10000;

}

void seq_find_min(int *res, int arr[ARRAY_SIZE])
{
	struct timeval Tstart, Tend;


	gettimeofday(&Tstart,NULL);	

	*res = arr[0];
	for(int i = 1 ; i < ARRAY_SIZE ; ++i)
	{
		if( *res < arr[i])
			*res = arr[i];
	}

	gettimeofday(&Tend,NULL);

	printf("[Sequential] ");
	getElapsedTime(Tstart,Tend);

}

int a[ARRAY_SIZE];
int result_seq,result_cuda1;

int main( int argc , char *argv[])
{
	init_array(a);

	seq_find_min(&result_seq,a);



}
