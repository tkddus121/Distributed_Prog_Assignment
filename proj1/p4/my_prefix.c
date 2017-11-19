#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

#define RANDOM_MAX 1000

void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
    Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
    Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
    Tend.tv_usec += (Tend.tv_sec*1000000);

    printf("Elapsed Time: %lf sec\n", Tend.tv_usec / 1000000.0);
}

int prefix_sum(int rank,int num,int *val)
{
	int arrow, step = 1;
	//Phase 1

//	printf("%d %d %d\n",rank,num,*val);


	for(step = 1; step < num;step *=2)
	{
		//send
		//
		if(rank + step < num)
		{
			MPI_Send(val, 1, MPI_INT, rank + step, 0, MPI_COMM_WORLD);
			//printf("send %d to %d\n",rank,rank+step);
		}
		if( rank > 0 && rank - step >= 0)
		{
			MPI_Status stat;
			int tmp;

			MPI_Recv(&tmp, 1, MPI_INT, rank - step, 0, MPI_COMM_WORLD, &stat);
			*val += tmp;
		}
	}


}


int main( int  argc, char **argv)
{

	int rank, num,size, i;
	int *dataset, data, ori_data;

	struct timeval Tstart, Tend; //time value

	// No MPI call before it. 
	// MPI start
	MPI_Init(&argc, &argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num );

	if(rank == 0)
	{
		printf("Get Ready.\n");
		// start time
		gettimeofday(&Tstart, NULL); 
	}

	srand((rank*3+1)*time(NULL));
    

	data = rand()%RANDOM_MAX;
	ori_data = data;
    prefix_sum(rank,num,&data);


	printf("[Process %d of %d. ori = %d, result = %d\n",rank, num , ori_data,data);  
	
	
	MPI_Barrier(MPI_COMM_WORLD);

	//print time.
	if (rank == 0)
	{
		gettimeofday(&Tend, NULL);
		getElapsedTime(Tstart, Tend);
	}

	MPI_Finalize();
	// No MPI after here

	return 0;

}
