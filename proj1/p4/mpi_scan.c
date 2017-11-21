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

int main( int  argc, char **argv)
{

	int rank, num,size, i;
	int *dataset, localdata, localrecv;

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

	//generate random value.
    localdata = rand()%RANDOM_MAX;
   // printf("[Process %d]: has data %d\n", rank, localdata);

    //MPI_Scan(sendbuf,recvbuf, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	//
    MPI_Scan(&localdata,&localrecv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //synchronize
    MPI_Barrier(MPI_COMM_WORLD);
   // printf("[Process %d]: has received data: %d \n", rank,localrecv);

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
