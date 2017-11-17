#include<stdio.h>
#include<mpi.h>

int main( int  argc, char **argv)
{

	int rank, num, i;

	// No MPI call before it. 
	// MPI start
	MPI_Init(&argc, &argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num );

	printf("Hello from process %i of %i\n",rank,num);

	MPI_Finalize();
	// No MPI after here

	return 0;

}
