#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define NUM_OF_ELEM 1024
#define MAX_VALUE 100

void sort(int* , int , MPI_Comm);

void Set_Sending_Size(int , int , int* , int* , int* , int* );
void Set_Receiving_Size(int , int , int* , int* , int* , int* );

void sequential_quick_sort(int* ,int ,int );

int set_chunk_size(int elements, int proc, int rank){
	return elements / proc + ((rank < elements % proc) ? 1 : 0);
}


void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
    Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
    Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
    Tend.tv_usec += (Tend.tv_sec*1000000);

    printf("Elapsed Time: %lf sec\n", Tend.tv_usec / 1000000.0);
}

////////////////////////////////////////
//
//

void quick_sort(int * begin, int local_size, MPI_Comm comm){

	// this is a commom size which is less than root size in the tree of processors
    int comm_size; 
    int rank; 
    MPI_Comm_size(comm, &comm_size); 
    MPI_Comm_rank(comm, &rank); 


    int bigger_size = comm_size + 2 * local_size; 
    int* internal_array = (int*)malloc(sizeof(int) * bigger_size);

	//printf("current size is %d common size is %d\n", local_size , comm_size);
	
    // use a new array to store the initial size of array on each processors
    int* initial_sizes = (int*)malloc(sizeof(int) * comm_size);
    MPI_Allgather(&local_size, 1, MPI_INT, &initial_sizes[0], 1, MPI_INT, comm); 
	int i=0; 
    for( i = 0; i < local_size; i++)
	{
        internal_array[i] = begin[i]; 
    }

	// call the sort() function, the sort function will call itself recursively
	sort(internal_array, local_size, comm);

	// sorted_sizes is the the local_size of all processors
    int* sorted_sizes = (int*)malloc(sizeof(int) * comm_size);
    MPI_Allgather(&local_size, 1, MPI_INT, &sorted_sizes[0], 1, MPI_INT, comm); 

    int* sendcounts = (int*)malloc(sizeof(int) * comm_size);
    int* senddispls = (int*)malloc(sizeof(int) * comm_size);
    int* recvcounts = (int*)malloc(sizeof(int) * comm_size);
    int* recvdispls = (int*)malloc(sizeof(int) * comm_size); 

    // calculate sendcounts, senddispls, recvcounts, recvdispls to send it to MPI_Alltoallv to perform merging
    Set_Sending_Size(comm_size, rank, sorted_sizes, initial_sizes, sendcounts, senddispls ); 
	Set_Receiving_Size(comm_size, rank, sorted_sizes, initial_sizes, recvcounts, recvdispls); 
    MPI_Alltoallv(internal_array, sendcounts, senddispls, MPI_INT, begin, recvcounts, recvdispls, MPI_INT, comm); 

	// free all allocated temporary arrays    
	MPI_Barrier(comm); 
    free(internal_array); free(initial_sizes); free(sorted_sizes); 
    free(sendcounts);	  free(senddispls);
    free(recvcounts);	  free(recvdispls); 
	return ;
}

void sort(int* begin, int local_size, MPI_Comm comm)
{
	int comm_size,rank;

	MPI_Comm_size(comm,&comm_size);
	MPI_Comm_rank(comm,&rank);


	if( comm_size <= 1 )
	{
		//if processor is 1, we perform  bubble sort.\
		int c = 0 ;
		int d = 0 ;
		int swap = 0 ;
		for (int c = 0 ; c < ( local_size - 1 ); c++)
		{
			for (d = 0 ; d < local_size - c - 1; d++)
			{
				if (begin[d] > begin[d+1]) 
				{
					swap       = begin[d];
					begin[d]   = begin[d+1];
					begin[d+1] = swap;
				}
			}
		}
	}
	else
	{
		int tmp = 0,i = 0;
		//find root processor in tree.
		int root_size = 0;
		MPI_Allreduce(&local_size, &root_size, 1, MPI_INT, MPI_SUM, comm);

		//create random pivot.

		int pivot_idx = 0;
		if( rank == 0 )
		{
			pivot_idx = rand()% root_size;
		}
	
		MPI_Bcast(&pivot_idx, 1, MPI_INT, 0 , comm);


		//find pivot rank
		//
		int pivot_holder ; 
		int x = root_size / comm_size ;
		int y = root_size % comm_size ;
		if( pivot_idx < (x + 1) * (y) )
		{
			pivot_holder = pivot_idx/(x + 1); 
		} 
		else
		{
			pivot_holder = ( (pivot_idx - (x + 1) * (y) )/(x) + (y) ); 
		}
		////////////////


		// the pivot_holder processor broadcast the value of pivot to all
		int pivot = 0; 
		if(rank == pivot_holder){
			int z ;			
			if( pivot_idx < (x + 1) * (y) )
			{
				z =  ( pivot_idx % (x + 1) ); 
			} else
			{
				z = ( (pivot_idx - (x + 1) * (y) ) % (x) ); 
			}			
			pivot = *(begin + z); 
		}
		MPI_Bcast(&pivot, 1, MPI_INT, pivot_holder, comm); 

		// now the array will split to two subarrays
		// the offset of splitting is calculated here
		// split by pivot.
		int off = 0; 
		for( i = 0; i < local_size; i++)
		{
			if( *(begin + i) <= pivot )
			{

                tmp = *(begin + off); 
                *(begin + off) = *(begin + i); 
                *(begin + i) = tmp; 
				off += 1; 
			}
		}
//////////////////////
	int subarray_size[2] = {0, 0}; 
		subarray_size[0] = off;
		subarray_size[1] = local_size - off; 

		// save the left child of the array and right child of array after splitting by pivot
        int* left_child = begin; 
        int* right_child = begin + off; 
       

         int* left_child_sizes = (int*)malloc(sizeof(int) * comm_size);
         int* right_child_sizes = (int*)malloc(sizeof(int) * comm_size);
         
         // gather the number of left_child_sizes and right_child_sizes together
         int* all_subarray_sizes = (int*)malloc(sizeof(int) * 2*comm_size);
         MPI_Allgather(subarray_size, 2, MPI_INT, all_subarray_sizes, 2, MPI_INT, comm); 
         for( i = 0; i < comm_size; i++){
            left_child_sizes[i] = all_subarray_sizes[2*i]; 
            right_child_sizes[i] = all_subarray_sizes[2*i + 1]; 
         }
         
		// calculate the number of elements in left side of array and right side of array in total 
        int all_left_sizes = 0; 
        int all_right_sizes = 0;
        for( i = 0; i < comm_size; i++){
            all_left_sizes += left_child_sizes[i]; 
            all_right_sizes += right_child_sizes[i]; 
        }


        int my_left_child_size = (all_left_sizes * comm_size)/root_size; 
        my_left_child_size = ( (my_left_child_size > 0) ? my_left_child_size : 1 ); 
        my_left_child_size = ( (my_left_child_size < comm_size) ? my_left_child_size : (my_left_child_size - 1) ); 
        int my_right_child_size = comm_size - my_left_child_size;  

		// calsulate the final size of array for each processor
        int* final_sizes = (int*)malloc(sizeof(int) * comm_size);
        for( i = 0; i < my_left_child_size; i++){
            final_sizes[i] = set_chunk_size(all_left_sizes, my_left_child_size, i); 
        }
        for( i = my_left_child_size; i < comm_size; i++){
            final_sizes[i] = 0; 
        }

        
        int new_size; 
        if(rank < my_left_child_size){
            new_size = set_chunk_size(all_left_sizes, my_left_child_size, rank); 
        } else{
            new_size = set_chunk_size(all_right_sizes, my_right_child_size, rank - my_left_child_size); 

        }

        int* new_array = (int*)malloc(sizeof(int) * new_size);
        
        int* sendcounts = (int*)malloc(sizeof(int) * comm_size);
        int* senddispls = (int*)malloc(sizeof(int) * comm_size);
        int* recvcounts = (int*)malloc(sizeof(int) * comm_size);
        int* recvdispls = (int*)malloc(sizeof(int) * comm_size);

		// calculate sendcounts, senddispls, recvcounts, recvdispls to send it to MPI_Alltoallv to perform merging
        Set_Sending_Size(comm_size, rank, left_child_sizes, final_sizes, sendcounts, senddispls); 
		Set_Receiving_Size(comm_size, rank, left_child_sizes, final_sizes, recvcounts, recvdispls); 
		//Sends data from all to all processes
        MPI_Alltoallv(left_child, sendcounts, senddispls, MPI_INT, new_array, recvcounts, recvdispls, MPI_INT, comm); 

        // set the size of all right_child
        for(i = 0; i < my_left_child_size; i++){
            final_sizes[i] = 0; 
        } 
        for( i = my_left_child_size; i < comm_size; i++){
            final_sizes[i] = set_chunk_size(all_right_sizes, my_right_child_size, i - my_left_child_size); 
        }

		// calculate sendcounts, senddispls, recvcounts, recvdispls to send it to MPI_Alltoallv to perform merging
        Set_Sending_Size(comm_size, rank, right_child_sizes, final_sizes, sendcounts, senddispls); 
		Set_Receiving_Size(comm_size, rank, right_child_sizes, final_sizes, recvcounts, recvdispls);
        //Sends data from all processes to all processes
		MPI_Alltoallv(right_child, sendcounts, senddispls, MPI_INT, new_array, recvcounts, recvdispls, MPI_INT, comm); 

        // we change our array according to the new arrays that is received after 
		// assigning values to left children and right children
        for( i = 0; i < new_size; i++){
            begin[i] = new_array[i]; 
        }
        
		// now we split the common space to smaller spaces
        int color = ( (rank < my_left_child_size) ? 0 : 1); 
        int key = rank; 
        MPI_Comm sub_comm; 
        MPI_Comm_split(comm, color, key, &sub_comm);  

		// free all the temporary variables
        MPI_Barrier(comm); 
        free(left_child_sizes); 
        free(right_child_sizes); 
        free(sendcounts); 
        free(senddispls);
        free(recvcounts); 
        free(recvdispls); 
        free(final_sizes); 
        free(all_subarray_sizes); 
        free(new_array); 

		// call sort recursively again, with less number of processors
        sort(begin, new_size , sub_comm); 


///////////////////////////

	}

}

void Set_Sending_Size(int comm_size, int myrank, int* initial_sizes, int* final_sizes, int* sendcounts, int* senddispls)
{

	int i = 0 ;
   
    for(i = 0; i < comm_size; i++)
	{
        sendcounts[i] = 0;  senddispls[i] = 0;
    }
	// calculate senders 
    int residual = initial_sizes[myrank],d; 
    int receiver = 0; 
    int right_elements = final_sizes[0]; 
    int total_residuals = 0; 
	
    for(i = 0; i < myrank; i++)
	{
        total_residuals += initial_sizes[i]; 
    }
    while(total_residuals > 0){
        if(total_residuals < right_elements)
		{
            d = total_residuals; 
            right_elements -= d;
        } else
		{
            d = right_elements; 
            receiver += 1; 
            right_elements = final_sizes[receiver];
        }
        total_residuals -= d; 
    }
    residual = initial_sizes[myrank]; 
    while(residual > 0)
	{
        if(residual < right_elements)
		{
            d = residual; 
            right_elements -= d; 
            sendcounts[receiver] = d; 
            senddispls[receiver] = initial_sizes[myrank] - residual; 
       
        } else
		{
            d = right_elements; 
            sendcounts[receiver] = d; 
            senddispls[receiver] = initial_sizes[myrank] - residual; 
            receiver += 1; 
            right_elements = final_sizes[receiver]; 
        } 
        residual -= d; 
    }
}


// calculate sendcounts, senddispls, recvcounts, recvdispls to send it to MPI_Alltoallv to perform merging
void Set_Receiving_Size(int comm_size, int myrank, int* initial_sizes, int* final_sizes, int* recvcounts, int* recvdispls){

	int i = 0 ; 
    for(i = 0; i < comm_size; i++)
	{
        recvcounts[i] = 0; 
        recvdispls[i] = 0; 
    }
    // calculate receiver 
    int residual = final_sizes[myrank], d; 
    int sender = 0; 
    int left_elements = initial_sizes[0];
    int total_residuals = 0; 
    for(i = 0; i < myrank; i++)
	{
        total_residuals += final_sizes[i]; 
    }
    while(total_residuals > 0)
	{
        if(total_residuals < left_elements)
		{
            d = total_residuals; 
            left_elements -= d; 
        } else
		{
            d = left_elements; 
            sender += 1; 
            left_elements = initial_sizes[sender]; 
        }
        total_residuals -= d; 
    }
    residual = final_sizes[myrank];
    while(residual > 0)
	{
        if(residual < left_elements)
		{
            d = residual; 
            left_elements -= d; 
            recvcounts[sender] = d; 
            recvdispls[sender] = final_sizes[myrank] - residual; 
        } else{
            d = left_elements; 
            recvcounts[sender] = d; 
            recvdispls[sender] = final_sizes[myrank] - residual; 
            sender += 1; 
            left_elements = initial_sizes[sender]; 
        }
        residual -= d; 
    }
}

int main( int  argc, char **argv)
{

	int rank, num, size, i;

	struct timeval Tstart, Tend; //time value

	// No MPI call before it. 
	// MPI start
	MPI_Init(&argc, &argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num );


	//all processors stop.
	MPI_Barrier (MPI_COMM_WORLD);
	
	if(rank == 0)
	{
		printf("Get Ready.\n");
		// start time
		gettimeofday(&Tstart, NULL); 
	}
// implement
	srand((rank*3+1)*time(NULL));
	//create local array and generate element

	int* local_elem;
	int local_size = NUM_OF_ELEM / num ;
	local_elem = (int*)malloc(sizeof(int) *local_size);

	for(i = 0 ; i < local_size ;i++)
	{
		local_elem[i] = rand() % MAX_VALUE;
	}

	if (num ==1){ 
		// if number of processors is only one, it means we wanted to perform a sequential quick sort
		printf ("\n a sequential quicksort is running now");
		//sequential_quick_sort(local_elem, 0 , local_size-1 ); 		
	}
	else {	
		// call quick sort to start sorting the elements
		//quick_sort(&local_elem[0], local_size, MPI_COMM_WORLD);
	}

	// all processors stop here 
    MPI_Barrier (MPI_COMM_WORLD);	

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

void sequential_quick_sort(int* x,int first,int last){

	int pivot, j, swap, i;

	if(first < last){
		// in this sequential quick sort pivoting is done 
		// based on selecting first element
		pivot = first;
		i = first;
		j = last;

		while(i<j){
			while( x[i] <= x[pivot] && i<last)
				i++;
			while( x[j] > x[pivot] )
				j--;
			if(i<j){
				swap = x[i];
				x[i] = x[j];
				x[j] = swap;
			}
		}

		swap = x[pivot];
		x[pivot] = x[j];
		x[j] = swap;
		
		// recursive call of quicksort
		sequential_quick_sort(x,first,j-1);
		sequential_quick_sort(x,j+1,last);

	}
}

