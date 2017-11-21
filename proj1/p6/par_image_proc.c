#include "ppm.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef struct {

	int offset;
	RGB rgb;
} PIXEL;

int rank, size, root = 0;

MPI_Datatype RGB_type, RGB_RowGroup_type, Pixel_type;


void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
    Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
    Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
    Tend.tv_usec += (Tend.tv_sec*1000000);

    printf("Elapsed Time: %lf sec\n", Tend.tv_usec / 1000000.0);
}

void derived_RGB_Type()
{
	RGB sample;
	int block_len[3] = {1,};
	MPI_Datatype types[3];
	MPI_Aint ad_rgb[3], off, origin;

	types[0] = types[1] = types[2] = MPI_UNSIGNED_CHAR;

	ad_rgb[0] = 0;
	MPI_Get_address(&sample.R, &origin);	
	MPI_Get_address(&sample.G, &origin);
	ad_rgb[1] = off - origin;
	MPI_Get_address(&sample.B, &origin);
	ad_rgb[2] = off - origin;

	// for structural type.
	MPI_Type_create_struct( 3, block_len, ad_rgb, types, &RGB_type);
	// 
	MPI_Type_commit(&RGB_type);


}
void derived_Pixel_Type()
{
	PIXEL sample;

	int block_len[2];
	
	
	MPI_Datatype types[2];
	MPI_Aint displ[2], off, origin;

	block_len[0] = block_len[1] = 1;

	types[0] = MPI_INT;
	types[1] = RGB_type;

	displ[0] = 0;
	MPI_Get_address(&sample.offset, &origin);
	MPI_Get_address(&sample.rgb, &off);
	displ[1] = off - origin;

	MPI_Type_create_struct(2, block_len, displ, types, &Pixel_type);
	MPI_Type_commit(&Pixel_type);
}
void scatter_img_block(RGB *sendbuf, int send_block_cnt,
					   int **p_sendcounts, int **p_displs,
					   RGB **p_recvbuf, int *p_recv_block_cnt, int block_size)
{
	int i;
	int idx = 0;
	if (rank == root) 
	{
		*p_sendcounts = malloc(sizeof(int) * size);
		*p_displs = malloc(sizeof(int) * size);


		for(i = 0; i < size; ++i) 
		{
			(*p_displs)[i] = idx;

			if (i < send_block_cnt % size)
				(*p_sendcounts)[i] = send_block_cnt / size + 1;
			else
				(*p_sendcounts)[i] = send_block_cnt / size;

			idx += (*p_sendcounts)[i];
		}
	}

	if (rank < send_block_cnt % size)
		*p_recv_block_cnt = send_block_cnt / size + 1;
	else
		*p_recv_block_cnt = send_block_cnt / size;

	*p_recvbuf = malloc(sizeof(RGB) * block_size * (*p_recv_block_cnt));
	MPI_Scatterv(sendbuf, 
			*p_sendcounts, 
			*p_displs, 
			RGB_RowGroup_type,	 
			*p_recvbuf, 
			*p_recv_block_cnt, 
			RGB_RowGroup_type, 
			root, MPI_COMM_WORLD);
}



int main(int argc, char *argv[])
{
	int width,height;
	struct timeval Tstart, Tend; //time value

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	PPMImage img;

	if( argc < 3)
	{
		printf("error\n");
		return 1;
	}

	gettimeofday(&Tstart, NULL); 

	//implement
	if (rank == root) {
		fnReadPPM(argv[1], &img);
		width = img.width;
		height = img.height;
	}

	MPI_Bcast(&width, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&height, 1, MPI_INT, root, MPI_COMM_WORLD);

	derived_RGB_Type();
	derived_Pixel_Type();
	MPI_Type_vector(1, width, 1, RGB_type, &RGB_RowGroup_type);
	MPI_Type_commit(&RGB_RowGroup_type);

	RGB *data = NULL;
	int data_block_cnt;
	int *counts = NULL, *displs = NULL;
	scatter_img_block(img.pixels, height,
					  &counts, &displs,
					  &data, &data_block_cnt, width);

	printf("[%d] received data block cnt = %d\n", rank, data_block_cnt);


	
	
	
	
	//
	gettimeofday(&Tend, NULL);
	getElapsedTime(Tstart, Tend);

	fnWritePPM(argv[2], &img);
	fnClosePPM(&img);

}
