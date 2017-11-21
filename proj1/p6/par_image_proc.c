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
void derivedPixelType()
{
	Pixel sample;

	int block_len[2];
	
	
	MPI_Datatype types[2];
	MPI_Aint displ[2], off, origin;

	block_len[0] = block_len[1] = 1;

	types[0] = MPI_INT;
	types[1] = RGB_Type;

	displ[0] = 0;
	MPI_Get_address(&sample.offset, &origin);
	MPI_Get_address(&sample.rgb, &off);
	displ[1] = off - origin;

	MPI_Type_create_struct(2, blklen, displ, types, &Pixel_Type);
	MPI_Type_commit(&Pixel_type);
}



int main(int argc, char *argv[])
{
	struct timeval Tstart, Tend; //time value

	if( argc < 3)
	{
		printf("error\n");
		return 1;
	}

	gettimeofday(&Tstart, NULL); 

	//implement

	
	
	
	
	//
	gettimeofday(&Tend, NULL);
	getElapsedTime(Tstart, Tend);

	fnWritePPM(argv[2], &img);
	fnClosePPM(&img);

}
