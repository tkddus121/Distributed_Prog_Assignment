#include "ppm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
    Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
    Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
    Tend.tv_usec += (Tend.tv_sec*1000000);

    printf("Elapsed Time: %lf sec\n", Tend.tv_usec / 1000000.0);
}

void img_flip_horizentally(PPMImage *img)
{

	int Height = img->height;

	for(int i = 0; i < Height ; ++i)
	{
		RGB *origin = &img->pixels[i * img->width];// H * W

		for(int j = 0 ; j < img->width /2 ; ++j)
		{
			RGB temp = origin[j];
			origin[j] = origin[img->width - j - 1];
			origin[img->width -j -1] = temp;
		}
	}

}

void img_grayscale(PPMImage *img)
{
	int H = img->height , W = img->width;

	for(int i = 0 ; i < H ; ++i)
	{
		for(int j = 0 ; j < W ;j++)
		{
			unsigned char avg = ( img->pixels[i * W + j].R
								+ img->pixels[i * W + j].G
								+ img->pixels[i * W + j].B)/3;

			img->pixels[i * W + j].R = avg;
			img->pixels[i * W + j].G = avg;
			img->pixels[i * W + j].B = avg;

		}

	}
	
	;

}


int main(int argc, char *argv[])
{
	struct timeval Tstart, Tend; //time value

	if( argc < 3)
	{
		printf("error\n");
		return 1;
	}

	PPMImage img;
	fnReadPPM(argv[1], &img);
	gettimeofday(&Tstart, NULL); 
	
	//implement

	img_flip_horizentally(&img);
	img_grayscale(&img);


	//

	gettimeofday(&Tend, NULL);
	getElapsedTime(Tstart, Tend);

	fnWritePPM(argv[2], &img);
	fnClosePPM(&img);
	
}
