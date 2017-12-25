#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
void getElapsedTime(struct timeval Tstart, struct timeval Tend);

int main(int argc, char* argv[] )
{
	int thread_count = strtol(argv[1], NULL, 10);
	char words[30000][20];

	int pal_num = 0;
	
	if (argc != 2) {
		printf("[Usage] %s <the number of threads>\n",argv[0]);
		return 1;
	}

//	int thread_cnt = atoi(argv[1]);
//

	FILE *fp;
	FILE *fd;
	if ( (fp = fopen("words.txt","r")) == NULL)
	{
		printf("file open error.\n");
		return 0;
	}

	fd = fopen("output.txt","w");
	int k = 0;

	struct timeval Tstart, Tend;

    gettimeofday(&Tstart, NULL);

	while (1)
	{
		fscanf(fp,"%s",words[k]);
		k++;
		if(feof(fp) != 0 )
			break;
	}

#pragma omp parallel for reduction(+:pal_num) schedule(dynamic) num_threads(thread_count)
	for(int i = 0 ; i < k+1 ;++i)
	{
		//palin test
		int pal_ck = 0;
		int len = strlen(words[i]);
		for(int j = 0 ; j < len/2;++j)
		{
			if(words[i][j] == words[i][len-j-1])
					;
			else
			{
				pal_ck = 1;
				break;
			}
		}
	
		if( pal_ck == 0)
		{
			fprintf(fd,"%s\n",words[i]);
			#pragma omp critical(output)
			pal_num++;
		}
		else
		{
			//palindromic test
			//back word
			for(int j = 0 ; j < len/2; ++j)
			{			
				char tmp;
				tmp = words[i][len-j-1];
				words[i][len-j-1] = words[i][j];
				words[i][j] = tmp;
			}

			//find
			//
			for (int j = 0 ; j < k ;j++)
			{
				if( j != i)
				{
					if( strcmp(words[i],words[j]) == 0 )
					{
						fprintf(fd,"%s\n",words[j]);

						#pragma omp critical(output)
						pal_num++;
					}
				}
			}
		}
	}

	fclose(fp);
	fclose(fd);

    gettimeofday(&Tend, NULL);	// 현재시간 구하기(측정할 소스부분 수행 후에 사용)
    getElapsedTime(Tstart, Tend);	// 정의된 함수로 시간 구하기

	printf("# of palindrome words : %d\n",pal_num);

	return 0;
}
void getElapsedTime(struct timeval Tstart, struct timeval Tend)
{
    Tend.tv_usec = Tend.tv_usec - Tstart.tv_usec;
    Tend.tv_sec  = Tend.tv_sec - Tstart.tv_sec;
    Tend.tv_usec += (Tend.tv_sec*1000000);

    printf("Elapsed Time: %lf sec\n", Tend.tv_usec / 1000000.0);
}



