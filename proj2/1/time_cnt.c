#include "time.h"
#include "time_cnt.h"
#include <sys/time.h>
#include <pthread.h>
#include <string.h>
#include <signal.h>
#include <stdio.h>

static int tod = 0;

static void timer_handler (int signum)
{
	__sync_fetch_and_add(&tod,1);
}

/* timer source at googling.*/
static void* timer_thread_main()
{
    struct sigaction sig;
    struct itimerval t;

    /* Install timer_handler as the signal handler for SIGVTALRM. */
    memset (&sig, 0, sizeof (sig));
    sig.sa_handler = &timer_handler;
    sigaction (SIGVTALRM, &sig, NULL);

    t.it_value.tv_sec = 1;
    t.it_value.tv_usec = 0;

    t.it_interval.tv_sec = 1;
    t.it_interval.tv_usec = 0;

    /* Start a virtual timer. It counts down whenever this process is executing. */
    setitimer (ITIMER_VIRTUAL, &t, NULL);

	/* Do busy work. */
    while(1);
    pthread_exit(0);
}

int timer_get_tod()
{
	int ret = __sync_fetch_and_add(&tod, 0);
    return ret;
}

int timer_start()
{
	//new pthread for timer.
	pthread_t tid; 
	pthread_attr_t attr;
	
	pthread_attr_init(&attr);
	pthread_create( &tid, &attr, timer_thread_main, NULL);
	return 0;
}

int timer_delay(int interval){
	int init_tod = timer_get_tod();
	int now_tod = init_tod;

	printf("Start delay by %d sec. \n", interval);

	for(;now_tod - init_tod < interval;)
		now_tod = timer_get_tod();

	printf("End   delay by %d sec. \n",interval);
	return now_tod;

}
