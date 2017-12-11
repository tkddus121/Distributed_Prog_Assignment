#ifndef _MY_RWL_H
#define _MY_RWL_H

#include <pthread.h>

typedef struct  pth_rwlock_t
{
	pthread_cond_t  conv_read;
	pthread_cond_t  conv_write;
	pthread_mutex_t mutex;

	int read_cnt;
	int write_cnt;
	int is_write;
}pth_rwlock_t;

int pth_rwlock_init( pth_rwlock_t* rwlock_t );
int pth_rwlock_destory(pth_rwlock_t* rwlock_t );
int pth_rwlock_rdlock(pth_rwlock_t* rwlock_t );
int pth_rwlock_wrlock(pth_rwlock_t* rwlock_t );
int pth_rwlock_unlock(pth_rwlock_t* rwlock_t  );

#endif
