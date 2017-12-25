#include "my_pth_rwl.h"

int pth_rwlock_init(pth_rwlock_t* rwlock_t )
{
	if( pthread_mutex_init(&rwlock_t->mutex, NULL) != 0)
	{
		return -1;
	}
	if( pthread_cond_init(&rwlock_t->conv_read, NULL) != 0)
	{
		pthread_mutex_destroy(&rwlock_t->mutex);
		return -1;
	}
	if(  pthread_cond_init(&rwlock_t->conv_write, NULL) != 0)
	{
		pthread_cond_destroy(&rwlock_t->conv_read);
		return -1;
	}

	rwlock_t->read_cnt = 0;
	rwlock_t->write_cnt = 0;
	rwlock_t->is_write = 0;
	return 0;
}

int pth_rwlock_rdlock(pth_rwlock_t* rwlock_t )
{
	pthread_mutex_lock(&rwlock_t->mutex);
	(rwlock_t->read_cnt)++;
	
	if( rwlock_t->is_write)
		pthread_cond_wait( &rwlock_t->conv_read, &rwlock_t->mutex);

	pthread_mutex_unlock(&rwlock_t->mutex);
	return 0;

}
int pth_rwlock_wrlock(pth_rwlock_t* rwlock_t )
{
	pthread_mutex_lock(&rwlock_t->mutex);
	(rwlock_t->write_cnt)++;
	
	if( rwlock_t->is_write || rwlock_t->read_cnt > 0 )
		pthread_cond_wait( &rwlock_t->conv_write, &rwlock_t->mutex);
	rwlock_t->is_write = 1;
	pthread_mutex_unlock(&rwlock_t->mutex);
	return 0;

}
int pth_rwlock_unlock(pth_rwlock_t* rwlock_t  )
{
	pthread_mutex_lock( &rwlock_t->mutex );

	if( !(rwlock_t->is_write) )
	{
		// writing..
		(rwlock_t->read_cnt)--;

		if( rwlock_t->read_cnt <= 0)
		{
			if( rwlock_t->write_cnt > 0)
			{
				rwlock_t->is_write = 1;
				pthread_cond_signal(&rwlock_t->conv_write);
			}	
		}
		else
		{
			pthread_cond_broadcast(&rwlock_t->conv_read);
		}
	}
	else
	{
		// not writing.
		(rwlock_t->write_cnt)--;
		
		if( rwlock_t->read_cnt > 0 )
		{
			
			rwlock_t->is_write = 0;
			pthread_cond_broadcast(&rwlock_t->conv_read);
		}
		else if( rwlock_t->write_cnt <= 0)
		{
			rwlock_t->is_write = 0;
		}
		else
		{
			pthread_cond_signal(&rwlock_t->conv_write);
		}

	}

	pthread_mutex_unlock(&rwlock_t->mutex);

	return 0;


}

int pth_rwlock_destory(pth_rwlock_t* rwlock_t )
{
	int fail = 0;
	fail = fail || pthread_mutex_destroy(&rwlock_t->mutex);
	fail = fail ||  pthread_cond_destroy(&rwlock_t->conv_read);
	fail = fail || pthread_cond_destroy(&rwlock_t->conv_write);

//	if( pthread_mutex_destroy(&rwlock_t->mutex) || 
//			pthread_cond_destroy(&rwlock_t->conv_read) || 
//			pthread_cond_destroy(&rwlock_t->conv_write) || 0)
	if(fail)
		return -1;
	else
		return 0;

}
