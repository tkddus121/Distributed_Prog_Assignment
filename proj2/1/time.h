/*
 * Please do not edit this file.
 * It was generated using rpcgen.
 */

#ifndef _TIME_H_RPCGEN
#define _TIME_H_RPCGEN

#include <rpc/rpc.h>


#ifdef __cplusplus
extern "C" {
#endif


#define TIME_PROG 0x31111111
#define TIME_VERS 1

#if defined(__STDC__) || defined(__cplusplus)
#define get_time 1
extern  int * get_time_1(void *, CLIENT *);
extern  int * get_time_1_svc(void *, struct svc_req *);
#define delay 2
extern  int * delay_1(int *, CLIENT *);
extern  int * delay_1_svc(int *, struct svc_req *);
extern int time_prog_1_freeresult (SVCXPRT *, xdrproc_t, caddr_t);

#else /* K&R C */
#define get_time 1
extern  int * get_time_1();
extern  int * get_time_1_svc();
#define delay 2
extern  int * delay_1();
extern  int * delay_1_svc();
extern int time_prog_1_freeresult ();
#endif /* K&R C */

#ifdef __cplusplus
}
#endif

#endif /* !_TIME_H_RPCGEN */