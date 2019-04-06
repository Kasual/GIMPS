/*$Id$*/
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2001-2006 Guillermo Ballester Valor, Klaus Kastens
 
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
 
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    Contact to the author:
    Guillermo Ballester Valor
    c/ Cordoba 19. 18151-Ogijares (Granada), Spain. 
    gbv@oxixares.com
*/
/* This file contains the pthread include code for glucas-yeafft */
/* this also includes some code to use with OPENMP and SUNMP     */

/* Typedefs for the use of pthreads */
#if defined(_OPENMP)
# include <omp.h>
#endif

#ifdef _PTHREADS

/*The barrier for multithreaded loops */
typedef struct
  {
    int   count;
    pthread_mutex_t barrier_mutex;
    pthread_cond_t barrier_cond;
    char  desc[32];
  }
y_barrier_t;

/* REIX */
void mutex_cleaner(void *);
void y_barrier(int, y_barrier_t *);

/*
#define	PTHREAD_TRACE	1
*/

#define	PTHREAD_JOIN	1
#if defined(PTHREAD_TRACE)
# define PTHREAD_trace(trace) trace
# define PTHREAD_trace2(trace) trace
#else
# define PTHREAD_trace(trace)
# define PTHREAD_trace2(trace)
#endif

/*
   This barrier differs from prior because and exit 
   after freeing 
*/
#define y_last_barrier( _num, _b)                 \
  y_barrier( _num, _b);				  \
  PTHREAD_trace(fprintf(stderr, "Auto canceling thread %p.\n", (void *) pthread_self());); \
  pthread_exit(PTHREAD_CANCELED);


/*The args struct to pass in glucas_loop*/

typedef struct
  {
    int   id;  /* The # of thread being passed */
    int   nt;  /* number of threads in the family */
    int   n1;  /* number of threads with work to do in pass 1 */
    int   n2;  /* number of threads with work to do in pass 2 */
    double *ws; /* The big array pointer, the work area */
    unsigned long N; /*Number of complex in big array*/
    int    ir;  /* pass2 threshold */
  }
y_gparm;


#endif /* _PTHREADS */


/************************************************************************
	  GLOBAL STATIC VARS 
************************************************************************/
#if defined(_PTHREADS)
extern y_barrier_t gbarrier1, gbarrier2, gbarrier_start, gbarrier_finished;
extern pthread_attr_t y_attr;
extern pthread_t *y_threads;
extern y_gparm *garg;
extern int Y_ERR_FLAG;
# if defined(linux)
extern y_ptr y_user;
# endif
#endif

#if defined(Y_USE_SSE2)
extern y_ptr *Y_XBJS, *Y_XBJS0;
#endif

extern y_ptr *Y_CARRIES, *Y_CARRIES0, *Y_TWX, Y_ERRS;
#if defined(SUM_CHECK)
extern y_ptr Y_SUMOUT, Y_SUMIN;
#endif
extern UL **Y_BJS;
extern int Y_NTHREADS, Y_ERR_FLAG, Y_THREADS_TERMINATE;
extern int Y_T1, Y_T2;
extern int *Y_2N0, *Y_2NN, *Y_DI, *Y_DJ, Y_BS, Y_IR, Y_NC;
extern int *Y_1NN, *Y_1NC;

/************************************************************************
	  FUNCTION PROTOTYPES  
************************************************************************/

int y_threads_init_pass1(int);

int y_threads_init_pass2(int);

int y_check_threshold_pass2(void);

void y_pass2_mp(y_ptr, int, int *, int *, int, int, int, int);


BIG_DOUBLE dit_carry_norm_dif_5_mp(y_ptr, y_ptr, y_ptr,
                                   UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_5_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_6_mp(y_ptr, y_ptr, y_ptr,
                                   UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_6_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_7_mp(y_ptr, y_ptr, y_ptr,
                                   UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_7_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_8_mp(y_ptr, y_ptr, y_ptr,
                                   UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_8_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_9_mp(y_ptr, y_ptr, y_ptr,
                                   UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_9_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_10_mp(y_ptr, y_ptr, y_ptr,
                                    UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_10_last_carries_mp(BIG_DOUBLE *, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_12_mp(y_ptr, y_ptr, y_ptr,
                                    UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_12_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_14_mp(y_ptr, y_ptr, y_ptr,
                                    UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_14_last_carries_mp(y_ptr, UL *, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_16_mp(y_ptr, y_ptr, y_ptr,
                                    UL *, UL, UL, UL, UL);

void dit_carry_norm_dif_16_last_carries_mp(y_ptr, UL *, UL, UL, int);

#if defined(Y_USE_SSE2)
BIG_DOUBLE dit_carry_norm_dif_5_mp_sse2(y_ptr, y_ptr, y_ptr, y_ptr,
                                        UL, UL, UL, UL);

void dit_carry_norm_dif_5_last_carries_mp_sse2(y_ptr, y_ptr, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_6_mp_sse2(y_ptr, y_ptr, y_ptr, y_ptr,
                                        UL, UL, UL, UL);

void dit_carry_norm_dif_6_last_carries_mp_sse2(y_ptr, y_ptr, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_7_mp_sse2(y_ptr, y_ptr, y_ptr, y_ptr,
                                        UL, UL, UL, UL);

void dit_carry_norm_dif_7_last_carries_mp_sse2(y_ptr, y_ptr, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_8_mp_sse2(y_ptr, y_ptr, y_ptr, y_ptr,
                                        UL, UL, UL, UL);

void dit_carry_norm_dif_8_last_carries_mp_sse2(y_ptr, y_ptr, UL, UL, int);

BIG_DOUBLE dit_carry_norm_dif_9_mp_sse2(y_ptr, y_ptr, y_ptr, y_ptr,
                                        UL, UL, UL, UL);

void dit_carry_norm_dif_9_last_carries_mp_sse2(y_ptr, y_ptr, UL, UL, int);

#endif /* Y_USE_SSE2 */


#if defined(_PTHREADS)

void y_create_lucas_threads(int, int, int);

void y_modify_lucas_threads(int, int, int);

void y_barrier_init( y_barrier_t *, char *);

void y_barrier_destroy(y_barrier_t *);

void y_init_attr( pthread_attr_t *);

void y_threads_init(void);

void y_threads_destroy(void);

void y_cancel_all_threads(void);

void y_join_and_freemem_all_threads(void);

void * glucas_loop( void *);

# if defined(linux)
double sec_elapsed_user_thread(void);
# endif
#endif

/*$Id$*/


