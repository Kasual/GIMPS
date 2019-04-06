/*$Id$*/
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006 Guillermo Ballester Valor, Klaus Kastens
 
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
/* This file contains the thread tasks */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#define NDEBUG1

/******************************************************************
 yeafft is the include file for FFT. Read it to obtain more 
 information about it 
******************************************************************/
#include "yeafft.h"
#include "glucas.h"

#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)
int *Y_2N0 = NULL, *Y_2NN = NULL, *Y_DI = NULL, *Y_DJ = NULL,
                            Y_BS, Y_IR, Y_NC, Y_T1, Y_T2;
int *Y_1NN = NULL, *Y_1NC = NULL, Y_ERR_FLAG, Y_THREADS_TERMINATE = 0;

#if defined(Y_USE_SSE2)
y_ptr *Y_XBJS, *Y_XBJS0;
#endif

y_ptr *Y_CARRIES = NULL, *Y_CARRIES0 = NULL, *Y_TWX = NULL, Y_ERRS = NULL;

#if defined(SUM_CHECK)
y_ptr  Y_SUMOUT = NULL, Y_SUMIN = NULL;
#endif

UL **Y_BJS = NULL;
int Y_NTHREADS = 0;
# if defined(linux)
y_ptr y_user = NULL;
# endif
#endif
# ifdef _PTHREADS
y_barrier_t gbarrier1, gbarrier2, gbarrier3, gbarrier_start, gbarrier_finished;
pthread_attr_t y_attr;
pthread_t *y_threads = NULL;
y_gparm *garg = NULL;


/************************************************************************
	  BARRIER
************************************************************************/
/*
   this function define a barrier. to free the 'b' barrier it is needed 
   '_num' threads reaching the barrier
*/

void y_barrier( int _num, y_barrier_t * _b)
{
  int ret;

  PTHREAD_trace(fprintf(stderr, "-Th: %p  In  y_barrier: %s\n", (void *) pthread_self(), (_b->desc)););

  if( (ret = pthread_mutex_lock(&( _b->barrier_mutex))) != 0 )
    {
      fprintf(stderr, "  : %p  Error %d while taking the lock\n", (void *) pthread_self(), ret);
    }

  _b->count += 1;

  PTHREAD_trace(fprintf(stderr, "   : %p  in  %s->count: %d  _num = %d\n", (void *) pthread_self(), _b->desc, _b->count, _num););

  if( _b->count >= _num )
    {
      _b->count = 0;
      if( (ret = pthread_cond_broadcast( &(_b->barrier_cond))) != 0 )
        {
          fprintf(stderr, "Error %d while waiting in cond_wait\n", ret);
        }
    }
  else
    {
      pthread_cleanup_push(mutex_cleaner, (void *)&(_b->barrier_mutex));
      do
        {
          if( (ret = pthread_cond_wait(&( _b->barrier_cond), &( _b->barrier_mutex))) != 0 )
            {
              fprintf(stderr, "Error %d while waiting in cond_wait\n", ret);
            }
        }
      while (_b->count != 0);

      pthread_cleanup_pop(0);
    }

  PTHREAD_trace(fprintf(stderr, "   : %p  out %s->count: %d  _num = %d\n", (void *) pthread_self(), _b->desc, _b->count, _num););

  if( (ret = pthread_mutex_unlock(&( _b->barrier_mutex))) != 0 )
    {
      fprintf(stderr, "  : %p  Error %d while releasing the lock\n", (void *) pthread_self(), ret);
    }

  PTHREAD_trace(fprintf(stderr, "-Th: %p  Out y_barrier: %s\n", (void *) pthread_self(), _b->desc););
}


/*This function inits the barrier */

void y_barrier_init( y_barrier_t *barr, char * desc)
{
  int ret;

  strncpy(barr->desc, desc, 31);

  PTHREAD_trace(fprintf(stderr, "#Init    y_barrier: %s\n", barr->desc););

  if( (ret = pthread_mutex_init (&( barr->barrier_mutex), NULL)) != 0 )
    {
      fprintf(stderr, "Error %d while initializing the mutex of barrier: %s\n", ret, barr->desc);
    }

  if( (ret = pthread_cond_init (&( barr->barrier_cond), NULL)) != 0 )
    {
      fprintf(stderr, "Error %d while initializing the condvar of barrier: %s\n", ret, barr->desc);
    }
  barr->count = 0;
}

/* destroy a defined barrier */
void y_barrier_destroy(y_barrier_t *barr)
{
  int ret;

  PTHREAD_trace(fprintf(stderr, "#Destroy  y_barrier: %s\n", barr->desc););

  if( (ret = pthread_mutex_destroy (&( barr->barrier_mutex))) != 0 )
    {
      fprintf(stderr, "Error %d while destroying the mutex of barrier: %s\n", ret, barr->desc);
    }

  if( (ret = pthread_cond_destroy (&( barr->barrier_cond))) != 0 )
    {
      fprintf(stderr, "Error %d while destroying the condvar of barrier: %s\n", ret, barr->desc);
    }
}


/* Inits the attr of threads */

void y_init_attr( pthread_attr_t *attr)
{
  int ret;

  PTHREAD_trace(fprintf(stderr, "#Init attributes of threads\n"););

  if( (ret = pthread_attr_init ( attr)) != 0 )
    {
      fprintf(stderr, "Error %d while initializing a pthread attribute\n", ret);
    }

  /* REIX */
  if( (ret = pthread_attr_setdetachstate (attr, PTHREAD_CREATE_JOINABLE)) != 0 )
    {
      fprintf(stderr, "Error %d while setting the detach state of a thread attribute\n", ret);
    }

  if( (ret = pthread_attr_setinheritsched (attr, PTHREAD_INHERIT_SCHED)) != 0 )
    {
      fprintf(stderr, "Error %d while setting the inherit sched of a thread attribute\n", ret);
    }
}

/* Inits all barrier needed in a normal L-L iteration cycle */
void y_threads_init(void)
{
  PTHREAD_trace(fprintf(stderr, "#Init    all barriers\n"););
  y_init_attr (&y_attr );
  y_barrier_init (&gbarrier1, "gbarrier1");
  y_barrier_init (&gbarrier2, "gbarrier2");
  y_barrier_init (&gbarrier3, "gbarrier3");
  y_barrier_init (&gbarrier_start, "gbarrier_start");
  y_barrier_init (&gbarrier_finished, "gbarrier_finished");
}

/* Destroys all barrier created */
void y_threads_destroy(void)
{
  PTHREAD_trace(fprintf(stderr, "#Destroy  all barriers\n"););
  y_barrier_destroy (&gbarrier1);
  y_barrier_destroy (&gbarrier2);
  y_barrier_destroy (&gbarrier3);
  y_barrier_destroy (&gbarrier_start);
  y_barrier_destroy (&gbarrier_finished);
  pthread_attr_destroy (&y_attr );
}

/* REIX */
void mutex_cleaner(void * arg)
{
  int ret;

  if( (ret = pthread_mutex_unlock((pthread_mutex_t *)arg)) != 0 )
    {
      fprintf(stderr, "Error %d while unlocking a mutex while in cond_wait\n", ret);
    }
}
#endif

#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)

/*
   This routine initializes some array needed in 
   pass 1 (forward and backward)
 
   input: nthr required number of threads
   return actually number of threads
 
   It allocates and fills arrays Y_1NC[] and Y_1NN[]
 
   Y_1NC[i] is the number of data to process  by thread i
   
   Y_1NN[i] is the index of first data to process by thread i 
 
*/
int y_threads_init_pass1(int nthr)
{
  int n0, i, ninc, nb, n = Y_LRIGHT[0];

  /* clean memory */
  if(Y_1NC != NULL)
    free(Y_1NC);
  if(Y_1NN != NULL)
    free(Y_1NN);

  /* allocate memory */
  Y_1NC = malloc(nthr * sizeof(int));
  Y_1NN = malloc(nthr * sizeof(int));

  /* see whether is possible to init 'nthr' threads

   After ynorm_* routines, where first forward DIF step
   is made, there are Y_PLAN[0] independent data blocks.

   If we were not able to split these Y_PLAN[0] blocks in small
   chunks, then If (Y_PLAN[0] < nthr) we would assign work 
   only for Y_PLAN[0] threads (which is bad)
   
   So, in pass 1, every thread has chunks of work with a number
   of data divisible by ninc = Y_LRIGHT[Y_PASS2].  All the array
   is splited as smoothly as possible with ninc-grained chunks of
   memory.

   else if nthr <= Y_PLAN[0] then assign Y_PLAN[0]/nthr blocks
   per thread plus some aditional block if needed.
  */
  /* There are Y_PLAN[0]*Y_PLAN[1] basic undivisible blocks of length
     Y_LRIGHT[2]
  */

  ninc = Y_LRIGHT[Y_PASS2];

  if ((Y_LRIGHT[0] / ninc) >= nthr)
    Y_T1 = nthr;
  else
    Y_T1 = Y_LRIGHT[0] / ninc;

  /* minimum number of blocks per thread*/
  nb = (Y_LRIGHT[0] / ninc) / nthr;
  nb *= ninc;

  for(n0 = 0, i = 0; n0 < Y_T1; n0++)
    {
      Y_1NC[n0]= nb;
      i += nb;
    }

  /*
     Are there still blocks to assign ?
     then add one more block per thread while possible 
  */
  /*if( ninc * nthr < n )*/
  if (i < n)
    {
      for(n0 = 0; (n0 < nthr) && (i < n); n0++)
        {
          Y_1NC[n0] += ninc;
          i += ninc;
        }
    }

  /* Finally fill the arrays */
  for( n0 = 0, i = 0; n0 < nthr; n0++)
    {
      Y_1NN[n0] = i;
      i += Y_1NC[n0];
      /*printf("%d %d\n",Y_1NN[n0],Y_1NC[n0]);*/
    }

  Y_T1 = nthr;
  return nthr;
}

/*
   This routine distribute pass 2 and dyadic mul 
   among the threads. 
 
   Because every data block (with size <= Y_MEM_THRESHOLD) 
   is independent, the forward pass 2, dyadic mul and
   backward pass 2 can be made by a thread. The problem 
   here is to distribute the blocks.
 
   The dyadic mul imposes to work with paired blocks
   Y_DI[i] is the lower block index of ith pair 
   Y_DJ[i] is the upper block index of ith pair
*/
int y_threads_init_pass2(int ir)
{
  int nd, nr, kk, ii, jj, bi, bj, ninc;
  Y_BS = Y_LRIGHT[ir]; /* the size of block */
  Y_NC = Y_NRADICES - 1 - ir; /* number of steps in pass 2 */
  Y_IR = ir;
  nd = Y_LRIGHT[0] / Y_BS; /* number of blocks */
  /*printf("NUmber of blocks for pass2 =%d\n",nd);*/
  /* alloc for indexes */
  if(Y_DI != NULL)
    free(Y_DI);
  if(Y_DJ != NULL)
    free(Y_DJ);
  Y_DI = malloc (nd * sizeof(int));
  Y_DJ = malloc (nd * sizeof(int));

  /* first pair */
  Y_DI[0] = 0;
  Y_DJ[0] = 0;

  /* now computes the following pair of indexes */
  for (nr = ir - 1, kk = 1; nr >= 0; nr--)
    {
      bi = (Y_LRIGHT[nr + 1] / Y_BS);
      bj = (Y_LRIGHT[nr] / Y_BS) - 1;
      for(ii = bi, jj = bj; ii < jj; ii++, jj--)
        {
          Y_DI[kk] = ii;
          Y_DJ[kk] = jj;
          kk++;
        }
      if(ii == jj)
        {
          Y_DI[kk] = ii;
          Y_DJ[kk] = ii;
          kk++;
        }
    }
  /* and now distribute the pairs in parallel sections */
  if(Y_2NN != NULL)
    free(Y_2NN);
  if(Y_2N0 != NULL)
    free(Y_2N0);
  Y_2NN = malloc(Y_NTHREADS * sizeof(int)); /* number of pairs in a thread */
  Y_2N0 = malloc(Y_NTHREADS * sizeof(int));/* Offset of first pair in a thread */
  ninc = kk / Y_NTHREADS; /* minimum number of pairs per thread */

  /* First assign ninc pairs per thread */
  for(jj = 0, ii = 0; jj < Y_NTHREADS; jj++)
    {
      Y_2NN[jj] = ninc;
      ii += ninc;
    }

  /* Now assing the ramaining pairs */
  if( ninc * Y_NTHREADS < kk )
    {
      for(jj = 0; (jj < Y_NTHREADS) && (ii < kk); jj++)
        {
          Y_2NN[jj] += 1;
          ii++;
        }
    }

  /* Fill the arrays of index */
  for( jj = 0, ii = 0; jj < Y_NTHREADS; jj++)
    {
      Y_2N0[jj] = ii;
      ii += Y_2NN[jj];
      /*printf("%d %d\n",Y_2N0[jj],Y_2NN[jj]);*/
    }

  /* return the number of threads with assigned work */
  if(kk < Y_NTHREADS)
    {
      Y_T2 = kk;
      return kk;
    }
  else
    {
      Y_T2 = Y_NTHREADS;
      return Y_NTHREADS;
    }
}
#endif

/*
   This routine create controls Posix threads
   Once all passes initialized with some parameters computed
   the threads are created.
 
   At the begining, they will be waiting at first barrier until a start
   signal is launched. Then every task will execute the code between 
   the defined barriers in the inner loop. 
   
*/
#if defined(_PTHREADS)

void * glucas_loop(void *arg)
{
  UL N, N2, n0;
  int i, ith, nc, nr, nt, ir, n1, n2, err_flag;
  UL *bj;
  y_ptr w, twx, carry;
  int ret;

  /*
  #if defined(Y_USE_SSE2)
    ysse2_ptr carrysse2;
  #endif
  */
  y_gparm  *p = (y_gparm *) arg;

  PTHREAD_trace(fprintf(stderr, "New Thread: %p in glucas_loop\n", (void *) pthread_self()););

  HACK_ALIGN_STACK_EVEN();
  nt = p ->nt;
  ith = p->id;
  N = p->N;
  w = p->ws;
  ir = p->ir;
  nc = Y_1NC[ith];
  nr = Y_1NN[ith];
  n1 = p->n1;
  n2 = p->n2;
  N2 = (N << 1);

  n0 = Y_LRIGHT[1] / Y_NTHREADS;
  /* make n0 divisble by Y_UPDATE */
  n0 = ((n0 >> (Y_SHIFT_UPDATE - 1)) << (Y_SHIFT_UPDATE - 1));

  twx = Y_TWX[ith];
  carry = Y_CARRIES[ith];
  /*
  #if defined(Y_USE_SSE2)
    carrysse2 = MM_CARRIES[ith];
  #endif
  */
  bj = Y_BJS[ith];
#ifndef __APPLE__

  if( (ret = pthread_setcancelstate ( PTHREAD_CANCEL_ENABLE, NULL)) != 0 )
    {
      fprintf(stderr, "Error %d while setting the Cancel State of a thread\n", ret);
    }

  if( (ret = pthread_setcanceltype ( PTHREAD_CANCEL_ASYNCHRONOUS, NULL)) != 0 )
    {
      fprintf(stderr, "Error %d while setting the Cancel Type of a thread\n", ret);
    }
#endif
  do
    {
      /* First barrier. We are all ready !         */
      /* and wait the signal to start the iteration*/
      PTHREAD_trace2(fprintf(stderr, "First Barrier: gbarrier_start in thread %p\n",
                             (void *) pthread_self()););
      y_barrier ((nt + 1), &gbarrier_start);

      pthread_testcancel();

      /* last dit, carry and norm, first dif       */
      err_flag = Y_ERR_FLAG;
      switch(Y_PLAN[0])
        {
        case(5):
#if !defined(Y_USE_SSE2)
                Y_ERRS[ith] = dit_carry_norm_dif_5_mp( w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
#else

                Y_ERRS[ith] = dit_carry_norm_dif_5_mp_sse2( w, twx, carry, Y_XBJS[ith], N2, n0,
                              ith, err_flag);
#endif

          break;
        case(6):
#if !defined(Y_USE_SSE2)
                Y_ERRS[ith] = dit_carry_norm_dif_6_mp( w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
#else

                Y_ERRS[ith] = dit_carry_norm_dif_6_mp_sse2( w, twx, carry, Y_XBJS[ith], N2, n0,
                              ith, err_flag);
#endif

          break;
        case(7):
#if !defined(Y_USE_SSE2)
                Y_ERRS[ith] = dit_carry_norm_dif_7_mp( w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
#else

                Y_ERRS[ith] = dit_carry_norm_dif_7_mp_sse2( w, twx, carry, Y_XBJS[ith], N2, n0,
                              ith, err_flag);
#endif

          break;
        case(8):
#if !defined(Y_USE_SSE2)
                Y_ERRS[ith] = dit_carry_norm_dif_8_mp( w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
#else

                Y_ERRS[ith] = dit_carry_norm_dif_8_mp_sse2( w, twx, carry, Y_XBJS[ith], N2, n0,
                              ith, err_flag);
#endif

          break;
        case(9):
#if !defined(Y_USE_SSE2)
                Y_ERRS[ith] = dit_carry_norm_dif_9_mp( w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
#else

                Y_ERRS[ith] = dit_carry_norm_dif_9_mp_sse2( w, twx, carry, Y_XBJS[ith], N2, n0,
                              ith, err_flag);
#endif

          break;
#if( Y_AVAL > 3) && defined(Y_MANY_REGISTERS)

        case(10):
                Y_ERRS[ith] = dit_carry_norm_dif_10_mp(w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
          break;
        case(12):
                Y_ERRS[ith] = dit_carry_norm_dif_12_mp(w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
          break;
        case(14):
                Y_ERRS[ith] = dit_carry_norm_dif_14_mp(w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
          break;
        case(16):
                Y_ERRS[ith] = dit_carry_norm_dif_16_mp(w, twx, carry, bj, N2, n0, ith,
                                                       err_flag);
          break;
#endif

        }

      /* wait the barrier to continue with last carry propagation */
      PTHREAD_trace2(fprintf(stderr, "Second Barrier: gbarrier1\n"););
      y_barrier( nt, &gbarrier1);

      /* last carry propagation */
      switch(Y_PLAN[0])
    {
        case(5):
#if !defined(Y_USE_SSE2)
                dit_carry_norm_dif_5_last_carries_mp( w, bj, N2, n0, ith);
#else

                dit_carry_norm_dif_5_last_carries_mp_sse2( w, Y_XBJS[ith], N2, n0, ith);
#endif

          break;
        case(6):
#if !defined(Y_USE_SSE2)
                dit_carry_norm_dif_6_last_carries_mp( w, bj, N2, n0, ith);
#else

                dit_carry_norm_dif_6_last_carries_mp_sse2( w, Y_XBJS[ith], N2, n0, ith);
#endif

          break;
        case(7):
#if !defined(Y_USE_SSE2)
                dit_carry_norm_dif_7_last_carries_mp( w, bj, N2, n0, ith);
#else

                dit_carry_norm_dif_7_last_carries_mp_sse2( w, Y_XBJS[ith], N2, n0, ith);
#endif

          break;
        case(8):
#if !defined(Y_USE_SSE2)
                dit_carry_norm_dif_8_last_carries_mp( w, bj, N2, n0, ith);
#else

                dit_carry_norm_dif_8_last_carries_mp_sse2( w, Y_XBJS[ith], N2, n0, ith);
#endif

          break;
        case(9):
#if !defined(Y_USE_SSE2)
                dit_carry_norm_dif_9_last_carries_mp( w, bj, N2, n0, ith);
#else

                dit_carry_norm_dif_9_last_carries_mp_sse2( w, Y_XBJS[ith], N2, n0, ith);
#endif

          break;
#if (Y_AVAL > 3) && defined(Y_MANY_REGISTERS)

        case(10):
                dit_carry_norm_dif_10_last_carries_mp( w, bj, N2, n0, ith);
          break;
        case(12):
                dit_carry_norm_dif_12_last_carries_mp( w, bj, N2, n0, ith);
          break;
        case(14):
                dit_carry_norm_dif_14_last_carries_mp( w, bj, N2, n0, ith);
          break;
        case(16):
                dit_carry_norm_dif_16_last_carries_mp( w, bj, N2, n0, ith);
          break;
#endif

        }

      /* wait the barrier to continue to pass1 forward DIF */
      PTHREAD_trace2(fprintf(stderr, "Third Barrier: gbarrier2\n"););
      y_barrier( nt, &gbarrier2);

      /* Set the error properly, this task is left to thread 0 */
      if(ith == 0)
    {
          Err = 0.0;
          if(err_flag)
            {
              for(i = 0; i < nt; i++)
                if(Y_ERRS[i] > Err)
                  Err = Y_ERRS[i];
            }
        }

      /* Substract two, this is left to thread Y_NPTHREADS - 1 */
      if(ith == (Y_NTHREADS - 1) && Y_SBIT)
        {
          switch(Y_PLAN[0])
            {
            case(5):
#if !defined(Y_USE_SSE2)
                    substract_two_5( w, N2);
#else

                    substract_two_5_sse2( w, N2);
#endif

              break;
            case(6):
#if !defined(Y_USE_SSE2)
                    substract_two_6( w, N2);
#else

                    substract_two_6_sse2( w, N2);
#endif

              break;
            case(7):
#if !defined(Y_USE_SSE2)
                    substract_two_7( w, N2);
#else

                    substract_two_7_sse2( w, N2);
#endif

              break;
            case(8):
#if !defined(Y_USE_SSE2)
                    substract_two_8( w, N2);
#else

                    substract_two_8_sse2( w, N2);
#endif

              break;
            case(9):
#if !defined(Y_USE_SSE2)
                    substract_two_9( w, N2);
#else

                    substract_two_9_sse2( w, N2);
#endif

              break;
#if (Y_AVAL > 3) && defined(Y_MANY_REGISTERS)

            case(10):
                    substract_two_10( w, N2);
              break;
            case(12):
                    substract_two_12( w, N2);
              break;
            case(14):
                    substract_two_14( w, N2);
              break;
            case(16):
                    substract_two_16( w, N2);
              break;
#endif

            }
#if defined(SUM_CHECK)
          /* Computes SumOut and SumIn */
          for(i = 0; i < nt; i++)
        {
              SumOut += Y_SUMOUT[i];
              SumIn += Y_SUMIN[i];
            }
#endif

        }
      /* Wait barrier to substract two */
      PTHREAD_trace2(fprintf(stderr, "Fourth Barrier: gbarrier3\n"););
      y_barrier( nt, &gbarrier3);

      if (Y_PASS2 > 1)
        {
          /* forward pass1 if this thread has work to do */
          if(ith < n1)
            y_fftf_pass1_lucas_mp (w, nc, nr, 1);


          /* wait the barrier to continue next step */
          PTHREAD_trace2(fprintf(stderr, "Fifth Barrier: gbarrier2\n"););
          y_barrier( nt, &gbarrier2);
        }

      if (Y_PASS2 > 2)
        {
          if(ith < n1)
            y_fftf_pass1_lucas_mp (w, nc, nr, 2);
          /*printf("ith=%d From second fftf %d\n",ith, 2);*/

          /* wait the barrier to continue next step */
          PTHREAD_trace2(fprintf(stderr, "Fifth-1 Barrier: gbarrier3\n"););
          y_barrier( nt, &gbarrier3);
        }

      if (Y_PASS2 > 3)
        {
          if(ith < n1)
            y_fftf_pass1_lucas_mp (w, nc, nr, 3);
          /*printf("ith=%d From second fftf %d\n",ith, 3);*/

          /* wait the barrier to continue next step */
          PTHREAD_trace2(fprintf(stderr, "Fifth-2 Barrier: gbarrier2\n"););
          y_barrier( nt, &gbarrier2);
        }

      /* Pass2 */
      if(ith < n2)
        y_pass2_mp( w, Y_BS, Y_DI, Y_DJ, Y_IR, Y_NC, Y_2N0[ith],
                    Y_2NN[ith]);

      /* wait the barrier to continue to pass1 backward DIT*/
      PTHREAD_trace2(fprintf(stderr, "Sixth Barrier: gbarrier1\n"););
      y_barrier( nt, &gbarrier1);

      if (Y_PASS2 > 3)
        {
          if(ith < n1)
            y_fftb_pass1_lucas_mp (w, nc, nr, 3);

          /* wait the barrier to continue next step */
          PTHREAD_trace2(fprintf(stderr, "seventh-2 Barrier: gbarrier2\n"););
          y_barrier( nt, &gbarrier2);
        }

      if (Y_PASS2 > 2)
        {
          if(ith < n1)
            y_fftb_pass1_lucas_mp (w, nc, nr, 2);

          /* wait the barrier to continue next step */
          PTHREAD_trace2(fprintf(stderr, "seventh-1 Barrier: gbarrier3\n"););
          y_barrier( nt, &gbarrier3);
        }

      if (Y_PASS2 > 1)
        {
          /* Backward pass1 if this thread has work to do */
          if(ith < n1)
            y_fftb_pass1_lucas_mp (w, nc, nr, 1);
        }

#if defined(linux)
      /*trick to compute the user time of all threads */
      y_user[ith] = sec_elapsed_user_thread();
#endif
      /* Barrier which tells to parent thread the iter has finished */
      PTHREAD_trace2(fprintf(stderr, "eighth Barrier: gbarrier_finished\n"););
      if(Y_THREADS_TERMINATE)
        {
          y_last_barrier( (nt + 1), &gbarrier_finished);
        }
      else
        {
          y_barrier( (nt + 1), &gbarrier_finished);
        }
    }
  while(1);
  return NULL;
}

/*
   This routine create the threads. 
   nr is the number of threads to open (all busy in pass0)
   n1 is the number of threads busy in pass1 
   n2 is the number of threads busy in pass2
   
   The best scenary is where nr==n1==n2
*/

void y_create_lucas_threads(int nr, int n1, int n2)
{
  int i,ir,ret;
  ir = Y_PASS2;
  for (i = 0; i < nr; i++)
    {
      garg[i].id = i;
      garg[i].nt = nr;
      garg[i].n1 = n1;
      garg[i].n2 = n2;
      garg[i].ws = Y_WORK;
      garg[i].N = Y_LRIGHT[0];
      garg[i].ir = ir;

      PTHREAD_trace(fprintf(stderr, "#Create thread\n"););

      if( (ret = pthread_create (&y_threads[i], &y_attr, glucas_loop, (void *)(garg + i))) != 0 )
        {
          fprintf(stderr, "Glucas: unable to create a new Posix Thread. Exiting.\n");
          fprintf(stderr, "Error: %d.\n", ret);
          exit(EXIT_FAILURE);
        }
    }
}

void y_modify_lucas_threads(int nr, int n1, int n2)
{
  int i, ir;
  ir = Y_PASS2;
  for (i = 0; i < nr; i++)
    {
      garg[i].id = i;
      garg[i].nt = nr;
      garg[i].n1 = n1;
      garg[i].n2 = n2;
      garg[i].ws = Y_WORK;
      garg[i].N = Y_LRIGHT[0];
      garg[i].ir = ir;
    }
}

/*
   This routine cancel and destroy the threads. It is used when finish a work,
   restarted or interrumpted
*/
void y_cancel_all_threads(void)
{
  int i, ret;
  if (y_threads == NULL)
    return; /* no work to do */
  PTHREAD_trace(fprintf(stderr, "Cancel all threads from thread %p !\n", (void *) pthread_self()););
  for(i = 0;i < Y_NTHREADS; i++)
    {
      PTHREAD_trace(fprintf(stderr, "#Cancel thread %d %p\n", i, (void *) y_threads[i]););
      if( (ret = pthread_cancel (y_threads[i])) != 0 )
        {
          fprintf(stderr, "Error \"%s\" while cancelling a thread\n", (char *) strerror(ret));
        }
    }

# if defined(PTHREAD_JOIN)
  PTHREAD_trace(fprintf(stderr, "#Join threads from thread %p !\n", (void *) pthread_self()););
  for(i = 0;i < Y_NTHREADS; i++)
    {
      PTHREAD_trace(fprintf(stderr, "#Join thread %d %p\n", i, (void *) y_threads[i]););
      if( (ret = pthread_join (y_threads[i], NULL)) != 0 )
        {
          fprintf(stderr, "Error \"%s\" while joining a thread\n", (char *) strerror(ret));
        }
      PTHREAD_trace(fprintf(stderr, "#Joined thread %d %p\n", i, (void *) y_threads[i]););

    }
  PTHREAD_trace(fprintf(stderr, "#After Join from thread %p !\n", (void *) pthread_self()););
# endif

  y_threads_destroy();
  free((char *) y_threads);
  y_threads = NULL;
}

void y_join_and_freemem_all_threads(void)
{
  int i, ret;
  if (y_threads == NULL)
    return; /* no work to do */
# if defined(PTHREAD_JOIN)

  PTHREAD_trace(fprintf(stderr, "#Join threads from thread %p !\n", (void *) pthread_self()););

  for(i = 0;i < Y_NTHREADS; i++)
    {
      PTHREAD_trace(fprintf(stderr, "#Join thread %d %p\n", i, (void *) y_threads[i]););
      if( (ret = pthread_join (y_threads[i], NULL)) != 0 )
        {
          fprintf(stderr, "Error \"%s\" while joining a thread\n", (char *) strerror(ret));
        }
    }
  PTHREAD_trace(fprintf(stderr, "#After Join from thread %p !\n", (void *) pthread_self()););
# endif

  y_threads_destroy();
  free((char *) y_threads);
  y_threads = NULL;
}

/*
   Linux does not computes properly the user time of the whole thread 
   family, so we need to do some trick 
*/
#  if defined(linux)
#   if (defined(HAVE_SYS_RESOURCE_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
#   include <sys/resource.h>
#  endif
double sec_elapsed_user_thread(void)
{
  double user_time=0.0;
#  if (defined(HAVE_GETRUSAGE) || !defined(HAVE_CONFIG_H)) && !defined(vms) && !defined(pccompiler) && !defined(macintosh)

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  user_time = (BIG_DOUBLE)(ru.ru_utime.tv_sec) +
              1e-6*(BIG_DOUBLE)(ru.ru_utime.tv_usec);
#  endif

  return user_time;
}

# endif
#endif
