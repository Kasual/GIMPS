/*$Id$*/
/*
   (c) 2000-2006 Guillermo Ballester Valor

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

   This routine includes all routines initializing data for 
   Lucas Lehmer test and DWT
*/
/* Include Files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef _OPENMP
# include <omp.h>
#endif

/* include needed to use yeafft package */
#include "yeafft.h"
#include "ydebug.h"

/* glucas include */
#include "glucas.h"

/*
   This routine computes the current Y_SBIT assuming the initial 
   shift random SHIFT bit is bit0 (at iteration #0)
 
   Y_SBIT = ((2 ^ iter) * bit0) mod q
   
   To be careful with the overflows with 32 bits integer we will use 
   doubles. This was the reason of infamous v.19 prime95 bug
*/
UL compute_sbit(UL bit0, UL iter, UL q)
{
  UL x, aux = iter, a = 1U << 31;

  do
    {
      a >>= 1;
    }
  while ((aux & a) == 0);
  /* 2^iter mod q */
  x = 2U;
  a >>= 1;
  do
    {
#if ULONG_MAX > 0xFFFFFFFF
      /* 64 bits integer arithmetic */
      x = (x * x) % q;
#else
      /* 32 bits integer arithmetic */
      modmul_32(x, x, x, q);
#endif

      if(a & aux)
        x <<= 1;
      if( x >= q)
        x -= q;
      a >>= 1;
    }
  while (a);

  /* bit0 * 2^j mod q */
#if ULONG_MAX > 0xFFFFFFFF

  /* 64 bits integer arithmetic */
  a = (x * bit0) % q;
#else
  /* 32 bits integer arithmetic */
  modmul_32(a, x, bit0, q);
#endif

  return a;
}

/*
   This routine computes the bj[index] based on b, index, and n 
   It is used to avoid overflows. It is not fast but accurace. 
 
   It computes (_i * _b) mode _n if i>0
*/
UL bj_comp( UL _i, UL _b, UL _n)
{
  UL j, aux = _b, acu = 0;
  for (j = 1U; j <= _i; j <<= 1)
    {
      if (j & _i)
        {
          acu += aux;
          if (acu >= _n)
            acu -= _n;
        }
      aux <<= 1;
      if (aux >= _n)
        aux -= _n;
    }
  return acu;
}


/* Sets the initial L[0] according to Y_SBIT or generate one*/
void set_L0(BIG_DOUBLE *x, UL q, UL n)
{
  UL ni, ib;
  Y_SBIT = Y_SBIT0;
  if(Y_SBIT0 == 0)
    {
      x[0] = 4.0;
      return;
    }
  Y_XBITS = (BIG_DOUBLE) q / (BIG_DOUBLE) n;
  ib = Y_SBIT + 2;
  if(ib >= Y_EXPONENT)
    ib -= Y_EXPONENT;
  ni = floor ((BIG_DOUBLE) ib / Y_XBITS);
  ib -= ceil (Y_XBITS * ni);
  x[ addr(ni) ] = (BIG_DOUBLE)(1U << ib);
}

/* Init the vars needed for Lucas-Lehmer test */
void init_lucas(UL q, UL n)
{
  UL j, qn, a, i;

#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)

  UL inc;
  int nr;
#endif

  if (pttp != NULL)
    Y_FREE ((char *) pttp);
  if (pttmp != NULL)
    Y_FREE ((char *) pttmp);
  pttp = ALLOC_DOUBLES(addr(n>>SHIFT_UPDATE));
  if(pttp == NULL)
    {
      fprintf(stderr, "init_lucas: No memory to allocate DWT data.\n");
      exit(EXIT_FAILURE);
    }
  two_to_phi = ALIGN_DOUBLES(pttp);
  pttmp = ALLOC_DOUBLES(addr(n>>SHIFT_UPDATE));
  if(pttmp == NULL)
    {
      fprintf(stderr,"init_lucas: No memory to allocate DWT data.\n");
      exit(EXIT_FAILURE);
    }
  two_to_minusphi = ALIGN_DOUBLES(pttmp);

  Y_XBITS = (BIG_DOUBLE)q / (BIG_DOUBLE)n;
  low = floor((exp(floor(Y_XBITS) * M_LN2)) + 0.5);
  high = low + low;
  lowinv = 1.0 / low;
  highinv = 1.0 / high;
  b = q % n;
  c = n - b;
  /*printf("n=%ld q=%ld b=%ld\n",n,q,b);*/

  two_to_phi[0] = low;
  two_to_minusphi[0] = 2.0 / (BIG_DOUBLE)(n);
  qn = (b << SHIFT_UPDATE) % n;

  for(i = 1, j = UPDATE; j < n; j += UPDATE, i++)
    {
      a = n - qn;
      two_to_phi[ addr(i) ] = exp(a * M_LN2 / n)*low;
      two_to_minusphi[ addr(i) ] = high / (two_to_phi[ addr(i) ] * n);
      qn += (b << SHIFT_UPDATE);
      qn %= n;
    }

  Hbig = exp(c * M_LN2 / n);
  Gbig = 1 / Hbig;
  Hsmall= 0.5 * Hbig;
  Gsmall= 2.0 * Gbig;

#if defined(Y_USE_SSE2)
  /* SSE2 constants to use in carry and norm phase */
  Y_MM_SET_PD (MM_bc[0], b, b );
  Y_MM_SET_PD (MM_bc[1], (-(long)c), b );
  Y_MM_SET_PD (MM_bc[2], b, (-(long)c) );
  Y_MM_SET_PD (MM_bc[3], (-(long)c), (-(long)c) );
  Y_MM_SET_PD (MM_c, ((double)c) - 0.4, ((double)c) - 0.4);
  Y_MM_SET_PD (MM_auxt[0], Gsmall, Gsmall);
  Y_MM_SET_PD (MM_auxt[1], Gbig, Gsmall);
  Y_MM_SET_PD (MM_auxt[2], Gsmall, Gbig);
  Y_MM_SET_PD (MM_auxt[3], Gbig, Gbig);
  Y_MM_SET_PD (MM_inv[0], lowinv, lowinv);
  Y_MM_SET_PD (MM_inv[1], highinv, lowinv);
  Y_MM_SET_PD (MM_inv[2], lowinv, highinv);
  Y_MM_SET_PD (MM_inv[3], highinv, highinv);
  Y_MM_SET_PD (MM_Hsmall, Hsmall, Hsmall);
#endif /* Y_USE_SSE2 */

#if defined(SUM_CHECK)
  /*
     this code, used to compute the limit of expected errors
     on comparing sum of inputs with outputs, is from 
     G. Wotlman PRIME95.

     We will need to assure than 

     abs (SumIn * SumIn - SumOut * 2 / n ) < limit
  */
  {
    double bits_per_double, total_bits, loglen;

    bits_per_double = (double) q / (double) n - 1.0;
    if (bits_per_double < 15.0)
      bits_per_double = 15.0;
    loglen = log ((double) n) / log (2.0);
    loglen *= 0.69;
    total_bits = bits_per_double * 2.0 + loglen * 2.0;
    ErrLimit = pow (2.0, total_bits - 49.08);
  }
  /*printf("ErrLimit = %lf\n",ErrLimit);
    ErrLimit = 100;*/
  /* this is the Sumout factor , 2 / n */

  SumNorm = (BIG_DOUBLE)2.0 / n;

#endif /* SUM_CHECK */

#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)
  /*
     This is added for MP threading of last_DIT-carry-norm-first_DIFF
     step
  */
  if(!Y_NTHREADS)
    {
# if defined(_OPENMP)
#  ifdef Y_NUM_THREADS
#   pragma omp parallel default(shared)

      {
        omp_set_dynamic (Y_NUM_THREADS );
        omp_set_num_threads (Y_NUM_THREADS);
        Y_NTHREADS = Y_NUM_THREADS;
      }
#  else
#   pragma omp parallel default(shared)

      {
        Y_NTHREADS = omp_get_num_procs();
        omp_set_dynamic (Y_NTHREADS );
        omp_set_num_threads (Y_NTHREADS);
      }
#  endif
# elif defined(_SUNMP)
      Y_NTHREADS = Y_NUM_THREADS;
# else

      if(Y_NTHREADS == 0)
        Y_NTHREADS = _PTHREADS;
# endif

    }

  if( Y_CARRIES == NULL)
    {
      Y_CARRIES = (y_ptr *) Y_ALLOC (Y_NTHREADS * sizeof(y_ptr));
      Y_CARRIES0 = (y_ptr *) Y_ALLOC (Y_NTHREADS * sizeof(y_ptr));
      for(i = 0; (int)i < Y_NTHREADS; i++)
        {
          Y_CARRIES[i] = NULL;
          Y_CARRIES0[i] = NULL;
        }
    }
  else
    {
      for(i = 0; (int)i < Y_NTHREADS; i++)
        {
          Y_FREE ((char *) Y_CARRIES0[i]);
          Y_CARRIES0[i] = NULL;
          Y_CARRIES[i] = NULL;
        }
    }

  if( Y_BJS == NULL)
    {
      Y_BJS = Y_ALLOC (Y_NTHREADS * sizeof(UL *));
      for(i = 0; (int)i < Y_NTHREADS; i++)
        Y_BJS[i] = NULL;
    }
  else
    {
      for(i = 0; (int)i < Y_NTHREADS; i++)
        {
          Y_FREE ((char *) Y_BJS[i]);
          Y_BJS[i] = NULL;
        }
    }

#if defined(Y_USE_SSE2)
  if( Y_XBJS == NULL)
    {
      Y_XBJS0 = Y_ALLOC (Y_NTHREADS * sizeof(y_ptr *));
      Y_XBJS = Y_ALLOC (Y_NTHREADS * sizeof(y_ptr *));
      for(i = 0; (int)i < Y_NTHREADS; i++)
        {
          Y_XBJS[i] = NULL;
          Y_XBJS0[i] = NULL;
        }
    }
  else
    {
      for(i = 0; (int)i < Y_NTHREADS; i++)
        {
          Y_FREE ((char *) Y_XBJS0[i]);
          Y_XBJS[i] = NULL;
          Y_XBJS0[i] = NULL;
        }
    }
#endif

  if( Y_ERRS == NULL)
    Y_ERRS = Y_ALLOC (Y_NTHREADS * sizeof(BIG_DOUBLE));

#if defined(SUM_CHECK)

  if( Y_SUMOUT == NULL)
    Y_SUMOUT = Y_ALLOC (Y_NTHREADS * sizeof(BIG_DOUBLE));
  if( Y_SUMIN == NULL)
    Y_SUMIN = Y_ALLOC (Y_NTHREADS * sizeof(BIG_DOUBLE));
#endif

  if( Y_TWX == NULL)
    Y_TWX = Y_ALLOC (Y_NTHREADS * sizeof(y_ptr));

  for (i = 0; (int)i < Y_NTHREADS; i++)
    {
      Y_CARRIES0[i] = (y_ptr) ALLOC_DOUBLES(Y_PLAN[0] + 2);
      Y_CARRIES[i] = ALIGN_DOUBLES(Y_CARRIES0[i]);

# if defined(Y_USE_SSE2)

      Y_XBJS0[i] = ALLOC_DOUBLES( Y_PLAN[0] + 2);
      Y_XBJS[i] = ALIGN_DOUBLES( Y_XBJS0[i]);
#endif

      Y_BJS[i] = Y_ALLOC (Y_PLAN[0] * sizeof(UL));
    }

  /* Init the Twiddle pointers */
  a = (Y_LRIGHT[1] << 1) / Y_NTHREADS;

  /* make a divisble by Y_UPDATE */
  a = ((a >> Y_SHIFT_UPDATE) << Y_SHIFT_UPDATE );
  Y_TWX[0] = Y_TWDB[Y_NRADICES -2];
  if(Y_NTHREADS > 1)
    {
      inc = (Y_POWERS[Y_PLAN[0]-1]);
      for (i = 1; (int)i < Y_NTHREADS; i++)
        {
          Y_TWX[i] = Y_TWDB[Y_NRADICES - 2] + (i * a - 2) * inc;
        }
    }
  /* Init the Bjs */
  /*
    a is the number of complex Y_PLAN[0]-FFTS to make by every thread:
    a=pad/Y_NTHREADS = Y_LRIGTH[1]/Y_NTHREADS 

    Then
    Y_BJS[i][j] = (2 * (a*i + j*Y_LRIGHT[1]) * b) % N;

    And if we suposse that 'a' is a power_of_two that divides N and 
    N/Y_PLAN, then:
    y_BJS[i][j] = ((2 * a *(i + j*Y_NTHREADS)) * b) % N =
                = (((i + j*Y_NTHREADS)*b) % (N/(2*a))) * (2*a)

  */
  for (i = 0; (int)i < Y_NTHREADS; i++)
    {
      for(j = 0; j < Y_PLAN[0]; j++)
        {
          /*Y_BJS[i][j] = bj_comp( a * i + 2 * j * Y_LRIGHT[1], b, n);*/
          modmul_32(Y_BJS[i][j], a * i + 2 * j * Y_LRIGHT[1], b, n);
          /*Y_BJS[i][j] = ((( i + j * Y_NTHREADS) * b) % (n / a)) * a;*/
#if defined(Y_USE_SSE2)

          Y_XBJS[i][j] = (BIG_DOUBLE) Y_BJS[i][j];
#endif

        }
    }
  Y_BJS[0][0] = n;

#if defined(Y_USE_SSE2)

  Y_XBJS[0][0] = (BIG_DOUBLE) n;
#endif

  /* also init pass1 and pass2 */
  nr = Y_PASS2;
  y_threads_init_pass1 (Y_NTHREADS);
  y_threads_init_pass2 (nr);
  Y_THREADS_TERMINATE = 0;

  /*Check whether the exponent is too small for threaded pass0 */
  if(Y_UPDATE > (2 * Y_LRIGHT[1] / Y_NTHREADS))
    {
      fprintf (stderr, "Glucas: Exponent too small for this threaded version.\n");
      exit (EXIT_FAILURE);
    }

  /* and create/modify the threads */
# if defined(_PTHREADS)
  if( garg != NULL)
    free((char *) garg);
  garg = (y_gparm *) Y_ALLOC (Y_NTHREADS * sizeof(y_gparm));
#  if defined(linux)
  /* user time for threads */
  if(y_user == NULL)
    y_user = Y_ALLOC (Y_NTHREADS * sizeof(double));
#  endif

  y_cancel_all_threads();

  y_threads_init();
  y_threads = (pthread_t *) Y_ALLOC (Y_NTHREADS * sizeof(pthread_t));
  if(y_threads == NULL)
    {
      fprintf(stderr,"Glucas: Could not allocate space for threads. Exiting\n");
      exit(EXIT_FAILURE);
    }
  y_create_lucas_threads (Y_NTHREADS, Y_T1, Y_T2);

# endif
#endif
}
