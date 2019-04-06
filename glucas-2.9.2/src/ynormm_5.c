/*$Id$*/
/*
   (c) 2003-2006 Guillermo Ballester Valor 

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
   
   
   It includes ideas from many Authors.
 
   This is the sse2 version of ynorm_5.c
 
   When it gets the control to perform the last backward DIT pass
   the memory are in fast scramble mode, i.e.
 
   d[0] = real part of first complex data;
   d[1] = real part of second complex data;
   d[2] = imag part of first complex data;
   d[3] = imag part of second complex data;
 
   this four double will be loaded respectively as Y__M128D t0r, t0i 
 
   t0r.d[0] <- data0[0]
   t0r.d[1] <- data0[2]
   t0i.d[0] <- data0[1]
   t0i.d[1] <- data0[3]
 
   idem for t1 pair  
 
   e = d + pad1
   t1r.d[0] <- data1[0]
   t1r.d[1] <- data1[2]
   t1i.d[0] <- data1[1]
   t1i.d[1] <- data1[3]
 
   A) So, actually we will perform two normal (non_SSE2) loops in just one.
 
   b) After last DIT pass, the data needs to be scrambled a bit prior to carry
   and normalization phase. This carry_and_norm will be as vectorized as
   possible. We call #i chain one of possible 0 to (n-1) radixn sequence of 
   carry and normalization data. It will join two chain of data when posible.
   For odd radices (5 and 9) last chain will be paired with a bogus 0.0 
   sequence. 
 
   Example of how it joins #0 and #1 chain: 
 
     Y_MM_MOV_INTER_PAIR( t0r, t1r) :
       t0r.d[0] <- t0r.d[0] ( data0[0] )
       t0r.d[1] <- t1r.d[0] ( data1[0] )
       t1r.d[0] <- t0r.d[1] ( data0[2] )
       t1r.d[1] <- t1r.d[1] ( data1[2] ) 
 
     Y_MM_MOV_INTER_PAIR( t0i, t1i) :
       t0i.d[0] <- t0i.d[0] ( data0[1] )
       t0i.d[1] <- t1i.d[0] ( data1[1] )
       t1i.d[0] <- t0i.d[1] ( data0[3] )
       t1i.d[1] <- t1i.d[1] ( data1[3] )
 
  B) carry_and_norm_vector #0 and #1 data in the following order: 
       t0r -> first double on #0 and #1 
       t0i -> second double on #0 and #1
       t1r -> third double on #0 and #1
       t1i -> fourth double on #0 and #1
       
first real part and then imag
   
  C) disjoin #0 and #1 data
     Y_MM_MOV_INTER_PAIR( t0r, t1r) :
       t0r.d[0] <- t0r.d[0]
       t0r.d[1] <- t1r.d[0]
       t1r.d[0] <- t0r.d[1]
       t1r.d[1] <- t1r.d[1]
 
     Y_MM_MOV_INTER_PAIR( t0i, t1i) :
       t0i.d[0] <- t0i.d[0]
       t0i.d[1] <- t1i.d[0]
       t1i.d[0] <- t0i.d[1]
       t1i.d[1] <- t1i.d[1]
 
  D) and then, the data are in the same scramble than after 
     last DIT pass. Do first DIT pass.
 
 
  In this SSE2 version, bj's will be stored as float SSE2. because we will
  ask about (bj >= c), we will actually ask for (bj >= (c - 0.49) ) to avoid 
  precision problems. We could round bj's every UPDATE passes.
 
*/


/* to avoid problems while writing code */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#include "glucas.h"
#include "gsetup.h"

#if defined (Y_USE_SSE2)
#include "ygensse2.h"
#include "fft5sse2.h"
#include "ynormsse2.h"


#define get_pads_5          \
  pad2 = (pad << 1);        \
  pad4 = (pad << 2);        \
  pad3 = pad + pad2;        \
  pd0 = x;                  \
  prefetch_data (pd0, 0);   \
  pd1 = x + addr(pad << 1); \
  prefetch_data (pd1, 0);   \
  pd2 = x + addr(pad2 << 1);\
  prefetch_data (pd2, 0);   \
  pd3 = x + addr(pad3 << 1);\
  prefetch_data (pd3, 0);   \
  pd4 = x + addr(pad4 << 1);\
  prefetch_data (pd4, 0);


#define init_bjs_5                       \
  bj0 = N;                               \
  bj1 = ( b % 5) << (Y_K + 1);           \
  bj2 = (( 2 * b) % 5) << (Y_K + 1);     \
  bj3 = (( 3 * b) % 5) << (Y_K + 1);     \
  bj4 = (( 4 * b) % 5) << (Y_K + 1);     \
  /* now it passes to SSE2 vars */       \
  Y_MM_SET_PD( bj02, bj0,  bj2);         \
  Y_MM_SET_PD( bj13, bj1,  bj3);         \
  Y_MM_SET_PD( bj4x, bj4,    0);

#define round_bjs_5                          \
  Y_MM_TWO_RINT_PD( bj02, bj13, bj02, bj13); \
  Y_MM_RINT_PD ( bj4x, bj4x);

#define  sse2_load_dwt_5(_i)                                            \
  {                                                                     \
    size_t j0, j1;                                                      \
    j0 = addr (_i);                                                     \
    j1 = addr (((pad2 >> (SHIFT_UPDATE - 1)) + _i));                    \
    Y_MM_SET_PD (ttmp02, two_to_minusphi[j0], two_to_minusphi[j1]);     \
    Y_MM_SET_PD (ttp02, two_to_phi[j0], two_to_phi[j1]);                \
    prefetch_data_trig (two_to_minusphi, j0 + Y_CACHE_LINE);            \
    prefetch_data_trig (two_to_minusphi, j1 + Y_CACHE_LINE);            \
    prefetch_data_trig (two_to_phi, j0 + Y_CACHE_LINE);                 \
    prefetch_data_trig (two_to_phi, j1 + Y_CACHE_LINE);                 \
    j0 = addr (((pad2 >> SHIFT_UPDATE) + _i));                          \
    j1 = addr (((pad3 >> (SHIFT_UPDATE - 1)) + _i));                    \
    Y_MM_SET_PD (ttmp13, two_to_minusphi[j0], two_to_minusphi[j1]);     \
    Y_MM_SET_PD (ttp13, two_to_phi[j0], two_to_phi[j1]);                \
    prefetch_data_trig (two_to_minusphi, j0 + Y_CACHE_LINE);            \
    prefetch_data_trig (two_to_minusphi, j1 + Y_CACHE_LINE);            \
    prefetch_data_trig (two_to_phi, j0 + Y_CACHE_LINE);                 \
    prefetch_data_trig (two_to_phi, j1 + Y_CACHE_LINE);                 \
    j0 = addr (((pad4 >> (SHIFT_UPDATE - 1)) + _i));                    \
    Y_MM_SET_PD (ttmp4x, two_to_minusphi[j0], 0.0);                     \
    Y_MM_SET_PD (ttp4x, two_to_phi[j0], 0.0);                           \
    prefetch_data_trig (two_to_minusphi, j0 + Y_CACHE_LINE);            \
    prefetch_data_trig (two_to_phi, j0 + Y_CACHE_LINE);                 \
  }

/* Interleaves before and after carry and norm phase */

#define inter_chain_5              \
    Y_MM_MOV_INTER_PAIR(t0r, t2r); \
    Y_MM_MOV_INTER_PAIR(t0i, t2i); \
    Y_MM_MOV_INTER_PAIR(t1r, t3r); \
    Y_MM_MOV_INTER_PAIR(t1i, t3i); \
    Y_MM_MOV_INTER_PAIR(t4r, t5r); \
    Y_MM_MOV_INTER_PAIR(t4i, t5i);

/* get the basic twiddle from memory and computes its powers */
#ifdef Y_MINIMUM
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# define init_twiddle_factors_sse2_5( _pd )      \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, _pd);\
   sse2_square(tw2, tw1);                        \
   sse2_square(tw4, tw2);                        \
   sse2_mul(tw3, tw2, tw1);                      \
   prefetch_p_trig( _pd );

#  define get_twiddle_factors_sse2_5                     \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);\
   sse2_square(tw2, tw1);                                \
   sse2_square(tw4, tw2);                                \
   sse2_mul(tw3, tw2, tw1);                              \
   prefetch_p_trig(px);


#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 8

# define init_twiddle_factors_sse2_5( _pd )           \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, _pd);     \
   Y_MM_LOAD_INTER_ONE_HALF(tw2##r, tw2##i, _pd + 2); \
   Y_MM_LOAD_INTER_ONE_HALF(tw3##r, tw3##i, _pd + 4); \
   Y_MM_LOAD_INTER_ONE_HALF(tw4##r, tw4##i, _pd + 6);

#  define get_twiddle_factors_sse2_5                               \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, px + 4, px + Y_STEP + 4);  \
   Y_MM_LOAD_INTER_PAIR(tw4##r, tw4##i, px + 6, px + Y_STEP + 6);

#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 4

# define init_twiddle_factors_sse2_5(_pd)        \
   Y_MM_LOAD_INTER_ONE_HALF(tw1r, tw1i, _pd);    \
   Y_MM_LOAD_INTER_ONE_HALF(tw3r, tw3i, _pd + 2);\
   sse2_divmul(tw2, tw4, tw3, tw1);              \
   prefetch_data_trig (_pd, Y_STEP);             \

#  define get_twiddle_factors_sse2_5                            \
   Y_MM_LOAD_INTER_PAIR (tw1r, tw1i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR (tw3r, tw3i, px + 2, px + Y_STEP + 2);  \
   sse2_divmul (tw2, tw4, tw3, tw1);                            \
   prefetch_data_trig (px, 2 * Y_STEP);                         \
   prefetch_data_trig (px, 3 * Y_STEP);

#endif

#if !defined(ALL_INTERLACED)

# define radixmm_5_twd_last_dit(_j)                          \
  sse2_data_to_local_pp_inter (t0, _j, pd0);                 \
  sse2_load_inter_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_inter_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  Y_MM_SETZERO_PD( t5r);                                     \
  Y_MM_SETZERO_PD( t5i);                                     \
  sse2_fft5 (B, t0, t1, t2, t3, t4);                         \
  inter_chain_5;



# define radixmm_5_first_dif(_j)                                               \
  inter_chain_5;                                                               \
  sse2_fft5_store_inter_p (_j, pd0, pd1, pd2, pd3, pd4, F, t0, t1, t2, t3, t4);



# define radixmm_5_first_dif_add(_d)          \
  inter_chain_5;                              \
  sse2_fft5 (F, t0, t1, t2, t3, t4);          \
  sse2_local_to_data_inter_add (_d,    0, t0);\
  sse2_local_to_data_inter_add (_d,  pad, t1);\
  sse2_local_to_data_inter_add (_d, pad2, t2);\
  sse2_local_to_data_inter_add (_d, pad3, t3);\
  sse2_local_to_data_inter_add (_d, pad4, t4);


# define radixmm_5_first_dif_add_shift(_d, _ofs)\
  sse2_fft5 (F, t0, t1, t2, t3, t4);           \
  sse2_local_to_data_inter_add (_d, _ofs, t0); \
  sse2_local_to_data_inter_add (_d,  pad, t1); \
  sse2_local_to_data_inter_add (_d, pad2, t2); \
  sse2_local_to_data_inter_add (_d, pad3, t3); \
  sse2_local_to_data_inter_add (_d, pad4, t4);

#else

# define radixmm_5_twd_last_dit(_j)                    \
  sse2_data_to_local_pp (t0, _j, pd0);                 \
  sse2_load_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  Y_MM_SETZERO_PD( t5r);                               \
  Y_MM_SETZERO_PD( t5i);                               \
  sse2_fft5 (B, t0, t1, t2, t3, t4);                   \
  inter_chain_5;



# define radixmm_5_first_dif(_j)                                        \
  inter_chain_5;                                                        \
  sse2_fft5_store_p (_j, pd0, pd1, pd2, pd3, pd4, F, t0, t1, t2, t3, t4);



# define radixmm_5_first_dif_add(_d)    \
  inter_chain_5;                        \
  sse2_fft5 (F, t0, t1, t2, t3, t4);    \
  sse2_local_to_data_add (_d,    0, t0);\
  sse2_local_to_data_add (_d,  pad, t1);\
  sse2_local_to_data_add (_d, pad2, t2);\
  sse2_local_to_data_add (_d, pad3, t3);\
  sse2_local_to_data_add (_d, pad4, t4);


# define radixmm_5_first_dif_add_shift(_d, _ofs)\
  sse2_fft5 (F, t0, t1, t2, t3, t4);           \
  sse2_local_to_data_add (_d, _ofs, t0);       \
  sse2_local_to_data_add (_d,  pad, t1);       \
  sse2_local_to_data_add (_d, pad2, t2);       \
  sse2_local_to_data_add (_d, pad3, t3);       \
  sse2_local_to_data_add (_d, pad4, t4);


#endif

#define sse2_sum_5( _Sum )        \
{                                 \
  Y__M128D aux0, aux1, aux2, aux3;\
  BIG_DOUBLE xsum[2];             \
  Y_MM_ADD_PD (aux0, t0r, t0i);   \
  Y_MM_ADD_PD (aux1, t1r, t1i);   \
  Y_MM_ADD_PD (aux2, t2r, t2i);   \
  Y_MM_ADD_PD (aux3, t3r, t3i);   \
  Y_MM_ADD_PD (aux0, aux0, t4r);  \
  Y_MM_ADD_PD (aux1, aux1, t4i);  \
  Y_MM_ADD_PD (aux2, aux2, t5r);  \
  Y_MM_ADD_PD (aux3, aux3, t5i);  \
  Y_MM_ADD_PD (aux0, aux0, aux1); \
  Y_MM_ADD_PD (aux2, aux2, aux3); \
  Y_MM_ADD_PD (aux0, aux0, aux2); \
  Y_MM_STORE_PD (xsum, aux0);     \
  _Sum += (xsum[0] + xsum[1]);    \
}


void substract_two_5_sse2(BIG_DOUBLE *x, UL N)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i,
  t5r, t5i;
  Y__M128D Y_XTWO[2];
  double Y_ALIGNED(16) xtwo[2];
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4, ofs, ni, who, ib, a;

  /* the bit to substract */
  ib = Y_SBIT + 1;
  if(ib == Y_EXPONENT)
    ib = 0;

  /* The element ni where bit 'ib' is */
  ni = floor((BIG_DOUBLE)ib / Y_XBITS);

  /* now ib is the bit within x[ni] */
  ib -= (UL)ceil(Y_XBITS * ni);

  /* The offset, in size of complex, of group

     the difference in this SSE2 version is now the ofset is forced to be
     even because of two complex are grouped in a pair (tr, ti)
  */
  ofs = (ni % (N / Y_PLAN[0])) >> 1;
  ofs -= (ofs & 1U);
  prefetch_data ( x, addr(ofs << 1));

  /* prepare pointers */
  pad2 = (pad << 1);
  pad4 = (pad << 2);
  pad3 = pad + pad2;
  pad += ofs;
  pad2 += ofs;
  prefetch_data ( x, addr(pad << 1));
  prefetch_data ( x, addr(pad2 << 1));
  pad3 += ofs;
  pad4 += ofs;
  prefetch_data ( x, addr(pad3 << 1));
  prefetch_data ( x, addr(pad4 << 1));

  /* What element of group is x[ni] */
  who = (ni * Y_PLAN[0]) / N;

  /* The shifted value of minustwo */
  xtwo[0] = -(BIG_DOUBLE)(1U << ib);

  /*computes (N - ni*b % N). We suposse up to 256 million exponents */
  if(ni)
    {
#if ULONG_MAX > 0xFFFFFFFF
      /* 64 bits integer machine */
      a = N - (ni * b) % N;
#else
      /* 32 bits integer machine */
      modmul_32 (a, b, ni, N);
      a = N - a;
#endif

    }
  else
    a = 0;

  /* Multiply by two_to_phi() */
  xtwo[0] *= exp(a * M_LN2 / N);
  xtwo[1] = 0.0;
  Y_MM_LOAD_PD( Y_XTWO[0], xtwo);

  /* Sumcheck */
#if defined(SUM_CHECK)

  SumIn += (xtwo[0] + xtwo[1]);
#endif

  xtwo[1] = xtwo[0];
  xtwo[0] = 0.0;
  Y_MM_LOAD_PD( Y_XTWO[1], xtwo);

  /* do some util work waiting the end of the above long latency ins*/
  Y_MM_SETZERO_PD( t0r );
  Y_MM_SETZERO_PD( t0i );
  Y_MM_SETZERO_PD( t1r );
  Y_MM_SETZERO_PD( t1i );
  Y_MM_SETZERO_PD( t2r );
  Y_MM_SETZERO_PD( t2i );
  Y_MM_SETZERO_PD( t3r );
  Y_MM_SETZERO_PD( t3i );
  Y_MM_SETZERO_PD( t4r );
  Y_MM_SETZERO_PD( t4i );
  Y_MM_SETZERO_PD( t5r );
  Y_MM_SETZERO_PD( t5i );

  switch (who)
    {
    case 0:
      switch (ni & 3U)
        {
        case 0:
          Y_MM_MOV_PD ( t0r , Y_XTWO[0]);
          break;
        case 1:
          Y_MM_MOV_PD ( t0i , Y_XTWO[0]);
          break;
        case 2:
          Y_MM_MOV_PD ( t0r , Y_XTWO[1]);
          break;
        case 3:
          Y_MM_MOV_PD ( t0i , Y_XTWO[1]);
          break;
        }
      break;
    case 1:
      switch (ni & 3U)
        {
        case 0:
          Y_MM_MOV_PD ( t1r , Y_XTWO[0]);
          break;
        case 1:
          Y_MM_MOV_PD ( t1i , Y_XTWO[0]);
          break;
        case 2:
          Y_MM_MOV_PD ( t1r , Y_XTWO[1]);
          break;
        case 3:
          Y_MM_MOV_PD ( t1i , Y_XTWO[1]);
          break;
        }
      break;
    case 2:
      switch (ni & 3U)
        {
        case 0:
          Y_MM_MOV_PD ( t2r , Y_XTWO[0]);
          break;
        case 1:
          Y_MM_MOV_PD ( t2i , Y_XTWO[0]);
          break;
        case 2:
          Y_MM_MOV_PD ( t2r , Y_XTWO[1]);
          break;
        case 3:
          Y_MM_MOV_PD ( t2i , Y_XTWO[1]);
          break;
        }
      break;
    case 3:
      switch (ni & 3U)
        {
        case 0:
          Y_MM_MOV_PD ( t3r , Y_XTWO[0]);
          break;
        case 1:
          Y_MM_MOV_PD ( t3i , Y_XTWO[0]);
          break;
        case 2:
          Y_MM_MOV_PD ( t3r , Y_XTWO[1]);
          break;
        case 3:
          Y_MM_MOV_PD ( t3i , Y_XTWO[1]);
          break;
        }
      break;
    case 4:
      switch (ni & 3U)
        {
        case 0:
          Y_MM_MOV_PD ( t4r , Y_XTWO[0]);
          break;
        case 1:
          Y_MM_MOV_PD ( t4i , Y_XTWO[0]);
          break;
        case 2:
          Y_MM_MOV_PD ( t4r , Y_XTWO[1]);
          break;
        case 3:
          Y_MM_MOV_PD ( t4i , Y_XTWO[1]);
          break;
        }
      break;
    }

  /* make first DIF and add to array */
  radixmm_5_first_dif_add_shift( x, ofs);

  /* New Y_SBIT and Y_SBJ*/
  Y_SBIT += Y_SBIT;
  if(Y_SBIT >= Y_EXPONENT)
    Y_SBIT -= Y_EXPONENT;
}

#if !defined(_OPENMP) && !defined(_SUNMP) && !defined(_PTHREADS)

void dit_carry_norm_dif_5_sse2 (BIG_DOUBLE *x, UL N ,UL err_flag)
{
  Y__M128D ttp02, ttp13, ttp4x, ttmp02, ttmp13, ttmp4x;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i;
  Y__M128D bj02, bj13, bj4x;
  Y__M128D carry02, carry13, carry4x, maxerr;
#if defined(SUM_CHECK)

  Y__M128D summ;
#endif

  UL bj0, bj1, bj2, bj3, bj4;
  y_ptr px, pd0, pd1, pd2, pd3, pd4;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4;
  UL i, j, k, l, ll;

  /* init carries */

  px = Y_TWDB[Y_NRADICES - 2];

  prefetch_data_trig( px, 0);

  Y_MM_SETZERO_PD (carry02);
  Y_MM_SETZERO_PD (carry13);
  Y_MM_SETZERO_PD (carry4x);
  Y_MM_SETZERO_PD (maxerr);
#if defined(SUM_CHECK)

  Y_MM_SETZERO_PD (summ);
#endif

  if(Y_SBIT == 0)
    Y_MM_SET_PD (carry02, -2.0, 0.0);

  get_pads_5;

  init_bjs_5;

  init_twiddle_factors_sse2_5 (px);

  px += Y_STEP;

  for (i = 0, j = 0; j < pad2; j += UPDATE, i++)
    {
      sse2_load_dwt_5 (i);
      /*round_bjs_5;*/
      for (k = 0; k < UPDATE; k += 4, px += 2 * Y_STEP)
        {
          l = (j + k) >> 1;
          ll = addr(j+k);

          radixmm_5_twd_last_dit( ll );

#if defined(SUM_CHECK)
          /*sse2_sum_5 (SumOut);*/
#endif

          if(err_flag)
            {
#if defined(Y_VECTORIZE)
              /* Intel 7.x compiler crash with this code */
              cplx_carry_norm_check_sse2_vector(t0r, t1r, 02, 13);
              cplx_carry_norm_check_sse2_vector(t0i, t1i, 02, 13);
              cplx_carry_norm_check_sse2_vector(t2r, t3r, 02, 13);
              cplx_carry_norm_check_sse2_vector(t2i, t3i, 02, 13);
#else

              cplx_carry_norm_check_sse2(t0r, 02);
              cplx_carry_norm_check_sse2(t0i, 02);
              cplx_carry_norm_check_sse2(t2r, 02);
              cplx_carry_norm_check_sse2(t2i, 02);
              cplx_carry_norm_check_sse2(t1r, 13);
              cplx_carry_norm_check_sse2(t1i, 13);
              cplx_carry_norm_check_sse2(t3r, 13);
              cplx_carry_norm_check_sse2(t3i, 13);
#endif /* Y_VECTORIZE */

              cplx_carry_norm_check_sse2_sd(t4r, 4x);
              cplx_carry_norm_check_sse2_sd(t4i, 4x);
              cplx_carry_norm_check_sse2_sd(t5r, 4x);
              cplx_carry_norm_check_sse2_last_odd_sd(t5i, 4x);
            }
          else
            {
#if defined(Y_VECTORIZE)
              /* Intel 7.x compiler crash with this code */
              cplx_carry_norm_sse2_vector(t0r, t1r, 02, 13);
              cplx_carry_norm_sse2_vector(t0i, t1i, 02, 13);
              cplx_carry_norm_sse2_vector(t2r, t3r, 02, 13);
              cplx_carry_norm_sse2_vector(t2i, t3i, 02, 13);
#else

              cplx_carry_norm_sse2(t0r, 02);
              cplx_carry_norm_sse2(t0i, 02);
              cplx_carry_norm_sse2(t2r, 02);
              cplx_carry_norm_sse2(t2i, 02);
              cplx_carry_norm_sse2(t1r, 13);
              cplx_carry_norm_sse2(t1i, 13);
              cplx_carry_norm_sse2(t3r, 13);
              cplx_carry_norm_sse2(t3i, 13);
#endif /* Y_VECTORIZE */

              cplx_carry_norm_sse2_sd(t4r, 4x);
              cplx_carry_norm_sse2_sd(t4i, 4x);
              cplx_carry_norm_sse2_sd(t5r, 4x);
              cplx_carry_norm_sse2_last_odd_sd(t5i, 4x);
            }
          radixmm_5_first_dif( ll );

#if defined(SUM_CHECK)

          sse2_sum (SumIn, t0r, t0i);
#endif

          get_twiddle_factors_sse2_5;
        }
    }

  /* adjust the last carries */
  {
    double Y_ALIGNED(16) aux[4], aux0;
    Y_MM_STORE_PD( aux, carry13);
    Y_MM_STORE_PD( aux + 2, carry4x);
    aux0   = aux[1];
    aux[1] =    0.0;
    aux[3] = aux[0];
    aux[0] =   aux0;
    Y_MM_LOAD_PD( t4r, aux);
    Y_MM_LOAD_PD( t0r, aux + 2);
  }
  Y_MM_MOV_PD( t1r, carry02);
  Y_MM_SETZERO_PD( carry13);
  Y_MM_SETZERO_PD( carry02);
  Y_MM_SETZERO_PD( carry4x);

  Y_MM_SETZERO_PD ( t0i );
  Y_MM_SETZERO_PD ( t2r );
  Y_MM_SETZERO_PD ( t2i );
  Y_MM_SETZERO_PD ( t1i );
  Y_MM_SETZERO_PD ( t3r );
  Y_MM_SETZERO_PD ( t3i );
  Y_MM_SETZERO_PD ( t4i );
  Y_MM_SETZERO_PD ( t5r );
  Y_MM_SETZERO_PD ( t5i );

  init_bjs_5;

  j = 0;
  sse2_load_dwt_5(j);

  cplx_carry_sse2 (t0r, 02);
  cplx_carry_sse2 (t0i, 02);
  cplx_carry_sse2 (t2r, 02);
  cplx_carry_sse2 (t2i, 02);
  prefetch_data (pd0, 0);
  prefetch_data (pd2, 0);
  cplx_carry_sse2 (t1r, 13);
  cplx_carry_sse2 (t1i, 13);
  cplx_carry_sse2 (t3r, 13);
  cplx_carry_sse2 (t3i, 13);
  prefetch_data (pd1, 0);
  prefetch_data (pd3, 0);
  cplx_carry_sse2 (t4r, 4x);
  cplx_carry_sse2 (t4i, 4x);
  cplx_carry_sse2 (t5r, 4x);
  cplx_carry_sse2 (t5i, 4x);
  prefetch_data (pd4, 0);

  /* add carry signal */
  radixmm_5_first_dif_add (x);

#if defined(SUM_CHECK)

  sse2_sum (SumIn, t0r, t0i);
  sse2_sum0 (SumOut, summ);
#endif

  if(err_flag)
    {
      double Y_ALIGNED(16) xmax[2];
      Y_MM_STORE_PD (xmax, maxerr);
      Err = (xmax[0] > xmax[1]) ? xmax[0] : xmax[1];
    }
  else
    Err = 0.0;

  if(Y_SBIT)
    substract_two_5_sse2 (x, N);
}


#else

/***********************************************************************/
/* These are the Multithread adds needed                               */
/***********************************************************************/
/*
To work with mutithread environment, data is cut into independent slims.
As example, assume we are working with N=256 complex data size and plan
is [4,....]. Also assume we want four threads.

The access patern is
c[0], c[64], c[128], c[192] work with it and store new values in same
place
c[1], c[65], c[129], c[193] in carry and normalization we work with them
as r[]
....................
c[63], c[127], c[191], c[255]


Then there are 64 groups of 8 complex data. Divide them in four slims and
assign a slim per thread

Thread 0 :
c[0],c[64]...
c[1],c[65]...
..........
c[15],c[79]...

Thread 1 :
c[16],c[80]...
c[17],c[81]...
...........
c[31],c[47]

Thread 2 :
c[32],c[96]...
..........

Thread 3:
c[48],c[112]...
..........
c[63],c[127]....

It will require, of course, a last carry adjust among threads.

_ofs is the offset, in complex size,  of every thread nthr
_ofs = _n0 * nthr;  0 <= nthr < max_num_of_threads
where _n0 is the basic length of a thread, and MUST be a
power of two
*/
# define get_pads_5_mp(_ofs)         \
  pad2 = (pad << 1);                 \
  pad4 = (pad << 2);                 \
  pad3 = pad + pad2;                 \
  pd0 = x + addr((_ofs << 1));       \
  prefetch_data( pd0, 0);            \
  pd1 = x + addr((pad  + _ofs) << 1);\
  prefetch_data( pd1, 0);            \
  pd2 = x + addr((pad2 + _ofs) << 1);\
  prefetch_data( pd2, 0);            \
  pd3 = x + addr((pad3 + _ofs) << 1);\
  prefetch_data( pd3, 0);            \
  pd4 = x + addr((pad4 + _ofs) << 1);\
  prefetch_data( pd4, 0);

/*
To compute the initial bj of every thread we use:

b(i,nthr)= (((nthr + i*(pad/n0))*b)% (N/(2*n0)) * 2 * n0;
where
pad/n0 must be integer and 0< (pad/n0) < max_num_of_threads
n0 must be a power of two.
It will be computed at first time. so we avoid to compute it every
iteration
*/

# define init_bjs_sse2_5_mp        \
  Y_MM_SET_PD( bj02, bp[0], bp[2]);\
  Y_MM_SET_PD( bj13, bp[1], bp[3]);\
  Y_MM_SET_PD( bj4x, bp[4], N);

/*
The twidle factors now are not (1.0,0.0) at the begining
so it should be computed

inc=(Y_POWERS[Y_PLAN[0]-1])<<1;
px = Y_TWDB[Y_NRADICES - 2] + _ofs * inc - 2;

It could be computed at a first iteration, so we save
some work
*/

# if defined(Y_MINIMUM)
#  define init_twiddle_factors_sse2_5_mp(_px)              \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, _px, _px + Y_STEP);\
   sse2_square (tw2, tw1);                                 \
   sse2_square (tw4, tw2);                                 \
   sse2_mul (tw3, tw2, tw1);                               \
   prefetch_p_trig (_px);

# elif defined(Y_MAXIMUM)

#  define init_twiddle_factors_sse2_5_mp(_px)                        \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, _px, _px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, _px + 2, _px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, _px + 4, _px + Y_STEP + 4);  \
   Y_MM_LOAD_INTER_PAIR(tw4##r, tw4##i, _px + 6, _px + Y_STEP + 6);  \
   prefetch_data_trig(_px, Y_STEP);                                  \
   prefetch_data_trig(_px, Y_STEP + Y_CACHE_LINE);

# else

#  define init_twiddle_factors_sse2_5_mp(_px)                   \
   Y_MM_LOAD_INTER_PAIR (tw1r, tw1i, _px, _px + Y_STEP);        \
   Y_MM_LOAD_INTER_PAIR (tw3r, tw3i, _px + 2, _px + Y_STEP + 2);\
   sse2_divmul (tw2, tw4, tw3, tw1);                            \
   prefetch_data_trig (_px, 2 * Y_STEP);                        \
   prefetch_data_trig (_px, 3 * Y_STEP);

# endif

# define save_carries_5_mp(_pcarr)          \
        Y_MM_STORE_PD( _pcarr , carry02);   \
        Y_MM_STORE_PD( _pcarr + 2, carry13);\
        Y_MM_STORE_PD( _pcarr + 4, carry4x);\
        _pcarr[5] = 0.0;

# if !defined(ALL_INTERLACED)
#  define radixmm_5_first_dif_add_mp        \
  inter_chain_5;                            \
  sse2_fft5 (F, t0, t1, t2, t3, t4);        \
  sse2_local_to_data_inter_add (pd0, 0, t0);\
  sse2_local_to_data_inter_add (pd1, 0, t1);\
  sse2_local_to_data_inter_add (pd2, 0, t2);\
  sse2_local_to_data_inter_add (pd3, 0, t3);\
  sse2_local_to_data_inter_add (pd4, 0, t4);

# else

#  define radixmm_5_first_dif_add_mp  \
  inter_chain_5;                      \
  sse2_fft5 (F, t0, t1, t2, t3, t4);  \
  sse2_local_to_data_add (pd0, 0, t0);\
  sse2_local_to_data_add (pd1, 0, t1);\
  sse2_local_to_data_add (pd2, 0, t2);\
  sse2_local_to_data_add (pd3, 0, t3);\
  sse2_local_to_data_add (pd4, 0, t4);

# endif

BIG_DOUBLE dit_carry_norm_dif_5_mp_sse2(y_ptr x, y_ptr ptx, y_ptr pcarr,
                                        y_ptr bp, UL N ,UL n0, UL nthr, UL err_flag)
{
  Y__M128D ttp02, ttp13, ttp4x, ttmp02, ttmp13, ttmp4x;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i;
  Y__M128D bj02, bj13, bj4x;
  Y__M128D carry02, carry13, carry4x, maxerr ;
  y_ptr pd0, pd1, pd2, pd3, pd4, px;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4;
  UL i, j, k, l, ll, ofs, lim;

  ASSERT_ALIGNED_DOUBLE();

  ofs = n0 * nthr;

# if defined(SUM_CHECK)

  Y_SUMOUT[nthr] = 0.0;
  Y_SUMIN[nthr] = 0.0;
# endif

  prefetch_data_trig( ptx, 0);

  /* Clean the carries */
  Y_MM_SETZERO_PD( carry02 );
  Y_MM_SETZERO_PD( carry13 );
  Y_MM_SETZERO_PD( carry4x );
  Y_MM_SETZERO_PD( maxerr  );

  /* get the pads */
  get_pads_5;

  /* get the bjs */
  init_bjs_sse2_5_mp;

  if(nthr)
    {
      /*
      init of other than first thread
      */
      init_twiddle_factors_sse2_5_mp (ptx);
      px = ptx + 2 * Y_STEP;
    }
  else
    {
      /*
      First thread. Init steps are the same than
      non-multithreaded case
      */
      if(Y_SBIT == 0)
        Y_MM_SET_PD(carry02, -2.0, 0.0);
      init_twiddle_factors_sse2_5 (ptx);
      px = ptx + Y_STEP;
    }

  if(nthr < (Y_NTHREADS - 1))
    lim = ((ofs + n0) << 1);
  else
    lim = (Y_LRIGHT[1] << 1);

  for (i = (ofs << 1) / UPDATE, j = (ofs << 1); j < lim ; j += UPDATE, i++)
    {
      sse2_load_dwt_5(i);

      for (k = 0; k < UPDATE; k += 4, px += (2 * Y_STEP))
        {
          l = (j + k) >> 1;
          ll = addr(j + k);
          radixmm_5_twd_last_dit( ll );

# ifdef SUM_CHECK

          sse2_sum_5 (Y_SUMOUT[nthr]);
# endif

          if(err_flag)
            {
#if !defined(__ICC) && defined(Y_VECTORIZE)
              /* Intel 7.x compiler crash with this code */
              cplx_carry_norm_check_sse2_vector(t0r, t1r, 02, 13);
              cplx_carry_norm_check_sse2_vector(t0i, t1i, 02, 13);
              cplx_carry_norm_check_sse2_vector(t2r, t3r, 02, 13);
              cplx_carry_norm_check_sse2_vector(t2i, t3i, 02, 13);
#else

cplx_carry_norm_check_sse2(t0r, 02);
cplx_carry_norm_check_sse2(t0i, 02);
cplx_carry_norm_check_sse2(t2r, 02);
cplx_carry_norm_check_sse2(t2i, 02);
cplx_carry_norm_check_sse2(t1r, 13);
cplx_carry_norm_check_sse2(t1i, 13);
cplx_carry_norm_check_sse2(t3r, 13);
cplx_carry_norm_check_sse2(t3i, 13);
#endif /* __ICC */

              cplx_carry_norm_check_sse2_sd(t4r, 4x);
              cplx_carry_norm_check_sse2_sd(t4i, 4x);
              cplx_carry_norm_check_sse2_sd(t5r, 4x);
              cplx_carry_norm_check_sse2_last_odd_sd(t5i, 4x);
            }
          else
            {
#if !defined(__ICC) && defined(Y_VECTORIZE)
              /* Intel 7.x compiler crash with this code */
              cplx_carry_norm_sse2_vector(t0r, t1r, 02, 13);
              cplx_carry_norm_sse2_vector(t0i, t1i, 02, 13);
              cplx_carry_norm_sse2_vector(t2r, t3r, 02, 13);
              cplx_carry_norm_sse2_vector(t2i, t3i, 02, 13);
#else

cplx_carry_norm_sse2(t0r, 02);
cplx_carry_norm_sse2(t0i, 02);
cplx_carry_norm_sse2(t2r, 02);
cplx_carry_norm_sse2(t2i, 02);
cplx_carry_norm_sse2(t1r, 13);
cplx_carry_norm_sse2(t1i, 13);
cplx_carry_norm_sse2(t3r, 13);
cplx_carry_norm_sse2(t3i, 13);
#endif /* __ICC */

              cplx_carry_norm_sse2_sd(t4r, 4x);
              cplx_carry_norm_sse2_sd(t4i, 4x);
              cplx_carry_norm_sse2_sd(t5r, 4x);
              cplx_carry_norm_sse2_last_odd_sd(t5i, 4x);
            }
          radixmm_5_first_dif( ll );

#if defined(SUM_CHECK)

          sse2_sum (Y_SUMIN[nthr], t0r, t0i);
#endif

          if( (k + j + 4 ) < lim )
            get_twiddle_factors_sse2_5;
        }
    }
  /* save the carries */
  save_carries_5_mp(pcarr);

  /* return the error */
  if(err_flag)
    {
      double Y_ALIGNED(16) xmax[2],aux;
      Y_MM_STORE_PD (xmax, maxerr);
      aux = (xmax[0] > xmax[1]) ? xmax[0] : xmax[1];
      return aux;
    }
  else
    {
      return 0.0;
    }
}


void dit_carry_norm_dif_5_last_carries_mp_sse2 (y_ptr x, y_ptr bp, UL N ,UL n0, int nthr)
{
  Y__M128D ttp02, ttp13, ttp4x, ttmp02, ttmp13, ttmp4x;
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i;
  Y__M128D bj02, bj13, bj4x;
  Y__M128D carry02, carry13, carry4x;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4;
  UL  j, ofs;

  ASSERT_ALIGNED_DOUBLE();

  Y_MM_SETZERO_PD( carry02 );
  Y_MM_SETZERO_PD( carry13 );
  Y_MM_SETZERO_PD( carry4x );
  Y_MM_SETZERO_PD ( t0i );
  Y_MM_SETZERO_PD ( t2r );
  Y_MM_SETZERO_PD ( t2i );
  Y_MM_SETZERO_PD ( t1i );
  Y_MM_SETZERO_PD ( t3r );
  Y_MM_SETZERO_PD ( t3i );
  Y_MM_SETZERO_PD ( t4i );
  Y_MM_SETZERO_PD ( t5r );
  Y_MM_SETZERO_PD ( t5i );

  /* Init the data */
  if(nthr)
    {
      Y_MM_LOAD_PD( t0r, Y_CARRIES[nthr - 1] + 0);
      Y_MM_LOAD_PD( t1r, Y_CARRIES[nthr - 1] + 2);
      Y_MM_LOAD_PD( t4r, Y_CARRIES[nthr - 1] + 4);
    }
  else
    {
      double Y_ALIGNED(16) aux[4];
      aux[0] = Y_CARRIES[Y_NTHREADS - 1][3];
      aux[1] = 0.0;
      aux[2] = Y_CARRIES[Y_NTHREADS - 1][4];
      aux[3] = Y_CARRIES[Y_NTHREADS - 1][2];
      Y_MM_LOAD_PD( t4r, aux);
      Y_MM_LOAD_PD( t0r, aux + 2);
      Y_MM_LOAD_PD( t1r, Y_CARRIES[Y_NTHREADS - 1] + 0);
    }

  init_bjs_sse2_5_mp;

  ofs = (n0 * nthr);

  j = (ofs << 1) / UPDATE;

  get_pads_5_mp (ofs);

  sse2_load_dwt_5 (j);

  cplx_carry_sse2 (t0r, 02);
  cplx_carry_sse2 (t0i, 02);
  cplx_carry_sse2 (t2r, 02);
  cplx_carry_sse2 (t2i, 02);
  cplx_carry_sse2 (t1r, 13);
  cplx_carry_sse2 (t1i, 13);
  cplx_carry_sse2 (t3r, 13);
  cplx_carry_sse2 (t3i, 13);
  cplx_carry_sse2 (t4r, 4x);
  cplx_carry_sse2 (t4i, 4x);
  cplx_carry_sse2 (t5r, 4x);
  cplx_carry_sse2 (t5i, 4x);

  /* add carry signal */

  radixmm_5_first_dif_add_mp;

#if defined(SUM_CHECK)

  sse2_sum (Y_SUMIN[nthr], t0r, t0i);
#endif
}

void dit_carry_norm_dif_5_sse2(BIG_DOUBLE *x, UL N ,UL err_flag)
{
  int i, n0;
  n0 = Y_LRIGHT[1] / Y_NTHREADS;

  /* make n0 divisble by Y_UPDATE */
  n0 = ((n0 >> (Y_SHIFT_UPDATE - 1)) << (Y_SHIFT_UPDATE - 1));


  /* The multithread loop */
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(i)
#  pragma MP taskloop shared(x, Y_ERRS, Y_TWX, Y_CARRIES, Y_XBJS, N, n0, err_flag)
# endif

  for(i = 0; i < Y_NTHREADS; i++)
    {
      Y_ERRS[i] = dit_carry_norm_dif_5_mp_sse2 ( x, Y_TWX[i], Y_CARRIES[i], Y_XBJS[i],
                  N, n0, i, err_flag);
    }

  /* maxerr */
  Err = 0.0;
  if(err_flag)
    {
      for(i = 0; i < Y_NTHREADS; i++)
        if(Y_ERRS[i] > Err)
          Err = Y_ERRS[i];
    }

  /* And finally, the last carry propagation */

# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(x, Y_XBJS, N, n0)
# endif
  for(i = 0; i < Y_NTHREADS; i++)
    {
      dit_carry_norm_dif_5_last_carries_mp_sse2( x, Y_XBJS[i],
          N, n0, i);
    }

  /* Substract two if shifts*/
  if(Y_SBIT)
    substract_two_5_sse2( x, N);

# if defined(SUM_CHECK)

  for (i = 0; i < Y_NTHREADS; i++)
    {
      SumOut += Y_SUMOUT[i];
      SumIn += Y_SUMIN[i];
    }
# endif

}

# endif /* openmp .... */

#else /* Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void void_carry_norm_5_sse2 (void)
{
  printf("The routine void_carry_norm_5_sse2 should never be called\n");
  exit(EXIT_FAILURE);
}

#endif /* Y_USE_SSE2 */
/*$Id$*/








