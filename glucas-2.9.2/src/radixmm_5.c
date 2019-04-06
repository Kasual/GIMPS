/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2003-2006  Guillermo Ballester Valor
 
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
/*
   About read/write mode 
   if y_rmode == 0 then it reads from normal (legacy mode)
      y_rmode == 1 then it reads from fast scramble mode
   if y_wmode == 0 then it writes in normal mode
      y_wmode == 1 then it writes in fast mode
*/

#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

#if defined(Y_USE_SSE2)
#include "ygensse2.h"
#include "fft5sse2.h"
#define NDEBUG1


/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* Prefetch the needed constants */
# define set_sse2_5_constants  \
prefetch_data ( &MM_F_1_5i, 0);\
prefetch_data ( &MM_F_2_5i, 0);\
prefetch_data ( &MM_FNM125, 0);


/* asign the correct values to every pad */

#define get_pads_5                     \
  jj = pad;                            \
  pd0 = d;                             \
  prefetch_data (pd0, addr((n0 << 1)));\
  pd1 = d + addr((jj << 1));           \
  jj += pad;                           \
  prefetch_data (pd1, addr((n0 << 1)));\
  pd2 = d + addr((jj << 1));           \
  jj += pad;                           \
  prefetch_data (pd2, addr((n0 << 1)));\
  pd3 = d + addr((jj << 1));           \
  jj += pad;                           \
  prefetch_data (pd3, addr((n0 << 1)));\
  pd4 = d + addr((jj << 1));           \
  prefetch_data (pd4, addr((n0 << 1)));\
  bigpad = (pad << 2) + pad;


#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2


/* get the basic twiddle from memory and computes its powers */

# define get_twiddle_factors_sse2_5                      \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);\
   sse2_square(tw2, tw1);                                \
   sse2_square(tw4, tw2);                                \
   sse2_mul(tw3, tw2, tw1);                              \
   prefetch_p_trig(px);

# define get_twiddle_factors_sse2_5_doubled     \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px); \
   sse2_square(tw2, tw1);                       \
   sse2_square(tw4, tw2);                       \
   sse2_mul(tw3, tw2, tw1);                     \
   prefetch_p_trig(px);

/* this version loads ones and trig factor */
#  define get_twiddle_factors_sse2_5_one        \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);\
   sse2_square(tw2, tw1);                       \
   sse2_square(tw4, tw2);                       \
   sse2_mul(tw3, tw2, tw1);                     \
   prefetch_p_trig(px);

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 8

# define get_twiddle_factors_sse2_5                                \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, px + 4, px + Y_STEP + 4);  \
   Y_MM_LOAD_INTER_PAIR(tw4##r, tw4##i, px + 6, px + Y_STEP + 6);

#  define get_twiddle_factors_sse2_5_doubled         \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);      \
   Y_MM_LOAD_INTER_DOUBLED(tw2##r, tw2##i, px + 2);  \
   Y_MM_LOAD_INTER_DOUBLED(tw3##r, tw3##i, px + 4);  \
   Y_MM_LOAD_INTER_DOUBLED(tw4##r, tw4##i, px + 6);

#  define get_twiddle_factors_sse2_5_one             \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);     \
   Y_MM_LOAD_INTER_ONE_HALF(tw2##r, tw2##i, px + 2); \
   Y_MM_LOAD_INTER_ONE_HALF(tw3##r, tw3##i, px + 4); \
   Y_MM_LOAD_INTER_ONE_HALF(tw4##r, tw4##i, px + 6);

#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 4


# define get_twiddle_factors_sse2_5                            \
   Y_MM_LOAD_INTER_PAIR(tw1r, tw1i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw3r, tw3i, px + 2, px + Y_STEP + 2);  \
   sse2_divmul (tw2, tw4, tw3, tw1);                           \
   prefetch_data_trig (px, 2 * Y_STEP);

#  define get_twiddle_factors_sse2_5_doubled        \
   Y_MM_LOAD_INTER_DOUBLED (tw1r, tw1i, px);        \
   Y_MM_LOAD_INTER_DOUBLED (tw3r, tw3i, px + 2);    \
   sse2_divmul (tw2, tw4, tw3, tw1);                \
   prefetch_data_trig (px, Y_STEP);

#  define get_twiddle_factors_sse2_5_one            \
   Y_MM_LOAD_INTER_ONE_HALF(tw1r, tw1i, px);        \
   Y_MM_LOAD_INTER_ONE_HALF(tw3r, tw3i, px + 2);    \
   sse2_divmul (tw2, tw4, tw3, tw1);                \
   prefetch_data_trig (px, Y_STEP);

#endif

/*
	Get data from memory, mul by twiddle trig factors,
	make a complex 9-length FFT and store the data to memory 
*/

/* read 0 write 0 */
# define radixmm0_5_twd_pp(_S, _j)                           \
  sse2_data_to_local_pp_inter (t0, _j, pd0);                 \
  sse2_load_inter_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_inter_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  sse2_fft5_store_inter_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

# define radixmm0_5_notwd_pp(_S, _j)        \
  sse2_data_to_local_pp_inter (t0, _j, pd0);\
  sse2_data_to_local_pp_inter (t1, _j, pd1);\
  sse2_data_to_local_pp_inter (t2, _j, pd2);\
  sse2_data_to_local_pp_inter (t3, _j, pd3);\
  sse2_data_to_local_pp_inter (t4, _j, pd4);\
  sse2_fft5_store_inter_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

/* read 0 write 1 */

# define radixmm1_5_twd_pp(_S, _j)                           \
  sse2_data_to_local_pp_inter (t0, _j, pd0);                 \
  sse2_load_inter_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_inter_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  sse2_fft5_store_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

# define radixmm1_5_notwd_pp(_S, _j)        \
  sse2_data_to_local_pp_inter (t0, _j, pd0);\
  sse2_data_to_local_pp_inter (t1, _j, pd1);\
  sse2_data_to_local_pp_inter (t2, _j, pd2);\
  sse2_data_to_local_pp_inter (t3, _j, pd3);\
  sse2_data_to_local_pp_inter (t4, _j, pd4);\
  sse2_fft5_store_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);


/* read 1 write 0 */

# define radixmm2_5_twd_pp(_S, _j)                     \
  sse2_data_to_local_pp (t0, _j, pd0);                 \
  sse2_load_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  sse2_fft5_store_inter_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

# define radixmm2_5_notwd_pp(_S, _j)  \
  sse2_data_to_local_pp (t0, _j, pd0);\
  sse2_data_to_local_pp (t1, _j, pd1);\
  sse2_data_to_local_pp (t2, _j, pd2);\
  sse2_data_to_local_pp (t3, _j, pd3);\
  sse2_data_to_local_pp (t4, _j, pd4);\
  sse2_fft5_store_inter_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

/* read 1 write 1 */
# define radixmm3_5_twd_pp(_S, _j)                     \
  sse2_data_to_local_pp (t0, _j, pd0);                 \
  sse2_load_mulmul_pp (t1, t2, tw1, tw2, _j, pd1, pd2);\
  sse2_load_mulmul_pp (t3, t4, tw3, tw4, _j, pd3, pd4);\
  sse2_fft5_store_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);

# define radixmm3_5_notwd_pp(_S, _j)  \
  sse2_data_to_local_pp (t0, _j, pd0);\
  sse2_data_to_local_pp (t1, _j, pd1);\
  sse2_data_to_local_pp (t2, _j, pd2);\
  sse2_data_to_local_pp (t3, _j, pd3);\
  sse2_data_to_local_pp (t4, _j, pd4);\
  sse2_fft5_store_p (_j, pd0, pd1, pd2, pd3, pd4, _S, t0, t1, t2, t3, t4);




/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r  WITHOUT TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = NULL in this no twiddle version. 
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radixmm0_5_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                         y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t nc, i, j, jj;

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm0_5_dif_notw n=%i n0=%i pad=%i \n",n, n0, pad);
#endif

  set_sse2_5_constants;
  get_pads_5;
  nc = n + n0;
  for(i = n0; i < nc; i += bigpad)
    {
      for(j = i; j < (pad + i);j += 2)
        {
          jj = addr(j << 1);
          radixmm0_5_notwd_pp (F, jj);
        }
    }
}

void radixmm1_5_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                         y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t nc, i, j, jj;

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm1_5_dif_notw n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  set_sse2_5_constants;
  get_pads_5;
  nc = n + n0;
  for(i = n0; i < nc; i += bigpad)
    {
      for(j = i; j < (pad + i);j += 2)
        {
          jj = addr(j << 1);
          radixmm1_5_notwd_pp (F, jj);
        }
    }
}

void radixmm2_5_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                         y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t nc, i, j, jj;

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm1_5_dif_notw n=%i n0=%i pad=%i \n",n, n0, pad);
#endif

  set_sse2_5_constants;
  get_pads_5;
  nc = n + n0;
  for(i = n0; i < nc; i += bigpad)
    {
      for(j = i; j < (pad + i);j += 2)
        {
          jj = addr(j << 1);
          radixmm2_5_notwd_pp (F, jj);
        }
    }
}

void radixmm3_5_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                         y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t nc, i, j, jj;

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm3_5_dif_notw n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  set_sse2_5_constants;
  get_pads_5;
  nc = n + n0;
  for(i = n0; i < nc; i += bigpad)
    {
      for(j = i; j < (pad + i);j += 2)
        {
          jj = addr(j << 1);
          radixmm3_5_notwd_pp (F, jj);
        }
    }
}


/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = Twidle Backward factors. Not scrambled.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radixmm0_5_dit (y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, tw1r, tw1i,
  tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t i, j, jj, nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm0_5_dit n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  nc = n0 + n;

  set_sse2_5_constants;
  get_pads_5;

  for(i = n0; i < nc; i += bigpad)
    {
      px = tw;
      get_twiddle_factors_sse2_5_one;
      jj = addr (i << 1);
      radixmm0_5_twd_pp (B, jj);
      j = i + 2;
      px += Y_STEP;
      for( ;j < (pad + i); j += 2, px += 2 * Y_STEP)
        {
          get_twiddle_factors_sse2_5;
          jj = addr (j << 1);
          radixmm0_5_twd_pp (B, jj);
        }
    }
}

void radixmm1_5_dit (y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, tw1r, tw1i,
  tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t i, j, jj, nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm1_5_dit n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  nc = n0 + n;

  set_sse2_5_constants;
  get_pads_5;

  for(i = n0; i < nc; i += bigpad)
    {
      px = tw;
      get_twiddle_factors_sse2_5_one;
      jj = addr (i << 1);
      radixmm1_5_twd_pp(B,jj);
      j = i + 2;
      px += Y_STEP;
      for( ;j < (pad + i); j += 2, px += 2 * Y_STEP)
        {
          get_twiddle_factors_sse2_5;
          jj = addr (j << 1);
          radixmm1_5_twd_pp (B, jj);
        }
    }
}

void radixmm2_5_dit (y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, tw1r, tw1i,
  tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t i, j, jj, nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm2_5_dit n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  nc = n0 + n;

  set_sse2_5_constants;
  get_pads_5;

  for(i = n0; i < nc; i += bigpad)
    {
      px = tw;
      get_twiddle_factors_sse2_5_one;
      jj = addr (i << 1);
      radixmm2_5_twd_pp(B,jj);
      j = i + 2;
      px += Y_STEP;
      for( ;j < (pad + i); j += 2, px += 2 * Y_STEP)
        {
          get_twiddle_factors_sse2_5;
          jj = addr (j << 1);
          radixmm2_5_twd_pp (B, jj);
        }
    }
}

void radixmm3_5_dit (y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, tw1r, tw1i,
  tw2r, tw2i, tw3r, tw3i, tw4r, tw4i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4;
  y_size_t i, j, jj, nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmm3_5_dit n=%i n0=%i pad=%i \n", n, n0, pad);
#endif

  nc = n0 + n;

  set_sse2_5_constants;
  get_pads_5;

  for(i = n0; i < nc; i += bigpad)
    {
      px = tw;
      get_twiddle_factors_sse2_5_one;
      jj = addr (i << 1);
      radixmm3_5_twd_pp(B,jj);
      j = i + 2;
      px += Y_STEP;
      for( ;j < (pad + i); j += 2, px += 2 * Y_STEP)
        {
          get_twiddle_factors_sse2_5;
          jj = addr (j << 1);
          radixmm3_5_twd_pp (B, jj);
        }
    }
}

#else /* Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radixmm_5_void( void )
{
  printf("The routine radixmm_5_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif /* Y_USE_SSE2 */
/*$Id$*/
