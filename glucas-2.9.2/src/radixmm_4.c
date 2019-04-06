/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2002-2006  Guillermo Ballester Valor
 
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
/* if not defined Y_USE_SSE2 then skip this file */
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

#if defined(Y_USE_SSE2)
# include "ygensse2.h"

# define NDEBUG1
# define Y_SAVE_TRIGS

/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

#define get_pads_4         \
  bigpad = (pad<<2);       \
  jj = pad;                \
  pd0 = d;                 \
  prefetch_data( pd0, 0);  \
  pd1 = d + addr((jj<<1)); \
  jj += pad;               \
  prefetch_data( pd1, 0);  \
  pd2 = d + addr((jj<<1)); \
  jj += pad;               \
  prefetch_data( pd2, 0);  \
  pd3 = d + addr((jj<<1)); \
  prefetch_data( pd3, 0);

#if defined(Y_MAXIMUM)

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

#  define get_twiddle_factors_sse2_4                               \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, px + 4, px + Y_STEP + 4);


#  define get_twiddle_factors_sse2_4_doubled         \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);      \
   Y_MM_LOAD_INTER_DOUBLED(tw2##r, tw2##i, px + 2);  \
   Y_MM_LOAD_INTER_DOUBLED(tw3##r, tw3##i, px + 4);

#  define get_twiddle_factors_sse2_4_one             \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);     \
   Y_MM_LOAD_INTER_ONE_HALF(tw2##r, tw2##i, px + 2); \
   Y_MM_LOAD_INTER_ONE_HALF(tw3##r, tw3##i, px + 4);

#else

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

#  define get_twiddle_factors_sse2_4                     \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);\
   sse2_square(tw2,tw1);                                 \
   sse2_mul(tw3,tw2,tw1);                                \
   prefetch_p_trig_nta(px);

#  define get_twiddle_factors_sse2_4_doubled             \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);          \
   sse2_square(tw2,tw1);                                 \
   sse2_mul(tw3,tw2,tw1);                                \
   prefetch_p_trig(px);

#  define get_twiddle_factors_sse2_4_one        \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);\
   sse2_square(tw2,tw1);                        \
   sse2_mul(tw3,tw2,tw1);                       \
   prefetch_p_trig(px);

#endif

/* read 0 write 0 */
# define radixmm0_4_twd_pp(_S,_j)                            \
   sse2_load_inter_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);       \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3);\
   sse2_addsub_store_inter_p(_j,pd0,pd2,t0,t2);              \
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd1,pd3,t1,t3);

# define radixmm0_4_notwd_pp(_S,_j)                          \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd2);              \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd1,pd3);              \
   sse2_addsub_store_inter_p(_j,pd0,pd2,t0,t2);              \
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd1,pd3,t1,t3);


/* read 0 write 1 */
# define radixmm1_4_twd_pp(_S,_j)                            \
   sse2_load_inter_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);       \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3);\
   sse2_addsub_store_p(_j,pd0,pd2,t0,t2);                    \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);


# define radixmm1_4_notwd_pp(_S,_j)                          \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd2);              \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd1,pd3);              \
   sse2_addsub_store_p(_j,pd0,pd2,t0,t2);                    \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);


/* read 1 write 0 */
# define radixmm2_4_twd_pp(_S,_j)                            \
   sse2_load_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);             \
   sse2_load_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3);      \
   sse2_addsub_store_inter_p(_j,pd0,pd2,t0,t2);              \
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd1,pd3,t1,t3);

# define radixmm2_4_notwd_pp(_S,_j)                          \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd2);                    \
   sse2_load_addsub_pp(t2,t3,_j,pd1,pd3);                    \
   sse2_addsub_store_inter_p(_j,pd0,pd2,t0,t2);              \
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd1,pd3,t1,t3);

/* read 1 write 1 */
# define radixmm3_4_twd_pp(_S,_j)                       \
   sse2_load_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);        \
   sse2_load_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3); \
   sse2_addsub_store_p(_j,pd0,pd2,t0,t2);               \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);

# define radixmm3_4_notwd_pp(_S,_j)                    \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd2);              \
   sse2_load_addsub_pp(t2,t3,_j,pd1,pd3);              \
   sse2_addsub_store_p(_j,pd0,pd2,t0,t2);              \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);


/* Mode 0: read =0   write= 0 */
void radixmm0_4_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_4_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm0_4_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm0_4_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm0_4_twd_pp(F,jj);
            }
        }
    }
}

void radixmm0_4_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_4_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;
  nc=n+n0;
  for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
    {
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radixmm0_4_notwd_pp(F,jj);
        }
    }
}


void radixmm0_4_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_4_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm0_4_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm0_4_twd_pp(B,jj);
        }
    }
# else
  /*
  Y_SAVE_TRIGS defined

  This optimization tries to use the trig factor load/or computed in a cycle
  in the following cycle (bigpad ahead).

  It is no aplicable when defined Y_ITANIUM, Y_PRELOAD_TRIGS, Y_PRELOAD_DATA
  or prefetch expensive
  */
  i=n0;
  i1=n0+bigpad;
  while(i1<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm0_4_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm0_4_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm0_4_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm0_4_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm0_4_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm0_4_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

/* Mode 0: read =0   write= 0 */
void radixmm1_4_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_4_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm1_4_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm1_4_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm1_4_twd_pp(F,jj);
            }
        }
    }
}

void radixmm1_4_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_4_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;
  nc=n+n0;
  for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
    {
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radixmm1_4_notwd_pp(F,jj);
        }
    }
}


void radixmm1_4_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_4_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm1_4_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm1_4_twd_pp(B,jj);
        }
    }
# else
  /*
  Y_SAVE_TRIGS defined

  This optimization tries to use the trig factor load/or computed in a cycle
  in the following cycle (bigpad ahead).

  It is no aplicable when defined Y_ITANIUM, Y_PRELOAD_TRIGS, Y_PRELOAD_DATA
  or prefetch expensive
  */
  i=n0;
  i1=n0+bigpad;
  while(i1<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm1_4_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm1_4_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm1_4_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm1_4_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm1_4_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm1_4_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}


/* Mode 0: read =0   write= 0 */
void radixmm2_4_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_4_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm2_4_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm2_4_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm2_4_twd_pp(F,jj);
            }
        }
    }
}

void radixmm2_4_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_4_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm2_4_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm2_4_twd_pp(B,jj);
        }
    }
# else
  /*
  Y_SAVE_TRIGS defined

  This optimization tries to use the trig factor load/or computed in a cycle
  in the following cycle (bigpad ahead).

  It is no aplicable when defined Y_ITANIUM, Y_PRELOAD_TRIGS, Y_PRELOAD_DATA
  or prefetch expensive
  */
  i=n0;
  i1=n0+bigpad;
  while(i1<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm2_4_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm2_4_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm2_4_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm2_4_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm2_4_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm2_4_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}


/* Mode 3: read =1   write= 1 */
void radixmm3_4_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_4_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm3_4_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm3_4_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_4_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm3_4_twd_pp(F,jj);
            }
        }
    }
}


void radixmm3_4_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm3_4_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  get_pads_4;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm3_4_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm3_4_twd_pp(B,jj);
        }
    }
# else
  /*
  Y_SAVE_TRIGS defined

  This optimization tries to use the trig factor load/or computed in a cycle
  in the following cycle (bigpad ahead).

  It is no aplicable when defined Y_ITANIUM, Y_PRELOAD_TRIGS, Y_PRELOAD_DATA
  or prefetch expensive
  */
  i=n0;
  i1=n0+bigpad;
  while(i1<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm3_4_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm3_4_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm3_4_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm3_4_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_4_one;
      jj=addr(i<<1);
      radixmm3_4_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_4;
          jj=addr(j<<1);
          radixmm3_4_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

#else /* defined Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radixmm_4_void( void )
{
  printf("The routine radixmm_4_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif
/* defined Y_USE_SSE2 */
