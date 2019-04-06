/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor
 
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

#define Y_INC_AHEAD 0

#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

/* if not defined Y_USE_SSE2 then skip this file */
#if defined(Y_USE_SSE2)
# include "ygensse2.h"

# define NDEBUG1
# define Y_SAVE_TRIGS

/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* Prefetch the needed constants */
# define set_sse2_8_constants  prefetch_data ( &MM_F_1_8r, 0);


# define PRINTMM(_op) printf(" %lf,%lf : %lf,%lf\n",_op##r.d[0],_op##r.d[1],_op##i.d[0],_op##i.d[1]);

/* asign the correct values to every pad */
#define get_pads_8                   \
  bigpad = (pad<<3);                 \
  jj = pad;                          \
  pd0 = d;                           \
  prefetch_data( pd0, addr((n0<<1)));\
  pd1 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd1, addr((n0<<1)));\
  pd2 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd2, addr((n0<<1)));\
  pd3 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd3, addr((n0<<1)));\
  pd4 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd4, addr((n0<<1)));\
  pd5 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd5, addr((n0<<1)));\
  pd6 = d + addr((jj<<1));           \
  jj += pad;                         \
  prefetch_data( pd6, addr((n0<<1)));\
  pd7 = d + addr((jj<<1));           \
  prefetch_data( pd7, addr((n0<<1)));

/* get the basic twiddle from memory and computes its powers */
# if defined(Y_MINIMUM)

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 2

#  define get_twiddle_factors_sse2_8                     \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);\
   sse2_square(tw2,tw1);                                 \
   sse2_square(tw4,tw2);                                 \
   sse2_divmul(tw3,tw5,tw4,tw1);                         \
   sse2_square(tw6,tw3);                                 \
   sse2_mul(tw7,tw4,tw3);                                \
   prefetch_p_trig_nta (px);

#  define get_twiddle_factors_sse2_8_doubled             \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);          \
   sse2_square(tw2,tw1);                                 \
   sse2_square(tw4,tw2);                                 \
   sse2_divmul(tw3,tw5,tw4,tw1);                         \
   sse2_square(tw6,tw3);                                 \
   sse2_mul(tw7,tw4,tw3);                                \
   prefetch_p_trig(px);

/* this version loads ones and trig factor */
#  define get_twiddle_factors_sse2_8_one        \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);\
   sse2_square(tw2,tw1);                        \
   sse2_square(tw4,tw2);                        \
   sse2_divmul(tw3,tw5,tw4,tw1);                \
   sse2_square(tw6,tw3);                        \
   sse2_mul(tw7,tw4,tw3);                       \
   prefetch_p_trig(px);

# elif defined(Y_MAXIMUM)

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 14

#  define get_twiddle_factors_sse2_8                               \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, px + 4, px + Y_STEP + 4);  \
   Y_MM_LOAD_INTER_PAIR(tw4##r, tw4##i, px + 6, px + Y_STEP + 6);  \
   Y_MM_LOAD_INTER_PAIR(tw5##r, tw5##i, px + 8, px + Y_STEP + 8);  \
   Y_MM_LOAD_INTER_PAIR(tw6##r, tw6##i, px + 10, px + Y_STEP + 10);\
   Y_MM_LOAD_INTER_PAIR(tw7##r, tw7##i, px + 12, px + Y_STEP + 12);

#  define get_twiddle_factors_sse2_8_doubled         \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);      \
   Y_MM_LOAD_INTER_DOUBLED(tw2##r, tw2##i, px + 2);  \
   Y_MM_LOAD_INTER_DOUBLED(tw3##r, tw3##i, px + 4);  \
   Y_MM_LOAD_INTER_DOUBLED(tw4##r, tw4##i, px + 6);  \
   Y_MM_LOAD_INTER_DOUBLED(tw5##r, tw5##i, px + 8);  \
   Y_MM_LOAD_INTER_DOUBLED(tw6##r, tw6##i, px + 10); \
   Y_MM_LOAD_INTER_DOUBLED(tw7##r, tw7##i, px + 12);

#  define get_twiddle_factors_sse2_8_one             \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);     \
   Y_MM_LOAD_INTER_ONE_HALF(tw2##r, tw2##i, px + 2); \
   Y_MM_LOAD_INTER_ONE_HALF(tw3##r, tw3##i, px + 4); \
   Y_MM_LOAD_INTER_ONE_HALF(tw4##r, tw4##i, px + 6); \
   Y_MM_LOAD_INTER_ONE_HALF(tw5##r, tw5##i, px + 8); \
   Y_MM_LOAD_INTER_ONE_HALF(tw6##r, tw6##i, px + 10);\
   Y_MM_LOAD_INTER_ONE_HALF(tw7##r, tw7##i, px + 12);

# else

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 6

#  define get_twiddle_factors_sse2_8                           \
   Y_MM_LOAD_INTER_PAIR(tw1r, tw1i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2r, tw2i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw5r, tw5i, px + 4, px + Y_STEP + 4);  \
   sse2_divmul(tw3,tw7,tw5,tw2);                               \
   sse2_divmul(tw4,tw6,tw5,tw1);                               \
   prefetch_data_trig_nta(px, 2 * Y_STEP);                     \
   prefetch_data_trig_nta(px, 3 * Y_STEP);

#  define get_twiddle_factors_sse2_8_doubled        \
   Y_MM_LOAD_INTER_DOUBLED(tw1r, tw1i, px);         \
   Y_MM_LOAD_INTER_DOUBLED(tw2r, tw2i, px + 2);     \
   Y_MM_LOAD_INTER_DOUBLED(tw5r, tw5i, px + 4);     \
   sse2_divmul(tw3,tw7,tw5,tw2);                    \
   sse2_divmul(tw4,tw6,tw5,tw1);                    \
   prefetch_data_trig (px,Y_STEP);                \
   prefetch_data_trig (px,Y_STEP + Y_CACHE_LINE);

#  define get_twiddle_factors_sse2_8_one            \
   Y_MM_LOAD_INTER_ONE_HALF(tw1r, tw1i, px);        \
   Y_MM_LOAD_INTER_ONE_HALF(tw2r, tw2i, px + 2);    \
   Y_MM_LOAD_INTER_ONE_HALF(tw5r, tw5i, px + 4);    \
   sse2_divmul(tw3,tw7,tw5,tw2);                    \
   sse2_divmul(tw4,tw6,tw5,tw1);                    \
   prefetch_data_trig(px,Y_STEP);                   \
   prefetch_data_trig(px,Y_STEP + Y_CACHE_LINE);
# endif

/* read 0 write 0 */
# define radixmm0_8_twd_pp(_S,_j)                             \
   sse2_load_inter_muladdsub_pp(t0,t1,tw4,_j,pd0, pd4);       \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw2,tw6,_j,pd2,pd6); \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_inter_mulmuladdsub_pp(t4,t5,tw1,tw5,_j,pd1,pd5); \
   sse2_load_inter_mulmuladdsub_pp(t6,t7,tw3,tw7,_j,pd3,pd7); \
   sse2_addsub(t4,t6);                                        \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_addsub_store_inter_p(_j,pd0,pd4,t0,t4);               \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd1,pd5,t1,t5);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd2,pd6,t2,t6);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd3,pd7,t3,t7);

# define radixmm0_8_notwd_pp(_S,_j)                           \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd4);               \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd2,pd6);               \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_inter_addsub_pp(t4,t5,_j,pd1,pd5);               \
   sse2_load_inter_addsub_pp(t6,t7,_j,pd3,pd7);               \
   sse2_addsub(t4,t6);                                        \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_addsub_store_inter_p(_j,pd0,pd4,t0,t4);               \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd1,pd5,t1,t5);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd2,pd6,t2,t6);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd3,pd7,t3,t7);

/* read 0 write 1 */
# define radixmm1_8_twd_pp(_S,_j)                             \
   sse2_load_inter_muladdsub_pp(t0,t1,tw4,_j,pd0, pd4);       \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw2,tw6,_j,pd2,pd6); \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_inter_mulmuladdsub_pp(t4,t5,tw1,tw5,_j,pd1,pd5); \
   sse2_load_inter_mulmuladdsub_pp(t6,t7,tw3,tw7,_j,pd3,pd7); \
   sse2_addsub(t4,t6);                                        \
   sse2_addsub_store_p(_j,pd0,pd4,t0,t4);                     \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);      \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

# define radixmm1_8_notwd_pp(_S,_j)                           \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd4);               \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd2,pd6);               \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_inter_addsub_pp(t4,t5,_j,pd1,pd5);               \
   sse2_load_inter_addsub_pp(t6,t7,_j,pd3,pd7);               \
   sse2_addsub(t4,t6);                                        \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_addsub_store_p(_j,pd0,pd4,t0,t4);                     \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);      \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

/* read 1 write 0 */
# define radixmm2_8_twd_pp(_S,_j)                             \
   sse2_load_muladdsub_pp(t0,t1,tw4,_j,pd0, pd4);             \
   sse2_load_mulmuladdsub_pp(t2,t3,tw2,tw6,_j,pd2,pd6);       \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_mulmuladdsub_pp(t4,t5,tw1,tw5,_j,pd1,pd5);       \
   sse2_load_mulmuladdsub_pp(t6,t7,tw3,tw7,_j,pd3,pd7);       \
   sse2_addsub(t4,t6);                                        \
   sse2_addsub_store_inter_p(_j,pd0,pd4,t0,t4);               \
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd2,pd6,t2,t6);\
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd1,pd5,t1,t5);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd3,pd7,t3,t7);

# define radixmm2_8_notwd_pp(_S,_j)                           \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd4);                     \
   sse2_load_addsub_pp(t2,t3,_j,pd2,pd6);                     \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_addsub_pp(t4,t5,_j,pd1,pd5);                     \
   sse2_load_addsub_pp(t6,t7,_j,pd3,pd7);                     \
   sse2_addsub(t4,t6);                                        \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_addsub_store_inter_p(_j,pd0,pd4,t0,t4);               \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd1,pd5,t1,t5);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd2,pd6,t2,t6);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd3,pd7,t3,t7);

/* read 1 write 1 */
# define radixmm3_8_twd_pp(_S,_j)                             \
   sse2_load_muladdsub_pp(t0,t1,tw4,_j,pd0, pd4);             \
   sse2_load_mulmuladdsub_pp(t2,t3,tw2,tw6,_j,pd2,pd6);       \
   sse2_addsub(t0,t2);                                        \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                         \
   sse2_load_mulmuladdsub_pp(t4,t5,tw1,tw5,_j,pd1,pd5);       \
   sse2_load_mulmuladdsub_pp(t6,t7,tw3,tw7,_j,pd3,pd7);       \
   sse2_addsub(t4,t6);                                        \
   sse2_addsub_store_p(_j,pd0,pd4,t0,t4);                     \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);      \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                         \
                                                              \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

# define radixmm3_8_notwd_pp(_S,_j)                     \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd4);               \
   sse2_load_addsub_pp(t2,t3,_j,pd2,pd6);               \
   sse2_addsub(t0,t2);                                  \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                   \
   sse2_load_addsub_pp(t4,t5,_j,pd1,pd5);               \
   sse2_load_addsub_pp(t6,t7,_j,pd3,pd7);               \
   sse2_addsub(t4,t6);                                  \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                   \
                                                        \
   sse2_addsub_store_p(_j,pd0,pd4,t0,t4);               \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);\
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);\
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);



/* SSE2 versioned */
/* Mode 0: read =0 (standard)  write= 0 (standard) */
void radixmm0_8_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_8_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm0_8_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm0_8_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm0_8_twd_pp(F,jj);
            }
        }
    }
}

void radixmm0_8_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_8_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;
  nc=n+n0;
  for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
    {
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radixmm0_8_notwd_pp(F,jj);
        }
    }
}

void radixmm0_8_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_8_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm0_8_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm0_8_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm0_8_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm0_8_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm0_8_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm0_8_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm0_8_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm0_8_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

/* Mode 1: read =0 (standard)  write= 1 (fast-sse2) */
void radixmm1_8_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_8_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm1_8_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm1_8_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm1_8_twd_pp(F,jj);
            }
        }
    }
}

void radixmm1_8_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i,
  t5r, t5i, t6r, t6i, t7r, t7i;
  y_size_t bigpad;
  y_ptr pd0, pd1, pd2, pd3, pd4, pd5, pd6, pd7;
  y_size_t i, j, jj, nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_8_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;
  nc = n + n0;
  for(i = n0; i < nc; i += bigpad)
    {
      for(j = i; j < (pad + i); j += 2)
        {
          jj = addr(j << 1);
          radixmm1_8_notwd_pp(F,jj);
        }
    }
}

void radixmm1_8_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_8_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm1_8_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm1_8_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm1_8_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm1_8_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm1_8_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm1_8_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm1_8_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm1_8_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

/* Mode 2: read =1 (fast-sse2)  write= 0 (standard) */
void radixmm2_8_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_8_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm2_8_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm2_8_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm2_8_twd_pp(F,jj);
            }
        }
    }
}


void radixmm2_8_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_8_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm2_8_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm2_8_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm2_8_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm2_8_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm2_8_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm2_8_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm2_8_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm2_8_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

/* Mode 3: read =1 (fast-sse2)  write= 1 (fast-sse2) */
void radixmm3_8_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm3_8_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm3_8_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm3_8_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_sse2_8_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              /*BEGIN_TIME;*/
              radixmm3_8_twd_pp(F,jj);
              /*END_TIME;*/
            }
        }
    }
}


void radixmm3_8_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm3_8_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_8_constants;
  get_pads_8;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm3_8_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm3_8_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm3_8_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm3_8_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          /*BEGIN_TIME;*/
          radixmm3_8_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm3_8_twd_pp(B,jj1);
          /*END_TIME;*/
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_8_one;
      jj=addr(i<<1);
      radixmm3_8_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_8;
          jj=addr(j<<1);
          radixmm3_8_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
}

#else /* defined Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radixmm_8_void( void )
{
  printf("The routine radixmm_8_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif
/* defined Y_USE_SSE2 */
