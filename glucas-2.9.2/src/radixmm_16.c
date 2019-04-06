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
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

/* if not defined Y_USE_SSE2 then skip this file */
#if defined(Y_USE_SSE2) && (Y_AVAL > 3)
# include "ygensse2.h"

# define NDEBUG1
# define Y_SAVE_TRIGS

/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* Set constant */

# define set_sse2_16_constants_F    \
  prefetch_data ( &MM_F_1_8r , 0 ); \
  prefetch_data ( &MM_F_1_16r, 0 ); \
  prefetch_data ( &MM_F_3_16r, 0 ); \
  prefetch_data ( &MM_F_5_16r, 0 ); \
  prefetch_data ( &MM_F_7_16r, 0 );

# define set_sse2_16_constants_B    \
  prefetch_data ( &MM_F_1_8r, 0 );  \
  prefetch_data ( &MM_B_1_16r, 0 ); \
  prefetch_data ( &MM_B_3_16r, 0 ); \
  prefetch_data ( &MM_B_5_16r, 0 ); \
  prefetch_data ( &MM_B_7_16r, 0 );

/* asign the correct values to every pad */
#define get_pads_16                   \
  bigpad = (pad<<4);                  \
  jj = pad;                           \
  pd0 = d;                            \
  prefetch_data( pd0, addr((n0<<1))); \
  pd1 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd1, addr((n0<<1))); \
  pd2 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd2, addr((n0<<1))); \
  pd3 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd3, addr((n0<<1))); \
  pd4 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd4, addr((n0<<1))); \
  pd5 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd5, addr((n0<<1))); \
  pd6 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd6, addr((n0<<1))); \
  pd7 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd7, addr((n0<<1))); \
  pd8 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd8, addr((n0<<1))); \
  pd9 = d + addr((jj<<1));            \
  jj += pad;                          \
  prefetch_data( pd9, addr((n0<<1))); \
  pd10 = d + addr((jj<<1));           \
  jj += pad;                          \
  prefetch_data( pd10, addr((n0<<1)));\
  pd11 = d + addr((jj<<1));           \
  jj += pad;                          \
  prefetch_data( pd11, addr((n0<<1)));\
  pd12 = d + addr((jj<<1));           \
  jj += pad;                          \
  prefetch_data( pd12, addr((n0<<1)));\
  pd13 = d + addr((jj<<1));           \
  jj += pad;                          \
  prefetch_data( pd13, addr((n0<<1)));\
  pd14 = d + addr((jj<<1));           \
  jj += pad;                          \
  prefetch_data( pd14, addr((n0<<1)));\
  pd15 = d + addr((jj<<1));           \
  prefetch_data( pd15, addr((n0<<1)));

/*
# define prefetch_data_16( _i) \
{                              \
  size_t off;                  \
  off = addr((_i << 1));       \
  prefetch_data( pd0, off);    \
  prefetch_data( pd1, off);    \
  prefetch_data( pd2, off);    \
  prefetch_data( pd3, off);    \
  prefetch_data( pd4, off);    \
  prefetch_data( pd5, off);    \
  prefetch_data( pd6, off);    \
  prefetch_data( pd7, off);    \
  prefetch_data( pd8, off);    \
  prefetch_data( pd9, off);    \
  prefetch_data( pd10, off);   \
  prefetch_data( pd11, off);   \
  prefetch_data( pd12, off);   \
  prefetch_data( pd13, off);   \
  prefetch_data( pd14, off);   \
  prefetch_data( pd15, off);   \
}
*/
# define prefetch_data_16( _i) /* */

/* the basic twiddle from memory and computes its powers */
# if defined(Y_MINIMUM)

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 2

#  define get_twiddle_factors_sse2_16                    \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);\
   sse2_square(tw2,tw1);                                 \
   sse2_square(tw4,tw2);                                 \
   sse2_square(tw8,tw4);                                 \
   sse2_divmul(tw3,tw5,tw4,tw1);                         \
   sse2_divmul(tw7,tw9,tw8,tw1);                         \
   sse2_divmul(tw6,tw10,tw8,tw2);                        \
   sse2_square(tw12,tw6);                                \
   sse2_divmul(tw11,tw13,tw12,tw1);                      \
   sse2_square(tw14,tw7);                                \
   sse2_mul(tw15,tw7,tw8);                               \
   prefetch_p_trig_nta(px);

#  define get_twiddle_factors_sse2_16_doubled            \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);          \
   sse2_square(tw2,tw1);                                 \
   sse2_square(tw4,tw2);                                 \
   sse2_square(tw8,tw4);                                 \
   sse2_divmul(tw3,tw5,tw4,tw1);                         \
   sse2_divmul(tw7,tw9,tw8,tw1);                         \
   sse2_divmul(tw6,tw10,tw8,tw2);                        \
   sse2_square(tw12,tw6);                                \
   sse2_divmul(tw11,tw13,tw12,tw1);                      \
   sse2_square(tw14,tw7);                                \
   sse2_mul(tw15,tw7,tw8);                               \
   prefetch_p_trig(px);

/* this version loads ones and trig factor */
#  define get_twiddle_factors_sse2_16_one        \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);\
   sse2_square(tw2,tw1);                        \
   sse2_square(tw4,tw2);                        \
   sse2_square(tw8,tw4);                        \
   sse2_divmul(tw3,tw5,tw4,tw1);                \
   sse2_divmul(tw7,tw9,tw8,tw1);                \
   sse2_divmul(tw6,tw10,tw8,tw2);               \
   sse2_square(tw12,tw6);                       \
   sse2_divmul(tw11,tw13,tw12,tw1);             \
   sse2_square(tw14,tw7);                       \
   sse2_mul(tw15,tw7,tw8);                      \
   prefetch_p_trig(px);

# elif defined(Y_MAXIMUM)

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 30

#  define get_twiddle_factors_sse2_16                               \
   Y_MM_LOAD_INTER_PAIR(tw1##r, tw1##i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2##r, tw2##i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw3##r, tw3##i, px + 4, px + Y_STEP + 4);  \
   Y_MM_LOAD_INTER_PAIR(tw4##r, tw4##i, px + 6, px + Y_STEP + 6);  \
   Y_MM_LOAD_INTER_PAIR(tw5##r, tw5##i, px + 8, px + Y_STEP + 8);  \
   Y_MM_LOAD_INTER_PAIR(tw6##r, tw6##i, px + 10, px + Y_STEP + 10);\
   Y_MM_LOAD_INTER_PAIR(tw7##r, tw7##i, px + 12, px + Y_STEP + 12);\
   Y_MM_LOAD_INTER_PAIR(tw8##r, tw8##i, px + 14, px + Y_STEP + 14);\
   Y_MM_LOAD_INTER_PAIR(tw9##r, tw9##i, px + 16, px + Y_STEP + 16);\
   Y_MM_LOAD_INTER_PAIR(tw10##r, tw10##i, px + 18, px + Y_STEP + 18);\
   Y_MM_LOAD_INTER_PAIR(tw11##r, tw11##i, px + 20, px + Y_STEP + 20);\
   Y_MM_LOAD_INTER_PAIR(tw12##r, tw12##i, px + 22, px + Y_STEP + 22);\
   Y_MM_LOAD_INTER_PAIR(tw13##r, tw13##i, px + 24, px + Y_STEP + 24);\
   Y_MM_LOAD_INTER_PAIR(tw14##r, tw14##i, px + 26, px + Y_STEP + 26);\
   Y_MM_LOAD_INTER_PAIR(tw15##r, tw15##i, px + 28, px + Y_STEP + 28);

#  define get_twiddle_factors_sse2_16_doubled          \
   Y_MM_LOAD_INTER_DOUBLED(tw1##r, tw1##i, px);       \
   Y_MM_LOAD_INTER_DOUBLED(tw2##r, tw2##i, px + 2);   \
   Y_MM_LOAD_INTER_DOUBLED(tw3##r, tw3##i, px + 4);   \
   Y_MM_LOAD_INTER_DOUBLED(tw4##r, tw4##i, px + 6);   \
   Y_MM_LOAD_INTER_DOUBLED(tw5##r, tw5##i, px + 8);   \
   Y_MM_LOAD_INTER_DOUBLED(tw6##r, tw6##i, px + 10);  \
   Y_MM_LOAD_INTER_DOUBLED(tw7##r, tw7##i, px + 12);  \
   Y_MM_LOAD_INTER_DOUBLED(tw8##r, tw8##i, px + 14);  \
   Y_MM_LOAD_INTER_DOUBLED(tw9##r, tw9##i, px + 16);  \
   Y_MM_LOAD_INTER_DOUBLED(tw10##r, tw10##i, px + 18);\
   Y_MM_LOAD_INTER_DOUBLED(tw11##r, tw11##i, px + 20);\
   Y_MM_LOAD_INTER_DOUBLED(tw12##r, tw12##i, px + 22);\
   Y_MM_LOAD_INTER_DOUBLED(tw13##r, tw13##i, px + 24);\
   Y_MM_LOAD_INTER_DOUBLED(tw14##r, tw14##i, px + 26);\
   Y_MM_LOAD_INTER_DOUBLED(tw15##r, tw15##i, px + 28);

#  define get_twiddle_factors_sse2_16_one               \
   Y_MM_LOAD_INTER_ONE_HALF(tw1##r, tw1##i, px);       \
   Y_MM_LOAD_INTER_ONE_HALF(tw2##r, tw2##i, px + 2);   \
   Y_MM_LOAD_INTER_ONE_HALF(tw3##r, tw3##i, px + 4);   \
   Y_MM_LOAD_INTER_ONE_HALF(tw4##r, tw4##i, px + 6);   \
   Y_MM_LOAD_INTER_ONE_HALF(tw5##r, tw5##i, px + 8);   \
   Y_MM_LOAD_INTER_ONE_HALF(tw6##r, tw6##i, px + 10);  \
   Y_MM_LOAD_INTER_ONE_HALF(tw7##r, tw7##i, px + 12);  \
   Y_MM_LOAD_INTER_ONE_HALF(tw8##r, tw8##i, px + 14);  \
   Y_MM_LOAD_INTER_ONE_HALF(tw9##r, tw9##i, px + 16);  \
   Y_MM_LOAD_INTER_ONE_HALF(tw10##r, tw10##i, px + 18);\
   Y_MM_LOAD_INTER_ONE_HALF(tw11##r, tw11##i, px + 20);\
   Y_MM_LOAD_INTER_ONE_HALF(tw12##r, tw12##i, px + 22);\
   Y_MM_LOAD_INTER_ONE_HALF(tw13##r, tw13##i, px + 24);\
   Y_MM_LOAD_INTER_ONE_HALF(tw14##r, tw14##i, px + 26);\
   Y_MM_LOAD_INTER_ONE_HALF(tw15##r, tw15##i, px + 28);

# else

#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 10


#  define get_twiddle_factors_sse2_16                          \
   Y_MM_LOAD_INTER_PAIR(tw1r, tw1i, px, px + Y_STEP);          \
   Y_MM_LOAD_INTER_PAIR(tw2r, tw2i, px + 2, px + Y_STEP + 2);  \
   Y_MM_LOAD_INTER_PAIR(tw4r, tw4i, px + 4, px + Y_STEP + 4);  \
   sse2_divmul(tw3,tw5,tw4,tw1);                               \
   Y_MM_LOAD_INTER_PAIR(tw8r, tw8i, px + 6, px + Y_STEP + 6);  \
   Y_MM_LOAD_INTER_PAIR(tw13r, tw13i, px + 8, px + Y_STEP + 8);\
   sse2_divmul(tw6,tw10,tw8,tw2);                              \
   sse2_divmul(tw7,tw9,tw8,tw1);                               \
   sse2_divmul(tw11,tw15,tw13,tw2);                            \
   sse2_divmul(tw12,tw14,tw13,tw1);                            \
   prefetch_data_trig_nta(px,2 * Y_STEP);                      \
   prefetch_data_trig_nta(px,3 * Y_STEP);

#  define get_twiddle_factors_sse2_16_doubled       \
   Y_MM_LOAD_INTER_DOUBLED(tw1r, tw1i, px);         \
   Y_MM_LOAD_INTER_DOUBLED(tw2r, tw2i, px + 2);     \
   Y_MM_LOAD_INTER_DOUBLED(tw4r, tw4i, px + 4);     \
   sse2_divmul(tw3,tw5,tw4,tw1);                    \
   Y_MM_LOAD_INTER_DOUBLED(tw8r, tw8i, px + 6);     \
   Y_MM_LOAD_INTER_DOUBLED(tw13r, tw13i, px + 8);   \
   sse2_divmul(tw6,tw10,tw8,tw2);                   \
   sse2_divmul(tw7,tw9,tw8,tw1);                    \
   sse2_divmul(tw11,tw15,tw13,tw2);                 \
   sse2_divmul(tw12,tw14,tw13,tw1);                 \
   prefetch_data_trig(px,Y_STEP);                   \
   prefetch_data_trig(px,Y_STEP + Y_CACHE_LINE);

#  define get_twiddle_factors_sse2_16_one           \
   Y_MM_LOAD_INTER_ONE_HALF(tw1r, tw1i, px);        \
   Y_MM_LOAD_INTER_ONE_HALF(tw2r, tw2i, px + 2);    \
   Y_MM_LOAD_INTER_ONE_HALF(tw4r, tw4i, px + 4);    \
   sse2_divmul(tw3,tw5,tw4,tw1);                    \
   Y_MM_LOAD_INTER_ONE_HALF(tw8r, tw8i, px + 6);    \
   Y_MM_LOAD_INTER_ONE_HALF(tw13r, tw13i, px + 8);  \
   sse2_divmul(tw6,tw10,tw8,tw2);                   \
   sse2_divmul(tw7,tw9,tw8,tw1);                    \
   sse2_divmul(tw11,tw15,tw13,tw2);                 \
   sse2_divmul(tw12,tw14,tw13,tw1);                 \
   prefetch_data_trig(px,Y_STEP);                   \
   prefetch_data_trig(px,Y_STEP + Y_CACHE_LINE);

# endif

/* read 0 write 0 */
# define radixmm0_16_twd_pp(_S,_j)                                 \
   sse2_load_inter_muladdsub_pp(t0,t1,tw8,_j,pd0,pd8);             \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw4,tw12,_j,pd4,pd12);    \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_inter_mulmuladdsub_pp(t4,t5,tw2,tw10,_j,pd2,pd10);    \
   sse2_load_inter_mulmuladdsub_pp(t6,t7,tw6,tw14,_j,pd6,pd14);    \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_inter_mulmuladdsub_pp(t8,t9,tw1,tw9,_j,pd1,pd9);      \
   sse2_load_inter_mulmuladdsub_pp(t10,t11,tw5,tw13,_j,pd5,pd13);  \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_inter_mulmuladdsub_pp(t12,t13,tw3,tw11,_j,pd3,pd11);  \
   sse2_load_inter_mulmuladdsub_pp(t14,t15,tw7,tw15,_j,pd7,pd15);  \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_inter_p(_j,pd0,pd8,t0,t8);                    \
   sse2_muladdsub_store_inter_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);  \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd2,pd10,t2,t10);   \
   sse2_muladdsub_store_inter_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd4,pd12,t4,t12);   \
   sse2_muladdsub_store_inter_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd6,pd14,t6,t14);   \
   sse2_muladdsub_store_inter_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

# define radixmm0_16_notwd_pp(_S,_j)                               \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd8);                    \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd4,pd12);                   \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_inter_addsub_pp(t4,t5,_j,pd2,pd10);                   \
   sse2_load_inter_addsub_pp(t6,t7,_j,pd6,pd14);                   \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_inter_addsub_pp(t8,t9,_j,pd1,pd9);                    \
   sse2_load_inter_addsub_pp(t10,t11,_j,pd5,pd13);                 \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_inter_addsub_pp(t12,t13,_j,pd3,pd11);                 \
   sse2_load_inter_addsub_pp(t14,t15,_j,pd7,pd15);                 \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_inter_p(_j,pd0,pd8,t0,t8);                    \
   sse2_muladdsub_store_inter_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);  \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd2,pd10,t2,t10);   \
   sse2_muladdsub_store_inter_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd4,pd12,t4,t12);   \
   sse2_muladdsub_store_inter_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd6,pd14,t6,t14);   \
   sse2_muladdsub_store_inter_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

/* read 0 write 1 */
# define radixmm1_16_twd_pp(_S,_j)                                 \
   sse2_load_inter_muladdsub_pp(t0,t1,tw8,_j,pd0,pd8);             \
   sse2_load_inter_mulmuladdsub_pp(t2,t3,tw4,tw12,_j,pd4,pd12);    \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_inter_mulmuladdsub_pp(t4,t5,tw2,tw10,_j,pd2,pd10);    \
   sse2_load_inter_mulmuladdsub_pp(t6,t7,tw6,tw14,_j,pd6,pd14);    \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_inter_mulmuladdsub_pp(t8,t9,tw1,tw9,_j,pd1,pd9);      \
   sse2_load_inter_mulmuladdsub_pp(t10,t11,tw5,tw13,_j,pd5,pd13);  \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_inter_mulmuladdsub_pp(t12,t13,tw3,tw11,_j,pd3,pd11);  \
   sse2_load_inter_mulmuladdsub_pp(t14,t15,tw7,tw15,_j,pd7,pd15);  \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_p(_j,pd0,pd8,t0,t8);                          \
   sse2_muladdsub_store_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);        \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);         \
   sse2_muladdsub_store_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);      \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);         \
   sse2_muladdsub_store_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);         \
   sse2_muladdsub_store_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

# define radixmm1_16_notwd_pp(_S,_j)                               \
   sse2_load_inter_addsub_pp(t0,t1,_j,pd0,pd8);                    \
   sse2_load_inter_addsub_pp(t2,t3,_j,pd4,pd12);                   \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_inter_addsub_pp(t4,t5,_j,pd2,pd10);                   \
   sse2_load_inter_addsub_pp(t6,t7,_j,pd6,pd14);                   \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_inter_addsub_pp(t8,t9,_j,pd1,pd9);                    \
   sse2_load_inter_addsub_pp(t10,t11,_j,pd5,pd13);                 \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_inter_addsub_pp(t12,t13,_j,pd3,pd11);                 \
   sse2_load_inter_addsub_pp(t14,t15,_j,pd7,pd15);                 \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_p(_j,pd0,pd8,t0,t8);                          \
   sse2_muladdsub_store_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);        \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);         \
   sse2_muladdsub_store_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);      \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);         \
   sse2_muladdsub_store_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);         \
   sse2_muladdsub_store_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

/* read 1 write 0 */
# define radixmm2_16_twd_pp(_S,_j)                                 \
   sse2_load_muladdsub_pp(t0,t1,tw8,_j,pd0,pd8);                   \
   sse2_load_mulmuladdsub_pp(t2,t3,tw4,tw12,_j,pd4,pd12);          \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_mulmuladdsub_pp(t4,t5,tw2,tw10,_j,pd2,pd10);          \
   sse2_load_mulmuladdsub_pp(t6,t7,tw6,tw14,_j,pd6,pd14);          \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_mulmuladdsub_pp(t8,t9,tw1,tw9,_j,pd1,pd9);            \
   sse2_load_mulmuladdsub_pp(t10,t11,tw5,tw13,_j,pd5,pd13);        \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_mulmuladdsub_pp(t12,t13,tw3,tw11,_j,pd3,pd11);        \
   sse2_load_mulmuladdsub_pp(t14,t15,tw7,tw15,_j,pd7,pd15);        \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_inter_p(_j,pd0,pd8,t0,t8);                    \
   sse2_muladdsub_store_inter_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);  \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd2,pd10,t2,t10);   \
   sse2_muladdsub_store_inter_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd4,pd12,t4,t12);   \
   sse2_muladdsub_store_inter_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd6,pd14,t6,t14);   \
   sse2_muladdsub_store_inter_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

# define radixmm2_16_notwd_pp(_S,_j)                               \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd8);                          \
   sse2_load_addsub_pp(t2,t3,_j,pd4,pd12);                         \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_addsub_pp(t4,t5,_j,pd2,pd10);                         \
   sse2_load_addsub_pp(t6,t7,_j,pd6,pd14);                         \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_addsub_pp(t8,t9,_j,pd1,pd9);                          \
   sse2_load_addsub_pp(t10,t11,_j,pd5,pd13);                       \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_addsub_pp(t12,t13,_j,pd3,pd11);                       \
   sse2_load_addsub_pp(t14,t15,_j,pd7,pd15);                       \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_inter_p(_j,pd0,pd8,t0,t8);                    \
   sse2_muladdsub_store_inter_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);  \
   sse2_mul_1_8_##_S##_addsub_store_inter_p(_j,pd2,pd10,t2,t10);   \
   sse2_muladdsub_store_inter_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);\
   sse2_mul_1_4_##_S##_addsub_store_inter_p(_j,pd4,pd12,t4,t12);   \
   sse2_muladdsub_store_inter_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);\
   sse2_mul_3_8_##_S##_addsub_store_inter_p(_j,pd6,pd14,t6,t14);   \
   sse2_muladdsub_store_inter_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);


/* read 1 write 1 */
# define radixmm3_16_twd_pp(_S,_j)                                 \
   sse2_load_muladdsub_pp(t0,t1,tw8,_j,pd0,pd8);                   \
   sse2_load_mulmuladdsub_pp(t2,t3,tw4,tw12,_j,pd4,pd12);          \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_mulmuladdsub_pp(t4,t5,tw2,tw10,_j,pd2,pd10);          \
   sse2_load_mulmuladdsub_pp(t6,t7,tw6,tw14,_j,pd6,pd14);          \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_mulmuladdsub_pp(t8,t9,tw1,tw9,_j,pd1,pd9);            \
   sse2_load_mulmuladdsub_pp(t10,t11,tw5,tw13,_j,pd5,pd13);        \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_mulmuladdsub_pp(t12,t13,tw3,tw11,_j,pd3,pd11);        \
   sse2_load_mulmuladdsub_pp(t14,t15,tw7,tw15,_j,pd7,pd15);        \
   sse2_addsub(t12,t14);                                           \
   sse2_addsub(t8,t12);                                            \
   sse2_addsub_store_p(_j,pd0,pd8,t0,t8);                          \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);         \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_muladdsub_store_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);        \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);         \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_muladdsub_store_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);      \
   sse2_muladdsub_store_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);         \
   sse2_muladdsub_store_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);

# define radixmm3_16_notwd_pp(_S,_j)                               \
   sse2_load_addsub_pp(t0,t1,_j,pd0,pd8);                          \
   sse2_load_addsub_pp(t2,t3,_j,pd4,pd12);                         \
   sse2_addsub(t0,t2);                                             \
   sse2_mul_1_4_##_S##_addsub(t1,t3);                              \
   sse2_load_addsub_pp(t4,t5,_j,pd2,pd10);                         \
   sse2_load_addsub_pp(t6,t7,_j,pd6,pd14);                         \
   sse2_addsub(t4,t6);                                             \
   sse2_mul_1_4_##_S##_addsub(t5,t7);                              \
   sse2_addsub(t0,t4);                                             \
   sse2_mul_1_8_##_S##_addsub(t1,t5);                              \
   sse2_mul_1_4_##_S##_addsub(t2,t6);                              \
   sse2_mul_3_8_##_S##_addsub(t3,t7);                              \
   sse2_load_addsub_pp(t8,t9,_j,pd1,pd9);                          \
   sse2_load_addsub_pp(t10,t11,_j,pd5,pd13);                       \
   sse2_addsub(t8,t10);                                            \
   sse2_mul_1_4_##_S##_addsub(t9,t11);                             \
   sse2_load_addsub_pp(t12,t13,_j,pd3,pd11);                       \
   sse2_load_addsub_pp(t14,t15,_j,pd7,pd15);                       \
   sse2_addsub(t12,t14);                                           \
   sse2_mul_1_4_##_S##_addsub(t13,t15);                            \
   sse2_addsub(t8,t12);                                            \
   sse2_mul_1_8_##_S##_addsub(t9,t13);                             \
   sse2_mul_1_4_##_S##_addsub(t10,t14);                            \
   sse2_mul_3_8_##_S##_addsub(t11,t15);                            \
	                                                           \
   sse2_addsub_store_p(_j,pd0,pd8,t0,t8);                          \
   sse2_muladdsub_store_p(_j,pd1,pd9,t1,t9,MM_##_S##_1_16);        \
   sse2_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);         \
   sse2_muladdsub_store_p(_j,pd3,pd11,t3,t11,MM_##_S##_3_16);      \
   sse2_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);         \
   sse2_muladdsub_store_p(_j,pd5,pd13,t5,t13,MM_##_S##_5_16);      \
   sse2_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);         \
   sse2_muladdsub_store_p(_j,pd7,pd15,t7,t15,MM_##_S##_7_16);


/* SSE2 versioned */
/* Mode 0: read =0 (standard)  write= 0 (standard) */
void radixmm0_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_16_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm0_16_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm0_16_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm0_16_twd_pp(F,jj);
            }
        }
    }
}

void radixmm0_16_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_16_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;
  nc=n+n0;
  for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
    {
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radixmm0_16_notwd_pp(F,jj);
        }
    }
}


void radixmm0_16_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm0_16_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_B;
  get_pads_16;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm0_16_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm0_16_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm0_16_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm0_16_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm0_16_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm0_16_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm0_16_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm0_16_twd_pp(B,jj);
        }
      i += bigpad;
    }

# endif
}

/* Mode 1: read =0 (standard)  write= 1 (fast-sse2) */
void radixmm1_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_16_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm1_16_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm1_16_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm1_16_twd_pp(F,jj);
            }
        }
    }
}

void radixmm1_16_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_16_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;
  nc=n+n0;
  for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
    {
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radixmm1_16_notwd_pp(F,jj);
        }
    }
}

void radixmm1_16_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm1_16_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_B;
  get_pads_16;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm1_16_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm1_16_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm1_16_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm1_16_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm1_16_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm1_16_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm1_16_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm1_16_twd_pp(B,jj);
        }
      i += bigpad;
    }

# endif
}

/* Mode 2: read =1 (fast-sse2)  write= 0 (standard) */
void radixmm2_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_16_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm2_16_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm2_16_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radixmm2_16_twd_pp(F,jj);
            }
        }
    }
}

void radixmm2_16_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm2_16_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_B;
  get_pads_16;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm2_16_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm2_16_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm2_16_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm2_16_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm2_16_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm2_16_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm2_16_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm2_16_twd_pp(B,jj);
        }
      i += bigpad;
    }

# endif
}

/* Mode 3: read =1 (fast-sse2)  write= 0 (fast-sse2) */
void radixmm3_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm3_16_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_F;
  get_pads_16;

  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radixmm3_16_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj = addr(j<<1);
              radixmm3_16_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          prefetch_data_16( i );
          get_twiddle_factors_sse2_16_doubled;
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              /*BEGIN_TIME;*/
              radixmm3_16_twd_pp(F,jj);
              /*END_TIME;*/
            }
        }
    }
}

void radixmm3_16_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  Y__M128D t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
  Y__M128D tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,pd14,pd15;
  y_size_t i,j,jj,nc;
# ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmm3_16_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  set_sse2_16_constants_B;
  get_pads_16;

  nc=n0+n;

# if !defined(Y_SAVE_TRIGS)

  for(i=n0;i<nc;i+=bigpad)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm3_16_twd_pp(B,jj);
      j=i+2;
      px += Y_STEP;
      for( ;j<(pad+i); j+=2, px += 2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm3_16_twd_pp(B,jj);
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
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm3_16_twd_pp(B,jj);
      j=i+2;
      jj1 = addr(i1<<1);
      radixmm3_16_twd_pp(B,jj1);
      j1 = i1 + 2;
      px += Y_STEP;
      for(;j<(pad+i);j+=2,j1+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm3_16_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radixmm3_16_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      px=tw;
      get_twiddle_factors_sse2_16_one;
      jj=addr(i<<1);
      radixmm3_16_twd_pp(B,jj);
      px += Y_STEP;
      j=i+2;
      for(;j<(pad+i);j+=2,px+=2*Y_STEP)
        {
          get_twiddle_factors_sse2_16;
          jj=addr(j<<1);
          radixmm3_16_twd_pp(B,jj);
        }
      i += bigpad;
    }

# endif
}

#else /* defined Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radixmm_16_void( void )
{
  printf("The routine radixmm_16_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif
