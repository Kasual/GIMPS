/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006 Guillermo Ballester Valor
 
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
/*#define Y_DEBUG*/
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define Y_SAVE_TRIGS
#define NDEBUG1
#ifdef Y_ITANIUM
# ifndef Y_PRELOAD_TRIGS
#  define Y_PRELOAD_TRIGS
# endif
# ifndef Y_PRELOAD_DATA
#  define Y_PRELOAD_DATA
# endif
#  ifdef Y_LONG_MACROS
#   undef Y_LONG_MACROS
#  endif
#endif

#define RADIX 16
/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* asign the correct values to every pad */


#define get_pads(_n0)                  \
  jj = pad;                            \
  pd0 = d;                             \
  prefetch_data( pd0, addr((_n0<<1))); \
  pd1 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd1, addr((_n0<<1))); \
  pd2 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd2, addr((_n0<<1))); \
  pd3 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd3, addr((_n0<<1))); \
  pd4 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd4, addr((_n0<<1))); \
  pd5 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd5, addr((_n0<<1))); \
  pd6 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd6, addr((_n0<<1))); \
  pd7 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd7, addr((_n0<<1))); \
  pd8 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd8, addr((_n0<<1))); \
  pd9 = d + addr((jj<<1));             \
  jj += pad;                           \
  prefetch_data( pd9, addr((_n0<<1))); \
  pd10 = d + addr((jj<<1));            \
  jj += pad;                           \
  prefetch_data( pd10, addr((_n0<<1)));\
  pd11 = d + addr((jj<<1));            \
  jj += pad;                           \
  prefetch_data( pd11, addr((_n0<<1)));\
  pd12 = d + addr((jj<<1));            \
  jj += pad;                           \
  prefetch_data( pd12, addr((_n0<<1)));\
  pd13 = d + addr((jj<<1));            \
  jj += pad;                           \
  prefetch_data( pd13, addr((_n0<<1)));\
  pd14 = d + addr((jj<<1));            \
  jj += pad;                           \
  prefetch_data( pd14, addr((_n0<<1)));\
  pd15 = d + addr((jj<<1));            \
  prefetch_data( pd15, addr((_n0<<1)));

#ifdef Y_PRELOAD_DATA
# define radix_16_init_preload(_tp,_pd,_index)                        \
  _tp##0r = *(_pd##0 + _index); _tp##0i = *(_pd##0 + _index + 1);     \
  _tp##8r = *(_pd##8 + _index); _tp##8i = *(_pd##8 + _index + 1);     \
  _tp##4r = *(_pd##4 + _index); _tp##4i = *(_pd##4 + _index + 1);     \
  _tp##12r = *(_pd##12 + _index); _tp##12i = *(_pd##12 + _index + 1); \
  _tp##2r = *(_pd##2 + _index); _tp##2i = *(_pd##2 + _index + 1);     \
  _tp##10r = *(_pd##10 + _index); _tp##10i = *(_pd##10 + _index + 1); \
  _tp##6r = *(_pd##6 + _index); _tp##6i = *(_pd##6 + _index + 1);     \
  _tp##14r = *(_pd##14 + _index); _tp##14i = *(_pd##14 + _index + 1); \
  _tp##1r = *(_pd##1 + _index); _tp##1i = *(_pd##1 + _index + 1);     \
  _tp##9r = *(_pd##9 + _index); _tp##9i = *(_pd##9 + _index + 1);     \
  _tp##5r = *(_pd##5 + _index); _tp##5i = *(_pd##5 + _index + 1);     \
  _tp##13r = *(_pd##13 + _index); _tp##13i = *(_pd##13 + _index + 1); \
  _tp##3r = *(_pd##3 + _index); _tp##3i = *(_pd##3 + _index + 1);     \
  _tp##11r = *(_pd##11 + _index); _tp##11i = *(_pd##11 + _index + 1); \
  _tp##7r = *(_pd##7 + _index); _tp##7i = *(_pd##7 + _index + 1);     \
  _tp##15r = *(_pd##15 + _index); _tp##15i = *(_pd##15 + _index + 1); \


#endif




/* get the basic twiddle from memory and computes its powers */
#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# include "yminimum.h"

# ifdef Y_PRELOAD_TRIGS

#  define trig_16_load_init(_px)              \
          tpw1r = *(_px );                    \
          tpw1i = *(_px + 1);

#  define trig_16_preload(_px)                \
          tpw1r = *(_px + Y_STEP);            \
          tpw1i = *(_px + Y_STEP + 1);

#  define get_twiddle_factors_16_preload      \
          tw1r = tpw1r; tw1i = tpw1i;         \
          cplx_trig_min_square(tw2,twp1);     \
          cplx_trig_min_square(tw4,tw2);      \
          cplx_trig_min_square(tw8,tw4);      \
          cplx_divmul(tw3,tw5,tw4,tpw1);      \
          cplx_divmul(tw7,tw9,tw8,tpw1);      \
          cplx_divmul(tw6,tw10,tw8,tw2);      \
          cplx_trig_min_square(tw12,tw6);     \
          cplx_divmul(tw11,tw13,tw12,tpw1);   \
          cplx_trig_min_square(tw14,tw7);     \
          cplx_trig_min_mul(tw15,tw7,tw8);    \
          trig_8_preload(px);
# else

#  define get_twiddle_factors_16              \
          cplx_trig_min_load(tw1,px);         \
          cplx_trig_min_square(tw2,tw1);      \
          cplx_trig_min_square(tw4,tw2);      \
          cplx_trig_min_square(tw8,tw4);      \
          cplx_divmul(tw3,tw5,tw4,tw1);       \
          cplx_divmul(tw7,tw9,tw8,tw1);       \
          cplx_divmul(tw6,tw10,tw8,tw2);      \
          cplx_trig_min_square(tw12,tw6);     \
          cplx_divmul(tw11,tw13,tw12,tw1);    \
          cplx_trig_min_square(tw14,tw7);     \
          cplx_trig_min_mul(tw15,tw7,tw8);    \
          prefetch_p_trig(px+2);

#   define get_twiddle_factors_16_inloop  get_twiddle_factors_16
#   define get_twiddle_factors_16_outloop get_twiddle_factors_16

# endif

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 30

# ifdef Y_PRELOAD_TRIGS
#  define trig_16_load_init(_px)                          \
          tpw1r = *(_px);      tpw1i = *(_px + 1);        \
          tpw2r = *(_px + 2);  tpw2i = *(_px + 3);        \
          tpw3r = *(_px + 4);  tpw3i = *(_px + 5);        \
          tpw4r = *(_px + 6);  tpw4i = *(_px + 7);        \
          tpw5r = *(_px + 8);  tpw5i = *(_px + 9);        \
          tpw6r = *(_px + 10);  tpw6i = *(_px + 11);      \
          tpw7r = *(_px + 12);  tpw7i = *(_px + 13);      \
          tpw8r = *(_px + 14);  tpw8i = *(_px + 15);      \
          tpw9r = *(_px + 16);  tpw9i = *(_px + 17);      \
          tpw10 = *(_px + 18);  tpw10i = *(_px + 19);     \
          tpw11 = *(_px + 20);  tpw11i = *(_px + 21);     \
          tpw12 = *(_px + 22);  tpw12i = *(_px + 23);      \
          tpw13 = *(_px + 24);  tpw13i = *(_px + 25);      \
          tpw14 = *(_px + 26);  tpw14i = *(_px + 27);      \
          tpw15 = *(_px + 28);  tpw15i = *(_px + 29);

#  define trig_16_preload(_px)                                              \
          tpw1r = *(_px + Y_STEP);      tpw1i = *(_px + Y_STEP + 1);        \
          tpw2r = *(_px + Y_STEP + 2);  tpw2i = *(_px + Y_STEP + 3);        \
          tpw3r = *(_px + Y_STEP + 4);  tpw3i = *(_px + Y_STEP + 5);        \
          tpw4r = *(_px + Y_STEP + 6);  tpw4i = *(_px + Y_STEP + 7);        \
          tpw5r = *(_px + Y_STEP + 8);  tpw5i = *(_px + Y_STEP + 9);        \
          tpw6r = *(_px + Y_STEP + 10);  tpw6i = *(_px + Y_STEP + 11);      \
          tpw7r = *(_px + Y_STEP + 12);  tpw7i = *(_px + Y_STEP + 13);      \
          tpw8r = *(_px + Y_STEP + 14);  tpw8i = *(_px + Y_STEP + 15);      \
          tpw9r = *(_px + Y_STEP + 16);  tpw9i = *(_px + Y_STEP + 17);      \
          tpw10 = *(_px + Y_STEP + 18);  tpw10i = *(_px + Y_STEP + 19);     \
          tpw11 = *(_px + Y_STEP + 20);  tpw11i = *(_px + Y_STEP + 21);     \
          tpw12 = *(_px + Y_STEP + 22);  tpw12i = *(_px + Y_STEP + 23);     \
          tpw13 = *(_px + Y_STEP + 24);  tpw13i = *(_px + Y_STEP + 25);     \
          tpw14 = *(_px + Y_STEP + 26);  tpw14i = *(_px + Y_STEP + 27);     \
          tpw15 = *(_px + Y_STEP + 28);  tpw15i = *(_px + Y_STEP + 29);     \


#  define get_twiddle_factors_16_preload      \
          tw1r = tpw1r; tw1i = tpw1i;         \
          tw2r = tpw2r; tw2i = tpw2i;         \
          tw3r = tpw3r; tw3i = tpw3i;         \
          tw4r = tpw4r; tw4i = tpw4i;         \
          tw5r = tpw5r; tw5i = tpw5i;         \
          tw6r = tpw6r; tw6i = tpw6i;         \
          tw7r = tpw7r; tw7i = tpw7i;         \
          tw8r = tpw8r; tw8i = tpw8i;         \
          tw9r = tpw9r; tw9i = tpw9i;         \
          tw10r = tpw10r; tw10i = tpw10i;     \
          tw11r = tpw11r; tw11i = tpw11i;     \
          tw12r = tpw12r; tw12i = tpw12i;     \
          tw13r = tpw13r; tw13i = tpw13i;     \
          tw14r = tpw14r; tw14i = tpw14i;     \
          tw15r = tpw15r; tw15i = tpw15i;     \
          tw16r = tpw16r; tw16i = tpw16i;     \
          trig_16_preload(px);

# else

#  if defined(Y_MAXIMUM) && (Y_TARGET == 0)

#   define get_twiddle_factors_16_inloop  /*  */

#   define get_twiddle_factors_16_outloop                  \
          prefetch_data_trig(px,Y_STEP);                   \
          prefetch_data_trig(px,Y_STEP + Y_CACHE_LINE);    \
          prefetch_data_trig(px,Y_STEP + 2*Y_CACHE_LINE);  \
          prefetch_data_trig(px,Y_STEP + 3*Y_CACHE_LINE);  \
          prefetch_data_trig(px,Y_STEP + 4*Y_CACHE_LINE);

#   define asim_prefetch_F_trig(_j)   /* */



#   define asim_prefetch_B_trig(_j)          \
          prefetch_data_trig(px,Y_STEP + _j*Y_CACHE_LINE);

#  else


#   define get_twiddle_factors_16                       \
          tw1r = *(px);      tw1i = *(px + 1);         \
          tw2r = *(px + 2);  tw2i = *(px + 3);         \
          tw3r = *(px + 4);  tw3i = *(px + 5);         \
          tw4r = *(px + 6);  tw4i = *(px + 7);         \
          tw5r = *(px + 8);  tw5i = *(px + 9);         \
          tw6r = *(px + 10);  tw6i = *(px + 11);       \
          tw7r = *(px + 12);  tw7i = *(px + 13);       \
          tw8r = *(px + 14);  tw8i = *(px + 15);       \
          tw9r = *(px + 16);  tw9i = *(px + 17);       \
          tw10r = *(px + 18);  tw10i = *(px + 19);     \
          tw11r = *(px + 20);  tw11i = *(px + 21);     \
          tw12r = *(px + 22);  tw12i = *(px + 23);     \
          tw13r = *(px + 24);  tw13i = *(px + 25);     \
          tw14r = *(px + 26);  tw14i = *(px + 27);     \
          tw15r = *(px + 28);  tw15i = *(px + 29);     \
          prefetch_data_trig(px,30);                   \
          prefetch_data_trig(px,30 + Y_CACHE_LINE);    \
          prefetch_data_trig(px,30 + 2*Y_CACHE_LINE);  \
          prefetch_data_trig(px,30 + 3*Y_CACHE_LINE);  \
          prefetch_data_trig(px,30 + 4*Y_CACHE_LINE);

#   define get_twiddle_factors_16_inloop  get_twiddle_factors_16
#   define get_twiddle_factors_16_outloop get_twiddle_factors_16

#  endif
# endif

#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 10

# ifdef Y_PRELOAD_TRIGS

#  define trig_16_load_init(_px)     \
          tpw1r = *(_px);            \
          tpw1i = *(_px + 1);        \
          tpw2r = *(_px + 2);        \
          tpw2i = *(_px + 3);        \
          tpw4r = *(_px + 4);        \
          tpw4i = *(_px + 5);        \
          tpw8r = *(_px + 6);        \
          tpw8i = *(_px + 7);        \
          tpw13r = *(_px + 8);       \
          tpw13i = *(_px + 9);       \


#  define trig_16_preload(_px)                \
          tpw1r = *(_px + Y_STEP);            \
          tpw1i = *(_px + Y_STEP + 1);        \
          tpw2r = *(_px + Y_STEP + 2);        \
          tpw2i = *(_px + Y_STEP + 3);        \
          tpw4r = *(_px + Y_STEP + 4);        \
          tpw4i = *(_px + Y_STEP + 5);        \
          tpw8r = *(_px + Y_STEP + 6);        \
          tpw8i = *(_px + Y_STEP + 7);        \
          tpw13r = *(_px + Y_STEP + 8);       \
          tpw13i = *(_px + Y_STEP + 9);       \


#  define get_twiddle_factors_16_preload      \
          tw1r = tpw1r; tw1i = tpw1i;         \
          tw2r = tpw2r; tw2i = tpw2i;         \
          tw4r = tpw4r; tw4i = tpw4i;         \
          tw8r = tpw8r; tw8i = tpw8i;         \
          tw13r = tpw13r; tw13i = tpw13i;     \
          cplx_divmul(tw3,tw5,tpw4,tpw1);     \
          cplx_divmul(tw6,tw10,tpw8,tpw2);    \
	  cplx_divmul(tw7,tw9,tpw8,tpw1);     \
	  cplx_divmul(tw11,tw15,tpw13,tpw2);  \
	  cplx_divmul(tw12,tw14,tpw13,tpw1);  \
          trig_16_preload(px);

# else

#  define get_twiddle_factors_16                   \
          cplx_data_trig_load_3(tw1,tw2,tw4,px);   \
          cplx_data_trig_load_2(tw8,tw13,px+6);    \
          cplx_divmul(tw3,tw5,tw4,tw1);            \
          cplx_divmul(tw6,tw10,tw8,tw2);           \
	  cplx_divmul(tw7,tw9,tw8,tw1);            \
	  cplx_divmul(tw11,tw15,tw13,tw2);         \
	  cplx_divmul(tw12,tw14,tw13,tw1);         \
          prefetch_data_trig(px,Y_STEP );          \
          prefetch_data_trig(px,Y_STEP + Y_CACHE_LINE);  \
          prefetch_data_trig(px,Y_STEP + 2*Y_CACHE_LINE);

#   define get_twiddle_factors_16_inloop  get_twiddle_factors_16
#   define get_twiddle_factors_16_outloop get_twiddle_factors_16

# endif
#endif

/*
	Get data from memory, mul by twiddle trig factors,
        make a complex 16-length FFT and store the data to memory 
*/

#if defined(Y_MAXIMUM) && (YTARGET == 0)

# define radix_16_twd_pp(_S,_j)                                   \
	  cplx_load_muladdsub_pp_ind(t0,t1,px,8,_j,pd0,pd8);      \
	  cplx_load_mulmuladdsub_pp_ind(t2,t3,px,4,12,_j,pd4,pd12);  \
          asim_prefetch_##_S##_trig(0);                              \
	  cplx_addsub(t0,t2);                                     \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                      \
	  cplx_load_mulmuladdsub_pp_ind(t4,t5,px,2,10,_j,pd2,pd10);  \
	  cplx_load_mulmuladdsub_pp_ind(t6,t7,px,6,14,_j,pd6,pd14);  \
          asim_prefetch_##_S##_trig(1);                           \
	  cplx_addsub(t4,t6);                                     \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                      \
	  cplx_addsub(t0,t4);                                     \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                      \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                      \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                      \
	  cplx_load_mulmuladdsub_pp_ind(t8,t9,px,1,9,_j,pd1,pd9);    \
	  cplx_load_mulmuladdsub_pp_ind(t10,t11,px,5,13,_j,pd5,pd13);\
          asim_prefetch_##_S##_trig(2);                             \
	  cplx_addsub(t8,t10);                                    \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                     \
	  cplx_load_mulmuladdsub_pp_ind(t12,t13,px,3,11,_j,pd3,pd11);\
	  cplx_load_mulmuladdsub_pp_ind(t14,t15,px,7,15,_j,pd7,pd15);\
          asim_prefetch_##_S##_trig(3);                           \
	  cplx_addsub(t12,t14);                                   \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                    \
	  cplx_addsub(t8,t12);                                    \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                     \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                    \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                    \
	                                                          \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                  \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);     \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10); \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);   \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12); \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);   \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14); \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);

#else


# define radix_16_twd_pp(_S,_j)                                   \
	  cplx_load_muladdsub_pp(t0,t1,tw8,_j,pd0,pd8);           \
	  cplx_load_mulmuladdsub_pp(t2,t3,tw4,tw12,_j,pd4,pd12);  \
	  cplx_addsub(t0,t2);                                     \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                      \
	  cplx_load_mulmuladdsub_pp(t4,t5,tw2,tw10,_j,pd2,pd10);  \
	  cplx_load_mulmuladdsub_pp(t6,t7,tw6,tw14,_j,pd6,pd14);  \
	  cplx_addsub(t4,t6);                                     \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                      \
	  cplx_addsub(t0,t4);                                     \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                      \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                      \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                      \
	  cplx_load_mulmuladdsub_pp(t8,t9,tw1,tw9,_j,pd1,pd9);    \
	  cplx_load_mulmuladdsub_pp(t10,t11,tw5,tw13,_j,pd5,pd13);\
	  cplx_addsub(t8,t10);                                    \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                     \
	  cplx_load_mulmuladdsub_pp(t12,t13,tw3,tw11,_j,pd3,pd11);\
	  cplx_load_mulmuladdsub_pp(t14,t15,tw7,tw15,_j,pd7,pd15);\
	  cplx_addsub(t12,t14);                                   \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                    \
	  cplx_addsub(t8,t12);                                    \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                     \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                    \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                    \
	                                                          \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                  \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);     \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10); \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);   \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12); \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);   \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14); \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);



#endif

# define radix_16_twd_pp_no_fetch(_S,_j)                                   \
	  cplx_load_muladdsub_pp_no_fetch(t0,t1,tw8,_j,pd0,pd8);           \
	  cplx_load_mulmuladdsub_pp_no_fetch(t2,t3,tw4,tw12,_j,pd4,pd12);  \
	  cplx_addsub(t0,t2);                                              \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                               \
	  cplx_load_mulmuladdsub_pp_no_fetch(t4,t5,tw2,tw10,_j,pd2,pd10);  \
	  cplx_load_mulmuladdsub_pp_no_fetch(t6,t7,tw6,tw14,_j,pd6,pd14);  \
	  cplx_addsub(t4,t6);                                              \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                               \
	  cplx_addsub(t0,t4);                                              \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                               \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                               \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                               \
	  cplx_load_mulmuladdsub_pp_no_fetch(t8,t9,tw1,tw9,_j,pd1,pd9);    \
	  cplx_load_mulmuladdsub_pp_no_fetch(t10,t11,tw5,tw13,_j,pd5,pd13);\
	  cplx_addsub(t8,t10);                                             \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                              \
	  cplx_load_mulmuladdsub_pp_no_fetch(t12,t13,tw3,tw11,_j,pd3,pd11);\
	  cplx_load_mulmuladdsub_pp_no_fetch(t14,t15,tw7,tw15,_j,pd7,pd15);\
	  cplx_addsub(t12,t14);                                            \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                             \
	  cplx_addsub(t8,t12);                                             \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                              \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                             \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                             \
	                                                                   \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                           \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);              \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);          \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);            \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);          \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);            \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);          \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);

# define radix_16_twd_preload(_S,_j,_jj)                                    \
	  cplx_preload_muladdsub(t0,t1,tw8,_jj,pd0,pd8,tp0,tp8);            \
	  cplx_preload_mulmuladdsub(t2,t3,tw4,tw12,_jj,pd4,pd12,tp4,tp12);  \
	  cplx_preload_mulmuladdsub(t4,t5,tw2,tw10,_jj,pd2,pd10,tp2,tp10);  \
	  cplx_preload_mulmuladdsub(t6,t7,tw6,tw14,_jj,pd6,pd14,tp6,tp14);  \
	  cplx_preload_mulmuladdsub(t8,t9,tw1,tw9,_jj,pd1,pd9,tp1,tp9);     \
	  cplx_preload_mulmuladdsub(t10,t11,tw5,tw13,_jj,pd5,pd13,tp5,tp13);\
	  cplx_preload_mulmuladdsub(t12,t13,tw3,tw11,_jj,pd3,pd11,tp3,tp11);\
	  cplx_preload_mulmuladdsub(t14,t15,tw7,tw15,_jj,pd7,pd15,tp7,tp15);\
	  cplx_addsub(t0,t2);                                               \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                                \
	  cplx_addsub(t4,t6);                                               \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                                \
	  cplx_addsub(t8,t10);                                              \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                               \
	  cplx_addsub(t12,t14);                                             \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                              \
	  cplx_addsub(t0,t4);                                               \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                                \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                                \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                                \
	  cplx_addsub(t8,t12);                                              \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                               \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                              \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                              \
	                                                                    \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                            \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);           \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);             \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);           \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);             \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);           \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);             \
          _j = _jj;


/*
	Get data from memory, make a complex 16-length FFT and
	store the data to memory 
*/
# define radix_16_notwd_pp(_S,_j)                                   \
	  cplx_load_addsub_pp(t0,t1,_j,pd0,pd8);                    \
	  cplx_load_addsub_pp(t2,t3,_j,pd4,pd12);                   \
	  cplx_addsub(t0,t2);                                       \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                        \
	  cplx_load_addsub_pp(t4,t5,_j,pd2,pd10);                   \
	  cplx_load_addsub_pp(t6,t7,_j,pd6,pd14);                   \
	  cplx_addsub(t4,t6);                                       \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                        \
	  cplx_addsub(t0,t4);                                       \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                        \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                        \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                        \
	  cplx_load_addsub_pp(t8,t9,_j,pd1,pd9);                    \
	  cplx_load_addsub_pp(t10,t11,_j,pd5,pd13);                 \
	  cplx_addsub(t8,t10);                                      \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                       \
	  cplx_load_addsub_pp(t12,t13,_j,pd3,pd11);                 \
	  cplx_load_addsub_pp(t14,t15,_j,pd7,pd15);                 \
	  cplx_addsub(t12,t14);                                     \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                      \
	  cplx_addsub(t8,t12);                                      \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                       \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                      \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                      \
	                                                            \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                    \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);       \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);   \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);     \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);   \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);     \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);   \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);

# define radix_16_notwd_pp_no_fetch(_S,_j)                          \
	  cplx_load_addsub_pp_no_fetch(t0,t1,_j,pd0,pd8);           \
	  cplx_load_addsub_pp_no_fetch(t2,t3,_j,pd4,pd12);          \
	  cplx_addsub(t0,t2);                                       \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                        \
	  cplx_load_addsub_pp_no_fetch(t4,t5,_j,pd2,pd10);          \
	  cplx_load_addsub_pp_no_fetch(t6,t7,_j,pd6,pd14);          \
	  cplx_addsub(t4,t6);                                       \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                        \
	  cplx_addsub(t0,t4);                                       \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                        \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                        \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                        \
	  cplx_load_addsub_pp_no_fetch(t8,t9,_j,pd1,pd9);           \
	  cplx_load_addsub_pp_no_fetch(t10,t11,_j,pd5,pd13);        \
	  cplx_addsub(t8,t10);                                      \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                       \
	  cplx_load_addsub_pp_no_fetch(t12,t13,_j,pd3,pd11);        \
	  cplx_load_addsub_pp_no_fetch(t14,t15,_j,pd7,pd15);        \
	  cplx_addsub(t12,t14);                                     \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                      \
	  cplx_addsub(t8,t12);                                      \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                       \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                      \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                      \
	                                                            \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                    \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);       \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);   \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);     \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);   \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);     \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);   \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);

# define radix_16_notwd_preload(_S,_j,_jj)                                  \
	  cplx_preload_addsub(t0,t1,_jj,pd0,pd8,tp0,tp8);                   \
	  cplx_preload_addsub(t2,t3,_jj,pd4,pd12,tp4,tp12);                 \
	  cplx_preload_addsub(t4,t5,_jj,pd2,pd10,tp2,tp10);                 \
	  cplx_preload_addsub(t6,t7,_jj,pd6,pd14,tp6,tp14);                 \
	  cplx_preload_addsub(t8,t9,_jj,pd1,pd9,tp1,tp9);                   \
	  cplx_preload_addsub(t10,t11,_jj,pd5,pd13,tp5,tp13);               \
	  cplx_preload_addsub(t12,t13,_jj,pd3,pd11,tp3,tp11);               \
	  cplx_preload_addsub(t14,t15,_jj,pd7,pd15,tp7,tp15);               \
	  cplx_addsub(t0,t2);                                               \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                                \
	  cplx_addsub(t4,t6);                                               \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                                \
	  cplx_addsub(t8,t10);                                              \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                               \
	  cplx_addsub(t12,t14);                                             \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                              \
	  cplx_addsub(t0,t4);                                               \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                                \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                                \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                                \
	  cplx_addsub(t8,t12);                                              \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                               \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                              \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                              \
	                                                                    \
	  cplx_addsub_store_p(_j,pd0,pd8,t0,t8);                            \
	  cplx_muladdsub_store_p(_j,pd1,pd9,t1,t9,_S##_1_16);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd2,pd10,t2,t10);           \
	  cplx_muladdsub_store_p(_j,pd3,pd11,t3,t11,_S##_3_16);             \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd4,pd12,t4,t12);           \
	  cplx_muladdsub_store_p(_j,pd5,pd13,t5,t13,_S##_5_16);             \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd6,pd14,t6,t14);           \
	  cplx_muladdsub_store_p(_j,pd7,pd15,t7,t15,_S##_7_16);             \
          _j = _jj;




#if (Y_AVAL > 3)

/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = Scrambled Twidle factors. Not padded.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
 
   This is a version modified to support multithreaded work. In the standard 
   version, it is supossed that n0 was divisible by bigpad, i.e. it should start
   at the beginnig of a block. 
 
   NOW n0 is supossed to be even a divisor of bigpad. 
   tw[] is supossed to point at the begining of a block
 
*************************************************************************/

void radixmp_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
#if !defined(Y_MAXIMUM) || (Y_TARGET != 0)

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
#endif
#ifdef Y_PRELOAD_DATA

  y_limb_t tp0r,tp0i,tp1r,tp1i,tp2r,tp2i,tp3r,tp3i,tp4r,tp4i,tp5r,tp5i,
  tp6r,tp6i,tp7r,tp7i,tp8r,tp8i,tp9r,tp9i,tp10r,tp10i,tp11r,tp11i,tp12r,
  tp12i,tp13r,tp13i,tp14r,tp14i,tp15r,tp15i;
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,
  tpw6r,tpw6i,tpw7r,tpw7i,tpw8r,tpw8i,tpw9r,tpw9i,tpw10r,tpw10i,
  tpw11r,tpw11i,tpw12r,tpw12i,tpw13r,tpw13i,tpw14r,tpw14i,tpw15r,tpw15i;
# elif defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw4r,tpw4i,tpw8r,tpw8i,tpw13r,tpw13i;
# endif
#endif

  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,
  pd14,pd15;
  y_size_t i,j,jj;
#ifdef Y_PRELOAD_DATA

  y_size_t kk;
#endif

  y_size_t done, mask, ioff, ibase, n00;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();
# ifndef NDEBUG1

  printf(" Radixmp_16_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  bigpad = pad << 4;
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad);
  ioff = (n0 & mask) >> 4;
  n00 = ibase * bigpad + ioff;
  get_pads(n00);

  /*printf("ibase = %d ioff= %d\n",ibase,ioff);*/

  /* If ibase == 0 first elements has twiddle factors equal to one */
  if(ibase == 0)
    {
#ifdef Y_PRELOAD_DATA
      radix_16_init_preload(tp,pd,0);
      kk=0;
#endif

      for(j = ioff; j < pad && done < n; j++)
        {
#ifdef Y_PRELOAD_DATA
          jj = addr((j+1)<<1);
          radix_16_notwd_preload(F,kk,jj);
          done += RADIX;
#else

          jj = addr(j<<1);
          radix_16_notwd_pp(F,jj);
          done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

          j++;
          jj += 2;
          radix_16_notwd_pp_no_fetch(F,jj);
          done += RADIX;
#  if Y_CACHE_LINE == 8

          j += 2;
          jj += 2;
          radix_16_notwd_pp_no_fetch(F,jj);
          done += RADIX;
          jj += 2;
          radix_16_notwd_pp_no_fetch(F,jj);
          done += RADIX;
#  endif
# endif
#endif

        }
      if (done == n)
        return;
      px = tw + Y_STEP;
#ifdef Y_PRELOAD_TRIGS

      trig_16_load_init(px);
#endif

      for(i = bigpad; done < n;i += bigpad, px += Y_STEP)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_16_preload;
#else

          get_twiddle_factors_16_outloop;
#endif
#ifdef Y_PRELOAD_DATA

          kk = addr(i<<1);
          radix_16_init_preload(tp,pd,kk);
          done += RADIX;
#endif

          for(j = i; j < (pad + i) && done < n; j++)
            {
#ifdef Y_PRELOAD_DATA
              jj = addr((j+1)<<1);
              radix_16_twd_preload(F,kk,jj);
              done += RADIX;
#else

              jj = addr(j<<1);
              radix_16_twd_pp(F,jj);
              done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  if Y_CACHE_LINE == 8

              j += 2;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  endif
# endif
#endif

            }
        }
    }
  else
    {
      px = tw + Y_STEP * ibase;
#ifdef Y_PRELOAD_TRIGS

      trig_16_load_init(px);
#endif

      for(i = ibase * bigpad; done < n; i += bigpad, px += Y_STEP)
        {
          /*BEGIN_TIME;*/
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_16_preload;
#else

          get_twiddle_factors_16_outloop;
#endif
#ifdef Y_PRELOAD_DATA

          kk = addr((i + ioff)<<1);
          radix_16_init_preload(tp, pd, kk);
#endif
          /*END_TIME;*/
          for(j = i + ioff; j < (pad + i) && done < n; j++)
            {
#ifdef Y_PRELOAD_DATA
              jj = addr((j+1)<<1);
              radix_16_twd_preload(F,kk,jj);
              done += RADIX;
#else

              jj = addr(j<<1);
              /*BEGIN_TIME;*/
              radix_16_twd_pp(F,jj);
              done += RADIX;
              /*END_TIME;*/
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  if Y_CACHE_LINE == 8

              j += 2;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
              jj += 2;
              radix_16_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  endif
# endif
#endif

            }
          ioff  = 0;
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
 
   This is a version modified to support multithreaded work. In the standard 
   version, it is supossed that n0 was divisible by bigpad, i.e. it should start
   at the beginnig of a block. 
 
   NOW n0 is supossed to be even a divisor of bigpad. 
   tw[] is supossed to point at the begining of a block
 
*************************************************************************/

void radixmp_16_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i;
#if !defined(Y_MAXIMUM) || (Y_TARGET != 0)

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,
  tpw6r,tpw6i,tpw7r,tpw7i,tpw8r,tpw8i,tpw9r,tpw9i,tpw10r,tpw10i,
  tpw11r,tpw11i,tpw12r,tpw12i,tpw13r,tpw13i,tpw14r,tpw14i,tpw15r,tpw15i;
# elif defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw4r,tpw4i,tpw8r,tpw8i,tpw13r,tpw13i;
# endif
#endif

  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,
  pd14,pd15;
  y_size_t i,j,jj;
#ifdef Y_SAVE_TRIGS
  /*y_size_t i1,j1,jj1;*/
#endif

  y_size_t done, mask, ioff, ibase, n00;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmp_16_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  bigpad = pad << 4;
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad) * bigpad;
  ioff = (n0 & mask) >> 4;
  n00 = ibase + ioff;

  get_pads(n00);


  /*printf("ibase = %d ioff= %d\n",ibase,ioff);*/

  for(i = ibase; done < n; i += bigpad)
    {
      if (ioff == 0)
        {
          jj = addr(i<<1);
          radix_16_notwd_pp(B,jj);
          done += RADIX;
          px = tw;
          j = i + 1;
#ifdef Y_PRELOAD_TRIGS

          trig_16_load_init(px);
#endif
#if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4
# ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
# else

          get_twiddle_factors_16_inloop;
# endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j++;
# if Y_CACHE_LINE == 8
#  ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
#  else

          get_twiddle_factors_16_inloop;
#  endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
#  ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
#  else

          get_twiddle_factors_16_inloop;
#  endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j += 2;
# endif
#endif

        }
      else
        {
          px = tw + (ioff - 1) * Y_STEP;
#ifdef Y_PRELOAD_TRIGS

          trig_16_load_init(px);
#endif

          j = i + ioff;
        }

      for(;j < (pad + i) && done < n; j++, px += Y_STEP)
        {
          /*BEGIN_TIME;*/
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_16_preload;
#else

          get_twiddle_factors_16_inloop;
#endif

          jj = addr(j<<1);
          radix_16_twd_pp(B,jj);
          done += RADIX;
#if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4

          px += Y_STEP;
          j++;
# ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
# else

          get_twiddle_factors_16_inloop;
# endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
# if Y_CACHE_LINE == 8

          px += Y_STEP;
#  ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
#  else

          get_twiddle_factors_16_inloop;
#  endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j += 2;
#  ifdef Y_PRELOAD_TRIGS

          get_twiddle_factors_16_preload;
#  else

          get_twiddle_factors_16_inloop;
#  endif

          jj += 2;
          radix_16_twd_pp_no_fetch(B,jj);
          done += RADIX;
# endif
#endif
          /*END_TIME;*/
        }
      ioff = 0;
    }
}

#else

void radixmp_16_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radixmp_16_dif , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_16.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

void radixmp_16_dit(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radixmp_16_dit , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_16.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

#endif
/*$Id$*/

















