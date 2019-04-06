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
#define NDEBUG1
#define Y_SAVE_TRIGS
#ifdef  Y_ITANIUM
#  include "yia64.h"
#  ifndef Y_PRELOAD_TRIGS
#   define Y_PRELOAD_TRIGS
#  endif
#  ifndef Y_PRELOAD_DATA
#   define Y_PRELOAD_DATA
#  endif
#  ifdef Y_LONG_MACROS
#   undef Y_LONG_MACROS
#  endif
#elif defined(Y_SSE2)
#  include "ysse2.h"
#  include "ygensse2.h"
#else
# ifdef Y_LONG_MACROS
#  include "ygeneric8.h"
# endif
#endif


/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* asign the correct values to every pad */

#define get_pads_8(_n0)                 \
  jj = pad;                             \
  pd0 = d;                              \
  prefetch_data( pd0, addr((_n0 << 1)));\
  pd1 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd1, addr((_n0 << 1)));\
  pd2 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd2, addr((_n0 << 1)));\
  pd3 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd3, addr((_n0 << 1)));\
  pd4 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd4, addr((_n0 << 1)));\
  pd5 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd5, addr((_n0 << 1)));\
  pd6 = d + addr((jj << 1));            \
  jj += pad;                            \
  prefetch_data( pd6, addr((_n0 << 1)));\
  pd7 = d + addr((jj << 1));            \
  prefetch_data( pd7, addr((_n0 << 1)));


#define get_pads_8_preload(_n0) \
  jj = pad + _n0;               \
  pd0 = d + addr(_n0 << 1);     \
  pd1 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd2 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd3 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd4 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd5 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd6 = d + addr((jj << 1));    \
  jj += pad;                    \
  pd7 = d + addr((jj << 1));

/* get the basic twiddle from memory and computes its powers */
#if defined(Y_MINIMUM)

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# include "yminimum.h"

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)              \
          tpw1r = *(_px );                   \
          tpw1i = *(_px + 1);

#  define trig_8_preload(_px)                \
          tpw1r = *(_px + Y_STEP);           \
          tpw1i = *(_px + Y_STEP + 1);

#  define get_twiddle_factors_8_preload      \
          tw1r = tpw1r; tw1i = tpw1i;        \
          trig_8_preload (px);               \
          cplx_trig_min_square (tw2, tw1);   \
          cplx_trig_min_square (tw4, tw2);   \
          cplx_divmul (tw3, tw5, tw4, tw1);  \
          cplx_trig_min_square (tw6, tw3);   \
          cplx_trig_min_mul (tw7, tw4, tw3);

#  define get_twiddle_factors_8_preload_consumed \
          t1r = tpw1r; t1i = tpw1i;              \
          trig_8_preload (px);                   \
          cplx_trig_min_square (t2, t1);         \
          cplx_trig_min_square (t4, t2);         \
          cplx_divmul (t3, t5, t4, t1);          \
          cplx_trig_min_square (t6, t3);         \
          cplx_trig_min_mul (t7, t4, t3);

# endif /* Y_PRELOAD_TRIGS */

# define get_twiddle_factors_8              \
          cplx_trig_min_load (tw1, px);     \
          cplx_trig_min_square (tw2, tw1);  \
          cplx_trig_min_square (tw4, tw2);  \
          cplx_divmul (tw3, tw5, tw4, tw1); \
          cplx_trig_min_square (tw6, tw3);  \
          cplx_trig_min_mul (tw7, tw4, tw3);\
          prefetch_p_trig (px);

#   define get_twiddle_factors_8_inloop  get_twiddle_factors_8
#   define get_twiddle_factors_8_outloop get_twiddle_factors_8

#elif defined(Y_MAXIMUM)

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 14

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)      \
          tpw1r = *(_px);            \
          tpw1i = *(_px + 1);        \
          tpw2r = *(_px + 2);        \
          tpw2i = *(_px + 3);        \
          tpw3r = *(_px + 4);        \
          tpw3i = *(_px + 5);        \
          tpw4r = *(_px + 6);        \
          tpw4i = *(_px + 7);        \
          tpw5r = *(_px + 8);        \
          tpw5i = *(_px + 9);        \
          tpw6r = *(_px + 10);       \
          tpw6i = *(_px + 11);       \
          tpw7r = *(_px + 12);       \
          tpw7i = *(_px + 13);       \


#  define trig_8_preload(_px)                 \
          tpw1r = *(_px + Y_STEP);            \
          tpw1i = *(_px + Y_STEP + 1);        \
          tpw2r = *(_px + Y_STEP + 2);        \
          tpw2i = *(_px + Y_STEP + 3);        \
          tpw3r = *(_px + Y_STEP + 4);        \
          tpw3i = *(_px + Y_STEP + 5);        \
          tpw4r = *(_px + Y_STEP + 6);        \
          tpw4i = *(_px + Y_STEP + 7);        \
          tpw5r = *(_px + Y_STEP + 8);        \
          tpw5i = *(_px + Y_STEP + 9);        \
          tpw6r = *(_px + Y_STEP + 10);       \
          tpw6i = *(_px + Y_STEP + 11);       \
          tpw7r = *(_px + Y_STEP + 12);       \
          tpw7i = *(_px + Y_STEP + 13);       \

#  define get_twiddle_factors_8_preload_consumed  \
   cplx_complete_twiddle_8_ia64 (t, tpw);         \
   trig_8_preload (px);

#  define get_twiddle_factors_8_preload           \
   cplx_complete_twiddle_8_ia64 (tw, tpw);        \
   trig_8_preload (px);

# endif /* Y_PRELOAD_TRIGS */

# if Y_TARGET == 0

#  define asim_prefetch_F_trig(_j)   /* */

#  define asim_prefetch_B_trig(_j)          \
          prefetch_data_trig_nta (px, 4 * Y_STEP + _j * Y_CACHE_LINE);

#  define get_twiddle_factors_8_inloop  /* */


#  define get_twiddle_factors_8_outloop                     \
          prefetch_data_trig (px, Y_STEP);                  \
          prefetch_data_trig (px, Y_STEP + Y_CACHE_LINE);   \
          prefetch_data_trig (px, Y_STEP + 2 * Y_CACHE_LINE);


# else /* Y_TARGET */
#  define get_twiddle_factors_8                              \
	  tw1r = px[0];  tw1i = px[1];                       \
          tw2r = px[2];  tw2i = px[3];                       \
          tw3r = px[4];  tw3i = px[5];                       \
          tw4r = px[6];  tw4i = px[7];                       \
          tw5r = px[8];  tw5i = px[9];                       \
          tw6r = px[10]; tw6i = px[11];                      \
          tw7r = px[12]; tw7i = px[13];                      \
          prefetch_data_trig (px, 2 * Y_STEP);               \
          prefetch_data_trig (px, 2 * Y_STEP + Y_CACHE_LINE);\
          prefetch_data_trig (px, 2 * Y_STEP + 2 * Y_CACHE_LINE);

#   define get_twiddle_factors_8_inloop  get_twiddle_factors_8
#   define get_twiddle_factors_8_outloop get_twiddle_factors_8


# endif


#else

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6


# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)      \
          tpw1r = *(_px);            \
          tpw1i = *(_px + 1);        \
          tpw2r = *(_px + 2);        \
          tpw2i = *(_px + 3);        \
          tpw5r = *(_px + 4);        \
          tpw5i = *(_px + 5);        \


#  define trig_8_preload(_px)                 \
          tpw1r = *(_px + Y_STEP);            \
          tpw1i = *(_px + Y_STEP + 1);        \
          tpw2r = *(_px + Y_STEP + 2);        \
          tpw2i = *(_px + Y_STEP + 3);        \
          tpw5r = *(_px + Y_STEP + 4);        \
          tpw5i = *(_px + Y_STEP + 5);        \

#  define get_twiddle_factors_8_preload_consumed  \
   cplx_complete_twiddle_8_ia64 (t, tpw);         \
   trig_8_preload (px);

#  define get_twiddle_factors_8_preload       \
   cplx_complete_twiddle_8_ia64 (tw, tpw);    \
   trig_8_preload (px);



# endif /* Y_PRELOAD_TRIGS */

# define get_twiddle_factors_8_inloop             \
	  tw1r = px[0];  tw1i = px[1];            \
          tw2r = px[2];  tw2i = px[3];            \
          tw5r = px[4];  tw5i = px[5];            \
	  cplx_divmul (tw3, tw7, tw5, tw2);       \
	  cplx_divmul (tw4, tw6, tw5, tw1);       \
          prefetch_data_trig_nta (px, 2 * Y_STEP);\
          prefetch_data_trig_nta (px, 2 * Y_STEP + Y_CACHE_LINE);

# define get_twiddle_factors_8_outloop            \
	  tw1r = px[0];  tw1i = px[1];            \
          tw2r = px[2];  tw2i = px[3];            \
          tw5r = px[4];  tw5i = px[5];            \
	  cplx_divmul (tw3, tw7, tw5, tw2);       \
	  cplx_divmul (tw4, tw6, tw5, tw1);       \
          prefetch_data_trig (px, Y_STEP);        \
          prefetch_data_trig (px, Y_STEP + Y_CACHE_LINE);


#endif /* Y_MAXIMUM || Y_MINIMUM */


/*
	Get data from memory, mul by twiddle trig factors,
	make a complex 8-length FFT and store the data to memory 
*/

#ifdef Y_ITANIUM

# define radix_8_twd_preload_F(_k)                                                 \
    t0r = s0r; t0i = s0i; ppd0r = ppdl0r; ppd0i = ppdl0i;                          \
    ppdl0r +=_k; ppdl0i += _k; s0r = *(ppdl0r); s0i = *(ppdl0i);                   \
    cplx_three_muls_fftn_preload (t2, t5, t7, 2, 5, 7, s, tw, ppd, ppdl, _k);      \
    cplx_four_muls_fftn_preload (t1, t4, t3, t6, 4, 1, 6, 3, s, tw, ppd, ppdl, _k);\
    cplx_fft8_ia64_store_F (ppd0, ppd1, ppd2, ppd3, ppd4, ppd5, ppd6, ppd7, t, 0, 1, 2, 3, 4, 5, 6, 7);

# define radix_8_twd_preload_B(_k)                                                                      \
    t0r = s0r; t0i = s0i; ppd0r = ppdl0r; ppd0i = ppdl0i;                                               \
    ppdl0r +=_k; ppdl0i += _k; s0r = *(ppdl0r); s0i =*(ppdl0i);                                         \
    cplx_three_muls_fftn_preload_consumed_trigs (t2, t5, t7, 2, 5, 7, s, t, ppd, ppdl, _k);             \
    cplx_four_muls_fftn_preload_consumed_trigs_paired (t1, t4, t3, t6, 1, 4, 3, 6, s, t, ppd, ppdl, _k);\
    cplx_fft8_ia64_store_B (ppd0, ppd1, ppd2, ppd3, ppd4, ppd5, ppd6, ppd7, t, 0, 1, 2, 3, 4, 5, 6, 7);


# define radix_8_notwd_preload_F(_inc) \
    cplx_fft8_preload_store_F (ppd0, ppd1, ppd2, ppd3, ppd4, ppd5, ppd6, ppd7, t, s, 0, 4, 2, 6, 1, 5, 3, 7, ppd, ppdl, _inc);

# define radix_8_notwd_preload_B(_inc) \
    cplx_fft8_preload_store_B (ppd0, ppd1, ppd2, ppd3, ppd4, ppd5, ppd6, ppd7, t, s, 0, 4, 2, 6, 1, 5, 3, 7, ppd, ppdl, _inc);

#endif

# ifdef Y_LONG_MACROS

#  define radix_8_twd_pp(_S,_j)                    \
          cplx_radix_8_twd_pp_pass1(t,tw,_j,pd);   \
          cplx_radix_8_pp_pass2_##_S(t);           \
          cplx_radix_8_pp_pass3_store_##_S(t,_j,pd);

#  define radix_8_twd_pp_no_fetch(_S,_j)                 \
          cplx_radix_8_twd_pp_pass1_no_fetch(t,tw,_j,pd);\
          cplx_radix_8_pp_pass2_##_S(t);                 \
          cplx_radix_8_pp_pass3_store_##_S(t,_j,pd);

#  define radix_8_notwd_pp(_S,_j)                  \
          cplx_radix_8_notwd_pp_pass1(t,_j,pd);    \
          cplx_radix_8_pp_pass2_##_S(t);           \
          cplx_radix_8_pp_pass3_store_##_S(t,_j,pd);

#  define radix_8_notwd_pp_no_fetch(_S,_j)              \
          cplx_radix_8_notwd_pp_pass1_no_fetch(t,_j,pd);\
          cplx_radix_8_pp_pass2_##_S(t);                \
          cplx_radix_8_pp_pass3_store_##_S(t,_j,pd);

# else

#  if defined(Y_MAXIMUM) && (Y_TARGET == 0)

#  define radix_8_twd_pp(_S,_j)                                \
	  cplx_load_muladdsub_pp_ind(t0,t1,px,4,_j,pd0, pd4);   \
	  cplx_load_mulmuladdsub_pp_ind(t2,t3,px,2,6,_j,pd2,pd6); \
          asim_prefetch_##_S##_trig(0);                        \
	  cplx_addsub(t0,t2);                                  \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                   \
	  cplx_load_mulmuladdsub_pp_ind(t4,t5,px,1,5,_j,pd1,pd5); \
	  cplx_load_mulmuladdsub_pp_ind(t6,t7,px,3,7,_j,pd3,pd7); \
          asim_prefetch_##_S##_trig(1);                        \
	  cplx_addsub(t4,t6);                                  \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                   \
                                                               \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);\
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);\
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

#  else

#  define radix_8_twd_pp(_S,_j)                                \
	  cplx_load_muladdsub_pp(t0,t1,tw4,_j,pd0, pd4);       \
	  cplx_load_mulmuladdsub_pp(t2,t3,tw2,tw6,_j,pd2,pd6); \
	  cplx_addsub(t0,t2);                                  \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                   \
	  cplx_load_mulmuladdsub_pp(t4,t5,tw1,tw5,_j,pd1,pd5); \
	  cplx_load_mulmuladdsub_pp(t6,t7,tw3,tw7,_j,pd3,pd7); \
	  cplx_addsub(t4,t6);                                  \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                   \
                                                               \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);\
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);\
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

# endif

#  define radix_8_twd_pp_no_fetch(_S,_j)                               \
	  cplx_load_muladdsub_pp_no_fetch(t0,t1,tw4,_j,pd0, pd4);      \
	  cplx_load_mulmuladdsub_pp_no_fetch(t2,t3,tw2,tw6,_j,pd2,pd6);\
	  cplx_addsub(t0,t2);                                          \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                           \
	  cplx_load_mulmuladdsub_pp_no_fetch(t4,t5,tw1,tw5,_j,pd1,pd5);\
	  cplx_load_mulmuladdsub_pp_no_fetch(t6,t7,tw3,tw7,_j,pd3,pd7);\
	  cplx_addsub(t4,t6);                                          \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                           \
                                                                       \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);                       \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);        \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);        \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

#  define radix_8_twd_preload(_S,_j,_jj)                               \
	  cplx_preload_muladdsub(t0,t1,tw4,_jj,pd0,pd4,tp0,tp4);       \
	  cplx_preload_mulmuladdsub(t2,t3,tw2,tw6,_jj,pd2,pd6,tp2,tp6);\
	  cplx_preload_mulmuladdsub(t4,t5,tw1,tw5,_jj,pd1,pd5,tp1,tp5);\
	  cplx_preload_mulmuladdsub(t6,t7,tw3,tw7,_jj,pd3,pd7,tp3,tp7);\
	  cplx_addsub(t0,t2);                                          \
          cplx_mul_1_4_##_S##_addsub(t1,t3);                           \
	  cplx_addsub(t4,t6);                                          \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                           \
                                                                       \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);                       \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);        \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);        \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);        \
          _j = _jj;

#  define radix_8_notwd_pp(_S,_j)                              \
	  cplx_load_addsub_pp(t0,t1,_j,pd0,pd4);               \
	  cplx_load_addsub_pp(t2,t3,_j,pd2,pd6);               \
	  cplx_addsub(t0,t2);                                  \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                   \
	  cplx_load_addsub_pp(t4,t5,_j,pd1,pd5);               \
	  cplx_load_addsub_pp(t6,t7,_j,pd3,pd7);               \
	  cplx_addsub(t4,t6);                                  \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                   \
                                                               \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);\
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);\
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

#  define radix_8_notwd_pp_no_fetch(_S,_j)                     \
	  cplx_load_addsub_pp_no_fetch(t0,t1,_j,pd0,pd4);      \
	  cplx_load_addsub_pp_no_fetch(t2,t3,_j,pd2,pd6);      \
	  cplx_addsub(t0,t2);                                  \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                   \
	  cplx_load_addsub_pp_no_fetch(t4,t5,_j,pd1,pd5);      \
	  cplx_load_addsub_pp_no_fetch(t6,t7,_j,pd3,pd7);      \
	  cplx_addsub(t4,t6);                                  \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                   \
                                                               \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);               \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5);\
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6);\
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7);

#  define radix_8_notwd_preload(_S,_j,_jj)                      \
	  cplx_preload_addsub(t0,t1,_jj,pd0,pd4,tp0,tp4);       \
	  cplx_preload_addsub(t2,t3,_jj,pd2,pd6,tp2,tp6);       \
	  cplx_preload_addsub(t4,t5,_jj,pd1,pd5,tp1,tp5);       \
	  cplx_preload_addsub(t6,t7,_jj,pd3,pd7,tp3,tp7);       \
	  cplx_addsub(t0,t2);                                   \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                    \
	  cplx_addsub(t4,t6);                                   \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                    \
                                                                \
	  cplx_addsub_store_p(_j,pd0,pd4,t0,t4);                \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd1,pd5,t1,t5); \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd2,pd6,t2,t6); \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd3,pd7,t3,t7); \
          _j = _jj;

# endif
/*#endif*/

#define RADIX 8

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
void radixmp_8_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
#if !defined(Y_MAXIMUM) || (Y_TARGET != 0)

  y_limb_t  tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
#endif
#ifdef Y_ITANIUM

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t fr=F_1_8r,fi=F_1_8i;
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
# elif defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,
  tpw6r,tpw6i,tpw7r,tpw7i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw5r,tpw5i;
# endif
#endif

  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
#ifdef Y_ITANIUM

  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i,ppd4r,ppd4i,
  ppd5r,ppd5i,ppd6r,ppd6i,ppd7r,ppd7i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i,ppdl4r,ppdl4i,
  ppdl5r,ppdl5i,ppdl6r,ppdl6i,ppdl7r,ppdl7i;
#endif

  y_size_t i,j,jj,nc;
#ifdef Y_PRELOAD_DATA

  y_size_t next,last,nexti;
#endif

  y_size_t done, mask, ioff, ibase, n00;


  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmp_8_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  bigpad = (pad << 3);

  nc = n + n0;
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad);
  ioff = (n0 & mask) >> 3;
  n00 = (ibase * bigpad + ioff);

#ifndef NDEBUG1

  printf(" ibase=%d ioff=%d n00=%d\n", ibase, ioff, n00);
#endif

#ifdef Y_PRELOAD_DATA

  get_pads_8_preload (n00);
#else

  get_pads_8 (n00);
#endif

#ifdef Y_PRELOAD_DATA

  cplx_radix_8_dif_dit_init(s,pd,ppdl);
  last = addr(n00 << 1);
#endif

  /* If ibase == 0 first elements has twiddle factors equal to one */
  if(ibase == 0)
    {
#ifdef Y_PRELOAD_DATA
      if (bigpad < nc )
        nexti = bigpad;
      else
        nexti = pad - 1;
#endif

      for(j = ioff; j < pad && done < n; j++)
        {
#ifdef Y_PRELOAD_DATA
          if((j + 1) < pad)
            next = j + 1;
          else
            next = nexti;
          next = addr(next<<1);
          jj= next - last;
          last = next;
          radix_8_notwd_preload_F(jj);
          done += RADIX;
#else

          jj = addr(j<<1);
          radix_8_notwd_pp(F,jj);
          done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

          j++;
          jj += 2;
          radix_8_notwd_pp_no_fetch(F,jj);
          done += RADIX;
#  if Y_CACHE_LINE == 8

          j += 2;
          jj += 2;
          radix_8_notwd_pp_no_fetch(F,jj);
          done += RADIX;
          jj+=2;
          radix_8_notwd_pp_no_fetch(F,jj);
          done += RADIX;
#  endif
# endif
#endif

        }
      /*if (done == n) printf("F\n");*/
      if( done == n)
        return;
      px = tw + Y_STEP;
#ifdef Y_PRELOAD_TRIGS

      trig_8_load_init(px);
#endif

      for(i = bigpad; done < n ;i += bigpad, px += Y_STEP)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_8_preload;
#else

          get_twiddle_factors_8_outloop;
#endif
#ifdef Y_PRELOAD_DATA

          if((i + bigpad) < nc)
            nexti = i + bigpad;
          else
            nexti= i + pad - 1 ;
          for(j = i; j < (pad + i) && done < n; j++)
            {
              if((j + 1) < (pad + i))
                next = j + 1;
              else
                next = nexti;
              next = addr(next<<1);
              jj= next - last;
              last = next;
              radix_8_twd_preload_F(jj);
              done += RADIX;
            }
#else
          for(j = i; j < (pad + i) && done < n; j++)
            {
              jj = addr(j<<1);
              radix_8_twd_pp(F,jj);
              done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  if Y_CACHE_LINE == 8

              j += 2;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  endif
# endif

            }
#endif

        }
    }
  else
    {
      px = tw + Y_STEP * ibase;
#ifdef Y_PRELOAD_TRIGS

      trig_8_load_init(px);
#endif

      for(i = ibase * bigpad; done < n; i += bigpad, px += Y_STEP)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_8_preload;
#else
          /*BEGIN_TIME;*/
          get_twiddle_factors_8_outloop;
          /*END_TIME;*/
#endif
#ifdef Y_PRELOAD_DATA

          if((i + bigpad) < nc)
            nexti = i + bigpad;
          else
            nexti= i + pad - 1;
          for(j = i + ioff; j < (pad + i) && done < n ; j++)
            {
              if((j + 1) < (pad + i))
                next = j + 1;
              else
                next = nexti;
              next = addr(next<<1);
              jj= next - last;
              last = next;
              radix_8_twd_preload_F(jj);
              done += RADIX;
            }
          ioff = 0;
#else

          for(j = i + ioff ; j < (pad + i) && done < n; j++)
            {
              jj = addr(j<<1);
              radix_8_twd_pp(F,jj);
              done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  if Y_CACHE_LINE == 8

              j += 2;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
              jj += 2;
              radix_8_twd_pp_no_fetch(F,jj);
              done += RADIX;
#  endif
# endif

            }
          ioff = 0;
#endif

        }
    }
  /*if (done == n) printf("F\n");*/
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

void radixmp_8_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
#ifdef Y_PRELOAD_DATA

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t fr=F_1_8r,fi=F_1_8i;
#elif !defined(Y_MAXIMUM) || (Y_TARGET != 0)

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i;
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
# elif defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,tpw6r,tpw6i,
  tpw7r,tpw7i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw5r,tpw5i;
# endif
#endif

  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7;
#ifdef Y_ITANIUM

  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i,ppd4r,ppd4i,
  ppd5r,ppd5i,ppd6r,ppd6i,ppd7r,ppd7i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i,ppdl4r,ppdl4i,
  ppdl5r,ppdl5i,ppdl6r,ppdl6i,ppdl7r,ppdl7i;
#endif

  y_size_t i,j,jj,nc;
#ifdef Y_SAVE_TRIGS
  /*y_size_t i1,j1,jj1;*/
#endif
#ifdef Y_PRELOAD_DATA

  y_size_t next, last, nexti;
#endif

  y_size_t done, mask, ioff, ibase, n00;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radixmp_8_dit n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  bigpad = (pad<<3);
  nc = n0 + n;
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad) * bigpad;
  ioff = (n0 & mask) >> 3;
  n00 = ibase + ioff;

#ifdef Y_PRELOAD_DATA

  get_pads_8_preload(n00);
#else

  get_pads_8(n00);
#endif


#ifdef Y_PRELOAD_DATA

  cplx_radix_8_dif_dit_init(s,pd,ppdl);
  last = addr(n00 << 1);
#endif

  for(i = ibase; done < n; i += bigpad)
    {
#ifdef Y_ITANIUM
      if((i + bigpad) < nc)
        nexti = i + bigpad;
      /*else nexti= i + pad - 1;*/
      else
        nexti= ((nc - i) >> 3) + i;
#endif

      if (ioff == 0)
        {
#ifdef Y_ITANIUM
          next = addr((i+1)<<1);
          jj= next - last;
          last = next;
          radix_8_notwd_preload_B(jj);
          done += RADIX;
          trig_8_load_init(tw);
          px = tw;
          j = i + 1;
          /*if((i + bigpad) < nc) nexti = i + bigpad;
            else nexti= i + pad - 1;*/
#else

          jj = addr(i<<1);
          radix_8_notwd_pp(B,jj);
          done += RADIX;
          px = tw;
          j = i + 1;
# if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4

          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j++;
#  if Y_CACHE_LINE == 8

          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j += 2;
#  endif
# endif
#endif

        }
      else
        {
          px = tw + (ioff - 1) * Y_STEP;
          j = i + ioff;
#ifdef Y_ITANIUM

          trig_8_load_init(px);
#endif

        }

      for(;j < (pad + i) && done < n; j++, px += Y_STEP)
        {
#ifdef Y_ITANIUM
          get_twiddle_factors_8_preload_consumed;
          if((j + 1) < (pad + i))
            next = j + 1;
          else
            next = nexti;
          next = addr(next << 1);
          jj= next - last;
          last = next;
          radix_8_twd_preload_B(jj);
          done += RADIX;
#else

          get_twiddle_factors_8_inloop;
          jj = addr(j << 1);
          radix_8_twd_pp(B,jj);
          done += RADIX;
# if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4

          px += Y_STEP;
          j++;
          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
#  if Y_CACHE_LINE == 8

          px+=Y_STEP;
          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
          px += Y_STEP;
          j += 2;
          get_twiddle_factors_8;
          jj += 2;
          radix_8_twd_pp_no_fetch(B,jj);
          done += RADIX;
#  endif
# endif
#endif

        }
      ioff = 0;
    }
  /*if (done == n) printf("B\n");*/
}

