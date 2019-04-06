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
/*
   This the subroutines are related with dyadic mul and square 
   It includes routines for a 8_last_radix reduction 
*/
#include <stdio.h>
#include <stdlib.h>

#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1
#ifdef Y_ITANIUM
# include "yia64.h"
# ifndef Y_PRELOAD_TRIGS
#  define Y_PRELOAD_TRIGS
# endif
# ifndef Y_PRELOAD_DATA
#  define Y_PRELOAD_DATA
# endif
# ifdef Y_LONG_MACROS
#  undef Y_LONG_MACROS
# endif
#else
# ifdef Y_LONG_MACROS
#  include "ygeneric8.h"
# endif
#endif

#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# include "yminimum.h"

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_tpw,_px,_j)               \
          _px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]); \
          _tpw##1r = *(_px );                         \
          _tpw##1i = *(_px + 1);

#  define trig_8_preload(_tpw,_px)              \
          _tpw##1r = *(_px + Y_STEP);           \
          _tpw##1i = *(_px + Y_STEP + 1);

#  define trig_8_preload_down(_tpw,_px)         \
          _tpw##1r = *(_px - Y_STEP);           \
          _tpw##1i = *(_px - Y_STEP + 1);

#  define get_twiddle_factors_8_preload(_px,_tw,_tpw)  \
          _tw##1r = _tpw##1r; _tw##1i = _tpw##1i;      \
          trig_8_preload(_tpw,_px);                    \
          cplx_trig_min_square(_tw##2,_tw##1);         \
          cplx_trig_min_square(_tw##4,_tw##2);         \
          cplx_divmul(_tw##3,_tw##5,_tw##4,_tw##1);    \
          cplx_trig_min_square(_tw##6,_tw##3);         \
          cplx_trig_min_mul(_tw##7,_tw##4,_tw##3);

#  define get_twiddle_factors_8_preload_down(_px,_tw,_tpw)  \
          _tw##1r = _tpw##1r; _tw##1i = _tpw##1i;           \
          trig_8_preload_down(_tpw,_px);                    \
          cplx_trig_min_square(_tw##2,_tw##1);              \
          cplx_trig_min_square(_tw##4,_tw##2);              \
          cplx_divmul(_tw##3,_tw##5,_tw##4,_tw##1);         \
          cplx_trig_min_square(_tw##6,_tw##3);              \
          cplx_trig_min_mul(_tw##7,_tw##4,_tw##3);

# endif

# define get_twiddle_factors_8_last(_j)            \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          cplx_trig_min_load(tw1,px);               \
          cplx_trig_min_square(tw2,tw1);            \
          cplx_trig_min_square(tw4,tw2);            \
          cplx_divmul(tw3,tw5,tw4,tw1);             \
          cplx_trig_min_square(tw6,tw3);            \
          cplx_trig_min_mul(tw7,tw4,tw3);           \
          prefetch_p_trig_nta(px);

# define get_twiddle_factors_8_last_down(_j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          cplx_trig_min_load(tw1,px);               \
          cplx_trig_min_square(tw2,tw1);            \
          cplx_trig_min_square(tw4,tw2);            \
          cplx_divmul(tw3,tw5,tw4,tw1);             \
          cplx_trig_min_square(tw6,tw3);            \
          cplx_trig_min_mul(tw7,tw4,tw3);           \
          prefetch_data_trig_nta(px,(-2));

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 14

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_tpw,_px,_j)                 \
          _px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);   \
          _tpw##1r = *(_px);                            \
          _tpw##1i = *(_px + 1);                        \
          _tpw##2r = *(_px + 2);                        \
          _tpw##2i = *(_px + 3);                        \
          _tpw##3r = *(_px + 4);                        \
          _tpw##3i = *(_px + 5);                        \
          _tpw##4r = *(_px + 6);                        \
          _tpw##4i = *(_px + 7);                        \
          _tpw##5r = *(_px + 8);                        \
          _tpw##5i = *(_px + 9);                        \
          _tpw##6r = *(_px + 10);                       \
          _tpw##6i = *(_px + 11);                       \
          _tpw##7r = *(_px + 12);                       \
          _tpw##7i = *(_px + 13);

#  define trig_8_preload(_tpw,_px)               \
          _tpw##1r = *(_px + Y_STEP);            \
          _tpw##1i = *(_px + Y_STEP + 1);        \
          _tpw##2r = *(_px + Y_STEP + 2);        \
          _tpw##2i = *(_px + Y_STEP + 3);        \
          _tpw##3r = *(_px + Y_STEP + 4);        \
          _tpw##3i = *(_px + Y_STEP + 5);        \
          _tpw##4r = *(_px + Y_STEP + 6);        \
          _tpw##4i = *(_px + Y_STEP + 7);        \
          _tpw##5r = *(_px + Y_STEP + 8);        \
          _tpw##5i = *(_px + Y_STEP + 9);        \
          _tpw##6r = *(_px + Y_STEP + 10);       \
          _tpw##6i = *(_px + Y_STEP + 11);       \
          _tpw##7r = *(_px + Y_STEP + 12);       \
          _tpw##7i = *(_px + Y_STEP + 13);

#  define trig_8_preload_down(_tpw,_px)          \
          _tpw##1r = *(_px - Y_STEP);            \
          _tpw##1i = *(_px - Y_STEP + 1);        \
          _tpw##2r = *(_px - Y_STEP + 2);        \
          _tpw##2i = *(_px - Y_STEP + 3);        \
          _tpw##3r = *(_px - Y_STEP + 4);        \
          _tpw##3i = *(_px - Y_STEP + 5);        \
          _tpw##4r = *(_px - Y_STEP + 6);        \
          _tpw##4i = *(_px - Y_STEP + 7);        \
          _tpw##5r = *(_px - Y_STEP + 8);        \
          _tpw##5i = *(_px - Y_STEP + 9);        \
          _tpw##6r = *(_px - Y_STEP + 10);       \
          _tpw##6i = *(_px - Y_STEP + 11);       \
          _tpw##7r = *(_px - Y_STEP + 12);       \
          _tpw##7i = *(_px - Y_STEP + 13);

#  define get_twiddle_factors_8_preload(_px,_tw,_tpw)     \
   cplx_complete_twiddle_8_ia64(_tw,_tpw);                \
   trig_8_preload(_tpw,_px);

#  define get_twiddle_factors_8_preload_down(_px,_tw,_tpw)\
   cplx_complete_twiddle_8_ia64(_tw,_tpw);                \
   trig_8_preload_down(_tpw,_px);

# endif

# if defined(Y_MAXIMUM) && (Y_TARGET == 0)

#  define get_twiddle_factors_8_last(_j)             \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          prefetch_data_trig_nta( px, 2*Y_STEP + Y_CACHE_LINE);\
          prefetch_data_trig_nta( px, 2*Y_STEP + 2*Y_CACHE_LINE);

#  define get_twiddle_factors_8_last_down(_j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          prefetch_data_trig_nta( px, Y_CACHE_LINE - 2*Y_STEP);\
          prefetch_data_trig_nta( px, Y_CACHE_LINE - 2*Y_STEP);

# else

#  define get_twiddle_factors_8_last(_j)             \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw1r=px[0];   tw1i=px[1];                 \
          tw2r=px[2];   tw2i=px[3];                 \
          tw3r=px[4];   tw3i=px[5];                 \
          tw4r=px[6];   tw4i=px[7];                 \
          tw5r=px[8];   tw5i=px[9];                 \
          tw6r=px[10];   tw6i=px[11];               \
          tw7r=px[12];   tw7i=px[13];               \
          prefetch_data_trig_nta( px, 14 + Y_CACHE_LINE);\
          prefetch_data_trig_nta( px, 14 + 2*Y_CACHE_LINE);

#  define get_twiddle_factors_8_last_down(_j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw1r=px[0];   tw1i=px[1];                 \
          tw2r=px[2];   tw2i=px[3];                 \
          tw3r=px[4];   tw3i=px[5];                 \
          tw4r=px[6];   tw4i=px[7];                 \
          tw5r=px[8];   tw5i=px[9];                 \
          tw6r=px[10];   tw6i=px[11];               \
          tw7r=px[12];   tw7i=px[13];               \
          prefetch_data_trig_nta( px, (-14) + Y_CACHE_LINE);\
          prefetch_data_trig_nta( px, (-14) + 2*Y_CACHE_LINE);

# endif

#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_tpw,_px,_j)                 \
          _px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);   \
          _tpw##1r = *(_px);                            \
          _tpw##1i = *(_px + 1);                        \
          _tpw##2r = *(_px + 2);                        \
          _tpw##2i = *(_px + 3);                        \
          _tpw##5r = *(_px + 4);                        \
          _tpw##5i = *(_px + 5);                        \


#  define trig_8_preload(_tpw,_px)               \
          _tpw##1r = *(_px + Y_STEP);            \
          _tpw##1i = *(_px + Y_STEP + 1);        \
          _tpw##2r = *(_px + Y_STEP + 2);        \
          _tpw##2i = *(_px + Y_STEP + 3);        \
          _tpw##5r = *(_px + Y_STEP + 4);        \
          _tpw##5i = *(_px + Y_STEP + 5);        \

#  define trig_8_preload_down(_tpw,_px)          \
          _tpw##1r = *(_px - Y_STEP);            \
          _tpw##1i = *(_px - Y_STEP + 1);        \
          _tpw##2r = *(_px - Y_STEP + 2);        \
          _tpw##2i = *(_px - Y_STEP + 3);        \
          _tpw##5r = *(_px - Y_STEP + 4);        \
          _tpw##5i = *(_px - Y_STEP + 5);        \

#  define get_twiddle_factors_8_preload(_px,_tw,_tpw)     \
   cplx_complete_twiddle_8_ia64(_tw,_tpw);                \
   trig_8_preload(_tpw,_px);

#  define get_twiddle_factors_8_preload_down(_px,_tw,_tpw)\
   cplx_complete_twiddle_8_ia64(_tw,_tpw);                \
   trig_8_preload_down(_tpw,_px);

# endif

# define get_twiddle_factors_8_last(_j)            \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw1r=px[0];   tw1i=px[1];                 \
          tw2r=px[2];   tw2i=px[3];                 \
          tw5r=px[4];   tw5i=px[5];                 \
          cplx_divmul(tw3,tw7,tw5,tw2);             \
          cplx_divmul(tw4,tw6,tw5,tw1);             \
          prefetch_data_trig_nta( px, 6);           \
          prefetch_data_trig_nta( px, 10);

# define get_twiddle_factors_8_last_down(_j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw1r=px[0];	  tw1i=px[1];               \
          tw2r=px[2];     tw2i=px[3];               \
          tw5r=px[4];     tw5i=px[5];               \
          cplx_divmul(tw3,tw7,tw5,tw2);             \
          cplx_divmul(tw4,tw6,tw5,tw1);             \
          prefetch_data_trig_nta( px, (-6));        \
          prefetch_data_trig_nta( px, (-2));


#endif

#ifdef Y_PRELOAD_DATA

# define radix_8_twd_last_preload(_t,_s,_ppd,_ppdl,_k)                     \
    _t##0r = _s##0r; _t##0i = _s##0i;                                      \
    _ppd##0r = _ppdl##0r; _ppd##0i = _ppdl##0i;                            \
    _ppdl##0r += _k; _ppdl##0i += _k;                                       \
    _s##0r= *(_ppdl##0r); _s##0i= *( _ppdl##0i);                           \
    cplx_three_muls_fftn_preload_consumed_trigs(_t##2,_t##5,_t##7,2,5,7,_s,_t,_ppd,_ppdl,_k);\
    cplx_four_muls_fftn_preload_consumed_trigs_paired(_t##1,_t##4,_t##3,_t##6,1,4,3,6,_s,_t,_ppd,_ppdl,_k);\
    cplx_fft8_ia64_F(_t,0,1,2,3,4,5,6,7);

# define radix_8_twd_last_preload_down(_t,_s,_ppd,_ppdl,_k)                \
    _t##0r = _s##0r; _t##0i = _s##0i;                                      \
    _ppd##0r = _ppdl##0r; _ppd##0i = _ppdl##0i;                            \
    _ppdl##0r -= _k; _ppdl##0i -= _k;                                       \
    _s##0r= *(_ppdl##0r); _s##0i= *( _ppdl##0i);                           \
    cplx_three_muls_fftn_preload_consumed_trigs_down(_t##2,_t##5,_t##7,2,5,7,_s,_t,_ppd,_ppdl,_k);\
    cplx_four_muls_fftn_preload_consumed_trigs_paired_down(_t##1,_t##4,_t##3,_t##6,1,4,3,6,_s,_t,_ppd,_ppdl,_k);\
    cplx_fft8_ia64_F(_t,0,1,2,3,4,5,6,7);


# define radix_8_notwd_first_preload(_p,_t) \
    cplx_fft8_ia64_store_B(_p##0,_p##1,_p##2,_p##3,_p##4,_p##5,_p##6,_p##7,_t,0,4,2,6,1,5,3,7);

#endif

#if defined(Y_LONG_MACROS)

#define radix_8_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
          cplx_radix_8_twd_pp_last_dif_pass1(_z,tw,_pd);\
          cplx_radix_8_pp_pass2_F(_z);\
          cplx_radix_8_pp_last_dif_pass3(_z);\

#define radix_8_twd_last_down(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
          cplx_radix_8_twd_pp_last_dif_pass1_down(_z,tw,_pd);\
          cplx_radix_8_pp_pass2_F(_z);\
          cplx_radix_8_pp_last_dif_pass3(_z);\

# define radix_8_notwd_first(_pd,_z,_j)\
          cplx_radix_8_pp_first_dit_pass1(_z);\
          cplx_radix_8_pp_first_dit_pass2(_z);\
          cplx_radix_8_pp_first_dit_pass3(_z,_pd);\

#else

# if defined(Y_MAXIMUM) && (Y_TARGET == 0)

# define radix_8_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
	  cplx_load_muladdsub_pp_ind( _z##0 , _z##1 , px, 4, 0, _pd, _pd + 8);\
          cplx_load_mulmuladdsub_pp_ind( _z##2,_z##3, px, 2, 6, 0,_pd + 4,_pd + 12);\
	  cplx_addsub( _z##0 ,_z##2 );\
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3);\
          cplx_load_mulmuladdsub_pp_ind( _z##4,_z##5, px, 1, 5, 0,_pd + 2,_pd + 10);\
          cplx_load_mulmuladdsub_pp_ind( _z##6,_z##7, px, 3, 7, 0,_pd + 6,_pd + 14);\
          prefetch_p( _pd + 8);\
          prefetch_p( _pd + 16);\
	  cplx_addsub( _z##4 ,_z##6 );\
	  cplx_mul_1_4_F_addsub( _z##5 , _z##7);\
	  \
	  cplx_addsub( _z##0, _z##4 );\
	  cplx_mul_1_8_F_addsub(_z##1 , _z##5 );\
	  cplx_mul_1_4_F_addsub(_z##2 , _z##6 );\
	  cplx_mul_3_8_F_addsub(_z##3 , _z##7 );\



# else


# define radix_8_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
	  cplx_load_muladdsub_pp( _z##0 , _z##1 , tw4, 0, _pd, _pd + 8);\
          cplx_load_mulmuladdsub_pp( _z##2,_z##3, tw2, tw6, 0,_pd + 4,_pd + 12);\
	  cplx_addsub( _z##0 ,_z##2 );\
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3);\
          cplx_load_mulmuladdsub_pp( _z##4,_z##5, tw1, tw5, 0,_pd + 2,_pd + 10);\
          cplx_load_mulmuladdsub_pp( _z##6,_z##7, tw3, tw7, 0,_pd + 6,_pd + 14);\
          prefetch_p( _pd + 8);\
          prefetch_p( _pd + 16);\
	  cplx_addsub( _z##4 ,_z##6 );\
	  cplx_mul_1_4_F_addsub( _z##5 , _z##7);\
	  \
	  cplx_addsub( _z##0, _z##4 );\
	  cplx_mul_1_8_F_addsub(_z##1 , _z##5 );\
	  cplx_mul_1_4_F_addsub(_z##2 , _z##6 );\
	  cplx_mul_3_8_F_addsub(_z##3 , _z##7 );\

# endif
# define radix_8_notwd_first(_pd,_z,_j)\
	  cplx_addsub( _z##0 , _z##4 );\
	  cplx_addsub( _z##2 , _z##6 );\
	  cplx_addsub( _z##0 , _z##2 );\
	  cplx_mul_1_4_B_addsub( _z##4 , _z##6 );\
	  cplx_addsub( _z##1 , _z##5 );\
	  cplx_addsub( _z##3 , _z##7 );\
	  cplx_addsub( _z##1 , _z##3 );\
	  cplx_mul_1_4_B_addsub( _z##5 , _z##7 );\
	  \
	  cplx_addsub_store_p( 0, _pd, _pd + 8, _z##0 , _z##1 );\
	  cplx_mul_1_8_B_addsub_store_p(0,_pd + 2,_pd + 10, _z##4 , _z##5 );\
	  cplx_mul_1_4_B_addsub_store_p(0,_pd + 4,_pd + 12, _z##2 , _z##3 );\
	  cplx_mul_3_8_B_addsub_store_p(0,_pd + 6,_pd + 14, _z##6 , _z##7 );


#endif

/*
#undef Y_PRELOAD_TRIGS
#undef Y_PRELOAD_DATA
#undef Y_ITANIUM
*/
void radix_8_dif_square_dit(y_ptr d, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i;
  y_limb_t fr=F_1_8r, fi=F_1_8i;
#ifdef Y_PRELOAD_DATA

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i;
#else

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i;
# if Y_PRELOAD_TRIGS

  y_limb_t uw1r,uw1i,uw2r,uw2i,uw3r,uw3i,uw4r,uw4i,uw5r,uw5i,uw6r,uw6i,
  uw7r,uw7i;
# endif
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
  y_limb_t upw1r,upw1i;
# elif defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,tpw6r,tpw6i,
  tpw7r,tpw7i;
  y_limb_t upw1r,upw1i,upw2r,upw2i,upw3r,upw3i,upw4r,upw4i,upw5r,upw5i,upw6r,upw6i,
  upw7r,upw7i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw5r,tpw5i;
  y_limb_t upw1r,upw1i,upw2r,upw2i,upw5r,upw5i;
# endif
#endif

  y_limb_t wkr,wki,aux;
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[7];
  int nr;
  y_ptr px;
#ifdef Y_PRELOAD_TRIGS

  y_ptr py;
#else

  y_ptr pt,pu;
#endif
#ifdef Y_PRELOAD_DATA

  y_ptr ppu0r,ppu0i,ppu1r,ppu1i,ppu2r,ppu2i,ppu3r,ppu3i,ppu4r,ppu4i,
  ppu5r,ppu5i,ppu6r,ppu6i,ppu7r,ppu7i;
  y_ptr ppul0r,ppul0i,ppul1r,ppul1i,ppul2r,ppul2i,ppul3r,ppul3i,ppul4r,ppul4i,
  ppul5r,ppul5i,ppul6r,ppul6i,ppul7r,ppul7i;
  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i,ppd4r,ppd4i,
  ppd5r,ppd5i,ppd6r,ppd6i,ppd7r,ppd7i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i,ppdl4r,ppdl4i,
  ppdl5r,ppdl5i,ppdl6r,ppdl6i,ppdl7r,ppdl7i;
#endif

  ASSERT_ALIGNED_DOUBLE();


  /*if(Y_PLAN[Y_NRADICES-1] != 8) return;*/
  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 8, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 8 
  */

  /* Cycle 0 block 0 */
  j=0;
#ifdef Y_PRELOAD_TRIGS

  trig_8_load_init(tpw,px,j);
  get_twiddle_factors_8_preload(px,t,tpw);
#else

  get_twiddle_factors_8_last(j);
#endif
#ifdef Y_PRELOAD_DATA

  cplx_radix_8_pad1_init(s,d,j,ppul);
  radix_8_twd_last_preload(t,s,ppu,ppul,16);
#else

  radix_8_twd_last(t,pt,d,j);
#endif
  /* dyadic mul for k=0 */

  square_nested_eq( t0 , 1.0 );
  square_nested_1_4( t2 , t6 , 1.0, 0.0);
  square_nested_eq( t4 , -1.0 );

  square_nested( t1 , t7 , fr, fi );
  square_nested_1_4( t3 , t5 , fr, fi );

#ifdef Y_PRELOAD_DATA

  radix_8_notwd_first_preload(ppu,t);
#else

  radix_8_notwd_first(pt,t,j);
#endif

  if(nc==0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-2-(int)(nc));nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>3;
      bj=(Y_LRIGHT[nr]>>3)-1;
#ifdef Y_PRELOAD_TRIGS

      trig_8_load_init(tpw,px,bi);
      trig_8_load_init(upw,py,bj);
#endif
#ifdef Y_PRELOAD_DATA

      cplx_radix_8_pad1_init(s,d,bi,ppul);
      cplx_radix_8_pad1_init(v,d,bj,ppdl);
#endif

      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_8_preload(px,t,tpw);
          wkr = t1r;
          wki = t1i;
          get_twiddle_factors_8_preload_down(py,u,upw);
          px += Y_STEP;
          py -= Y_STEP;
#else

          get_twiddle_factors_8_last_down(jj);
          l=jj<<3;
#endif
#ifdef Y_PRELOAD_DATA

          l = addr(((ii+1)<<4)) - addr((ii<<4));
          radix_8_twd_last_preload(t,s,ppu,ppul,l);
          l = addr(((jj)<<4)) - addr(((jj-1)<<4));
          radix_8_twd_last_preload_down(u,v,ppd,ppdl,l);
#else
# ifdef Y_LONG_MACROS

          radix_8_twd_last_down(u,pu,d,l);
# else

          radix_8_twd_last(u,pu,d,l);
# endif

          get_twiddle_factors_8_last(ii);
          l=ii<<3;
          radix_8_twd_last(t,pt,d,l);
# if defined(Y_MAXIMUM) && (Y_TARGET == 0)

          wkr=px[0];
          wki=px[1];
# else

          wkr=tw1r;
          wki=tw1i;
# endif
#endif
#ifndef Y_ITANIUM

          square_nested(t0, u7, wkr, wki);
          square_nested_1_4(t2, u5, wkr, wki);
          square_nested_1_2(t4, u3, wkr, wki);
          square_nested_3_4(t6, u1, wkr, wki);
#else

          square_nested_0_0_1_4(t0,u7,t2,u5,wkr, wki);
          square_nested_1_2_3_4(t4,u3,t6,u1,wkr, wki);
#endif

          aux=(wkr+wki)*fr;
          wki=(wki-wkr)*fr;
          wkr=aux;
#ifndef Y_ITANIUM

          square_nested(t1, u6, wkr, wki);
          square_nested_1_4(t3, u4, wkr, wki);
          square_nested_1_2(t5, u2, wkr, wki);
          square_nested_3_4(t7, u0, wkr, wki);
#else

          square_nested_0_0_1_4(t1,u6,t3,u4,wkr, wki);
          square_nested_1_2_3_4(t5,u2,t7,u0,wkr, wki);
#endif

#ifdef Y_PRELOAD_DATA

          radix_8_notwd_first_preload(ppu,t);
          radix_8_notwd_first_preload(ppd,u);
#else

          radix_8_notwd_first(pt,t,l);
          l=jj<<3;
          radix_8_notwd_first(pu,u,l);
#endif

        }
      if(ii==jj)
        {
#ifdef Y_PRELOAD_TRIGS
          trig_8_load_init(tpw,px,ii);
          get_twiddle_factors_8_preload(px,t,tpw);
          wkr = t1r;
          wki = t1i;
#else

          get_twiddle_factors_8_last(ii);
#endif
#ifdef Y_PRELOAD_DATA

          cplx_radix_8_pad1_init(s,d,ii,ppul);
          l = 16;
          radix_8_twd_last_preload(t,s,ppu,ppul,l);
#else

          l=ii<<3;
          radix_8_twd_last(t,pt,d,l);
# if defined(Y_MAXIMUM) && (Y_TARGET == 0)

          wkr=px[0];
          wki=px[1];
# else

          wkr=tw1r;
          wki=tw1i;
# endif
#endif

          square_nested(t0, t7, wkr, wki);
          square_nested_1_4(t2, t5, wkr,wki);
          aux=(wkr+wki)*fr;
          wki=(wki-wkr)*fr;
          wkr=aux;
          square_nested(t1, t6, wkr, wki);
          square_nested_1_4(t3, t4, wkr, wki);

#ifdef Y_PRELOAD_DATA

          radix_8_notwd_first_preload(ppu,t);
#else

          radix_8_notwd_first(pt,t,l);
#endif

        }
    }
}

/*
#define Y_PRELOAD_TRIGS
#define Y_PRELOAD_DATA
*/
void radix_8_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                  y_size_t jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i;
  y_limb_t fr=F_1_8r;
#ifdef Y_PRELOAD_DATA

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i;
  y_limb_t fi=F_1_8i;
#else

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i;
# if Y_PRELOAD_TRIGS

  y_limb_t uw1r,uw1i,uw2r,uw2i,uw3r,uw3i,uw4r,uw4i,uw5r,uw5i,uw6r,uw6i,
  uw7r,uw7i;
# endif
#endif
#ifdef Y_PRELOAD_TRIGS
# if defined(Y_MINIMUM)

  y_limb_t tpw1r,tpw1i;
  y_limb_t upw1r,upw1i;
# elif defined(Y_MAXIMUM)

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw3r,tpw3i,tpw4r,tpw4i,tpw5r,tpw5i,tpw6r,tpw6i,
  tpw7r,tpw7i;
  y_limb_t upw1r,upw1i,upw2r,upw2i,upw3r,upw3i,upw4r,upw4i,upw5r,upw5i,upw6r,upw6i,
  upw7r,upw7i;
# else

  y_limb_t tpw1r,tpw1i,tpw2r,tpw2i,tpw5r,tpw5i;
  y_limb_t upw1r,upw1i,upw2r,upw2i,upw5r,upw5i;
# endif
#endif

  y_limb_t wkr,wki,aux;
  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[7];
  y_ptr px;
#ifdef Y_PRELOAD_TRIGS

  y_ptr py;
#else

  y_ptr pt,pu;
#endif
#ifdef Y_PRELOAD_DATA

  y_ptr ppu0r,ppu0i,ppu1r,ppu1i,ppu2r,ppu2i,ppu3r,ppu3i,ppu4r,ppu4i,
  ppu5r,ppu5i,ppu6r,ppu6i,ppu7r,ppu7i;
  y_ptr ppul0r,ppul0i,ppul1r,ppul1i,ppul2r,ppul2i,ppul3r,ppul3i,ppul4r,ppul4i,
  ppul5r,ppul5i,ppul6r,ppul6i,ppul7r,ppul7i;
  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i,ppd4r,ppd4i,
  ppd5r,ppd5i,ppd6r,ppd6i,ppd7r,ppd7i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i,ppdl4r,ppdl4i,
  ppdl5r,ppdl5i,ppdl6r,ppdl6i,ppdl7r,ppdl7i;
#endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

  /* if(Y_PLAN[Y_NRADICES-1] != 8) return; */

  nb=bs>>3;
  bi=(bs*ii)>>3;
  bj=((bs*jj)>>3)+nb-1;
#ifdef Y_PRELOAD_TRIGS

  trig_8_load_init(tpw,px,bi);
  trig_8_load_init(upw,py,bj);
#endif
#ifdef Y_PRELOAD_DATA

  cplx_radix_8_pad1_init(s,d,bi,ppul);
  cplx_radix_8_pad1_init(v,d,bj,ppdl);
#endif
  /*printf("bi=%i , bj=%i, nb=%i\n",bi,bj,nb);*/
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
#ifdef Y_PRELOAD_TRIGS
      get_twiddle_factors_8_preload(px,t,tpw);
      wkr = t1r;
      wki = t1i;
      get_twiddle_factors_8_preload_down(py,u,upw);
      px += Y_STEP;
      py -= Y_STEP;
#else

      get_twiddle_factors_8_last_down(j);
      l=j<<3;
#endif
#ifdef Y_PRELOAD_DATA

      l = addr(((i+1)<<4)) - addr((i<<4));
      radix_8_twd_last_preload(t,s,ppu,ppul,l);
      l = addr(((j)<<4)) - addr(((j-1)<<4));
      radix_8_twd_last_preload_down(u,v,ppd,ppdl,l);
#else
# ifdef Y_LONG_MACROS

      radix_8_twd_last_down(u,pu,d,l);
# else

      radix_8_twd_last(u,pu,d,l);
# endif

      get_twiddle_factors_8_last(i);
#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

      wkr=px[0];
      wki=px[1];
#else

      wkr=tw1r;
      wki=tw1i;
#endif

      l=i<<3;
      radix_8_twd_last(t,pt,d,l);
#endif
#ifndef Y_ITANIUM

      square_nested(t0, u7, wkr, wki);
      square_nested_1_4(t2, u5, wkr, wki);
      square_nested_1_2(t4, u3, wkr, wki);
      square_nested_3_4(t6, u1, wkr, wki);
#else

      square_nested_0_0_1_4(t0,u7,t2,u5,wkr, wki);
      square_nested_1_2_3_4(t4,u3,t6,u1,wkr, wki);
#endif

      aux=(wkr+wki)*fr;
      wki=(wki-wkr)*fr;
      wkr=aux;

#ifndef Y_ITANIUM

      square_nested(t1, u6, wkr, wki);
      square_nested_1_4(t3, u4, wkr, wki);
      square_nested_1_2(t5, u2, wkr, wki);
      square_nested_3_4(t7, u0, wkr, wki);
#else

      square_nested_0_0_1_4(t1,u6,t3,u4,wkr, wki);
      square_nested_1_2_3_4(t5,u2,t7,u0,wkr, wki);
#endif
#ifdef Y_PRELOAD_DATA

      radix_8_notwd_first_preload(ppu,t);
      radix_8_notwd_first_preload(ppd,u);
#else

      radix_8_notwd_first(pt,t,l);
      l=j<<3;
      radix_8_notwd_first(pu,u,l);
#endif

    }
  if((i==j) && (ni < nb))
    {
#ifdef Y_PRELOAD_TRIGS
      trig_8_load_init(tpw,px,i);
      get_twiddle_factors_8_preload(px,t,tpw);
      wkr = t1r;
      wki = t1i;
#else

      get_twiddle_factors_8_last(i);
#endif
#ifdef Y_PRELOAD_DATA

      cplx_radix_8_pad1_init(s,d,i,ppul);
      l = 16;
      radix_8_twd_last_preload(t,s,ppu,ppul,l);
#else
#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

      wkr=px[0];
      wki=px[1];
#else

      wkr=tw1r;
      wki=tw1i;
#endif

      l=i<<3;
      radix_8_twd_last(t,pt,d,l);
#endif

      square_nested(t0, t7, wkr, wki);
      square_nested_1_4(t2, t5, wkr,wki);
      aux=(wkr+wki)*fr;
      wki=(wki-wkr)*fr;
      wkr=aux;
      square_nested(t1, t6, wkr, wki);
      square_nested_1_4(t3, t4, wkr, wki);

#ifdef Y_PRELOAD_DATA

      radix_8_notwd_first_preload(ppu,t);
#else

      radix_8_notwd_first(pt,t,l);
#endif

    }
}


void radix_8_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i;
  y_limb_t wkr,wki;
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[7];
  int nr;
  y_ptr px,pt,pu,ps,pv;

  ASSERT_ALIGNED_DOUBLE();

  /*if(Y_PLAN[Y_NRADICES-1] != 8) return;*/

  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 8, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 8 
  */

  /* Cycle 0 block 0 */
  j=0;
  get_twiddle_factors_8_last(j);
  radix_8_twd_last(s,ps,d1,j);
  radix_8_twd_last(t,pt,d2,j);
  /* dyadic mul for k=0 */
  conv_nested( s0, s0, t0, t0, 1.0, 0.0 );
  conv_nested_1_4( s2, s6, t2, t6, 1.0, 0.0);
  conv_nested_1_2( s4, s4, t4, t4, 1.0, 0.0 );

  wkr=F_1_8r;
  wki=F_1_8i;
  conv_nested( s1, s7, t1, t7, wkr, wki );
  conv_nested_1_4( s3, s5, t3, t5, wkr, wki );

  radix_8_notwd_first(ps,s,j);
  if(nc == 0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>3;
      bj=(Y_LRIGHT[nr]>>3)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_8_last_down(jj);
          l=jj<<3;
#ifdef Y_LONG_MACROS

          radix_8_twd_last_down(u,pu,d1,l);
          radix_8_twd_last_down(v,pv,d2,l);
#else

          radix_8_twd_last(u,pu,d1,l);
          radix_8_twd_last(v,pv,d2,l);
#endif

          get_twiddle_factors_8_last(ii);
          l=ii<<3;
          radix_8_twd_last(s,ps,d1,l);
          radix_8_twd_last(t,pt,d2,l);

#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

          conv_nested( s0, u7, t0, v7, px[0], px[1]);
          conv_nested_1_4( s2, u5, t2, v5, px[0], px[1]);
          conv_nested_1_2( s4, u3, t4, v3, px[0], px[1]);
          conv_nested_3_4( s6, u1, t6, v1, px[0], px[1]);

          wkr=(px[0] + px[1])*F_1_8r;
          wki=(px[1] - px[0])*F_1_8r;
#else

          conv_nested( s0, u7, t0, v7, tw1r, tw1i);
          conv_nested_1_4( s2, u5, t2, v5, tw1r, tw1i);
          conv_nested_1_2( s4, u3, t4, v3, tw1r, tw1i);
          conv_nested_3_4( s6, u1, t6, v1, tw1r, tw1i);

          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
#endif

          conv_nested( s1, u6, t1, v6, wkr, wki);
          conv_nested_1_4( s3, u4, t3, v4, wkr, wki);
          conv_nested_1_2( s5, u2, t5, v2, wkr, wki);
          conv_nested_3_4( s7, u0, t7, v0, wkr, wki);

          radix_8_notwd_first(ps,s,l);
          l=jj<<3;
          radix_8_notwd_first(pu,u,l);
        }
      if(ii==jj)
        {
          get_twiddle_factors_8_last(ii);
          l=ii<<3;
          radix_8_twd_last(s,ps,d1,l);
          radix_8_twd_last(t,pt,d2,l);
#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

          conv_nested(s0, s7, t0, t7, px[0], px[1]);
          conv_nested_1_4(s2, s5, t2, t5, px[0], px[1]);
          wkr=(px[0] + px[1])*F_1_8r;
          wki=(px[1] - px[0])*F_1_8r;
#else

          conv_nested(s0, s7, t0, t7, tw1r, tw1i);
          conv_nested_1_4(s2, s5, t2, t5, tw1r, tw1i);
          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
#endif

          conv_nested(s1, s6, t1, t6, wkr, wki);
          conv_nested_1_4(s3, s4, t3, t4, wkr, wki);
          radix_8_notwd_first(ps,s,l);
        }
    }
}

void radix_8_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                               y_size_t  jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i;
  y_limb_t wkr,wki;
  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[7];
  y_ptr px,pu,pt,ps,pv;

  ASSERT_ALIGNED_DOUBLE();


  /*  if(Y_PLAN[Y_NRADICES-1] != 8) return;*/

  nb=bs>>3;
  bi=(bs*ii)>>3;
  bj=((bs*jj)>>3)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      get_twiddle_factors_8_last_down(j);
      l=j<<3;
#ifdef Y_LONG_MACROS

      radix_8_twd_last_down(u,pu,d1,l);
      radix_8_twd_last_down(v,pv,d2,l);
#else

      radix_8_twd_last(u,pu,d1,l);
      radix_8_twd_last(v,pv,d2,l);
#endif

      get_twiddle_factors_8_last(i);
      l=i<<3;
      radix_8_twd_last(s,ps,d1,l);
      radix_8_twd_last(t,pt,d2,l);

#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

      conv_nested( s0, u7, t0, v7, px[0], px[1]);
      conv_nested_1_4( s2, u5, t2, v5, px[0], px[1]);
      conv_nested_1_2( s4, u3, t4, v3, px[0], px[1]);
      conv_nested_3_4( s6, u1, t6, v1, px[0], px[1]);

      wkr=(px[0] + px[1])*F_1_8r;
      wki=(px[1] - px[0])*F_1_8r;
#else

      conv_nested( s0, u7, t0, v7, tw1r, tw1i);
      conv_nested_1_4( s2, u5, t2, v5, tw1r, tw1i);
      conv_nested_1_2( s4, u3, t4, v3, tw1r, tw1i);
      conv_nested_3_4( s6, u1, t6, v1, tw1r, tw1i);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
#endif

      conv_nested( s1, u6, t1, v6, wkr, wki);
      conv_nested_1_4( s3, u4, t3, v4, wkr, wki);
      conv_nested_1_2( s5, u2, t5, v2, wkr, wki);
      conv_nested_3_4( s7, u0, t7, v0, wkr, wki);

      radix_8_notwd_first(ps,s,l);
      l=j<<3;
      radix_8_notwd_first(pu,u,l);
    }
  if((i==j)&&(ni<nb))
    {
      get_twiddle_factors_8_last(i);
      l=i<<3;
      radix_8_twd_last(s,ps,d1,l);
      radix_8_twd_last(t,pt,d2,l);
#if defined(Y_MAXIMUM) && (Y_TARGET == 0)

      conv_nested(s0, s7, t0, t7, px[0], px[1]);
      conv_nested_1_4(s2, s5, t2, t5, px[0], px[1]);

      wkr=(px[0] + px[1])*F_1_8r;
      wki=(px[1] - px[0])*F_1_8r;
#else

      conv_nested(s0, s7, t0, t7, tw1r, tw1i);
      conv_nested_1_4(s2, s5, t2, t5, tw1r, tw1i);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
#endif

      conv_nested(s1, s6, t1, t6, wkr, wki);
      conv_nested_1_4(s3, s4, t3, t4, wkr, wki);

      radix_8_notwd_first(ps,s,l);
    }
}

/*$Id$*/








