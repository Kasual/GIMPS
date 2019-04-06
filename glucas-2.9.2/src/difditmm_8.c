/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2002-2006 Guillermo Ballester Valor
 
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
 
 
   This is the version using SSE2 extensions. The plan is as follows
 
   a) 
     Excluding the first block. Last DIF pass is made by in pairs, 
     i.e., the two 8-complex blocks will be executed in parallel. 
 
     It loads memory data in legacy mode and returns the intermediate data
     also in legacy mode: a complex will be represented by an Y__M128D 
     data. 
 
   b) 
     The dyadic square phase is made also by pairs. It takes complex_data 
     as Y_M128D data, makes the dyadic work and return again data in this  
     format. Named 'txr' data will be the old 't' data and 'txi' the old 'u'  
 
   c) 
     Again, first DIT pass is made by pairs. It takes complex_data as 
     Y_M128D data and stores the data in memory in legacy format.
 
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

#if defined(Y_USE_SSE2)
# include "ygensse2.h"

# define NDEBUG1

# define set_sse2_8_constants         \
  Y_MM_UNPACKLO_PD (MM_AUX, MM_F_1_8r, MM_F_1_8i);

/*
# define set_sse2_8_constants         \
  Y_MM_SET1_PD ( MM_F_1_8r, F_1_8r);  \
  Y_MM_SET1_PD ( MM_F_1_8i, F_1_8i);  \
  Y_MM_SET1_PD ( MM_1_0, 1.0);        \
  Y_MM_SET1_PD ( MM_0_5, 0.5);        \
  Y_MM_UNPACKLO_PD (MM_AUX, MM_F_1_8r, MM_F_1_8i); 
*/

# if defined(Y_MINIMUM)
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 2

#  include "yminimum.h"

#  define sse2_get_twiddle_factors_8_last(_tw, _i, _j)      \
          px = &( Y_TWDF[Y_NRADICES-2][( _i * inc) << 1 ] );\
          py = &( Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ] );\
          Y_MM_LOAD_INTER_PAIR(_tw##1r, _tw##1i, px, py);   \
          sse2_square( _tw##2, _tw##1);                     \
          sse2_square( _tw##4, _tw##2);                     \
          sse2_divmul( _tw##3, _tw##5, _tw##4, _tw##1);     \
          sse2_square( _tw##6, _tw##3);                     \
          sse2_mul( _tw##7, _tw##4, _tw##3);                \
          prefetch_p_trig_nta(px);                          \
          prefetch_data_trig_nta(py,(-2));

#  define get_twiddle_factors_8_last(_tw, _j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          cplx_trig_min_load(_tw##1,px);            \
          cplx_trig_min_square(_tw##2,_tw##1);      \
          cplx_trig_min_square(_tw##4,_tw##2);      \
          cplx_divmul(_tw##3,_tw##5,_tw##4,_tw##1); \
          cplx_trig_min_square(_tw##6,_tw##3);      \
          cplx_trig_min_mul(_tw##7,_tw##4,_tw##3);  \
          prefetch_p_trig_nta(px);

#  define get_twiddle_factors_8_last_down(_tw, _j)  \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          cplx_trig_min_load(_tw##1,px);            \
          cplx_trig_min_square(_tw##2,_tw##1);      \
          cplx_trig_min_square(_tw##4,_tw##2);      \
          cplx_divmul(_tw##3,_tw##5,_tw##4,_tw##1); \
          cplx_trig_min_square(_tw##6,_tw##3);      \
          cplx_trig_min_mul(_tw##7,_tw##4,_tw##3);  \
          prefetch_data_trig_nta(px,(-2));


# elif defined(Y_MAXIMUM)
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 14


#  define sse2_get_twiddle_factors_8_last(_tw, _i, _j)              \
          px = &( Y_TWDF[Y_NRADICES-2][( _i * inc) << 1 ] );        \
          py = &( Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ] );        \
          Y_MM_LOAD_INTER_PAIR(_tw##1r, _tw##1i, px, py);           \
          Y_MM_LOAD_INTER_PAIR(_tw##2r, _tw##2i, px + 2, py + 2);   \
          Y_MM_LOAD_INTER_PAIR(_tw##3r, _tw##3i, px + 4, py + 4);   \
          Y_MM_LOAD_INTER_PAIR(_tw##4r, _tw##4i, px + 6, py + 6);   \
          Y_MM_LOAD_INTER_PAIR(_tw##5r, _tw##5i, px + 8, py + 8);   \
          Y_MM_LOAD_INTER_PAIR(_tw##6r, _tw##6i, px + 10, py + 10); \
          Y_MM_LOAD_INTER_PAIR(_tw##7r, _tw##7i, px + 12, py + 12); \
          prefetch_data_trig_nta( px, 14 + Y_CACHE_LINE);           \
          prefetch_data_trig_nta( px, 14 + 2*Y_CACHE_LINE);         \
          prefetch_data_trig_nta( py, (-14) + Y_CACHE_LINE);        \
          prefetch_data_trig_nta( py, (-14) + 2*Y_CACHE_LINE);


#  define get_twiddle_factors_8_last(_tw, _j)         \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);  \
	  _tw##1r=px[0];   _tw##1i=px[1];             \
          _tw##2r=px[2];   _tw##2i=px[3];             \
          _tw##3r=px[4];   _tw##3i=px[5];             \
          _tw##4r=px[6];   _tw##4i=px[7];             \
          _tw##5r=px[8];   _tw##5i=px[9];             \
          _tw##6r=px[10];  _tw##6i=px[11];            \
          _tw##7r=px[12];  _tw##7i=px[13];            \
          prefetch_data_trig_nta( px, 14 + Y_CACHE_LINE); \
          prefetch_data_trig_nta( px, 14 + 2*Y_CACHE_LINE);

#  define get_twiddle_factors_8_last_down(_tw, _j)      \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);    \
	  _tw##1r=px[0];   _tw##1i=px[1];               \
          _tw##2r=px[2];   _tw##2i=px[3];               \
          _tw##3r=px[4];   _tw##3i=px[5];               \
          _tw##4r=px[6];   _tw##4i=px[7];               \
          _tw##5r=px[8];   _tw##5i=px[9];               \
          _tw##6r=px[10];  _tw##6i=px[11];              \
          _tw##7r=px[12];  _tw##7i=px[13];              \
          prefetch_data_trig_nta( px, (-14) + Y_CACHE_LINE);\
          prefetch_data_trig_nta( px, (-14) + 2*Y_CACHE_LINE);

# else
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 6

#  define sse2_get_twiddle_factors_8_last(_tw, _i, _j)            \
          px = &( Y_TWDF[Y_NRADICES-2][( _i * inc) << 1 ] );      \
          py = &( Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ] );      \
          Y_MM_LOAD_INTER_PAIR(_tw##1r, _tw##1i, px, py);         \
          Y_MM_LOAD_INTER_PAIR(_tw##2r, _tw##2i, px + 2, py + 2); \
          Y_MM_LOAD_INTER_PAIR(_tw##5r, _tw##5i, px + 4, py + 4); \
          sse2_divmul( _tw##3, _tw##7, _tw##5, _tw##2);           \
          sse2_divmul( _tw##4, _tw##6, _tw##5, _tw##1);           \
          prefetch_data_trig_nta( px, Y_STEP);                    \
          prefetch_data_trig_nta( px, Y_STEP + Y_CACHE_LINE);     \
          prefetch_data_trig_nta( py, (-6));                      \
          prefetch_data_trig_nta( py, (-2));

#  define get_twiddle_factors_8_last(_tw, _j)       \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  _tw##1r=px[0];   _tw##1i=px[1];           \
          _tw##2r=px[2];   _tw##2i=px[3];           \
          _tw##5r=px[4];   _tw##5i=px[5];           \
          cplx_divmul(_tw##3,_tw##7,_tw##5,_tw##2); \
          cplx_divmul(_tw##4,_tw##6,_tw##5,_tw##1); \
          prefetch_data_trig_nta( px, Y_STEP);      \
          prefetch_data_trig_nta( px, Y_STEP + Y_CACHE_LINE);

#  define get_twiddle_factors_8_last_down(_tw, _j)  \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  _tw##1r=px[0];   _tw##1i=px[1];           \
          _tw##2r=px[2];   _tw##2i=px[3];           \
          _tw##5r=px[4];   _tw##5i=px[5];           \
          cplx_divmul(_tw##3,_tw##7,_tw##5,_tw##2); \
          cplx_divmul(_tw##4,_tw##6,_tw##5,_tw##1); \
          prefetch_data_trig( px, (-6));            \
          prefetch_data_trig( px, (-2));

# endif



#define sse2_mul_1_8_B_addsub_store_inter_no_next(_pd, _pu, _i, _j, _t0, _t1)\
{                                                                            \
   Y__M128D __a, __b;                                                        \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i );                                       \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r );                                       \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );                                      \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );                                      \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                        \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                        \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                        \
   Y_MM_STORE_INTER_PAIR( _pd + _j, _pu + _j, _t1##r, _t1##i);               \
   Y_MM_STORE_INTER_PAIR( _pd + _i, _pu + _i, _t0##r, _t0##i);               \
}


#define sse2_mul_3_8_B_addsub_store_inter_no_next(_pd, _pu, _i, _j, _t0, _t1)\
{                                                                            \
   Y__M128D __a, __b;                                                        \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i );                                       \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r );                                       \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );                                      \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );                                      \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                        \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                        \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                        \
   Y_MM_STORE_INTER_PAIR( _pd + _j, _pu + _j, _t1##r, _t1##i);               \
   Y_MM_STORE_INTER_PAIR( _pd + _i, _pu + _i, _t0##r, _t0##i);               \
}


# define sse2_radix_8_twd_last(_z, _tw, _pt, _pu, _d, _i, _j)                                  \
          _pt = _d + addr( (_i) << 1);                                                         \
          _pu = _d + addr( (_j) << 1);                                                         \
	  sse2_load_inter_muladdsub_no_next( _z##0, _z##1, _tw##4, 0, 8, _pt, _pu);            \
	  sse2_load_inter_mulmuladdsub_no_next( _z##2, _z##3, _tw##2, _tw##6, 4, 12, _pt, _pu);\
	  sse2_addsub( _z##0 ,_z##2 );                                                         \
	  sse2_mul_1_4_F_addsub( _z##1 , _z##3);                                               \
	  sse2_load_inter_mulmuladdsub_no_next( _z##4, _z##5, _tw##1, _tw##5, 2, 10, _pt, _pu);\
	  sse2_load_inter_mulmuladdsub_no_next( _z##6, _z##7, _tw##3, _tw##7, 6, 14, _pt, _pu);\
	  sse2_addsub( _z##4 ,_z##6 );                                                         \
	  sse2_mul_1_4_F_addsub( _z##5 , _z##7);                                               \
	                                                                                       \
	  sse2_addsub( _z##0, _z##4 );                                                         \
          Y_MM_INTER_PAIR (_z##0r, _z##0i);                                                    \
          Y_MM_INTER_PAIR (_z##4r, _z##4i);                                                    \
	  sse2_mul_1_8_F_addsub(_z##1 , _z##5 );                                               \
          Y_MM_INTER_PAIR (_z##1r, _z##1i);                                                    \
          Y_MM_INTER_PAIR (_z##5r, _z##5i);                                                    \
	  sse2_mul_1_4_F_addsub(_z##2 , _z##6 );                                               \
          Y_MM_INTER_PAIR (_z##2r, _z##2i);                                                    \
          Y_MM_INTER_PAIR (_z##6r, _z##6i);                                                    \
	  sse2_mul_3_8_F_addsub(_z##3 , _z##7 );                                               \
          Y_MM_INTER_PAIR (_z##3r, _z##3i);                                                    \
          Y_MM_INTER_PAIR (_z##7r, _z##7i);


# define radix_8_twd_last(_z,_tw,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
	  cplx_load_muladdsub_pp( _z##0 , _z##1 , _tw##4, 0, _pd, _pd + 8);\
          cplx_load_mulmuladdsub_pp( _z##2,_z##3, _tw##2, _tw##6, 0,_pd + 4,_pd + 12);\
	  cplx_addsub( _z##0 ,_z##2 );\
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3);\
          cplx_load_mulmuladdsub_pp( _z##4,_z##5, _tw##1, _tw##5, 0,_pd + 2,_pd + 10);\
          cplx_load_mulmuladdsub_pp( _z##6,_z##7, _tw##3, _tw##7, 0,_pd + 6,_pd + 14);\
          prefetch_p( _pd + 8);\
          prefetch_p( _pd + 16);\
	  cplx_addsub( _z##4 ,_z##6 );\
	  cplx_mul_1_4_F_addsub( _z##5 , _z##7);\
	  \
	  cplx_addsub( _z##0, _z##4 );\
	  cplx_mul_1_8_F_addsub(_z##1 , _z##5 );\
	  cplx_mul_1_4_F_addsub(_z##2 , _z##6 );\
	  cplx_mul_3_8_F_addsub(_z##3 , _z##7 );


# define sse2_radix_8_notwd_first(_pd, _pu, _z)                                       \
          Y_MM_INTER_PAIR (_z##0r, _z##0i);                                           \
          Y_MM_INTER_PAIR (_z##4r, _z##4i);                                           \
	  sse2_addsub( _z##0 , _z##4 );                                               \
          Y_MM_INTER_PAIR (_z##2r, _z##2i);                                           \
          Y_MM_INTER_PAIR (_z##6r, _z##6i);                                           \
	  sse2_addsub( _z##2 , _z##6 );                                               \
	  sse2_addsub( _z##0 , _z##2 );                                               \
	  sse2_mul_1_4_B_addsub( _z##4 , _z##6 );                                     \
          Y_MM_INTER_PAIR (_z##1r, _z##1i);                                           \
          Y_MM_INTER_PAIR (_z##5r, _z##5i);                                           \
	  sse2_addsub( _z##1 , _z##5 );                                               \
          Y_MM_INTER_PAIR (_z##3r, _z##3i);                                           \
          Y_MM_INTER_PAIR (_z##7r, _z##7i);                                           \
	  sse2_addsub( _z##3 , _z##7 );                                               \
	  sse2_addsub( _z##1 , _z##3 );                                               \
          sse2_mul_1_4_B_addsub( _z##5 , _z##7 );                                     \
                                                                                      \
          sse2_addsub_store_inter_no_next( _pd, _pu, 0, 8, _z##0, _z##1);             \
          sse2_mul_1_8_B_addsub_store_inter_no_next( _pd, _pu, 2, 10, _z##4, _z##5);  \
          sse2_mul_1_4_B_addsub_store_inter_no_next( _pd, _pu, 4, 12, _z##2, _z##3);  \
          sse2_mul_3_8_B_addsub_store_inter_no_next( _pd, _pu, 6, 14, _z##6, _z##7);


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
	  cplx_addsub_store_p(0 , _pd, _pd + 8, _z##0 , _z##1 );\
	  cplx_mul_1_8_B_addsub_store_p(0,_pd + 2,_pd + 10, _z##4 , _z##5 );\
	  cplx_mul_1_4_B_addsub_store_p(0,_pd + 4,_pd + 12, _z##2 , _z##3 );\
	  cplx_mul_3_8_B_addsub_store_p(0,_pd + 6,_pd + 14, _z##6 , _z##7 );

/* debug macros, only in emulation mode */
#define sse2_to_double(_t, _s) \
          _t##r = _s.d[0];      \
          _t##i = _s.d[1];

#define double_to_sse2(_t, _s) \
          _t.d[0] = _s##r;      \
          _t.d[1] = _s##i;


/*
#undef Y_PRELOAD_TRIGS
#undef Y_PRELOAD_DATA
#undef Y_ITANIUM
*/

/*
   this function only does the work for first block0 
   It does not use sse2 extensions
*/
void radix_8_dif_square_dit_block0(y_ptr d)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i;
  y_limb_t fr=F_1_8r, fi=F_1_8i;
  y_size_t j, inc = Y_POWERS[7];
  y_ptr pt, px;

  /* Cycle 0 block 0 */
  j=0;
  /* TO DO:
     we can simplify this because we don't need twidle factors 
     in block 0 (are ones)
  */
  get_twiddle_factors_8_last(tw, j);
  radix_8_twd_last(t,tw,pt,d,j);

  /* dyadic mul for k=0 */

  square_nested_eq( t0 , 1.0 );
  square_nested_1_4( t2 , t6 , 1.0, 0.0);
  square_nested_eq( t4 , -1.0 );

  square_nested( t1 , t7 , fr, fi );
  square_nested_1_4( t3 , t5 , fr, fi );

  radix_8_notwd_first(pt,t,j);

}

void radix_8_dif_square_dit_blockE(y_ptr d, y_size_t ii)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,wkr,wki,aux;
  y_limb_t fr=F_1_8r;
  y_size_t l, inc = Y_POWERS[7];
  y_ptr pt, px;

  get_twiddle_factors_8_last(tw, ii);
  l=ii<<3;
  radix_8_twd_last(t,tw,pt,d,l);

  wkr=tw1r;
  wki=tw1i;

  square_nested(t0, t7, wkr, wki);
  square_nested_1_4(t2, t5, wkr,wki);

  aux=(wkr+wki)*fr;
  wki=(wki-wkr)*fr;
  wkr=aux;

  square_nested(t1, t6, wkr, wki);
  square_nested_1_4(t3, t4, wkr, wki);

  radix_8_notwd_first(pt,t,l);
}

void radixmm_8_dif_square_dit(y_ptr d, y_size_t nc)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i,
  t4r, t4i, t5r, t5i, t6r, t6i, t7r, t7i;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i,
  tw5r, tw5i, tw6r, tw6i, tw7r, tw7i, MM_AUX;
  Y__M128D gkr, gki;

  y_size_t bi, bj, jj, ii, li, lj, inc = Y_POWERS[7];
  int nr;
  y_ptr px, py;
  y_ptr pt, pu;

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

  /* Cycle 0 block 0. In this version, left to
     radix_8_dif_square_dit_block0()  */

  radix_8_dif_square_dit_block0 ( d );

  if (nc == 0)
    return;

  /* set some constants */
  set_sse2_8_constants;

  for (nr = Y_NRADICES-2; nr > ((int)(Y_NRADICES) - 2 - (int)(nc)); nr-- )
    {
      bi = (Y_LRIGHT[nr+1]) >> 3;
      bj = (Y_LRIGHT[nr] >> 3) - 1;
      for(ii = bi, jj = bj; ii < jj; ii++, jj--)
        {
          sse2_get_twiddle_factors_8_last (tw, ii, jj);

          Y_MM_UNPACKLO_PD (gkr, tw1r, tw1i);
          Y_MM_COMPLEX_MUL (gki, gkr, MM_AUX);
          Y_MM_INTER_PAIR (gkr, gki);

          li = ii << 3;
          lj = jj << 3;

          sse2_radix_8_twd_last (t, tw, pt, pu, d, li, lj);

          sse2_square_nested (t0r, t7i, t1r, t6i, gkr, gki);
          sse2_square_nested_1_4 (t2r, t5i, t3r, t4i, gkr, gki);
          sse2_square_nested_1_2 (t4r, t3i, t5r, t2i, gkr, gki);
          sse2_square_nested_3_4 (t6r, t1i, t7r, t0i, gkr, gki);

          sse2_radix_8_notwd_first ( pt, pu, t);
        }
      if(ii == jj)
        {
          radix_8_dif_square_dit_blockE ( d , ii);
        }
    }
}

/*
#define Y_PRELOAD_TRIGS
#define Y_PRELOAD_DATA
*/

void radixmm_8_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                    y_size_t jj)
{
  Y__M128D t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i,
  t4r, t4i, t5r, t5i, t6r, t6i, t7r, t7i;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i,
  tw5r, tw5i, tw6r, tw6i, tw7r, tw7i, MM_AUX;
  Y__M128D gkr, gki;

  y_size_t bi, bj, j, i, li, lj, nb, ni,  inc = Y_POWERS[7];

  y_ptr px, py;
  y_ptr pt, pu;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

  /* if(Y_PLAN[Y_NRADICES-1] != 8) return; */

  nb=bs>>3;
  bi=(bs*ii)>>3;
  bj=((bs*jj)>>3)+nb-1;

  /* set some constants */
  set_sse2_8_constants;

  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      sse2_get_twiddle_factors_8_last (tw, i, j);

      Y_MM_UNPACKLO_PD (gkr, tw1r, tw1i);
      Y_MM_COMPLEX_MUL (gki, gkr, MM_AUX);
      Y_MM_INTER_PAIR (gkr, gki);

      li = i << 3;
      lj = j << 3;

      sse2_radix_8_twd_last (t, tw, pt, pu, d, li, lj);

      sse2_square_nested (t0r, t7i, t1r, t6i, gkr, gki);
      sse2_square_nested_1_4 (t2r, t5i, t3r, t4i, gkr, gki);
      sse2_square_nested_1_2 (t4r, t3i, t5r, t2i, gkr, gki);
      sse2_square_nested_3_4 (t6r, t1i, t7r, t0i, gkr, gki);

      sse2_radix_8_notwd_first ( pt, pu, t);
    }
  if((i == j) && (ni < nb))
    {
      radix_8_dif_square_dit_blockE ( d , i);
    }
}

#if defined(_KK)
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
#endif /* _kk */
#else /* defined Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radix_8_square_void( void )
{
  printf("The routine radix_8_square_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif /* Y_USE_SSE2 */

/*$Id$*/








