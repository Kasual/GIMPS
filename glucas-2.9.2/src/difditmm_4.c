/*$Id$*/
/*  This file is a part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2003-2006 Guillermo Ballester Valor
 
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
   It includes routines for a 4_last_radix reduction 
   
   This version includes SSE2 code 
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"

#if defined(Y_USE_SSE2)
# include "ygensse2.h"

#define NDEBUG1

# define set_sse2_4_constants         \
  Y_MM_SET1_PD ( MM_1_0, 1.0);        \
  Y_MM_SET1_PD ( MM_0_5, 0.5);        \
  Y_MM_SET_SD (MM_AUX, (-1.0));


#if defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

#  define sse2_get_twiddle_factors_4_last(_tw, _i, _j)            \
          px = &( Y_TWDF[Y_NRADICES-2][( _i * inc) << 1 ] );      \
          py = &( Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ] );      \
          Y_MM_LOAD_INTER_PAIR(_tw##1r, _tw##1i, px, py);         \
          Y_MM_LOAD_INTER_PAIR(_tw##2r, _tw##2i, px + 2, py + 2); \
          Y_MM_LOAD_INTER_PAIR(_tw##3r, _tw##3i, px + 4, py + 4); \
          prefetch_p_trig_nta(px);                                \
          prefetch_p_trig_nta(px + Y_CACHE_LINE);                 \
          prefetch_data_trig_nta(py, (-6));                       \
          prefetch_data_trig_nta(py, (-6 + Y_CACHE_LINE));


# define get_twiddle_factors_4_last(_tw, _j)        \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          _tw##1r=px[0];      _tw##1i=px[1];        \
          _tw##2r=px[2];      _tw##2i=px[3];        \
          _tw##3r=px[4];      _tw##3i=px[5];        \
          prefetch_p_trig_nta(px);                  \
          prefetch_p_trig_nta(px + Y_CACHE_LINE);

# define get_twiddle_factors_4_last_down(_j)        \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          _tw##1r=px[0];      _tw##1i=px[1];        \
          _tw##2r=px[2];      _tw##2i=px[3];        \
          _tw##3r=px[4];      _tw##3i=px[5];        \
          prefetch_data_trig_nta(px, (-6));         \
          prefetch_data_trig_nta(px, (-6 + Y_CACHE_LINE));


#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

#  define sse2_get_twiddle_factors_4_last(_tw, _i, _j)        \
          px = &( Y_TWDF[Y_NRADICES-2][( _i * inc) << 1 ] );  \
          py = &( Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ] );  \
          Y_MM_LOAD_INTER_PAIR(_tw##1r, _tw##1i, px, py);     \
          sse2_square( _tw##2, _tw##1);                       \
          sse2_mul (_tw##3, _tw##2, _tw##1);                  \
          prefetch_p_trig_nta(px);                            \
          prefetch_data_trig_nta(py,(-2));



# define get_twiddle_factors_4_last(_tw, _j)            \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);    \
	  _tw##2r=(px[0]+px[1])*(px[0]-px[1]);          \
	  _tw##2i=2.0*px[0]*px[1];                      \
	  _tw##1r=px[0];                                \
	  _tw##1i=px[1];                                \
	  _tw##3r=_tw##1r * _tw##2r - _tw##1i * _tw##2i;\
	  _tw##3i=_tw##1r * _tw##2i + _tw##1i * _tw##2r;\
          prefetch_p_trig_nta(px);

# define get_twiddle_factors_4_last_down(_j)            \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);    \
	  _tw##2r=(px[0]+px[1])*(px[0]-px[1]);          \
	  _tw##2i=2.0*px[0]*px[1];                      \
	  _tw##1r=px[0];                                \
	  _tw##1i=px[1];                                \
	  _tw##3r=_tw##1r * _tw##2r - _tw##1i * _tw##2i;\
	  _tw##3i=_tw##1r * _tw##2i + _tw##1i * _tw##2r;\
          prefetch_data_trig_nta(px,(-2));

#endif

# define sse2_radix_4_twd_last(_z, _tw, _pt, _pu, _d, _i, _j)                                  \
          _pt = _d + addr( (_i) << 1);                                                         \
          _pu = _d + addr( (_j) << 1);                                                         \
	  sse2_load_inter_muladdsub_no_next( _z##0, _z##1, _tw##2, 0, 4, _pt, _pu);            \
	  sse2_load_inter_mulmuladdsub_no_next( _z##2, _z##3, _tw##1, _tw##3, 2, 6, _pt, _pu); \
          prefetch_p(_pt + 4);                                                                 \
          prefetch_p(_pu + 4);                                                                 \
	  sse2_addsub( _z##0 ,_z##2 );                                                         \
          Y_MM_INTER_PAIR (_z##0r, _z##0i);                                                    \
          Y_MM_INTER_PAIR (_z##2r, _z##2i);                                                    \
	  sse2_mul_1_4_F_addsub( _z##1 , _z##3);                                               \
          Y_MM_INTER_PAIR (_z##1r, _z##1i);                                                    \
          Y_MM_INTER_PAIR (_z##3r, _z##3i);                                                    \


# define radix_4_twd_last(_z, _tw, _pd, _d, _j)                                        \
          _pd = _d + addr((_j)<<1);                                                    \
	  cplx_load_muladdsub_pp( _z##0 , _z##1 ,_tw##2, 0, _pd, _pd + 4);             \
	  cplx_load_mulmuladdsub_pp(_z##2 ,_z##3 ,_tw##1, _tw##3, 0, _pd + 2, _pd + 6);\
	                                                                               \
          prefetch_p(_pd+4);                                                           \
	  cplx_addsub( _z##0 ,_z##2 );                                                 \
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3);



# define sse2_radix_4_notwd_first(_pd, _pu, _z)                                       \
          Y_MM_INTER_PAIR (_z##0r, _z##0i);                                           \
          Y_MM_INTER_PAIR (_z##2r, _z##2i);                                           \
	  sse2_addsub( _z##0 , _z##2 );                                               \
          Y_MM_INTER_PAIR (_z##1r, _z##1i);                                           \
          Y_MM_INTER_PAIR (_z##3r, _z##3i);                                           \
	  sse2_addsub( _z##1 , _z##3 );                                               \
                                                                                      \
          sse2_addsub_store_inter_no_next( _pd, _pu, 0, 4, _z##0, _z##1);             \
          sse2_mul_1_4_B_addsub_store_inter_no_next( _pd, _pu, 2, 6, _z##2, _z##3);

# define radix_4_notwd_first(_pd,_z)                                       \
	  cplx_addsub( _z##0 , _z##2 );                                    \
	  cplx_addsub( _z##1 , _z##3 );                                    \
	                                                                   \
	  cplx_addsub_store_p(0, _pd,_pd + 4, _z##0 , _z##1 );             \
	  cplx_mul_1_4_B_addsub_store_p(0,_pd + 2,_pd + 6, _z##2 , _z##3 );



void radix_4_dif_square_dit_block0(y_ptr d)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;

  y_size_t j,inc=Y_POWERS[3];

  y_ptr px,pt;

  /* Cycle 0 block 0 */
  j=0;
  get_twiddle_factors_4_last(tw, j);
  radix_4_twd_last(t, tw, pt, d, j);
  /* dyadic mul for k=0 */

  square_nested_eq( t0 , 1.0 );

  square_nested_1_4( t1 , t3 ,1.0, 0.0);

  square_nested_eq( t2 , -1.0 );
  radix_4_notwd_first(pt,t);
}

void radix_4_dif_square_dit_blockE(y_ptr d, y_size_t ii)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;

  y_size_t l,inc=Y_POWERS[3];

  y_ptr px,pt;

  get_twiddle_factors_4_last(tw, ii);
  l=ii<<2;
  radix_4_twd_last(t, tw, pt, d, l);
  square_nested(t0, t3, tw1r, tw1i);
  square_nested_1_4(t1, t2, tw1r, tw1i);

  radix_4_notwd_first(pt,t);
}


void radixmm_4_dif_square_dit(y_ptr d, y_size_t nc)
{
  Y__M128D  t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i;
  Y__M128D MM_1_0, MM_0_5, MM_AUX, gkr, gki;
  y_size_t bi, bj, jj, ii, li, lj, inc=Y_POWERS[3];
  int nr;
  y_ptr px, py, pt, pu;

  ASSERT_ALIGNED_DOUBLE();

  if(Y_PLAN[Y_NRADICES-1] != 4)
    return;
  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 4, so we have
    to compute which blocks we need
    
    The first cycle, asociated   Y__M128D  t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i;
  Y__M128D MM_1_0, MM_0_5, MM_AUX, gkr, gki;
  y_size_t bi, bj, jj, ii, li, lj, inc=Y_POWERS[3];
  int nr;
  y_ptr px, py, pt, pu;
  with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 4 
  */

  radix_4_dif_square_dit_block0 ( d );

  if(nc==0)
    return;

  /* set some constants */
  set_sse2_4_constants;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-2-(int)nc);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>2;
      bj=(Y_LRIGHT[nr]>>2)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          sse2_get_twiddle_factors_4_last (tw, ii, jj);

          Y_MM_UNPACKLO_PD (gkr, tw1r, tw1i);
          Y_MM_COMPLEX_MUL (gki, gkr, MM_AUX);
          Y_MM_INTER_PAIR (gkr, gki);

          li = ii << 2;
          lj = jj << 2;

          sse2_radix_4_twd_last (t, tw, pt, pu, d, li, lj);

          sse2_square_nested (t0r, t3i, t2r, t1i, gkr, gki);
          sse2_square_nested_1_4 (t1r, t2i, t3r, t0i, gkr, gki);

          sse2_radix_4_notwd_first ( pt, pu, t);

        }
      if(ii==jj)
        {
          radix_4_dif_square_dit_blockE ( d , ii);
        }
    }
}


void radixmm_4_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                    y_size_t jj)
{
  Y__M128D  t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i;
  Y__M128D tw1r, tw1i, tw2r, tw2i, tw3r, tw3i;
  Y__M128D MM_1_0, MM_0_5, MM_AUX, gkr, gki;
  y_size_t nb, ni, bi, bj, j, i, li, lj, inc=Y_POWERS[3];

  y_ptr px, py, pt, pu;

  ASSERT_ALIGNED_DOUBLE();


  if(Y_PLAN[Y_NRADICES-1] != 4)
    return;
  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 4, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 4
  */

  /* set some constants */
  set_sse2_4_constants;

  nb=bs>>2;
  bi=(bs*ii)>>2;
  bj=((bs*jj)>>2)+nb-1;

  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      sse2_get_twiddle_factors_4_last (tw, i, j);

      Y_MM_UNPACKLO_PD (gkr, tw1r, tw1i);
      Y_MM_COMPLEX_MUL (gki, gkr, MM_AUX);
      Y_MM_INTER_PAIR (gkr, gki);

      li = i << 2;
      lj = j << 2;

      sse2_radix_4_twd_last (t, tw, pt, pu, d, li, lj);

      sse2_square_nested (t0r, t3i, t2r, t1i, gkr, gki);
      sse2_square_nested_1_4 (t1r, t2i, t3r, t0i, gkr, gki);

      sse2_radix_4_notwd_first ( pt, pu, t);
    }
  if((i==j) && (ni < nb))
    {
      radix_4_dif_square_dit_blockE ( d , i);
    }
}

#else /* defined Y_USE_SSE2 */

/* This line avoids errors on compilers expecting something to do in a file */
void radix_4_square_void( void )
{
  printf("The routine radix_4_square_void should never be called\n");
  exit(EXIT_FAILURE);
}

#endif /* Y_USE_SSE2 */
/*$Id$*/




