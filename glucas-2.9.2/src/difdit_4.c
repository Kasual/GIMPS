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
   It includes routines for a 4_last_radix reduction 
 
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#ifdef Y_ITANIUM
# include "yia64.h"
# ifdef Y_LONG_MACROS
#  undef Y_LONG_MACROS
# endif
#endif
#define NDEBUG1
#ifdef Y_LONG_MACROS
# include "ygeneric4.h"
#endif

#ifdef Y_ITANIUM

#define init_pointers_4( _pd ,_d )         \
  _pd##0r = _d;      _pd##0i = _d + 1;     \
  _pd##1r = _d + 2;  _pd##1i = _d + 3;     \
  _pd##2r = _d + 4;  _pd##2i = _d + 5;     \
  _pd##3r = _d + 6;  _pd##3i = _d + 7


#define init_trig_pointers(_j,_k)               \
  pxtr = &(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
  pxdr = &(Y_TWDF[Y_NRADICES-2][(_k * inc)<<1]);\
  pxti = pxtr + 1;                              \
  pxdi = pxdr + 1

#endif

#if defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

# define get_twiddle_factors_4_last(_j)             \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          tw1r=px[0];      tw1i=px[1];              \
          tw2r=px[2];      tw2i=px[3];              \
          tw3r=px[4];      tw3i=px[5];              \
          prefetch_p_trig_nta(px);                  \
          prefetch_p_trig_nta(px + Y_CACHE_LINE);

# define get_twiddle_factors_4_last_down(_j)        \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          tw1r=px[0];      tw1i=px[1];              \
          tw2r=px[2];      tw2i=px[3];              \
          tw3r=px[4];      tw3i=px[5];              \
          prefetch_data_trig_nta(px, (-6));         \
          prefetch_data_trig_nta(px, (-6 + Y_CACHE_LINE));


#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2
# define get_twiddle_factors_4_last(_j)             \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw2r=(px[0]+px[1])*(px[0]-px[1]);         \
	  tw2i=2.0*px[0]*px[1];                     \
	  tw1r=px[0];                               \
	  tw1i=px[1];                               \
	  tw3r=tw1r*tw2r - tw1i*tw2i;               \
	  tw3i=tw1r*tw2i + tw1i*tw2r;               \
          prefetch_p_trig_nta(px);

# define get_twiddle_factors_4_last_down(_j)        \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw2r=(px[0]+px[1])*(px[0]-px[1]);         \
	  tw2i=2.0*px[0]*px[1];                     \
	  tw1r=px[0];                               \
	  tw1i=px[1];                               \
	  tw3r=tw1r*tw2r - tw1i*tw2i;               \
	  tw3i=tw1r*tw2i + tw1i*tw2r;               \
          prefetch_data_trig_nta(px,(-2));
#endif
#ifdef Y_LONG_MACROS

#define radix_4_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
          cplx_radix_4_pp_last_dif_passes(_z,tw,_pd);

#define radix_4_twd_last_down(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
          cplx_radix_4_pp_last_dif_passes_down(_z,tw,_pd);


#else
# define radix_4_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
	  cplx_load_muladdsub_pp( _z##0 , _z##1 ,tw2, 0, _pd, _pd + 4);\
	  cplx_load_mulmuladdsub_pp(_z##2 ,_z##3 ,tw1,tw3,0,_pd + 2,_pd + 6);\
	  \
          prefetch_p(_pd+4);\
	  cplx_addsub( _z##0 ,_z##2 );\
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3);\

#endif

#ifdef Y_LONG_MACROS

# define radix_4_notwd_first(_pd,_z,_j)\
          cplx_radix_4_pp_first_dit_passes(_z,_pd);

#else

# define radix_4_notwd_first(_pd,_z,_j)\
	  cplx_addsub( _z##0 , _z##2 );\
	  cplx_addsub( _z##1 , _z##3 );\
	  \
	  cplx_addsub_store_p(0, _pd,_pd + 4, _z##0 , _z##1 );\
	  cplx_mul_1_4_B_addsub_store_p(0,_pd + 2,_pd + 6, _z##2 , _z##3 );

#endif

void radix_4_dif_square_dit(y_ptr d, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
#ifdef Y_ITANIUM

  y_limb_t tp0r,tp0i,tp1r,tp1i,tp2r,tp2i,tp3r,tp3i;
  y_limb_t sp0r,sp0i,sp1r,sp1i,sp2r,sp2i,sp3r,sp3i;
  y_limb_t dw1r,dw1i,dw2r,dw2i,dw3r,dw3i,tpw1r,tpw1i;
#endif

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
#ifdef Y_ITANIUM

  y_ptr p0r,p0i,p1r,p1i,p2r,p2i,p3r,p3i;
  y_ptr pt0r,pt0i,pt1r,pt1i,pt2r,pt2i,pt3r,pt3i;
  y_ptr pd0r,pd0i,pd1r,pd1i,pd2r,pd2i,pd3r,pd3i;
  y_ptr ptl0r,ptl0i,ptl1r,ptl1i,ptl2r,ptl2i,ptl3r,ptl3i;
  y_ptr pdl0r,pdl0i,pdl1r,pdl1i,pdl2r,pdl2i,pdl3r,pdl3i;
  y_ptr pxtr,pxti,pxdr,pxdi;
#endif

  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[3];
  int nr;
  y_ptr px,pt;
#ifndef Y_ITANIUM

  y_ptr pu;
#endif

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


  /* Cycle 0 block 0 */
  j=0;
  get_twiddle_factors_4_last(j);
  radix_4_twd_last(t,pt,d,j);
  /* dyadic mul for k=0 */

  square_nested_eq( t0 , 1.0 );

  square_nested_1_4( t1 , t3 ,1.0, 0.0);

  square_nested_eq( t2 , -1.0 );
  radix_4_notwd_first(pt,t,j);
  /*printf(" %i --> 0 0 0 \n",nc);*/

  if(nc==0)
    return;

#ifdef Y_ITANIUM

  init_pointers_4(p,d);
#endif

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-2-(int)nc);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>2;
      bj=(Y_LRIGHT[nr]>>2)-1;
      /*printf("%i %i -->\n",bi,bj);*/
#ifdef Y_ITANIUM

      init_trig_pointers(bi,bj);
      last_dif_4_before_loop_preload(bi,bj);
#endif

      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
#ifndef Y_ITANIUM
          get_twiddle_factors_4_last_down(jj);
          l=jj<<2;
# ifdef Y_LONG_MACROS

          radix_4_twd_last_down(u,pu,d,l);
# else

          radix_4_twd_last(u,pu,d,l);
# endif

          get_twiddle_factors_4_last(ii);
          l=ii<<2;
          radix_4_twd_last(t,pt,d,l);

          square_nested(t0, u3, tw1r, tw1i);
          square_nested_1_4(t1, u2, tw1r, tw1i);
          square_nested_1_2(t2, u1, tw1r, tw1i);
          square_nested_3_4(t3, u0, tw1r, tw1i);

          radix_4_notwd_first(pt,t,l);
          l=jj<<2;
          radix_4_notwd_first(pu,u,l);

#else

          big_macro_last_dif_4_preload(ii,jj);
          square_nested_0_0_1_4(t0,u3,t1,u2,tw1r, tw1i);
          square_nested_1_2_3_4(t2,u1,t3,u0,tw1r, tw1i);
          big_macro_first_dit_4_preload();

#endif

        }
      if(ii==jj)
        {
          get_twiddle_factors_4_last(ii);
          l=ii<<2;
          radix_4_twd_last(t,pt,d,l);
          square_nested(t0, t3, tw1r, tw1i);
          square_nested_1_4(t1, t2, tw1r, tw1i);

          radix_4_notwd_first(pt,t,l);
        }
    }
}


void radix_4_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[3];
  int nr;
  y_ptr px,pt,pu,ps,pv;

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

  /* Cycle 0 block 0 */
  j=0;
  get_twiddle_factors_4_last(j);
  radix_4_twd_last(s,ps,d1,0);
  radix_4_twd_last(t,pt,d2,0);
  /* dyadic mul for k=0 */
  conv_nested( s0, s0, t0, t0, 1.0, 0.0 );
  conv_nested_1_4( s1, s3, t1, t3, 1.0, 0.0 );
  conv_nested_1_2( s2, s2, t2, t2, 1.0, 0.0 );

  radix_4_notwd_first(ps,s,j);
  if(nc == 0)
    return;
  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>2;
      bj=(Y_LRIGHT[nr]>>2)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_4_last_down(jj);
          l=jj<<2;
#ifdef Y_LONG_MACROS

          radix_4_twd_last_down(u,pu,d1,l);
          radix_4_twd_last_down(v,pv,d2,l);
#else

          radix_4_twd_last(u,pu,d1,l);
          radix_4_twd_last(v,pv,d2,l);
#endif

          get_twiddle_factors_4_last(ii);
          l=ii<<2;
          radix_4_twd_last(s,ps,d1,l);
          radix_4_twd_last(t,pt,d2,l);

          conv_nested( s0, u3, t0, v3, tw1r, tw1i);
          conv_nested_1_4( s1, u2, t1, v2, tw1r, tw1i);
          conv_nested_1_2( s2, u1, t2, v1, tw1r, tw1i);
          conv_nested_3_4( s3, u0, t3, v0, tw1r, tw1i);

          radix_4_notwd_first(ps,s,l);
          l=jj<<2;
          radix_4_notwd_first(pu,u,l);
        }
      if(ii==jj)
        {
          get_twiddle_factors_4_last(ii);
          l=ii<<2;
          radix_4_twd_last(s,ps,d1,l);
          radix_4_twd_last(t,pt,d2,l);
          conv_nested(s0, s3, t0, t3, tw1r, tw1i);
          conv_nested_1_4(s1, s2, t1, t2, tw1r, tw1i);
          radix_4_notwd_first(ps,s,l);
        }
    }
}

void radix_4_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                               y_size_t  jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  y_size_t i,j,bi,bj,ni,l,nb,inc=Y_POWERS[3];
  y_ptr px,pu,ps,pt,pv;

  ASSERT_ALIGNED_DOUBLE();


  if(Y_PLAN[Y_NRADICES-1] != 4)
    return;

  nb=bs>>2;
  bi=(bs*ii)>>2;
  bj=((bs*jj)>>2)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb) && (i<j);i++,j--,ni++)
    {
      get_twiddle_factors_4_last_down(j);
      l=j<<2;
#ifdef Y_LONG_MACROS

      radix_4_twd_last_down(u,pu,d1,l);
      radix_4_twd_last_down(v,pv,d2,l);
#else

      radix_4_twd_last(u,pu,d1,l);
      radix_4_twd_last(v,pv,d2,l);
#endif

      get_twiddle_factors_4_last(i);
      l=i<<2;
      radix_4_twd_last(s,ps,d1,l);
      radix_4_twd_last(t,pt,d2,l);

      conv_nested( s0, u3, t0, v3, tw1r, tw1i);
      conv_nested_1_4( s1, u2, t1, v2, tw1r, tw1i);
      conv_nested_1_2( s2, u1, t2, v1, tw1r, tw1i);
      conv_nested_3_4( s3, u0, t3, v0, tw1r, tw1i);

      radix_4_notwd_first(ps,s,l);
      l=j<<2;
      radix_4_notwd_first(pu,u,l);
    }
  if((i==j)&&(ni<nb))
    {
      get_twiddle_factors_4_last(i);
      l=i<<2;
      radix_4_twd_last(s,ps,d1,l);
      radix_4_twd_last(t,pt,d2,l);
      conv_nested(s0, s3, t0, t3, tw1r, tw1i);
      conv_nested_1_4(s1, s2, t1, t2, tw1r, tw1i);
      radix_4_notwd_first(ps,s,l);
    }
}

void radix_4_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                  y_size_t jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
#ifdef Y_ITANIUM

  y_limb_t tp0r,tp0i,tp1r,tp1i,tp2r,tp2i,tp3r,tp3i;
  y_limb_t sp0r,sp0i,sp1r,sp1i,sp2r,sp2i,sp3r,sp3i;
  y_limb_t dw1r,dw1i,dw2r,dw2i,dw3r,dw3i,tpw1r,tpw1i;
#endif

  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
#ifdef Y_ITANIUM

  y_ptr p0r,p0i,p1r,p1i,p2r,p2i,p3r,p3i;
  y_ptr pt0r,pt0i,pt1r,pt1i,pt2r,pt2i,pt3r,pt3i;
  y_ptr pd0r,pd0i,pd1r,pd1i,pd2r,pd2i,pd3r,pd3i;
  y_ptr ptl0r,ptl0i,ptl1r,ptl1i,ptl2r,ptl2i,ptl3r,ptl3i;
  y_ptr pdl0r,pdl0i,pdl1r,pdl1i,pdl2r,pdl2i,pdl3r,pdl3i;
  y_ptr pxtr,pxti,pxdr,pxdi;
#endif

  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[3];
  y_ptr px,pt;
#ifndef Y_ITANIUM

  y_ptr pu;
#endif

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

  nb=bs>>2;
  bi=(bs*ii)>>2;
  bj=((bs*jj)>>2)+nb-1;

#ifdef Y_ITANIUM

  init_pointers_4(p,d);
  init_trig_pointers(bi,bj);
  last_dif_4_before_loop_preload(bi,bj);
#endif

  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
#ifndef Y_ITANIUM
      get_twiddle_factors_4_last_down(j);
      l=j<<2;
# ifdef Y_LONG_MACROS

      radix_4_twd_last_down(u,pu,d,l);
# else

      radix_4_twd_last(u,pu,d,l);
# endif

      get_twiddle_factors_4_last(i);
      l=i<<2;
      radix_4_twd_last(t,pt,d,l);
      square_nested(t0, u3, tw1r, tw1i);
      square_nested_1_4(t1, u2, tw1r, tw1i);
      square_nested_1_2(t2, u1, tw1r, tw1i);
      square_nested_3_4(t3, u0, tw1r, tw1i);
      radix_4_notwd_first(pt,t,l);
      l=j<<2;
      radix_4_notwd_first(pu,u,l);
#else

      big_macro_last_dif_4_preload(i,j);
      square_nested_0_0_1_4(t0,u3,t1,u2,tw1r, tw1i);
      square_nested_1_2_3_4(t2,u1,t3,u0,tw1r, tw1i);
      big_macro_first_dit_4_preload();
#endif

    }
  if((i==j) && (ni < nb))
    {
      get_twiddle_factors_4_last(i);
      l=i<<2;
      radix_4_twd_last(t,pt,d,l);

      square_nested(t0, t3, tw1r, tw1i);
      square_nested_1_4(t1, t2, tw1r, tw1i);

      radix_4_notwd_first(pt,t,l);
    }
}
/*$Id$*/




