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
   It includes routines for a 2_last_radix reduction 
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1


#define get_twiddle_factors_2_last(_j) \
	  px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
	  tw1r= px[0];\
	  tw1i= px[1];\



#define radix_2_twd_last(_z,_d,_j)\
	  cplx_data_to_local( _z##0 ,_d, _j);\
	  cplx_data_to_local( _z##1 , _d,(_j + 1));\
	  cplx_muladdsub( _z##0 , _z##1 , tw1);\



#define radix_2_notwd_first(_d,_z,_j)\
	  cplx_addsub( _z##0 , _z##1 );\
	  cplx_local_to_data( _d, _j , _z##0 );\
	  cplx_local_to_data( _d,(_j + 1), _z##1 );\



void radix_2_dif_square_dit(y_ptr d, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_limb_t u0r,u0i,u1r,u1i;
  y_limb_t tw1r,tw1i;
  /*y_limb_t wkr,wki;*/
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[1];
  int nr;
  y_ptr px;

  ASSERT_ALIGNED_DOUBLE();


  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 4, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 2 
  */


  /* Cycle 0 block 0 */
  j=0;
  get_twiddle_factors_2_last(j);
  radix_2_twd_last(t,d,j);
  /* dyadic mul for k=0 */
  /*wkr=1.0; wki=0.0;*/
  square_nested_eq( t0 , 1.0);
  /*wkr=-1.0; wki=0.0;*/
  square_nested_eq( t1 , -1.0 );
  radix_2_notwd_first(d,t,j);
  /*printf(" %i --> 0 0 0 \n",nc);*/

  if(nc==0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-2-(int)nc);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>1;
      bj=(Y_LRIGHT[nr]>>1)-1;
      /*printf("%i %i -->\n",bi,bj);*/
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_2_last(jj);
          l=jj<<1;
          radix_2_twd_last(u,d,l);
          get_twiddle_factors_2_last(ii);
          l=ii<<1;
          radix_2_twd_last(t,d,l);
          square_nested(t0, u1, tw1r,tw1i);
          /*wkr=-tw1r; wki=-tw1i;*/
          square_nested(t1, u0, -tw1r, -tw1i);
          radix_2_notwd_first(d,t,l);
          l=jj<<1;
          radix_2_notwd_first(d,u,l);
          /*printf(" %i %i %i \n",nr,ii,jj);*/
        }
      if(ii==jj)
        {
          get_twiddle_factors_2_last(ii);
          l=ii<<1;
          radix_2_twd_last(t,d,l);
          square_nested(t0, t1, tw1r, tw1i);
          radix_2_notwd_first(d,t,l);
          /*printf(" %i %i %i \n",nr,ii,jj);*/
        }
    }
}


void radix_2_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_limb_t u0r,u0i,u1r,u1i;
  y_limb_t s0r,s0i,s1r,s1i;
  y_limb_t v0r,v0i,v1r,v1i;
  y_limb_t tw1r,tw1i;
  /*y_limb_t wkr,wki;*/
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[1];
  int nr;
  y_ptr px;

  ASSERT_ALIGNED_DOUBLE();


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
  get_twiddle_factors_2_last(j);
  radix_2_twd_last(s,d1,0);
  radix_2_twd_last(t,d2,0);
  /* dyadic mul for k=0 */
  /*wkr=1.0; wki= 0.0;*/
  conv_nested( s0, s0, t0, t0, 1.0, 0.0);
  /*wkr=-1.0;*/
  conv_nested( s1, s1, t1, t1, -1.0, 0.0);
  radix_2_notwd_first(d1,s,j);
  if(nc == 0)
    return;
  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>1;
      bj=(Y_LRIGHT[nr]>>1)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_2_last(jj);
          l=jj<<1;
          radix_2_twd_last(u,d1,l);
          radix_2_twd_last(v,d2,l);
          get_twiddle_factors_2_last(ii);
          l=ii<<1;
          radix_2_twd_last(s,d1,l);
          radix_2_twd_last(t,d2,l);
          conv_nested( s0, u1, t0, v1, tw1r, tw1i);
          /*wkr=-tw1r; wki= -tw1i;*/
          conv_nested( s1, u0, t1, v0, -tw1r, -tw1i);
          radix_2_notwd_first(d1,s,l);
          l=jj<<1;
          radix_2_notwd_first(d1,u,l);
        }
      if(ii==jj)
        {
          get_twiddle_factors_2_last(ii);
          l=ii<<1;
          radix_2_twd_last(s,d1,l);
          radix_2_twd_last(t,d2,l);
          conv_nested(s0, s1, t0, t1, tw1r, tw1i);
          radix_2_notwd_first(d1,s,l);
        }
    }
}

void radix_2_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                               y_size_t  jj)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_limb_t u0r,u0i,u1r,u1i;
  y_limb_t s0r,s0i,s1r,s1i;
  y_limb_t v0r,v0i,v1r,v1i;
  y_limb_t tw1r,tw1i;
  /*y_limb_t wkr,wki;*/
  y_size_t i,j,bi,bj,ni,l,nb,inc=Y_POWERS[1];
  y_ptr px;

  ASSERT_ALIGNED_DOUBLE();


  nb=bs>>1;
  bi=(bs*ii)>>1;
  bj=((bs*jj)>>1)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb) && (i<j);i++,j--,ni++)
    {
      get_twiddle_factors_2_last(j);
      l=j<<1;
      radix_2_twd_last(u,d1,l);
      radix_2_twd_last(v,d2,l);
      get_twiddle_factors_2_last(i);
      l=i<<1;
      radix_2_twd_last(s,d1,l);
      radix_2_twd_last(t,d2,l);
      conv_nested( s0, u1, t0, v1, tw1r, tw1i);
      /*wkr=-tw1r; wki= -tw1i;*/
      conv_nested( s1, u0, t1, v0, -tw1r, -tw1i);
      radix_2_notwd_first(d1,s,l);
      l=j<<1;
      radix_2_notwd_first(d1,u,l);
    }
  if((i==j)&&(ni<nb))
    {
      get_twiddle_factors_2_last(i);
      l=i<<1;
      radix_2_twd_last(s,d1,l);
      radix_2_twd_last(t,d2,l);
      conv_nested(s0, s1, t0, t1, tw1r, tw1i);
      radix_2_notwd_first(d1,s,l);
    }
}

void radix_2_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                  y_size_t jj)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_limb_t u0r,u0i,u1r,u1i;
  y_limb_t tw1r,tw1i;
  /*y_limb_t wkr,wki;*/
  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[1];
  y_ptr px;

  ASSERT_ALIGNED_DOUBLE();


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

  nb=bs>>1;
  bi=(bs*ii)>>1;
  bj=((bs*jj)>>1)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      get_twiddle_factors_2_last(j);
      l=j<<1;
      radix_2_twd_last(u,d,l);
      get_twiddle_factors_2_last(i);
      l=i<<1;
      radix_2_twd_last(t,d,l);

      square_nested(t0, u1, tw1r, tw1i);
      square_nested_1_2(t1, u0, tw1r, tw1i);

      radix_2_notwd_first(d,t,l);
      l=j<<1;
      radix_2_notwd_first(d,u,l);
      /*printf("%i %i\n",i,j);*/
    }
  if((i==j)&&(ni<nb))
    {
      get_twiddle_factors_2_last(i);
      l=i<<1;
      radix_2_twd_last(t,d,l);
      square_nested(t0, t1, tw1r, tw1i);
      radix_2_notwd_first(d,t,l);
    }
}
/*$Id$*/


