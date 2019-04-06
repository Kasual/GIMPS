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
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1

#define get_pad_2 \
  pd0= d;\
  pd1= d + addr((pad<<1));


#define radix_2_twd_pp(_j)\
      cplx_load_muladdsub_pp(t0,t1,a,_j,pd0,pd1);\
      cplx_local_to_data_p(pd1,t1);\
      cplx_local_to_data_p(pd0,t0);

#define radix_2_notwd_pp(_j)\
      cplx_load_addsub_pp(t0,t1,_j,pd0,pd1);\
      cplx_local_to_data_p(pd1,t1);\
      cplx_local_to_data_p(pd0,t0);


/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw[] = Scrambled Twidle factors. Not padded.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
 
   NOTE: All data lengths and pads are in size_of_complex units.
*************************************************************************/
void radix_2_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,ar,ai;
  y_ptr px;
  y_size_t i,j,jj,nc,nl;
  y_ptr pd0,pd1;


#ifndef NDEBUG1

  printf(" Radix_2_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pad_2;
  if(n0==0)
    {
      for(j=0;j<pad;j++)
        {
          jj=addr(j<<1);
          radix_2_notwd_pp(jj);
        }
      nc=pad<<1;
      if(n==nc)
        return;
      px=tw+2;
      for(i=nc;i<n;i+=nc,px+=2)
        {
          ar=(*px);
          ai=(*(px+1));
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_2_twd_pp(jj);
            }
        }
    }
  else
    {
      nc=pad<<1;
      nl=n0+n;
      px=tw;
      for(i=n0;i<nl;i+=nc,px+=2)
        {
          ar=(*px);
          ai=(*(px+1));
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_2_twd_pp(jj);
            }
        }
    }
}

/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r NO TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw = NULL in this  no  twiddle version.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radix_2_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_size_t i,j,nc,jj;
  y_ptr pd0,pd1;


#ifndef NDEBUG1

  printf(" Radix_2_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pad_2;
  if(tw != NULL)
    return;
  nc=pad<<1;
  for(i=n0;i<(n+n0);i+=nc)
    {
      for(j=i;j<(pad+i);j++)
        {
          jj=addr(j<<1);
          radix_2_notwd_pp(jj);
        }
    }
}


/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw[] = Twidle Backward factors. Not scrambled.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/
void radix_2_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,ar,ai;
  y_size_t i,j,nc,jj;
  y_ptr pd0,pd1;
  y_ptr px;


#ifndef NDEBUG1

  printf(" Radix_2_dit n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pad_2;
  nc=pad<<1;
  for(i=n0;i<(n+n0);i+=nc)
    {
      jj=addr(i<<1);
      radix_2_notwd_pp(jj);
      for(j=i+1,px=tw;j<(i+pad);j++,px+=2)
        {
          ar=(*px);
          ai=(*(px+1));
          jj=addr(j<<1);
          radix_2_twd_pp(jj);
        }
    }
}

/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r. NO TWIDDLE 
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw = NULL in this no twiddle version.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/
void radix_2_dit_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i;
  y_size_t i,j,jj,nc;
  y_ptr pd0,pd1;

#ifndef NDEBUG1

  printf(" Radix_2_dit_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pad_2;
  if(tw != NULL)
    return;
  nc=pad<<1;
  for(i=n0;i<(n+n0);i+=nc)
    {
      for(j=i;j<(i+pad);j++)
        {
          jj=addr(j<<1);
          radix_2_notwd_pp(jj);
        }
    }
}
/*$Id$*/










