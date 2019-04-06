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
#include "fft3.h"
#define NDEBUG1


/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* asign the correct values to every pad */

#define get_pads           \
  bigpad=(pad<<2)+(pad<<1);\
  jj=pad;                  \
  pd0= d;                  \
  pd1= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd2= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd3= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd4= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd5= d + addr((jj<<1));



#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

/* get the basic twiddle from memory and computes its powers */
# include "yminimum.h"

# define get_twiddle_factors_6          \
          cplx_trig_min_load(tw1,px);   \
          cplx_trig_min_square(tw2,tw1);\
          cplx_trig_min_square(tw4,tw2);\
          cplx_divmul(tw3,tw5,tw4,tw1);

#elif defined(Y_MAXIMUM)

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 10

# define get_twiddle_factors_6         \
        tw1r= px[0]; tw1i= px[1];      \
        tw2r= px[2]; tw2i= px[3];      \
        tw3r= px[4]; tw3i= px[5];      \
        tw4r= px[6]; tw4i= px[7];      \
        tw5r= px[8]; tw5i= px[9];      \
        prefetch_data_trig(px, Y_STEP);\
        prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);\



#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 4

# define get_twiddle_factors_6          \
	  tw1r= px[0];                  \
	  tw1i= px[1];                  \
	  tw2r= (tw1r+tw1i)*(tw1r-tw1i);\
	  tw2i= 2.0*tw1r*tw1i;          \
	  tw4r= px[2];                  \
	  tw4i= px[3];                  \
	  cplx_divmul(tw3,tw5,tw4,tw1); \
          prefetch_p_trig(px);


#endif



/*
	Get data from memory, mul by twiddle trig factors,
	make a complex 6-length FFT and store the data to memory 
*/

# define radix_6_twd_pp(_S,_j)                               \
	  cplx_data_to_local_pp(t0,_j,pd0);                  \
	  cplx_load_mulmul_pp(t1,t2,tw2,tw4,_j,pd2,pd4);     \
	  cplx_fft3(_S,t0,t1,t2);                            \
	  cplx_load_mul_pp(t3,tw1,_j,pd1);                   \
	  cplx_load_mulmul_pp(t4,t5,tw3,tw5,_j,pd3,pd5);     \
	  cplx_fft3(_S,t3,t4,t5);                            \
	                                                     \
	  cplx_addsub_store_p(_j,pd0,pd3,t0,t3);             \
	  cplx_muladdsub_store_p(_j,pd1,pd4,t1,t4, _S##_1_6);\
	  cplx_muladdsub_store_p(_j,pd2,pd5,t2,t5, _S##_1_3);


/*
	Get data from memory, make a complex 6-length FFT and
	store the data to memory 
*/

# define radix_6_notwd_pp(_S,_j)                             \
	  cplx_data_to_local_pp(t0,_j,pd0);                  \
	  cplx_data_to_local_pp(t1,_j,pd2);                  \
	  cplx_data_to_local_pp(t2,_j,pd4);                  \
	  cplx_fft3(_S,t0,t1,t2);                            \
	  cplx_data_to_local_pp(t3,_j,pd1);                  \
	  cplx_data_to_local_pp(t4,_j,pd3);                  \
	  cplx_data_to_local_pp(t5,_j,pd5);                  \
	  cplx_fft3(_S,t3,t4,t5);                            \
	                                                     \
	  cplx_addsub_store_p(_j,pd0,pd3,t0,t3);             \
	  cplx_muladdsub_store_p(_j,pd1,pd4,t1,t4, _S##_1_6);\
	  cplx_muladdsub_store_p(_j,pd2,pd5,t2,t5, _S##_1_3);


#ifndef Y_SAVE
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
*************************************************************************/

void radix_6_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,tw1r,tw1i,
  tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5;
  y_size_t i,j,jj,nc;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();
#ifndef NDEBUG1

  printf(" Radix_6_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pads;
  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j++)
        {
          /*BEGIN_TIME;*/
          jj=addr(j<<1);
          radix_6_notwd_pp(F,jj);
          /*END_TIME;*/
        }
      if(bigpad==n)
        return;
      px=tw+Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_6;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_6_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_6;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_6_twd_pp(F,jj);
            }
        }
    }
}
#endif

/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r  WITHOUT TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = NULL in this no twiddle version. 
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radix_6_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5;
  y_size_t nc,i,j,jj;
  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_6_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  if(tw != NULL)
    return;
  get_pads;
  nc=n+n0;
  for(i=n0;i<nc;i+=bigpad)
    {
      for(j=i;j<pad+i;j++)
        {
          jj=addr(j<<1);
          radix_6_notwd_pp(F,jj);
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
*************************************************************************/

void radix_6_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,tw1r,tw1i,
  tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5;
  y_size_t i,j,jj,nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_12_dit n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  nc=n0+n;
  get_pads;

  for(i=n0;i<nc;i+=bigpad)
    {
      jj=addr(i<<1);
      radix_6_notwd_pp(B,jj);
      for(j=i+1,px=tw;j<(pad+i);j++,px+=Y_STEP)
        {
          get_twiddle_factors_6;
          jj=addr(j<<1);
          radix_6_twd_pp(B,jj);
        }
    }
}
#ifndef Y_SAVE
/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r NO TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
	   element of the whole data array. 
   tw[] = NULL in this no twiddle version. 
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radix_6_dit_notw(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5;
  y_size_t nc,i,j,jj;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_6_dit_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  if (tw != NULL)
    return;
  nc=n0+n;

  get_pads;

  for(i=n0;i<nc;i+=bigpad)
    {
      for(j=i;j<(pad+i);j++)
        {
          jj=addr(j<<1);
          radix_6_notwd_pp(B,jj);
        }
    }
}
#endif
/*$Id$*/
















