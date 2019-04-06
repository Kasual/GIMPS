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
#include "fft7.h"
#define NDEBUG1


/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* asign the correct values to every pad */

#define get_pads         \
  bigpad=(pad<<3) - pad; \
  jj=pad;                \
  pd0= d;                \
  pd1= d + addr((jj<<1));\
  jj+=pad;               \
  pd2= d + addr((jj<<1));\
  jj+=pad;               \
  pd3= d + addr((jj<<1));\
  jj+=pad;               \
  pd4= d + addr((jj<<1));\
  jj+=pad;               \
  pd5= d + addr((jj<<1));\
  jj+=pad;               \
  pd6= d + addr((jj<<1));


#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

/* get the basic twiddle from memory and computes its powers */
# include "yminimum.h"
# define get_twiddle_factors_7          \
          cplx_trig_min_load(tw1,px);   \
          cplx_trig_min_square(tw2,tw1);\
          cplx_trig_min_square(tw4,tw2);\
          cplx_divmul(tw3,tw5,tw4,tw1); \
          cplx_trig_min_square(tw6,tw3);



#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 12

# define get_twiddle_factors_7         \
        tw1r= px[0]; tw1i= px[1];      \
        tw2r= px[2]; tw2i= px[3];      \
        tw3r= px[4]; tw3i= px[5];      \
        tw4r= px[6]; tw4i= px[7];      \
        tw5r= px[8]; tw5i= px[9];      \
        tw6r= px[10]; tw6i= px[11];    \
        prefetch_data_trig(px, Y_STEP);\
        prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);\


#else

# define get_twiddle_factors_7              \
	  tw1r= px[0];                      \
	  tw1i= px[1];                      \
	  tw2r= (tw1r+tw1i)*(tw1r-tw1i);    \
	  tw2i= 2.0*tw1r*tw1i;              \
	  tw4r= px[2];                      \
	  tw4i= px[3];                      \
	  tw6r= px[4];                      \
	  tw6i= px[5];                      \
	  cplx_divmul(tw3,tw5,tw4,tw1);     \
          prefetch_data(px,2);              \
          prefetch_data(px,2 + Y_CACHE_LINE);



# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6
#endif


/*
	Get data from memory, mul by twiddle trig factors,
	make a complex 7-length FFT and store the data to memory 
*/

# define radix_7_twd_pp(_S,_j)                          \
	  cplx_data_to_local_pp(t0,_j,pd0);             \
	  cplx_load_mulmul_pp(t1,t2,tw1,tw2,_j,pd1,pd2);\
	  cplx_load_mulmul_pp(t3,t4,tw3,tw4,_j,pd3,pd4);\
	  cplx_load_mulmul_pp(t5,t6,tw5,tw6,_j,pd5,pd6);\
	  cplx_fft7_store_p(_j,pd0,pd1,pd2,pd3,pd4,pd5,pd6,_S,t0,t1,t2,t3,t4,t5,t6);

/*
	Get data from memory, make a complex 7-length FFT and
	store the data to memory 
*/
# define radix_7_notwd_pp(_S,_j)           \
	  cplx_data_to_local_pp(t0,_j,pd0);\
	  cplx_data_to_local_pp(t1,_j,pd1);\
	  cplx_data_to_local_pp(t2,_j,pd2);\
	  cplx_data_to_local_pp(t3,_j,pd3);\
	  cplx_data_to_local_pp(t4,_j,pd4);\
	  cplx_data_to_local_pp(t5,_j,pd5);\
	  cplx_data_to_local_pp(t6,_j,pd6);\
	  cplx_fft7_store_p(_j,pd0,pd1,pd2,pd3,pd4,pd5,pd6,_S,t0,t1,t2,t3,t4,t5,t6);

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
void radix_7_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,
  tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6;
  y_size_t i,j,jj,nc;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();
#ifndef NDEBUG1

  printf(" Radix_7_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  get_pads;
  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j++)
        {
          jj=addr(j<<1);
          radix_7_notwd_pp(F,jj);
        }
      if(bigpad==n)
        return;
      px=tw+Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_7;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_7_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_7;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_7_twd_pp(F,jj);
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
void radix_7_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6;
  y_size_t nc,i,j,jj;
  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_7_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  if(tw != NULL)
    return;
  get_pads;
  nc=n+n0;
  for(i=n0;i<nc;i+=bigpad)
    {
      for(j=i;j<pad+i;j++)
        {
          /*BEGIN_TIME;*/
          jj=addr(j<<1);
          radix_7_notwd_pp(F,jj);
          /*END_TIME;*/
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
void radix_7_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,
  tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6;
  y_size_t i,j,jj,nc;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_7_dit n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  nc=n0+n;
  get_pads;

  for(i=n0;i<nc;i+=bigpad)
    {
      jj=addr(i<<1);
      radix_7_notwd_pp(B,jj);
      for(j=i+1,px=tw;j<(pad+i);j++,px+=Y_STEP)
        {
          get_twiddle_factors_7;
          jj=addr(j<<1);
          radix_7_twd_pp(B,jj);
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
void radix_7_dit_notw(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i;
  y_size_t bigpad;
  y_size_t pd0,pd1,pd2,pd3,pd4,pd5,pd6;
  y_size_t nc,i,j,jj;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_7_dit_notw n=%i n0=%i pad=%i \n",n,n0,pad);
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
          radix_7_notwd_pp(B,jj);
        }
    }
}
#endif
/*$Id$*/
