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
  bigpad=(pad<<3)+(pad<<2);\
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
  pd5= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd6= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd7= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd8= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd9= d + addr((jj<<1));  \
  jj+=pad;                 \
  pd10= d + addr((jj<<1)); \
  jj+=pad;                 \
  pd11= d + addr((jj<<1));


#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2


/* get the basic twiddle from memory and computes its powers */
# include "yminimum.h"
# define get_twiddle_factors_12          \
          cplx_trig_min_load(tw1,px);    \
          cplx_trig_min_square(tw2,tw1); \
          cplx_trig_min_square(tw4,tw2); \
          cplx_trig_min_square(tw8,tw4); \
          cplx_divmul(tw3,tw5,tw4,tw1);  \
          cplx_divmul(tw7,tw9,tw8,tw1);  \
          cplx_divmul(tw6,tw10,tw8,tw2); \
          cplx_trig_min_mul(tw11,tw8,tw3);

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 22

# define get_twiddle_factors_12        \
        tw1r= px[0]; tw1i= px[1];      \
        tw2r= px[2]; tw2i= px[3];      \
        tw3r= px[4]; tw3i= px[5];      \
        tw4r= px[6]; tw4i= px[7];      \
        tw5r= px[8]; tw5i= px[9];      \
        tw6r= px[10]; tw6i= px[11];    \
        tw7r= px[12]; tw7i= px[13];    \
        tw8r= px[14]; tw8i= px[15];    \
        tw9r= px[16]; tw9i= px[17];    \
        tw10r= px[18]; tw10i= px[19];    \
        tw11r= px[20]; tw11i= px[21];    \
        prefetch_data_trig(px, Y_STEP);\
        prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);\
        prefetch_data_trig(px, Y_STEP + 2*Y_CACHE_LINE);

#else

# define get_twiddle_factors_12             \
          tw1r= px[0];                      \
          tw1i= px[1];                      \
	  tw2r= px[2];                      \
	  tw2i= px[3];                      \
	  tw4r= (tw2r +tw2i)*(tw2r - tw2i); \
	  tw4i= 2.0*tw2r*tw2i;              \
	  tw3r= px[4];                      \
	  tw3i= px[5];                      \
	  tw8r= px[6];                      \
	  tw8i= px[7];                      \
	  cplx_divmul(tw5,tw11,tw8,tw3);    \
	  cplx_divmul(tw6,tw10,tw8,tw2);    \
	  cplx_divmul(tw7,tw9,tw8,tw1);     \
          prefetch_data(px,4);              \
          prefetch_data(px,4 + Y_CACHE_LINE);



# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 8
#endif


/*
	Get data from memory, mul by twiddle trig factors,
	make a complex 12-length FFT and store the data to memory 
*/
# define radix_12_twd_pp(_S,_j)                                \
	  cplx_data_to_local_pp(t0,_j,pd0);                    \
	  cplx_load_mulmul_pp(t1,t2,tw4,tw8,_j,pd4,pd8);       \
	  cplx_fft3(_S,t0,t1,t2);                              \
	  cplx_load_mul_pp(t3,tw2,_j,pd2);                     \
	  cplx_load_mulmul_pp(t4,t5,tw6,tw10,_j,pd6,pd10);     \
	  cplx_fft3(_S,t3,t4,t5);                              \
	  cplx_addsub(t0,t3);                                  \
	  cplx_muladdsub(t1,t4, _S##_1_6);                     \
	  cplx_muladdsub(t2,t5, _S##_1_3);                     \
	  cplx_load_mul_pp(t6,tw1,_j,pd1);                     \
	  cplx_load_mulmul_pp(t7,t8,tw5,tw9,_j,pd5,pd9);       \
	  cplx_fft3(_S,t6,t7,t8);                              \
	  cplx_load_mul_pp(t9,tw3,_j,pd3);                     \
	  cplx_load_mulmul_pp(t10,t11,tw7,tw11,_j,pd7,pd11);   \
	  cplx_fft3(_S,t9,t10,t11);                            \
	  cplx_addsub(t6,t9);                                  \
	  cplx_muladdsub(t7,t10, _S##_1_6);                    \
	  cplx_muladdsub(t8,t11, _S##_1_3);                    \
	                                                       \
	  cplx_addsub_store_p(_j,pd0,pd6,t0,t6);               \
	  cplx_muladdsub_store_p(_j,pd1,pd7,t1,t7, _S##_1_12); \
	  cplx_muladdsub_store_p(_j,pd2,pd8,t2,t8, _S##_1_6);  \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd3,pd9,t3,t9);\
	  cplx_muladdsub_store_p(_j,pd4,pd10,t4,t10, _S##_1_3);\
	  cplx_muladdsub_store_p(_j,pd5,pd11,t5,t11,_S##_5_12);


/*
	Get data from memory, make a complex 12-length FFT and
	store the data to memory 
*/
# define radix_12_notwd_pp(_S,_j)                              \
	  cplx_data_to_local_pp(t0,_j,pd0);                    \
	  cplx_data_to_local_pp(t1,_j,pd4);                    \
	  cplx_data_to_local_pp(t2,_j,pd8);                    \
	  cplx_fft3(_S,t0,t1,t2);                              \
	  cplx_data_to_local_pp(t3,_j,pd2);                    \
	  cplx_data_to_local_pp(t4,_j,pd6);                    \
	  cplx_data_to_local_pp(t5,_j,pd10);                   \
	  cplx_fft3(_S,t3,t4,t5);                              \
	  cplx_addsub(t0,t3);                                  \
	  cplx_muladdsub(t1,t4, _S##_1_6);                     \
	  cplx_muladdsub(t2,t5, _S##_1_3);                     \
	  cplx_data_to_local_pp(t6,_j,pd1);                    \
	  cplx_data_to_local_pp(t7,_j,pd5);                    \
	  cplx_data_to_local_pp(t8,_j,pd9);                    \
	  cplx_fft3(_S,t6,t7,t8);                              \
	  cplx_data_to_local_pp(t9,_j,pd3);                    \
	  cplx_data_to_local_pp(t10,_j,pd7);                   \
	  cplx_data_to_local_pp(t11,_j,pd11);                  \
	  cplx_fft3(_S,t9,t10,t11);                            \
	  cplx_addsub(t6,t9);                                  \
	  cplx_muladdsub(t7,t10, _S##_1_6);                    \
	  cplx_muladdsub(t8,t11, _S##_1_3);                    \
	                                                       \
	  cplx_addsub_store_p(_j,pd0,pd6,t0,t6);               \
	  cplx_muladdsub_store_p(_j,pd1,pd7,t1,t7, _S##_1_12); \
	  cplx_muladdsub_store_p(_j,pd2,pd8,t2,t8, _S##_1_6);  \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd3,pd9,t3,t9);\
	  cplx_muladdsub_store_p(_j,pd4,pd10,t4,t10, _S##_1_3);\
	  cplx_muladdsub_store_p(_j,pd5,pd11,t5,t11,_S##_5_12);


#if (Y_AVAL > 3) && defined(Y_MANY_REGISTERS)

# ifndef Y_SAVE
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

void radix_12_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,
  tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,
  tw11r,tw11i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  y_size_t i,j,jj,nc;

  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();
#  ifndef NDEBUG1

  printf(" Radix_12_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#  endif

  get_pads;
  /* If n0==0 first elements has twiddle factors equal to one */
  if(n0==0)
    {
      for(j=0;j<pad;j++)
        {
          /*BEGIN_TIME;*/
          jj=addr(j<<1);
          radix_12_notwd_pp(F,jj);
          /*END_TIME;*/
        }
      if(bigpad==n)
        return;
      px=tw+Y_STEP;
      for(i=bigpad;i<n;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_12;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_12_twd_pp(F,jj);
            }
        }
    }
  else
    {
      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad,px+=Y_STEP)
        {
          get_twiddle_factors_12;
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_12_twd_pp(F,jj);
            }
        }
    }
}
# endif

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

void radix_12_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                       y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  y_size_t nc,i,j,jj;
  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

#ifndef NDEBUG1

  printf(" Radix_12_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
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
          radix_12_notwd_pp(F,jj);
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

void radix_12_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,
  tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,
  tw11r,tw11i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
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
      radix_12_notwd_pp(B,jj);
      for(j=i+1,px=tw;j<(pad+i);j++,px+=Y_STEP)
        {
          get_twiddle_factors_12;
          jj=addr(j<<1);
          radix_12_twd_pp(B,jj);
        }
    }
}

# ifndef Y_SAVE
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

void radix_12_dit_notw(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0,
                       y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  y_size_t nc,i,j,jj;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radix_12_dit_notw n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  if (tw != NULL)
    return;
  nc=n0+n;

  get_pads;

  for(i=n0;i<nc;i+=bigpad)
    {
      for(j=i;j<(pad+i);j++)
        {
          jj=addr(j<<1);
          radix_12_notwd_pp(B,jj);
        }
    }
}
# endif

#else

void radix_12_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radix_12_dif , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_12.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

void radix_12_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                       y_size_t pad)
{
  printf("Called radix_12_dif_notwd , d=%p, tw=%p, n=%d, n=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_12.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

void radix_12_dit(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radix_12_dit , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_12.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

void radix_12_dit_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                       y_size_t pad)
{
  printf("Called radix_12_dit_notwd , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_12.c with -DY_AVAL=4 flag\n");
  exit(-1);
}

#endif
/*$Id$*/
















