/*
$Id$
*/
/*
   This file is a part of 
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
#define NDEBUG1
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"

/*#define SEE_PLAN*/

#if defined(Y_NEW_PADDING)
y_size_t Y_SHIFT0, Y_SHIFT1, Y_SHIFT2, Y_SHIFT3, Y_SHIFT4;
#endif

/************  A lot of global definitions *********************
Most of them are set when calling to yea_init() procedure: 
Y_TWDF the array of pointers to twidle factors for forward transform. 
	Example, TWDF[i] points to twidle factors for pass i+1 on 
        forward transform.
Y_TWDB idem to Y_TWDF for backward transform.
 
F_x_yr   real(G^(-x/y)). The real part for radix x forward.
F_x_yi   The imaginary part of powers for forward
B_x_yr   real(G^(x/y)). The real part for backward transform.
B_x_yi   The imaginary part for backward transform.
 
Y_LENGTH the FFT length (in size_of_complex units). 
Y_Q, Y_K . Y_LENGTH= TP_Q*2^(TP_K)  Here, TP_Q=5,6,7 or 8 
Y_NRADICES number of radix reduction passes in the FFT.
Y_PLAN the plan for radix reduction. Y_PLAN[i] is the i+1 radix
	reduction on forward transform and Y_NRADICES-1-i on backward.
Y_LLEFT[i] the product of Y_PLAN elements less or equal than i.  
Y_LRIGHT[i] the product of Y_PLAN elelments great or equal than i.
**************************************************************************/
y_limb_t F_1_4r, F_1_4i, F_1_8r, F_1_8i, F_3_8r, F_3_8i, F_1_16r, F_1_16i,
F_3_16r, F_3_16i, F_5_16r, F_5_16i, F_7_16r, F_7_16i;
y_limb_t  F_1_32r, F_1_32i, F_3_32r, F_3_32i, F_5_32r, F_5_32i, F_7_32r,
F_7_32i, F_9_32r, F_9_32i, F_11_32r, F_11_32i, F_13_32r, F_13_32i, F_15_32r,
F_15_32i;
y_limb_t F_1_3r, F_1_3i, FM150, F_1_6r, F_1_6i, F_1_12r, F_1_12i, F_5_12r,
F_5_12i;
y_limb_t F_1_9r, F_1_9i, F_2_9r, F_2_9i, F_4_9r, F_4_9i;
y_limb_t F_1_5r, F_1_5i, F_2_5r, F_2_5i;
y_limb_t FN1_5r, FN1_5i, FNM125, FN2_5r, FN2_5i;
y_limb_t F_1_10r, F_1_10i, F_3_10r, F_3_10i;
y_limb_t F_1_7r, F_1_7i, F_2_7r, F_2_7i, F_3_7r, F_3_7i;
y_limb_t FN1_7r, FN1_7i, FN2_7r, FN2_7i, FN3_7r, FN3_7i, FN4_7r, FN4_7i;
y_limb_t F_1_14r, F_1_14i, F_3_14r, F_3_14i, F_5_14r, F_5_14i;
y_limb_t B_1_4r, B_1_4i, B_1_8r, B_1_8i, B_3_8r, B_3_8i, B_1_16r, B_1_16i,
B_3_16r, B_3_16i, B_5_16r, B_5_16i, B_7_16r, B_7_16i;
y_limb_t  B_1_32r, B_1_32i, B_3_32r, B_3_32i, B_5_32r, B_5_32i, B_7_32r,
B_7_32i, B_9_32r, B_9_32i, B_11_32r, B_11_32i, B_13_32r, B_13_32i, B_15_32r,
B_15_32i;
y_limb_t B_1_3r, B_1_3i, BM150, B_1_6r, B_1_6i, B_1_12r, B_1_12i, B_5_12r,
B_5_12i;
y_limb_t B_1_9r, B_1_9i, B_2_9r, B_2_9i, B_4_9r, B_4_9i;
y_limb_t B_1_5r, B_1_5i, B_2_5r, B_2_5i;
y_limb_t BN1_5r, BN1_5i, BNM125, BN2_5r, BN2_5i;
y_limb_t B_1_10r, B_1_10i, B_3_10r, B_3_10i;
y_limb_t B_1_7r, B_1_7i, B_2_7r, B_2_7i, B_3_7r, B_3_7i;
y_limb_t BN1_7r, BN1_7i, BN2_7r, BN2_7i, BN3_7r, BN3_7i, BN4_7r, BN4_7i;
y_limb_t B_1_14r, B_1_14i, B_3_14r, B_3_14i, B_5_14r, B_5_14i;
y_size_t *Y_PLAN = NULL, *Y_LRIGHT = NULL, *Y_LLEFT = NULL;
y_ptr *Y_TWDF = NULL, *Y_TWDB = NULL, Y_DYADIC = NULL, *Y_PF = NULL,
                                *Y_PB = NULL, Y_D=NULL, Y_WORK;
y_size_t Y_LENGTH = 0, Y_Q = 0, Y_K = 0, Y_NRADICES = 0, Y_OLDRADICES = 0;

#if defined(Y_USE_SSE2)
Y__M128D MM_1_0, MM_0_5, MM_0_0, MM_Y0LO, MM_Y0HI, MM_MTWO;

# if !defined(Y_SIMULATE_SSE2)
Y__M128D MM_YCHS, MM_YABS;
# endif /* !Y_SIMULATE_SSE2 */

Y__M128D MM_F_1_4r, MM_F_1_4i, MM_F_1_8r, MM_F_1_8i, MM_F_3_8r,
MM_F_3_8i, MM_F_1_16r, MM_F_1_16i, MM_F_3_16r, MM_F_3_16i, MM_F_5_16r,
MM_F_5_16i, MM_F_7_16r, MM_F_7_16i;

Y__M128D MM_B_1_4r, MM_B_1_4i, MM_B_1_8r, MM_B_1_8i, MM_B_3_8r,
MM_B_3_8i, MM_B_1_16r, MM_B_1_16i, MM_B_3_16r, MM_B_3_16i, MM_B_5_16r,
MM_B_5_16i, MM_B_7_16r, MM_B_7_16i;

Y__M128D MM_F_1_3r, MM_F_1_3i, MM_FM150, MM_F_1_6r, MM_F_1_6i,
MM_F_1_12r, MM_F_1_12i, MM_F_5_12r, MM_F_5_12i;

Y__M128D MM_F_1_9r, MM_F_1_9i, MM_F_2_9r, MM_F_2_9i, MM_F_4_9r,
MM_F_4_9i;

Y__M128D MM_B_1_3r, MM_B_1_3i, MM_BM150, MM_B_1_6r, MM_B_1_6i,
MM_B_1_12r, MM_B_1_12i, MM_B_5_12r, MM_B_5_12i;

Y__M128D MM_B_1_9r, MM_B_1_9i, MM_B_2_9r, MM_B_2_9i, MM_B_4_9r,
MM_B_4_9i;

Y__M128D MM_FNM125, MM_BNM125;

Y__M128D MM_F_1_5r, MM_F_1_5i, MM_F_2_5r, MM_F_2_5i;

Y__M128D MM_FN1_5r, MM_FN1_5i, MM_FN2_5r, MM_FN2_5i;

Y__M128D MM_F_1_10r, MM_F_1_10i, MM_F_3_10r, MM_F_3_10i;

Y__M128D MM_B_1_5r, MM_B_1_5i, MM_B_2_5r, MM_B_2_5i;

Y__M128D MM_BN1_5r, MM_BN1_5i, MM_BN2_5r, MM_BN2_5i;

Y__M128D MM_B_1_10r, MM_B_1_10i, MM_B_3_10r, MM_B_3_10i;

Y__M128D MM_F_1_7r, MM_F_1_7i, MM_F_2_7r, MM_F_2_7i, MM_F_3_7r,
MM_F_3_7i;

Y__M128D MM_FN1_7r, MM_FN1_7i, MM_FN2_7r, MM_FN2_7i, MM_FN3_7r,
MM_FN3_7i, MM_FN4_7r, MM_FN4_7i;

Y__M128D MM_F_1_14r, MM_F_1_14i, MM_F_3_14r, MM_F_3_14i, MM_F_5_14r, MM_F_5_14i;

Y__M128D MM_B_1_7r, MM_B_1_7i, MM_B_2_7r, MM_B_2_7i, MM_B_3_7r,
MM_B_3_7i;

Y__M128D MM_BN1_7r, MM_BN1_7i, MM_BN2_7r, MM_BN2_7i, MM_BN3_7r,
MM_BN3_7i, MM_BN4_7r, MM_BN4_7i;

Y__M128D MM_B_1_14r, MM_B_1_14i, MM_B_3_14r, MM_B_3_14i, MM_B_5_14r, MM_B_5_14i;

/* This macro is used to set the SSE2 version of most used double constant*/
# define SSE2_VERS( _var)              \
Y_MM_SET1_PD( MM_##_var##r , _var##r); \
Y_MM_SET1_PD( MM_##_var##i , _var##i);

#endif /* Y_USE_SSE2 */

/*number of trig factors powers computed for every radix */
#if defined(Y_MINIMUM)
y_size_t Y_POWERS[32]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#elif defined(Y_MAXIMUM)
y_size_t Y_POWERS[32]={1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                       16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                       30, 31};
#else
y_size_t Y_POWERS[32]={1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 1, 4, 1, 5, 1, 5,
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9};
#endif

/* This function performs the ith root of one, i.e w^i=1 mod p */
void y_root( y_ptr r, y_ptr i, y_size_t num, y_size_t den, int ifor)
{
  *r = cos ((BIG_DOUBLE)(num) / ((BIG_DOUBLE)(den)) * TWOPI);
  *i = sin ((BIG_DOUBLE)(num) / ((BIG_DOUBLE)(den)) * TWOPI);
  if(ifor)
    *i = (-(*i));
}

/*
  This routine inits the basic powers for radix r. The basic powers we need
  are:
  radix 2: none.
  radix 3: w^(1/3), w^(2/3);
  radix 4: w^(1/4).
  radix 5: w^(1/5),w^(2/5),w^(3/5),w^(4/5)
  radix 6: w^(1/6),w^(2/6)
  radix 7: w^(1/7),w^(2/7),w^(3/7),w^(4/7),w^(5/7),w^(6/7)
  radix 8: w^(1/8),w^(2/8),w^(3/8)
  radix 10: w^(1/10),w^(2/10),w^(3/10),w^(4/10)
  radix 12: w^(1/12),w^(2/12)...w^(5/12)
  radix 14: w^(1/14),w^(2/14)...w^(6/14)
  radix 16: w^(1/16),w^(2/16)...w^(7/16)
  radix 32: w^(1/32),w^(32/16)...w^(15/32)
*/
void y_initradix(void)
{
  y_limb_t re=0.0,im=0.0;

  FM150 = -1.5;
  BM150 = -1.5;
  FNM125 = -1.25;
  BNM125 = -1.25;

#if defined(Y_USE_SSE2)

  Y_MM_SET1_PD( MM_0_0, 0.0);
  Y_MM_SET1_PD( MM_0_5, 0.5);
  Y_MM_SET1_PD( MM_1_0, 1.0);
  Y_MM_SET1_PD( MM_MTWO, -2.0);
  Y_MM_SET1_PD( MM_FM150, -1.5);
  Y_MM_SET1_PD( MM_BM150, -1.5);
  Y_MM_SET1_PD( MM_FNM125, -1.25);
  Y_MM_SET1_PD( MM_BNM125, -1.25);

# if !defined(Y_SIMULATE_SSE2)
  /*
     an xor to MM_YCHS will change the sign 
     an and to MM_YABS will compute the abs value
     an and to MM_Y0HI will change d[1] part to 0.0
  */
  {
    Y__M128D aux;
    Y_MM_SUB_PD( aux, MM_0_0, MM_1_0);  /* aux = -1.0 */
    Y_MM_XOR_PD( MM_YCHS, MM_1_0, aux); /* sign bit to 1 */
    Y_MM_CMPGT_PD( MM_YABS, MM_1_0, MM_0_0); /* aux = 0xff...ff */
    Y_MM_UNPACKHI_PD( MM_Y0HI, MM_YABS, MM_0_0); /* y0hi = ( 0xffff..., 0x0000) */
    Y_MM_UNPACKHI_PD( MM_Y0LO, MM_0_0, MM_YABS); /* y0hi = ( 0xffff..., 0x0000) */
    Y_MM_XOR_PD( MM_YABS, MM_YABS, MM_YCHS); /* y_chs = 0x7ff..f */
  }
# endif
#endif

  /* radix 3*/
  y_root(&re, &im, 1, 3, 0);
  F_1_3r = re;
  F_1_3i = -im;
  B_1_3r = re;
  B_1_3i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_3);
  SSE2_VERS (B_1_3);
#endif

  /* radix 6*/
  y_root(&re, &im, 1, 6, 0);
  F_1_6r = re;
  F_1_6i = -im;
  B_1_6r = re;
  B_1_6i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_6);
  SSE2_VERS (B_1_6);
#endif

  /* radix 9*/
  y_root(&re, &im, 1, 9, 0);
  F_1_9r = re;
  F_1_9i = -im;
  B_1_9r = re;
  B_1_9i = im;
  y_root(&re, &im, 2, 9, 0);
  F_2_9r = re;
  F_2_9i = -im;
  B_2_9r = re;
  B_2_9i = im;
  y_root(&re, &im, 4, 9, 0);
  F_4_9r = re;
  F_4_9i = -im;
  B_4_9r = re;
  B_4_9i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_9);
  SSE2_VERS (B_1_9);
  SSE2_VERS (F_2_9);
  SSE2_VERS (B_2_9);
  SSE2_VERS (F_4_9);
  SSE2_VERS (B_4_9);
#endif


  /* radix 12*/
  y_root(&re, &im, 1, 12,0);
  F_1_12r = re;
  F_1_12i = -im;
  B_1_12r = re;
  B_1_12i = im;
  y_root(&re, &im, 5, 12, 0);
  F_5_12r = re;
  F_5_12i = -im;
  B_5_12r = re;
  B_5_12i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_12);
  SSE2_VERS (B_1_12);
  SSE2_VERS (F_5_12);
  SSE2_VERS (B_5_12);
#endif

  /* radix 4*/
  F_1_4r = 0.0;
  F_1_4i = -1.0;
  B_1_4r = 0.0;
  B_1_4i = 1.0;

  /* radix 8 */
  F_1_8r = SQRTHALF;
  F_1_8i = -SQRTHALF;
  F_3_8r = -SQRTHALF;
  F_3_8i = -SQRTHALF;
  B_1_8r = SQRTHALF;
  B_1_8i = SQRTHALF;
  B_3_8r = -SQRTHALF;
  B_3_8i = SQRTHALF;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_8);
  SSE2_VERS (B_1_8);
  SSE2_VERS (F_3_8);
  SSE2_VERS (B_3_8);
#endif

  /* radix 16 */
  y_root(&re, &im, 1, 16, 0);
  F_1_16r = re;
  F_1_16i = -im;
  B_1_16r = re;
  B_1_16i = im;
  y_root(&re, &im, 3, 16, 0);
  F_3_16r = re;
  F_3_16i = -im;
  B_3_16r = re;
  B_3_16i = im;
  y_root(&re, &im, 5, 16, 0);
  F_5_16r = re;
  F_5_16i = -im;
  B_5_16r = re;
  B_5_16i = im;
  y_root(&re, &im, 7, 16, 0);
  F_7_16r = re;
  F_7_16i = -im;
  B_7_16r = re;
  B_7_16i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_16);
  SSE2_VERS (B_1_16);
  SSE2_VERS (F_3_16);
  SSE2_VERS (B_3_16);
  SSE2_VERS (F_5_16);
  SSE2_VERS (B_5_16);
  SSE2_VERS (F_7_16);
  SSE2_VERS (B_7_16);
#endif


  /* radix 32 */
  y_root(&re, &im, 1, 32, 0);
  F_1_32r = re;
  F_1_32i = -im;
  B_1_32r = re;
  B_1_32i = im;
  y_root(&re, &im, 3, 32, 0);
  F_3_32r = re;
  F_3_32i = -im;
  B_3_32r = re;
  B_3_32i = im;
  y_root(&re, &im, 5, 32, 0);
  F_5_32r = re;
  F_5_32i = -im;
  B_5_32r = re;
  B_5_32i = im;
  y_root(&re, &im, 7, 32, 0);
  F_7_32r = re;
  F_7_32i = -im;
  B_7_32r = re;
  B_7_32i = im;
  y_root(&re, &im, 9, 32, 0);
  F_9_32r = re;
  F_9_32i = -im;
  B_9_32r = re;
  B_9_32i = im;
  y_root(&re, &im, 11, 32, 0);
  F_11_32r = re;
  F_11_32i = -im;
  B_11_32r = re;
  B_11_32i = im;
  y_root(&re, &im, 13, 32, 0);
  F_13_32r = re;
  F_13_32i = -im;
  B_13_32r = re;
  B_13_32i = im;
  y_root(&re, &im, 15, 32, 0);
  F_15_32r = re;
  F_15_32i = -im;
  B_15_32r = re;
  B_15_32i = im;


  /* radix 5 */
  y_root(&re, &im, 1, 5, 0);
  F_1_5r = re;
  F_1_5i = -im;
  B_1_5r = re;
  B_1_5i = im;
  y_root(&re,&im,2,5,0);
  F_2_5r = re;
  F_2_5i = -im;
  B_2_5r = re;
  B_2_5i = im;

  /* for Nussbaumer fft-5 */
  FN1_5r = -1.25;
  FN1_5i = F_1_5i + F_2_5i;
  FN2_5r = 0.5 * (F_1_5r - F_2_5r);
  FN2_5i = F_1_5i - F_2_5i;
  BN1_5r = -1.25;
  BN1_5i = B_1_5i + B_2_5i;
  BN2_5r = FN2_5r;
  BN2_5i = B_1_5i - B_2_5i;

  /* radix 10 */
  y_root(&re, &im, 1, 10,0);
  F_1_10r = re;
  F_1_10i = -im;
  B_1_10r = re;
  B_1_10i = im;
  y_root(&re, &im, 3, 10, 0);
  F_3_10r = re;
  F_3_10i = -im;
  B_3_10r = re;
  B_3_10i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_5);
  SSE2_VERS (B_1_5);
  SSE2_VERS (F_2_5);
  SSE2_VERS (B_2_5);
  SSE2_VERS (FN1_5);
  SSE2_VERS (BN1_5);
  SSE2_VERS (FN2_5);
  SSE2_VERS (BN2_5);
  SSE2_VERS (F_1_10);
  SSE2_VERS (B_1_10);
  SSE2_VERS (F_3_10);
  SSE2_VERS (B_3_10);
#endif

  /* radix 7*/
  y_root(&re, &im, 1, 7, 0);
  F_1_7r = re;
  F_1_7i = -im;
  B_1_7r = re;
  B_1_7i = im;
  y_root(&re, &im, 2, 7, 0);
  F_2_7r = re;
  F_2_7i = -im;
  B_2_7r = re;
  B_2_7i = im;
  y_root(&re, &im, 3, 7, 0);
  F_3_7r = re;
  F_3_7i = -im;
  B_3_7r = re;
  B_3_7i = im;

  /* For Nussbaumer fft-7 */
  FN1_7r = (F_1_7r + F_2_7r + B_3_7r) / 3.0;
  FN1_7i = (F_1_7i + F_2_7i + B_3_7i) / 3.0;
  FN2_7r = (2.0 * F_1_7r - F_2_7r - B_3_7r) / 3.0;
  FN2_7i = (2.0 * F_1_7i - F_2_7i - B_3_7i) / 3.0;
  FN3_7r = (F_1_7r - 2.0 * F_2_7r + B_3_7r) / 3.0;
  FN3_7i = (F_1_7i - 2.0 * F_2_7i + B_3_7i) / 3.0;
  FN4_7r = (F_1_7r + F_2_7r - 2.0 * B_3_7r) / 3.0;
  FN4_7i = (F_1_7i + F_2_7i - 2.0 * B_3_7i) / 3.0;

  BN1_7r = (B_1_7r + B_2_7r + F_3_7r) / 3.0;
  BN1_7i = (B_1_7i + B_2_7i + F_3_7i) / 3.0;
  BN2_7r = (2.0 * B_1_7r - B_2_7r - F_3_7r) / 3.0;
  BN2_7i = (2.0 * B_1_7i - B_2_7i - F_3_7i) / 3.0;
  BN3_7r = (B_1_7r - 2.0 * B_2_7r + F_3_7r) / 3.0;
  BN3_7i = (B_1_7i - 2.0 * B_2_7i + F_3_7i) / 3.0;
  BN4_7r = (B_1_7r + B_2_7r - 2.0 * F_3_7r) / 3.0;
  BN4_7i = (B_1_7i + B_2_7i - 2.0 * F_3_7i) / 3.0;

  /* radix 14 */
  y_root(&re, &im, 1, 14, 0);
  F_1_14r = re;
  F_1_14i = -im;
  B_1_14r = re;
  B_1_14i = im;
  y_root(&re, &im, 3, 14, 0);
  F_3_14r = re;
  F_3_14i = -im;
  B_3_14r = re;
  B_3_14i = im;
  y_root(&re, &im, 5, 14, 0);
  F_5_14r = re;
  F_5_14i = -im;
  B_5_14r = re;
  B_5_14i = im;

#if defined(Y_USE_SSE2)

  SSE2_VERS (F_1_7);
  SSE2_VERS (B_1_7);
  SSE2_VERS (F_2_7);
  SSE2_VERS (B_2_7);
  SSE2_VERS (F_3_7);
  SSE2_VERS (B_3_7);
  SSE2_VERS (FN1_7);
  SSE2_VERS (BN1_7);
  SSE2_VERS (FN2_7);
  SSE2_VERS (BN2_7);
  SSE2_VERS (FN3_7);
  SSE2_VERS (BN3_7);
  SSE2_VERS (FN4_7);
  SSE2_VERS (BN4_7);
  SSE2_VERS (F_1_14);
  SSE2_VERS (B_1_14);
  SSE2_VERS (F_3_14);
  SSE2_VERS (B_3_14);
  SSE2_VERS (F_5_14);
  SSE2_VERS (B_5_14);
#endif

}

/*
   This routine checks the step where pass 2 begins 
   It return the step 'ir' where Y_LRIGHT[i] <= Y_MEM_THRESHOLD
*/
int y_check_threshold_pass2(void)
{
  int ir = 1;
  if(Y_NRADICES <= 2)
    return ir;
  while ( Y_LRIGHT[ir] > Y_MEM_THRESHOLD)
    ir++;
  /*printf("pass 2 at ir=%d\n",ir);*/
  return ir;
}


/* Mr  Proper. Cleans the allocated memory for trig. factors. */
void Mr_Proper(void)
{
  y_size_t j;
  if(Y_PF != NULL)
    {
      for (j = 0; j < Y_OLDRADICES - 1; j++)
        free( (char *) Y_PF[j]);
      free( (char *) Y_PF);
      free( (char *) Y_TWDF);
    }
  if(Y_PB != NULL)
    {
      for (j = 0; j < Y_OLDRADICES - 1; j++)
        free( (char *) Y_PB[j]);
      free( (char *) Y_PB);
      free( (char *) Y_TWDB);
    }
}

/*
   This init the numerical factors. For every prime, numerical factors are 
   divided in two parts. First arrays, TP_NUMF[]  is for forward transform, 
   and has the bit reverse numerical factors. Second arrays TP_NUMB is for 
   backward transform, and has the numerical factors in normal place.
 
*/
void y_initnumfactg(void)
{
  y_limb_t tr = 0.0, ti = 0.0;
  y_ptr paux, ptar;
  y_size_t d, j, j1, np, nn, nrow, nr, i, inc;

  /*
     Allocates memory space for twidle forward and backward FFT
     trig factors 
  */

  if(Y_PF != NULL)
    {
      for (j = 0; j < Y_OLDRADICES - 1; j++)
        free( (char *) Y_PF[j]);
      free( (char *) Y_PF);
      free( (char *) Y_TWDF);
    }
  Y_PF = (y_ptr *) malloc((Y_NRADICES - 1) * sizeof(y_ptr));
  Y_TWDF = (y_ptr *) malloc((Y_NRADICES - 1) * sizeof(y_ptr));
  if ( Y_PF == NULL || Y_TWDF == NULL)
    {
      fprintf(stderr,"initnumfactg: Memory exhausted\n");
      exit(EXIT_FAILURE);
    }

  for (j = 0; j < Y_NRADICES - 1; j++)
    {
      d = Y_POWERS[Y_PLAN[j + 1] - 1];
      Y_PF[j] = ALLOC_DOUBLES((2 * Y_LLEFT[j] * d));
      if( Y_PF[j] == NULL)
        {
          fprintf(stderr,"initnumfactg: No memory to alloc trig. factors\n");
          exit(EXIT_FAILURE);
        }
      Y_TWDF[j] = ALIGN_DOUBLES(Y_PF[j]);
    }
#ifndef NDEBUG1
  fprintf(stderr,"Allocated space for forward trig factors\n");
#endif

  if(Y_PB != NULL)
    {
      for (j = 0; j < Y_OLDRADICES - 1; j++)
        free( (char *) Y_PB[j]);
      free( (char *) Y_PB);
      free( (char *) Y_TWDB);
    }
  Y_PB = (y_ptr *) malloc((Y_NRADICES - 1) * sizeof(y_ptr));
  Y_TWDB= (y_ptr *) malloc((Y_NRADICES - 1) * sizeof(y_ptr));
  if ( Y_PB == NULL || Y_TWDB == NULL)
    {
      fprintf(stderr,"initnumfactg: Memory exhausted\n");
      exit(EXIT_FAILURE);
    }

  for (j = 0; j < Y_NRADICES - 1; j++)
    {
      d = Y_POWERS[Y_PLAN[Y_NRADICES - 2 - j] - 1];
      Y_PB[j] = ALLOC_DOUBLES((2 * Y_LRIGHT[Y_NRADICES - 1 - j] * d));
      if( Y_PB[j] == NULL)
        {
          fprintf(stderr, "initnumfactg: No memory to alloc trig. factors\n");
          exit(EXIT_FAILURE);
        }
      Y_TWDB[j] = ALIGN_DOUBLES(Y_PB[j]);
    }
#ifndef NDEBUG1
  fprintf(stderr,"Allocated space for backward trig factors\n");
#endif

  /******************* TWIDDLE FACTORS FOR BACKWARD DIT *************/
  for(j=1;j<=Y_NRADICES-1;j++)
    {
      d=Y_LRIGHT[Y_NRADICES-j-1];
      i=Y_POWERS[Y_PLAN[Y_NRADICES-1-j]-1];
      for (j1=1;j1<Y_LRIGHT[Y_NRADICES-j];j1++)
        {
#if defined(Y_MAXIMUM)
          y_size_t ik,ij1;
          for (ik=1,ij1=0;ik<Y_PLAN[Y_NRADICES-j-1];ik+=1,ij1+=2)
            {
              y_root(&tr,&ti,ik*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+ij1]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+ij1+1]=ti;
            }

#else
          switch (Y_PLAN[Y_NRADICES-j-1])
            {
# ifndef Y_MINIMUM
            case (8):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,5*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              break;
            case (16):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,4*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              y_root(&tr,&ti,8*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+6]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+7]=ti;
              y_root(&tr,&ti,13*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+8]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+9]=ti;
              break;
            case (5):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,3*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              break;
            case (6):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,4*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              break;
            case (7):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,4*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,6*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              break;
            case (9):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,6*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              break;
            case (10):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,7*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              break;
            case (12):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,3*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              y_root(&tr,&ti,8*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+6]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+7]=ti;
              break;
            case (14):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,3*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              y_root(&tr,&ti,4*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+6]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+7]=ti;
              y_root(&tr,&ti,9*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+8]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+9]=ti;
              break;
#  if Y_AVAL > 4

            case (32):
                    y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
              y_root(&tr,&ti,2*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+2]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+3]=ti;
              y_root(&tr,&ti,3*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+4]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+5]=ti;
              y_root(&tr,&ti,4*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+6]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+7]=ti;
              y_root(&tr,&ti,5*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+8]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+9]=ti;
              y_root(&tr,&ti,6*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+10]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+11]=ti;
              y_root(&tr,&ti,7*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+12]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+13]=ti;
              y_root(&tr,&ti,15*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+14]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+15]=ti;
              y_root(&tr,&ti,27*j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)+16]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+17]=ti;
              break;
#  endif
# endif

            default:
              y_root(&tr,&ti,j1,d,0);
              Y_TWDB[j-1][((j1-1)*i<<1)]=tr;
              Y_TWDB[j-1][((j1-1)*i<<1)+1]=ti;
            }
#endif
#ifndef NDEBUG1
          fprintf(stderr,"%i %i  %10.5f %10.5f \n",j,j1,tr,ti);
#endif

        }
    }

#ifndef NDEBUG1
  printf(" Backward twidle factors computed\n");
#endif
  /************* TWIDLE FACTORS FOR FORWARD DIF **********************/
  /* It has been hard to me to find an efficient algorithm to compute
     the scrambled twiddle factors needed in this in-place-decimation- 
     in-frecuency (and in dyadic mul phase).

     Suposse the radix reduction plan is  (R0,R1,R2...R(N-1),RN). 
     In Forward transform we first make an R0 reduction, then R1... and  
     the last RN. On backward transform we begin with RN and finish with 
     R0. Only R0 uses to be a non-power-of-two radix. 
     
     The twidle factors for pass 0, reduction R0 are trivial one
     The twidle factor for pass 1, reduction R1 are (we only write
     the first of R1-1 powers needed for unit W)
  ROW 1) (0,1, ....(R0-1))/R0*R1
     The twidle factors for pass 2, reduction R2 are computed taken this
     basic  rows:
  ROW 1) (0,1,....(R1-1))/(R1*R2)
  ROW 2) (0,1,....(R0-1))/(R0*R1*R2)

     The twidle factors for pass 3, reduction R3 are computed using   
  ROW 1) (0,1,....(R2-1))/(R2*R3)
  ROW 2) (0,1,....(R1-1))/(R1*R2*R3)
  ROW 3) (0,1,....(R0-1))/(R0*R1*R2*R3)
     
     and so on, the last  
  ROW 1) (0,1,....(R(N-1) - 1))/(R(N-1)*RN)
  ROW 2) (0,1,....(R(N-2) - 1))/(R(N-2)*R(N-1)*RN)
  ROW 3) (0,1,....(R(N-3) - 1))/(R(N-3)*R(N-2)*R(N-1)*RN)
  .....
  ROW N-1) (0,1,....(R0 - 1))/(R0*R1*R2*....*RN)
     With this  basic rows, we construct a big row (it can be very long)
     1) Start moving ROW 1 to big row. (lenth R(N-1))
     2) Multiply all element in big row sequentaly by the first element of
  ROW 2.Then by the second ....Finally we will get a R(N-1)*R(N-2) big
  row. I'm talking about complex multiplication, it is the same than 
  addition of exponents written above.
     3) Repeat step 2 but now multiplying elements of big row by the elements  
  on ROW 3,  then ROW 4 up to ROW N-1. 
     At the end, big rows will became the rows of basic powers for twidle  
     factors in the passes. We will obtain the other powers trivialy if 
     needed.
  */
  /*
  np is the radix pass, from 1 to Y_NRADICES-1
  pass 0 has no  twiddle factors
  */
  for(np=1;np<Y_NRADICES;np++)
    {
      /* first twiddle factors are ones */
      inc=Y_POWERS[Y_PLAN[np]-1]*2;

      for(nn=0;nn<inc;nn+=2)
        {
          Y_TWDF[np-1][nn]=0.0;
          Y_TWDF[np-1][nn+1]=(y_limb_t)Y_LENGTH;
        }
      nn=1;
      for(nrow=1;nrow<=np;nrow++)
        {
          /* first we need to find the denominator on row  */
#ifndef NDEBUG1
          printf("Computing row %i of pass %i...\n",nrow,np);
#endif

          d=1;
          for(j=np-nrow;j<=np;j++)
            d*=Y_PLAN[j];

          paux=Y_TWDF[np-1];
          ptar= paux+nn*inc;
          for (nr=1;nr<Y_PLAN[np-nrow];nr++)
            {
              /*y_root(&tr,&ti,nr,d,1);*/
              for(i=0;i<nn;i++,ptar+=inc)
                {
#if defined(Y_MAXIMUM)
                  y_size_t ik,ikinc;
                  for(ik=1,ikinc=0;ik<Y_PLAN[np];ik++,ikinc+=2)
                    {
                      ptar[ikinc]=paux[i*inc + ikinc]+(y_limb_t)(ik*(nr*Y_LENGTH/d));
                      ptar[ikinc+1]=(y_limb_t)(Y_LENGTH);
                    }
#else
                  switch(Y_PLAN[np])
                    {
# ifndef Y_MINIMUM
                    case(8):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(5*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(5):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(3*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(6):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(4*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(7):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(4*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(6*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(9):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(6*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(10):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(7*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(12):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(3*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      ptar[6]=paux[i*inc+6]+(y_limb_t)(8*(nr*Y_LENGTH/d));
                      ptar[7]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(14):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(3*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      ptar[6]=paux[i*inc+6]+(y_limb_t)(4*(nr*Y_LENGTH/d));
                      ptar[7]=(y_limb_t)(Y_LENGTH);
                      ptar[8]=paux[i*inc+8]+(y_limb_t)(9*(nr*Y_LENGTH/d));
                      ptar[9]=(y_limb_t)(Y_LENGTH);
                      break;
                    case(16):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(4*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      ptar[6]=paux[i*inc+6]+(y_limb_t)(8*(nr*Y_LENGTH/d));
                      ptar[7]=(y_limb_t)(Y_LENGTH);
                      ptar[8]=paux[i*inc+8]+(y_limb_t)(13*(nr*Y_LENGTH/d));
                      ptar[9]=(y_limb_t)(Y_LENGTH);
                      break;
#  if Y_AVAL > 4

                    case(32):
                            ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      ptar[2]=paux[i*inc+2]+(y_limb_t)(2*(nr*Y_LENGTH/d));
                      ptar[3]=(y_limb_t)(Y_LENGTH);
                      ptar[4]=paux[i*inc+4]+(y_limb_t)(3*(nr*Y_LENGTH/d));
                      ptar[5]=(y_limb_t)(Y_LENGTH);
                      ptar[6]=paux[i*inc+6]+(y_limb_t)(4*(nr*Y_LENGTH/d));
                      ptar[7]=(y_limb_t)(Y_LENGTH);
                      ptar[8]=paux[i*inc+8]+(y_limb_t)(5*(nr*Y_LENGTH/d));
                      ptar[9]=(y_limb_t)(Y_LENGTH);
                      ptar[10]=paux[i*inc+10]+(y_limb_t)(6*(nr*Y_LENGTH/d));
                      ptar[11]=(y_limb_t)(Y_LENGTH);
                      ptar[12]=paux[i*inc+12]+(y_limb_t)(7*(nr*Y_LENGTH/d));
                      ptar[13]=(y_limb_t)(Y_LENGTH);
                      ptar[14]=paux[i*inc+14]+(y_limb_t)(15*(nr*Y_LENGTH/d));
                      ptar[15]=(y_limb_t)(Y_LENGTH);
                      ptar[16]=paux[i*inc+16]+(y_limb_t)(27*(nr*Y_LENGTH/d));
                      ptar[17]=(y_limb_t)(Y_LENGTH);
                      break;
#  endif
# endif

                    default:
                      ptar[0]=paux[i*inc]+(y_limb_t)(nr*Y_LENGTH/d);
                      ptar[1]=(y_limb_t)(Y_LENGTH);
                      break;
                    }
                  /*printf("%i -->  %10.6f %10.6f \n",np,ptar[0],ptar[1]);*/
#endif

                }
              /*paux=Y_TWDF[np-1];*/
            }
          nn*=Y_PLAN[np-nrow];
        }
      /*now, compute the trig factors */
      paux=Y_TWDF[np-1];
      for(i=0;i<inc*Y_LLEFT[np-1];i+=2)
        {
          j=floor(paux[i]+0.5);
          /*printf(" %i \n",j);*/
          y_root(&tr,&ti,j,Y_LENGTH,1);
          paux[i]=tr;
          paux[i+1]=ti;
        }
    }
#ifndef NDEBUG1
  printf(" Forward twidle factors computed\n");
  for(i=0;i<Y_NRADICES-1;i++)
    {
      inc=Y_POWERS[Y_PLAN[i+1]-1]*2;
      for(j=0;j<Y_LLEFT[i]*inc;j+=2)
        {
          printf(" %i %i %10.5f %10.5f \n",i,j,Y_TWDF[i][j],Y_TWDF[i][j+1]);
        }
    }
#endif

  /*************** TRIG FACTORS FOR DYADIC NESTED MUL ***************/
  if( Y_D != NULL)
    free( (char *) Y_D);
  Y_D = ALLOC_DOUBLES(Y_PLAN[Y_NRADICES-1]*2);
  Y_DYADIC = ALIGN_DOUBLES(Y_D);
  if (Y_D == NULL || Y_DYADIC == NULL)
    {
      fprintf(stderr,"initnumfactg: Memory exhausted\n");
      exit(EXIT_FAILURE);
    }
  for (j = 0, i = 0; j < Y_LENGTH;j += Y_LLEFT[Y_NRADICES - 2],i++)
    {
      y_root(&tr,&ti,j,Y_LENGTH,1);
      Y_DYADIC[i<<1]=tr;
      Y_DYADIC[(i<<1)+1]=ti;
      /*printf(" %10.5f %10.5f \n",tr,ti);*/
    }
#ifndef NDEBUG1
  printf(" Dyadic mul factors computed\n");
#endif

}


/********************************************************************
This routine:
    1) Returns the length for a needed n
    2) Set the pointers to the correct data acording to n
    3) Make the sequence of radix reduction plan, and... 
    4) Call to init numeric twidle factors.
 
INPUT:
	n = is the minimum FFT needed length (in terms of BIG_DOUBLE).
            In the caller program it has been derived from Y_BITS (the bits 
	    of a limb) and the overall bits.
 
OUTPUT:
	returns the optimal length=q*2^k. (in terms of COMPLEX size). 
	Also, those q and k has been stored as globals Y_Q and Y_K
	
GLOBALS: 
	Y_AVAL= is the maximum power_of_two-radix avalaible. I .e. 
	      if we set TP_AVAL=4, we have radix 2,4,8, and 16 
 
NOTES: we suposse n<2^31.  
*****************************************************************/

y_size_t y_init(y_size_t n)
{
  y_size_t i = (1L << (sizeof(y_size_t) * 8 - 2)), j,
               k = sizeof(y_size_t) * 8 - 1, old = Y_LENGTH;

  /*
     If we've called this routine before to work with the same length
     then all work is done so we return now 
  */
  if((n >> 1) == old)
    return old;

  while (i >= n )
    {
      k--;
      i >>= 1;
    } ;
  /* now k is the smallest power of two great or equal than n */
  j = 1 << (k - 3);
  for(i = 8; i > 4; i--)
    {
      if(j * i < n)
        break;
    }
  i++;
  if((i == 5) && (k > 3) && ((j >> 1) * 9 >= n))
    {
      i = 9;
      j >>= 1;
      k--;
    }


  /* Second chance to skip */
  if((( i * j) >> 1) == old )
    return( old );

  Y_LENGTH = (i * j) >> 1;
#ifndef NDEBUG1

  printf (" %i \n", Y_LENGTH);
#endif
  /* select the right pointers to data for selected length=q*2^k */
  Y_Q = i;
  Y_K = k - 4;
  if(Y_K > k)
    {
      fprintf (stderr, "y_init: Error assigning fft plan.\n");
      exit (EXIT_FAILURE);
    }

  Y_NRADICES = Y_K / Y_AVAL + 2;
  j= Y_K % Y_AVAL;
#if Y_AVAL == 3
  /* This cond #if is a patch to avoid the use of difdit_16 for x86 */
  if(j == 0)
    Y_NRADICES--;
#endif

#if defined(Y_MANY_REGISTERS) && (Y_AVAL > 3)

  if((j == 1) && (Y_Q != 9))
    {
      Y_Q <<= 1;
      Y_K--;
      Y_NRADICES--;
    }
  else if (Y_Q != 9)
    {
      Y_Q <<= 1;
      Y_K--;
    }
#endif

  if(Y_NRADICES < 2)
    {
      fprintf (stderr, "\ny_init: FFT length too small for yeafft.\n");
      exit (EXIT_FAILURE);
    }

  if (Y_PLAN != NULL)
    {
      free ((char *) Y_PLAN);
      free ((char *) Y_LLEFT);
      free((char *)Y_LRIGHT);
      /*
      #if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)
      if(Y_2N0 !=NULL) free((char *)Y_2N0);
      if(Y_2NN !=NULL) free((char *)Y_2NN);
      if(Y_DI != NULL) free((char *)Y_DI);
      if(Y_DJ != NULL) free((char *)Y_DJ);
      Y_2N0 = NULL; Y_2NN=NULL; Y_DI=NULL; Y_DJ=NULL;
      # if defined(_PTHREADS)
      y_threads_destroy();
      if(y_threads != NULL) free((char *)y_threads);
      # endif
      #endif
      */
    }
  Y_PLAN = (y_size_t *)(malloc(20 * sizeof(y_size_t)));
  Y_LLEFT = (y_size_t *)(malloc(20 * sizeof(y_size_t)));
  Y_LRIGHT = (y_size_t *)(malloc(20 * sizeof(y_size_t)));
  if( Y_PLAN == NULL || Y_LLEFT == NULL || Y_LRIGHT == NULL)
    {
      fprintf(stderr,"y_init: Memory exhausted\n");
      exit(EXIT_FAILURE);
    }
  /*
  #if defined(_PTHREADS)
  y_threads_init();
  #endif
  */
  Y_PLAN[0] = Y_Q;

  for(i = 1; i < Y_NRADICES; i++)
    Y_PLAN[i] = 1;
  j = Y_NRADICES - 1;
  for(i = 0; i < Y_K; i++)
    Y_PLAN[i % j + 1] <<= 1;

  /* Compute Y_LLEFT and Y_LRIGHT */
  Y_LLEFT[0] = Y_PLAN[0];
  Y_LRIGHT[Y_NRADICES - 1] = Y_PLAN[Y_NRADICES-1];
  for(i = 1; i < Y_NRADICES; i++)
    {
      Y_LLEFT[i] = Y_LLEFT[i - 1] * Y_PLAN[i];
      Y_LRIGHT[Y_NRADICES - 1 - i] = Y_LRIGHT[Y_NRADICES - i] *
                                     Y_PLAN[Y_NRADICES - 1 - i];
    }

#ifndef NDEBUG1
  /* To be sure the radix-plan is correct */
  i = 1;
  for (j = 0;j < Y_NRADICES; j++)
    {
      i *= Y_PLAN[j];
      printf(" %i %i %i %i\n",Y_PLAN[j],Y_LRIGHT[j],Y_LLEFT[j],i);
    }
  assert (i == Y_LENGTH);
#endif

#if defined(Y_NEW_PADDING)

  i = 1;
  Y_SHIFT0 = 1;
  switch (Y_PLAN[Y_NRADICES - 1])
    {
    case 4:
      i += 2;
      break;
    case 8:
      i += 3;
      break;
    case 16:
      i += 4;
      break;
    }
  if( Y_NRADICES > 1)
    {
      switch (Y_PLAN[Y_NRADICES - 2])
        {
        case 4:
          i += 2;
          break;
        case 8:
          i += 3;
          break;
        case 16:
          i += 4;
          break;
        case 32:
          i += 5;
          break;
        }
    }
  Y_SHIFT1 = i;
  if( Y_NRADICES > 2)
    {
      switch (Y_PLAN[Y_NRADICES - 3])
        {
        case 4:
          i += 2;
          break;
        case 8:
          i += 3;
          break;
        case 16:
          i += 4;
          break;
        case 32:
          i += 5;
          break;
        }
    }
  Y_SHIFT2 = i;
  if( Y_NRADICES > 3)
    {
      switch (Y_PLAN[Y_NRADICES - 4])
        {
        case 4:
          i += 2;
          break;
        case 8:
          i += 3;
          break;
        case 16:
          i += 4;
          break;
        case 32:
          i += 5;
          break;
        }
    }
  Y_SHIFT3 = i;
  if( Y_NRADICES > 4)
    {
      switch (Y_PLAN[Y_NRADICES - 5])
        {
        case 4:
          i += 2;
          break;
        case 8:
          i += 3;
          break;
        case 16:
          i += 4;
          break;
        case 32:
          i += 5;
          break;
        }
    }
  Y_SHIFT4 = i;
#endif

#ifdef SEE_PLAN

  for (j = 0; j < Y_NRADICES; j++)
    {
      printf (" %i,", Y_PLAN[j]);
    }
  printf ("\n");
  /* this is for a pentium 3, it is referred as power of two */
# define L_LINES 9
# define L_MASK ((1L << L_LINES) - 1)
  /* Number of doubles in a cache line as power of two */
# define L_CACHE_LINE 2

  for (j = 0; j < Y_NRADICES; j++)
    {
      printf ("Cache lines pad, try to avoid 0 values!\n");
      for(k = 1; k <= Y_PLAN[j]; k++)
        {
          i = addr((k * Y_LRIGHT[j] / Y_PLAN[j]) << 1);
          i >>= L_CACHE_LINE;
          i &= L_MASK;
          printf (" %0X,",i);
        }
      printf ("\n");
    }
#endif

  /* New oportunity to skip the heavy routines */
  if(old == Y_LENGTH)
    return(Y_LENGTH);

  /* Set Y_PASS2 */
  Y_PASS2 = y_check_threshold_pass2();

  /* init the numeric factors */
  y_initnumfactg ();

  Y_OLDRADICES = Y_NRADICES;

  /* init the radix routines */
  y_initradix ();
  return Y_LENGTH;
}
/*$Id$*/
