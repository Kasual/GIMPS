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
#define Y_SAVE_TRIGS
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
#elif defined(Y_VECTORIZE2)
# include "ygenvect.h"
# ifdef Y_LONG_MACROS
#  undef Y_LONG_MACROS
# endif
#else
# ifdef Y_LONG_MACROS
#  include "ygeneric4.h"
# endif
#endif



/* common blocks  code */

#define get_pads_4         \
  bigpad = (pad<<2);       \
  jj = pad;                \
  pd0 = d;                 \
  prefetch_data( pd0, 0);  \
  pd1 = d + addr((jj<<1)); \
  jj += pad;               \
  prefetch_data( pd1, 0);  \
  pd2 = d + addr((jj<<1)); \
  jj += pad;               \
  prefetch_data( pd2, 0);  \
  pd3 = d + addr((jj<<1)); \
  prefetch_data( pd3, 0);


#define get_pads_4_preload(_n0) \
  bigpad=(pad<<2);              \
  jj= pad + _n0;                \
  pd0= d + addr(_n0<<1);        \
  pd1= d + addr((jj<<1));       \
  jj+=pad;                      \
  pd2= d + addr((jj<<1));       \
  jj+=pad;                      \
  pd3= d + addr((jj<<1));


#if defined(Y_MAXIMUM)

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

# define get_twiddle_factors_4             \
	  tw1r=px[0];                      \
	  tw1i=px[1];                      \
	  tw2r=px[2];                      \
	  tw2i=px[3];                      \
	  tw3r=px[4];                      \
	  tw3i=px[5];                      \
          prefetch_p_trig(px);

# ifdef Y_PRELOAD_TRIGS
#  define trig_4_load_init(_px)              \
          tpw1r = 1.0;                       \
          tpw1i = 0.0;                       \
          tpw2r = 1.0;                       \
          tpw2i = 0.0;                       \
          tpw3r = 1.0;                       \
          tpw3i = 0.0;                       \
          upw1r = *(_px);                    \
          upw1i = *(_px + 1);                \
          upw2r = *(_px + 2);                \
          upw2i = *(_px + 3);                \
          upw3r = *(_px + 4);                \
          upw3i = *(_px + 5);

#  define trig_4_preload(_px)                \
          tpw1r = *(_px);                    \
          tpw1i = *(_px + 1);                \
          tpw2r = *(_px + 2);                \
          tpw2i = *(_px + 3);                \
          tpw3r = *(_px + 4);                \
          tpw3i = *(_px + 5);                \
          upw1r = *(_px + 6);                \
          upw1i = *(_px + 7);                \
          upw2r = *(_px + 8);                \
          upw2i = *(_px + 9);                \
          upw3r = *(_px + 10);               \
          upw3i = *(_px + 11);

#  define get_twiddle_factors_4_preload      \
          tw1r = tpw1r;  tw1i = tpw1i;       \
          tw2r = tpw2r;  tw2i = tpw2i;       \
          tw3r = tpw3r;  tw3i = tpw3i;       \
          uw1r = upw1r;  uw1i = upw1i;       \
          uw2r = upw2r;  uw2i = upw2i;       \
          uw3r = upw3r;  uw3i = upw3i;       \
          trig_4_preload(px);

#  define trig_4_load_dif_init(_px)          \
          tpw1r = *(_px );                   \
          tpw1i = *(_px + 1);                \
          tpw2r = *(_px + 2);                \
          tpw2i = *(_px + 3);                \
          tpw3r = *(_px + 4);                \
          tpw3i = *(_px + 5);                \

#  define trig_4_dif_preload(_px)            \
          tpw1r = *(_px + 6);                \
          tpw1i = *(_px + 7);                \
          tpw2r = *(_px + 8);                \
          tpw2i = *(_px + 9);                \
          tpw3r = *(_px + 10);               \
          tpw3i = *(_px + 11);

#  define get_twiddle_factors_4_dif_preload  \
          tw1r = tpw1r; tw1i = tpw1i;        \
          tw2r = tpw2r; tw2i = tpw2i;        \
          tw3r = tpw3r; tw3i = tpw3i;        \
          trig_4_dif_preload(px);

# endif

#else

# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# define get_twiddle_factors_4             \
	  tw2r=(px[0]+px[1])*(px[0]-px[1]);\
	  tw2i=2.0*px[0]*px[1];            \
	  tw1r=px[0];                      \
	  tw1i=px[1];                      \
	  tw3r=px[0]*tw2r - px[1]*tw2i;    \
	  tw3i=px[0]*tw2i + px[1]*tw2r;    \
          prefetch_p_trig(px);


/* Macros used frecuently */

# ifdef Y_PRELOAD_TRIGS
#  define trig_4_load_init(_px)              \
          tpw1r = 1.0;                       \
          tpw1i = 0.0;                       \
          upw1r = *(_px);                    \
          upw1i = *(_px + 1);

#  define trig_4_preload(_px)                \
          tpw1r = *(_px);                    \
          tpw1i = *(_px + 1);                \
          upw1r = *(_px + 2);                \
          upw1i = *(_px + 3);

#  define get_twiddle_factors_4_preload      \
{                                            \
          y_limb_t a0,b0,a1,b1;              \
          tw2i = tpw1r * tpw1i;              \
          uw2i = upw1r * upw1i;              \
          a0 = tpw1i * tpw1i;                \
          a1 = upw1i * upw1i;                \
          tw1r = tpw1r; tw1i = tpw1i;        \
          uw1r = upw1r; uw1i = upw1i;        \
          tw2i += tw2i;                      \
          uw2i += uw2i;                      \
          tw2r = tw1r * tw1r - a0;           \
          uw2r = uw1r * uw1r - a1;           \
          trig_4_preload(px);                \
          a0 = tw1i * tw2i;                  \
          b0 = tw1r * tw2i;                  \
          a1 = uw1i * uw2i;                  \
          b1 = uw1r * uw2i;                  \
          tw3r = tw1r * tw2r - a0;           \
          tw3i = tw1i * tw2r + b0;           \
          uw3r = uw1r * uw2r - a1;           \
          uw3i = uw1i * uw2r + b1;           \
}

#  define trig_4_load_dif_init(_px)          \
          tpw1r = *(_px );                   \
          tpw1i = *(_px + 1);

#  define trig_4_dif_preload(_px)            \
          tpw1r = *(_px + 2);                \
          tpw1i = *(_px + 3);

#  define get_twiddle_factors_4_dif_preload  \
{                                            \
          y_limb_t a0,b0;                    \
          a0 = tpw1i * tpw1i;                \
          tw2i = tpw1r * tpw1i;              \
          tw1r = tpw1r; tw1i = tpw1i;        \
          tw2r = tpw1r * tpw1r - a0;         \
          tw2i += tw2i;                      \
          trig_4_dif_preload(px);            \
          a0 = tw1i * tw2i;                  \
          b0 = tw1r * tw2i;                  \
          tw3r = tw1r * tw2r - a0;           \
          tw3i = tw1i * tw2r + b0;           \
}

# endif
#endif

#ifdef Y_PRELOAD_DATA

# define radix_4_twd_preload_B(_k)                                        \
   t0r = s0r; t0i = s0i; ppd0r = ppdl0r; ppd0i = ppdl0i;                  \
   ppdl0r +=_k; ppdl0i += _k; s0r= *(ppdl0r); s0i=*(ppdl0i);              \
   u0r = v0r; u0i = v0i; qqd0r = qqdl0r; qqd0i = qqdl0i;                  \
   qqdl0r +=_k; qqdl0i += _k; v0r= *(qqdl0r); v0i=*(qqdl0i);              \
   cplx_three_muls_fftn_preload(t1,t2,t3,2,1,3,s,tw,ppd,ppdl,_k);         \
   cplx_three_muls_fftn_preload(u1,u2,u3,2,1,3,v,uw,qqd,qqdl,_k);         \
   cplx_two_fft4_preload_store_B(ppd0,ppd1,ppd2,ppd3,qqd0,qqd1,qqd2,qqd3,t,u,0,1,2,3);

# define radix_4_notwd_preload_B(_k)                                      \
  cplx_two_fft4_notwd_preload_store_B(ppd0,ppd1,ppd2,ppd3,qqd0,qqd1,qqd2,qqd3,t,s,u,v,0,2,1,3,ppd,ppdl,qqd,qqdl,_k);

# define radix_4_twd_preload_F(_k)                                      \
   t0r = s0r; t0i = s0i; ppd0r = ppdl0r; ppd0i = ppdl0i;                  \
   ppdl0r +=_k; ppdl0i += _k; s0r= *(ppdl0r); s0i=*(ppdl0i);              \
   u0r = v0r; u0i = v0i; qqd0r = qqdl0r; qqd0i = qqdl0i;                  \
   qqdl0r +=_k; qqdl0i += _k; v0r= *(qqdl0r); v0i=*(qqdl0i);              \
   cplx_three_muls_fftn_preload(t1,t2,t3,2,1,3,s,tw,ppd,ppdl,_k);         \
   cplx_three_muls_fftn_preload(u1,u2,u3,2,1,3,v,tw,qqd,qqdl,_k);         \
   cplx_two_fft4_preload_store_F(ppd0,ppd1,ppd2,ppd3,qqd0,qqd1,qqd2,qqd3,t,u,0,1,2,3);

# define radix_4_notwd_preload_F(_k)                                      \
  cplx_two_fft4_notwd_preload_store_F(ppd0,ppd1,ppd2,ppd3,qqd0,qqd1,qqd2,qqd3,t,s,u,v,0,2,1,3,ppd,ppdl,qqd,qqdl,_k);

#endif

#ifdef Y_VECTORIZE2

#if defined(Y_MAXIMUM)

# define get_twiddle_factors_4_vector      \
          tw1r = *(px);                    \
          tw1i = *(px + 1);                \
          tw2r = *(px + 2);                \
          tw2i = *(px + 3);                \
          tw3r = *(px + 4);                \
          tw3i = *(px + 5);                \
          uw1r = *(px + 6);                \
          uw1i = *(px + 7);                \
          uw2r = *(px + 8);                \
          uw2i = *(px + 9);                \
          uw3r = *(px + 10);               \
          uw3i = *(px + 11);               \
	  prefetch_p_trig(px);             \
	  prefetch_p_trig(px+Y_CACHE_LINE);

#else
# define get_twiddle_factors_4_vector      \
	  tw2r=(px[0]+px[1])*(px[0]-px[1]);\
          uw2r=(px[2]+px[3])*(px[2]-px[3]);\
	  tw2i=2.0*px[0]*px[1];            \
	  uw2i=2.0*px[2]*px[3];            \
	  tw1r=px[0];                      \
	  tw1i=px[1];                      \
	  uw1r=px[2];                      \
	  uw1i=px[3];                      \
	  tw3r=px[0]*tw2r - px[1]*tw2i;    \
	  uw3r=px[2]*uw2r - px[3]*uw2i;    \
	  tw3i=px[0]*tw2i + px[1]*tw2r;    \
	  uw3i=px[2]*uw2i + px[3]*uw2r;    \
          prefetch_p_trig(px);
#endif

# define radix_4_twd_pp_vector_F(_j)                                         \
	  cplx_load_muladdsub_pp_vector_F(t0,t1,u0,u1,tw2,_j,pd0,pd2);       \
	  cplx_load_mulmuladdsub_pp_vector_F(t2,t3,u2,u3,tw1,tw3,_j,pd1,pd3);\
	                                                                     \
	  cplx_addsub_store_p_vector(_j,pd0,pd2,t0,t2,u0,u2);                \
	  cplx_mul_1_4_F_addsub_store_p_vector(_j,pd1,pd3,t1,t3,u1,u3);

# define radix_4_twd_pp_vector_B(_j)                                         \
	  cplx_load_muladdsub_pp_vector_B(t0,t1,u0,u1,tw2,uw2,_j,pd0,pd2);   \
	  cplx_load_mulmuladdsub_pp_vector_B(t2,t3,u2,u3,tw1,tw3,uw1,uw3,_j,pd1,pd3);\
	                                                                     \
	  cplx_addsub_store_p_vector(_j,pd0,pd2,t0,t2,u0,u2);                \
	  cplx_mul_1_4_B_addsub_store_p_vector(_j,pd1,pd3,t1,t3,u1,u3);


# define radix_4_notwd_pp_vector(_S,_j)                                      \
	  cplx_load_addsub_pp_vector(t0,t1,u0,u1,_j,pd0,pd2);                \
	  cplx_load_addsub_pp_vector(t2,t3,u2,u3,_j,pd1,pd3);                \
	  \
	  cplx_addsub_store_p_vector(_j,pd0,pd2,t0,t2,u0,u2);                \
	  cplx_mul_1_4_##_S##_addsub_store_p_vector(_j,pd1,pd3,t1,t3,u1,u3);

#endif


#ifdef Y_LONG_MACROS

# define radix_4_twd_pp(_S,_j)                        \
          cplx_radix_4_twd_pp_passes_##_S(t,tw,pd,_j);

# define radix_4_twd_pp_no_fetch(_S,_j)                       \
          cplx_radix_4_twd_pp_passes_no_fetch_##_S(t,tw,pd,_j);

# define radix_4_notwd_pp(_S,_j)                    \
          cplx_radix_4_notwd_pp_passes_##_S(t,pd,_j);

# define radix_4_notwd_pp_no_fetch(_S,_j)                    \
          cplx_radix_4_notwd_pp_passes_no_fetch_##_S(t,pd,_j);

#else

# define radix_4_twd_pp(_S,_j)                                \
	  cplx_load_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);       \
	  cplx_load_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3);\
	                                                      \
	  cplx_addsub_store_p(_j,pd0,pd2,t0,t2);              \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);

# define radix_4_twd_pp_no_fetch(_S,_j)                                \
	  cplx_load_muladdsub_pp_no_fetch(t0,t1,tw2,_j,pd0,pd2);       \
	  cplx_load_mulmuladdsub_pp_no_fetch(t2,t3,tw1,tw3,_j,pd1,pd3);\
	  \
	  cplx_addsub_store_p(_j,pd0,pd2,t0,t2);               \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);

# define radix_4_notwd_pp(_S,_j)                               \
	  cplx_load_addsub_pp(t0,t1,_j,pd0,pd2);               \
	  cplx_load_addsub_pp(t2,t3,_j,pd1,pd3);               \
	                                                       \
	  cplx_addsub_store_p(_j,pd0,pd2,t0,t2);               \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);

# define radix_4_notwd_pp_no_fetch(_S,_j)                      \
	  cplx_load_addsub_pp_no_fetch(t0,t1,_j,pd0,pd2);      \
	  cplx_load_addsub_pp_no_fetch(t2,t3,_j,pd1,pd3);      \
	                                                       \
	  cplx_addsub_store_p(_j,pd0,pd2,t0,t2);               \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd1,pd3,t1,t3);


#endif

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
*************************************************************************/

void radix_4_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
#if defined Y_PRELOAD_TRIGS

  y_limb_t tpw1r,tpw1i;
# if defined(Y_MAXIMUM)

  y_limb_t tpw2r,tpw2i,tpw3r,tpw3i;
# endif
#endif
#if defined(Y_VECTORIZE2) || defined(Y_PRELOAD_DATA)

  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
#endif
#ifdef Y_PRELOAD_DATA

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i;
  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i;
  y_ptr qqd0r,qqd0i,qqd1r,qqd1i,qqd2r,qqd2i,qqd3r,qqd3i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i;
  y_ptr qqdl0r,qqdl0i,qqdl1r,qqdl1i,qqdl2r,qqdl2i,qqdl3r,qqdl3i;
#endif

  y_ptr px;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t bigpad;
  y_size_t i,j,jj,nc;
#ifdef Y_PRELOAD_DATA

  y_size_t next,last,nexti;
#endif

  ASSERT_ALIGNED_DOUBLE();


#ifndef NDEBUG1

  printf(" Radix_4_dif n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

#ifdef Y_PRELOAD_DATA

  get_pads_4_preload(n0);
#else

  get_pads_4;
#endif

#ifdef Y_PRELOAD_DATA

  cplx_two_radix_4_dif_dit_init(s,ppdl,v,qqdl,pd);
  last = addr(n0<<1);
#endif

  /*
     If n0==0 first elements has twiddle factors equal to one.
     then we try to save this twiddle factors muls.
  */
  if(n0==0)
    {
#ifdef Y_PRELOAD_DATA
      if (bigpad<n)
        nexti = bigpad;
      else
        nexti = pad -2;
      for(j=0;j<pad;j+=2)
        {
          if((j+2)<pad)
            next = j+2;
          else
            next = nexti;
          next = addr(next<<1);
          jj= next - last;
          last = next;
          radix_4_notwd_preload_F(jj);
        }
#elif defined(Y_VECTORIZE2)
      for(j=0;j<pad;j+=2)
        {
          jj=addr(j<<1);
          radix_4_notwd_pp_vector(F,jj);
        }
#else
      for(j=0;j<pad;j++)
        {
          jj=addr(j<<1);
          radix_4_notwd_pp(F,jj);
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

          j++;
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
#  if Y_CACHE_LINE == 8

          j+=2;
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
#  endif
# endif

        }
#endif
      if(bigpad==n)
        return;
      px=tw + Y_STEP;
#ifdef Y_PRELOAD_TRIGS

      trig_4_load_dif_init(px);
#endif

      for(i=bigpad;i<n;i+=bigpad, px += Y_STEP)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_4_dif_preload;
#else

          get_twiddle_factors_4;
#endif
#ifdef Y_PRELOAD_DATA

          if((i+bigpad)<n)
            nexti = i + bigpad;
          else
            nexti= i + pad - 2 ;
          for(j=i;j<(pad+i);j+=2)
            {
              if((j+2)<(pad+i))
                next = j+2;
              else
                next=nexti;
              next = addr(next<<1);
              jj= next - last;
              last = next;
              radix_4_twd_preload_F(jj);
            }
#elif defined(Y_VECTORIZE2)
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radix_4_twd_pp_vector_F(jj);
            }
#else
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_4_twd_pp(F,jj);
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
#  if Y_CACHE_LINE == 8

              j+=2;
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
#  endif
# endif

            }
#endif

        }
    }
  else
    {
#ifdef Y_PRELOAD_TRIGS
      trig_4_load_dif_init(tw);
#endif

      nc=n+n0;
      for(i=n0,px=tw;i<nc;i+=bigpad, px += Y_STEP)
        {
#ifdef Y_PRELOAD_TRIGS
          get_twiddle_factors_4_dif_preload;
#else

          get_twiddle_factors_4;
#endif
#ifdef Y_PRELOAD_DATA

          if((i+bigpad)<nc)
            nexti =i + bigpad;
          else
            nexti= i + pad - 2;
          for(j=i;j<(pad+i);j+=2)
            {
              if((j+2)<(pad+i))
                next=j+2;
              else
                next=nexti;
              next = addr(next<<1);
              jj= next - last;
              last = next;
              radix_4_twd_preload_F(jj);
            }
#elif defined(Y_VECTORIZE2)
          for(j=i;j<(pad+i);j+=2)
            {
              jj=addr(j<<1);
              radix_4_twd_pp_vector_F(jj);
            }
#else
          for(j=i;j<(pad+i);j++)
            {
              jj=addr(j<<1);
              radix_4_twd_pp(F,jj);
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

              j++;
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
#  if Y_CACHE_LINE == 8

              j+=2;
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
              jj+=2;
              radix_4_twd_pp_no_fetch(F,jj);
#  endif
# endif

            }
#endif

        }
    }
}

/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r  WITHOUT TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw[] = NULL in this no twiddle version. 
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radix_4_dif_notw(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
#ifdef Y_VECTORIZE2

  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i;
#endif

  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;

  ASSERT_ALIGNED_DOUBLE();


#ifndef NDEBUG1

  printf(" Radix_4_dif_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  if(tw != NULL)
    return;

  get_pads_4;

  nc=n+n0;
  for(i=n0;i<nc;i+=bigpad)
    {
#ifdef Y_VECTORIZE2
      for(j=i;j<(pad+i);j+=2)
        {
          jj=addr(j<<1);
          radix_4_notwd_pp_vector(F,jj);
        }
#else
      for(j=i;j<(pad+i);j++)
        {
          jj=addr(j<<1);
          radix_4_notwd_pp(F,jj);
# if defined(Y_PREFETCH_EXPENSIVE) && (Y_CACHE_LINE >= 4)

          j++;
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
#  if Y_CACHE_LINE == 8

          j+=2;
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
          jj+=2;
          radix_4_notwd_pp_no_fetch(F,jj);
#  endif
# endif

        }
#endif

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

void radix_4_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0,y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
#if defined(Y_VECTORIZE2) || defined(Y_PRELOAD_TRIGS)

  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,uw1r,uw1i,uw2r,uw2i,uw3r,uw3i;
#endif
#ifdef Y_PRELOAD_TRIGS

  y_limb_t tpw1r,tpw1i;
  y_limb_t upw1r,upw1i;
# if defined(Y_MAXIMUM)

  y_limb_t tpw2r,tpw2i,tpw3r,tpw3i;
  y_limb_t upw2r,upw2i,upw3r,upw3i;
# endif
#endif
#ifdef Y_PRELOAD_DATA

  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i;
  y_ptr ppd0r,ppd0i,ppd1r,ppd1i,ppd2r,ppd2i,ppd3r,ppd3i;
  y_ptr qqd0r,qqd0i,qqd1r,qqd1i,qqd2r,qqd2i,qqd3r,qqd3i;
  y_ptr ppdl0r,ppdl0i,ppdl1r,ppdl1i,ppdl2r,ppdl2i,ppdl3r,ppdl3i;
  y_ptr qqdl0r,qqdl0i,qqdl1r,qqdl1i,qqdl2r,qqdl2i,qqdl3r,qqdl3i;
#endif

  y_ptr px;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t bigpad;
  y_size_t i,j,jj,nc;
#ifdef Y_SAVE_TRIGS

  y_size_t i1,j1,jj1;
#endif
#ifdef Y_PRELOAD_DATA

  y_size_t next,last,nexti;
#endif

  ASSERT_ALIGNED_DOUBLE();


#ifndef NDEBUG1

  printf(" Radix_4_dit n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  nc=n0+n;

#ifdef Y_PRELOAD_DATA

  get_pads_4_preload(n0);
#else

  get_pads_4;
#endif

#if !defined(Y_SAVE_TRIGS) || defined(Y_ITANIUM) || defined(Y_PRELOAD_TRIGS) || defined(Y_PRELOAD_DATA) || defined(Y_PREFETCH_EXPENSIVE)
# if defined(Y_ITANIUM)

  cplx_two_radix_4_dif_dit_init(s,ppdl,v,qqdl,pd);
  last = addr(n0<<1);
  for(i=n0;i<nc;i+=bigpad)
    {
      if((i+bigpad)<nc)
        nexti = i+ bigpad;
      else
        nexti = i + pad - 2;
      trig_4_load_init(tw);
      px = tw + Y_STEP;
      for(j=i;j<(pad+i);j+=2, px += 2*Y_STEP)
        {
          if((j+2)<(pad+i))
            next = j+2;
          else
            next = nexti;
          next = addr(next<<1);
          jj= next - last;
          last = next;
          get_twiddle_factors_4_preload;
          radix_4_twd_preload_B(jj);
        }
    }
# elif defined(Y_VECTORIZE2)
  for(i=n0;i<nc;i+=bigpad)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
      get_twiddle_factors_4;
      jj+=2;
      radix_4_twd_pp(B,jj);
      px += Y_STEP;
      j++;
      for(;j<(pad+i);j+=2,px += 2*Y_STEP)
        {
          jj=addr(j<<1);
          get_twiddle_factors_4_vector;
          radix_4_twd_pp_vector_B(jj);
        }
    }

# else
  for(i=n0;i<nc;i+=bigpad)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
#  if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4

      get_twiddle_factors_4;
      jj+=2;
      radix_4_twd_pp_no_fetch(B,jj);
      px += Y_STEP;
      j++;
#  if Y_CACHE_LINE == 8

      get_twiddle_factors_4;
      jj+=2;
      radix_4_twd_pp_no_fetch(B,jj);
      px += Y_STEP;
      get_twiddle_factors_4;
      jj+=2;
      radix_4_twd_pp_no_fetch(B,jj);
      px += Y_STEP;
      j+=2;
#   endif
#  endif

      for(;j<(pad+i);j++, px += Y_STEP)
        {
          jj=addr(j<<1);
          get_twiddle_factors_4;
          radix_4_twd_pp(B,jj);
#  if defined(Y_PREFETCH_EXPENSIVE) && Y_CACHE_LINE >= 4

          px += Y_STEP;
          j++;
          get_twiddle_factors_4;
          jj+=2;
          radix_4_twd_pp_no_fetch(B,jj);
#   if Y_CACHE_LINE == 8

          px += Y_STEP;
          get_twiddle_factors_4;
          jj+=2;
          radix_4_twd_pp_no_fetch(B,jj);
          px += Y_STEP;
          j+=2;
          get_twiddle_factors_4;
          jj+=2;
          radix_4_twd_pp_no_fetch(B,jj);
#   endif
#  endif

        }
    }
# endif
#else
  /*
  Y_SAVE_TRIGS defined

  This optimization tries to use the trig factor load/or computed in a cycle
  in the following cycle (bigpad ahead).

  It is no aplicable when defined Y_ITANIUM, Y_PRELOAD_TRIGS, Y_PRELOAD_DATA
  or Y_PREFETCH_EXPENSIVE
  */
  i=n0;
  i1=n0+bigpad;
# if defined(Y_VECTORIZE2)

  while(i1<nc)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
      jj1 = addr(i1<<1);
      radix_4_notwd_pp(B,jj1);
      j1 = i1 + 1;
      get_twiddle_factors_4;
      jj += 2;
      radix_4_twd_pp(B,jj);
      px += Y_STEP;
      j++;
      jj1 += 2;
      radix_4_twd_pp(B,jj1);
      j1++;
      for(;j<(pad+i);j+=2,j1+=2,px += 2*Y_STEP)
        {
          jj=addr(j<<1);
          get_twiddle_factors_4_vector;
          radix_4_twd_pp_vector_B(jj);
          jj1=addr(j1<<1);
          radix_4_twd_pp_vector_B(jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }
  while (i<nc)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
      get_twiddle_factors_4;
      jj+=2;
      radix_4_twd_pp(B,jj);
      px += Y_STEP;
      j++;
      for(;j<(pad+i);j+=2,px += 2*Y_STEP)
        {
          jj=addr(j<<1);
          get_twiddle_factors_4_vector;
          radix_4_twd_pp_vector_B(jj);
        }
      i += bigpad;
    }
# else
  while(i1<nc)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
      jj1 = addr(i1<<1);
      radix_4_notwd_pp(B,jj1);
      j1 = i1 + 1;
      for(;j<(pad+i);j++,j1++,px+=Y_STEP)
        {
          get_twiddle_factors_4;
          jj=addr(j<<1);
          radix_4_twd_pp(B,jj);
          jj1=addr(j1<<1);
          radix_4_twd_pp(B,jj1);
        }
      i+=(bigpad<<1);
      i1+=(bigpad<<1);
    }

  while(i<nc)
    {
      jj=addr(i<<1);
      radix_4_notwd_pp(B,jj);
      px=tw;
      j=i+1;
      for(;j<(pad+i);j++,px+=Y_STEP)
        {
          get_twiddle_factors_4;
          jj=addr(j<<1);
          radix_4_twd_pp(B,jj);
        }
      i += bigpad;
    }
# endif
#endif

}

/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r NO TWIDDLE
   INPUTS:
   d[] =all the data. Because of padding, mp_w[0] must be the first
           element of the whole data array. 
   tw[] = NULL in this no twiddle version. 
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radix_4_dit_notw(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0,
                      y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3;
  y_size_t i,j,jj,nc;

  ASSERT_ALIGNED_DOUBLE();


#ifndef NDEBUG1

  printf(" Radix_4_dit_notw n=%i n0=%i pad=%i \n",n,n0,pad);
#endif

  if(tw!= NULL)
    return;
  nc=n0+n;

  get_pads_4;

  for(i=n0;i<nc;i+=bigpad)
    {
      for(j=i;j<(pad+i);j++)
        {
          jj=addr(j<<1);
          radix_4_notwd_pp(B,jj);
        }
    }
}
/*$Id$*/
