/*$Id$*/
/*
    yealucas.c. An interface to use YEAFFT in a Lucas Lehmer test
 
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NDEBUG1
/******************************************************************
 yeafft is the include file for FFT. Read it to obtain more 
 information about it 
******************************************************************/
#include "yeafft.h"
#include "glucas.h"
void y_fftf_pass1_lucas_first(y_ptr w, y_size_t n)
{
  int n0=0;
  y_size_t pad;
  y_ptr tw;

  HACK_ALIGN_STACK_EVEN();
#ifndef NDEBUG1

  printf("y_fftf_pass1_lucas_first size %i \n",n);
#endif

  if(Y_NRADICES == 1)
    return;

  /* set next Y_SBIT */
  Y_SBIT += Y_SBIT;
  if(Y_SBIT >= Y_EXPONENT )
    Y_SBIT -= Y_EXPONENT;

  pad=Y_LRIGHT[1];
  tw=NULL;
  switch (Y_PLAN[0])
    {
    case(2):
            radix_2_dif_notw (w, tw, n, n0, pad);
      break;
    case(4):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_4_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_4_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_4_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(8):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_8_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_8_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_8_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(5):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_5_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_5_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_5_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(6):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_6_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_6_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_6_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(7):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_7_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_7_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_7_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(9):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_9_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_9_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_9_dif_notw (w, tw, n, n0, pad);
#endif

      break;
#if Y_AVAL > 3

    case(16):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_16_dif_notw (w, tw, n, n0, pad);
# else

            radixmm1_16_dif_notw (w, tw, n, n0, pad);
# endif
#else

            radix_16_dif_notw (w, tw, n, n0, pad);
#endif

      break;
    case(10):
            radix_10_dif_notw (w, tw, n, n0, pad);
      break;
    case(12):
            radix_12_dif_notw (w, tw, n, n0, pad);
      break;
    case(14):
            radix_14_dif_notw (w, tw, n, n0, pad);
      break;
#endif

    default:
      fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
               Y_PLAN[0]);
      exit(1);
    }
}

int y_fftf_pass1_lucas(y_ptr w, y_size_t n, y_size_t n0)
{
  int ir = 1, inc;
  y_size_t pad;
#if defined(Y_USE_SSE2)

  int r_mode, w_mode, rw;
#endif

  y_ptr tw;

  HACK_ALIGN_STACK_EVEN();

#ifndef NDEBUG1

  printf ("y_fftf_pass1_lucas size %i \n", n);
#endif

  if(Y_NRADICES <= 2)
    return ir;

  while ( Y_LRIGHT[ir] > Y_MEM_THRESHOLD)
    {
#if defined(Y_USE_SSE2)
# if defined(ALL_INTERLACED)
      /*w_mode = 0; r_mode = 0;*/
      w_mode = (ir == (Y_NRADICES - 2)) ? 0 : 1;
      r_mode = (ir == 1 )? 0 : 1;
      rw = w_mode + 2 * r_mode;
# else

      rw = 0; /* legacy mode */
# endif
#endif

      pad = Y_LRIGHT[ir + 1];
      inc = Y_POWERS[Y_PLAN[ir] - 1] * 2;
      tw = &(Y_TWDF[ir - 1][inc * n0 / Y_LRIGHT[ir]]);
      switch (Y_PLAN[ir])
        {
        case(8):
#if defined(Y_USE_SSE2)
                switch(rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_8_dif (w, tw, n, n0, pad);
# else

                      radixmm2_8_dif (w, tw, n, n0, pad);
# endif

                break;
              case(1):
# if !defined(ALL_INTERLACED)
                      radixmm1_8_dif (w, tw, n, n0, pad);
# else

                      radixmm3_8_dif (w, tw, n, n0, pad);
# endif

                break;
              case(2):
                      radixmm2_8_dif (w, tw, n, n0, pad);
                break;
              case(3):
                      radixmm3_8_dif (w, tw, n, n0, pad);
                break;
              }
#else
                radix_8_dif (w, tw, n, n0, pad);
#endif

          break;
#if Y_AVAL > 3

        case(16):
#if defined(Y_USE_SSE2)
                switch (rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_16_dif (w, tw, n, n0, pad);
# else

                      radixmm2_16_dif (w, tw, n, n0, pad);
# endif

                break;
              case(1):
# if !defined(ALL_INTERLACED)
                      radixmm1_16_dif (w, tw, n, n0, pad);
# else

                      radixmm3_16_dif (w, tw, n, n0, pad);
# endif

                break;
              case(2):
                      radixmm2_16_dif (w, tw, n, n0, pad);
                break;
              case(3):
                      radixmm3_16_dif (w, tw, n, n0, pad);
                break;
              }
#else
                radix_16_dif (w, tw, n, n0, pad);
#endif

          break;
# if Y_AVAL > 4

        case(32):
                radix_32_dif (w, tw, n, n0, pad);
          break;
# endif
#endif

        case(4):
#if defined(Y_USE_SSE2)
                switch (rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_4_dif (w, tw, n, n0, pad);
# else

                      radixmm2_4_dif (w, tw, n, n0, pad);
# endif

                break;
              case(1):
# if !defined(ALL_INTERLACED)
                      radixmm1_4_dif (w, tw, n, n0, pad);
# else

                      radixmm3_4_dif (w, tw, n, n0, pad);
# endif

                break;
              case(2):
                      radixmm2_4_dif (w, tw, n, n0, pad);
                break;
              case(3):
                      radixmm3_4_dif (w, tw, n, n0, pad);
                break;
              }
#else
                radix_4_dif (w, tw, n, n0, pad);
#endif

          break;
        case(2):
                radix_2_dif (w, tw, n, n0, pad);
          break;
        default:
          fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir++;
    }
  return ir;
}


void y_fftb_pass1_lucas(y_ptr w, y_size_t n, y_size_t n0, int irlast)
{
  int ir = irlast - 1;
  y_size_t pad;
#if defined(Y_USE_SSE2)

  int r_mode, w_mode, rw;
#endif

  y_ptr tw = NULL;

#ifndef NDEBUG1

  printf ("y_fftb_pass1_lucas size %i irlast %i\n", n, irlast);
#endif

  HACK_ALIGN_STACK_EVEN();
  while (ir >= 1)
    {
#if defined(Y_USE_SSE2)
# if defined(ALL_INTERLACED)
      /*w_mode = 0; r_mode = 0;*/
      r_mode = (ir == (Y_NRADICES - 2)) ? 0 : 1;
      w_mode = (ir == 1 ) ? 0 : 1;
      rw = w_mode + 2 * r_mode;
# else

      rw = 0;
# endif
#endif

      pad = Y_LRIGHT[ir + 1];
      tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
      switch (Y_PLAN[ir])
        {
        case(8):
#if defined(Y_USE_SSE2)
                switch (rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_8_dit (w, tw, n, n0, pad);
# else

                      radixmm1_8_dit (w, tw, n, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_8_dit (w, tw, n, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_8_dit (w, tw, n, n0, pad);
# else

                      radixmm3_8_dit (w, tw, n, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_8_dit (w, tw, n, n0, pad);
                break;
              }
#else
                radix_8_dit (w, tw, n, n0, pad);
#endif

          break;
#if Y_AVAL > 3

        case(16):
#if defined(Y_USE_SSE2)
                switch (rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_16_dit (w, tw, n, n0, pad);
# else

                      radixmm1_16_dit (w, tw, n, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_16_dit (w, tw, n, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_16_dit (w, tw, n, n0, pad);
# else

                      radixmm3_16_dit (w, tw, n, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_16_dit (w, tw, n, n0, pad);
                break;
              }
#else
                radix_16_dit (w, tw, n, n0, pad);
#endif

          break;
# if Y_AVAL > 4

        case(32):
                radix_32_dit (w, tw, n, n0, pad);
          break;
# endif
#endif

        case(4):
#if defined(Y_USE_SSE2)
                switch (rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_4_dit (w, tw, n, n0, pad);
# else

                      radixmm1_4_dit (w, tw, n, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_4_dit (w, tw, n, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_4_dit (w, tw, n, n0, pad);
# else

                      radixmm3_4_dit (w, tw, n, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_4_dit (w, tw, n, n0, pad);
                break;
              }
#else
                radix_4_dit (w, tw, n, n0, pad);
#endif

          break;
        case(2):
                radix_2_dit (w, tw, n, n0, pad);
          break;
        default:
          fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir--;
    }
}

void lucas_dit_carry_norm_dif(y_ptr w, UL N, UL err_flag)
{
  y_size_t n0 = 0, pad = Y_LRIGHT[1], n = (y_size_t)(N) >> 1;
  y_ptr tw = &(Y_TWDB[Y_NRADICES - 2][0]);

#ifndef NDEBUG1

  printf ("lucas_dit_carry_norm_dif  Y_PLAN[0]=%i \n", Y_PLAN[0]);
#endif

  HACK_ALIGN_STACK_EVEN();

  switch(Y_PLAN[0])
    {
    case(2):
            radix_2_dit (w, tw, n, n0, pad);
      HACK_ALIGN_STACK_EVEN();
      normalize (w, N,err_flag);
      tw = NULL;
      HACK_ALIGN_STACK_EVEN();
      radix_2_dif_notw (w, tw, n, n0, pad);
      break;
    case(4):
            dit_carry_norm_dif_4 (w, N, err_flag);
      break;
    case(5):
#if defined(Y_USE_SSE2)
            dit_carry_norm_dif_5_sse2 (w, N, err_flag);
#else

            dit_carry_norm_dif_5 (w, N, err_flag);
#endif

      break;
    case(6):
#if defined(Y_USE_SSE2)
            dit_carry_norm_dif_6_sse2 (w, N, err_flag);
#else

            dit_carry_norm_dif_6 (w, N, err_flag);
#endif

      break;
    case(7):
#if defined(Y_USE_SSE2)
            dit_carry_norm_dif_7_sse2 (w, N, err_flag);
#else

            dit_carry_norm_dif_7 (w, N, err_flag);
#endif

      break;
    case(8):
#if defined(Y_USE_SSE2)
            dit_carry_norm_dif_8_sse2 (w, N, err_flag);
#else

            dit_carry_norm_dif_8 (w, N, err_flag);
#endif

      break;
    case(9):
#if defined(Y_USE_SSE2)
            dit_carry_norm_dif_9_sse2 (w, N, err_flag);
#else

            dit_carry_norm_dif_9 (w, N, err_flag);
#endif

      break;
#if Y_AVAL > 3

    case(10):
            dit_carry_norm_dif_10 (w, N, err_flag);
      break;
    case(12):
            dit_carry_norm_dif_12 (w, N, err_flag);
      break;
    case(14):
            dit_carry_norm_dif_14 (w, N, err_flag);
      break;
    case(16):
            dit_carry_norm_dif_16 (w, N, err_flag);
      break;
#endif

    }
}

void y_fftb_pass1_lucas_last(y_ptr w, y_size_t n)
{
  y_size_t n0 = 0,pad = Y_LRIGHT[1];
  y_ptr tw = &(Y_TWDB[Y_NRADICES - 2][0]);

  HACK_ALIGN_STACK_EVEN();

#ifndef NDEBUG1

  printf ("y_fftb_pass1_lucas_last  Y_PLAN[0]=%i \n", Y_PLAN[0]);
#endif

  switch(Y_PLAN[0])
    {
    case(2):
            radix_2_dit (w, tw, n, n0, pad);
      break;
    case(4):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_4_dit (w, tw, n, n0, pad);
# else

            radixmm2_4_dit (w, tw, n, n0, pad);
# endif
#else

            radix_4_dit (w, tw, n, n0, pad);
#endif

      break;
    case(5):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_5_dit (w, tw, n, n0, pad);
# else

            radixmm2_5_dit (w, tw, n, n0, pad);
# endif
#else

            radix_5_dit (w, tw, n, n0, pad);
#endif

      break;
    case(6):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_6_dit (w, tw, n, n0, pad);
# else

            radixmm2_6_dit (w, tw, n, n0, pad);
# endif
#else

            radix_6_dit (w, tw, n, n0, pad);
#endif

      break;
    case(7):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_7_dit (w, tw, n, n0, pad);
# else

            radixmm2_7_dit (w, tw, n, n0, pad);
# endif
#else

            radix_7_dit (w, tw, n, n0, pad);
#endif

      break;
    case(8):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_8_dit (w, tw, n, n0, pad);
# else

            radixmm2_8_dit (w, tw, n, n0, pad);
# endif
#else

            radix_8_dit (w, tw, n, n0, pad);
#endif

      break;
    case(9):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_9_dit (w, tw, n, n0, pad);
# else

            radixmm2_9_dit (w, tw, n, n0, pad);
# endif
#else

            radix_9_dit (w, tw, n, n0, pad);
#endif

      break;
#if Y_AVAL > 3

    case(10):
            radix_10_dit (w, tw, n, n0, pad);
      break;
    case(12):
            radix_12_dit (w, tw, n, n0, pad);
      break;
    case(14):
            radix_14_dit (w, tw, n, n0, pad);
      break;
    case(16):
#if defined(Y_USE_SSE2)
# if !defined(ALL_INTERLACED)
            radixmm0_16_dit (w, tw, n, n0, pad);
# else

            radixmm2_16_dit (w, tw, n, n0, pad);
# endif
#else

            radix_16_dit (w, tw, n, n0, pad);
#endif

      break;
#endif

    }
}


void lucas_first_auto_convolution(y_ptr w, y_size_t n)
{
  int ir;

  y_fftf_pass1_lucas_first (w, n);
  ir= y_fftf_pass1_lucas (w, n, 0);
  y_fftf_squar_fftb_pass2 (w, ir);
  y_fftb_pass1_lucas (w, n, 0, ir);
}

void lucas_auto_convolution(y_ptr w, y_size_t n)
{
  int ir;

#if defined(_OPENMP) || defined(_SUNMP)

  int n0, nthr;

  nthr = Y_NTHREADS;

  if (Y_PASS2 > 1)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared) lastprivate(ir)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(ir)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          ir = y_fftf_pass1_lucas_mp(w,Y_1NC[n0],Y_1NN[n0], 1);
        }
    }
  if (Y_PASS2 > 2)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared) lastprivate(ir)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(ir)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          ir = y_fftf_pass1_lucas_mp(w,Y_1NC[n0],Y_1NN[n0], 2);
        }

    }

  if (Y_PASS2 > 3)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared) lastprivate(ir)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(ir)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          ir = y_fftf_pass1_lucas_mp(w,Y_1NC[n0],Y_1NN[n0], 3);
        }

    }

  y_fftf_squar_fftb_pass2 (w, Y_PASS2);

  if (Y_PASS2 > 3)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          y_fftb_pass1_lucas_mp(w, Y_1NC[n0], Y_1NN[n0], 3);
        }

    }

  if (Y_PASS2 > 2)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          y_fftb_pass1_lucas_mp(w, Y_1NC[n0], Y_1NN[n0], 2);
        }

    }

  if (Y_PASS2 > 1)
    {
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(w, Y_1NC, Y_1NN)
# endif
      for (n0 = 0; n0 < Y_T1; n0++)
        {
          y_fftb_pass1_lucas_mp(w, Y_1NC[n0], Y_1NN[n0], 1);
        }
    }
#else

  ir = y_fftf_pass1_lucas (w, n, 0);
  y_fftf_squar_fftb_pass2 (w, ir);
  y_fftb_pass1_lucas (w, n, 0, ir);

#endif


}

void lucas_last_auto_convolution(y_ptr w, y_size_t n)
{
  int ir;

  ir = y_fftf_pass1_lucas (w, n, 0);
  y_fftf_squar_fftb_pass2 (w, ir);
  y_fftb_pass1 (w, n, ir);
}
/*$Id$*/
