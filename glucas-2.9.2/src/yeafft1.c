/*$Id$*/
/*
    This file is a part of 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ON_WORK
#define NDEBUG1
/******************************************************************
 yeafft is the include file for FFT. Read it to obtain more 
 information about it 
******************************************************************/
#include "yeafft.h"

/*
   Forward Fast Fourier Transfrom of w[] array of size n.
   This is an in-place decimation in frecuency.
*/
void y_fftf( y_ptr w, y_size_t n)
{
  int ir0, im, j, n0 = 0, nt;
  y_size_t ir = 0, pad, inc;
  y_ptr tw = NULL;
  HACK_ALIGN_STACK_ODD();
  /* first pass 0 */
  while ( Y_LRIGHT[ir] > Y_MEM_THRESHOLD)
    {
      n0 = 0;
      if(ir < (Y_NRADICES - 1))
        pad = Y_LRIGHT[ir + 1];
      else
        pad = 1;
      if(ir > 0)
        tw = &(Y_TWDF[ir - 1][0]);
      else
        tw = NULL;
      switch (Y_PLAN[ir])
        {
        case(8):
                if (ir == 0)
                  radix_8_dif_notw (w, tw, n, n0, pad);
            else
              radix_8_dif (w, tw, n, n0, pad);
          break;
        case(2):
                if (ir == 0)
                  radix_2_dif_notw (w, tw, n, n0, pad);
            else
              radix_2_dif (w, tw, n, n0, pad);
          break;
        case(4):
                if (ir == 0)
                  radix_4_dif_notw (w, tw, n, n0, pad);
            else
              radix_4_dif (w, tw, n, n0, pad);
          break;
#ifdef Y_SAVE

        case(5):
                radix_5_dif_notw (w, tw, n, n0, pad);
          break;
        case(6):
                radix_6_dif_notw (w, tw, n, n0, pad);
          break;
        case(7):
                radix_7_dif_notw (w, tw, n, n0, pad);
          break;
        case(9):
                radix_9_dif_notw (w, tw, n, n0, pad);
          break;
#else

        case(5):
                if (ir == 0)
                  radix_5_dif_notw (w, tw, n, n0, pad);
            else
              radix_5_dif (w, tw, n, n0, pad);
          break;
        case(6):
                if (ir == 0)
                  radix_6_dif_notw (w, tw, n, n0, pad);
            else
              radix_6_dif (w, tw, n, n0, pad);
          break;
        case(7):
                if (ir == 0)
                  radix_7_dif_notw (w, tw, n, n0, pad);
            else
              radix_7_dif(w, tw, n, n0, pad);
          break;
        case(9):
                if (ir == 0)
                  radix_9_dif_notw (w, tw, n, n0, pad);
            else
              radix_9_dif (w, tw, n, n0, pad);
          break;
#endif
#if Y_AVAL > 3

        case(16):
                if (ir == 0)
                  radix_16_dif_notw (w, tw, n, n0, pad);
            else
              radix_16_dif (w, tw, n, n0, pad);
          break;
# if Y_AVAL > 4

        case(32):
                if (ir == 0)
                  radix_32_dif_notw (w, tw, n, n0, pad);
            else
              radix_32_dif (w, tw, n, n0, pad);
          break;
# endif
# ifdef Y_SAVE

        case(10):
                radix_10_dif_notw (w, tw, n, n0, pad);
          break;
        case(12):
                radix_12_dif_notw (w, tw, n, n0, pad);
          break;
        case(14):
                radix_14_dif_notw (w, tw, n, n0, pad);
          break;
# else

        case(10):
                if (ir == 0)
                  radix_10_dif_notw (w, tw, n, n0, pad);
            else
              radix_10_dif (w, tw, n, n0, pad);
          break;
        case(12):
                if (ir == 0)
                  radix_12_dif_notw (w, tw, n, n0, pad);
            else
              radix_12_dif (w, tw, n, n0, pad);
          break;
        case(14):
                if (ir == 0)
                  radix_14_dif_notw (w, tw, n, n0, pad);
            else
              radix_14_dif (w, tw, n, n0, pad);
          break;
# endif
#endif

        default:
          fprintf (stderr,"Subroutine radix_%i_dif has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir++;
    }
  ir0 = ir;
  nt = Y_LRIGHT[ir];
  im = (n / nt);
  for(j = 0; j < im; j++, n0 += nt,ir = ir0)
{
      do
        {
          if(ir < (Y_NRADICES - 1))
            pad = Y_LRIGHT[ir + 1];
          else
            pad = 1;
          inc = Y_POWERS[Y_PLAN[ir] - 1] * 2;
          if(ir > 0)
            tw = &(Y_TWDF[ir - 1][inc * n0 / Y_LRIGHT[ir]]);
          else
            tw = NULL;
          switch (Y_PLAN[ir])
            {
            case(8):
                    if (ir == 0)
                      radix_8_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_8_dif (w, tw, nt, n0, pad);
              break;
            case(2):
                    if (ir == 0)
                      radix_2_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_2_dif (w, tw, nt, n0, pad);
              break;
            case(4):
                    if (ir == 0)
                      radix_4_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_4_dif (w, tw, nt, n0, pad);
              break;
#ifdef Y_SAVE

            case(5):
                    radix_5_dif_notw (w, tw, nt, n0, pad);
              break;
            case(6):
                    radix_6_dif_notw (w, tw, nt, n0, pad);
              break;
            case(7):
                    radix_7_dif_notw (w, tw, nt, n0, pad);
              break;
            case(9):
                    radix_9_dif_notw (w, tw, nt, n0, pad);
              break;
#else

            case(5):
                    if (ir == 0)
                      radix_5_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_5_dif (w, tw, nt, n0, pad);
              break;
            case(6):
                    if (ir == 0)
                      radix_6_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_6_dif (w, tw, nt, n0, pad);
              break;
            case(7):
                    if (ir == 0)
                      radix_7_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_7_dif (w, tw, nt, n0, pad);
              break;
            case(9):
                    if (ir == 0)
                      radix_9_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_9_dif (w, tw, nt, n0, pad);
              break;
#endif
#if Y_AVAL > 3

            case(16):
                    if (ir == 0)
                      radix_16_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_16_dif (w, tw, nt, n0, pad);
              break;
# if Y_AVAL > 4

            case(32):
                    if (ir == 0)
                      radix_32_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_32_dif (w, tw, nt, n0, pad);
              break;
# endif
# ifdef Y_SAVE

            case(10):
                    radix_10_dif_notw (w, tw, nt, n0, pad);
              break;
            case(12):
                    radix_12_dif_notw (w, tw, nt, n0, pad);
              break;
            case(14):
                    radix_14_dif_notw (w, tw, nt, n0, pad);
              break;
# else

            case(10):
                    if (ir == 0)
                      radix_10_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_10_dif (w, tw, nt, n0, pad);
              break;
            case(12):
                    if (ir == 0)
                      radix_12_dif_notw (w, tw, nt, n0, pad);
                else
                  radix_12_dif (w, tw, nt, n0, pad);
              break;
            case(14):
                    if (ir == 0)
                      radix_14_dif_notw(w, tw, nt, n0, pad);
                else
                  radix_14_dif (w, tw, nt, n0, pad);
              break;
# endif
#endif

            default:
              fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
                       Y_PLAN[ir]);
              exit(1);
            }
          ir++;
        }
      while (ir < Y_NRADICES);
    }
}

/*
   Backward Fast Fourier Transfrom of w[] array of size n.
   This is an in-place decimation in time.
*/
void y_fftb( y_ptr w, y_size_t n)
{
  int jm, j, n0 =0 , im, nt, ir = 0;
  y_size_t pad;
  y_ptr tw = NULL;
  HACK_ALIGN_STACK_ODD();

  jm = 0;
  while(Y_LRIGHT[jm] > Y_MEM_THRESHOLD)
    {
      jm++;
    };

  nt = Y_LRIGHT[jm];
  im = (n / nt);
  for(j = 0; j < im; j++, n0 += nt)
    {
      ir = Y_NRADICES - 1;
      do
        {
          if(ir < (int)(Y_NRADICES - 1))
            {
              pad = Y_LRIGHT[ir + 1];
              tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
            }
          else
            {
              pad = 1;
              tw = NULL;
            }
          switch (Y_PLAN[ir])
            {
            case(8):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_8_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_8_dit (w, tw, nt, n0, pad);
              break;
            case(2):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_2_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_2_dit (w, tw, nt, n0, pad);
              break;
            case(4):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_4_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_4_dit (w, tw, nt, n0, pad);
              break;
#ifdef Y_SAVE

            case(5):
                    radix_5_dit (w, tw, nt, n0, pad);
              break;
            case(6):
                    radix_6_dit (w, tw, nt, n0, pad);
              break;
            case(7):
                    radix_7_dit (w, tw, nt, n0, pad);
              break;
            case(9):
                    radix_9_dit (w, tw, nt, n0, pad);
              break;
#else

            case(5):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_5_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_5_dit (w, tw, nt, n0, pad);
              break;
            case(6):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_6_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_6_dit(w, tw, nt, n0, pad);
              break;
            case(7):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_7_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_7_dit (w, tw, nt, n0, pad);
              break;
            case(9):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_9_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_9_dit (w, tw, nt, n0, pad);
              break;
#endif
#if Y_AVAL > 3

            case(16):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_16_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_16_dit (w, tw, nt, n0, pad);
              break;

# if Y_AVAL > 4

            case(32):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_32_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_32_dit (w, tw, nt, n0, pad);
              break;
# endif
#ifdef Y_SAVE

            case(10):
                    radix_10_dit (w, tw, nt, n0, pad);
              break;
            case(12):
                    radix_12_dit (w, tw, nt, n0, pad);
              break;
            case(14):
                    radix_14_dit (w, tw, nt, n0, pad);
              break;
# else

            case(10):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_10_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_10_dit (w, tw, nt, n0, pad);
              break;
            case(12):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_12_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_12_dit (w, tw, nt, n0, pad);
              break;
            case(14):
                    if (ir == (int)(Y_NRADICES - 1))
                      radix_14_dit_notw (w, tw, nt, n0, pad);
                else
                  radix_14_dit (w, tw, nt, n0, pad);
              break;
# endif
#endif

            default:
              fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
                       Y_PLAN[ir]);
              exit(1);
            }
          ir--;
        }
      while (ir >= jm);
    }
  while (ir >= 0)
{
      n0 = 0;
      if(ir < (int)(Y_NRADICES - 1))
        {
          pad = Y_LRIGHT[ir+1];
          tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
        }
      else
        {
          pad = 1;
          tw = NULL;
        }
      switch (Y_PLAN[ir])
        {
        case(8):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_8_dit_notw( w, tw, n, n0, pad);
            else
              radix_8_dit( w, tw, n, n0, pad);
          break;
        case(2):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_2_dit_notw (w, tw, n, n0, pad);
            else
              radix_2_dit (w, tw, n, n0, pad);
          break;
        case(4):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_4_dit_notw (w, tw, n, n0, pad);
            else
              radix_4_dit (w, tw, n, n0, pad);
          break;
#ifdef Y_SAVE

        case(5):
                radix_5_dit (w, tw, n, n0, pad);
          break;
        case(6):
                radix_6_dit (w, tw, n, n0, pad);
          break;
        case(7):
                radix_7_dit (w, tw, n, n0, pad);
          break;
        case(9):
                radix_9_dit (w, tw, n, n0, pad);
          break;
#else

        case(5):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_5_dit_notw (w, tw, n, n0, pad);
            else
              radix_5_dit (w, tw, n, n0, pad);
          break;
        case(6):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_6_dit_notw (w, tw, n, n0, pad);
            else
              radix_6_dit (w, tw, n, n0, pad);
          break;
        case(7):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_7_dit_notw (w, tw, n, n0, pad);
            else
              radix_7_dit (w, tw, n, n0, pad);
          break;
        case(9):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_9_dit_notw (w, tw, n, n0, pad);
            else
              radix_9_dit(w, tw, n, n0, pad);
          break;
#endif
#if Y_AVAL > 3

        case(16):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_16_dit_notw(w, tw, n, n0, pad);
            else
              radix_16_dit(w, tw, n, n0, pad);
          break;

# if Y_AVAL > 4

        case(32):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_32_dit_notw (w, tw, n, n0, pad);
            else
              radix_32_dit (w, tw, n, n0, pad);
          break;
# endif

#ifdef Y_SAVE

        case(10):
                radix_10_dit (w, tw, n, n0, pad);
          break;
        case(12):
                radix_12_dit (w, tw, n, n0, pad);
          break;
        case(14):
                radix_14_dit (w, tw, n, n0, pad);
          break;
#else

        case(10):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_10_dit_notw (w, tw, n, n0, pad);
            else
              radix_10_dit (w, tw, n, n0, pad);
          break;
        case(12):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_12_dit_notw (w, tw, n, n0, pad);
            else
              radix_12_dit (w, tw, n, n0, pad);
          break;
        case(14):
                if (ir == (int)(Y_NRADICES - 1))
                  radix_14_dit_notw (w, tw, n, n0, pad);
            else
              radix_14_dit (w, tw, n, n0, pad);
          break;
# endif
#endif

        default:
          fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir--;
    }
}


int y_fftf_pass1(y_ptr w, y_size_t n)
{
  int ir = 0, n0 = 0;
  y_size_t pad;
  y_ptr tw;

  HACK_ALIGN_STACK_EVEN();
#ifndef NDEBUG1

  printf("y_fftf_pass1 size %i \n",n);
#endif

  if(Y_NRADICES == 1)
    return ir;

  while ( Y_LRIGHT[ir] > Y_MEM_THRESHOLD)
    {

      pad = Y_LRIGHT[ir + 1];
      if(ir > 0)
        tw = &(Y_TWDF[ir - 1][0]);
      else
        tw = NULL;
      switch (Y_PLAN[ir])
        {
        case(8):
                if (ir == 0)
                  radix_8_dif_notw (w, tw, n, n0, pad);
            else
              radix_8_dif (w, tw, n, n0, pad);
          break;
        case(4):
                if (ir == 0)
                  radix_4_dif_notw (w, tw, n, n0, pad);
            else
              radix_4_dif (w, tw, n, n0, pad);
          break;
#if Y_AVAL > 3

        case(16):
                if (ir == 0)
                  radix_16_dif_notw (w, tw, n, n0, pad);
            else
              radix_16_dif (w, tw, n, n0, pad);
          break;
# if Y_AVAL > 4

        case(32):
                if (ir == 0)
                  radix_32_dif_notw (w, tw, n, n0, pad);
            else
              radix_32_dif (w, tw, n, n0, pad);
          break;
# endif

# ifdef Y_SAVE

        case(10):
                radix_10_dif_notw (w, tw, n, n0, pad);
          break;
        case(12):
                radix_12_dif_notw (w, tw, n, n0, pad);
          break;
        case(14):
                radix_14_dif_notw (w, tw, n, n0, pad);
          break;
# else

        case(10):
                if (ir == 0)
                  radix_10_dif_notw (w, tw, n, n0, pad);
            else
              radix_10_dif (w, tw, n, n0, pad);
          break;
        case(12):
                if (ir == 0)
                  radix_12_dif_notw (w, tw, n, n0, pad);
            else
              radix_12_dif (w, tw, n, n0, pad);
          break;
        case(14):
                if (ir == 0)
                  radix_14_dif_notw (w, tw, n, n0, pad);
            else
              radix_14_dif (w, tw, n, n0, pad);
          break;
# endif
#endif

        case(2):
                if (ir == 0)
                  radix_2_dif_notw (w, tw, n, n0, pad);
            else
              radix_2_dif (w, tw, n, n0, pad);
          break;
#ifdef Y_SAVE

        case(5):
                radix_5_dif_notw (w, tw, n, n0, pad);
          break;
        case(6):
                radix_6_dif_notw (w, tw, n, n0, pad);
          break;
        case(7):
                radix_7_dif_notw (w, tw, n, n0, pad);
          break;
        case(9):
                radix_9_dif_notw (w, tw, n, n0, pad);
          break;
#else

        case(5):
                if (ir == 0)
                  radix_5_dif_notw (w, tw, n, n0, pad);
            else
              radix_5_dif (w, tw, n, n0, pad);
          break;
        case(6):
                if (ir == 0)
                  radix_6_dif_notw (w, tw, n, n0, pad);
            else
              radix_6_dif (w, tw, n, n0, pad);
          break;
        case(7):
                if (ir == 0)
                  radix_7_dif_notw (w, tw, n, n0, pad);
            else
              radix_7_dif (w, tw, n, n0, pad);
          break;
        case(9):
                if (ir == 0)
                  radix_9_dif_notw (w, tw, n, n0, pad);
            else
              radix_9_dif (w, tw, n, n0, pad);
          break;
#endif

        default:
          fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir++;
    }
  return ir;
}


void y_fftb_pass1(y_ptr w, y_size_t n, int irlast)
{
  int ir = irlast - 1, n0 = 0;
  y_size_t pad;
  y_ptr tw = NULL;

#ifndef NDEBUG1

  printf ("y_fftb_pass1 size %i irlast %i\n", n, irlast);
#endif

  HACK_ALIGN_STACK_EVEN();


  while (ir >= 0)
    {
      n0 = 0;
      pad = Y_LRIGHT[ir+1];
      tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
      switch (Y_PLAN[ir])
        {
        case(8):
                radix_8_dit (w, tw, n, n0, pad);
          break;
#if Y_AVAL > 3

        case(16):
                radix_16_dit (w, tw, n, n0, pad);
          break;
# if Y_AVAL > 4

        case(32):
                radix_32_dit (w, tw, n, n0, pad);
          break;
# endif
#endif

        case(5):
                radix_5_dit (w, tw, n, n0, pad);
          break;
        case(6):
                radix_6_dit (w, tw, n, n0, pad);
          break;
        case(7):
                radix_7_dit (w, tw, n, n0, pad);
          break;
        case(9):
                radix_9_dit (w, tw, n, n0, pad);
          break;
        case(10):
                radix_10_dit (w, tw, n, n0, pad);
          break;
        case(12):
                radix_12_dit (w, tw, n, n0, pad);
          break;
        case(14):
                radix_14_dit (w, tw, n, n0, pad);
          break;
        case(2):
                radix_2_dit (w, tw, n, n0, pad);
          break;
        case(4):
                radix_4_dit (w, tw, n, n0, pad);
          break;
        default:
          fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir--;
    }
}

void y_fftf_pass2(y_ptr w, y_size_t nt, y_size_t n0, int irfirst)
{
  int ir = irfirst;
  y_size_t pad, inc;
#if defined(Y_USE_SSE2)

  int r_mode, w_mode, rw;
#endif

  y_ptr tw;

  HACK_ALIGN_STACK_EVEN();
#ifndef NDEBUG1

  printf ("y_fftf_pass2 size %i offset %i irfirst %i\n", nt, n0, irfirst);
#endif

  if(irfirst == (int)(Y_NRADICES - 1))
    return; /* pass2 has no work to done */

  do
    {
#if defined(Y_USE_SSE2)
      /*r_mode = 0; w_mode = 0;*/
# if defined(ALL_INTERLACED)
      r_mode = (ir == 1)? 0 : 1;
      w_mode = (ir == (Y_NRADICES -2))? 0 : 1;
      rw = w_mode + 2 * r_mode;
# else

      rw = 0; /* legacy mode */
# endif
#endif

      pad = Y_LRIGHT[ir + 1];
      inc = Y_POWERS[Y_PLAN[ir] - 1] * 2;
      if(ir > 0)
        tw = &(Y_TWDF[ir - 1][inc * n0 / Y_LRIGHT[ir]]);
      else
        tw = NULL;
      switch (Y_PLAN[ir])
        {
        case(8):
                if (ir == 0)
                  radix_8_dif_notw (w, tw, nt, n0, pad);
#if defined(Y_USE_SSE2)

            else
              switch (rw)
            {
                case(0):
# if !defined(ALL_INTERLACED)
                        radixmm0_8_dif (w, tw, nt, n0, pad);
# else

                        radixmm2_8_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(1):
# if !defined(ALL_INTERLACED)
                        radixmm1_8_dif (w, tw, nt, n0, pad);
# else

                        radixmm3_8_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(2):
                        radixmm2_8_dif (w, tw, nt, n0, pad);
                  break;
                case(3):
                        radixmm3_8_dif (w, tw, nt, n0, pad);
                  break;
                }
#else
            else
              radix_8_dif (w, tw, nt, n0, pad);
#endif

          break;
        case(4):
                if (ir == 0)
                  radix_4_dif_notw (w, tw, nt, n0, pad);
#if defined(Y_USE_SSE2)

            else
              switch(rw)
            {
                case(0):
# if !defined(ALL_INTERLACED)
                        radixmm0_4_dif (w, tw, nt, n0, pad);
# else

                        radixmm2_4_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(1):
# if !defined(ALL_INTERLACED)
                        radixmm1_4_dif (w, tw, nt, n0, pad);
# else

                        radixmm3_4_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(2):
                        radixmm2_4_dif (w, tw, nt, n0, pad);
                  break;
                case(3):
                        radixmm3_4_dif (w, tw, nt, n0, pad);
                  break;
                }
#else
            else
              radix_4_dif (w, tw, nt, n0, pad);
#endif

          break;
#if Y_AVAL > 3

        case(16):
                if (ir == 0)
                  radix_16_dif_notw (w, tw, nt, n0, pad);
#if defined(Y_USE_SSE2)

            else
              switch(rw)
            {
                case(0):
# if !defined(ALL_INTERLACED)
                        radixmm0_16_dif (w, tw, nt, n0, pad);
# else

                        radixmm2_16_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(1):
# if !defined(ALL_INTERLACED)
                        radixmm1_16_dif (w, tw, nt, n0, pad);
# else

                        radixmm3_16_dif (w, tw, nt, n0, pad);
# endif

                  break;
                case(2):
                        radixmm2_16_dif (w, tw, nt, n0, pad);
                  break;
                case(3):
                        radixmm3_16_dif (w, tw, nt, n0, pad);
                  break;
                }
#else
            else
              radix_16_dif (w, tw, nt, n0, pad);
#endif

          break;
# if Y_AVAL > 4

        case(32):
                if (ir == 0)
                  radix_32_dif_notw (w, tw, nt, n0, pad);
            else
              radix_32_dif (w, tw, nt, n0, pad);
          break;
# endif
# ifdef Y_SAVE

        case(10):
                radix_10_dif_notw (w, tw, nt, n0, pad);
          break;
        case(12):
                radix_12_dif_notw (w, tw, nt, n0, pad);
          break;
        case(14):
                radix_14_dif_notw (w, tw, nt, n0, pad);
          break;
# else

        case(10):
                if (ir == 0)
                  radix_10_dif_notw (w, tw, nt, n0, pad);
            else
              radix_10_dif (w, tw, nt, n0, pad);
          break;
        case(12):
                if (ir == 0)
                  radix_12_dif_notw (w, tw, nt, n0, pad);
            else
              radix_12_dif (w, tw, nt, n0, pad);
          break;
        case(14):
                if (ir == 0)
                  radix_14_dif_notw (w, tw, nt, n0, pad);
            else
              radix_14_dif (w, tw, nt, n0, pad);
          break;
# endif
#endif
#ifdef Y_SAVE

        case(5):
                radix_5_dif_notw (w, tw, nt, n0, pad);
          break;
        case(6):
                radix_6_dif_notw (w, tw, nt, n0, pad);
          break;
        case(7):
                radix_7_dif_notw (w, tw, nt, n0, pad);
          break;
        case(9):
                radix_9_dif_notw (w, tw, nt, n0, pad);
          break;
#else

        case(5):
                if (ir == 0)
                  radix_5_dif_notw (w, tw, nt, n0, pad);
            else
              radix_5_dif (w, tw, nt, n0, pad);
          break;
        case(6):
                if (ir == 0)
                  radix_6_dif_notw (w, tw, nt, n0, pad);
            else
              radix_6_dif (w, tw, nt, n0, pad);
          break;
        case(7):
                if (ir == 0)
                  radix_7_dif_notw (w, tw, nt, n0, pad);
            else
              radix_7_dif (w, tw, nt, n0, pad);
          break;
        case(9):
                if (ir == 0)
                  radix_9_dif_notw (w, tw, nt, n0, pad);
            else
              radix_9_dif (w, tw, nt, n0, pad);
          break;
#endif

        case(2):
                if (ir == 0)
                  radix_2_dif_notw (w, tw, nt, n0, pad);
            else
              radix_2_dif (w, tw, nt, n0, pad);
          break;
        default:
          fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
                   Y_PLAN[ir]);
          exit(1);
        }
      ir++;
    }
  while (ir < (int)(Y_NRADICES - 1));
}

void y_fftb_pass2( y_ptr w, y_size_t nt, y_size_t n0, int irlast)
{
  int ir;
  y_size_t pad;
#if defined(Y_USE_SSE2)

  int r_mode, w_mode, rw;
#endif

  y_ptr tw = NULL;
  HACK_ALIGN_STACK_EVEN();
#ifndef NDEBUG1

  printf ("y_fftb_pass2 size %i offset %i irlast %i\n", nt, n0, irlast);
#endif

  if (irlast == (int)(Y_NRADICES - 1))
    return;  /* work already done */

  if (Y_NRADICES == 1)
    return;
  ir = Y_NRADICES - 2;
  do
    {
#if defined(Y_USE_SSE2)
      /*r_mode = 0; w_mode = 0;*/
# if defined(ALL_INTERLACED)
      r_mode = (ir == (Y_NRADICES - 2))? 0 : 1;
      w_mode = (ir == 1) ? 0 : 1;
      rw = w_mode + 2 * r_mode;
# else

      rw = 0; /* legacy mode */
# endif
#endif

      pad = Y_LRIGHT[ir + 1];
      tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
      switch (Y_PLAN[ir])
        {
        case(8):
#if defined(Y_USE_SSE2)
                switch(rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_8_dit (w, tw, nt, n0, pad);
# else

                      radixmm1_8_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_8_dit (w, tw, nt, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_8_dit (w, tw, nt, n0, pad);
# else

                      radixmm3_8_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_8_dit (w, tw, nt, n0, pad);
                break;
              }
#else
                radix_8_dit (w, tw, nt, n0, pad);
#endif

          break;
        case(4):
#if defined(Y_USE_SSE2)
                switch(rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_4_dit (w, tw, nt, n0, pad);
# else

                      radixmm1_4_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_4_dit (w, tw, nt, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_4_dit (w, tw, nt, n0, pad);
# else

                      radixmm3_4_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_4_dit (w, tw, nt, n0, pad);
                break;
              }
#else
                radix_4_dit (w, tw, nt, n0, pad);
#endif

          break;
#if Y_AVAL > 3

        case(16):
# if defined(Y_USE_SSE2)
                switch(rw)
              {
              case(0):
# if !defined(ALL_INTERLACED)
                      radixmm0_16_dit (w, tw, nt, n0, pad);
# else

                      radixmm1_16_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(1):
                      radixmm1_16_dit (w, tw, nt, n0, pad);
                break;
              case(2):
# if !defined(ALL_INTERLACED)
                      radixmm2_16_dit (w, tw, nt, n0, pad);
# else

                      radixmm3_16_dit (w, tw, nt, n0, pad);
# endif

                break;
              case(3):
                      radixmm3_16_dit (w, tw, nt, n0, pad);
                break;
              }
# else
                radix_16_dit (w, tw, nt, n0, pad);
# endif

          break;
# if Y_AVAL > 4

        case(32):
                radix_32_dit (w, tw, nt, n0, pad);
          break;
# endif

        case(10):
                radix_10_dit (w, tw, nt, n0, pad);
          break;
        case(12):
                radix_12_dit (w, tw, nt, n0, pad);
          break;
        case(14):
                radix_14_dit (w, tw, nt, n0, pad);
          break;
#endif

        case(5):
                radix_5_dit (w, tw, nt, n0, pad);
          break;
        case(6):
                radix_6_dit (w, tw, nt, n0, pad);
          break;
        case(7):
                radix_7_dit (w, tw, nt, n0, pad);
          break;
        case(9):
                radix_9_dit (w, tw, nt, n0, pad);
          break;
        case(2):
                radix_2_dit (w, tw, nt, n0, pad);
          break;
        default:
          fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
                   Y_PLAN[ir]);
          exit (1);
        }
      ir--;
    }
  while (ir >= irlast);
}

void y_fftf_mul_fftb(y_ptr w1, y_ptr w2, y_size_t nc)
{
  HACK_ALIGN_STACK_EVEN();
  switch (Y_PLAN[Y_NRADICES - 1])
    {
    case (8):
            radix_8_dif_mul_dit (w1, w2 , nc);
      break;
#if Y_AVAL > 3

    case (16):
            radix_16_dif_mul_dit (w1, w2 , nc);
      break;
#endif
#if Y_AVAL > 4

    case (32):
            radix_32_dif_mul_dit (w1, w2 , nc);
      break;
#endif

    case (2):
            radix_2_dif_mul_dit (w1, w2 , nc);
      break;
    case (4):
            radix_4_dif_mul_dit (w1, w2 , nc);
      break;
    default:
      fprintf (stderr, "Subroutine radix_%i_dif_mul_dit has not been writen yet\n",
               Y_PLAN[Y_NRADICES - 1]);
      exit(1);
    }
}

void y_fftf_square_fftb(y_ptr w, y_size_t nc)
{
  HACK_ALIGN_STACK_ODD();
  switch (Y_PLAN[Y_NRADICES - 1])
    {
#if defined(Y_USE_SSE2)
    case (8):
            radixmm_8_dif_square_dit (w, nc);
      break;
#else

    case (8):
            radix_8_dif_square_dit (w, nc);
      break;
#endif
#if Y_AVAL >3

    case (16):
            radix_16_dif_square_dit (w, nc);
      break;
#endif
#if Y_AVAL >4

    case (32):
            radix_32_dif_square_dit (w, nc);
      break;
#endif

    case (2):
            radix_2_dif_square_dit (w, nc);
      break;
#if defined(Y_USE_SSE2)

    case (4):
            radixmm_4_dif_square_dit (w, nc);
      break;
#else

    case (4):
            radix_4_dif_square_dit (w, nc);
      break;
#endif

    default:
      fprintf (stderr, "Subroutine radix_%i_dif_square_dit has not been writen yet\n",
               Y_PLAN[Y_NRADICES-1]);
      exit(1);
    }
}

void y_fftf_mul_fftb_block(y_ptr w1, y_ptr w2, y_size_t bs, y_size_t i,
                           y_size_t j)
{
  HACK_ALIGN_STACK_EVEN();
  switch (Y_PLAN[Y_NRADICES - 1])
    {
    case (8):
            radix_8_dif_mul_dit_block (w1, w2, bs, i, j);
      break;
#if Y_AVAL > 3

    case (16):
            radix_16_dif_mul_dit_block (w1, w2, bs, i, j);
      break;
#endif
#if Y_AVAL > 4

    case (32):
            radix_32_dif_mul_dit_block (w1, w2, bs, i, j);
      break;
#endif

    case (2):
            radix_2_dif_mul_dit_block (w1, w2, bs, i, j);
      break;
    case (4):
            radix_4_dif_mul_dit_block (w1, w2, bs, i, j);
      break;
    default:
      fprintf (stderr, "Subroutine radix_%i_fftf_mul_fftb_block has not been writen yet\n",
               Y_PLAN[Y_NRADICES - 1]);
      exit(1);
    }
}



void y_fftf_square_fftb_block(y_ptr w, y_size_t bs, y_size_t i, y_size_t j)
{
  HACK_ALIGN_STACK_ODD();
  switch (Y_PLAN[Y_NRADICES - 1])
    {
#if defined(Y_USE_SSE2)
    case (8):
            radixmm_8_dif_square_dit_block (w, bs, i, j);
      break;
#else

    case (8):
            radix_8_dif_square_dit_block (w, bs, i, j);
      break;
#endif
#if Y_AVAL > 3

    case (16):
            radix_16_dif_square_dit_block (w, bs, i, j);
      break;
#endif
#if Y_AVAL > 4

    case (32):
            radix_32_dif_square_dit_block (w, bs, i, j);
      break;
#endif

    case (2):
            radix_2_dif_square_dit_block (w, bs, i, j);
      break;
#if defined(Y_USE_SSE2)

    case (4):
            radixmm_4_dif_square_dit_block (w, bs, i, j);
      break;
#else

    case (4):
            radix_4_dif_square_dit_block (w, bs, i, j);
      break;
#endif

    default:
      fprintf (stderr, "Subroutine radix_%i_fftf_mul_fftb_block has not been writen yet\n",
               Y_PLAN[Y_NRADICES - 1]);
      exit(1);
    }
}

void y_fftf_mul_fftb_pass2(y_ptr w1, y_ptr w2, int irfirst )
{
  int ii, jj, bs, bi, bj, nr;
  y_size_t nc;

  bs = Y_LRIGHT[irfirst];
  nc = Y_NRADICES - 1 - irfirst;

  /* first block */
  y_fftf_pass2 (w2, bs, 0, irfirst);
  y_fftf_pass2 (w1, bs, 0, irfirst);
  y_fftf_mul_fftb (w1, w2, nc);
  y_fftb_pass2 (w1, bs, 0, irfirst);
  if( irfirst == 0)
    return;

  for (nr = irfirst - 1; nr >= 0; nr--)
    {
      bi = (Y_LRIGHT[nr + 1] / bs);
      bj = (Y_LRIGHT[nr] / bs) - 1;
      for(ii = bi, jj = bj; ii < jj; ii++, jj--)
        {
          y_fftf_pass2 (w2, bs, ii * bs, irfirst);
          y_fftf_pass2 (w2, bs, jj * bs, irfirst);
          y_fftf_pass2 (w1, bs, ii * bs, irfirst);
          y_fftf_pass2 (w1, bs, jj * bs, irfirst);
          y_fftf_mul_fftb_block (w1, w2, bs, ii, jj);
          y_fftb_pass2 (w1, bs, ii * bs, irfirst);
          y_fftb_pass2 (w1, bs, jj * bs, irfirst);
        }
      if(ii == jj)
        {
          y_fftf_pass2 (w2, bs, ii * bs, irfirst);
          y_fftf_pass2 (w1, bs, ii * bs, irfirst);
          y_fftf_mul_fftb_block (w1, w2, bs, ii, ii);
          y_fftb_pass2 (w1, bs, ii * bs, irfirst);
        }
    }
}


#if (!defined(_OPENMP) && !defined(_SUNMP) && !defined(_PTHREADS))

void y_fftf_squar_fftb_pass2(y_ptr w, int irfirst )
{
  int ii, jj, bi, bj, bs, nr;
  y_size_t nc;

  bs = Y_LRIGHT[irfirst];
  nc = Y_NRADICES - 1 - irfirst;


  y_fftf_pass2 (w, bs, 0, irfirst);
  y_fftf_square_fftb (w,nc);
  y_fftb_pass2 (w, bs, 0, irfirst);

  if( irfirst == 0 )
    return;

  for (nr = irfirst - 1; nr >= 0; nr--)
    {
      bi = (Y_LRIGHT[nr + 1] / bs);
      bj = (Y_LRIGHT[nr] / bs) - 1;
      for(ii = bi, jj = bj; ii < jj; ii++, jj--)
        {
          y_fftf_pass2 (w, bs, ii * bs, irfirst);
          y_fftf_pass2 (w, bs, jj * bs, irfirst);
          y_fftf_square_fftb_block (w, bs, ii, jj);
          y_fftb_pass2 (w, bs, ii * bs, irfirst);
          y_fftb_pass2 (w, bs, jj * bs, irfirst);
        }
      if(ii == jj)
        {
          y_fftf_pass2 (w, bs, ii * bs, irfirst);
          y_fftf_square_fftb_block (w, bs, ii, ii);
          y_fftb_pass2 (w, bs, ii * bs, irfirst);
        }
    }
}

#else

void y_pass2_mp(y_ptr w, int bs, int *di, int *dj, int nr, int nc, int ifir,
                int ni)
{
  int i;
  for(i = ifir; i < ifir + ni; i++)
    {
      if(i == 0)
        {
          y_fftf_pass2 (w, bs, 0, nr);
          y_fftf_square_fftb (w, nc);
          y_fftb_pass2 (w, bs, 0, nr);
        }
      else if(di[i] == dj[i])
        {
          y_fftf_pass2 (w, bs, di[i] * bs, nr);
          y_fftf_square_fftb_block (w, bs, di[i], di[i]);
          y_fftb_pass2 (w, bs, di[i] * bs, nr);
        }
      else
        {
          y_fftf_pass2 (w, bs, di[i] * bs, nr);
          y_fftf_pass2 (w, bs, dj[i] * bs, nr);
          y_fftf_square_fftb_block (w, bs, di[i], dj[i]);
          y_fftb_pass2 (w, bs, di[i] * bs, nr);
          y_fftb_pass2 (w, bs, dj[i] * bs, nr);
        }
    }

}


/* Openmp multithread version */
void y_fftf_squar_fftb_pass2(y_ptr w, int irfirst )
{
  int ii, bs;
  y_size_t nc;

  bs = Y_LRIGHT[irfirst];
  nc = Y_NRADICES - 1 - irfirst;

  /* if irfirst == 0 there is no multithread */
  if( irfirst == 0)
    {
      y_fftf_pass2 (w, bs, 0, 0);
      y_fftf_square_fftb (w,nc);
      y_fftb_pass2 (w, bs, 0, 0);
      return;
    }

  /*  if(Y_DI == NULL) y_threads_init_pass2(irfirst);*/

#if defined(_OPENMP)
# pragma omp parallel for default(shared)
#elif defined(_SUNMP)
# pragma MP taskloop maxcpus(Y_NUM_THREADS)
# pragma MP taskloop shared(w, Y_BS, Y_DI, Y_DJ, Y_IR, Y_NC, Y_2N0, Y_2NN)
#endif
  for (ii = 0; ii < Y_T2; ii++)
    {
      y_pass2_mp( w, Y_BS, Y_DI, Y_DJ, Y_IR, Y_NC, Y_2N0[ii], Y_2NN[ii]);
    }
}

#endif

void y_convolution(y_ptr w1, y_ptr w2 , y_size_t n)
{
  int ir;

  ir = y_fftf_pass1 (w1, n);
  y_fftf_pass1 (w2, n);
  y_fftf_mul_fftb_pass2 (w1, w2, ir);
  y_fftb_pass1 (w1, n, ir);
}

void y_auto_convolution(y_ptr w, y_size_t n)
{
  int ir;

  ir = y_fftf_pass1 (w, n);
  y_fftf_squar_fftb_pass2 (w, ir);
  y_fftb_pass1 (w, n, ir);
}

/*$Id$*/

