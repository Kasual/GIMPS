/*$Id$*/
/*
    yealucasmp.c. An interface to use YEAFFT in a Lucas Lehmer test
 
    Copyright (C) 2004-2006  Guillermo Ballester Valor
 
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

/*
   version for multithreaded pass1 
   For Pass1 with only one radix reduction FFT step, there is no problem.
 
   But if there two or more FFT steps then we have to be careful. As example, 
   assuming two radix-4 steps, a thread could be use the memory marked by
   X, other the Y and other the Z
 
   XXXY XXXY XXXY XXXY YYYY YYYY YYYY YYYY YYZZ YYZZ YYZZ YYZZ .....
 
   It ilustrate the fact that for the inner passes there can be many 
   non-contiguous blocks. It may occurr for ir > 1.
   
   THE CONCLUSION is we only can make a single radix pass. If PASS1 has 
   more than 1 radix reduction we need to make several steps (a radix reduction
   per step)
*/
int y_fftf_pass1_lucas_mp(y_ptr w, y_size_t n, y_size_t n0, int ir)
{
  y_size_t pad;
  y_ptr tw;

  HACK_ALIGN_STACK_EVEN();

#ifndef NDEBUG1

  printf ("y_fftf_pass1_lucas_mp size %i \n", n);
#endif

  if(Y_NRADICES <= 2)
    return ir;
  if (ir >= Y_PASS2)
    return ir;

  pad = Y_LRIGHT[ir + 1];
  tw = &(Y_TWDF[ir - 1][0]);
  switch (Y_PLAN[ir])
    {
    case(8):
            radixmp_8_dif (w, tw, n, n0, pad);
      break;
#if Y_AVAL > 3

    case(16):
            radixmp_16_dif (w, tw, n, n0, pad);
      break;
# if Y_AVAL > 4

    case(32):
            radixmp_32_dif (w, tw, n, n0, pad);
      break;
# endif
#endif

    case(4):
            radixmp_4_dif (w, tw, n, n0, pad);
      break;
    default:
      fprintf (stderr, "Subroutine radix_%i_dif has not been writen yet\n",
               Y_PLAN[ir]);
      exit(1);
    }
  return ir;
}



/* version for multithreaded pass1 */
void y_fftb_pass1_lucas_mp(y_ptr w, y_size_t n, y_size_t n0, int ir)
{
  y_size_t pad;
  y_ptr tw = NULL;

#ifndef NDEBUG1

  printf ("y_fftb_pass1_lucas_mp size %i irlast %i\n", n, ir);
#endif

  HACK_ALIGN_STACK_EVEN();

  pad = Y_LRIGHT[ir + 1];
  tw = &(Y_TWDB[Y_NRADICES - 2 - ir][0]);
  switch (Y_PLAN[ir])
    {
    case(8):
            radixmp_8_dit (w, tw, n, n0, pad);
      break;
#if Y_AVAL > 3

    case(16):
            radixmp_16_dit (w, tw, n, n0, pad);
      break;
# if Y_AVAL > 4

    case(32):
            radixmp_32_dit (w, tw, n, n0, pad);
      break;
# endif
#endif

    case(4):
            radixmp_4_dit (w, tw, n, n0, pad);
      break;
    default:
      fprintf (stderr, "Subroutine radix_%i_dit has not been writen yet\n",
               Y_PLAN[ir]);
      exit(1);
    }
}


