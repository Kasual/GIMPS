/*$Id$*/
/*  This file is a part of 
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor, Klaus Kastens

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
#include <memory.h>
#include <math.h>
#include "yeafft.h"
#include "glucas.h"

#if defined(TRICKY_ROUND)
# define TEST_LIMIT 29
/*# define DEBUG_ROUND*/
# if(defined(__INTEL_COMPILER) && defined(__ICC))
#  define TEST_RINT(_x) ((long long) (_x))
# else
#  define TEST_RINT(_x) (((_x) + A ) - B)

/*
#  define TEST_RINT(x) ({ double __arg ;\
   __asm__ volatile ("frndint \n" : "=t" (__arg) : "0" (x));\
   __arg;})
*/
# endif

BIG_DOUBLE tricky(void)
{
  BIG_DOUBLE A=3.0, B=3.0, x, y, z, r1, r2;
  long int i, j, l, fail;

  l = 8 * sizeof(BIG_DOUBLE);
  i = 0;
  fail = 0;
  x = 1.3;
  x = 3.3;
  r1 = floor(x + 0.5);
  do
    {
      A *= 2.0;
      B *= 2.0;
      i++;
      r2 = TEST_RINT(x);
#ifdef DEBUG_ROUND	      
      printf("i=%d A=%lf B=%lf r1=%lf r2=%lf\n", i, A, 
	     B, r1, r2);
#endif      
    }while ((r1 != r2) && (i < l));

  /* 
     problems with new gcc-ssa compiler and ia32 , it is 
     better to  set directly the values 
  */
# if (defined(__i386__) || defined(__i486__) || defined(__i586__) || defined(__i686__))
  A = ((3.0 *(1L << 31)) * (1L << 31));
  B = A;
#endif


  
  if(i < l)
    {
      x=1.0; y=1.263885; z=.723547123;
      /*for(j = 1; j < (2L << TEST_LIMIT); j <<= 1)
	{*/
	  for (x = -0.498; x < 0.5; x += 0.002)
	    {
	      r1 = floor(x + 0.5);
	      r2 = TEST_RINT(x);
#ifdef DEBUG_ROUND	      
	      printf("test 1: x=%lf r1=%lf r2=%lf\n", x, 
		     r1, r2);
#endif
	      if(r1 != r2) fail = 1;
	      r1 = floor(-x + 0.5);
	      r2 = TEST_RINT(-x);
#ifdef DEBUG_ROUND	      
	      printf("test 2: -x=%lf r1=%lf r2=%lf\n", -x, 
		     r1, r2);
#endif
	      if (r1 != r2) fail = 1;
	      r1 = floor( x * y + 0.5);
	      r2 = TEST_RINT(x * y);
#ifdef DEBUG_ROUND	      
	      printf("test 3: x+y=%lf r1=%lf r2=%lf\n", x + y, 
		     r1, r2);
#endif
	      if (r1 != r2) fail = 1;
	      r1 = floor(x * y + z + 0.5);
	      r2 = TEST_RINT ( x * y + z);
#ifdef DEBUG_ROUND	      
	      printf("test 4: x+y+z=%lf r1=%lf r2=%lf\n", x + y + z, 
		     r1, r2);
#endif
	      if (r1 !=  r2) fail = 1;
	      if (fail) break;
	    }
	  /*}*/
    }
  else fail = 1;
    if(fail)
    {
      printf("Sorry, I've not been able to find the tricky number :-(\n");
      printf("\n You should recompile this program undefining TRICKY_ROUND\n");
      printf("in file round.h\n");
      exit(1);
    }

# if defined(Y_USE_SSE2)
    /* now lets try a small check for sse2 constant */
    {
      double Y_ALIGNED(16) AA= 3.0,x0 = 1.4, x1 = 0.6, *rr, r0[4], *r1, y0[2] = {1.0, 1.0};
      Y__M128D MM_A, MM_X0, MM_X1;
      int ii=0;
      /*
	memaling() does not work for intel CC !?, and because of align problem
	in Y_MM_STORE_PD macro, we need to do the follwing awful trick
      */
      /*rr = ALLOC_DOUBLES( 4 );
	r0 = (double *) ALIGN_DOUBLES( rr );*/
      

      r1 = r0 + 2;
      do 
	{
	  AA *= 2.0;
	  Y_MM_SET1_PD( MM_A, AA);
	  Y_MM_SET1_PD( MM_X0, x0);
	  Y_MM_SET1_PD( MM_X1, x1);

	  /* here we prove the most used macro for SSE2 round */
#  if defined(Y_SIMULATE_SSE2)
	  Y_MM_ADD_PD( MM_X0, MM_X0, MM_A);              
	  Y_MM_ADD_PD( MM_X1, MM_X1, MM_A);             
	  Y_MM_SUB_PD( MM_X0, MM_X0, MM_A);              
	  Y_MM_SUB_PD( MM_X1, MM_X1, MM_A);
#  else
	  __asm__ volatile ("addpd %4, %0 \n"          
			    "        addpd %4, %1 \n"           
			    "        subpd %4, %0 \n"           
			    "        subpd %4, %1"              
			    : "=&x" (MM_X0) , "=&x" (MM_X1) :     
			    "0" (MM_X0), "1" (MM_X1), "x"(MM_A));

#  endif /* Y_SIMULATE_SSE2 */	  
	  Y_MM_STORE_PD( &r0[0], MM_X0);
	  Y_MM_STORE_PD( &r0[2], MM_X1);
	  ii++;
	} while ((r0[0] != r1[0]) && (ii < l));
      if (r0[0] != y0[0] || r0[1] != y0[1] || r1[0] != y0[0] || 
	  r1[1] != y0[1])
	{
	  printf("Sorry, SSE2-vector round trick constant search fails\n");
	  exit(1);
	}
      /*free (rr);*/
      /* set SSE2 vector round constant when success */
      Y_MM_MOV_PD( MM_bigA, MM_A);
      Y_MM_MOV_PD( MM_bigB, MM_A);
    }
# endif /* Y_USE_SSE2 */
  return (A);
}

#endif /* TRICKY_ROUND */
/*$Id$*/

