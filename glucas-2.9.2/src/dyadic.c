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
/*
   This the subroutines are related with dyadic mul and square 
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1

void y_squar(y_ptr w1)
{
  y_size_t j;
  y_limb_t t0r,t0i,tm0r,tm0i,wkr,wki,a0r,a0i,a1r,a1i;
  y_size_t mj,nr,jj,ii,l=Y_PLAN[Y_NRADICES-1],
                         inc=Y_POWERS[Y_PLAN[Y_NRADICES-1]-1];

  /* dyadic mul for k=0 */
  t0r=w1[0];
  t0i=w1[1];
  wkr=1.0;
  wki=0.0;
  square_nested_eq( t0 , 1.0 );
  w1[0]=t0r;
  w1[1]=t0i;
  /* now, for every radix reduction make a dyadic mul cycle */
  for (nr=Y_NRADICES-1;nr<Y_NRADICES;nr--)
    {
      if(nr==(Y_NRADICES-1))
        {
          j=1;
          jj=1;
        }
      else
        {
          j=Y_LRIGHT[nr+1];
          jj=0;
        }
      mj=Y_LRIGHT[nr]-1;
      ii=j/Y_PLAN[Y_NRADICES-1];
      cplx_mem_to_local(a1,Y_TWDF[Y_NRADICES-2],ii*inc);
      for(;j<=mj;j++,mj--)
        {
          cplx_mem_to_local(a0,Y_DYADIC,jj);
          cplx_mul(wk,a0,a1);
          cplx_data_to_local(t0,w1,j);
          cplx_data_to_local(tm0,w1,mj);
          square_nested(t0,tm0,wkr,wki);
          cplx_local_to_data(w1,j,t0);
          cplx_local_to_data(w1,mj,tm0);
          jj++;
          if(jj>=l)
            {
              jj-=l;
              ii++;
              cplx_mem_to_local(a1,Y_TWDF[Y_NRADICES-2],ii*inc);
            }
        }
#ifndef NDEBUG1
      printf("Made dyadic mul for pass %i \n",nr);
#endif

    }
}

void y_conv(y_ptr w1, y_ptr w2)
{
  y_limb_t t0r,t0i,tm0r,tm0i,t1r,t1i,tm1r,tm1i,wkr,wki,a0r,a0i,a1r,a1i;
  y_size_t j,mj,nr,jj,ii,l=Y_PLAN[Y_NRADICES-1],
                           inc=Y_POWERS[Y_PLAN[Y_NRADICES-1]-1];

  /* dyadic mul for k=0 */
  t0r=w1[0];
  t0i=w1[1];
  t1r=w2[0];
  t1i=w2[1];
  /*wkr=1.0; wki=0.0;*/
  conv_nested( t0 , t0 , t1 , t1 , 1.0, 0.0 );
  w1[0]=t0r;
  w1[1]=t0i;
  /* now, for every radix reduction make a dyadic mul cycle */
  for (nr=Y_NRADICES-1;nr<Y_NRADICES;nr--)
    {
      if(nr==(Y_NRADICES-1))
        {
          j=1;
          jj=1;
        }
      else
        {
          j=Y_LRIGHT[nr+1];
          jj=0;
        }
      mj=Y_LRIGHT[nr]-1;
      ii=j/Y_PLAN[Y_NRADICES-1];
      cplx_mem_to_local(a1,Y_TWDF[Y_NRADICES-2],ii*inc);
      for(;j<=mj;j++,mj--)
        {
          cplx_mem_to_local(a0,Y_DYADIC,jj);
          cplx_mul(wk,a0,a1);
          cplx_data_to_local(t0,w1,j);
          cplx_data_to_local(tm0,w1,mj);
          cplx_data_to_local(t1,w2,j);
          cplx_data_to_local(tm1,w2,mj);
          conv_nested(t0,tm0,t1,tm1,wkr,wki);
          cplx_local_to_data(w1,j,t0);
          cplx_local_to_data(w1,mj,tm0);
          jj++;
          if(jj>=l)
            {
              jj-=l;
              ii++;
              cplx_mem_to_local(a1,Y_TWDF[Y_NRADICES-2],ii*inc);
            }
        }
#ifndef NDEBUG1
      printf("Made dyadic mul for pass %i \n",nr);
#endif

    }
}
/*$Id$*/









