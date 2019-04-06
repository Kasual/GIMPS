/*$Id$*/
/*
   (c) 2000-2006 Guillermo Ballester Valor 
   
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
   
   It includes ideas from many Authors.
 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "yeafft.h"
#include "ydebug.h"
#include "glucas.h"

#define Y_BITSINT (sizeof(int)*8-1)

/*
   This routine substrac two from an array in a random shift bit
   contest, i. e. substract the bit (Y_SBIT + 1). 
 
   This is different from routines substract_two_n() because here
   all FFTs passes have been made.
 
   NOTE: IT SUPOSSES THIS ROUTINE IS CALLED BEFORE NORMALIZATION,
   SO WE HAVE TO TRANSFORM WITH 1/TWO_TO_MINUS_PHI 
*/
void substract_two(BIG_DOUBLE *x, UL N)
{
  BIG_DOUBLE xtwo;
  UL ib, ni, a;
  /* the bit to substract */
  ib = Y_SBIT + 1;
  if(ib == Y_EXPONENT)
    ib = 0;
  /* The element ni where bit 'ib' is */
  ni=floor((BIG_DOUBLE)ib/Y_XBITS);
  /* now ib is the bit within x[ni] */
  ib -= (UL)ceil(Y_XBITS*ni);
  /* The shifted value of two */
  xtwo = (BIG_DOUBLE)(1<<ib);

  /*computes (N - ni*b % N). We suposse up to 2^30 exponents */
  if(ni)
    {
#if ULONG_MAX > 0xFFFFFFFF
      /* 64 bits integer machine */
      a = N - (ni * b) % N;
#else
      /* 32 bits integer machine */
      modmul_32 (a, ni, b, N);
      a = N - a;
#endif

    }
  else
    a=0;
  /* Multiply by 1/two_to_minusphi() */
  xtwo *= (0.5 * exp( a * M_LN2 / N) * N);
  /* finally substract two */
  x[addr(ni)] -= xtwo;
}

void normalize(BIG_DOUBLE *x, UL N, UL err_flag)
{
  BIG_DOUBLE hiinv=highinv, loinv=lowinv;
  BIG_DOUBLE temp0,tempErr;
  BIG_DOUBLE maxerr=0.0,err=0.0,ttmpSmall=Gsmall,ttmpBig=Gbig,
                                          ttmp,ttp,ttmpTemp;
  BIG_DOUBLE carry,ttpSmall=Hsmall;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  UL i,j,k,lastloop,bj,offset,newOffset;


  ASSERT_ALIGNED_DOUBLE();

  if(Y_SBIT==0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two( x, N);
      carry = 0.0;
    }
  lastloop = N - UPDATE;

  bj = N;
  offset=0;

  for (i=0,j=0; j<N; j+=UPDATE,i++)
    {
      ttmp = two_to_minusphi[addr(i)];
      ttp = two_to_phi[addr(i)];

      for (k=1; k<UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr((j + k));
          tempErr = RINT( temp0*ttmp );
          err = fabs(temp0*ttmp-tempErr);
          if (err>maxerr)
            maxerr=err;
          temp0 = tempErr + carry;
          if (bj >= c)
            {
              ttp +=ttp;
              temp0 *= hiinv;
              carry = RINT(temp0);
              bj-=c;
              ttmp *=ttmpBig;
            }
          else
            {
              temp0 *= loinv;
              carry = RINT(temp0);
              bj+=b;
              ttmp *=ttmpSmall;
            }
          x[offset] = (temp0-carry) * ttp;
          ttp *= ttpSmall;
          offset=newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j == lastloop)
        bj = 0;
      tempErr = RINT(temp0*ttmp);
      err = fabs(temp0*ttmp-tempErr);
      if (err>maxerr)
        maxerr=err;
      temp0 = tempErr + carry;
      if (bj >= c)
        {
          temp0 *= hiinv;
          ttp +=ttp;
          carry = RINT(temp0);
          bj-=c;
          ttmp *=ttmpBig;
        }
      else
        {
          temp0 *= loinv;
          carry = RINT(temp0);
          bj+=b;
          ttmp *=ttmpSmall;
        }
      x[offset] = (temp0-carry) * ttp;
      ttp *= ttpSmall;
      offset=newOffset;
    }
  bj = N;
  ttmp = two_to_minusphi[0];
  ttp = two_to_phi[0];

  k=0;
  while(carry != 0)
    {
      ttmpTemp = x[addr(k)]*ttmp*N*2;
      temp0 = (ttmpTemp + carry);
      if (bj >= c)
        {
          bj-=c;
          ttp += ttp;
          temp0 *= hiinv;
          ttmp *= ttmpBig;
          carry = RINT(temp0);
        }
      else
        {
          bj +=b;
          temp0 *= loinv;
          ttmp *= ttmpSmall;
          carry = RINT(temp0);
        }
      x[addr(k)] = (temp0 -carry)* ttp;
      ttp *= ttpSmall;
      k++;
    }
  if(err_flag)
    Err= maxerr;
}

void first_normalize(BIG_DOUBLE *x, UL N)
{
  BIG_DOUBLE hiinv = highinv, loinv = lowinv;
  BIG_DOUBLE temp0;
  BIG_DOUBLE ttp, ttmpTemp;
  BIG_DOUBLE carry, ttpSmall = Hsmall;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A = bigA, B = bigB;
#endif

  UL i, j, k, lastloop, bj, offset, newOffset;


  /*ASSERT_ALIGNED_DOUBLE();*/

  carry = 0.0;

  lastloop = N - UPDATE;

  bj = N;
  offset = 0;

  for (i = 0,j = 0; j < N; j += UPDATE, i++)
    {
      ttp = two_to_phi[addr(i)];

      for (k = 1; k < UPDATE; k++)
        {
          temp0 = x[offset] + carry;
          newOffset = addr((j + k));
          if (bj >= c)
            {
              ttp += ttp;
              temp0 *= hiinv;
              carry = RINT(temp0);
              bj -= c;
            }
          else
            {
              temp0 *= loinv;
              carry = RINT(temp0);
              bj += b;
            }
          x[offset] = (temp0-carry) * ttp;
          ttp *= ttpSmall;
          offset = newOffset;
        }
      temp0 = x[offset] + carry;
      newOffset = addr((j + UPDATE));
      if (j == lastloop)
        bj = 0;
      if (bj >= c)
        {
          temp0 *= hiinv;
          ttp += ttp;
          carry = RINT(temp0);
          bj -= c;
        }
      else
        {
          temp0 *= loinv;
          carry = RINT(temp0);
          bj += b;
        }
      x[offset] = (temp0 - carry) * ttp;
      ttp *= ttpSmall;
      offset = newOffset;
    }

  bj = N;

  ttp = two_to_phi[0];

  k=0;
  while(carry != 0)
    {
      ttmpTemp = x[addr(k)] * low / ttp;
      temp0 = (ttmpTemp + carry);
      if (bj >= c)
        {
          bj -= c;
          ttp += ttp;
          temp0 *= hiinv;
          carry = RINT(temp0);
        }
      else
        {
          bj += b;
          temp0 *= loinv;
          carry = RINT(temp0);
        }
      x[addr(k)] = (temp0 - carry) * ttp;
      ttp *= ttpSmall;
      k++;
    }
}



void last_normalize(BIG_DOUBLE *x, UL N, UL err_flag)
{
  BIG_DOUBLE hi = high, hiinv = highinv, lo = low, loinv = lowinv;
  BIG_DOUBLE temp0, tempErr;
  BIG_DOUBLE maxerr = 0.0, err = 0.0, ttmpSmall = Gsmall, ttmpBig = Gbig,
                                 ttmp;
  BIG_DOUBLE carry;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A = bigA,B = bigB;
#endif

  UL i, j, k, lastloop, bj, offset, newOffset;


  ASSERT_ALIGNED_DOUBLE();

  if(Y_SBIT == 0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two (x, N);
      carry = 0.0;
    }

  lastloop = N - UPDATE;

  bj = N;
  offset = 0;

  for (j = 0,i = 0; j < N; j += UPDATE, i++)
    {
      ttmp = two_to_minusphi[addr(i)];

      for (k=1; k < UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr((j + k));
          tempErr = RINT(temp0 * ttmp);
          err = fabs(temp0 * ttmp - tempErr);
          if (err > maxerr)
            maxerr = err;
          temp0 = tempErr + carry;
          if ( bj >= c)
            {
              temp0 *= hiinv;
              carry = RINT(temp0);
              bj -= c;
              ttmp *= ttmpBig;
              if(carry > temp0)
                carry -= 1.0;
              x[offset] = (temp0 - carry) * hi;
            }
          else
            {
              temp0 *= loinv;
              carry = RINT(temp0);
              if(carry > temp0)
                carry -= 1.0;
              bj += b;
              ttmp *= ttmpSmall;
              x[offset] = (temp0 - carry) * lo;
            }
          offset = newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j == lastloop)
        bj = 0;
      tempErr = RINT(temp0 * ttmp);
      err = fabs(temp0 * ttmp - tempErr);
      if (err > maxerr)
        maxerr = err;
      temp0 = tempErr + carry;
      if (bj >= c)
        {
          temp0 *= hiinv;
          carry = RINT(temp0);
          if(carry > temp0)
            carry-=1.0;
          bj-=c;
          ttmp *=ttmpBig;
          x[offset] = (temp0-carry) * hi;
        }
      else
        {
          temp0 *= loinv;
          carry = RINT(temp0);
          if(carry > temp0)
            carry -= 1.0;
          bj += b;
          ttmp *= ttmpSmall;
          x[offset] = (temp0 - carry) * lo;
        }
      offset = newOffset;
    }
  bj = N;
  k = 0;
  while(carry != 0 && k < (N - 1))
    {
      temp0 = (x[addr(k)] + carry);
      if (bj >= c)
        {
          bj -= c;
          temp0 *= hiinv;
          carry = RINT(temp0);
          if(carry > temp0)
            carry -= 1.0;
          x[addr(k)] = (temp0 - carry) * hi;
        }
      else
        {
          bj += b;
          temp0 *= loinv;
          carry = RINT(temp0);
          if(carry > temp0)
            carry -= 1.0;
          x[addr(k)] = (temp0 - carry) * lo;
        }
      k++;
    }
  if (carry)
    {
      temp0 = (x[addr(k)] + carry);
      temp0 *= loinv;
      carry = RINT(temp0);
      if(carry > temp0)
        carry -= 1.0;
      x[addr(k)] = (temp0 - carry) * lo;
    }
  if(err_flag)
    Err = maxerr;
  /*balancetosdrep(x,N);*/
}

void y_normalize(BIG_DOUBLE *x,UL N, UL err_flag )
{
  BIG_DOUBLE hi=high, hiinv=highinv, lo=low, loinv=lowinv;
  BIG_DOUBLE temp0, tempErr;
  BIG_DOUBLE maxerr=0.0, err=0.0, ttmpSmall = Gsmall, ttmpBig = Gbig,
                             ttmp, ttp, ttmpTemp;
  BIG_DOUBLE carry, ttpSmall = Hsmall, ttpBig = Hbig;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A = bigA, B = bigB;
#endif

  UL i, j, k, lastloop, offset, newOffset, bj;

  ASSERT_ALIGNED_DOUBLE();

  if(Y_SBIT==0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two( x, N);
      carry = 0.0;
    }

  lastloop = N-UPDATE;

  bj=N;
  bj=0;
  offset=0;

  for (i=0,j=0; j<N; j+=UPDATE,i++)
    {
      ttmp = two_to_minusphi[addr(i)];
      ttp = two_to_phi[addr(i)];

      for (k=1; k<UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr(j + k);
          tempErr = RINT( temp0*ttmp );
          if (err_flag)
            {
              err = fabs(temp0*ttmp-tempErr);
              if (err>maxerr)
                maxerr=err;
            }
          temp0 = tempErr + carry;
          if (bj >= c)
            {
              bj-=c;
              temp0 *= hiinv;
              ttp += ttp;
            }
          else
            {
              bj+=b;
              temp0 *= loinv;
              ttmp += ttmp;
            }
          carry = RINT(temp0);
          ttmp *= ttmpBig;
          x[offset] = (temp0-carry) * ttp ;
          ttp *= ttpSmall;
          offset=newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j==lastloop)
        bj = 0;
      tempErr = RINT(temp0*ttmp);
      if (err_flag)
        {
          err = fabs(temp0*ttmp-tempErr);
          if (err>maxerr)
            maxerr=err;
        }
      temp0 = tempErr + carry;
      if (bj>=c)
        {
          bj-=c;
          temp0 *= hiinv;
          ttp += ttp;
        }
      else
        {
          bj+=b;
          temp0 *= loinv;
          ttmp += ttmp;
        }
      carry = RINT(temp0);
      ttmp *= ttmpBig;
      x[offset] = (temp0-carry) * ttp;
      ttp *= ttpSmall;
      offset=newOffset;
    }
  bj = N;
  ttmp = two_to_minusphi[0];
  ttp = two_to_phi[0];

  k=0;
  while(carry != 0.0)
    {
      ttmpTemp = x[addr(k)]*ttmp*N*2;
      temp0 = (ttmpTemp + carry);
      if (bj >= c)
        {
          bj  -= c;
          temp0 *= hiinv;
          ttmp *= ttmpBig;
          carry = RINT(temp0);
          x[addr(k)] = (temp0 - carry) * ttp * hi;
          ttp *= ttpBig;
        }
      else
        {
          bj += b;
          temp0 *= loinv;
          ttmp *= ttmpSmall;
          carry = RINT(temp0);
          x[addr(k)] = (temp0 -carry) * ttp * lo;
          ttp *= ttpSmall;
        }
      k++;
    }
  if(err_flag)
    Err=maxerr;
}


void y2_normalize(BIG_DOUBLE *x,UL N, UL err_flag )
{
  UL i, j, k, lastloop, offset, newOffset, bj;
  BIG_DOUBLE temp0,tempErr;
  double Y_ALIGNED(16) w[2], winv[2], wttp[2], wttmp[2];
  BIG_DOUBLE maxerr=0.0,err=0.0,ttmp,ttp,ttmpTemp,carry;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A = bigA, B = bigB;
#endif

  int isbig;

  w[0]=low;
  w[1]=high;
  winv[0]=lowinv;
  winv[1]=highinv;
  wttp[0]=Hsmall;
  wttp[1]=Hbig;
  wttmp[0]=Gsmall;
  wttmp[1]=Gbig;

  if(Y_SBIT==0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two( x, N);
      carry = 0.0;
    }

  ASSERT_ALIGNED_DOUBLE();
  lastloop = N-UPDATE;

  bj=N;
  isbig = 1;
  offset = 0;

  for (i=0,j=0; j<N; j+=UPDATE,i++)
    {
      ttmp = two_to_minusphi[addr(i)];
      ttp = two_to_phi[addr(i)];

      for (k=1; k<UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr(j + k);
          tempErr = RINT( temp0*ttmp );
          if (err_flag)
            {
              err = fabs(temp0*ttmp-tempErr);
              if (err>maxerr)
                maxerr=err;
            }
          temp0 = tempErr + carry;
          bj -= c;
          temp0 *= winv[isbig];
          bj += (bj>>Y_BITSINT) & N;
          carry = RINT(temp0);
          ttmp *= wttmp[isbig];
          x[offset] = (temp0-carry) * (ttp * w[isbig]);
          ttp *= wttp[isbig];
          isbig = 1+((bj-(int)c)>>Y_BITSINT);
          offset=newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j==lastloop)
        isbig = 0;
      tempErr = RINT(temp0*ttmp);
      if (err_flag)
        {
          err = fabs(temp0*ttmp-tempErr);
          if (err>maxerr)
            maxerr=err;
        }
      temp0 = tempErr + carry;
      bj -= c;
      temp0 *= winv[isbig];
      bj += (bj>>Y_BITSINT) & N;
      carry = RINT(temp0);
      ttmp *= wttmp[isbig];
      x[offset] = (temp0-carry) * (ttp * w[isbig]);
      ttp *= wttp[isbig];
      isbig = 1+((bj-(int)c)>>Y_BITSINT);
      offset=newOffset;
    }
  bj = N;
  ttmp = two_to_minusphi[0];
  ttp = two_to_phi[0];
  isbig=1;
  k=0;
  while(carry != 0)
    {
      ttmpTemp = x[addr(k)]*ttmp*N*2;
      bj += b;
      temp0 = (ttmpTemp + carry);
      if(bj >= N)
        bj-=N;
      temp0 *= winv[isbig];
      ttmp *= wttmp[isbig];
      carry = RINT(temp0);
      x[addr(k)] = (temp0 - carry)* ttp * w[isbig];
      ttp *= wttp[isbig];
      isbig = (bj>=c);
      k++;
    }
  if(err_flag)
    Err=maxerr;
}


void y_last_normalize(BIG_DOUBLE *x,UL N,UL err_flag )
{
  BIG_DOUBLE hi = high, hiinv = highinv, lo = low, loinv = lowinv;
  BIG_DOUBLE temp0, tempErr;
  BIG_DOUBLE maxerr = 0.0,err = 0.0, ttmpSmall = Gsmall, ttmpBig = Gbig,
                                ttmp;
  BIG_DOUBLE carry;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  UL i, j, k, lastloop, bj, offset, newOffset;
  UL size0;


  if(Y_SBIT==0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two( x, N);
      carry = 0.0;
    }

  lastloop = N-UPDATE;

  bj=N;
  size0 = 1;
  offset=0;

  for (j=0,i=0; j<N; j+=UPDATE,i++)
    {
      ttmp = two_to_minusphi[addr(i)];

      for (k=1; k<UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr((j + k));
          tempErr = RINT(temp0*ttmp);
          if (err_flag)
            {
              err = fabs(temp0*ttmp-tempErr);
              if (err>maxerr)
                maxerr=err;
            }
          temp0 = tempErr + carry;
          if (bj>=c)
            {
              temp0 *= hiinv;
              carry = RINT(temp0 - 0.5);
              bj-=c;
              ttmp *=ttmpBig;
              x[offset] = (temp0-carry) * hi;
            }
          else
            {
              temp0 *= loinv;
              carry = RINT(temp0 - 0.5);
              bj+=b;
              ttmp *=ttmpSmall;
              x[offset] = (temp0-carry) * lo;
            }
          offset=newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j==lastloop)
        bj=0;
      tempErr = RINT(temp0*ttmp);
      if (err_flag)
        {
          err = fabs(temp0*ttmp-tempErr);
          if (err>maxerr)
            maxerr=err;
        }
      temp0 = tempErr + carry;
      if (bj>=c)
        {
          temp0 *= hiinv;
          carry = RINT(temp0 - 0.5);
          bj-=c;
          ttmp *=ttmpBig;
          x[offset] = (temp0-carry) * hi;
        }
      else
        {
          temp0 *= loinv;
          carry = RINT(temp0 - 0.5);
          bj+=b;
          ttmp *=ttmpSmall;
          x[offset] = (temp0-carry) * lo;
        }
      offset=newOffset;
    }
  bj = N;
  k=0;
  while(carry != 0)
    {
      size0 = (bj>=c);
      bj += b;
      temp0 = (x[addr(k)] + carry);
      if(bj >= N)
        bj-=N;
      if (size0)
        {
          temp0 *= hiinv;
          carry = RINT(temp0 - 0.5);
          x[addr(k)] = (temp0 -carry) * hi;
        }
      else
        {
          temp0 *= loinv;
          carry = RINT(temp0 - 0.5);
          x[addr(k)] = (temp0 -carry)* lo;
        }
      k++;
    }
  if(err_flag)
    Err=maxerr;
}

void y3_normalize(BIG_DOUBLE *x,UL N, UL err_flag )
{
  UL i, j, k, bj, lastloop, offset, newOffset;
  BIG_DOUBLE hi=high, hiinv=highinv, lo=low, loinv=lowinv;
  BIG_DOUBLE temp0,tempErr;
  BIG_DOUBLE maxerr=0.0, err=0.0, ttmpSmall=Gsmall, ttmpBig=Gbig, ttmp, ttp, ttmpTemp;
  double Y_ALIGNED(16) tt[2];
  BIG_DOUBLE carry,ttpSmall=Hsmall,ttpBig=Hbig;
#if defined(TRICKY_ROUND)

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  UL size0;
  int issmall;

  if(Y_SBIT==0)
    carry = - 2.0; /* this is the -2 of the LL x*x - 2 */
  else
    {
      substract_two( x, N);
      carry = 0.0;
    }

  lastloop = N-UPDATE;

  bj=N;
  size0 = 1;
  offset=0;

  for (i=0,j=0; j<N; j+=UPDATE,i++)
    {
      tt[1] = two_to_minusphi[addr(i)];
      tt[0] = two_to_phi[addr(i)] * lo;

      for (k=1; k<UPDATE; k++)
        {
          temp0 = x[offset];
          newOffset = addr(j + k);
          tempErr = RINT( temp0*tt[1] );
          if (err_flag)
            {
              err = fabs(temp0*tt[1]-tempErr);
              if (err>maxerr)
                maxerr=err;
            }
          temp0 = tempErr + carry;
          bj -= c;
          issmall=(bj>>Y_BITSINT);
          temp0 *= hiinv;
          bj += issmall & N;
          issmall &= 1;
          temp0 += issmall*temp0;
          tt[issmall] += tt[issmall];
          carry = RINT(temp0);
          tt[1] *= ttmpBig;
          x[offset] = (temp0-carry) * tt[0] ;
          tt[0] *= ttpSmall;
          offset=newOffset;
        }
      temp0 = x[offset];
      newOffset = addr((j + UPDATE));
      if (j==lastloop)
        bj = 0;
      tempErr = RINT(temp0*tt[1]);
      if (err_flag)
        {
          err = fabs(temp0*tt[1]-tempErr);
          if (err>maxerr)
            maxerr=err;
        }
      temp0 = tempErr + carry;
      if (bj>=c)
        {
          temp0 *= hiinv;
          bj-=c;
          tt[0] += tt[0];
        }
      else
        {
          temp0 *= loinv;
          bj+=b;
          tt[1] += tt[1];
        }
      carry = RINT(temp0);
      tt[1] *=ttmpBig;
      x[offset] = (temp0-carry) * tt[0];
      tt[0] *= ttpSmall;
      offset=newOffset;
    }
  bj = N;
  ttmp = two_to_minusphi[0];
  ttp = two_to_phi[0];

  k=0;
  while(carry != 0)
    {
      ttmpTemp = x[addr(k)]*ttmp*N*2;
      size0 = (bj>=c);
      bj += b;
      temp0 = (ttmpTemp + carry);
      if(bj >= N)
        bj-=N;
      if (size0)
        {
          temp0 *= hiinv;
          ttmp *= ttmpBig;
          carry = RINT(temp0);
          x[addr(k)] = (temp0 - carry) * ttp * hi;
          ttp *= ttpBig;
        }
      else
        {
          temp0 *= loinv;
          ttmp *= ttmpSmall;
          carry = RINT(temp0);
          x[addr(k)] = (temp0 -carry) * ttp * lo;
          ttp *= ttpSmall;
        }
      k++;
    }
  if(err_flag)
    Err=maxerr;
}
/*$Id$*/






