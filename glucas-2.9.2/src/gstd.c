/*$Id$*/
/*
    gstd.c. An collection of routines to transform from DWT-float residues
    to an all-integer format (an viceversa)
 
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
   This is a file including the routines to trnasform residues from
   float to integer
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
/* Patch suggested by B.J.Beesley */
# if defined SUN_V9_GCC
#  undef ULONG_MAX
#  define ULONG_MAX 0xFFFFFFFF
# endif
#include "yeafft.h"
#include "glucas.h"

#define BITS_PER_UL (sizeof(UL)*8)

/* pass from real in *x to standar int residue in *tar */
/* It returns the number of UL written in tar          */
int write_standar(UL q, UL n, BIG_DOUBLE *x, UL *tar, UL *sumcheck)
{
  BIG_DOUBLE a,carry=0.0;
#ifdef TRICKY_ROUND

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  UL  bj=n,aux;
  int nb,ns,j=0,k=0,jbit=0,ibit;

  ns= (int)(q / n);
  nb= ns + 1;
  tar[0]= (UL)0;
  do
    {
      a = x[addr(j)] + carry;
      ibit = jbit;
      if (bj >= c)
        {
          bj -= c;
          if (a < 0.0 )
            {
              a += high;
              carry = -1.0;
            }
          else
            carry=0.0;
          aux =(UL) RINT(a);
          jbit += nb;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          tar[k] += (aux << ibit);
          if(jbit < ibit)
            {
              add_sumcheck( sumcheck , tar[k]);
              tar[++k] = (aux >> (nb-jbit));
            }
        }
      else
        {
          bj += b;
          if (a < 0.0)
            {
              a += low;
              carry = -1.0;
            }
          else
            carry=0.0;

          aux =(UL) RINT(a);
          jbit += ns;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          tar[k] += (aux << ibit);
          if(jbit < ibit)
            {
              add_sumcheck( sumcheck , tar[k]);
              tar[++k] = (aux >> (ns-jbit));
            }
        }
      j++;
      if(j == ((int)n-1) )
        bj=0;
    }
  while ((UL)j<n);
  if(jbit != 0)
    {
      add_sumcheck( sumcheck , tar[k]);
      k++;
    }
  if (carry == 0.0)
    return k;
  j=0;
  add_sumcheck (sumcheck, ~1U );
  do
    {
      aux=tar[j];
      tar[j]--;
      j++;
    }
  while (aux == 0);
  return k;
}

/*
   This is almost a copy of write_standar() but it does not write an
   integer residue. It only computes the integer residue 
   mod (2 BITS_PER_UL -1). 
   It assumes shift = 0
*/
int residue_sumcheck(UL q, UL n, BIG_DOUBLE *x, UL *sumcheck)
{
  BIG_DOUBLE a,carry=0.0;
#ifdef TRICKY_ROUND

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  UL  bj=n,aux,auxb;
  int nb,ns,j=0,k=0,jbit=0,ibit;

  ns= (int)(q / n);
  nb= ns + 1;
  auxb = (UL)0;
  do
    {
      a = x[addr(j)] + carry;
      ibit = jbit;
      if (bj >= c)
        {
          bj -= c;
          if (a < 0.0 )
            {
              a += high;
              carry = -1.0;
            }
          else
            carry=0.0;
          aux =(UL) RINT(a);
          jbit += nb;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          auxb += (aux << ibit);
          if(jbit < ibit)
            {
              add_sumcheck( sumcheck , auxb);
              auxb = (aux >> (nb-jbit));
              k++;
            }
        }
      else
        {
          bj += b;
          if (a < 0.0)
            {
              a += low;
              carry = -1.0;
            }
          else
            carry=0.0;

          aux =(UL) RINT(a);
          jbit += ns;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          auxb += (aux << ibit);
          if(jbit < ibit)
            {
              add_sumcheck( sumcheck , auxb);
              auxb = (aux >> (ns-jbit));
              k++;
            }
        }
      j++;
      if(j == ((int)n-1) )
        bj=0;
    }
  while ((UL)j<n);
  if(jbit != 0)
    {
      add_sumcheck( sumcheck , auxb);
      k++;
    }
  if (carry == 0.0)
    return auxb;
  j=0;
  add_sumcheck (sumcheck, ~1U );
  return k;
}


/*
   It computes the residue mod 2^64-1
   If BITS_PER_UL == 32 , *sumcheck has the low half (less significant)
   and the returned value are the high half (most significant bits)
 
   In a 64 bits, both *sumcheck and the returned value are the same 
*/
UL residue64_sumcheck_shift(UL q, UL n, BIG_DOUBLE *x, UL *sumcheck)
{
  UL i,j,k=Y_SBIT,aux,aux2,sum=0;
  unsigned char b[8];

  *sumcheck = 0;
  i= ((q-1) / 64U) + 1U;

  for (j=0;j<(i-1);j++,k+=64)
    {
      residue_char(x, q, n, b, 64, k);
#if ULONG_MAX > 0xFFFFFFFF

      aux = (UL)b[7];
      aux += ((UL)b[6])<<8;
      aux += ((UL)b[5])<<16;
      aux += ((UL)b[4])<<24;
      aux += ((UL)b[3])<<32;
      aux += ((UL)b[2])<<40;
      aux += ((UL)b[1])<<48;
      aux += ((UL)b[0])<<56;
      add_sumcheck( sumcheck, aux);
#else

      aux = (UL)b[7];
      aux += ((UL)b[6])<<8;
      aux += ((UL)b[5])<<16;
      aux += ((UL)b[4])<<24;
      aux2 = (UL)b[3];
      aux2 += ((UL)b[2])<<8;
      aux2 += ((UL)b[1])<<16;
      aux2 += ((UL)b[0])<<24;
      add_sumcheck2( &sum, sumcheck, aux2, aux);
#endif

    }
  for (j=0;j<8;j++)
    b[j]=0;
  i= q % 64U;
  if(i)
    {
      residue_char(x, q, n, b, i, k);
#if ULONG_MAX > 0xFFFFFFFF

      aux = (UL)b[7];
      aux += ((UL)b[6])<<8;
      aux += ((UL)b[5])<<16;
      aux += ((UL)b[4])<<24;
      aux += ((UL)b[3])<<32;
      aux += ((UL)b[2])<<40;
      aux += ((UL)b[1])<<48;
      aux += ((UL)b[0])<<56;
      add_sumcheck( sumcheck, aux);
#else

      aux = (UL)b[7];
      aux += ((UL)b[6])<<8;
      aux += ((UL)b[5])<<16;
      aux += ((UL)b[4])<<24;
      aux2 = (UL)b[3];
      aux2 += ((UL)b[2])<<8;
      aux2 += ((UL)b[1])<<16;
      aux2 += ((UL)b[0])<<24;
      add_sumcheck2( &sum, sumcheck, aux2, aux);
#endif

    }
#if ULONG_MAX > 0xFFFFFFFF
  return *sumcheck;
#else

  return sum;
#endif
}


/*        pass from standar int residue in *tar to real in *x          */
/*        It returns the number of x's readed                          */
int read_standar(UL q, UL n, BIG_DOUBLE *x, UL *tar,UL *sumcheck)
{
  BIG_DOUBLE hi,lo;
  UL  a,bj=n,aux,masklow,maskbig,bl,sl,bb,cc;
  int j=0,k=0,jbit=0,ibit,nb,ns,icarry=0;

  ns= (int)( q/n);
  nb= ns + 1;
  masklow = ((UL)1 << ns) - (UL)1;
  maskbig = ((UL)1 << nb) - (UL)1;
  bl= (UL)1<<(nb-1);
  sl= (UL)1<<(ns-1);
  bb= q % n;
  cc= n - bb;
  hi= (BIG_DOUBLE)((UL)1 << nb);
  lo= (BIG_DOUBLE)((UL)1 << ns);

  a=tar[0];
  add_sumcheck(sumcheck, a);
  /* *sumcheck += a; */
  do
    {
      ibit=jbit;
      if (bj >= cc)
        {
          bj -= cc;
          jbit += nb;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          aux = (a >> ibit) & maskbig;
          if(jbit < ibit)
            {
              a=tar[++k];
              add_sumcheck(sumcheck, a);
              /* *sumcheck += a; */
              aux |= ((maskbig >> (nb - jbit)) & a)<< (nb - jbit);
            }
          aux += (UL)icarry;
          icarry = (aux >= bl);
          if (icarry)
            x[addr(j)]= (BIG_DOUBLE)aux - hi;
          else
            x[addr(j)]= (BIG_DOUBLE) aux;
          j++;
        }
      else
        {
          bj += bb;
          jbit += ns;
          if( jbit >= (int)BITS_PER_UL)
            jbit -= BITS_PER_UL;
          aux = (a >> ibit) & masklow;
          if(jbit < ibit)
            {
              a=tar[++k];
              add_sumcheck(sumcheck, a);
              /* *sumcheck += a; */
              aux |= ((masklow >> (ns - jbit)) & a)<< (ns - jbit);
            }
          aux += (UL)icarry;
          icarry= (aux >= sl);
          if (icarry)
            x[addr(j)]= (BIG_DOUBLE)aux - lo;
          else
            x[addr(j)]= (BIG_DOUBLE) aux;
          j++;
        }
      if(j == ((int)n-1) )
        bj=0;
    }
  while ( j< (int)n);
  if (icarry)
    x[0]+=(BIG_DOUBLE) icarry;
  return j;
}






