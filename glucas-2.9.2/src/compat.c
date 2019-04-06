/*$Id$*/
/*
    compat.c. An implementation of Interchangeable Mersenne Residue Format.
    This file includes a collection of needed routines
 
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
#include "stdlib.h"
#include "stdio.h"
#include "limits.h"

/* Patch suggested by B.J.Beesley */
# if defined(SUN_V9_GCC) || defined (__sparc)
#  undef ULONG_MAX
#  define ULONG_MAX 0xFFFFFFFF
# endif
#include "yeafft.h"
#include "glucas.h"

/* this routine writes an unsigned long to be read by any other system
	the writting is byte oriented:
 
	byte 0 : bits 0 to 7
	byte 1 : bits 8 to 15
	byte 2 : bits 16 to 23
	byte 3 : bits 24 to 31
	byte 4 : bits 32 to 39         (0 for 32 bits UL)
	byte 5 : bits 40 to 47         (0 for 32 bits UL)
	byte 6 : bits 48 to 55         (0 for 32 bits UL)
	byte 7 : bits 56 to 63         (0 for 32 bits UL)
 
*/
size_t fwrite_UL_pad(UL u, FILE * chkpnt)
{
  UL mask=(UL)255;
  size_t retval;
  unsigned char  aux[8];


#if ULONG_MAX > 0xFFFFFFFF
  /* case 64 bits */
  aux[0]= (unsigned char) (u & mask);
  aux[1]= (unsigned char) ((u>>8) & mask);
  aux[2]= (unsigned char) ((u>>16) & mask);
  aux[3]= (unsigned char) ((u>>24) & mask);
  aux[4]= (unsigned char) ((u>>32) & mask);
  aux[5]= (unsigned char) ((u>>40) & mask);
  aux[6]= (unsigned char) ((u>>48) & mask);
  aux[7]= (unsigned char) ((u>>56) & mask);
#else
  /* case 32 bits */
  aux[0]= (unsigned char) (u & mask);
  aux[1]= (unsigned char) ((u>>8) & mask);
  aux[2]= (unsigned char) ((u>>16) & mask);
  aux[3]= (unsigned char) ((u>>24) & mask);
  aux[4]= 0;
  aux[5]= 0;
  aux[6]= 0;
  aux[7]= 0;
#endif

  retval = fwrite (&aux[0],1,8, chkpnt);
  if (retval == 8 || retval == 1)
    return 1;
  return retval;
}

size_t fwrite_UL(UL u, FILE * chkpnt)
{
  UL mask=(UL)255;
  size_t retval;
#if ULONG_MAX > 0xFFFFFFFF
  /* case 64 bits */
  unsigned char  aux[8];
  aux[0]= (unsigned char) (u & mask);
  aux[1]= (unsigned char) ((u>>8) & mask);
  aux[2]= (unsigned char) ((u>>16) & mask);
  aux[3]= (unsigned char) ((u>>24) & mask);
  aux[4]= (unsigned char) ((u>>32) & mask);
  aux[5]= (unsigned char) ((u>>40) & mask);
  aux[6]= (unsigned char) ((u>>48) & mask);
  aux[7]= (unsigned char) ((u>>56) & mask);
  retval = fwrite (&aux[0],1,8, chkpnt);
  if (retval == 8 || retval == 1)
    return 1;
  return retval;
#else
  /* case 32 bits */
  unsigned char  aux[4];
  aux[0]= (unsigned char) (u & mask);
  aux[1]= (unsigned char) ((u>>8) & mask);
  aux[2]= (unsigned char) ((u>>16) & mask);
  aux[3]= (unsigned char) ((u>>24) & mask);
  retval = fwrite (&aux[0],1,4, chkpnt);
  if (retval == 4 || retval == 1)
    return 1;
  return retval;
#endif
}

/* writes two 32 bits words. Lower is u0 */
size_t fwrite_twin(UL u0, UL u1, FILE * chkpnt)
{
  UL mask=(UL)255;
  size_t retval;
  unsigned char  aux[8];
  aux[0]= (unsigned char) (u0 & mask);
  aux[1]= (unsigned char) ((u0>>8) & mask);
  aux[2]= (unsigned char) ((u0>>16) & mask);
  aux[3]= (unsigned char) ((u0>>24) & mask);
  aux[4]= (unsigned char) (u1 & mask);
  aux[5]= (unsigned char) ((u1>>8) & mask);
  aux[6]= (unsigned char) ((u1>>16) & mask);
  aux[7]= (unsigned char) ((u1>>24) & mask);
  retval = fwrite (&aux[0],1,8, chkpnt);
  if (retval == 8 || retval == 1)
    return 1;
  return retval;
}



size_t fread_UL_pad(UL * u, FILE * chkpnt)
{
  unsigned char  aux[8];
  UL u1,u2,u3;
  size_t retval;

  /* this routine reads a universal UL */
  retval = fread( &aux[0], sizeof(unsigned char), 8, chkpnt);
  if( (retval != 1) && (retval != 8*sizeof(unsigned char)))
    return retval;
  *u = (UL) aux[0];
  u1 = ((UL) aux[1])<<8;
  u2 = ((UL) aux[2])<<16;
  u3 = ((UL) aux[3])<<24;
#if ULONG_MAX > 0xFFFFFFFF

  *u += (UL) aux[4]<<32;
  u1 += ((UL) aux[5])<<40;
  u2 += ((UL) aux[6])<<48;
  u3 += ((UL) aux[7])<<56;
#endif

  *u += (u1 + u2 + u3);
  return 1;
}

size_t fread_UL(UL * u, FILE * chkpnt)
{
  UL u1,u2,u3;
  size_t retval;
#if ULONG_MAX > 0xFFFFFFFF

  unsigned char  aux[8];

  retval = fread( &aux[0], sizeof(unsigned char), 8, chkpnt);
  if( (retval != 1) && (retval != 8*sizeof(unsigned char)))
    return retval;
  *u = (UL) aux[0];
  u1 = ((UL) aux[1])<<8;
  u2 = ((UL) aux[2])<<16;
  u3 = ((UL) aux[3])<<24;
  *u += (UL) aux[4]<<32;
  u1 += ((UL) aux[5])<<40;
  u2 += ((UL) aux[6])<<48;
  u3 += ((UL) aux[7])<<56;
#else

  unsigned char  aux[4];

  retval = fread( &aux[0], sizeof(unsigned char), 4, chkpnt);
  if( (retval != 1) && (retval != 4*sizeof(unsigned char)))
    return retval;
  *u = (UL) aux[0];
  u1 = ((UL) aux[1])<<8;
  u2 = ((UL) aux[2])<<16;
  u3 = ((UL) aux[3])<<24;
#endif

  *u += (u1 + u2 + u3);
  return 1;
}


size_t fread_twin(UL *v0, UL *v1, FILE * chkpnt)
{
  unsigned char  aux[8];
  UL u1,u2,u3;
  size_t retval;

  retval = fread( &aux[0], sizeof(unsigned char), 8, chkpnt);
  if( (retval != 1) && (retval != 8*sizeof(unsigned char)))
    return retval;
  *v0 = (UL) aux[0];
  u1 = ((UL) aux[1])<<8;
  u2 = ((UL) aux[2])<<16;
  u3 = ((UL) aux[3])<<24;
  *v0 +=( u1 + u2 + u3);
  *v1 = (UL) aux[4];
  u1 = ((UL) aux[5])<<8;
  u2 = ((UL) aux[6])<<16;
  u3 = ((UL) aux[7])<<24;
  *v1 += (u1 + u2 + u3);
  return 1;
}

/* this routine writes the mantissa in integer form. It returns the
   number of 8-bytes blocks wrote */
#if ULONG_MAX > 0xFFFFFFFF
/* case 64 bits */
size_t fwrite_UL_array(UL *u, size_t nitems, FILE * chkpnt)
{
  UL mask=255;
  int i;
  size_t retval;
  unsigned char  aux[8];

  for(i=0;i<nitems;i++)
    {
      aux[0]= (unsigned char) (u[i] & mask);
      aux[1]= (unsigned char) ((u[i]>>8) & mask);
      aux[2]= (unsigned char) ((u[i]>>16) & mask);
      aux[3]= (unsigned char) ((u[i]>>24) & mask);
      aux[4]= (unsigned char) ((u[i]>>32) & mask);
      aux[5]= (unsigned char) ((u[i]>>40) & mask);
      aux[6]= (unsigned char) ((u[i]>>48) & mask);
      aux[7]= (unsigned char) ((u[i]>>56) & mask);
      retval= fwrite (&aux[0],1,8, chkpnt);
      if (retval != 8)
        return retval;
    }
  return nitems;
}
#else

/* case 32 bits */
size_t fwrite_UL_array(UL *u, size_t nitems, FILE * chkpnt)
{
  UL mask=255;
  size_t i, retval;
  unsigned char  aux[4];

  for(i = 0; i < nitems; i++)
    {
      aux[0]= (unsigned char) (u[i] & mask);
      aux[1]= (unsigned char) ((u[i]>>8) & mask);
      aux[2]= (unsigned char) ((u[i]>>16) & mask);
      aux[3]= (unsigned char) ((u[i]>>24) & mask);
      retval= fwrite (&aux[0],1,4, chkpnt);
      if (retval != 4)
        return retval;
    }
  /* to complete a 8-bytes block */
  if (nitems & 1)
    {
      aux[0]= 0;
      aux[1]= 0;
      aux[2]= 0;
      aux[3]= 0;
      retval= fwrite (&aux[0],1,4, chkpnt);
      if (retval != 4)
        return retval;
      return ((nitems + 1)>>1);
    }
  return (nitems >> 1);
}
#endif

#if ULONG_MAX > 0xFFFFFFFF
size_t fread_UL_array(UL * u, size_t nitems, FILE * chkpnt)
{
  unsigned char  aux[8];
  UL u0,u1,u2,u3;
  int i;
  size_t retval;

  for(i=0;i<nitems;i++)
    {
      retval = fread( &aux[0], sizeof(unsigned char), 8, chkpnt);
      if( (retval != 8) && (retval != 1))
        return retval;
      u0 = (UL) aux[0];
      u1 = ((UL) aux[1])<<8;
      u2 = ((UL) aux[2])<<16;
      u3 = ((UL) aux[3])<<24;
      u0 += ((UL) aux[4])<<32;
      u1 += ((UL) aux[5])<<40;
      u2 += ((UL) aux[6])<<48;
      u3 += ((UL) aux[7])<<56;
      u[i] = (u0 + u1 + u2 + u3);
    }
  return nitems;
}
#else
size_t fread_UL_array(UL * u, size_t nitems, FILE * chkpnt)
{
  unsigned char  aux[4];
  UL u0,u1,u2,u3;
  size_t i, retval;

  for(i = 0; i < nitems; i++)
    {
      retval = fread( &aux[0], 1, 4, chkpnt);
      if( retval != 4)
        return retval;
      u0 = (UL) aux[0];
      u1 = ((UL) aux[1])<<8;
      u2 = ((UL) aux[2])<<16;
      u3 = ((UL) aux[3])<<24;
      u[i] = (u0 + u1 + u2 + u3);
    }
  if(nitems & 1) /* complete 8-byte blocks */
    {
      retval = fread( &aux[0], 1, 4, chkpnt);
      if( retval != 4)
        return retval;
      return ((nitems+1)>>1);
    }
  return (nitems >> 1);
}
#endif

void add_sumcheck(UL * total, UL partial)
{
  *total += partial;
  if (*total < partial)
    *total = *total +1;
}

/* This is to add mod (2^2*BITS_PER_UL-1) */
void add_sumcheck2(UL *total_high,UL *total_low, UL high, UL low)
{
  UL cl,ch;
  *total_low += low;
  cl=(*total_low < low);
  *total_high += high;
  ch=(*total_high < high);
  if(cl)
    {
      *total_high += cl;
      ch += (*total_high == 0);
    }
  if(ch)
    {
      *total_low += ch;
      cl = (*total_low == 0);
      *total_high += cl;
    }
}


UL sumcheck32( UL sumcheck)
{
#if ULONG_MAX > 0xFFFFFFFF
  UL aux1,aux2;
  aux1 = (sumcheck >> 32);
  aux2 = (sumcheck & 0xFFFFFFFF);
  aux1+=aux2;
  if(aux1 > 0xFFFFFFFF)
    aux1 -= 0xFFFFFFFF;
  return aux1;
#else

  return sumcheck;
#endif
}












