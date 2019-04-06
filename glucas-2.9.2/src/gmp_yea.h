/*$Id$*/
/*
  YEAFFT. A library to make real convolutions using Fast Fourier
  Transforms. 
  Copyright (C) 2000-2006 Guillermo Ballester Valor,
  
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

/* These are utils definition to use YEAFFT with GMP */


#ifdef USE_GMP
# include "gmp.h"
#endif

#ifdef FFT_CHECK
# if defined(HAVE_LIMITS_H) || !defined(HAVE_CONFIG_H)
#  include <limits.h>
#  if defined(SUN_V9_GCC) || defined(__sparc)
#   undef ULONG_MAX
#   define ULONG_MAX 0xFFFFFFFF
#  endif
#  if ULONG_MAX == 0xFFFFFFFF
#   define W_TYPE_SIZE 32
#  else
#   define W_TYPE_SIZE 64
#  endif
# endif
# if defined(SIZEOF_UNSIGNED_LONG_INT)
#  ifdef W_TYPE_SIZE
#   undef W_TYPE_SIZE
#  endif
#  if SIZEOF_UNSIGNED_LONG_INT == 4
#   define W_TYPE_SIZE 32
#  elif SIZEOF_UNSIGNED_LONG_INT == 8
#   define W_TYPE_SIZE 64
#  endif
# endif
#endif
typedef unsigned long int UWtype;
typedef unsigned int UHWtype;
typedef unsigned long int UDWtype;
typedef unsigned char UQItype;
typedef int SItype;
typedef unsigned int USItype;
typedef long int DItype;
typedef unsigned long int UDItype;

/* Some utils definitions */
#ifndef USE_GMP
typedef unsigned long mp_limb_t;
typedef unsigned long * mp_ptr;
#endif
#define BYTES_PER_LIMB sizeof(mp_limb_t)
#define BITS_PER_LIMB (BYTES_PER_LIMB*8)

/*$Id$*/
