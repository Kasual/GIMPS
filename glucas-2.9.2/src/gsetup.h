/*$Id$*/
/*
 
   (c) 2000-2006 Guillermo Ballester Valor, Klaus Kastens
 
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
   
   This file has derived from the W. Edgington one with modifications
   from the authors above
 
*/

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#if ((defined(__MWERKS__) && defined(__INTEL__)) || defined(_MSC_VER) || defined(__WATCOMC__) || defined(__TURBOC__)) && !defined(HAVE_CONFIG_H)
# define pccompiler 1
#endif


# define PROTO(x) x

typedef unsigned short US;

/* NBBY is Number of Bits per BYte, usually 8 */
#ifndef NBBY
# define NBBY 8
#endif


#ifndef HAVE_CONFIG_H
# define BYTES_PER_UL sizeof(unsigned long)
# define BYTES_PER_BIG_DOUBLE sizeof(double)
#else
# define BYTES_PER_UL SIZEOF_UNSIGNED_LONG_INT
# define BYTES_PER_BIG_DOUBLE sizeof(double)
#endif
#define BITS_PER_UL (BYTES_PER_UL*NBBY)

/*#if !defined(vms) */
# if defined(__sgi) && !defined(__sgi__)
/* SGI's cc uses __sgi; gcc on SGI's uses __sgi__ */
#  define __sgi__ (__sgi)
/* avoid -xansi's inlining that fails to set errno correctly in some cases */
#  undef __INLINE_INTRINSICS
# endif

# include <stdlib.h>


# include <math.h>

# ifdef sun
#  include "prof.h"
# endif

# ifdef hpux
/* some HPUX compilers in ANSI mode only #define __hpux, so use it */
#  ifndef __hpux
#   define __hpux (hpux)
#  endif
# endif

# include <stdio.h>
# include <ctype.h>
# include <assert.h>
# include <string.h>
# include <signal.h>
# include <errno.h>
/*#endif*/

#if !defined(__APPLE__) && !defined(__FreeBSD__) && !defined(__OpenBSD__) && !defined(pc7300) && !defined(mips) && !defined(linux) && !defined(_AIX)
# if !defined(__ultrix) && !defined(macintosh) && !defined(__hpux) && !defined(pccompiler) && !defined(__osf__) && !defined(__sun)
typedef int handler;
#  define return_handler return(0)
# endif
#endif
#ifndef return_handler
typedef void handler;
# define return_handler return
#endif


#ifndef __STDC__
# define volatile /* */
#endif

extern volatile sig_atomic_t terminate;	/* Flag: have we gotten a SIGTERM ? */

#ifndef errno
extern int errno;
#endif

#if defined(IBMRTAIX) || defined(pc7300)
# define SIGCHLD SIGCLD
#endif

extern void setup PROTO((int nonice));

#ifndef vms

# if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(linux) || defined(__ultrix) || defined(_AIX) || defined(__hpux) || defined(macintosh) || (defined(__MWERKS__) && defined(__INTEL__)) || defined(_MSC_VER) || defined(__osf__) || defined(__sun)
extern handler term_handler PROTO((int huh));
# else
extern handler term_handler PROTO((void));
# endif
#endif














