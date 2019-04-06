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
*/
/*
   This file has derived from the W. Edgington one with modifications
   from the authors above
 
*/
#include "gsetup.h"

/* Nice value, default = 40. It is changed by -N nice option */
int Nice = 40;

#ifndef vms
volatile sig_atomic_t terminate = 0; /* Flag: have we gotten a SIGTERM ? */
#endif

#if (defined(__hpux) || defined(mips)) && defined(TCPIP)
int h_errno = 0;
#endif

/* only used to make the process's priority as low as possible */
#if !defined(macintosh) && !defined(pccompiler)
extern int nice PROTO((int priority_increment));
#endif

#ifdef __sgi__
#include <limits.h>
#include <sys/types.h>
#include <sys/prctl.h>
#include <sys/schedctl.h>
#endif

void setup (int nonice)
{
  /*setlinebuf(stdout);
  setlinebuf(stderr);*/
  if (!nonice)
    {
#ifdef __sgi__
      (void) schedctl (NDPRI, 0, NDPLOMIN);
#endif
#if !defined(macintosh) && !defined(pccompiler)

      (void) nice (Nice);
#endif

    }
#ifdef SIGTERM
  (void) signal (SIGTERM, term_handler);
#endif
#ifdef SIGINT

  (void) signal (SIGINT, term_handler);
#endif
#ifdef SIGXCPU

  (void) signal (SIGXCPU, term_handler);
#endif
#ifdef SIGHUP

  if(signal (SIGHUP,SIG_IGN) != SIG_IGN)
    (void) signal (SIGHUP, term_handler);
#endif
#ifdef SIGPIPE

  (void) signal (SIGPIPE, SIG_IGN);
#endif
#ifdef linux
  /*!!really ought to examine the signal and let it happen if it's an overflow, NaN, etc.,*/
  /*!! but it's almost always an underflow or round off warning rather than a problem like overflow or NaN */
  (void) signal (SIGFPE, SIG_IGN);
#endif

  return;
}

#if defined(__MWERKS__)
#pragma warn_unusedarg off
#endif

#if !defined(vms)
/*ARGSUSED*/ /* they aren't but that's because only some OS's want an argument here */
# if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(linux) || defined(__ultrix) || defined(_AIX) || defined(__hpux) || defined(macintosh) || (defined(__MWERKS__) && defined(__INTEL__)) || defined(_MSC_VER) || defined(__osf__) || defined(__sun)
handler term_handler(int huh)
# else
handler term_handler PROTO((void))
# endif
{
# ifdef SIGTERM
  (void) signal (SIGTERM, SIG_IGN);
# endif
# ifdef SIGINT

  (void) signal (SIGINT, SIG_IGN);
# endif
# ifdef SIGXCPU

  (void) signal (SIGXCPU, SIG_IGN);
# endif
# ifdef SIGHUP

  (void) signal (SIGHUP, SIG_IGN);
# endif

  terminate = 1;
  return_handler;
}
#endif

#if defined(__MWERKS__)
#pragma warn_unusedarg reset
#endif
