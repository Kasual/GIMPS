/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2005 Guillermo Ballester Valor
 
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
/* file includes macros to tune the code properly for 8086 family
  procesors with gnu gcc compiler */
/*#define Y_DEBUG*/

#if (defined(__GNUC__) || defined (__ICC)) && (defined(__i386__) || defined(__i486__) || defined (__i586__) || defined(__i686__)) && defined(Y_DEBUG)

/*# include "longlong.h"*/
# include <assert.h>

/* from longlong.h */
#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  __asm__ ("addl %5,%1\n\tadcl %3,%0"					\
	   : "=r" ((USItype)(sh)), "=&r" ((USItype)(sl))		\
	   : "%0" ((USItype)(ah)), "g" ((USItype)(bh)),			\
	     "%1" ((USItype)(al)), "g" ((USItype)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  __asm__ ("subl %5,%1\n\tsbbl %3,%0"					\
	   : "=r" ((USItype)(sh)), "=&r" ((USItype)(sl))		\
	   : "0" ((USItype)(ah)), "g" ((USItype)(bh)),			\
	     "1" ((USItype)(al)), "g" ((USItype)(bl)))

# define DEBUG_VARS  mp_limb_t he, le, hb, lb;

# define DEBUG_GLOBAL_VARS extern unsigned int TP_TIMEH, TP_TIMEL, TP_COUNT, \
	TP_TIMEMAX, TP_TIMEMIN, TP_COUNTMIN;

# define RDTSC(_h, _l) {\
	__asm__ __volatile__ (".byte 0x0F; .byte 0x31"\
	 : "=a" (_l),"=d" (_h) : );}

# define BEGIN_TIME  RDTSC(hb, lb);

# define END_TIME RDTSC(he, le); \
	sub_ddmmss(he, le, he, le, hb, lb);\
        if(le < TP_TIMEMIN) {TP_TIMEMIN = le; TP_COUNTMIN= TP_COUNT;}\
        if(le > TP_TIMEMAX) TP_TIMEMAX = le;\
        if(le < (TP_TIMEMIN << 3) && le < 20000) {\
	add_ssaaaa (TP_TIMEH, TP_TIMEL, he, le, TP_TIMEH, TP_TIMEL);\
	TP_COUNT++;}


DEBUG_GLOBAL_VARS;

/* This is from GNU/FFTW package, a trick to see whether doubles on stack
 are aligned */
# define ASSERT_ALIGNED_DOUBLE() {                       \
     double __foo;                                       \
     assert((((long) &__foo) & 0x7)==0);                 \
}\


# define see_align_stack printf(" %8.8X \n",(int)(&t0r));

#elif (defined(__GNUC__) || defined (__ICC)) && (defined(Y_AMD64)) && defined(Y_DEBUG)

/* from longlong.h */
#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  __asm__ ("addl %5,%1\n\tadcl %3,%0"					\
	   : "=R" ((USItype)(sh)), "=&R" ((USItype)(sl))		\
	   : "%0" ((USItype)(ah)), "g" ((USItype)(bh)),			\
	     "%1" ((USItype)(al)), "g" ((USItype)(bl)))

#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  __asm__ ("subl %5,%1\n\tsbbl %3,%0"					\
	   : "=R" ((USItype)(sh)), "=&R" ((USItype)(sl))		\
	   : "0" ((USItype)(ah)), "g" ((USItype)(bh)),			\
	     "1" ((USItype)(al)), "g" ((USItype)(bl)))


# define DEBUG_VARS  USItype he, le, hb, lb;

# define DEBUG_GLOBAL_VARS extern unsigned int TP_TIMEH, TP_TIMEL, TP_COUNT, \
	TP_TIMEMAX, TP_TIMEMIN, TP_COUNTMIN;

# define RDTSC(_h, _l) {\
	__asm__ __volatile__ (".byte 0x0F; .byte 0x31"\
	 : "=a" (_l),"=d" (_h) : );}

# define BEGIN_TIME  RDTSC(hb, lb);

# define END_TIME RDTSC(he, le); \
	sub_ddmmss(he, le, he, le, hb, lb);\
        if(le < TP_TIMEMIN) {TP_TIMEMIN = le; TP_COUNTMIN= TP_COUNT;}\
        if(le > TP_TIMEMAX) TP_TIMEMAX = le;\
        if(le < (TP_TIMEMIN << 3) && le < 20000) {\
	add_ssaaaa (TP_TIMEH, TP_TIMEL, he, le, TP_TIMEH, TP_TIMEL);\
	TP_COUNT++;}


DEBUG_GLOBAL_VARS;

/* This is from GNU/FFTW package, a trick to see whether doubles on stack
are aligned */
# define ASSERT_ALIGNED_DOUBLE() /* */

#else
# define DEBUG_VARS /* */
# define BEGIN_TIME /* */
# define END_TIME  /* */
# define ASSERT_ALIGNED_DOUBLE() /* */
#endif

/*$Id$*/
