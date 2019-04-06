/* $Id$ */
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006 Guillermo Ballester Valor, Klaus kastens,
                  Alexander Kruppa
 
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
/* This include file has the prefetch hints for some specific targets
 
   The targets defined are:
 
   0: Generic. No prefetch. 100% pure C-code.
 
   FOR GNU/GCC compiler, or other accepting assembler extensions
   
 
   1: Pentium-plain. Pentiummx. Pentium-pro. No prefetch. A lot of assembler
      macros. GNU/gcc syntax
   11: Pentium 3. Prefetch. A lot of assembler macros. GNU/gcc syntax
   12: Athlon. Prefetch. A lot of assembler macros. GNU/gcc syntax
   13: Pentium 4. Prefetch. A lot of assembler macros. GNU/gcc syntax
   16: Pentium 3. Prefetch. Only two lines of assembler. GNU/gcc syntax
   17: Athlon. Prefetch. Only two lines of assembler. GNU/gcc syntax
   21: PowerPC 601. Prefetch. Only two lines of assembler. 
       GNU/gcc or Metrowerks Codewarrior syntax
   23: PowerPC 604e, 750. Prefetch. Only two lines of assembler. 
       GNU/gcc or Metrowerks Codewarrior syntax
   24: PowerPC 970 aka G5. GNU/gcc or Metrowerks Codewarrior syntax
   31: Alpha ev56, ev6, ev67. Prefetch. Only two lines of assembler. 
       Compaq-c syntax
   32: Alpha ev56, ev6, ev67. Prefetch. Only two lines of assembler. 
       GNU/gcc syntax
   41: Ultrasparc v9: Prefetch. Only two lines of assembler. GNU/gcc syntax
   51: Itanium. Intel IA-64. Prefetch. Only two lines of assembler. 
       GNU/gcc syntax
   61: PA-RISC 2.0: Prefetch. Only two lines of assembler. GNU/gcc syntax
*/

/* some routines can include fast loops and we may need to increase the
   prefetch ahead value in Y_INC_AHEAD doubles */

#ifndef Y_INC_AHEAD
# define Y_INC_AHEAD 0
#endif

/* the cache line is 32 bytes by default */
#ifndef Y_CACHE_LINE
# define Y_CACHE_LINE 4
#endif

#define Y_PREFETCH_AHEAD (Y_CACHE_LINE + Y_INC_AHEAD)
#if Y_TARGET == 0

/* Prefetch is not defined for generic code */
# if defined (HAVE_CONFIG_H) &&  (HAVE_DECL___BUILTIN_PREFETCH == 1)
#  define Y_PREFETCH
#  define prefetch_data(_array,_k) __builtin_prefetch(_array + _k, 1, 3)

#  define prefetch_p(_pd) __builtin_prefetch(_pd + Y_PREFETCH_AHEAD, 1, 3)

#  define prefetch_data_trig(_array,_k) __builtin_prefetch(_array + _k, 0, 3)

#  define prefetch_p_trig(_pd) __builtin_prefetch(_pd + Y_PREFETCH_AHEAD, 0, 3)

#  define prefetch_data_trig_nta(_array,_k) __builtin_prefetch(_array + _k, 0, 0)

#  define prefetch_p_trig_nta(_pd) __builtin_prefetch(_pd + Y_PREFETCH_AHEAD, 0, 0)

# elif defined(__ICC)

#  include <xmmintrin.h>

#  define Y_PREFETCH

#  define prefetch_data(_array,_k) _mm_prefetch((char *)(_array + _k), _MM_HINT_T2)

#  define prefetch_p(_pd) _mm_prefetch((char *)(_pd + Y_PREFETCH_AHEAD), _MM_HINT_T2)

#  define prefetch_data_trig_nta(_array,_k) _mm_prefetch((char *)(_array + _k), _MM_HINT_NTA)

#  define prefetch_p_trig_nta(_pd) _mm_prefetch((char *)(_pd + Y_PREFETCH_AHEAD), _MM_HINT_NTA)

#  define prefetch_data_trig(_array,_k) _mm_prefetch((char *)(_array + _k), _MM_HINT_T2)

#  define prefetch_p_trig(_pd) _mm_prefetch((char *)(_pd + Y_PREFETCH_AHEAD), _MM_HINT_T2)

# else

#  define prefetch_data(_array,_k) /* */

#  define prefetch_p(_pd) /* */

#  define prefetch_data_trig(_array,_k) /* */

#  define prefetch_p_trig(_pd) /* */

#  define prefetch_data_trig_nta(_array,_k) /* */

#  define prefetch_p_trig_nta(_pd) /* */

# endif

/* Pentium 3 family with lot of assembler*/
#elif (Y_TARGET == 11)

# define Y_PREFETCH

# ifndef USE_ASM

#  define prefetch_data(_array, _k) \
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p(_pd)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# else /* USE_ASM */

#  define prefetch_data(_array,_k) /* */

#  define prefetch_p(_pd) /* */

#  define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# endif /* USE_ASM */


/* Athlons with a lot of assembler lines */
#elif (Y_TARGET == 12)
# define Y_PREFETCH
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8
# ifndef USE_ASM

#  define prefetch_data(_array,_k)\
    __asm__ volatile ("prefetchw (%0)": :"r"(_array + _k));

#  define prefetch_p(_pd)\
    __asm__ volatile ("prefetchw (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# else /* USE_ASM */

#  define prefetch_data(_array,_k) /* */

#  define prefetch_p(_pd) /* */

#  define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# endif /* USE_ASM */


/* Pentium 4 family with lot of assembler*/
#elif (Y_TARGET == 13)
# define Y_PREFETCH
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 16
# ifndef USE_ASM
#  define prefetch_data(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p(_pd)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# else /* USE_ASM */

#  define prefetch_data(_array,_k) /* */

#  define prefetch_p(_pd) /* */

#  define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

#  define prefetch_p_trig(_pd)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

#  define prefetch_p_trig_nta(_pd)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# endif /* USE_ASM */


/* Pentium 3 with only prefetch hints in assembler */
#elif (Y_TARGET == 16)
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 4

# define prefetch_data(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

# define prefetch_p(_pd)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

# define prefetch_p_trig(_pd)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

# define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));


/* Athlons with only prefetch hints in assembler */
#elif (Y_TARGET == 17)
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8

# define prefetch_data(_array,_k)\
    __asm__ volatile ("prefetchw (%0)": :"r"(_array + _k));

# define prefetch_p(_pd)\
    __asm__ volatile ("prefetchw (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_array + _k));

# define prefetch_p_trig(_pd)\
    __asm__ volatile ("prefetcht2 (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_array + _k));

# define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("prefetchnta (%0)": :"r"(_pd + Y_PREFETCH_AHEAD));


/* PowerPC 601, 604e, 750 */
#elif (Y_TARGET == 21) || (Y_TARGET == 23) || (Y_TARGET == 24)

# undef Y_CACHE_LINE
# if Y_TARGET == 21
#  define Y_CACHE_LINE 8
# elif Y_TARGET == 23
#  define Y_CACHE_LINE 4
# elif Y_TARGET == 24
#  define Y_CACHE_LINE 16
# endif
# define Y_PREFETCH

# ifdef __GNUC__
#  define prefetch_data(_array,_k)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_array + _k), "r"(0));

#  define prefetch_p(_pd)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_pd), "r"(Y_PREFETCH_AHEAD * sizeof(double)));

#  define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_array + _k), "r"(0));

#  define prefetch_p_trig(_pd)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_pd), "r"(Y_PREFETCH_AHEAD * sizeof(double)));

#  define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_array + _k), "r"(0));

#  define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("dcbt %0,%1": :"%b"(_pd), "r"(Y_PREFETCH_AHEAD * sizeof(double)));

# endif

# ifdef __MWERKS__
#  define prefetch_data(_array,_k) __dcbt(_array, _k);
#  define prefetch_p(_pd) __dcbt(_pd, Y_PREFETCH_AHEAD * sizeof(double));
#  define prefetch_data_trig(_array, _k) __dcbt(_array, _k);
#  define prefetch_p_trig(_pd) __dcbt(_pd, Y_PREFETCH_AHEAD * sizeof(double));
#  define prefetch_data_trig_nta(_array,_k) __dcbt(_array, _k);
#  define prefetch_p_trig_nta(_pd) __dcbt(_pd, Y_PREFETCH_AHEAD * sizeof(double));
# endif


/* Alpha ev56 and above : Compaq-c*/
#elif (Y_TARGET == 31)
# include <c_asm.h>

# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8
# define Y_PREFETCH

# define prefetch_data(_array,_k) asm("lds %f31,0(%a0)",(_array + _k));

# define prefetch_p(_pd) asm("lds %f31,64(%a0)",(_pd + Y_INC_AHEAD));

# define prefetch_data_trig(_array,_k) asm("ldf %f31,0(%a0)",(_array + _k));

# define prefetch_p_trig(_pd) asm("ldf %f31,64(%a0)",(_pd + Y_INC_AHEAD));

# define prefetch_data_trig_nta(_array,_k) asm("ldf %f31,0(%a0)",(_array + _k));

# define prefetch_p_trig_nta(_pd) asm("ldf %f31,64(%a0)",(_pd + Y_INC_AHEAD));

/* Alpha ev56 and above */
#elif (Y_TARGET == 32)

# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8
# define Y_PREFETCH

# define prefetch_data(_array,_k)\
   __asm__ volatile ("lds $f31,0(%0)": :"r"(_array + _k));

# define prefetch_p(_pd)\
   __asm__ volatile ("lds $f31,0(%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig(_array,_k)\
   __asm__ volatile ("ldf $f31,0(%0)": :"r"(_array + _k));

# define prefetch_p_trig(_pd)\
   __asm__ volatile ("ldf $f31,0(%0)": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig_nta(_array,_k)\
   __asm__ volatile ("ldf $f31,0(%0)": :"r"(_array + _k));

# define prefetch_p_trig_nta(_pd)\
   __asm__ volatile ("ldf $f31,0(%0)": :"r"(_pd + Y_PREFETCH_AHEAD));


/* Ultrasparc v9 */
#elif (Y_TARGET == 41)

# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8
# define Y_PREFETCH

# if defined(__sun) && defined(__sparcv9)
#  include <sun_prefetch.h>
#  define prefetch_data(_array,_k)\
       sparc_prefetch_write_once (_array + _k);
#  define prefetch_p(_pd) \
       sparc_prefetch_write_once (_pd + Y_PREFETCH_AHEAD);
#  define prefetch_data_trig(_array,_k)\
       sparc_prefetch_read_once (_array + _k);
#  define prefetch_p_trig(_pd)\
       sparc_prefetch_read_once (_pd + Y_PREFETCH_AHEAD);
#  define prefetch_data_trig_nta(_array,_k)\
       sparc_prefetch_read_once (_array + _k);
#  define prefetch_p_trig_nta(_pd)\
       sparc_prefetch_read_once (_pd + Y_PREFETCH_AHEAD);
# else

#  define prefetch_data(_array,_k)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_array + _k));

#  define prefetch_p(_pd)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_array+_k));

#  define prefetch_p_trig(_pd)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_pd + Y_PREFETCH_AHEAD));

#  define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_array+_k));

#  define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("prefetch [%0],3": :"g"(_pd + Y_PREFETCH_AHEAD));

/* Preload better than prefetch?
#  define prefetch_data(_array,_k)\
    __asm__ volatile ("ldd [%0],%%f30": :"g"(_array+_k):"f30");

#  define prefetch_p(_pd)\
    __asm__ volatile ("ldd [%0+16],%%f30": :"g"(_pd):"f30");
*/

# endif

/* Itanium IA-64 */
#elif (Y_TARGET == 51)
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8
# define Y_PREFETCH

# define prefetch_data(_array,_k)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_array+_k));

# define prefetch_p(_pd)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_array+_k));

# define prefetch_p_trig(_pd)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_array+_k));

# define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("lfetch.nta [%0]": :"r"(_pd + Y_PREFETCH_AHEAD));


/* PA-RISC 2.0
According to "PA-RISC 2.0 Architecture", Table 6-10:
ldd - Prefetch cache line for write (not supported by binutils 2.12)
ldw - Prefetch cache line for read */
#elif (Y_TARGET == 61)
# undef Y_CACHE_LINE
# define Y_CACHE_LINE 8

# define prefetch_data(_array,_k)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_array + _k));

# define prefetch_p(_pd)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_array + _k));

# define prefetch_p_trig(_pd)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_pd + Y_PREFETCH_AHEAD));

# define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_array + _k));

# define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("ldw 0(%0), %%r0": :"r"(_pd + Y_PREFETCH_AHEAD));


/* Any experimental */
# elif (Y_TARGET == 99)

# undef Y_CACHE_LINE
# define Y_CACHE_LINE 2
# define Y_PREFETCH

# define prefetch_data(_array,_k)\
    __asm__ volatile ("ldd [%0],%%f30": :"g"(_array+_k):"f30");

# define prefetch_p(_pd)\
    __asm__ volatile ("ldd [%0+16],%%f30": :"g"(_pd):"f30");

# define prefetch_data_trig(_array,_k)\
    __asm__ volatile ("ldd [%0],%%f30": :"g"(_array+_k):"f30");

# define prefetch_p_trig(_pd)\
    __asm__ volatile ("ldd [%0+16],%%f30": :"g"(_pd):"f30");

# define prefetch_data_trig_nta(_array,_k)\
    __asm__ volatile ("ldd [%0],%%f30": :"g"(_array+_k):"f30");

# define prefetch_p_trig_nta(_pd)\
    __asm__ volatile ("ldd [%0+16],%%f30": :"g"(_pd):"f30");


#else

# define prefetch_data(_array,_k) /* */

# define prefetch_p(_pd) /* */

# define prefetch_data_trig(_array,_k) /* */

# define prefetch_p_trig(_pd) /* */

# define prefetch_data_trig_nta(_array,_k) /* */

# define prefetch_p_trig_nta(_pd) /* */

#endif

/* $Id$ */
