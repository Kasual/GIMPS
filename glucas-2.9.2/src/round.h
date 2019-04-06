/*$Id$*/
/*  round.h. An collection of routines to make fast selftest
 
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
/********  The TRICK to round to nearest ******************************
   This plays with the internal hardward round to nearest when we add a
small number to a very big number, called bigA number. For intel fpu
this is 3*2^62. Then you have to substract the same number (named different 
to confuse the clever compilers). This makes rint() very fast.
 
If you want to use this, define TRICKY_ROUND. It is recommended to do so.
When init the FFT package on first call, it detects whether there is a 
valid trick constant and assign then properly 
***********************************************************************/
#define TRICKY_ROUND
#ifdef TRICKY_ROUND
# if (defined(__INTEL_COMPILER) && defined(__ICC))
#  define RINT(x) ((long long) (x))
# else

#  define RINT(x) (((x) + A ) - B)
/*
#  define RINT(x) ({ double __arg ;\
   __asm__ volatile ("frndint \n" : "=t" (__arg) : "0" (x));\
   __arg;})
*/
# endif

extern BIG_DOUBLE bigA, bigB;

# if defined(Y_USE_SSE2)
extern Y__M128D MM_bigA, MM_bigB;

/* these are short small assembler macros valid for both GCC and Intel
   compiler */

#  if defined(Y_SIMULATE_SSE2)
#   define Y_MM_RINT_PART0_PD( _res, _op) Y_MM_ADD_PD(_res, _op, MM_bigA)

#   define Y_MM_RINT_PART1_PD( _res, _op) Y_MM_SUB_PD(_res, _op, MM_bigB)

#   if defined(__GNUC__)
#    define Y_MM_RINT_PD( _res, _op) _res.d[0] = ((_op.d[0] + bigA) - bigB);  \
     _res.d[1] = ((_op.d[1] + bigA) - bigB);

#    define Y_MM_RINT_SD( _res, _op) _res.d[0] = ((_op.d[0] + bigA) - bigB);

#   elif defined(__ICC)
#    define Y_MM_RINT_PD( _res, _op) _res.d[0] = (long long)(_op.d[0]);  \
     _res.d[1] = (long long)(_op.d[1]);

#    define Y_MM_RINT_SD( _res, _op) _res.d[0] = (long long)(_op.d[0]);

#   endif

#   define Y_MM_TWO_RINT_PD( _res0, _res1, _op0, _op1) \
   _res0.d[0] = ((_op0.d[0] + bigA) - bigB);         \
   _res0.d[1] = ((_op0.d[1] + bigA) - bigB);         \
   _res1.d[0] = ((_op1.d[0] + bigA) - bigB);         \
   _res1.d[1] = ((_op1.d[1] + bigA) - bigB);


#  else /* Y_SIMULATE_SSE2 */
/* First part, add bigA */

#   define Y_MM_RINT_PART0_PD( _res, _op)   \
 __asm__ volatile ("movapd %1, %0 \n"       \
"        addpd %2, %0"                      \
: "=&x" (_res) : "x" (_op), "m"(MM_bigA))


/* Second part, substract bigA */

#   define Y_MM_RINT_PART1_PD( _res)                 \
 __asm__ volatile ("subpd %1, %0"                    \
                : "+&x" (_res) : "x"(MM_bigA))


/* divided in two parts */

#   define Y_MM_RINT2_PD( _res, _op)\
    Y_MM_RINT_PART0_PD( _res, _op); \
    Y_MM_SUB_PD( _res, _res, MM_bigA);

#   if (defined(Y_AMD64) || defined(Y_PENTIUM4)) && defined(__GNUC__)

#    define Y_MM_RINT_PD( _res, _op)                     \
     __asm__ volatile ("addpd %2, %0\n"                  \
                   "        subpd %2, %0\n"              \
                : "=&x" (_res) : "0" (_op), "X"(MM_bigA))

#    define Y_MM_RINT_SD( _res, _op)                     \
     __asm__ volatile ("addsd %2, %0 \n"                 \
                   "        subsd %2, %0"                \
                : "=&x" (_res) : "0" (_op), "X"(MM_bigA))

#    define Y_MM_RINT_EQ_PD( _res)                       \
     __asm__ volatile ("addpd %1, %0\n"                  \
                   "        subpd %1, %0\n"              \
                : "+&x" (_res) : "X"(MM_bigA))


#   else

#    define Y_MM_RINT_PD( _res, _op)                     \
     __asm__ volatile ("movapd %1, %0 \n"                \
                   "        addpd %2, %0\n"              \
                   "        subpd %2, %0\n"              \
                : "=&x" (_res) : "x" (_op), "m"(MM_bigA))

#    define Y_MM_RINT_SD( _res, _op)                     \
     __asm__ volatile ("movapd %1, %0 \n"                \
                   "        addsd %2, %0 \n"             \
                   "        subpd %2, %0"                \
                : "=&x" (_res) : "x" (_op), "m"(MM_bigA))

#    define Y_MM_RINT_EQ_PD( _res)                       \
     __asm__ volatile (" addpd %1, %0\n"                 \
                   "        subpd %1, %0\n"              \
                : "+&x" (_res) : "m"(MM_bigA))


#   endif /* Y_AMD64 */



#   if defined(__ICC)
/* This is required because of bad support of Intel C++ compiler
on GNU-style __asm__ */

#   define Y_MM_TWO_RINT_PD( _res0, _res1, _op0, _op1) \
 __asm__ volatile (" movapd MM_bigA, %%xmm7 \n"        \
                   "        movapd %2, %0 \n"          \
                   "        movapd %3, %1 \n"          \
                   "        addpd %%xmm7, %0 \n"       \
                   "        addpd %%xmm7, %1 \n"       \
                   "        subpd %%xmm7, %0 \n"       \
                   "        subpd %%xmm7, %1"          \
                   : "=&x" (_res0) , "=&x" (_res1) :   \
		   "x" (_op0), "x" (_op1) : "xmm7" )

#   else /* __ICC */

#    if defined(Y_AMD64)

#     define Y_MM_TWO_RINT_PART1_PD( _res0, _res1, _op0, _op1) \
       __asm__ volatile ("addpd %4, %0 \n"           \
                   "        addpd %4, %1"           \
                   : "=&x" (_res0) , "=&x" (_res1) :   \
		   "0" (_op0), "1" (_op1), "X"(MM_bigA))

#     define Y_MM_TWO_RINT_PART2_PD( _res0, _res1) \
       __asm__ volatile ("subpd %2, %0 \n"           \
                   "        subpd %2, %1"              \
                   : "+&x" (_res0) , "+&x" (_res1) :   \
		   "X"(MM_bigA))

#     define Y_MM_TWO_RINT_PD( _res0, _res1, _op0, _op1) \
      Y_MM_TWO_RINT_PART1_PD (_res0, _res1, _op0, _op1); \
      Y_MM_TWO_RINT_PART2_PD (_res0, _res1)

/*
#     define Y_MM_TWO_RINT_PD( _res0, _res1, _op0, _op1) \
       __asm__ volatile ("addpd %4, %0 \n"           \
                   "        addpd %4, %1 \n"           \
                   "        subpd %4, %0 \n"           \
                   "        subpd %4, %1"              \
                   : "=&x" (_res0) , "=&x" (_res1) :   \
		   "0" (_op0), "1" (_op1), "X"(MM_bigA))
*/
#     define Y_MM_TWO_RINT_EQ_PD( _res0, _res1)        \
       __asm__ volatile ("addpd %2, %0 \n"             \
                   "        addpd %2, %1 \n"           \
                   "        subpd %2, %0 \n"           \
                   "        subpd %2, %1"              \
                   : "+&x" (_res0) , "+&x" (_res1) :   \
		   "X"(MM_bigA))

#    else

#     define Y_MM_TWO_RINT_PD( _res0, _res1, _op0, _op1) \
       __asm__ volatile ("movapd %2, %0 \n"              \
                   "        movapd %3, %1 \n"            \
                   "        addpd %4, %0 \n"             \
                   "        addpd %4, %1 \n"             \
                   "        subpd %4, %0 \n"             \
                   "        subpd %4, %1"                \
                   : "=&x" (_res0) , "=&x" (_res1) :     \
		   "x" (_op0), "x" (_op1), "m"(MM_bigA))

#     define Y_MM_TWO_RINT_EQ_PD( _res00, _res11) \
      Y_MM_TWO_RINT_PD(_res00, _res11, _res00, _res11)

#    endif /* Y_AMD64 */
#   endif /* __ICC */
#  endif /* Y_SIMULATE_SSE2 */
# endif /* Y_USE_SSE2 */
#else
# define RINT(x) floor((x) + 0.5)
#endif /*TRICKY_ROUND */
