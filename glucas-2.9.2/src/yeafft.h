/*$Id$*/
/*
  YEAFFT. A library to make real convolutions using Fast Fourier
  Transforms. 
  Copyright (C) 2000-2006 Guillermo Ballester Valor, Klaus Kastens
  
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

/**********************************************
   The use of OpenMP or _SUNMP is over PTHREADS
************************************************/
#if defined(_OPENMP) || defined(_SUNMP)
# ifdef _PTHREADS
#  undef _PTHREADS
# endif
#endif
#if defined(_PTHREADS)
# include <pthread.h>
#endif


/* the standar includes */
#include <stdlib.h>
#include <assert.h>
#include "gmp_yea.h"

/* Includes autogenerate yeafft flags */
#ifdef Y_DEFINES_FILE
# include "yeafft_defines.h"
#endif
/*
   This is to set we are using YEAFFT 
*/
#define YEAFFT
/*#define Y_NEW_PADDING*/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

/*
  This is to set the target machine.
  We have to write a optimized macro file for every target machine
  at the moment, we only have a generic machine   
  Would you like to write a tuned one for you machine ?
*/
#ifndef Y_TARGET
# define Y_TARGET 0
#endif

/* Set it to perform sum check */
/*#define SUM_CHECK*/

/*
  When Y_TARGET is 1 and the compiler is GNU gcc then we can use inline 
  assembler macros, and gain about 30% performance. These macros only runs with
  Y_BLOCKSIZE = 1024 and Y_SHIFT = 8
 
  There is an aditional target 11, for pentium3 and athlon with some hints to 
  prefetch data it gains about 15% or more extra performance.
*/

#if defined(__GNUC__)
# if (Y_TARGET == 1)
#  if !defined(Y_BLOCKSIZE)
#   define Y_BLOCKSIZE 1024
#  endif
#  if !defined(Y_SHIFT)
#   define Y_SHIFT 8
#  endif
#  if ((Y_BLOCKSIZE == 1024) && (Y_SHIFT == 8))
#    define USE_ASM
#  endif
# elif (Y_TARGET == 11)
#  if !defined(Y_BLOCKSIZE)
#   define Y_BLOCKSIZE 1024
#  endif
#  if !defined(Y_SHIFT)
#   define Y_SHIFT 8
#  endif
#  if ((Y_BLOCKSIZE == 1024) && (Y_SHIFT == 8))
#    define USE_ASM
#  endif
# elif (Y_TARGET == 12)
#  if !defined(Y_BLOCKSIZE)
#   define Y_BLOCKSIZE 1024
#  endif
#  if !defined(Y_SHIFT)
#   define Y_SHIFT 8
#  endif
#  if ((Y_BLOCKSIZE == 1024) && (Y_SHIFT == 8))
#    define USE_ASM
#  endif
# elif (Y_TARGET == 13)
#  if !defined(Y_BLOCKSIZE)
#   define Y_BLOCKSIZE 1024
#  endif
#  if !defined(Y_SHIFT)
#   define Y_SHIFT 8
#  endif
#  if ((Y_BLOCKSIZE == 1024) && (Y_SHIFT == 8))
#    define USE_ASM
#  endif

# endif
#endif


/*#define Y_MANY_REGISTERS*/

/*
   The max. radix reduction avalable (in powers of two) 
   if Y_AVAL == 3 then we have up to radix-8 reduction
   if Y_AVAL == 4 then we have up to radix-16 reduction
   if Y_AVAL == 5 then we have up to radix-32 reduction in middle passes
                and up to 16 on first and last radix-reduction
*/
#ifndef Y_AVAL
# define Y_AVAL 3
#endif

/*
   In this early version, the non_power_of_two radices are always
   the first radix reduction in plan, so some routines for these files
   are not used. To save space and time they will not compiled if we
   define Y_SAVE
*/
#define Y_SAVE


/*
 * FIXME: Which platforms/compilers benefit/support  using
 * 16 Byte algned double?
 */
#if defined(__GNUC__) && ( __GNUC__ < 4) && !(defined(__powerpc__) || defined(__ppc__))
# define Y_ALIGNED(_jj) __attribute__ ((aligned(_jj))) 
typedef double BIG_DOUBLE Y_ALIGNED(16);
typedef double y_limb_t Y_ALIGNED(16);
#elif defined(__GNUC__) && ( __GNUC__ == 4) && !(defined(__powerpc__) || defined(__ppc__))
# define Y_ALIGNED(_jj) __attribute__ ((aligned(_jj))) 
typedef double BIG_DOUBLE;
typedef double y_limb_t Y_ALIGNED(16);
#else
# define Y_ALIGNED(_j) /* */
typedef double BIG_DOUBLE;
typedef double y_limb_t;
#endif

/*
   This are some type used.
*/

typedef unsigned long UL;
typedef BIG_DOUBLE* y_ptr;
typedef unsigned int y_size_t;
struct y_complex
  {
    BIG_DOUBLE r;
    BIG_DOUBLE i;
  };

/*
  basic macros for sse2 extensions
 
  ALL_INTERLACED is defined when both read and write data in 
  last_dit_carry_norm_first_dif routines are assumed already 
  interlaced, so there is no need to interlace data on input/output
  This could save a lot of SSE2 moves.
*/

/* SSE2 is not used for more than 4 THREADS
#if defined(Y_USE_SSE2) && ((Y_NUM_THREADS > 4) || (_PTHREADS > 4))
# undef Y_USE_SSE2
#endif
*/

#if defined(Y_USE_SSE2)
# include "ysse2.h"
# if !defined(_OPENMP) && !defined(_PTHREADS) && !defined(_SUNMP)
/* for pass1 multithreaded we don't use SSE2, so we need legacy
    mode in pass0 and pass2 */
#  define ALL_INTERLACED
# endif
#endif

/* some util definitions */

#define TWOPI (BIG_DOUBLE)(2*3.1415926535897932384626433)
#define SQRTHALF (BIG_DOUBLE)(0.707106781186547524400844362104)
#define SQRT2 (BIG_DOUBLE)(1.414213562373095048801688724209)


/*
Y_BITS is the number of bits we can put in an element of FFT 
In a 64 bits float IEEE-359 compliant machine uses to be 19 or 20
 
#ifndef Y_BITS
# define Y_BITS 20
#endif
#define Y_MASQ ((1L<<Y_BITS)-1L)
#define Y_LIMB (1L<<Y_BITS)
#define Y_MAX (1L<<(Y_BITS-1))
*/

/****************** IMPORTANT MEMORY MANAGEMENT *********************
This definitions are used to manage memory in a eficcient way
Y_BLOCKSIZE: Is the number of complex data elements stored in memory
	      without holes.
Y_SHIFT: (Y_BLOCKSIZE >>Y_SHIFT is the extension of padding holes between 
	     blocksizes.
Y_MEM_THRESHOLD The threshold to change the scheme in FFT. Typical is 
	  The L2 cache size divided by 4*sizeof(double).
Y_SHIFT2: Is optional defined, it adds aditional padding spaces to avoid
          overlapping cache lines. 
*********************************************************************/
#define Y_PADDING
#ifdef Y_PADDING
# ifndef Y_PADDING_LEVEL
#  define Y_PADDING_LEVEL 1
# endif
# ifndef Y_BLOCKSIZE
#  define Y_BLOCKSIZE 512
# endif
# ifndef Y_SHIFT
#  define Y_SHIFT (UL)8
# endif
# define Y_MASK (~(UL)(Y_BLOCKSIZE - 1))
# if defined(Y_SHIFT2)
#  if Y_SHIFT2 < 4
#   undef Y_SHIFT2
#   define Y_SHIFT2 4
#  endif
#  define addr(_j) ((_j) + (((_j) & Y_MASK )>>Y_SHIFT) + (((_j)>>Y_SHIFT2)<<1))
# else
#  if defined(Y_NEW_PADDING)
#   if Y_PADDING_LEVEL == 1
#    define addr(_j) ((_j) + (((_j) >> Y_SHIFT1 )<< Y_SHIFT0))

#   elif Y_PADDING_LEVEL == 2
#    define addr(_j) ((_j) + (((_j) >> Y_SHIFT2 )<< Y_SHIFT0) +\
      (((_j) >> Y_SHIFT1 )<< 1))

#   elif Y_PADDING_LEVEL == 3
#    define addr(_j) ((_j) + (((_j) >> Y_SHIFT3 )<< Y_SHIFT0) +\
      (((_j) >> Y_SHIFT1 )<< 1)  +\
      (((_j) >> Y_SHIFT2 )<< 1))

#   elif Y_PADDING_LEVEL == 4
#    define addr(_j) ((_j) + (((_j) >> Y_SHIFT4 )<< Y_SHIFT0) +\
      (((_j) >> Y_SHIFT1 )<< 1)  +\
      (((_j) >> Y_SHIFT2 )<< 1)  +\
      (((_j) >> Y_SHIFT3 )<< 1))
#   endif
#  else
#   if Y_PADDING_LEVEL == 1
#    define addr(_j) ((_j) + (((_j) & Y_MASK )>>(Y_SHIFT )))

#   elif Y_PADDING_LEVEL == 2
#    define addr(_j) ((_j) + (((_j) & Y_MASK )>>(Y_SHIFT )) +\
      (((_j) & ((Y_MASK) << Y_AVAL))>>(Y_SHIFT + Y_AVAL)))

#   elif Y_PADDING_LEVEL == 3
#    define addr(_j) ((_j) + (((_j) & Y_MASK )>>(Y_SHIFT )) +\
      (((_j) & ((Y_MASK) << Y_AVAL))>>(Y_SHIFT + Y_AVAL)) +\
      (((_j) & ((Y_MASK) << 2*Y_AVAL))>>(Y_SHIFT + 2*Y_AVAL)))

#   elif Y_PADDING_LEVEL == 4
#    define addr(_j) ((_j) + (((_j) & Y_MASK )>>(Y_SHIFT )) +\
     (((_j) & ((Y_MASK) << Y_AVAL))>>(Y_SHIFT + Y_AVAL)) +\
     (((_j) & ((Y_MASK) << 2*Y_AVAL))>>(Y_SHIFT + 2*Y_AVAL))+\
     (((_j) & ((Y_MASK) << 3*Y_AVAL))>>(Y_SHIFT + 3*Y_AVAL)))
#   endif
#  endif
# endif
# define ac_r(_j) (addr(_j<<1))
# define ac_i(_j) (addr((_j<<1)+1))
#else
# define addr(_j) (_j)
# define ac_r(_j) (_j<<1)
# define ac_i(_j) ((_j<<1)+1)
# define Y_SHIFT 0
# define Y_MASK 0
#endif

#ifndef Y_MEM_THRESHOLD
# define Y_MEM_THRESHOLD 2048
#endif

/*
   The amount of precomputed trig factors is important. If we are limited
   with memory we can define Y_MINIMUM, in this manner it will only be 
   precomputed the basic power for every radix reduction, but 
   - The precision can drop dangerously and perhaps the program we will 
   need to redefine AVERBITS and so a larger FFT.
   - The cost of computing must be higher than access to memory.
 
   For an IEEE double with a 53 bits mantissa is not recommended define 
   Y_MINIMUM. For intel x86, the use of internal CPU 64 bits-mantissa can avoid
   the lack of precision. Anyway, is slower.
*/
/*#define Y_MINIMUM*/

/*
   To align memory to 64.  This is useful on intel x86 family
   To avoid core dumps, the memory is allocated in other pointer than
   the pointer used to data references. We must de-allocate memory using
   the pointer we used in alloc.
 */
/*
#if (defined(__i386__) || defined(__i486__) || defined(__i586__) || defined (__i686__) || defined(__sparc))
# define Y_ALIGN
#else
# ifdef Y_ALIGN
#   undef Y_ALIGN
# endif
#endif
*/

/*
# if !defined(__FreeBSD__) && !defined(__OpenBSD__)
#  define Y_ALLOC( _x ) (memalign ( 16 , (_x)))
#  define Y_FREE( _x ) (free((_x)))
# else
#  define Y_ALLOC(_x) (malloc((_x)))
#  define Y_FREE(_x) (free((_x)))
# endif
*/
#  define Y_ALLOC(_x) (malloc((_x)))
#  define Y_FREE(_x) (free((_x)))

#define Y_ALIGN
#ifdef Y_ALIGN
# define ALLOC_DOUBLES(_n) ((BIG_DOUBLE *)malloc((_n)*sizeof(BIG_DOUBLE)+128))
# define CALLOC_DOUBLES(_n) ((BIG_DOUBLE *)calloc((size_t)((_n) + 16),sizeof(BIG_DOUBLE)))
# define ALIGN_DOUBLES(_p) ((BIG_DOUBLE *)(((long int)(_p) | 63L)+1))
#else
# define ALLOC_DOUBLES(_n) ((BIG_DOUBLE *)malloc((_n)*sizeof(BIG_DOUBLE)))
# define CALLOC_DOUBLES(_n) ((BIG_DOUBLE *)calloc(_n,sizeof(BIG_DOUBLE)))
# define ALIGN_DOUBLES(_p) (_p)
#endif

/*
 On pentiums, alignment on stack is very important 
 This uses the GNU gcc  __builtin_alloca function to align doubles properly
 This is taken from GNU/FFTW package 
*/

#if defined(__GNUC__)
# if (defined(__i386__) || defined(__i486__) || defined(__i586__) || defined(__i686__)) && !defined(NO_HACK_ALIGN)

/* For some extrange reasons, a bunch of GCC old versions does not set properly
   esp, doubling the execution time :( This is an ugly patch for it */

#  define HACK_CHECK_ESP() { \
if (!(((long) (__builtin_alloca(0))) & 0x7)) BAD_STACK = 1;    \
}

# define HACK_ALIGN_STACK_EVEN()                         \
    if (BAD_STACK)                                        \
    {                                                     \
	if (!(((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4); \
    }                                                                     \
    else                                                                  \
    {                                                                     \
	if ((((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);  \
    }

# define HACK_ALIGN_STACK_ODD()                           \
    if (BAD_STACK)                                        \
    {                                                     \
	if ((((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4); \
    }                                                                     \
    else                                                                  \
    {                                                                     \
	if (!(((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);  \
    }

# else
#  define HACK_CHECK_ESP() /* */
#  define HACK_ALIGN_STACK_EVEN() /* */
#  define HACK_ALIGN_STACK_ODD() /* */
# endif
#else
# define HACK_CHECK_ESP() /* */
# define HACK_ALIGN_STACK_EVEN() /* */
# define HACK_ALIGN_STACK_ODD() /* */
#endif


/******************************************************************
 This interface includes the nine's proof. We are pretty sure the 
 result is correct if it passes the proof. It wastes very few cpu 
 cycles, so is a good thing to left this option active  
 
 BE_CAREFUL!, the use of YEAFFT package with MERS package to 
 perform Lucas-Lehmer test is incompatible with FFT_CHECK. This 
 proof is only useful when calling from GMP/GNU rutines 
 In this case, you will need to make yeafft with -DUSE_GMP option.
*******************************************************************/
#define DO_NOT_USE_FFT_CHECK
#ifndef DO_NOT_USE_FFT_CHECK
# define FFT_CHECK
#endif

#ifdef FFT_CHECK
typedef struct
  {
    mp_limb_t _bph; /* high word of mod (2^(bits_per_limb)+1) it has
    		     the value 0 or 1 */
    mp_limb_t _bpl; /* low word of mod (2^(bits_per_limb)+1) */
    mp_limb_t _bm;  /* mod (2^(bits_per_limb)-1) */
    mp_limb_t _sp;  /* mod (2^(bits_per_double)+1) */
    mp_limb_t _sm;  /* mod (2^(bits_per_double)-1) */
    mp_limb_t _v;   /* void, to alignment */
  }
mod_test;
#endif

#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)
# include "gthreads.h"
#endif


/************  A lot of global definitions *********************
Most of them are set when calling to yea_init() procedure: 
Y_TWDF: the array of pointers to twidle trig. factors for forward transform. 
	Example, TWDF[i] points to twidle factors for pass i+1 on 
	forward transform.
Y_TWDB: idem to Y_TWDF for backward transform.
Y_DYADIC: Trig. factors used in dyadic mul. 
Y_PF, Y_PB:  The pointers used to allocate trig. factors 
Y_D:    The pointer used to allocate dyadic trig. factors.
F_x_yr:   real(G^(-x/y)). The real part for radix x forward.
F_x_yi:   The imaginary part of powers for forward
B_x_yr:   real(G^(x/y)). The real part for backward transform.
B_x_yi:   The imaginary part for backward transform.
 
Y_LENGTH the FFT length (in size_of_complex units). 
Y_Q, Y_K . Y_LENGTH= TP_Q*2^(TP_K)  Here, TP_Q=5,6,7 or 8 
Y_NRADICES number of radix reduction passes in the FFT.
Y_PLAN the plan for radix reduction. Y_PLAN[i] is the i+1 radix
	reduction on forward transform and Y_NRADICES-1-i on backward.
Y_LLEFT[i] the product of Y_PLAN elements less or equal than i.  
Y_LRIGHT[i] the product of Y_PLAN elelments great or equal than i.
**************************************************************************/
/* Pseudo-complex version */
extern y_limb_t F_1_4r,F_1_4i,F_1_8r,F_1_8i,F_3_8r,F_3_8i,F_1_16r,F_1_16i,
  F_3_16r,F_3_16i,F_5_16r,F_5_16i,F_7_16r,F_7_16i;
extern y_limb_t  F_1_32r,F_1_32i,F_3_32r,F_3_32i,F_5_32r,F_5_32i,F_7_32r,
  F_7_32i,F_9_32r,F_9_32i,F_11_32r,F_11_32i,F_13_32r,F_13_32i,F_15_32r,
  F_15_32i;
extern y_limb_t F_1_3r,F_1_3i,FM150,F_1_6r,F_1_6i,F_1_12r,F_1_12i,F_5_12r,F_5_12i;
extern y_limb_t F_1_9r,F_1_9i,F_2_9r,F_2_9i,F_4_9r,F_4_9i;
extern y_limb_t F_1_5r,F_1_5i,F_2_5r,F_2_5i;
extern y_limb_t FN1_5r, FN1_5i, FNM125, FN2_5r, FN2_5i;
extern y_limb_t F_1_10r,F_1_10i,F_3_10r,F_3_10i;
extern y_limb_t F_1_7r,F_1_7i,F_2_7r,F_2_7i,F_3_7r,F_3_7i;
extern y_limb_t FN1_7r,FN1_7i,FN2_7r,FN2_7i,FN3_7r,FN3_7i,FN4_7r,FN4_7i;
extern y_limb_t F_1_14r,F_1_14i,F_3_14r,F_3_14i,F_5_14r,F_5_14i;
extern y_limb_t B_1_4r,B_1_4i,B_1_8r,B_1_8i,B_3_8r,B_3_8i,B_1_16r,B_1_16i,
  B_3_16r,B_3_16i,B_5_16r,B_5_16i,B_7_16r,B_7_16i;
extern y_limb_t  B_1_32r,B_1_32i,B_3_32r,B_3_32i,B_5_32r,B_5_32i,B_7_32r,
  B_7_32i,B_9_32r,B_9_32i,B_11_32r,B_11_32i,B_13_32r,B_13_32i,B_15_32r,
  B_15_32i;
extern y_limb_t B_1_3r,B_1_3i,BM150,B_1_6r,B_1_6i,B_1_12r,B_1_12i,B_5_12r,B_5_12i;
extern y_limb_t B_1_9r,B_1_9i,B_2_9r,B_2_9i,B_4_9r,B_4_9i;
extern y_limb_t B_1_5r,B_1_5i,B_2_5r,B_2_5i;
extern y_limb_t BN1_5r, BN1_5i, BNM125, BN2_5r, BN2_5i;
extern y_limb_t B_1_10r,B_1_10i,B_3_10r,B_3_10i;
extern y_limb_t B_1_7r,B_1_7i,B_2_7r,B_2_7i,B_3_7r,B_3_7i;
extern y_limb_t BN1_7r,BN1_7i,BN2_7r,BN2_7i,BN3_7r,BN3_7i,BN4_7r,BN4_7i;
extern y_limb_t B_1_14r,B_1_14i,B_3_14r,B_3_14i,B_5_14r,B_5_14i;
extern y_size_t *Y_PLAN,*Y_LRIGHT,*Y_LLEFT;
extern y_ptr *Y_TWDF,*Y_TWDB,Y_DYADIC,Y_WORK;
extern y_size_t Y_LENGTH,Y_Q,Y_K,Y_NRADICES,Y_OLDRADICES;
extern y_size_t Y_POWERS[];
extern int Y_FFTERR;
extern int Y_PASS2;
#if defined(__GNUC__)
extern int BAD_STACK;
#endif

#if defined(Y_NEW_PADDING)
extern  y_size_t Y_SHIFT0, Y_SHIFT1, Y_SHIFT2, Y_SHIFT3, Y_SHIFT4;
#endif
/************  A lot of global definitions *********************
 
**************** SSE2 versions *****************************************/
#if defined(Y_USE_SSE2)
extern Y__M128D MM_1_0, MM_0_5, MM_Y0LO, MM_Y0HI, MM_MTWO;

# if !defined(Y_SIMULATE_SSE2)
/*
    an xor to MM_YCHS will change the sign 
    an and to MM_YABS will compute the abs value
    an and to MM_Y0LO set low part to zero
    an and to MM_Y0HI set HI part to zero
*/
extern Y__M128D MM_YCHS, MM_YABS;
# endif

extern Y__M128D MM_F_1_4r, MM_F_1_4i, MM_F_1_8r, MM_F_1_8i, MM_F_3_8r,
  MM_F_3_8i, MM_F_1_16r, MM_F_1_16i, MM_F_3_16r, MM_F_3_16i, MM_F_5_16r,
  MM_F_5_16i, MM_F_7_16r, MM_F_7_16i;

extern Y__M128D MM_B_1_4r, MM_B_1_4i, MM_B_1_8r, MM_B_1_8i, MM_B_3_8r,
  MM_B_3_8i, MM_B_1_16r, MM_B_1_16i, MM_B_3_16r, MM_B_3_16i, MM_B_5_16r,
  MM_B_5_16i, MM_B_7_16r, MM_B_7_16i;

extern Y__M128D MM_F_1_3r, MM_F_1_3i, MM_FM150, MM_F_1_6r, MM_F_1_6i,
  MM_F_1_12r, MM_F_1_12i, MM_F_5_12r, MM_F_5_12i;

extern Y__M128D MM_F_1_9r, MM_F_1_9i, MM_F_2_9r, MM_F_2_9i, MM_F_4_9r,
  MM_F_4_9i;

extern Y__M128D MM_B_1_3r, MM_B_1_3i, MM_BM150, MM_B_1_6r, MM_B_1_6i,
  MM_B_1_12r, MM_B_1_12i, MM_B_5_12r, MM_B_5_12i;

extern Y__M128D MM_B_1_9r, MM_B_1_9i, MM_B_2_9r, MM_B_2_9i, MM_B_4_9r,
  MM_B_4_9i;

extern Y__M128D MM_FNM125, MM_BNM125;

extern Y__M128D MM_F_1_5r, MM_F_1_5i, MM_F_2_5r, MM_F_2_5i;

extern Y__M128D MM_FN1_5r, MM_FN1_5i, MM_FN2_5r, MM_FN2_5i;

extern Y__M128D MM_F_1_10r, MM_F_1_10i, MM_F_3_10r, MM_F_3_10i;

extern Y__M128D MM_B_1_5r, MM_B_1_5i, MM_B_2_5r, MM_B_2_5i;

extern Y__M128D MM_BN1_5r, MM_BN1_5i, MM_BN2_5r, MM_BN2_5i;

extern Y__M128D MM_B_1_10r, MM_B_1_10i, MM_B_3_10r, MM_B_3_10i;

extern Y__M128D MM_F_1_7r, MM_F_1_7i, MM_F_2_7r, MM_F_2_7i, MM_F_3_7r,
  MM_F_3_7i;

extern Y__M128D MM_FN1_7r, MM_FN1_7i, MM_FN2_7r, MM_FN2_7i, MM_FN3_7r,
  MM_FN3_7i, MM_FN4_7r, MM_FN4_7i;

extern Y__M128D MM_F_1_14r, MM_F_1_14i, MM_F_3_14r, MM_F_3_14i,
  MM_F_5_14r, MM_F_5_14i;

extern Y__M128D MM_B_1_7r, MM_B_1_7i, MM_B_2_7r, MM_B_2_7i, MM_B_3_7r,
  MM_B_3_7i;

extern Y__M128D MM_BN1_7r, MM_BN1_7i, MM_BN2_7r, MM_BN2_7i, MM_BN3_7r,
  MM_BN3_7i, MM_BN4_7r, MM_BN4_7i;

extern Y__M128D MM_B_1_14r, MM_B_1_14i, MM_B_3_14r, MM_B_3_14i,
  MM_B_5_14r, MM_B_5_14i;
#endif

/**************** and now,  the routines ********************************/
void y_root(y_ptr ,y_ptr , y_size_t, y_size_t, int );

void y_initradix(void);

void y_initnumfactg(void);

y_size_t y_init( y_size_t);

void Mr_Proper( void );

y_limb_t tricky(void);

#ifdef FFT_CHECK
mp_limb_t y_ninemod(mp_ptr, y_size_t);

int y_fntcheck(mp_limb_t ,mp_limb_t , mp_limb_t);

int y_fntcheck2(mp_limb_t ,mp_limb_t , mp_limb_t);

int y_fntcheck3(mod_test *, mod_test *, mod_test *);

int y_fntcheck4(mp_limb_t ,mp_limb_t , mp_limb_t);

int y_mod_tests(mod_test *, mod_test *, mod_test *);
#endif

int y_check_threshold_pass2(void);

y_ptr get_fftworkspace(y_size_t, y_ptr *);

void y_limb_to_workspace( y_ptr, mp_ptr, y_size_t, y_size_t);

void y_norm_and_carry(mp_ptr ,y_ptr, y_size_t ,y_size_t);

#ifdef FFT_CHECK
void y_limb_to_workspace_four_tests( y_ptr, mp_ptr, y_size_t, y_size_t,
                                     mod_test *);

void y_norm_and_carry_four_tests(mp_ptr ,y_ptr, y_size_t ,y_size_t,mod_test *);
#endif
void y_fftf( y_ptr ,y_size_t);

void y_fftb( y_ptr ,y_size_t);

void y_convfft(y_ptr ,y_ptr, y_size_t);

void y_squarfft(y_ptr ,y_size_t);

void y_conv(y_ptr, y_ptr);

void y_squar(y_ptr);

mp_limb_t y_mul(mp_ptr, mp_ptr, y_size_t, mp_ptr, y_size_t);

int y_mpn_mul_fft_full(mp_ptr, mp_ptr, y_size_t, mp_ptr, y_size_t);

int y_mpn_mul_fft_full_sqr(mp_ptr, mp_ptr, y_size_t);

void radix_2_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_2_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_2_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_2_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_4_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_4_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_4_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_4_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_8_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_8_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_8_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_8_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_16_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_16_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_16_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_16_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_32_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmp_32_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_12_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_12_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_12_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_12_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_9_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_9_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_9_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_9_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_9_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_9_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_9_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_9_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_9_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_9_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_9_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_9_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_10_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_10_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_10_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_10_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_14_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_14_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_14_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_14_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_5_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_5_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_5_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_5_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_5_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_5_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_5_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_5_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_5_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_5_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_5_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_5_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_6_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_6_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_6_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_6_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_6_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_6_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_6_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_6_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_6_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_6_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_6_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_6_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_7_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_7_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_7_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_7_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_7_dif_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_7_dif(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_7_dit_notw(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_7_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm0_7_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm1_7_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm2_7_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm3_7_dit(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm_4_void( void);

void radixmm_5_void( void);

void radixmm_6_void( void);

void radixmm_7_void( void);

void radixmm_8_void( void);

void radixmm_9_void( void);

void radixmm_16_void( void);

void radix_4_square_void( void );

void radix_8_square_void( void );

/* routines for yeafft1.c */

void y_convolution(y_ptr, y_ptr, y_size_t);

void y_auto_convolution(y_ptr, y_size_t);

int y_fftf_pass1(y_ptr , y_size_t);

void y_fftb_pass1(y_ptr , y_size_t, int);

void y_fftf_pass2(y_ptr , y_size_t , y_size_t , int );

void y_fftb_pass2(y_ptr , y_size_t , y_size_t , int );

void y_fftf_mul_fftb_pass2(y_ptr , y_ptr , int );

void y_fftf_squar_fftb_pass2(y_ptr , int );

void y_fftf_mul_fftb(y_ptr , y_ptr , y_size_t);

void y_fftf_mul_fftb_block(y_ptr , y_ptr , y_size_t , y_size_t, y_size_t);

void y_fftf_square_fftb(y_ptr , y_size_t);

void y_fftf_square_fftb_block(y_ptr,  y_size_t , y_size_t, y_size_t);

void radix_2_dif_mul_dit(y_ptr, y_ptr, y_size_t);

void radix_2_dif_square_dit(y_ptr, y_size_t);

void radix_2_dif_mul_dit_block(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_2_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dif_mul_dit(y_ptr, y_ptr, y_size_t);

void radix_4_dif_square_dit(y_ptr, y_size_t);

void radixmm_4_dif_square_dit(y_ptr, y_size_t);

void radix_4_dif_mul_dit_block(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_4_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm_4_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dif_mul_dit(y_ptr, y_ptr, y_size_t);

void radix_8_dif_square_dit(y_ptr, y_size_t);

void radixmm_8_dif_square_dit(y_ptr, y_size_t);

void radix_8_dif_mul_dit_block(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_8_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radixmm_8_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dif_mul_dit(y_ptr, y_ptr, y_size_t);

void radix_16_dif_square_dit(y_ptr, y_size_t);

void radix_16_dif_mul_dit_block(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_16_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dif_mul_dit(y_ptr, y_ptr, y_size_t);

void radix_32_dif_square_dit(y_ptr, y_size_t);

void radix_32_dif_mul_dit_block(y_ptr, y_ptr, y_size_t, y_size_t, y_size_t);

void radix_32_dif_square_dit_block(y_ptr, y_size_t, y_size_t, y_size_t);

void fget_twiddle_factors_32_last(struct y_complex *, y_size_t , y_size_t );

y_ptr fradix_32_twd_last(struct y_complex * , y_limb_t  * , y_size_t , struct y_complex  * );

void fradix_32_notwd_first(struct y_complex * ,y_ptr );

/*$Id$*/







