/* $Id$ */
/***
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2002-2006  Guillermo Ballester Valor
 
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
***/

/*
   This is the include file to use x-86 SSE2 extensions in YEAFFT-Glucas 
   routines 
   
   ABOUT HOW TO USE SSE2 extensions:
 
   the sse2-code tries to use as amount of normal code as possible. To 
   achieve this, it tries to parallelize the operations made in single
   mode to sse2 form.
 
   Memory access:
 
   In the normal code the elements in memory are stored:
 
   position in memory          element
   p                           real part of i complex
   p + 1                       imag part of i complex
   p + 2                       real part of i+1 complex
   p + 3                       imag part of i+1 complex
 
   As in the normal code the load/store operations are paired assuming real and
   imag part of a complex, in the sse2 code also are load/stored two __md128
   registers.
 
   When the sse2 code loads from normal (legacy) stored memory, it loads using
   '_INTER_' versioned macros. When stores to normal form it then also stores 
   with this versioned macros.
 
   In load/store macros we can include one or two pointers to the double array.
   If one, 'p' it refers to the lowest pointer memory and it assumes the next 
   complex to deal with is the next 'p + 2'. This macros are 'NEXT' versioned
   if two, every pointer links to real part of a complex. Anyway, the two 
   loaded complex in 'INTER' macros are interleaved in sse2 types as follows
 
   __md128 a, b;
   double *p1,*p2;
   a = (*p1, *p2)
   b = (*(p1 + 1),*(p2 + 1))
   
   When stores in _INTER_ format:
  
   *(p1) = a.d[0]; *(p1 + 1) = b.d[0];
   *(p2) = a.d[1]; *(p2 + 1) = b.d[1];
   
   If we want save some interleaved process, we can store/load in a new
   sse2 format
 
   position in memory          element
   p                           real part of i complex
   p + 1                       real part of i+1 complex
   p + 2                       imag part of i complex
   p + 3                       imag part of i+1 complex
   
   so, the loads and store from/to memory are easier and faster. This uses
   versioned macros without INTER.
 
   When loading
 
    __md128 a, b;
   double *p1,*p2;
   a = (*p1, *(p1 + 1))
   b = (*p2, *(p2 + 1))
  
   when storing
 
   *(p1) = a.d[0]; *(p1 + 1) = a.d[1];
   *(p2) = b.d[0]; *(p2 + 1) = b.d[1];
   
*/

#if !(defined(Y_PENTIUM4) || defined(Y_AMD64)) || defined(__GNUC__) && !((( __GNUC__ == 3) && ( __GNUC_MINOR__ >= 3 )) || ( __GNUC__ > 3)) && !defined(__ICC)
# define Y_SIMULATE_SSE2
#else
# if !defined(__GNUC__) && (!defined(__ICC) || !(defined(Y_PENTIUM4) || defined(Y_AMD64)))
#  define Y_SIMULATE_SSE2
# endif
#endif

#if defined(Y_SIMULATE_SSE2)
/*
   this macros are only to test correctness, no speed, of the code
   unfortunately I still have no access to a machine with real SSE2 
   extensions, so I have to simulate it.
 
   Updated on Apr-2003: 
   We already have pentium4 with SSE2 available, but still is 
   interesting an emulation.
*/
/*
#if (( __GNUC__ == 3) && ( __GNUC_MINOR__ >= 3 )) || ( __GNUC__ > 3)
 
typedef int Y__M128D __atribute__ ((mode(V2DF)));
*/

typedef struct
  {
    double Y_ALIGNED(16) d[2];
  }
Y__M128D __attribute__ ((aligned(16)));

# define Y_MM_LOAD_PD( _res, __p) \
   _res.d[0] = *(__p) ;           \
   _res.d[1] = *(__p + 1);

# define Y_MM_LOAD1_PD( _res, __p)\
   _res.d[0] = *(__p);            \
   _res.d[1] = *(__p);

# define Y_MM_STORE_PD(__p, _op)  \
   *(__p)  = _op.d[0];            \
   *(__p + 1) = _op.d[1];

# define Y_MM_SET_SD( _res, _op) \
  _res.d[0] = _op;               \
  _res.d[1] = 0.0;

# define Y_MM_SET_PD( _res, _op0, _op1) \
  _res.d[0] = _op0;        \
  _res.d[1] = _op1;

# define Y_MM_SET1_PD( _res, _op) \
  _res.d[0] = _op;                \
  _res.d[1] = _op;

# define Y_MM_SETZERO_PD( _res) \
  _res.d[0] = 0.0;              \
  _res.d[1] = 0.0;

# define Y_MM_MOV_PD( _res, _op) \
   _res.d[0] = _op.d[0];         \
   _res.d[1] = _op.d[1];

# define Y_MM_UNPACKHI_PD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[1];                    \
   _res.d[1] = _op2.d[1];

# define Y_MM_UNPACKLO_PD( _res, _op1, _op2) \
   _res.d[1] = _op2.d[0];                    \
   _res.d[0] = _op1.d[0];

# define Y_MM_ADD_PD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] + _op2.d[0];   \
   _res.d[1] = _op1.d[1] + _op2.d[1];

# define Y_MM_ADD_SD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] + _op2.d[0];

# define Y_MM_SUB_PD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] - _op2.d[0];   \
   _res.d[1] = _op1.d[1] - _op2.d[1];

# define Y_MM_SUB_SD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] - _op2.d[0];

# define Y_MM_MUL_PD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] * _op2.d[0];   \
   _res.d[1] = _op1.d[1] * _op2.d[1];

# define Y_MM_MUL_SD( _res, _op1, _op2) \
   _res.d[0] = _op1.d[0] * _op2.d[0];

# define Y_MM_SWAP_PD( _res) \
   {                         \
     BIG_DOUBLE _aux;        \
     _aux = _res.d[0];       \
     _res.d[0] = _res.d[1];  \
     _res.d[1] = _aux;       \
   }

# define Y_MM_ABS_PD( _res, _op1)\
   _res.d[0] = fabs( _op1.d[0] );\
   _res.d[1] = fabs( _op1.d[1] );

# define Y_MM_CHS_PD( _res, _op1)\
   _res.d[0] = -(_op1.d[0]);\
   _res.d[1] = -(_op1.d[1]);

# define Y_MM_MAX_PD( _res, _op1, _op2)                          \
   _res.d[0] = (_op1.d[0] >= _op2.d[0]) ? _op1.d[0] : _op2.d[0]; \
   _res.d[1] = (_op1.d[1] >= _op2.d[1]) ? _op1.d[1] : _op2.d[1];

# define Y_MM_MAX_SD( _res, _op1, _op2)                          \
   _res.d[0] = (_op1.d[0] >= _op2.d[0]) ? _op1.d[0] : _op2.d[0];

/* The following emulated macros have awful performance */
# define Y_MM_AND_PD( _res, _op1, _op2)                        \
  {                                                            \
     char *paux, *paux1, *paux2;                               \
     size_t ix;                                                \
     paux = (char *)(&_res.d[0]);                              \
     paux1 = (char *)(&_op1.d[0]);                             \
     paux2 = (char *)(&_op2.d[0]);                             \
     for( ix = 0; ix < sizeof(double) ; ix ++)                 \
        paux[ix] = paux1[ix] & paux2[ix];                      \
     paux = (char *)(&_res.d[1]);                              \
     paux1 = (char *)(&_op1.d[1]);                             \
     paux2 = (char *)(&_op2.d[1]);                             \
     for( ix = 0; ix < sizeof(double) ; ix ++)                 \
        paux[ix] = paux1[ix] & paux2[ix];                      \
  }



# define Y_MM_CMPGT_PD(_res, _op1, _op2)           \
{                                                  \
  if( (_op1.d[0] > _op2.d[0]) )                    \
    {                                              \
       char *paux;                                 \
       size_t ix;                                  \
       paux = (char *)(&_res.d[0]);                \
       for( ix = 0; ix < sizeof(double) ; ix ++)   \
          paux[ix] = 0xff;                         \
    }                                              \
  else _res.d[0] = 0.0;                            \
  if( (_op1.d[1] > _op2.d[1]) )                    \
    {                                              \
       char *paux;                                 \
       size_t ix;                                  \
       paux = (char *)(&_res.d[1]);                \
       for( ix = 0; ix < sizeof(double) ; ix++)    \
          paux[ix] = 0xff;                         \
    }                                              \
  else _res.d[1] = 0.0;                            \
}

# define Y_MM_CMPGT_SD(_res, _op1, _op2)           \
{                                                  \
  if( (_op1.d[0] > _op2.d[0]) )                    \
    {                                              \
       char *paux;                                 \
       size_t ix;                                  \
       paux = (char *)(&_res.d[0]);                \
       for( ix = 0; ix < sizeof(double) ; ix ++)   \
          paux[ix] = 0xff;                         \
    }                                              \
  else _res.d[0] = 0.0;                            \
}


# define Y_MM_MOVEMASK_PD( _res, _op)\
  _res = 0;                          \
  if (_op.d[0] != 0.0) _res += 1;    \
  if (_op.d[1] != 0.0) _res += 2;


#else

# if defined(__ICC)
#  include <emmintrin.h>
# elif defined(__GNUC__)
#  include <xmmintrin.h>
# endif

typedef __m128d Y__M128D /*__attribute__ ((aligned (16)))*/;

# define Y_MM_LOAD_PD( _res, _p) _res = _mm_load_pd(_p)

# define Y_MM_LOAD1_PD( _res, _p) _res = _mm_load1_pd(_p)

# define Y_MM_SET_SD( _res, _op) _res = _mm_set_sd( (BIG_DOUBLE) _op)

# define Y_MM_SET_PD( _res, _op0, _op1) _res = _mm_setr_pd( (BIG_DOUBLE) _op0, (BIG_DOUBLE) _op1)

# define Y_MM_SET1_PD( _res, _op) _res = _mm_set1_pd( (BIG_DOUBLE)_op )

# define Y_MM_SETZERO_PD( _res) _res = _mm_setzero_pd(  )

# define Y_MM_STORE_PD( _p, _op) _mm_store_pd( _p, _op)

# define Y_MM_UNPACKHI_PD( _res, _op1, _op2) _res = _mm_unpackhi_pd( _op1, _op2)

# define Y_MM_UNPACKLO_PD( _res, _op1, _op2) _res = _mm_unpacklo_pd( _op1, _op2)

# define Y_MM_ADD_PD( _res, _op1, _op2) _res = _mm_add_pd( _op1, _op2)

# define Y_MM_ADD_SD( _res, _op1, _op2) _res = _mm_add_sd( _op1, _op2)

# define Y_MM_SUB_PD( _res, _op1, _op2) _res = _mm_sub_pd( _op1, _op2)

# define Y_MM_SUB_SD( _res, _op1, _op2) _res = _mm_sub_sd( _op1, _op2)

# define Y_MM_MUL_PD( _res, _op1, _op2) _res = _mm_mul_pd( _op1, _op2)

# define Y_MM_MUL_SD( _res, _op1, _op2) _res = _mm_mul_sd( _op1, _op2)

# define Y_MM_MOV_PD( _res, _op1) _res = _op1

# define Y_MM_SWAP_PD( _res ) _res = _mm_shuffle_pd( _res, _res, 0x01)

# define Y_MM_AND_PD( _res, _op1, _op2) _res = _mm_and_pd( _op1, _op2)

# define Y_MM_ANDNOT_PD( _res, _op1, _op2) _res = _mm_andnot_pd( _op1, _op2)

# define Y_MM_OR_PD( _res, _op1, _op2) _res = _mm_or_pd( _op1, _op2)

# define Y_MM_XOR_PD( _res, _op1, _op2) _res = _mm_xor_pd( _op1, _op2)

# define Y_MM_ABS_PD( _res, _op1) _res = _mm_and_pd( _op1, MM_YABS)

# define Y_MM_CHS_PD( _res, _op1) _res = _mm_xor_pd( _op1, MM_YCHS)

# define Y_MM_CHS_SD( _res, _op1) _res = _mm_xor_sd( _op1, MM_YCHS)

# define Y_MM_CMPGT_PD( _res, _op1, _op2) _res = _mm_cmpgt_pd( _op1, _op2)

# define Y_MM_CMPGT_SD( _res, _op1, _op2) _res = _mm_cmpgt_sd( _op1, _op2)

# define Y_MM_MAX_PD( _res, _op1, _op2) _res = _mm_max_pd( _op1, _op2)

# define Y_MM_MAX_SD( _res, _op1, _op2) _res = _mm_max_sd( _op1, _op2)

# define Y_MM_MOVEMASK_PD( _i, _op) _i = _mm_movemask_pd( _op )

#endif

typedef Y__M128D * ysse2_ptr;


/*
   This macro loads two complex in main array to two Y__MD128 
   _dest0 and _dest1. _dest0 will contain the real part and _dest1 the 
   imaginary ones
 
   _pd0 and _pd1 are (double *) and are pointing to real part (lower memory)
   to the respective complex in main array
*/


#define Y_MM_LOAD_INTER_PAIR( _dest0, _dest1, __pd0, __pd1) \
   {                                                        \
      Y__M128D __aux0,__aux1;                               \
      Y_MM_LOAD_PD ( __aux0, __pd0);                        \
      Y_MM_LOAD_PD ( __aux1, __pd1);                        \
      Y_MM_UNPACKLO_PD( _dest0, __aux0, __aux1);            \
      Y_MM_UNPACKHI_PD( _dest1, __aux0, __aux1);            \
   }


#define Y_MM_MOV_INTER_PAIR( _dest0, _dest1)     \
   {                                             \
      Y__M128D __aux0,__aux1;                    \
      Y_MM_MOV_PD ( __aux0, _dest0);             \
      Y_MM_MOV_PD ( __aux1, _dest1);             \
      Y_MM_UNPACKLO_PD( _dest0, __aux0, __aux1); \
      Y_MM_UNPACKHI_PD( _dest1, __aux0, __aux1); \
   }

/*
#define Y_MM_LOAD_INTER_PAIR( _dest0, _dest1, __pd0, __pd1)  \
   Y_MM_LOAD_PD ( _dest0, __pd0);                            \
   Y_MM_LOAD_PD ( _dest1, __pd0);                            \
   Y_MM_UNPACKLO_PD( _dest0, _dest0, *((Y__M128D *)(__pd1)));\
   Y_MM_UNPACKHI_PD( _dest1, _dest1, *((Y__M128D *)(__pd1)));
*/

#define Y_MM_LOAD_INTER_TWO_PAIR( _dest0, _dest1, _dest2, _dest3, __pd0, __pd1, __pd2, __pd3) \
   {                                                             \
      Y__M128D __aux0, __aux1;                                   \
      Y_MM_LOAD_PD ( _dest0, __pd0);               /* #1-6   */  \
      Y_MM_LOAD_PD ( _dest2, __pd2);               /* #2-7   */  \
      Y_MM_LOAD_PD ( __aux0, __pd1);               /* #3-8   */  \
      Y_MM_LOAD_PD ( __aux1, __pd3);               /* #4-9   */  \
      Y_MM_MOV_PD( _dest1, _dest0);      /*STALL*/ /* #7-12  */  \
      Y_MM_MOV_PD( _dest3, _dest2);                /* #8-13  */  \
      Y_MM_UNPACKLO_PD( _dest0, _dest0, __aux0);   /* #9-14  */  \
      Y_MM_UNPACKLO_PD( _dest2, _dest2, __aux1);   /* #10-15 */  \
      Y_MM_UNPACKHI_PD( _dest1, _dest1, __aux0);   /* #13-16 */  \
      Y_MM_UNPACKHI_PD( _dest3, _dest3, __aux1);   /* #14-19 */  \
   }

#define Y_MM_LOAD_INTER_TWO_PAIR_NEXT( _dest0, _dest1, _dest2, _dest3, __pd0, __pd2) \
   {                                                             \
      Y__M128D __aux0, __aux1;                                   \
      Y_MM_LOAD_PD ( _dest0, __pd0);               /* #1-6   */  \
      Y_MM_LOAD_PD ( __aux0, __pd0 + 2);           /* #2-7   */  \
      Y_MM_LOAD_PD ( _dest2, __pd2);               /* #3-8   */  \
      Y_MM_LOAD_PD ( __aux1, __pd2 + 2);           /* #4-9   */  \
      Y_MM_MOV_PD( _dest1, _dest0);                /* #7-12  */  \
      Y_MM_MOV_PD( _dest3, _dest2);                /* #8-13  */  \
      Y_MM_UNPACKLO_PD( _dest0, _dest0, __aux0);   /* #9-14  */  \
      Y_MM_UNPACKLO_PD( _dest2, _dest2, __aux1);   /* #10-15 */  \
      Y_MM_UNPACKHI_PD( _dest1, _dest1, __aux0);   /* #13-16 */  \
      Y_MM_UNPACKHI_PD( _dest3, _dest3, __aux1);   /* #14-19 */  \
   }

/* the difference with the above macro is here there is no interleave step */
#define Y_MM_LOAD_PAIR( _dest0, _dest1, __pd0, __pd1) \
   Y_MM_LOAD_PD ( _dest0, __pd0);                     \
   Y_MM_LOAD_PD ( _dest1, __pd1);

/* This macro is to interchange format */
#define Y_MM_INTER_PAIR( _op0, _op1)           \
   {                                           \
      Y__M128D __aux0,__aux1;                  \
      Y_MM_MOV_PD( __aux0, _op0);              \
      Y_MM_MOV_PD( __aux1, _op1);              \
      Y_MM_UNPACKLO_PD( _op0, __aux0, __aux1); \
      Y_MM_UNPACKHI_PD( _op1, __aux0, __aux1); \
   }


/*
   This load assuming first is the unity. This macro
   is used when loading trig factors, when first is
   always the complex unit (1.0,0.0)
*/
#define Y_MM_LOAD_INTER_ONE_HALF( _dest0, _dest1, _pd0)\
   {                                                  \
      Y__M128D __aux0, __aux1;                        \
      Y_MM_LOAD_PD ( __aux1, _pd0);                   \
      Y_MM_SET_SD ( __aux0, 1.0);                     \
      Y_MM_UNPACKLO_PD( _dest0, __aux0, __aux1);      \
      Y_MM_UNPACKHI_PD( _dest1, __aux0, __aux1);      \
   }

/*
  This is used to load trig factors in forward phase 
 */
#define Y_MM_LOAD_INTER_DOUBLED( _dest0, _dest1, __pd0)   \
{                                                         \
    Y_MM_LOAD1_PD( _dest0, __pd0);                        \
    Y_MM_LOAD1_PD( _dest1, __pd0 + 1);                    \
}

/* this version assumes the _pd1 = _pd0 + 2  */

#define Y_MM_LOAD_INTER_PAIR_NEXT( _dest0, _dest1, __pd)  \
   {                                                      \
      Y__M128D __aux0,__aux1;                             \
      Y_MM_LOAD_PD ( __aux0, __pd);                       \
      Y_MM_LOAD_PD ( __aux1, (__pd + 2));                 \
      Y_MM_UNPACKLO_PD( _dest0, __aux0, __aux1);          \
      Y_MM_UNPACKHI_PD( _dest1, __aux0, __aux1);          \
   }

/*
#define Y_MM_LOAD_INTER_PAIR_NEXT( _dest0, _dest1, __pd)        \
   Y_MM_LOAD_PD ( _dest0, __pd);                                \
   Y_MM_LOAD_PD ( _dest1, __pd);                                \
   Y_MM_UNPACKLO_PD( _dest0, _dest0, *((Y__M128D *)(__pd + 2)));\
   Y_MM_UNPACKHI_PD( _dest1, _dest1, *((Y__M128D *)(__pd + 2)));
*/

#define Y_MM_LOAD_PAIR_NEXT( _dest0, _dest1, __pd) \
   Y_MM_LOAD_PD ( _dest0, __pd );                  \
   Y_MM_LOAD_PD ( _dest1, (__pd + 2));


/*
   This macro is the reverse of Y_LOAD_INTER_PAIR. 
   it stores two complex in main array to pointers 
   _pd0 and _pd1, extracting real and imaginary parts
   from _sour0 and _sour1. _sour0 it supossed to include 
   the real parts, and _sour1 the imaginary parts.
*/
#define Y_MM_STORE_INTER_PAIR( __pd0, __pd1, _sour0, _sour1) \
   {                                                         \
      Y__M128D __aux0,__aux1;                                \
      Y_MM_UNPACKLO_PD( __aux0, _sour0, _sour1);             \
      Y_MM_UNPACKHI_PD( __aux1, _sour0, _sour1);             \
      Y_MM_STORE_PD( __pd0, __aux0);                         \
      Y_MM_STORE_PD( __pd1, __aux1);                         \
   }

/* the difference with the above macro is here there is no interleave step */
#define Y_MM_STORE_PAIR( __pd0, __pd1, _sour0, _sour1) \
   Y_MM_STORE_PD( __pd0, _sour0);                      \
   Y_MM_STORE_PD( __pd1, _sour1)

#define Y_MM_STORE_INTER_PAIR_NEXT(__pd, _sour0, _sour1)   \
   {                                                       \
      Y__M128D __aux0,__aux1;                              \
      Y_MM_UNPACKLO_PD( __aux0, _sour0, _sour1);           \
      Y_MM_UNPACKHI_PD( __aux1, _sour0, _sour1);           \
      Y_MM_STORE_PD( __pd , __aux0);                       \
      Y_MM_STORE_PD((__pd + 2), __aux1);                   \
   }

#define Y_MM_STORE_INTER_LO(__pd, _sour0, _sour1)   \
   {                                                \
      Y__M128D __aux0;                              \
      Y_MM_UNPACKLO_PD( __aux0, _sour0, _sour1);    \
      Y_MM_STORE_PD( __pd , __aux0);                \
   }

#define Y_MM_STORE_INTER_HI(__pd, _sour0, _sour1)   \
   {                                                \
      Y__M128D __aux0;                              \
      Y_MM_UNPACKHI_PD( __aux0, _sour0, _sour1);    \
      Y_MM_STORE_PD( __pd , __aux0);                \
   }



#define Y_MM_STORE_PAIR_NEXT( __pd,_sour0, _sour1)    \
   Y_MM_STORE_PD( __pd, _sour0);                      \
   Y_MM_STORE_PD((__pd + 2), _sour1);

#define Y_MM_COMPLEX_MUL( __pd, _a, _b)                                  \
   {                                                                     \
     Y__M128D __aux0, __aux1;                                            \
     Y_MM_MOV_PD ( __aux0, _b); /* aux0 :(br,bi)*/                       \
     Y_MM_SWAP_PD ( __aux0 );   /* aux0 :(bi,br)*/                       \
     Y_MM_MOV_PD ( __aux1, _a); /* aux1 :(ar,ai)*/                       \
     Y_MM_MUL_PD ( __aux0, __aux0, _a); /* aux0 :(ar*bi,ai*br)*/         \
     Y_MM_MUL_PD ( __aux1, __aux1, _b); /* aux1 :(ar*br,ai*bi)*/         \
     Y_MM_MOV_PD ( __pd, __aux0);    /* _pd :(ar*bi,ai*br)  */           \
     Y_MM_SWAP_PD ( __aux0);         /* aux0 :(ai*br,ar*bi) */           \
     Y_MM_ADD_SD ( __pd, __pd, __aux0);/* _pd :(ar*bi+ai*br,ai*br)*/     \
     Y_MM_MOV_PD ( __aux0, __aux1); /* aux0 :(ar*br, ai*bi)*/            \
     Y_MM_SWAP_PD ( __aux1 );      /* aux1 :(ai*bi, ar*br)*/             \
     Y_MM_SUB_SD ( __aux0, __aux0, __aux1); /*aux0 :(ar*br-ai*bi,ai*bi)*/\
     Y_MM_UNPACKLO_PD ( __pd, __aux0, __pd);                             \
   }

#define Y_MM_COMPLEX_DIV( __pd, _sour0, _sour1)                          \
   {                                                                     \
     Y__M128D __aux0, __aux1;                                            \
     Y_MM_MOV_PD ( __aux0, _sour1);  /* aux0 :(br,bi)*/                  \
     Y_MM_SWAP_PD ( __aux0 );        /* aux0 :(bi,br)*/                  \
     Y_MM_MOV_PD ( __aux1, _sour0); /* aux1 :(ar,ai)*/                   \
     Y_MM_MUL_PD ( __aux0, __aux0, _sour0);  /* aux0 :(ar*bi,ai*br)*/    \
     Y_MM_MUL_PD ( __aux1, __aux1, _sour1);  /* aux1 :(ar*br,ai*bi)*/    \
     Y_MM_MOV_PD ( __pd, __aux0);    /* _pd :(ar*bi,ai*br)  */           \
     Y_MM_SWAP_PD ( __pd);           /* _pd :(ai*br,ar*bi) */            \
     Y_MM_SUB_SD ( __pd, __pd, __aux0);/* _pd :(ai*br-ar*bi,ar*bi)*/     \
     Y_MM_MOV_PD ( __aux0, __aux1); /* aux0 :(ar*br, ai*bi)*/            \
     Y_MM_SWAP_PD ( __aux1 );      /* aux1 :(ai*bi, ar*br)*/             \
     Y_MM_ADD_SD ( __aux0, __aux0, __aux1); /*aux0 :(ar*br+ai*bi,ai*bi)*/\
     Y_MM_UNPACKLO_PD ( __pd, __aux0, __pd);                             \
   }

# define Y_MM_PRINTF_PD(_op)                      \
  {                                               \
     y_ptr aux0, aux;                             \
      aux0 = ALLOC_DOUBLES( 2 );                  \
      aux = (double *) ALIGN_DOUBLES( aux0 );     \
      Y_MM_STORE_PD (aux, _op);                   \
      printf("(%lf,%lf)\n",aux[0],aux[1]);        \
      free (aux0);                                \
  }

/* This routine sums both parts of _op1 and _op2 and sums it to _Sum */
#define sse2_sum(_Sum, _op0, _op1)  \
  {                                 \
    Y__M128D _aux0;                 \
    BIG_DOUBLE _xaux[2];            \
    Y_MM_ADD_PD (_aux0, _op0, _op1);\
    Y_MM_STORE_PD (_xaux, _aux0);   \
    _Sum += (_xaux[0] + _xaux[1]);  \
  }

/* This routine sums both parts of _op1 and sums it to _Sum */
#define sse2_sum0(_Sum, _op0)     \
  {                               \
    BIG_DOUBLE _xaux[2];          \
    Y_MM_STORE_PD (_xaux, _op0);  \
    _Sum += (_xaux[0] + _xaux[1]);\
  }
