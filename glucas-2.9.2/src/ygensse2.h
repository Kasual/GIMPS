/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor
 
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
/* This file is the sse2 version of macros in ygeneric */
/*#define OLD_MACROS */
#if !defined(__ICC)
# define ALT_ADDSUB
#endif
#define ASM_SSE2
#if defined(ASM_SSE2) && !defined(__ICC)
/*# define ASM_ADDSUB*/
/*# define ASM_MULMUL */
/*# define ASM_DIVDIV */
/*# define ASM_DIVMUL */
/*# define ASM_MULADDSUB */
/*# define ASM_MULMULADDSUB */
#endif
/* to read some trigonomeric factors */
#define sse2_data_trig_load_inter_2(_t0,_t1,_p)                    \
   Y_MM_LOAD_INTER_PAIR( _t0##r, _t0##i, _p , _p + Y_STEP );       \
   Y_MM_LOAD_INTER_PAIR( _t1##r, _t1##i, _p + 2, _p + 2 + Y_STEP )

#define sse2_data_trig_load_inter_3(_t0, _t1, _t2, _p)             \
   Y_MM_LOAD_INTER_PAIR( _t0##r, _t0##i, _p , _p + Y_STEP );       \
   Y_MM_LOAD_INTER_PAIR( _t1##r, _t1##i, _p + 2, _p + 2 + Y_STEP );\
   Y_MM_LOAD_INTER_PAIR( _t2##r, _t2##i, _p + 4, _p + 4 + Y_STEP )

/* To load from real padded array of data as sse2 */
#define sse2_data_to_local_inter(_t,_array,_k)                          \
        Y_MM_LOAD_INTER_PAIR_NEXT( _t##r, _t##i, _array + addr((_k)<<1))

#define sse2_data_to_local(_t,_array,_k)                          \
        Y_MM_LOAD_PAIR_NEXT( _t##r, _t##i, _array + addr((_k)<<1))

#define sse2_data_to_local_p_inter(_t,_pd,_array,_k)                     \
	_pd = _array + addr(( _k)<<1);                                   \
        Y_MM_LOAD_INTER_PAIR_NEXT( _t##r, _t##i, _array + addr((_k)<<1));\
        prefetch_p(_pd)

#define sse2_data_to_local_p(_t,_pd,_array,_k)                     \
	_pd = _array + addr(( _k)<<1);                             \
        Y_MM_LOAD_PAIR_NEXT( _t##r, _t##i, _array + addr((_k)<<1));\
        prefetch_p(_pd)

#define sse2_data_to_local_pp_inter(_t,_index,_pd)             \
        Y_MM_LOAD_INTER_PAIR_NEXT( _t##r, _t##i, _pd + _index);\
        prefetch_p( _pd + _index)

#define sse2_data_to_local_pp(_t,_index,_pd)             \
        Y_MM_LOAD_PAIR_NEXT( _t##r, _t##i, _pd + _index);\
        prefetch_p( _pd + _index)

/* To store on read padded array from local data as sse2 */

#define sse2_local_to_data_inter(_array,_k,_t)                                \
        Y_MM_STORE_PAIR_INTER_NEXT( (_array + addr(( _k )<<1)) , _t##r, _t##i)

#define sse2_local_to_data(_array,_k,_t)                                \
        Y_MM_STORE_PAIR_NEXT( (_array + addr(( _k )<<1)) , _t##r, _t##i)

#define sse2_local_to_data_p_inter(_pd,_t)              \
        Y_MM_STORE_PAIR_INTER_NEXT( _pd , _t##r, _t##i);\
        prefetch_p(_pd)

#define sse2_local_to_data_p(_pd,_t)              \
        Y_MM_STORE_PAIR_NEXT( _pd , _t##r, _t##i);\
        prefetch_p(_pd)

/* To load from real non-padded array of data as a complex array */
#define sse2_mem_to_local_inter(_t,_array,_k)                         \
        Y_MM_LOAD_PAIR_INTER_NEXT(_t##r, _t##i, (_array + ( _k )<<1) )

#define sse2_mem_to_local(_t,_array,_k)                          \
        Y_MM_LOAD_PAIR_NEXT(_t##r, _t##i, (_array + ( _k )<<1) )

#define sse2_mem_to_local_p_inter(_t,_pd)             \
        Y_MM_LOAD_PAIR_INTER_NEXT(_t##r, _t##i, _pd );\
         prefetch_p (_pd)

#define sse2_mem_to_local_p(_t,_pd)             \
        Y_MM_LOAD_PAIR_NEXT(_t##r, _t##i, _pd );\
         prefetch_p (_pd)

/*********** PSEUDO_COMPLEX_SSE2 GENERAL AUXILIAR MACROS **************/

#define sse2_add(_r, _op0, _op1)       \
  Y_MM_ADD_PD( _r##r, _op0##r, _op1##r);\
  Y_MM_ADD_PD( _r##i, _op0##i, _op1##i);

#define sse2_sub(_r, _op0, _op1)       \
  Y_MM_SUB_PD( _r##r, _op0##r, _op1##r);\
  Y_MM_SUB_PD( _r##i, _op0##i, _op1##i);

/* here __t = __t * __f */
#define sse2_mul0(__t, __f)             \
{                                       \
   Y__M128D __a,__b;                    \
   Y_MM_MUL_PD(__a, __t##r, __f##i);    \
   Y_MM_MUL_PD(__t##r, __t##r, __f##r); \
   Y_MM_MUL_PD(__b, __t##i, __f##i);    \
   Y_MM_MUL_PD(__t##i, __t##i, __f##r); \
   Y_MM_SUB_PD(__t##r, __t##r, __b);    \
   Y_MM_ADD_PD(__t##i, __t##i, __a);    \
}

#define sse2_mul_1_4_F( __t)     \
{                                \
  Y__M128D __a;                  \
  Y_MM_MOV_PD(__a, __t##r);      \
  Y_MM_MOV_PD(__t##r, __t##i);   \
  Y_MM_CHS_PD(__t##i, __a);      \
}

#define sse2_mul_1_4_B( __t)     \
{                                \
  Y__M128D __a;                  \
  Y_MM_MOV_PD(__a, __t##i);      \
  Y_MM_MOV_PD(__t##i, __t##r);   \
  Y_MM_CHS_PD(__t##r, __a);      \
}

#define sse2_mul_1_8_F( __t)              \
{                                         \
  Y__M128D __a;                           \
  Y_MM_MUL_PD( __t##r, __t##r, MM_F_1_8r);\
  Y_MM_MUL_PD( __t##i, __t##i, MM_F_1_8r);\
  Y_MM_MOV_PD( __a, __t##r);              \
  Y_MM_ADD_PD( __t##r, __t##r, __t##i);   \
  Y_MM_SUB_PD( __t##i, __t##i, __a);      \
}

#define sse2_mul_1_8_B( __t)              \
{                                         \
  Y__M128D __a;                           \
  Y_MM_MUL_PD( __t##r, __t##r, MM_F_1_8r);\
  Y_MM_MUL_PD( __t##i, __t##i, MM_F_1_8r);\
  Y_MM_MOV_PD( __a, __t##r);              \
  Y_MM_SUB_PD( __t##r, __t##r, __t##i);   \
  Y_MM_ADD_PD( __t##i, __t##i, __a);      \
}

#define sse2_mul_3_8_F( __t)              \
{                                         \
  Y__M128D __a;                           \
  Y_MM_MUL_PD( __t##r, __t##r, MM_F_1_8i);\
  Y_MM_MUL_PD( __t##i, __t##i, MM_F_1_8i);\
  Y_MM_MOV_PD( __a, __t##r);              \
  Y_MM_SUB_PD( __t##r, __t##r, __t##i);   \
  Y_MM_ADD_PD( __t##i, __t##i, __a);      \
}

#define sse2_mul_3_8_B( __t)              \
{                                         \
  Y__M128D __a;                           \
  Y_MM_MUL_PD( __t##r, __t##r, MM_F_1_8i);\
  Y_MM_MUL_PD( __t##i, __t##i, MM_F_1_8i);\
  Y_MM_MOV_PD( __a, __t##r);              \
  Y_MM_ADD_PD( __t##r, __t##r, __t##i);   \
  Y_MM_SUB_PD( __t##i, __t##i, __a);      \
}


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/

#if defined(ASM_ADDSUB)

# define sse2_asm_addsub(__t0, __t1)       \
__asm__ volatile ("movapd %0, %%xmm6   \n" \
"        movapd %1, %%xmm7             \n" \
"        subpd %2, %%xmm6              \n" \
"        subpd %3, %%xmm7              \n" \
"        addpd %2, %0                  \n" \
"        addpd %3, %1                  \n" \
"        movapd %%xmm6, %2             \n" \
"        movapd %%xmm7, %3             \n" \
: "+&x" (__t0##r), "+&x" (__t0##i), "+&x" (__t1##r), "+&x" (__t1##i)  \
: : "xmm6", "xmm7");

# define sse2_addsub(_t0,_t1) \
  sse2_asm_addsub (_t0, _t1);

# define sse2_load_inter_addsub(_t0, _t1, _array, _k0, _k1)           \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, _array + addr(_k1<<1));  \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));  \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_addsub(_t0, _t1, _array, _k0, _k1)           \
 Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, _array + addr(_k1<<1));  \
 Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));  \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_inter_addsub_p(_t0, _t1, _pd0, _pd1, _array, _k0, _k1) \
 _pd1 = _array + addr(_k1<<1);                                           \
 _pd0 = _array + addr(_k0<<1);                                           \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, _pd1);                       \
 prefetch_p (_pd1);                                                      \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _pd0);                       \
 prefetch_p (_pd0);                                                      \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_addsub_p(_t0, _t1, _pd0, _pd1, _array, _k0, _k1) \
 _pd1 = _array + addr(_k1<<1);                                     \
 _pd0 = _array + addr(_k0<<1);                                     \
 Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, _pd1);                       \
 prefetch_p (_pd1);                                                \
 Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _pd0);                       \
 prefetch_p (_pd0);                                                \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_inter_addsub_pp(_t0, _t1, _index, _pd0, _pd1)\
 Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, (_pd1 + _index));  \
 prefetch_p (_pd1 + _index);                                   \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i,(_pd0 + _index));   \
 prefetch_p (_pd0 + _index);                                   \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_addsub_pp(_t0, _t1, _index, _pd0, _pd1)\
 Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, (_pd1 + _index));  \
 prefetch_p (_pd1 + _index);                             \
 Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, (_pd0 + _index));  \
 prefetch_p (_pd0 + _index);                             \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_inter_addsub_pp_no_fetch(_t0, _t1, _index, _pd0, _pd1)\
 Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, (_pd1 + _index));  \
 Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, (_pd0 + _index));  \
 sse2_asm_addsub (_t0, _t1);

# define sse2_load_addsub_pp_no_fetch(_t0, _t1, _index, _pd0, _pd1)\
 Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, (_pd1 + _index));  \
 Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, (_pd0 + _index));  \
 sse2_asm_addsub (_t0, _t1);

# define sse2_addsub_store_inter(_k0, _k1, _t0, _t1)       \
 sse2_asm_addsub (_t0, _t1);                              \
 Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i); \
 Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \

# define sse2_addsub_store(_k0, _k1, _t0, _t1)                 \
 sse2_asm_addsub (_t0, _t1);                                  \
 Y_MM_STORE_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i); \
 Y_MM_STORE_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i);

# define sse2_addsub_store_inter_p(_index, _pd0, _pd1, _t0, _t1) \
 sse2_asm_addsub (_t0, _t1);                                    \
 Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);    \
 Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);

# define sse2_addsub_store_p(_index, _pd0, _pd1, _t0, _t1) \
 sse2_asm_addsub (_t0, _t1);                              \
 Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);    \
 Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);

# define sse2_addsub_store_inter_no_next(_pd, _pu, _i, _j, _t0, _t1)  \
 sse2_asm_addsub (_t0, _t1);                                         \
 Y_MM_STORE_INTER_PAIR( _pd + _j, _pu + _j, _t1##r, _t1##i);         \
 Y_MM_STORE_INTER_PAIR( _pd + _i, _pu + _i, _t0##r, _t0##i);

#else /* ASM_ADDSUB */

#if defined ALT_ADDSUB
# define sse2_addsub(__t0,__t1)                \
   Y_MM_ADD_PD( __t0##r, __t0##r, __t1##r);    \
   Y_MM_MUL_PD( __t1##r, __t1##r, MM_MTWO);    \
   Y_MM_ADD_PD( __t0##i, __t0##i, __t1##i);    \
   Y_MM_MUL_PD( __t1##i, __t1##i, MM_MTWO);    \
   Y_MM_ADD_PD( __t1##r, __t1##r, __t0##r);    \
   Y_MM_ADD_PD( __t1##i, __t1##i, __t0##i)

#else

# define sse2_addsub(_t0,_t1)           \
{                                       \
   Y__M128D _t0r, _t0i;                 \
   Y_MM_MOV_PD( _t0r, _t0##r);          \
   Y_MM_MOV_PD( _t0i, _t0##i);          \
   Y_MM_SUB_PD( _t0r, _t0r, _t1##r);    \
   Y_MM_SUB_PD( _t0i, _t0i, _t1##i);    \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);\
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);\
   Y_MM_MOV_PD( _t1##r, _t0r);          \
   Y_MM_MOV_PD( _t1##i, _t0i);          \
}

#endif

#if defined(ALT_ADDSUB)
# define sse2_load_inter_addsub(_t0,_t1,_array,_k0,_k1)              \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));\
   Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, _array + addr(_k1<<1));\
   sse2_addsub (_t0, _t1)

# define sse2_load_addsub(_t0,_t1,_array,_k0,_k1)              \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));\
   Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, _array + addr(_k1<<1));\
   sse2_addsub (_t0, _t1)

# define sse2_load_inter_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1)   \
   _pd0 = _array + addr(_k0<<1);                                     \
   _pd1 = _array + addr(_k1<<1);                                     \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _pd0);                 \
   prefetch_p (_pd0);                                                \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i, _pd1);                 \
   prefetch_p (_pd1);                                                \
   sse2_addsub (_t0, _t1)

# define sse2_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1)   \
   _pd0 = _array + addr(_k0<<1);                                \
   _pd1 = _array + addr(_k1<<1);                                \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _pd0);                  \
   prefetch_p (_pd0);                                           \
   Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i, _pd1);                  \
   prefetch_p (_pd1);                                           \
   sse2_addsub (_t0, _t1)

# define sse2_load_inter_addsub_pp(_t0,_t1,_index,_pd0,_pd1)  \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i,(_pd0 + _index));\
   prefetch_p (_pd0 + _index);                                \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t1##r, _t1##i,(_pd1 + _index));\
   prefetch_p (_pd1 + _index);                                \
   sse2_addsub (_t0, _t1)

# define sse2_load_addsub_pp(_t0,_t1,_index,_pd0,_pd1)  \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i,(_pd0 + _index));\
   prefetch_p (_pd0 + _index);                          \
   Y_MM_LOAD_PAIR_NEXT( _t1##r, _t1##i,(_pd1 + _index));\
   prefetch_p (_pd1 + _index);                          \
   sse2_addsub (_t0, _t1)

#else
# define sse2_load_inter_addsub(_t0,_t1,_array,_k0,_k1)               \
{                                                                    \
   Y__M128D __a, __b;                                                \
   Y_MM_LOAD_INTER_PAIR_NEXT( __a, __b, _array + addr(_k1<<1));      \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));\
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                \
}

# define sse2_load_addsub(_t0,_t1,_array,_k0,_k1)               \
{                                                              \
   Y__M128D __a, __b;                                          \
   Y_MM_LOAD_PAIR_NEXT( __a, __b, _array + addr(_k1<<1));      \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _array + addr(_k0<<1));\
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                          \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                          \
}


# define sse2_load_inter_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1)   \
{                                                                    \
   Y__M128D __a, __b;                                                \
   _pd1 = _array + addr(_k1<<1);                                     \
   _pd0 = _array + addr(_k0<<1);                                     \
   Y_MM_LOAD_INTER_PAIR_NEXT( __a, __b, _pd1);                       \
   prefetch_p (_pd1);                                                \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _pd0);                 \
   prefetch_p (_pd0);                                                \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                \
}

# define sse2_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1)   \
{                                                              \
   Y__M128D __a, __b;                                          \
   _pd1 = _array + addr(_k1<<1);                               \
   _pd0 = _array + addr(_k0<<1);                               \
   Y_MM_LOAD_PAIR_NEXT( __a, __b, _pd1);                       \
   prefetch_p (_pd1);                                          \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _pd0);                 \
   prefetch_p (_pd0);                                          \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                          \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                          \
}

# define sse2_load_inter_addsub_pp(_t0,_t1,_index,_pd0,_pd1)   \
{                                                             \
   Y__M128D __a, __b;                                         \
   Y_MM_LOAD_INTER_PAIR_NEXT( __a, __b, (_pd1 + _index));     \
   prefetch_p (_pd1 + _index);                                \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i,(_pd0 + _index));\
   prefetch_p (_pd0 + _index);                                \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                         \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                         \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                         \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                         \
}

# define sse2_load_addsub_pp(_t0,_t1,_index,_pd0,_pd1)   \
{                                                       \
   Y__M128D __a, __b;                                   \
   Y_MM_LOAD_PAIR_NEXT( __a, __b, (_pd1 + _index));     \
   prefetch_p (_pd1 + _index);                          \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i,(_pd0 + _index));\
   prefetch_p (_pd0 + _index);                          \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                   \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                   \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                   \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                   \
}
#endif

#define sse2_load_inter_addsub_pp_no_fetch(_t0,_t1,_index,_pd0,_pd1)   \
{                                                             \
   Y__M128D __a, __b;                                         \
   Y_MM_LOAD_INTER_PAIR_NEXT( __a, __b, _pd1 + _index);       \
   Y_MM_LOAD_INTER_PAIR_NEXT( _t0##r, _t0##i, _pd0 + _index); \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                         \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                         \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                         \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                         \
}

#define sse2_load_addsub_pp_no_fetch(_t0,_t1,_index,_pd0,_pd1)   \
{                                                       \
   Y__M128D __a, __b;                                   \
   Y_MM_LOAD_PAIR_NEXT( __a, __b, _pd1 + _index);       \
   Y_MM_LOAD_PAIR_NEXT( _t0##r, _t0##i, _pd0 + _index); \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                   \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                   \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                   \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                   \
}

#if defined(ALT_ADDSUB)
# define sse2_addsub_store_inter(_k0,_k1,_t0,_t1)                     \
   sse2_addsub (_t0,_t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i);

# define sse2_addsub_store(_k0,_k1,_t0,_t1)                     \
   sse2_addsub (_t0,_t1);                                       \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1)   \
   sse2_addsub (_t0,_t1);                                     \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_addsub_store_p(_index,_pd0,_pd1,_t0,_t1)   \
   sse2_addsub (_t0,_t1);                               \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#else

# define sse2_addsub_store_inter(_k0,_k1,_t0,_t1)                     \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_SUB_PD( __a, _t0##r, _t1##r);                                 \
   Y_MM_SUB_PD( __b, _t0##i, _t1##i);                                 \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);                              \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);                              \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k1<<1), __a, __b);       \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
}

# define sse2_addsub_store(_k0,_k1,_t0,_t1)                      \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_SUB_PD( __a, _t0##r, _t1##r);                           \
   Y_MM_SUB_PD( __b, _t0##i, _t1##i);                           \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);                        \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);                        \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k1<<1), __a, __b);       \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
}

#define sse2_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1)   \
{                                                             \
   Y__M128D __a, __b;                                         \
   Y_MM_SUB_PD( __a, _t0##r, _t1##r);                         \
   Y_MM_SUB_PD( __b, _t0##i, _t1##i);                         \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);                      \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);                      \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, __a, __b);      \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
}

#define sse2_addsub_store_p(_index,_pd0,_pd1,_t0,_t1)   \
{                                                       \
   Y__M128D __a, __b;                                   \
   Y_MM_SUB_PD( __a, _t0##r, _t1##r);                   \
   Y_MM_SUB_PD( __b, _t0##i, _t1##i);                   \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);                \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);                \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, __a, __b);      \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
}
#endif

#define sse2_addsub_store_inter_no_next( _pd, _pu, _i, _j, _t0, _t1)   \
{                                                                      \
   Y__M128D __a, __b;                                                  \
   Y_MM_SUB_PD( __a, _t0##r, _t1##r);                                  \
   Y_MM_SUB_PD( __b, _t0##i, _t1##i);                                  \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##r);                               \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##i);                               \
   Y_MM_STORE_INTER_PAIR( _pd + _j, _pu + _j, __a, __b);               \
   Y_MM_STORE_INTER_PAIR( _pd + _i, _pu + _i, _t0##r, _t0##i);         \
}

#endif /* ASM_ADDSUB */

/**********************************************************************
   cplx_mulmul macro:
		_a= _t0 * _f0;
		_b= _t1 * _f1;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/

#if defined(ASM_MULMUL)

# define sse2_asm_mulmul(__t0, __t1, __f0, __f1)       \
__asm__ volatile ("movapd %0, %%xmm5             \n"   \
"        mulpd %5, %%xmm5                        \n"   \
"        mulpd %4, %0                            \n"   \
"        movapd %1, %%xmm6                       \n"   \
"        mulpd  %5, %%xmm6                       \n"   \
"        mulpd  %4, %1                           \n"   \
"        movapd  %2, %%xmm7                      \n"   \
"        mulpd  %7, %%xmm7                       \n"   \
"        subpd  %%xmm6, %0                       \n"   \
"        mulpd  %6, %2                           \n"   \
"        movapd  %3, %%xmm6                      \n"   \
"        mulpd  %7, %%xmm6                       \n"   \
"        addpd  %%xmm5, %1                       \n"   \
"        mulpd  %6, %3                           \n"   \
"        subpd  %%xmm6, %2                       \n"   \
"        addpd  %%xmm7, %3                       \n"   \
: "+&x" (__t0##r), "+&x" (__t0##i), "+&x" (__t1##r), "+&x" (__t1##i)  \
: "m" (__f0##r), "m" (__f0##i), "m" (__f1##r), "m" (__f1##i) \
: "xmm5", "xmm6", "xmm7" )


#define sse2_mulmul(_t0, _t1, _f0, _f1)        \
 sse2_asm_mulmul( _t0, _t1, _f0, _f1);

#define sse2_load_inter_mulmul_p(_t0, _t1, _pd0, _pd1, _f0, _f1, _array, _k0, _k1) \
   _pd0 = _array + addr((_k0)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0);              \
   prefetch_p( _pd0);                                            \
   _pd1 = _array + addr((_k1)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1);              \
   prefetch_p( _pd1);                                            \
   sse2_asm_mulmul( _t0, _t1, _f0, _f1);

#define sse2_load_mulmul_p(_t0, _t1, _pd0, _pd1, _f0, _f1, _array, _k0, _k1) \
   _pd0 = _array + addr((_k0)<<1);                               \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0);                    \
   prefetch_p( _pd0);                                            \
   _pd1 = _array + addr((_k1)<<1);                               \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1);                    \
   prefetch_p( _pd1);                                            \
   sse2_asm_mulmul( _t0, _t1, _f0, _f1);

#define sse2_load_inter_mulmul_pp(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);     \
   prefetch_p( _pd0 + _index);                                   \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);     \
   prefetch_p( _pd1 + _index);                                   \
   sse2_asm_mulmul( _t0, _t1, _f0, _f1);

#define sse2_load_mulmul_pp(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);              \
   prefetch_p( _pd0 + _index);                                      \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);              \
   prefetch_p( _pd1 + _index);                                      \
   sse2_asm_mulmul( _t0, _t1, _f0, _f1);


#else /* ASM_MULMUL */

#define sse2_mulmul(_t0,_t1,_f0,_f1)        \
   sse2_mul0(_t0, _f0);                     \
   sse2_mul0(_t1, _f1);

#define sse2_load_inter_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
   _pd0 = _array + addr((_k0)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0);              \
   prefetch_p( _pd0);                                            \
   sse2_mul0(_t0, _f0);                                          \
   _pd1 = _array + addr((_k1)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1);              \
   prefetch_p( _pd1);                                            \
   sse2_mul0(_t1, _f1);

#define sse2_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
   _pd0 = _array + addr((_k0)<<1);                                   \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0);                        \
   prefetch_p( _pd0);                                                \
   sse2_mul0(_t0, _f0);                                              \
   _pd1 = _array + addr((_k1)<<1);                                   \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1);                        \
   prefetch_p( _pd1);                                                \
   sse2_mul0(_t1, _f1);

/*
#define sse2_mulmul(_t0,_t1,_f0,_f1)        \
{                                           \
   Y__M128D __a,__b,__c,__d;                \
   Y_MM_MUL_PD( __a, _t0##r, _f0##r);       \
   Y_MM_MUL_PD( __b, _t0##i, _f0##i);       \
   Y_MM_MUL_PD( __c, _t1##r, _f1##r);       \
   Y_MM_MUL_PD( __d, _t1##i, _f1##i);       \
   Y_MM_SUB_PD( __a, __a, __b);             \
   Y_MM_SUB_PD( __c, __c, __d);             \
   Y_MM_MUL_PD( __b, _t0##r, _f0##i);       \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);       \
   Y_MM_MUL_PD( _t0##i, _t0##i, _f0##r);    \
   Y_MM_MUL_PD( _t1##i, _t1##i, _f1##r);    \
   Y_MM_MOV_PD( _t0##r, __a);               \
   Y_MM_MOV_PD( _t1##r, __c);               \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);       \
   Y_MM_ADD_PD( _t1##i, _t1##i, __d);       \
}

#define sse2_load_inter_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
{                                                                \
   Y__M128D __a, __b, __c, __d;                                  \
   _pd0 = _array + addr((_k0)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0);              \
   prefetch_p( _pd0);                                            \
   Y_MM_MUL_PD( __a, _t0##i, _f0##i);                            \
   Y_MM_MUL_PD( __b, _t0##r, _f0##i);                            \
   Y_MM_MUL_PD( _t0##r, _t0##r, _f0##r);                         \
   Y_MM_MUL_PD( _t0##i, _t0##i, _f0##r);                         \
   _pd1 = _array + addr((_k1)<<1);                               \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1);              \
   prefetch_p( _pd1);                                            \
   Y_MM_MUL_PD( __c, _t1##i, _f1##i);                            \
   Y_MM_SUB_PD( _t0##r, _t0##r, __a);                            \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);                            \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                            \
   Y_MM_MUL_PD( _t1##r, _t1##r, _f1##r);                         \
   Y_MM_MUL_PD( _t1##i, _t1##i, _f1##r);                         \
   Y_MM_SUB_PD( _t1##r, _t1##r, __c);                            \
   Y_MM_ADD_PD( _t1##i, _t1##i, __d);                            \
}

#define sse2_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
{                                                                \
   Y__M128D __a, __b, __c, __d;                                  \
   _pd0 = _array + addr((_k0)<<1);                               \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0);                    \
   prefetch_p( _pd0);                                            \
   Y_MM_MUL_PD( __a, _t0##i, _f0##i);                            \
   Y_MM_MUL_PD( __b, _t0##r, _f0##i);                            \
   Y_MM_MUL_PD( _t0##r, _t0##r, _f0##r);                         \
   Y_MM_MUL_PD( _t0##i, _t0##i, _f0##r);                         \
   _pd1 = _array + addr((_k1)<<1);                               \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1);                    \
   prefetch_p( _pd1);                                            \
   Y_MM_MUL_PD( __c, _t1##i, _f1##i);                            \
   Y_MM_SUB_PD( _t0##r, _t0##r, __a);                            \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);                            \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                            \
   Y_MM_MUL_PD( _t1##r, _t1##r, _f1##r);                         \
   Y_MM_MUL_PD( _t1##i, _t1##i, _f1##r);                         \
   Y_MM_SUB_PD( _t1##r, _t1##r, __c);                            \
   Y_MM_ADD_PD( _t1##i, _t1##i, __d);                            \
}
*/
/*
#define sse2_load_inter_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
{                                                                \
   Y__M128D __a, __b, __c, __d;                                  \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);     \
   prefetch_p( _pd0 + _index);                                   \
   Y_MM_MUL_PD( __a, _t0##i, _f0##i);                            \
   Y_MM_MUL_PD( __b, _t0##r, _f0##i);                            \
   Y_MM_MUL_PD( _t0##r, _t0##r, _f0##r);                         \
   Y_MM_MUL_PD( _t0##i, _t0##i, _f0##r);                         \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);     \
   prefetch_p( _pd1 + _index);                                   \
   Y_MM_MUL_PD( __c, _t1##i, _f1##i);                            \
   Y_MM_SUB_PD( _t0##r, _t0##r, __a);                            \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);                            \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                            \
   Y_MM_MUL_PD( _t1##r, _t1##r, _f1##r);                         \
   Y_MM_MUL_PD( _t1##i, _t1##i, _f1##r);                         \
   Y_MM_SUB_PD( _t1##r, _t1##r, __c);                            \
   Y_MM_ADD_PD( _t1##i, _t1##i, __d);                            \
}
*/

#define sse2_load_inter_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);        \
   prefetch_p( _pd0 + _index);                                      \
   sse2_mul0(_t0, _f0);                                             \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);        \
   prefetch_p( _pd1 + _index);                                      \
   sse2_mul0(_t1, _f1);

#define sse2_load_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);        \
   prefetch_p( _pd0 + _index);                                \
   sse2_mul0(_t0, _f0);                                       \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);        \
   prefetch_p( _pd1 + _index);                                \
   sse2_mul0(_t1, _f1);

/*
#define sse2_load_inter_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
{                                                             \
   Y__M128D __a, __b, __c, __d, _t0r, _t0i, _t1r, _t1i;       \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0r, _t0i, _pd0 + _index);            \
   prefetch_p( _pd0 + _index);                                \
   Y_MM_MUL_PD( __a, _t0i, _f0##i);                           \
   Y_MM_MUL_PD( __b, _t0r, _f0##i);                           \
   Y_MM_MUL_PD( _t0r, _t0r, _f0##r);                          \
   Y_MM_MUL_PD( _t0i, _t0i, _f0##r);                          \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1r, _t1i, _pd1 + _index);            \
   prefetch_p( _pd1 + _index);                                \
   Y_MM_MUL_PD( __c, _t1i, _f1##i);                           \
   Y_MM_SUB_PD( _t0r, _t0r, __a);                             \
   Y_MM_MUL_PD( __d, _t1r, _f1##i);                           \
   Y_MM_ADD_PD( _t0i, _t0i, __b);                             \
   Y_MM_MUL_PD( _t1r, _t1r, _f1##r);                          \
   Y_MM_MUL_PD( _t1i, _t1i, _f1##r);                          \
   Y_MM_SUB_PD( _t1r, _t1r, __c);                             \
   Y_MM_ADD_PD( _t1i, _t1i, __d);                             \
   Y_MM_MOV_PD( _t0##r, _t0r);                                \
   Y_MM_MOV_PD( _t0##i, _t0i);                                \
   Y_MM_MOV_PD( _t1##r, _t1r);                                \
   Y_MM_MOV_PD( _t1##i, _t1i);                                \
}


#define sse2_load_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
{                                                             \
   Y__M128D __a, __b, __c, __d, _t0r, _t0i, _t1r, _t1i;       \
   Y_MM_LOAD_PAIR_NEXT(_t0r, _t0i, _pd0 + _index);            \
   prefetch_p( _pd0 + _index);                                \
   Y_MM_MUL_PD( __a, _t0i, _f0##i);                           \
   Y_MM_MUL_PD( __b, _t0r, _f0##i);                           \
   Y_MM_MUL_PD( _t0r, _t0r, _f0##r);                          \
   Y_MM_MUL_PD( _t0i, _t0i, _f0##r);                          \
   Y_MM_LOAD_PAIR_NEXT(_t1r, _t1i, _pd1 + _index);            \
   prefetch_p( _pd1 + _index);                                \
   Y_MM_MUL_PD( __c, _t1i, _f1##i);                           \
   Y_MM_SUB_PD( _t0r, _t0r, __a);                             \
   Y_MM_MUL_PD( __d, _t1r, _f1##i);                           \
   Y_MM_ADD_PD( _t0i, _t0i, __b);                             \
   Y_MM_MUL_PD( _t1r, _t1r, _f1##r);                          \
   Y_MM_MUL_PD( _t1i, _t1i, _f1##r);                          \
   Y_MM_SUB_PD( _t1r, _t1r, __c);                             \
   Y_MM_ADD_PD( _t1i, _t1i, __d);                             \
   Y_MM_MOV_PD( _t0##r, _t0r);                                \
   Y_MM_MOV_PD( _t0##i, _t0i);                                \
   Y_MM_MOV_PD( _t1##r, _t1r);                                \
   Y_MM_MOV_PD( _t1##i, _t1i);                                \
}
*/

#endif /* ASM_MULMUL */

#if defined(ASM_DIVDIV)

# define sse2_asm_divdiv(__t0, __t1, __f0, __f1)       \
__asm__ volatile ("movapd %0, %%xmm5             \n"   \
"        mulpd %5, %%xmm5                        \n"   \
"        mulpd %4, %0                            \n"   \
"        movapd %1, %%xmm6                       \n"   \
"        mulpd  %5, %%xmm6                       \n"   \
"        mulpd  %4, %1                           \n"   \
"        movapd  %2, %%xmm7                      \n"   \
"        mulpd  %7, %%xmm7                       \n"   \
"        addpd  %%xmm6, %0                       \n"   \
"        mulpd  %6, %2                           \n"   \
"        movapd  %3, %%xmm6                      \n"   \
"        mulpd  %7, %%xmm6                       \n"   \
"        subpd  %%xmm5, %1                       \n"   \
"        mulpd  %6, %3                           \n"   \
"        addpd  %%xmm6, %2                       \n"   \
"        subpd  %%xmm7, %3                       \n"   \
: "+&x" (__t0##r), "+&x" (__t0##i), "+&x" (__t1##r), "+&x" (__t1##i)  \
: "m" (__f0##r), "m" (__f0##i), "m" (__f1##r), "m" (__f1##i) \
: "xmm5", "xmm6", "xmm7" )

#define sse2_divdiv(_t0, _t1, _f0, _f1)        \
 sse2_asm_divdiv (_t0, _t1, _f0, _f1);

#else /* ASM_DIVDIV */

#define sse2_divdiv(_t0,_t1,_f0,_f1)        \
{                                           \
   Y__M128D __a,__b,__c,__d;                \
   Y_MM_MUL_PD( __a, _t0##r, _f0##r);       \
   Y_MM_MUL_PD( __b, _t0##i, _f0##i);       \
   Y_MM_MUL_PD( __c, _t1##r, _f1##r);       \
   Y_MM_MUL_PD( __d, _t1##i, _f1##i);       \
   Y_MM_ADD_PD( __a, __a, __b);             \
   Y_MM_ADD_PD( __c, __c, __d);             \
   Y_MM_MUL_PD( __b, _t0##r, _f0##i);       \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);       \
   Y_MM_MUL_PD( _t0##i, _t0##i, _f0##r);    \
   Y_MM_MUL_PD( _t1##i, _t1##i, _f1##r);    \
   Y_MM_MOV_PD( _t0##r, __a);               \
   Y_MM_MOV_PD( _t1##r, __c);               \
   Y_MM_SUB_PD( _t0##i, _t0##i, __b);       \
   Y_MM_SUB_PD( _t1##i, _t1##i, __d);       \
}

#endif /* ASM_DIVDIV */

#if defined(ASM_DIVMUL)

#define sse2_asm_divmul( __t0, __t1, __s0, __s1) \
__asm__ volatile (" movapd %4, %0            \n" \
"        movapd %5, %1                       \n" \
"        movapd %5, %%xmm6                   \n" \
"        movapd %4, %%xmm7                   \n" \
"        mulpd %7, %%xmm6                    \n" \
"        mulpd %7, %%xmm7                    \n" \
"        mulpd %6, %0                        \n" \
"        mulpd %6, %1                        \n" \
"        movapd %0, %2                       \n" \
"        movapd %1, %3                       \n" \
"        subpd %%xmm6, %2                    \n" \
"        addpd %%xmm6, %0                    \n" \
"        addpd %%xmm7, %3                    \n" \
"        subpd %%xmm7, %1                    \n" \
: "=x" (__t0##r), "=x" (__t0##i), "=x" (__t1##r), "=x" (__t1##i)  \
: "m" (__s0##r), "m" (__s0##i), "m" (__s1##r), "m" (__s1##i) \
: "xmm6", "xmm7" )

#define sse2_divmul( _t0, _t1, _s0, _s1) \
 sse2_asm_divmul( _t0, _t1, _s0, _s1);

#else /* ASM_DIVMUL */

#define sse2_divmul( _t0, _t1, _s0, _s1) \
{                                        \
   Y__M128D __x0,__x1;                   \
   Y_MM_MOV_PD( _t0##r, _s0##r);         \
   Y_MM_MOV_PD( __x0, _s0##i);           \
   Y_MM_MOV_PD( _t0##i, _s0##i);         \
   Y_MM_MOV_PD( __x1,  _s0##r);          \
   Y_MM_MUL_PD( _t0##r, _t0##r, _s1##r); \
   Y_MM_MUL_PD( __x0, __x0, _s1##i);     \
   Y_MM_MUL_PD( _t0##i, _t0##i, _s1##r); \
   Y_MM_MOV_PD( _t1##r, _t0##r);         \
   Y_MM_MOV_PD( _t1##i, _t0##i);         \
   Y_MM_MUL_PD( __x1, __x1, _s1##i);     \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0 );  \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0 );  \
   Y_MM_ADD_PD( _t1##i, _t1##i, __x1 );  \
   Y_MM_SUB_PD( _t0##i, _t0##i, __x1 );  \
}

#endif /* ASM_DIVMUL */

/**********************************************************************
   cplx_load_mulmuladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/
#if defined(ASM_MULADDSUB)

# define sse2_asm_muladdsub(__t0, __t1, __f)\
   __asm__ volatile ("movapd %2, %%xmm6 \n" \
"        movapd %3, %%xmm7      \n"         \
"        mulpd %5, %%xmm6       \n"         \
"        mulpd %5, %%xmm7       \n"         \
"        mulpd %4, %2           \n"         \
"        mulpd %4, %3           \n"         \
"        subpd %%xmm7, %2       \n"         \
"        addpd %%xmm6, %3       \n"         \
"        movapd %0, %%xmm6      \n"         \
"        movapd %1, %%xmm7      \n"         \
"        subpd %2, %%xmm6       \n"         \
"        subpd %3, %%xmm7       \n"         \
"        addpd %2, %0           \n"         \
"        addpd %3, %1           \n"         \
"        movapd %%xmm6, %2      \n"         \
"        movapd %%xmm7, %3      \n"         \
  : "+&x" (__t0##r), "+&x" (__t0##i), "+&x" (__t1##r), "+&x" (__t1##i) \
  :"m" (__f##r) , "m" (__f##i): "xmm6",  "xmm7");

# define sse2_muladdsub(_t0, _t1, _f)      \
  sse2_asm_muladdsub (_t0, _t1, _f);

#define sse2_load_inter_muladdsub(_t0, _t1, _f, _array, _k0, _k1) \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));\
  prefetch_p(_array + addr(_k1<<1));                              \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));\
  prefetch_p(_array + addr(_k0<<1));                              \
  sse2_asm_muladdsub(_t0, _t1, _f);                               \

#define sse2_load_muladdsub(_t0, _t1, _f, _array, _k0, _k1) \
  Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));\
  prefetch_p(_array + addr(_k1<<1));                        \
  Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));\
  prefetch_p(_array + addr(_k0<<1));                        \
  sse2_asm_muladdsub(_t0, _t1, _f);                         \

#define sse2_load_inter_muladdsub_p(_t0, _t1, _pd0, _pd1, _f, _array, _k0, _k1)\
   _pd1 = _array + addr((_k1)<<1);                                             \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1);                            \
   prefetch_p(_pd1);                                                           \
   _pd0 = _array + addr((_k0)<<1);                                             \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0);                            \
   prefetch_p(_pd0);                                                           \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_muladdsub_p(_t0, _t1, _pd0, _pd1, _f, _array, _k0, _k1)\
   _pd1 = _array + addr((_k1)<<1);                                       \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1);                            \
   prefetch_p(_pd1);                                                     \
   _pd0 = _array + addr((_k0)<<1);                                       \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0);                            \
   prefetch_p(_pd0);                                                     \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_inter_muladdsub_pp(_t0, _t1, _f, _index, _pd0, _pd1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);           \
   prefetch_p( _pd1 + _index);                                         \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);           \
   prefetch_p( _pd0 + _index);                                         \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_muladdsub_pp(_t0, _t1, _f, _index, _pd0, _pd1) \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);           \
   prefetch_p( _pd1 + _index);                                   \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);           \
   prefetch_p( _pd0 + _index);                                   \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_inter_muladdsub_pp_no_fetch(_t0, _t1, _f, _index, _pd0, _pd1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                    \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                    \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_muladdsub_pp_no_fetch(_t0, _t1, _f, _index, _pd0, _pd1) \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                    \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                    \
   sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_load_inter_muladdsub_no_next(_t0, _t1, _f, _ii, _jj, _pd0, _pd1) \
  Y_MM_LOAD_INTER_PAIR (_t1##r, _t1##i, (_pd0 + _jj), (_pd1 + _jj));          \
  prefetch_p ( _pd1 + _jj);                                                   \
  prefetch_p ( _pd0 + _jj);                                                   \
  Y_MM_LOAD_INTER_PAIR (_t0##r, _t0##i, (_pd0 + _ii), (_pd1 + _ii));          \
  prefetch_p ( _pd1 + _ii);                                                   \
  prefetch_p ( _pd0 + _ii);                                                   \
  sse2_asm_muladdsub(_t0, _t1, _f);

#define sse2_muladdsub_store_inter(_array, _k0, _k1, _t0, _t1, _f)    \
   sse2_asm_muladdsub (_t0, _t1, _f);                                 \
   Y_MM_STORE_INTER_PAIR_NEXT (_array + addr(_k0<<1), _t0##r, _t0##i);\
   Y_MM_STORE_INTER_PAIR_NEXT (_array + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_muladdsub_store(_array, _k0, _k1, _t0, _t1, _f)    \
   sse2_asm_muladdsub (_t0, _t1, _f);                           \
   Y_MM_STORE_PAIR_NEXT (_array + addr(_k0<<1), _t0##r, _t0##i);\
   Y_MM_STORE_PAIR_NEXT (_array + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_muladdsub_store_inter_p(_index, _pd0, _pd1, _t0, _t1, _f) \
   sse2_asm_muladdsub (_t0, _t1, _f);                                  \
   Y_MM_STORE_INTER_PAIR_NEXT (_pd0 + _index, _t0##r, _t0##i);         \
   Y_MM_STORE_INTER_PAIR_NEXT (_pd1 + _index, _t1##r, _t1##i);

#define sse2_muladdsub_store_p(_index, _pd0, _pd1, _t0, _t1, _f) \
   sse2_asm_muladdsub (_t0, _t1, _f);                            \
   Y_MM_STORE_PAIR_NEXT (_pd0 + _index, _t0##r, _t0##i);         \
   Y_MM_STORE_PAIR_NEXT (_pd1 + _index, _t1##r, _t1##i);

#else /* ASM_MULADDSUB */

# define sse2_muladdsub(_t0,_t1,_f)      \
  sse2_mul0(_t1, _f);                    \
  sse2_addsub(_t0, _t1);

#define sse2_load_inter_muladdsub(_t0,_t1,_f,_array,_k0,_k1)       \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));\
   prefetch_p(_array + addr(_k1<<1));                              \
   sse2_mul0(_t1, _f);                                             \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));\
   prefetch_p(_array + addr(_k0<<1));                              \
   sse2_addsub(_t0, _t1);

#define sse2_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1)       \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));\
   prefetch_p(_array + addr(_k1<<1));                        \
   sse2_mul0(_t1, _f);                                       \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));\
   prefetch_p(_array + addr(_k0<<1));                        \
   sse2_addsub(_t0, _t1);

#define sse2_load_inter_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
   _pd1= _array + addr((_k1)<<1);                                        \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1);                      \
   prefetch_p(_pd1);                                                     \
   sse2_mul0(_t1, _f);                                                   \
   _pd0= _array + addr((_k1)<<1);                                        \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0);                      \
   prefetch_p(_pd0);                                                     \
   sse2_addsub(_t0, _t1);

/*
# define sse2_muladdsub(_t0,_t1,_f)      \
{                                        \
   Y__M128D __ar,__ai;                   \
   Y_MM_MUL_PD(__ar, _t1##r, _f##i);     \
   Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);   \
   Y_MM_MUL_PD(__ai, _t1##i, _f##i);     \
   Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);   \
   Y_MM_SUB_PD(_t1##r, _t1##r, __ai);    \
   Y_MM_ADD_PD(_t1##i, _t1##i, __ar);    \
   sse2_addsub(_t0,_t1);                 \
}   


#define sse2_load_inter_muladdsub(_t0,_t1,_f,_array,_k0,_k1)   \
{                                                              \
   Y__M128D __ar,__ai,__br,__bi,__cr,__ci;                     \
   Y_MM_LOAD_INTER_PAIR_NEXT(__ai, __bi, array + addr(_k1<<1));\
   prefetch_p(_array + addr(_k1<<1));                          \
   Y_MM_MUL_PD(__ar, __ai, _f##r);                             \
   Y_MM_MUL_PD(__ai, __ai, _f##i);                             \
   Y_MM_MUL_PD(__br, __bi, _f##i);                             \
   Y_MM_MUL_PD(__bi, __bi, _f##r);                             \
   Y_MM_LOAD_INTER_PAIR_NEXT(__cr, __ci, array + addr(_k0<<1));\
   prefetch_p(_array + addr(_k0<<1));                          \
   Y_MM_SUB_PD(__ar, __ar, __br);                              \
   Y_MM_ADD_PD(__ai, __ai, __bi);                              \
   Y_MM_ADD_PD(_t0##r, __cr, __ar);                            \
   Y_MM_SUB_PD(_t1##r, __cr, __ar);                            \
   Y_MM_ADD_PD(_t0##i, __ci, __ai);                            \
   Y_MM_SUB_PD(_t1##i, __ci, __ai);                            \
}


#define sse2_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1)   \
{                                                        \
   Y__M128D __ar,__ai,__br,__bi,__cr,__ci;               \
   Y_MM_LOAD_PAIR_NEXT(__ai, __bi, array + addr(_k1<<1));\
   prefetch_p(_array + addr(_k1<<1));                    \
   Y_MM_MUL_PD(__br, __bi, _f##i);                       \
   Y_MM_MUL_PD(__bi, __bi, _f##r);                       \
   Y_MM_MUL_PD(__ar, __ai, _f##r);                       \
   Y_MM_MUL_PD(__ai, __ai, _f##i);                       \
   Y_MM_LOAD_PAIR_NEXT(__cr, __ci, array + addr(_k0<<1));\
   prefetch_p(_array + addr(_k0<<1));                    \
   Y_MM_SUB_PD(__ar, __ar, __br);                        \
   Y_MM_ADD_PD(__ai, __ai, __bi);                        \
   Y_MM_ADD_PD(_t0##r, __cr, __ar);                      \
   Y_MM_SUB_PD(_t1##r, __cr, __ar);                      \
   Y_MM_ADD_PD(_t0##i, __ci, __ai);                      \
   Y_MM_SUB_PD(_t1##i, __ci, __ai);                      \
}


#define sse2_load_inter_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
{                                                                        \
   Y__M128D __ar,__ai,__br,__bi,__cr,__ci;                               \
   _pd1= _array + addr((_k1)<<1);                                        \
   Y_MM_LOAD_INTER_PAIR_NEXT(__ai, __bi, _pd1);                          \
   prefetch_p(_pd1);                                                     \
   Y_MM_MUL_PD(__ar, __ai, _f##r);                                       \
   Y_MM_MUL_PD(__ai, __ai, _f##i);                                       \
   Y_MM_MUL_PD(__br, __bi, _f##i);                                       \
   Y_MM_MUL_PD(__bi, __bi, _f##r);                                       \
   _pd0= _array + addr((_k1)<<1);                                        \
   Y_MM_LOAD_INTER_PAIR_NEXT(__cr, __ci, _pd0);                          \
   prefetch_p(_pd0);                                                     \
   Y_MM_SUB_PD(__ar, __ar, __br);                                        \
   Y_MM_ADD_PD(__ai, __ai, __bi);                                        \
   Y_MM_ADD_PD(_t0##r, __cr, __ar);                                      \
   Y_MM_SUB_PD(_t1##r, __cr, __ar);                                      \
   Y_MM_ADD_PD(_t0##i, __ci, __ai);                                      \
   Y_MM_SUB_PD(_t1##i, __ci, __ai);                                      \
}

#define sse2_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
{                                                                  \
   Y__M128D __ar,__ai,__br,__bi,__cr,__ci;                         \
   _pd1= _array + addr((_k1)<<1);                                  \
   Y_MM_LOAD_PAIR_NEXT(__ai, __bi, _pd1);                          \
   prefetch_p(_pd1);                                               \
   Y_MM_MUL_PD(__ar, __ai, _f##r);                                 \
   Y_MM_MUL_PD(__ai, __ai, _f##i);                                 \
   Y_MM_MUL_PD(__br, __bi, _f##i);                                 \
   Y_MM_MUL_PD(__bi, __bi, _f##r);                                 \
   _pd0= _array + addr((_k1)<<1);                                  \
   Y_MM_LOAD_PAIR_NEXT(__cr, __ci, _pd0);                          \
   prefetch_p(_pd0);                                               \
   Y_MM_SUB_PD(__ar, __ar, __br);                                  \
   Y_MM_ADD_PD(__ai, __ai, __bi);                                  \
   Y_MM_ADD_PD(_t0##r, __cr, __ar);                                \
   Y_MM_SUB_PD(_t1##r, __cr, __ar);                                \
   Y_MM_ADD_PD(_t0##i, __ci, __ai);                                \
   Y_MM_SUB_PD(_t1##i, __ci, __ai);                                \
}
*/
/*
#define sse2_load_inter_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1)\
{                                                                \
  Y__M128D __x0,__x1,__x2;                                       \
  Y_MM_LOAD_INTER_PAIR_NEXT(__x0, __x1, _pd1 + _index);          \
  prefetch_p( _pd1 + _index);                                    \
  Y_MM_MOV_PD(_t1##r, __x0);                                     \
  Y_MM_MOV_PD(__x2, _f##r);                                      \
  Y_MM_MOV_PD(_t1##i, _f##i);                                    \
  Y_MM_MUL_PD(_t1##r, _t1##r, __x2);                             \
  Y_MM_MUL_PD(__x2, __x2, __x1);                                 \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);      \
  prefetch_p( _pd0 + _index);                                    \
  Y_MM_MUL_PD(__x1, __x1, _t1##i);                               \
  Y_MM_MUL_PD(_t1##i, _t1##i, __x0);                             \
  Y_MM_SUB_PD(_t1##r, _t1##r, __x1);                             \
  Y_MM_MOV_PD(__x0, _t0##r);                                     \
  Y_MM_ADD_PD(_t1##i, _t1##i, __x2);                             \
  Y_MM_MOV_PD(__x1, _t0##i);                                     \
  Y_MM_SUB_PD(__x0, __x0, _t1##r);                               \
  Y_MM_ADD_PD(_t0##r, _t0##r, _t1##r);                           \
  Y_MM_SUB_PD(__x1, __x1, _t1##i);                               \
  Y_MM_ADD_PD(_t0##i, _t0##i, _t1##i);                           \
  Y_MM_MOV_PD(_t1##r, __x0);                                     \
  Y_MM_MOV_PD(_t1##i, __x1);                                     \
}
*/

#define sse2_load_inter_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1)\
{                                                                \
  Y__M128D __ar, __ai;                                           \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);      \
  prefetch_p( _pd1 + _index);                                    \
  Y_MM_MUL_PD(__ar, _t1##r, _f##i);                              \
  Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                            \
  Y_MM_MUL_PD(__ai, _t1##i, _f##i);                              \
  Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                            \
  Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                             \
  Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                             \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);      \
  prefetch_p( _pd0 + _index);                                    \
  sse2_addsub( _t0, _t1);                                        \
}

#define sse2_load_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1)\
{                                                          \
  Y__M128D __ar, __ai;                                     \
  Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);      \
  prefetch_p( _pd1 + _index);                              \
  Y_MM_MUL_PD(__ar, _t1##r, _f##i);                        \
  Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                      \
  Y_MM_MUL_PD(__ai, _t1##i, _f##i);                        \
  Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                      \
  Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                       \
  Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                       \
  Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);      \
  prefetch_p( _pd0 + _index);                              \
  sse2_addsub( _t0, _t1);                                  \
}

#define sse2_load_inter_muladdsub_pp_no_fetch(_t0,_t1,_f,_index,_pd0,_pd1)\
{                                                                \
  Y__M128D __ar, __ai;                                           \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);      \
  Y_MM_MUL_PD(__ar, _t1##r, _f##i);                              \
  Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                            \
  Y_MM_MUL_PD(__ai, _t1##i, _f##i);                              \
  Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                            \
  Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                             \
  Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                             \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);      \
  sse2_addsub( _t0, _t1);                                        \
}

#define sse2_load_muladdsub_pp_no_fetch(_t0,_t1,_f,_index,_pd0,_pd1)\
{                                                          \
  Y__M128D __ar, __ai;                                     \
  Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);      \
  Y_MM_MUL_PD(__ar, _t1##r, _f##i);                        \
  Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                      \
  Y_MM_MUL_PD(__ai, _t1##i, _f##i);                        \
  Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                      \
  Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                       \
  Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                       \
  Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);      \
  sse2_addsub( _t0, _t1);                                  \
}

#define sse2_load_inter_muladdsub_no_next(_t0, _t1, _f, _ii, _jj, _pd0, _pd1) \
{                                                                             \
  Y__M128D __ar, __ai;                                                        \
  Y_MM_LOAD_INTER_PAIR (_t1##r, _t1##i, (_pd0 + _jj), (_pd1 + _jj));          \
  prefetch_p ( _pd1 + _jj);                                                   \
  prefetch_p ( _pd0 + _jj);                                                   \
  Y_MM_MUL_PD(__ar, _t1##r, _f##i);                                           \
  Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                                         \
  Y_MM_MUL_PD(__ai, _t1##i, _f##i);                                           \
  Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                                         \
  Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                                          \
  Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                                          \
  Y_MM_LOAD_INTER_PAIR (_t0##r, _t0##i, (_pd0 + _ii), (_pd1 + _ii));          \
  prefetch_p ( _pd1 + _ii);                                                   \
  prefetch_p ( _pd0 + _ii);                                                   \
  sse2_addsub( _t0, _t1);                                                     \
}

#define sse2_muladdsub_store_inter(_array,_k0,_k1,_t0,_t1,_f)         \
   sse2_mul0(_t1, _f);                                                \
   sse2_addsub(_t0, _t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
   Y_MM_STORE_INTER_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_muladdsub_store(_array,_k0,_k1,_t0,_t1,_f)         \
   sse2_mul0(_t1, _f);                                          \
   sse2_addsub(_t0, _t1);                                       \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k0<<1), _t0##r, _t0##i); \
   Y_MM_STORE_PAIR_NEXT(_array + addr(_k1<<1), _t1##r, _t1##i);


#define sse2_muladdsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1,_f) \
{                                                                 \
   Y__M128D __ar,__ai;                                            \
   Y_MM_MUL_PD(__ar, _t1##r, _f##i);                              \
   Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                            \
   Y_MM_MUL_PD(__ai, _t1##i, _f##i);                              \
   Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                            \
   Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                             \
   Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                             \
   sse2_addsub (_t0, _t1);                                        \
   Y_MM_STORE_INTER_PAIR_NEXT(_pd0 + _index, _t0##r, _t0##i);     \
   Y_MM_STORE_INTER_PAIR_NEXT(_pd1 + _index, _t1##r, _t1##i);     \
}

#define sse2_muladdsub_store_p(_index,_pd0,_pd1,_t0,_t1,_f) \
{                                                                 \
   Y__M128D __ar,__ai;                                            \
   Y_MM_MUL_PD(__ar, _t1##r, _f##i);                              \
   Y_MM_MUL_PD(_t1##r, _t1##r, _f##r);                            \
   Y_MM_MUL_PD(__ai, _t1##i, _f##i);                              \
   Y_MM_MUL_PD(_t1##i, _t1##i, _f##r);                            \
   Y_MM_SUB_PD(_t1##r, _t1##r, __ai);                             \
   Y_MM_ADD_PD(_t1##i, _t1##i, __ar);                             \
   sse2_addsub (_t0, _t1);                                        \
   Y_MM_STORE_PAIR_NEXT(_pd0 + _index, _t0##r, _t0##i);     \
   Y_MM_STORE_PAIR_NEXT(_pd1 + _index, _t1##r, _t1##i);     \
}

#endif /* ASM_MULADDSUB */

/**********************************************************************
   cplx_load_mulmuladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/

#if defined(ASM_MULMULADDSUB)

# define sse2_asm_mulmuladdsub(__t0, __t1, __f0, __f1) \
__asm__ volatile ("movapd %4, %%xmm4             \n"   \
"        movapd %0, %%xmm5                       \n"   \
"        movapd %1, %%xmm6                       \n"   \
"        movapd %2, %%xmm7                       \n"   \
"        mulpd %%xmm4, %0                        \n"   \
"        mulpd %%xmm4, %1                        \n"   \
"        mulpd %5, %%xmm5                        \n"   \
"        mulpd %5, %%xmm6                        \n"   \
"        mulpd %7, %%xmm7                        \n"   \
"        subpd %%xmm6, %0                        \n"   \
"        movapd  %3, %%xmm6                      \n"   \
"        mulpd  %6, %2                           \n"   \
"        mulpd  %7, %%xmm6                       \n"   \
"        addpd  %%xmm5, %1                       \n"   \
"        mulpd  %6, %3                           \n"   \
"        subpd  %2, %%xmm6                       \n"   \
"        movapd  %0, %2                          \n"   \
"        addpd  %3, %%xmm7                       \n"   \
"        movapd  %1, %3                          \n"   \
"        subpd  %%xmm6, %0                       \n"   \
"        addpd  %%xmm6, %2                       \n"   \
"        addpd  %%xmm7, %1                       \n"   \
"        subpd  %%xmm7, %3                       \n"   \
: "+x" (__t0##r), "+x" (__t0##i), "+x" (__t1##r), "+x" (__t1##i)  \
: "m" (__f0##r), "m" (__f0##i), "m" (__f1##r), "m" (__f1##i) \
: "xmm4", "xmm5", "xmm6", "xmm7" )

# define sse2_mulmuladdsub(_t0, _t1, _f0, _f1) \
  sse2_asm_mulmuladdsub (_t0, _t1, _f0, _f1);


# define sse2_load_inter_mulmuladdsub(_t0, _t1, _f0, _f1, _array, _k0, _k1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));         \
   prefetch_p(&_array[addr(_k0<<1)]);                                       \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));         \
   prefetch_p(&_array[addr(_k1<<1)]);                                       \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_mulmuladdsub(_t0, _t1, _f0, _f1, _array, _k0, _k1) \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));         \
   prefetch_p(&_array[addr(_k0<<1)]);                                 \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));         \
   prefetch_p(&_array[addr(_k1<<1)]);                                 \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_inter_mulmuladdsub_pp(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
  Y_MM_LOAD_INTER_PAIR_NEXT (_t0##r, _t0##i, _pd0 + _index);                     \
  prefetch_p( _pd0 + _index);                                                    \
  Y_MM_LOAD_INTER_PAIR_NEXT (_t1##r, _t1##i, _pd1 + _index);                     \
  prefetch_p( _pd1 + _index);                                                    \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_mulmuladdsub_pp(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
  Y_MM_LOAD_PAIR_NEXT (_t0##r, _t0##i, _pd0 + _index);                     \
  prefetch_p( _pd0 + _index);                                              \
  Y_MM_LOAD_PAIR_NEXT (_t1##r, _t1##i, _pd1 + _index);                     \
  prefetch_p( _pd1 + _index);                                              \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_inter_mulmuladdsub_pp_no_fetch(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
  Y_MM_LOAD_INTER_PAIR_NEXT (_t0##r, _t0##i, _pd0 + _index);                              \
  Y_MM_LOAD_INTER_PAIR_NEXT (_t1##r, _t1##i, _pd1 + _index);                              \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_mulmuladdsub_pp_no_fetch(_t0, _t1, _f0, _f1, _index, _pd0, _pd1) \
  Y_MM_LOAD_PAIR_NEXT (_t0##r, _t0##i, _pd0 + _index);                              \
  Y_MM_LOAD_PAIR_NEXT (_t1##r, _t1##i, _pd1 + _index);                              \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

# define sse2_load_inter_mulmuladdsub_no_next(_t0, _t1, _f0, _f1, _ii, _jj, _pd0, _pd1)\
  Y_MM_LOAD_INTER_PAIR (_t0##r, _t0##i, (_pd0 + _ii), (_pd1 + _ii));                   \
  prefetch_p ( _pd1 + _ii);                                                            \
  prefetch_p ( _pd0 + _ii);                                                            \
  Y_MM_LOAD_INTER_PAIR (_t1##r, _t1##i, (_pd0 + _jj), (_pd1 + _jj));                   \
  prefetch_p ( _pd1 + _jj);                                                            \
  prefetch_p ( _pd0 + _jj);                                                            \
  sse2_asm_mulmuladdsub ( _t0, _t1, _f0, _f1);

#else /* ASM_MULMULADDSUB */

# define sse2_mulmuladdsub(_t0,_t1,_f0,_f1) \
   sse2_mul0(_t0, _f0);                     \
   sse2_mul1(_t1, _f1);                     \
   sse2_addsub(_t0, _t1);

/*
# define sse2_mulmuladdsub(_t0,_t1,_f0,_f1) \
{                                           \
   Y__M128D __a,__b,__c,__d,__e,__f;        \
   Y_MM_MUL_PD( __a, _t0##r, _f0##r);       \
   Y_MM_MUL_PD( __b, _t0##i, _f0##i);       \
   Y_MM_MUL_PD( __c, _t1##r, _f1##r);       \
   Y_MM_MUL_PD( __d, _t1##i, _f1##i);       \
   Y_MM_SUB_PD( __a, __a, __b);             \
   Y_MM_MUL_PD( __e, _t0##i, _f0##r);       \
   Y_MM_MUL_PD( __f, _t0##r, _f0##i);       \
   Y_MM_SUB_PD( __c, __c, __d);             \
   Y_MM_MUL_PD( __b, _t1##i, _f1##r);       \
   Y_MM_MUL_PD( __d, _t1##r, _f1##i);       \
   Y_MM_ADD_PD( __e, __e, __f);             \
   Y_MM_ADD_PD( _t0##r, __a, __c);          \
   Y_MM_ADD_PD( __b, __b, __d);             \
   Y_MM_SUB_PD( _t1##r, __a, __c);          \
   Y_MM_ADD_PD( _t0##i, __e, __b);          \
   Y_MM_SUB_PD( _t1##i, __e, __b);          \
}
*/

#define sse2_load_inter_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1) \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));  \
   prefetch_p(&_array[addr(_k0<<1)]);                                \
   sse2_mul0(_t0, _f0);                                              \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));  \
   prefetch_p(&_array[addr(_k1<<1)]);                                \
   sse2_mul0(_t1, _f1);                                              \
   sse2_addsub (_t0, _t1);

#define sse2_load_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1) \
   Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, array + addr(_k0<<1));  \
   prefetch_p(&_array[addr(_k0<<1)]);                          \
   sse2_mul0(_t0, _f0);                                        \
   Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, array + addr(_k1<<1));  \
   prefetch_p(&_array[addr(_k1<<1)]);                          \
   sse2_mul0(_t1, _f1);                                        \
   sse2_addsub(_t0, _t1);


# define sse2_load_inter_mulmuladdsub_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                \
  prefetch_p(_pd0 + _index);                                               \
  sse2_mul0 (_t0, _f0);                                                    \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                \
  prefetch_p(_pd1 + _index);                                               \
  sse2_mul0 (_t1, _f1);                                                    \
  sse2_addsub (_t0, _t1);


# define sse2_load_mulmuladdsub_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
  Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                \
  prefetch_p(_pd0 + _index);                                         \
  sse2_mul0 (_t0, _f0);                                              \
  Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                \
  prefetch_p(_pd1 + _index);                                         \
  sse2_mul0 (_t1, _f1);                                              \
  sse2_addsub (_t0, _t1);

# define sse2_load_inter_mulmuladdsub_pp_no_fetch(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                \
  sse2_mul0 (_t0, _f0);                                                    \
  Y_MM_LOAD_INTER_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                \
  sse2_mul0 (_t1, _f1);                                                    \
  sse2_addsub (_t0, _t1);


# define sse2_load_mulmuladdsub_pp_no_fetch(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
  Y_MM_LOAD_PAIR_NEXT(_t0##r, _t0##i, _pd0 + _index);                \
  sse2_mul0 (_t0, _f0);                                              \
  Y_MM_LOAD_PAIR_NEXT(_t1##r, _t1##i, _pd1 + _index);                \
  sse2_mul0 (_t1, _f1);                                              \
  sse2_addsub (_t0, _t1);

# define sse2_load_inter_mulmuladdsub_no_next(_t0, _t1, _f0, _f1, _ii, _jj, _pd0, _pd1)\
  Y_MM_LOAD_INTER_PAIR (_t0##r, _t0##i, (_pd0 + _ii), (_pd1 + _ii));                   \
  prefetch_p ( _pd1 + _ii);                                                            \
  prefetch_p ( _pd0 + _ii);                                                            \
  sse2_mul0 (_t0, _f0);                                                                \
  Y_MM_LOAD_INTER_PAIR (_t1##r, _t1##i, (_pd0 + _jj), (_pd1 + _jj));                   \
  prefetch_p ( _pd1 + _jj);                                                            \
  prefetch_p ( _pd0 + _jj);                                                            \
  sse2_mul0 (_t1, _f1);                                                                \
  sse2_addsub (_t0, _t1);

#endif /* ASM_MULMULADDSUB */


#define sse2_mul( _t, _s0, _s1)          \
{                                        \
   Y__M128D __a,__b,__c,__d;             \
   Y_MM_MUL_PD( __a, _s0##r, _s1##r);    \
   Y_MM_MUL_PD( __b, _s0##i, _s1##i);    \
   Y_MM_MUL_PD( __c, _s0##i, _s1##r);    \
   Y_MM_MUL_PD( __d, _s0##r, _s1##i);    \
   Y_MM_SUB_PD( _t##r, __a, __b );       \
   Y_MM_ADD_PD( _t##i, __c, __d );       \
}


#define sse2_load_inter_mul_p(_t,_pd,_f,_array,_k)               \
{                                                                \
   Y__M128D __a, __b;                                            \
   _pd = _array + addr((_k)<<1);                                 \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t##r, _t##i, _pd);                 \
   Y_MM_MUL_PD( __a, _t##i, _f##i);                              \
   Y_MM_MUL_PD( __b, _t##r, _f##i);                              \
   Y_MM_MUL_PD( _t##r, _t##r, _f##r);                            \
   Y_MM_MUL_PD( _t##i, _t##i, _f##r);                            \
   Y_MM_SUB_PD( _t##r, _t##r, __a);                              \
   Y_MM_ADD_PD( _t##i, _t##i, __b);                              \
}

#define sse2_load_mul_p(_t,_pd,_f,_array,_k)               \
{                                                          \
   Y__M128D __a, __b;                                      \
   _pd = _array + addr((_k)<<1);                           \
   Y_MM_LOAD_PAIR_NEXT(_t##r, _t##i, _pd);                 \
   Y_MM_MUL_PD( __a, _t##i, _f##i);                        \
   Y_MM_MUL_PD( __b, _t##r, _f##i);                        \
   Y_MM_MUL_PD( _t##r, _t##r, _f##r);                      \
   Y_MM_MUL_PD( _t##i, _t##i, _f##r);                      \
   Y_MM_SUB_PD( _t##r, _t##r, __a);                        \
   Y_MM_ADD_PD( _t##i, _t##i, __b);                        \
}


#define sse2_load_mul_pp(_t,_f,_index,_pd)                 \
{                                                          \
   Y__M128D __a, __b;                                      \
   Y_MM_LOAD_PAIR_NEXT(_t##r, _t##i, _pd + _index);        \
   prefetch_p( _pd + _index);                              \
   Y_MM_MUL_PD( __a, _t##i, _f##i);                        \
   Y_MM_MUL_PD( __b, _t##r, _f##i);                        \
   Y_MM_MUL_PD( _t##r, _t##r, _f##r);                      \
   Y_MM_MUL_PD( _t##i, _t##i, _f##r);                      \
   Y_MM_SUB_PD( _t##r, _t##r, __a);                        \
   Y_MM_ADD_PD( _t##i, _t##i, __b);                        \
}

#define sse2_load_inter_mul_pp(_t,_f,_index,_pd)                 \
{                                                                \
   Y__M128D __a, __b;                                            \
   Y_MM_LOAD_INTER_PAIR_NEXT(_t##r, _t##i, _pd + _index);        \
   prefetch_p( _pd + _index);                                    \
   Y_MM_MUL_PD( __a, _t##i, _f##i);                              \
   Y_MM_MUL_PD( __b, _t##r, _f##i);                              \
   Y_MM_MUL_PD( _t##r, _t##r, _f##r);                            \
   Y_MM_MUL_PD( _t##i, _t##i, _f##r);                            \
   Y_MM_SUB_PD( _t##r, _t##r, __a);                              \
   Y_MM_ADD_PD( _t##i, _t##i, __b);                              \
}


#define sse2_div( _t, _s0, _s1)          \
{                                        \
   Y__M128D __a,__b,__c,__d;             \
   Y_MM_MUL_PD( __a, _s0##r, _s1##r);    \
   Y_MM_MUL_PD( __b, _s0##i, _s1##i);    \
   Y_MM_MUL_PD( __c, _s0##i, _s1##r);    \
   Y_MM_MUL_PD( __d, _s0##r, _s1##i);    \
   Y_MM_ADD_PD( _t##r, __a, __b );       \
   Y_MM_SUB_PD( _t##i, __c, __d );       \
}

#define sse2_square( _t, _s)          \
{                                     \
    Y__M128D __a, __b;                \
    Y_MM_MUL_PD( __a, _s##r, _s##r);  \
    Y_MM_MUL_PD( __b, _s##i, _s##i);  \
    Y_MM_MUL_PD( _t##i, _s##i, _s##r);\
    Y_MM_SUB_PD( _t##r, __a, __b);    \
    Y_MM_ADD_PD( _t##i, _t##i, _t##i);\
}

/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#if defined(OLD_MACROS)
# define sse2_mul_1_4_F_addsub(_t0,_t1)  \
{                                        \
   Y__M128D __ar, __ai, __a;             \
   Y_MM_MOV_PD( __ar, _t0##r);           \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##i); \
   Y_MM_MOV_PD( __ai, _t0##i);           \
   Y_MM_SUB_PD( _t0##i, _t0##i, _t1##r); \
   Y_MM_MOV_PD( __a, _t1##r);            \
   Y_MM_SUB_PD( _t1##r, __ar, _t1##i);   \
   Y_MM_ADD_PD( _t1##i, __ai, __a);      \
}

#define sse2_mul_1_4_F_addsub_store_inter(_d,_k0,_k1,_t0,_t1) \
{                                                             \
   Y__M128D __c, __d;                                         \
   Y_MM_SUB_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##i );                     \
   Y_MM_SUB_PD( _t0##i, _t0##i, _t1##r );                     \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k1<<1), __c, __d);  \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k0<<1), _t0##r, _t0##i);  \
}

#define sse2_mul_1_4_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __c, __d;                                         \
   Y_MM_SUB_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##i );                     \
   Y_MM_SUB_PD( _t0##i, _t0##i, _t1##r );                     \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, __c, __d);      \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
}

#define sse2_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b, __c, __d;                               \
   Y_MM_ADD_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_SUB_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_SUB_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, __a, __b);            \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, __c, __d);            \
}
#else


#define sse2_mul_1_4_F_addsub(_t0,_t1)   \
   sse2_mul_1_4_F(_t1);                  \
   sse2_addsub (_t0, _t1);

#define sse2_mul_1_4_F_addsub_store_inter(_d,_k0,_k1,_t0,_t1)       \
   sse2_mul_1_4_F(_t1);                                             \
   sse2_addsub (_t0, _t1);                                          \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k0<<1), _t0##r, _t0##i);  \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_mul_1_4_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_4_F(_t1);                                               \
   sse2_addsub (_t0, _t1);                                            \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_4_F(_t1);                                         \
   sse2_addsub (_t0, _t1);                                      \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

/*
#define sse2_mul_1_4_F_addsub(_t0,_t1)   \
{                                        \
   Y__M128D __x0, __x1;                  \
   Y_MM_MOV_PD( __x0, _t1##r);           \
   Y_MM_MOV_PD( __x1, _t1##i);           \
   Y_MM_MOV_PD( _t1##r, _t0##r);         \
   Y_MM_MOV_PD( _t1##i, _t0##i);         \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x1);   \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x1);   \
   Y_MM_ADD_PD( _t1##i, _t1##i, __x0);   \
   Y_MM_SUB_PD( _t0##i, _t0##i, __x0);   \
}


#define sse2_mul_1_4_F_addsub_store_inter(_d,_k0,_k1,_t0,_t1)       \
{                                                                   \
   Y__M128D __x0, __x1;                                             \
   Y_MM_MOV_PD( __x0, _t1##r);                                      \
   Y_MM_MOV_PD( __x1, _t1##i);                                      \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                    \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                    \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x1);                              \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x1);                              \
   Y_MM_ADD_PD( _t1##i, _t1##i, __x0);                              \
   Y_MM_SUB_PD( _t0##i, _t0##i, __x0);                              \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k0<<1), _t0##r, _t0##i);  \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k1<<1), _t1##r, _t1##i);  \
}

#define sse2_mul_1_4_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b;                                         \
   Y_MM_SUB_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##i );                     \
   Y_MM_SUB_PD( _t0##i, _t0##i, _t1##r );                     \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, __a, __b);      \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);\
}

#define sse2_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b;                                         \
   Y_MM_SUB_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( _t0##r, _t0##r, _t1##i );                     \
   Y_MM_SUB_PD( _t0##i, _t0##i, _t1##r );                     \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, __a, __b);            \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);      \
}
*/

#endif
/***********************************************************************
   cplx_mul_1_4_B macro:
		_t = _t1 * G^(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#if defined(OLD_MACROS)
# define sse2_mul_1_4_B_addsub(_t0,_t1)  \
{                                        \
   Y__M128D __ar, __ai, __a;             \
   Y_MM_MOV_PD( __ar, _t0##r);           \
   Y_MM_SUB_PD( _t0##r, _t0##r, _t1##i); \
   Y_MM_MOV_PD( __ai, _t0##i);           \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##r); \
   Y_MM_MOV_PD( __a, _t1##r);            \
   Y_MM_ADD_PD( _t1##r, __ar, _t1##i);   \
   Y_MM_SUB_PD( _t1##i, __ai, __a);      \
}

#define sse2_mul_1_4_B_addsub_store_inter(_d,_k0,_k1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b, __c, __d;                               \
   Y_MM_SUB_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_SUB_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr (_k0<<1), __a, __b); \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr (_k1<<1), __c, __d); \
}


#define sse2_mul_1_4_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b, __c, __d;                               \
   Y_MM_SUB_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_SUB_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, __a, __b);      \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, __c, __d);      \
}

#define sse2_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                             \
   Y__M128D __a, __b, __c, __d;                               \
   Y_MM_SUB_PD( __a, _t0##r, _t1##i );                        \
   Y_MM_ADD_PD( __b, _t0##i, _t1##r );                        \
   Y_MM_ADD_PD( __c, _t0##r, _t1##i );                        \
   Y_MM_SUB_PD( __d, _t0##i, _t1##r );                        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, __a, __b);            \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, __c, __d);            \
}

#else


#define sse2_mul_1_4_B_addsub(_t0,_t1)   \
   sse2_mul_1_4_B(_t1);                  \
   sse2_addsub (_t0, _t1);

#define sse2_mul_1_4_B_addsub_store_inter(_d,_k0,_k1,_t0,_t1)       \
   sse2_mul_1_4_B(_t1);                                             \
   sse2_addsub (_t0, _t1);                                          \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k0<<1), _t0##r, _t0##i);  \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k1<<1), _t1##r, _t1##i);

#define sse2_mul_1_4_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_4_B(_t1);                                               \
   sse2_addsub (_t0, _t1);                                            \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_4_B(_t1);                                         \
   sse2_addsub (_t0, _t1);                                      \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

/*
#define sse2_mul_1_4_B_addsub(_t0,_t1)   \
{                                        \
   Y__M128D __x0, __x1;                  \
   Y_MM_MOV_PD( __x0, _t1##r);           \
   Y_MM_MOV_PD( __x1, _t1##i);           \
   Y_MM_MOV_PD( _t1##r, _t0##r);         \
   Y_MM_MOV_PD( _t1##i, _t0##i);         \
   Y_MM_ADD_PD( _t1##r, _t1##r, __x1);   \
   Y_MM_SUB_PD( _t0##r, _t0##r, __x1);   \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x0);   \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x0);   \
}


#define sse2_mul_1_4_B_addsub_store_inter(_d,_k0,_k1,_t0,_t1)       \
{                                                                   \
   Y__M128D __x0, __x1;                                             \
   Y_MM_MOV_PD( __x0, _t1##r);                                      \
   Y_MM_MOV_PD( __x1, _t1##i);                                      \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                    \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                    \
   Y_MM_ADD_PD( _t1##r, _t1##r, __x1);                              \
   Y_MM_SUB_PD( _t0##r, _t0##r, __x1);                              \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x0);                              \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x0);                              \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k0<<1), _t0##r, _t0##i);  \
   Y_MM_STORE_INTER_PAIR_NEXT( _d + addr(_k1<<1), _t1##r, _t1##i);  \
}

#define sse2_mul_1_4_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_ADD_PD( __a, _t0##r, _t1##i );                                \
   Y_MM_SUB_PD( __b, _t0##i, _t1##r );                                \
   Y_MM_SUB_PD( _t0##r, _t0##r, _t1##i );                             \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##r );                             \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, __a, __b);              \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_ADD_PD( __a, _t0##r, _t1##i );                          \
   Y_MM_SUB_PD( __b, _t0##i, _t1##r );                          \
   Y_MM_SUB_PD( _t0##r, _t0##r, _t1##i );                       \
   Y_MM_ADD_PD( _t0##i, _t0##i, _t1##r );                       \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, __a, __b);              \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}
*/
#endif

#define sse2_mul_1_4_B_addsub_store_inter_no_next(_pd, _pu, _i, _j, _t0, _t1)\
{                                                                            \
   Y__M128D __x0, __x1;                                                      \
   Y_MM_MOV_PD( __x0, _t1##r);                                               \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                             \
   Y_MM_MOV_PD( __x1, _t1##i);                                               \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                             \
   Y_MM_ADD_PD( _t1##r, _t1##r, __x1);                                       \
   Y_MM_SUB_PD( _t0##r, _t0##r, __x1);                                       \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x0);                                       \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x0);                                       \
   Y_MM_STORE_INTER_PAIR( _pd + _j, _pu + _j, _t1##r, _t1##i);               \
   Y_MM_STORE_INTER_PAIR( _pd + _i, _pu + _i, _t0##r, _t0##i);               \
}

/***********************************************************************
   cplx_mul_1_8_F macro:
		_t = _t1 * G^-(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/

#if defined(OLD_MACROS)
# define sse2_mul_1_8_F_addsub(_t0,_t1)\
{                                      \
   Y__M128D __a, __b;                  \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i ); \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r ); \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r); \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r); \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);  \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);  \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);  \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);  \
}
#define sse2_mul_1_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i );                                \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r );                                \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );                               \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );                               \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                 \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                 \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                 \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                 \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i );                          \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r );                          \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );                         \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );                         \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                           \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                           \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                           \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                           \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#else

# define sse2_mul_1_8_F_addsub(_t0,_t1)    \
   sse2_mul_1_8_F( _t1);                   \
   sse2_addsub(_t0, _t1);

#define sse2_mul_1_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_8_F( _t1);                                              \
   sse2_addsub(_t0, _t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_8_F( _t1);                                        \
   sse2_addsub(_t0, _t1);                                       \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

/*
# define sse2_mul_1_8_F_addsub(_t0,_t1)    \
{                                          \
   Y__M128D __x1, __x0;                    \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);\
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);\
   Y_MM_MOV_PD( __x0, _t1##r);             \
   Y_MM_MOV_PD( __x1, _t1##i);             \
   Y_MM_MOV_PD( _t1##i, _t0##i);           \
   Y_MM_ADD_PD( __x0, __x0, __x1);         \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );      \
   Y_MM_MOV_PD( _t1##r, _t0##r);           \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);     \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);     \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);     \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);     \
}

#define sse2_mul_1_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __x1, __x0;                                               \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);                           \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);                           \
   Y_MM_MOV_PD( __x0, _t1##r);                                        \
   Y_MM_MOV_PD( __x1, _t1##i);                                        \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                      \
   Y_MM_ADD_PD( __x0, __x0, __x1);                                    \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );                                 \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                      \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                                \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                                \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __x1, __x0;                                         \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);                     \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);                     \
   Y_MM_MOV_PD( __x0, _t1##r);                                  \
   Y_MM_MOV_PD( __x1, _t1##i);                                  \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                \
   Y_MM_ADD_PD( __x0, __x0, __x1);                              \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );                           \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                          \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                          \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}
*/
#endif


/***********************************************************************
   cplx_mul_1_8_B macro:
		_t = _t1 * G^(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#if defined(OLD_MACROS)
# define sse2_mul_1_8_B_addsub(_t0,_t1)\
{                                      \
   Y__M128D __a, __b;                  \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i ); \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r ); \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );\
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );\
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);  \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);  \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);  \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);  \
}

#define sse2_mul_1_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i );                                \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r );                                \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );                               \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );                               \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                 \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                 \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                 \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                 \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i );                          \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r );                          \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8r );                         \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8r );                         \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                           \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                           \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                           \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                           \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#else


# define sse2_mul_1_8_B_addsub(_t0,_t1)    \
   sse2_mul_1_8_B( _t1);                   \
   sse2_addsub(_t0, _t1);

#define sse2_mul_1_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_8_B( _t1);                                              \
   sse2_addsub(_t0, _t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_1_8_B( _t1);                                        \
   sse2_addsub(_t0, _t1);                                       \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

/*
# define sse2_mul_1_8_B_addsub(_t0,_t1)    \
{                                          \
   Y__M128D __x1, __x0;                    \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);\
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);\
   Y_MM_MOV_PD( __x0, _t1##r);             \
   Y_MM_MOV_PD( __x1, _t1##i);             \
   Y_MM_MOV_PD( _t1##i, _t0##i);           \
   Y_MM_SUB_PD( __x0, __x0, __x1);         \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );      \
   Y_MM_MOV_PD( _t1##r, _t0##r);           \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);     \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);     \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);     \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);     \
}

#define sse2_mul_1_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __x1, __x0;                                               \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);                           \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);                           \
   Y_MM_MOV_PD( __x0, _t1##r);                                        \
   Y_MM_MOV_PD( __x1, _t1##i);                                        \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                      \
   Y_MM_SUB_PD( __x0, __x0, __x1);                                    \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );                                 \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                      \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                                \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                                \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __x1, __x0;                                         \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8r);                     \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8r);                     \
   Y_MM_MOV_PD( __x0, _t1##r);                                  \
   Y_MM_MOV_PD( __x1, _t1##i);                                  \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                \
   Y_MM_SUB_PD( __x0, __x0, __x1);                              \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );                           \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                          \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                          \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}
*/
#endif

/***********************************************************************
   cplx_mul_3_8_F macro:
		_t = _t1 * G^-(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#if defined(OLD_MACROS)
# define sse2_mul_3_8_F_addsub(_t0,_t1)\
{                                      \
   Y__M128D __a, __b;                  \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i ); \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r ); \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i); \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i); \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);  \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);  \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);  \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);  \
}


#define sse2_mul_3_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i );                                \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r );                                \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );                               \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );                               \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                 \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                 \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                 \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                 \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_SUB_PD( __a , _t1##r,_t1##i );                          \
   Y_MM_ADD_PD( __b , _t1##i,_t1##r );                          \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );                         \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );                         \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                           \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                           \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                           \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                           \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#else

# define sse2_mul_3_8_F_addsub(_t0,_t1)    \
   sse2_mul_3_8_F( _t1);                   \
   sse2_addsub(_t0, _t1);

#define sse2_mul_3_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_3_8_F( _t1);                                              \
   sse2_addsub(_t0, _t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_3_8_F( _t1);                                        \
   sse2_addsub(_t0, _t1);                                       \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

/*
# define sse2_mul_3_8_F_addsub(_t0,_t1)    \
{                                          \
   Y__M128D __x1, __x0;                    \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);\
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);\
   Y_MM_MOV_PD( __x0, _t1##r);             \
   Y_MM_MOV_PD( __x1, _t1##i);             \
   Y_MM_MOV_PD( _t1##i, _t0##i);           \
   Y_MM_SUB_PD( __x0, __x0, __x1);         \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );      \
   Y_MM_MOV_PD( _t1##r, _t0##r);           \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);     \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);     \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);     \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);     \
}

#define sse2_mul_3_8_F_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __x1, __x0;                                               \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);                           \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);                           \
   Y_MM_MOV_PD( __x0, _t1##r);                                        \
   Y_MM_MOV_PD( __x1, _t1##i);                                        \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                      \
   Y_MM_SUB_PD( __x0, __x0, __x1);                                    \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );                                 \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                      \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                                \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                                \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __x1, __x0;                                         \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);                     \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);                     \
   Y_MM_MOV_PD( __x0, _t1##r);                                  \
   Y_MM_MOV_PD( __x1, _t1##i);                                  \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                \
   Y_MM_SUB_PD( __x0, __x0, __x1);                              \
   Y_MM_ADD_PD( __x1, __x1, _t1##r );                           \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                          \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                          \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}
*/
#endif
/***********************************************************************
   cplx_mul_3_8_B macro:
		_t = _t1 * G^(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#if defined(OLD_MACROS)
# define sse2_mul_3_8_B_addsub(_t0,_t1)\
{                                      \
   Y__M128D __a, __b;                  \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i ); \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r ); \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );\
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );\
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);  \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);  \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);  \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);  \
}

#define sse2_mul_3_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __a, __b;                                                 \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i );                                \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r );                                \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );                               \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );                               \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                                 \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                                 \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                                 \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                                 \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __a, __b;                                           \
   Y_MM_ADD_PD( __a , _t1##r,_t1##i );                          \
   Y_MM_SUB_PD( __b , _t1##i,_t1##r );                          \
   Y_MM_MUL_PD( __a , __a, MM_F_1_8i );                         \
   Y_MM_MUL_PD( __b , __b, MM_F_1_8i );                         \
   Y_MM_SUB_PD( _t1##r, _t0##r, __a);                           \
   Y_MM_ADD_PD( _t0##r, _t0##r, __a);                           \
   Y_MM_SUB_PD( _t1##i, _t0##i, __b);                           \
   Y_MM_ADD_PD( _t0##i, _t0##i, __b);                           \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#else

/*
# define sse2_mul_3_8_B_addsub(_t0,_t1)    \
   sse2_mul_3_8_B( _t1);                   \
   sse2_addsub(_t0, _t1);

#define sse2_mul_3_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_3_8_B( _t1);                                              \
   sse2_addsub(_t0, _t1);                                             \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

#define sse2_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
   sse2_mul_3_8_B( _t1);                                        \
   sse2_addsub(_t0, _t1);                                       \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);

*/
# define sse2_mul_3_8_B_addsub(_t0,_t1)    \
{                                          \
   Y__M128D __x1, __x0;                    \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);\
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);\
   Y_MM_MOV_PD( __x0, _t1##r);             \
   Y_MM_MOV_PD( __x1, _t1##i);             \
   Y_MM_MOV_PD( _t1##i, _t0##i);           \
   Y_MM_ADD_PD( __x0, __x0, __x1);         \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );      \
   Y_MM_MOV_PD( _t1##r, _t0##r);           \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);     \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);     \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);     \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);     \
}

#define sse2_mul_3_8_B_addsub_store_inter_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                                     \
   Y__M128D __x1, __x0;                                               \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);                           \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);                           \
   Y_MM_MOV_PD( __x0, _t1##r);                                        \
   Y_MM_MOV_PD( __x1, _t1##i);                                        \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                      \
   Y_MM_ADD_PD( __x0, __x0, __x1);                                    \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );                                 \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                      \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                                \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                                \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                                \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                                \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#define sse2_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
{                                                               \
   Y__M128D __x1, __x0;                                         \
   Y_MM_MUL_PD( _t1##r, _t1##r, MM_F_1_8i);                     \
   Y_MM_MUL_PD( _t1##i, _t1##i, MM_F_1_8i);                     \
   Y_MM_MOV_PD( __x0, _t1##r);                                  \
   Y_MM_MOV_PD( __x1, _t1##i);                                  \
   Y_MM_MOV_PD( _t1##i, _t0##i);                                \
   Y_MM_ADD_PD( __x0, __x0, __x1);                              \
   Y_MM_SUB_PD( __x1, __x1, _t1##r );                           \
   Y_MM_MOV_PD( _t1##r, _t0##r);                                \
   Y_MM_SUB_PD( _t1##i, _t1##i, __x1);                          \
   Y_MM_ADD_PD( _t0##i, _t0##i, __x1);                          \
   Y_MM_SUB_PD( _t1##r, _t1##r, __x0);                          \
   Y_MM_ADD_PD( _t0##r, _t0##r, __x0);                          \
   Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _t1##r, _t1##i);        \
   Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _t0##r, _t0##i);        \
}

#endif

/*********************************************************************
This is the dyadic square for nested-compls representation
 xk =(xk + xmk* ) * (xk + xmk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(-k) * (xk - xmk* ) * (xk - xmk* );
 xkm =(xkm + xk* ) * (xmk + xk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(k) * (xmk - xk* ) * (xmk - xk* );
 where all are pseudo-complex vars. 
**********************************************************************/
/* 18 fadds  12  Fmuls
#define square_nested(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= (((gkr)+ 1.0) * _r2 - (gki) * _r1);\
   _r1= (((gkr)+ 1.0) * _r1 + (gki) * _r2);\
   xk##r -= _r0;\
   xmk##r -= _r0;\
   xk##i -= _r1;\
   xmk##i += _r1;\
}
*/

/*
   On calling, _xkr and _xmkr are the first time variables when 
   calling in legacy mode, _xki and _xmki are the arguments when
   calling a second time in legacy mode 
 
   Once interlaced, *.d[0] will correspond to first call and *.d[1]
   to second.
 
  _gkr y _gki have to be interleaved, i.e, _gkr.d[0] and _gki.d[0] 
  are the gkr and gki on first call in legacy mode
  _gkr.d[1] and _gki.d[1] are gkr and gki on second call
 
  */
/*
#define square_nested(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= (((gkr)+ 1.0) * _r2 - (gki) * _r1);\
   _r1= (((gkr)+ 1.0) * _r1 + (gki) * _r2);\
   xk##r -= _r0;\
   xmk##r -= _r0;\
   xk##i -= _r1;\
   xmk##i += _r1;\
}
*/

#define sse2_square_nested(_xkr, _xmkr, _xki, _xmki, _gkr, _gki)    \
{                                                                   \
   Y__M128D __x0, __x1, __x2, __x3, __x4;                           \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #0          */             \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #4          */             \
   Y_MM_MOV_PD(__x0, _xkr);           /* #5          */             \
   Y_MM_MOV_PD(__x1, _xki);           /* #6          */             \
   Y_MM_SUB_PD(__x0, __x0, _xmkr);    /* #8-11       */             \
   Y_MM_MOV_PD(__x2, _xki);           /* #8-13       */             \
   Y_MM_MOV_PD(__x3, _xkr);           /* #9-14       */             \
   Y_MM_ADD_PD(__x1, __x1, _xmki);    /* #10-13      */             \
   Y_MM_MUL_PD(__x0, __x0, MM_0_5);   /* #12-17   __x0 = _r0   */   \
   Y_MM_SUB_PD(_xkr, _xkr, _xki);     /* #14-19      */             \
   Y_MM_MUL_PD(__x1, __x1, MM_0_5);   /* #16-19   __x1 = _r1   */   \
   Y_MM_MUL_PD(_xki, _xki, __x3);     /* #18-23    xki *= xkr  */   \
   Y_MM_MOV_PD(__x4, __x0);           /* #18-23      */             \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #20-23      */             \
   Y_MM_MOV_PD(__x2, __x0);           /* #21-26      */             \
   Y_MM_SUB_PD(__x4, __x4, __x1);     /* #22-25      */             \
   Y_MM_ADD_PD(__x2, __x2, __x1);     /* #24-27      */             \
   Y_MM_MUL_PD(_xkr, _xkr, __x3);     /* #26-31   _xkr = _r2   */   \
   Y_MM_MOV_PD(__x3, _xmkr);          /* #27-32      */             \
   Y_MM_MUL_PD(__x1, __x1, __x0);     /* #28-33   _r1 *= _r0   */   \
   Y_MM_MOV_PD(__x0, _xmki);          /* #29-34      */             \
   Y_MM_MUL_PD(__x2, __x2, __x4);     /* #30-35   _r2 = __x2   */   \
   Y_MM_MOV_PD(__x4, _gkr);           /* #31-36   ** */             \
   Y_MM_SUB_PD(_xmkr, _xmkr, _xmki);  /* #32-34      */             \
   Y_MM_MUL_PD(_xmki, _xmki, __x3);   /* #34-39  xmki *= xmkr  */   \
   Y_MM_ADD_PD(__x3, __x3, __x0);     /* #36-39      */             \
   Y_MM_ADD_PD(__x1, __x1, __x1);     /* #38-41   _r1 += r1   */    \
   Y_MM_ADD_PD(__x4, __x4, MM_1_0);   /* #40-43   ** */             \
   Y_MM_MOV_PD(__x0, __x2);           /* #40-45      */             \
   Y_MM_MUL_PD(_xmkr, _xmkr, __x3);   /* #42-47   xmkr = _r3 = .. */\
   Y_MM_MOV_PD(__x3, __x1);           /* #42-47      */             \
   Y_MM_ADD_PD(_xki, _xki, _xki);     /* #44-47   xki += kxi   */   \
   Y_MM_MUL_PD(__x1, __x1, _gki);     /* #46-51   ** */             \
   Y_MM_ADD_PD(_xmki, _xmki, _xmki);  /* #        xmki += xmki   */ \
   Y_MM_MUL_PD(__x0, __x0, __x4);     /* #48-53   _r2*(gkr + 1)   */\
   Y_MM_MUL_PD(__x2, __x2, _gki);     /* #50-55   ** */             \
   Y_MM_MUL_PD(__x3, __x3, __x4);     /* #52-57   _r1*(gkr + 1)   */\
   Y_MM_SUB_PD(__x0, __x0, __x1);     /* #54-57   ** _r0 */         \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #58-61   ** _r1 */         \
   Y_MM_SUB_PD(_xkr, _xkr, __x0);     /* #60-63   ** */             \
   Y_MM_SUB_PD(_xmkr, _xmkr, __x0);   /* #62-65   ** */             \
   Y_MM_SUB_PD(_xki, _xki, __x3);     /* #64-67   ** */             \
   Y_MM_ADD_PD(_xmki, _xmki, __x3);   /* #66-69   ** */             \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #67         */             \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #71         */             \
}

/*
#define square_nested_eq(xk,gkr) \
{\
   y_limb_t _y0r=xk##r,_y0i=xk##i;\
   xk##r =(_y0r * _y0r + (gkr) * _y0i * _y0i);\
   xk##i = 2.0 * _y0r * _y0i;\
}
*/
#define sse2_nested_eq(__xkr, __xki, __gkr)                      \
{                                                                \
   Y__M128D __x0, __x1, __x2, __x3, __x4;                        \
   Y_MM_INTER_PAIR(__xkr, __xki);     /* #0          */          \
   Y_MM_MOV_PD( __x0, __xkr);                                    \
   Y_MM_MOV_PD( __x1, __xki);                                    \
   Y_MM_MUL_PD( __xkr, __xkr, __xkr);                            \
   Y_MM_MUL_PD( __x1, __x1, __x1);                               \
   Y_MM_MUL_PD( __xki, __xki, __x0);                             \
   Y_MM_MUL_PD( __x1, __x1, __gkr);                              \
   Y_MM_ADD_PD( __xki, __xki, __xki);                            \
   Y_MM_ADD_PD( __xkr, __xkr, __x1);                             \
   Y_MM_INTER_PAIR(__xkr, __xki);                                \
}

/* 
#define square_nested_1_4(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gki + 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki + 1.0) * _r1 - gkr * _r2);\
   xk##r -= _r0;\
   xmk##r -= _r0;\
   xk##i -= _r1;\
   xmk##i += _r1;\
}
 */
#define sse2_square_nested_1_4(_xkr, _xmkr, _xki, _xmki, _gkr, _gki) \
{                                                                \
   Y__M128D __x0, __x1, __x2, __x3, __x4;                        \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #0          */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #4          */          \
   Y_MM_MOV_PD(__x0, _xkr);           /* #5          */          \
   Y_MM_MOV_PD(__x1, _xki);           /* #6          */          \
   Y_MM_SUB_PD(__x0, __x0, _xmkr);    /* #8-11       */          \
   Y_MM_MOV_PD(__x2, _xki);           /* #8-13       */          \
   Y_MM_MOV_PD(__x3, _xkr);           /* #9-14       */          \
   Y_MM_ADD_PD(__x1, __x1, _xmki);    /* #10-13      */          \
   Y_MM_MUL_PD(__x0, __x0, MM_0_5);   /* #12-17      */          \
   Y_MM_SUB_PD(_xkr, _xkr, _xki);     /* #16-19      */          \
   Y_MM_MUL_PD(__x1, __x1, MM_0_5);   /* #14-19      */          \
   Y_MM_MUL_PD(_xki, _xki, __x3);     /* #18-23      */          \
   Y_MM_MOV_PD(__x4, __x0);           /* #18-23      */          \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #20-23      */          \
   Y_MM_MOV_PD(__x2, __x0);           /* #21-26      */          \
   Y_MM_SUB_PD(__x4, __x4, __x1);     /* #22-25      */          \
   Y_MM_ADD_PD(__x2, __x2, __x1);     /* #24-27      */          \
   Y_MM_MUL_PD(_xkr, _xkr, __x3);     /* #26-31      */          \
   Y_MM_MOV_PD(__x3, _xmkr);          /* #27-32      */          \
   Y_MM_MUL_PD(__x1, __x1, __x0);     /* #28-33      */          \
   Y_MM_MOV_PD(__x0, _xmki);          /* #29-34      */          \
   Y_MM_MUL_PD(__x2, __x2, __x4);     /* #30-35      */          \
   Y_MM_MOV_PD(__x4, _gki);           /* #31-36   ** */          \
   Y_MM_SUB_PD(_xmkr, _xmkr, _xmki);  /* #32-34      */          \
   Y_MM_MUL_PD(_xmki, _xmki, __x3);   /* #34-39      */          \
   Y_MM_ADD_PD(__x3, __x3, __x0);     /* #36-39      */          \
   Y_MM_ADD_PD(__x1, __x1, __x1);     /* #38-41      */          \
   Y_MM_ADD_PD(__x4, __x4, MM_1_0);   /* #40-43   ** */          \
   Y_MM_MOV_PD(__x0, __x2);           /* #40-45      */          \
   Y_MM_MUL_PD(_xmkr, _xmkr, __x3);   /* #42-47      */          \
   Y_MM_MOV_PD(__x3, __x1);           /* #42-47      */          \
   Y_MM_ADD_PD(_xki, _xki, _xki);     /* #44-47      */          \
   Y_MM_MUL_PD(__x1, __x1, _gkr);     /* #46-51   ** */          \
   Y_MM_ADD_PD(_xmki, _xmki, _xmki);  /* #           */          \
   Y_MM_MUL_PD(__x0, __x0, __x4);     /* #48-53      */          \
   Y_MM_MUL_PD(__x2, __x2, _gkr);     /* #50-55   ** */          \
   Y_MM_MUL_PD(__x3, __x3, __x4);     /* #52-57      */          \
   Y_MM_ADD_PD(__x0, __x0, __x1);     /* #54-57   ** */          \
   Y_MM_SUB_PD(__x3, __x3, __x2);     /* #58-61   ** */          \
   Y_MM_SUB_PD(_xkr, _xkr, __x0);     /* #60-63   ** */          \
   Y_MM_SUB_PD(_xmkr, _xmkr, __x0);   /* #62-65   ** */          \
   Y_MM_SUB_PD(_xki, _xki, __x3);     /* #64-67   ** */          \
   Y_MM_ADD_PD(_xmki, _xmki, __x3);   /* #66-69   ** */          \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #67         */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #71         */          \
}


/*
#define square_nested_1_2(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gkr - 1.0) * _r2 - gki * _r1);\
   _r1= ((gkr - 1.0) * _r1 + gki * _r2);\
   xk##r += _r0;\
   xmk##r += _r0;\
   xk##i += _r1;\
   xmk##i -= _r1;\
}
*/
#define sse2_square_nested_1_2(_xkr, _xmkr, _xki, _xmki, _gkr, _gki) \
{                                                                \
   Y__M128D __x0, __x1, __x2, __x3, __x4;                        \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #0          */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #4          */          \
   Y_MM_MOV_PD(__x0, _xkr);           /* #5          */          \
   Y_MM_MOV_PD(__x1, _xki);           /* #6          */          \
   Y_MM_SUB_PD(__x0, __x0, _xmkr);    /* #8-11       */          \
   Y_MM_MOV_PD(__x2, _xki);           /* #8-13       */          \
   Y_MM_MOV_PD(__x3, _xkr);           /* #9-14       */          \
   Y_MM_ADD_PD(__x1, __x1, _xmki);    /* #10-13      */          \
   Y_MM_MUL_PD(__x0, __x0, MM_0_5);   /* #12-17      */          \
   Y_MM_SUB_PD(_xkr, _xkr, _xki);     /* #16-19      */          \
   Y_MM_MUL_PD(__x1, __x1, MM_0_5);   /* #14-19      */          \
   Y_MM_MUL_PD(_xki, _xki, __x3);     /* #18-23      */          \
   Y_MM_MOV_PD(__x4, __x0);           /* #18-23      */          \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #20-23      */          \
   Y_MM_MOV_PD(__x2, __x0);           /* #21-26      */          \
   Y_MM_SUB_PD(__x4, __x4, __x1);     /* #22-25      */          \
   Y_MM_ADD_PD(__x2, __x2, __x1);     /* #24-27      */          \
   Y_MM_MUL_PD(_xkr, _xkr, __x3);     /* #26-31      */          \
   Y_MM_MOV_PD(__x3, _xmkr);          /* #27-32      */          \
   Y_MM_MUL_PD(__x1, __x1, __x0);     /* #28-33      */          \
   Y_MM_MOV_PD(__x0, _xmki);          /* #29-34      */          \
   Y_MM_MUL_PD(__x2, __x2, __x4);     /* #30-35      */          \
   Y_MM_MOV_PD(__x4, _gkr);           /* #31-36   ** */          \
   Y_MM_SUB_PD(_xmkr, _xmkr, _xmki);  /* #32-34      */          \
   Y_MM_MUL_PD(_xmki, _xmki, __x3);   /* #34-39      */          \
   Y_MM_ADD_PD(__x3, __x3, __x0);     /* #36-39      */          \
   Y_MM_ADD_PD(__x1, __x1, __x1);     /* #38-41      */          \
   Y_MM_SUB_PD(__x4, __x4, MM_1_0);   /* #40-43   ** */          \
   Y_MM_MOV_PD(__x0, __x2);           /* #40-45      */          \
   Y_MM_MUL_PD(_xmkr, _xmkr, __x3);   /* #42-47      */          \
   Y_MM_MOV_PD(__x3, __x1);           /* #42-47      */          \
   Y_MM_ADD_PD(_xki, _xki, _xki);     /* #44-47      */          \
   Y_MM_MUL_PD(__x1, __x1, _gki);     /* #46-51   ** */          \
   Y_MM_ADD_PD(_xmki, _xmki, _xmki);  /* #56-59      */          \
   Y_MM_MUL_PD(__x0, __x0, __x4);     /* #48-53      */          \
   Y_MM_MUL_PD(__x2, __x2, _gki);     /* #50-55   ** */          \
   Y_MM_MUL_PD(__x3, __x3, __x4);     /* #52-57      */          \
   Y_MM_SUB_PD(__x0, __x0, __x1);     /* #54-57   ** */          \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #58-61   ** */          \
   Y_MM_ADD_PD(_xkr, _xkr, __x0);     /* #60-63   ** */          \
   Y_MM_ADD_PD(_xmkr, _xmkr, __x0);   /* #62-65   ** */          \
   Y_MM_ADD_PD(_xki, _xki, __x3);     /* #64-67   ** */          \
   Y_MM_SUB_PD(_xmki, _xmki, __x3);   /* #66-69   ** */          \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #67         */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #71         */          \
}

/*
#define square_nested_3_4(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gki - 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki - 1.0) * _r1 - gkr * _r2);\
   xk##r += _r0;\
   xmk##r += _r0;\
   xk##i += _r1;\
   xmk##i -= _r1;\
}
*/
#define sse2_square_nested_3_4(_xkr, _xmkr, _xki, _xmki, _gkr, _gki) \
{                                                                \
   Y__M128D __x0, __x1, __x2, __x3, __x4;                        \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #0          */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #4          */          \
   Y_MM_MOV_PD(__x0, _xkr);           /* #5          */          \
   Y_MM_MOV_PD(__x1, _xki);           /* #6          */          \
   Y_MM_SUB_PD(__x0, __x0, _xmkr);    /* #8-11       */          \
   Y_MM_MOV_PD(__x2, _xki);           /* #8-13       */          \
   Y_MM_MOV_PD(__x3, _xkr);           /* #9-14       */          \
   Y_MM_ADD_PD(__x1, __x1, _xmki);    /* #10-13      */          \
   Y_MM_MUL_PD(__x0, __x0, MM_0_5);   /* #12-17      */          \
   Y_MM_SUB_PD(_xkr, _xkr, _xki);     /* #16-19      */          \
   Y_MM_MUL_PD(__x1, __x1, MM_0_5);   /* #14-19      */          \
   Y_MM_MUL_PD(_xki, _xki, __x3);     /* #18-23      */          \
   Y_MM_MOV_PD(__x4, __x0);           /* #18-23      */          \
   Y_MM_ADD_PD(__x3, __x3, __x2);     /* #20-23      */          \
   Y_MM_MOV_PD(__x2, __x0);           /* #21-26      */          \
   Y_MM_SUB_PD(__x4, __x4, __x1);     /* #22-25      */          \
   Y_MM_ADD_PD(__x2, __x2, __x1);     /* #24-27      */          \
   Y_MM_MUL_PD(_xkr, _xkr, __x3);     /* #26-31      */          \
   Y_MM_MOV_PD(__x3, _xmkr);          /* #27-32      */          \
   Y_MM_MUL_PD(__x1, __x1, __x0);     /* #28-33      */          \
   Y_MM_MOV_PD(__x0, _xmki);          /* #29-34      */          \
   Y_MM_MUL_PD(__x2, __x2, __x4);     /* #30-35      */          \
   Y_MM_MOV_PD(__x4, _gki);           /* #31-36   ** */          \
   Y_MM_SUB_PD(_xmkr, _xmkr, _xmki);  /* #32-34      */          \
   Y_MM_MUL_PD(_xmki, _xmki, __x3);   /* #34-39      */          \
   Y_MM_ADD_PD(__x3, __x3, __x0);     /* #36-39      */          \
   Y_MM_ADD_PD(__x1, __x1, __x1);     /* #38-41      */          \
   Y_MM_SUB_PD(__x4, __x4, MM_1_0);   /* #40-43   ** */          \
   Y_MM_MOV_PD(__x0, __x2);           /* #40-45      */          \
   Y_MM_MUL_PD(_xmkr, _xmkr, __x3);   /* #42-47      */          \
   Y_MM_MOV_PD(__x3, __x1);           /* #42-47      */          \
   Y_MM_ADD_PD(_xki, _xki, _xki);     /* #44-47      */          \
   Y_MM_MUL_PD(__x1, __x1, _gkr);     /* #46-51   ** */          \
   Y_MM_ADD_PD(_xmki, _xmki, _xmki);  /* #56-59      */          \
   Y_MM_MUL_PD(__x0, __x0, __x4);     /* #48-53      */          \
   Y_MM_MUL_PD(__x2, __x2, _gkr);     /* #50-55   ** */          \
   Y_MM_MUL_PD(__x3, __x3, __x4);     /* #52-57      */          \
   Y_MM_ADD_PD(__x0, __x0, __x1);     /* #54-57   ** */          \
   Y_MM_SUB_PD(__x3, __x3, __x2);     /* #58-61   ** */          \
   Y_MM_ADD_PD(_xkr, _xkr, __x0);     /* #60-63   ** */          \
   Y_MM_ADD_PD(_xmkr, _xmkr, __x0);   /* #62-65   ** */          \
   Y_MM_ADD_PD(_xki, _xki, __x3);     /* #64-67   ** */          \
   Y_MM_SUB_PD(_xmki, _xmki, __x3);   /* #66-69   ** */          \
   Y_MM_INTER_PAIR(_xkr, _xki);       /* #67         */          \
   Y_MM_INTER_PAIR(_xmkr, _xmki);     /* #71         */          \
}

/*
#define conv_nested(xk,xmk,yk,ymk,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;\
   _y0i=(xk##i + xmk##i)*0.25;\
   _y1r=(yk##r - ymk##r);\
   _y1i=(yk##i + ymk##i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);\
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);\
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);\
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);\
   _y3r= ((1.0 +(gkr)) * _y2r - (gki) * _y2i);\
   _y3i= ((1.0 +(gkr)) * _y2i + (gki) * _y2r);\
   xk##r = _y0r - _y3r;\
   xk##i = _y0i - _y3i;\
   xmk##r = _y1r - _y3r;\
   xmk##i = _y1i + _y3i;\
}
*/
