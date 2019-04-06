/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor, Klaus Kastens
 
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
/* This are macros to be used with pseudo complex FFT                       */

/* Include all common macros */
#include "yx86comm.h"

/***************  SINGLE LOAD AND STORE MACROS ***************/


#ifndef USE_ASM
/* To load from real padded array of data as a complex array */
#define cplx_data_to_local(_t,_array,_k)                      \
	_t##r = _array[addr(( _k )<<1)];                      \
	_t##i = _array[addr(( _k )<<1)+1];

#else

#define cplx_data_to_local(_t,_array,_k) \
	__asm__ volatile("        \n"                  \
"movl %3,%%ebx           #        \n"                  \
"movl %%ebx,%%eax        #1  U    \n"                  \
"andl $-512,%%ebx        #1  V    \n"                  \
"shrl $7,%%ebx           #2  U    \n"                  \
"addl %%eax,%%eax	 #2  V    \n"                  \
"addl %%eax,%%ebx	 #3  U    \n"                  \
"fldl (%%esi,%%ebx,8)             \n"                  \
"fldl 8(%%esi,%%ebx,8)            \n"                  \
"prefetchnta 64(%%esi,%%ebx,8)    \n"                  \
"fxch %%st(1)                     \n"                  \
"fstpl %0                         \n"                  \
"fstpl %1                         \n"                  \
: :"m"(_t##r),"m"(_t##i),"S"(_array),"m"(_k):          \
 "ax","bx","memory");

#endif


#ifndef USE_ASM

#define cplx_data_to_local_p(_t,_pd,_array,_k)         \
	_pd = _array + addr(( _k)<<1);                 \
	_t##r = *(_pd );                               \
	_t##i = *( _pd + 1);                           \
        prefetch_p(_pd);

#else

# define cplx_data_to_local_p(_t,_pd,_array,_k) \
	__asm__ volatile("                 \n"  \
"movl %%ebx,%%eax        #1  U             \n"  \
"andl $-512,%%ebx        #1  V             \n"  \
"shrl $7,%%ebx           #2  U             \n"  \
"addl %%eax,%%eax	 #2  V             \n"  \
"addl %%eax,%%ebx	 #3  U             \n"  \
"fldl (%%esi,%%ebx,8)    #4  (if cache hit)\n"  \
"fldl 8(%%esi,%%ebx,8)   #5                \n"  \
"prefetchnta 64(%%esi,%%ebx,8)             \n"  \
"leal (%%esi,%%ebx,8),%%ebx                \n"  \
"fstpl %2                #7                \n"  \
"fstpl %1                #8                \n"  \
:"=&b"(_pd) :"m"(_t##r),"m"(_t##i),"S"(_array),"0"(_k):\
	"ax","memory");

#endif

#ifndef USE_ASM

#define cplx_data_to_local_pp(_t,_index,_pd) \
	_t##r = *( _pd + _index);            \
	_t##i = *( _pd + _index + 1);        \
        prefetch_p ( _pd + _index);          \

#else

# define cplx_data_to_local_pp(_t,_index,_pd)   \
	__asm__ volatile("                   \n" \
"fldl (%%esi,%%ebx,8)    #1  (if cache hit)  \n" \
"fldl 8(%%esi,%%ebx,8)   #2                  \n" \
"prefetchnta 64(%%esi,%%ebx,8)               \n" \
"fxch %%st(1)                                \n" \
"fstpl %0                #7                  \n" \
"fstpl %1                #8                  \n" \
: :"m"(_t##r),"m"(_t##i),"S"(_pd),"b"(_index):"memory");

#endif

#ifndef USE_ASM
/* to read some trigonomeric factors */
#define cplx_data_trig_load_2(_t0,_t1,_p) \
   _t0##r = *( _p );                      \
   _t0##i = *( _p + 1);                   \
   _t1##r = *( _p + 2);                   \
   _t1##i = *( _p + 3);

#define cplx_data_trig_load_3(_t0,_t1,_t2,_p) \
   _t0##r = *( _p );                          \
   _t0##i = *( _p + 1);                       \
   _t1##r = *( _p + 2);                       \
   _t1##i = *( _p + 3);                       \
   _t2##r = *( _p + 4);                       \
   _t2##i = *( _p + 5);

#define cplx_data_trig_load_4(_t0,_t1,_t2,_t3,_p) \
   _t0##r = *( _p );                              \
   _t0##i = *( _p + 1);                           \
   _t1##r = *( _p + 2);                           \
   _t1##i = *( _p + 3);                           \
   _t2##r = *( _p + 4);                           \
   _t2##i = *( _p + 5);                           \
   _t3##r = *( _p + 6);                           \
   _t3##i = *( _p + 7);

#else

#define cplx_data_trig_load_2(_t0,_t1,_p)    \
   __asm__ volatile( "                    \n"  \
"fldl (%%esi)         #1   t0r            \n"  \
"fldl 8(%%esi)        #2   t0i,t0r        \n"  \
"fldl 16(%%esi)       #3   t1r,t0i,t0r    \n"  \
"fldl 24(%%esi)       #4   t1i,t1r,t0i,t0r\n"  \
"fxch %%st(3)         #5   t0r,t1r,t0i,t1i\n"  \
"fstpl %1             #6   t1r,t0i,t1i    \n"  \
"fstpl %3             #7   t0i,t1i        \n"  \
"fstpl %2             #8   t1i            \n"  \
"fstpl %4             #9                  \n"  \
"prefetchnta 64(%%esi)                    \n"  \
: :"S" (_p), "m" ( _t0##r), "m" (_t0##i),      \
"m" (_t1##r),"m"(_t1##i)  :"st","memory");

#define cplx_data_trig_load_3(_t0,_t1,_t2,_p)           \
   __asm__ volatile("                              \n"  \
"fldl (%%esi)         #1   t0r                     \n"  \
"fldl 8(%%esi)        #2   t0i,t0r                 \n"  \
"fldl 16(%%esi)       #3   t1r,t0i,t0r             \n"  \
"fldl 24(%%esi)       #4   t1i,t1r,t0i,t0r         \n"  \
"fldl 32(%%esi)       #3   t2r,t1i,t1r,t0i,t0r     \n"  \
"fldl 40(%%esi)       #4   t2i,t2r,t1i,t1r,t0i,t0r \n"  \
"fxch %%st(5)         #5   t0r,t2r,t1i,t1r,t0i,t2i \n"  \
"fstpl %1             #6   t2r,t1i,t1r,t0i,t2i     \n"  \
"fstpl %5             #7   t1i,t1r,t0i,t2i         \n"  \
"fstpl %4             #8   t1r,t0i,t2i             \n"  \
"fstpl %3             #9   t0i,t2i                 \n"  \
"fstpl %2             #10  t2i                     \n"  \
"fstpl %6             #11                          \n"  \
"prefetchnta 128(%%esi)                             \n"  \
"prefetchnta 192(%%esi)                             \n"  \
: :"S" (_p), "m" ( _t0##r), "m" (_t0##i), "m" (_t1##r), \
 "m"(_t1##i),"m" ( _t2##r), "m" (_t2##i):"st","memory");
#endif


/**************** PSEUDO_COMPLEX GENERAL AUXILIAR MACROS **************/

/**********************************************************************
   cplx_load_muladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
 where all the variables are pseudo-complex and _t is temporal aux var.
**********************************************************************/
#ifndef USE_ASM

#define cplx_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1)    \
	{                                                 \
		y_limb_t _ar,_ai,_br,_bi;                 \
		_br= _array[addr(_k1<<1)];                \
		_bi= _array[addr(_k1<<1)+1];              \
                prefetch_p(&_array[addr(_k1<<1)]);        \
		_ar = (_br) * (_f##r) - (_bi) * (_f##i);  \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r);  \
		_t0##r = _array[addr(_k0<<1)] + _ar;      \
		_t1##r = _array[addr(_k0<<1)] - _ar;      \
		_t0##i = _array[addr(_k0<<1)+1] + _ai;    \
		_t1##i = _array[addr(_k0<<1)+1] - _ai;    \
                prefetch_p(&_array[addr(_k0<<1)]);        \
	}

#else

# define cplx_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1)       \
__asm__ volatile("                                   \n"      \
"movl %3 ,%%ebx                                      \n"      \
"movl %4 ,%%edx                                      \n"      \
"movl %%ebx,%%eax             # 1  U                 \n"      \
"andl $-512,%%ebx             # 1  V                 \n"      \
"shrl $7,%%ebx                # 2  U                 \n"      \
"addl %%eax,%%eax	     # 2  V                  \n"      \
"addl %%eax,%%ebx	     # 3  U                  \n"      \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache\n"     \
"fldl 8(%%esi,%%ebx,8)        # 5 t1i,t1r            \n"      \
"prefetchnta 64(%%esi,%%ebx,8)                       \n"      \
"fldl %0     		     # 6 fr,t1i,t1r          \n"      \
"fmul %%st(1)                 # 7-11 IR,t1i,t1r      \n"      \
"fldl %1                      # 8    fi,IR,t1i,t1r   \n"      \
"fmul %%st(3)                 # 9-13 RI,IR,t1i,t1r   \n"      \
"fldl %0                      # 10   fr,RI,IR,t1i,t1r\n"      \
"fmulp %%st,%%st(4)           # 11-15   RI,IR,t1i,RR \n"      \
"fldl %1			     # 12   fi,RI,ir,t1i,RR\n"\
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR   \n"      \
"movl %%edx,%%eax             # 14  U                \n"      \
"andl $-512,%%edx             # 14  V                \n"      \
"shrl $7,%%edx                # 15  U                \n"      \
"addl %%eax,%%eax	     # 15  V                 \n"      \
"addl %%eax,%%edx	     # 16  U                 \n"      \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr    \n"      \
"fxch %%st(2)                 # 17     rr,II,NT1I    \n"      \
"fldl (%%esi,%%edx,8)         # 18     nt0r,rr,ii,NT1I\n"     \
"fxch %%st(1)                 # 18     rr,nt0r,ii,NT1I\n"     \
"fsubp %%st,%%st(2)           # 19-23  nt0r,NT1R,NT1I\n"      \
"fldl 8(%%esi,%%edx,8)        # 20     nt0i,nt0r,NT1R,NT1I\n" \
"fsub %%st(3)                 # 21-25  T1I,nt0r,NT1R,nt1i\n"  \
"fxch %%st(3)                 # 21     nt1i,nt0r,NT1R,T1I\n"  \
"faddl 8(%%esi,%%edx,8)       # 22-26  T0I,nt0r,nt1r,T1I\n"   \
"fxch %%st(1)                 # 22     nt0r,T0I,nt1r,T1I\n"   \
"fsub %%st(2)                 # 23-27  T1R,T0I,nt1r,TI1\n"    \
"fxch %%st(2)                 # 23     nt1r,T0I,T1R,T1I\n"    \
"faddl (%%esi,%%edx,8)        # 24-28  T0R,T0I,T1R,t1i\n"     \
"fxch %%st(3)                 # 25     t1i,T0I,T1R,T0R\n"     \
"prefetchnta 64(%%esi,%%edx,8)                        \n"     \
: : "m"(_f##r),"m"(_f##i),"S"(_array),                        \
"m"(_k1),"m"(_k0) :"ax","bx","dx");                           \
      __asm__ volatile("                              \n"     \
"fstpl  %3           # 27                              \n"    \
"fstpl  %1           # 28                              \n"    \
"fstpl  %2           # 29                              \n"    \
"fstpl  %0           # 30                              \n"    \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):   \
"st","memory");

#endif

#ifndef USE_ASM

#define cplx_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
	{                                                          \
		y_limb_t _ar,_ai,_br,_bi;                          \
		_pd1= _array + addr((_k1)<<1);                     \
		_br= *( _pd1 );                                    \
		_bi= *( _pd1 + 1);                                 \
                prefetch_p(_pd1);                                  \
		_ar = (_br) * (_f##r) - (_bi) * (_f##i);           \
		_pd0 = _array + addr((_k0)<<1);                    \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r);           \
		_t0##r = *( _pd0 ) + _ar;                          \
		_t1##r = *( _pd0 ) - _ar;                          \
		_t0##i = *( _pd0 + 1) + _ai;                       \
		_t1##i = *( _pd0 + 1) - _ai;                       \
                prefetch_p(_pd0);                                  \
	}

#else

#define cplx_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
__asm__ volatile("                                            \n"  \
"movl %%ebx,%%eax             # 1  U                          \n"  \
"andl $-512,%%ebx             # 1  V                          \n"  \
"shrl $7,%%ebx                # 2  U                          \n"  \
"addl %%eax,%%eax	     # 2  V                           \n"  \
"addl %%eax,%%ebx	     # 3  U                           \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache        \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5 t1i,t1r                     \n"  \
"leal (%%esi,%%ebx,8),%%ebx                                   \n"  \
"prefetchnta 64(%%ebx)                                        \n"  \
"fldl %2     		     # 6 fr,t1i,t1r                   \n"  \
"fmul %%st(1)                 # 7-11 IR,t1i,t1r               \n"  \
"fldl %3                      # 8    fi,IR,t1i,t1r            \n"  \
"fmul %%st(3)                 # 9-13 RI,IR,t1i,t1r            \n"  \
"fldl %2                      # 10   fr,RI,IR,t1i,t1r         \n"  \
"fmulp %%st,%%st(4)           # 11-15   RI,IR,t1i,RR          \n"  \
"fldl %3			     # 12   fi,RI,ir,t1i,RR   \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR            \n"  \
"movl %%edx,%%eax             # 14  U                         \n"  \
"andl $-512,%%edx             # 14  V                         \n"  \
"shrl $7,%%edx                # 15  U                         \n"  \
"addl %%eax,%%eax	     # 15  V                          \n"  \
"addl %%eax,%%edx	     # 16  U                          \n"  \
"leal (%%esi,%%edx,8),%%edx                                   \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr             \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I             \n"  \
"fldl (%%edx)                 # 18     nt0r,rr,ii,NT1I        \n"  \
"fxch %%st(1)                 # 18     rr,nt0r,ii,NT1I        \n"  \
"fsubp %%st,%%st(2)           # 19-23  nt0r,NT1R,NT1I         \n"  \
"fldl 8(%%edx)                # 20     nt0i,nt0r,NT1R,NT1I    \n"  \
"fsub %%st(3)                 # 21-25  T1I,nt0r,NT1R,nt1i     \n"  \
"fxch %%st(3)                 # 21     nt1i,nt0r,NT1R,T1I     \n"  \
"faddl 8(%%edx)               # 22-26  T0I,nt0r,nt1r,T1I      \n"  \
"fxch %%st(1)                 # 22     nt0r,T0I,nt1r,T1I      \n"  \
"prefetchnta 64(%%edx)                                        \n"  \
"fsub %%st(2)                 # 23-27  T1R,T0I,nt1r,T1I       \n"  \
"fxch %%st(2)                 # 23     nt1r,T0I,T1R,T1I       \n"  \
"faddl (%%edx)                # 24-28  T0R,T0I,T1R,t1i        \n"  \
"fxch %%st(3)                 # 25     t1i,T0I,T1R,T0R        \n"  \
: "=&b"(_pd1) ,"=&d" (_pd0): "m"(_f##r),"m"(_f##i),"S"(_array),    \
"0"(_k1),"1"(_k0) :"ax");                                          \
      __asm__ volatile("                                      \n"  \
"fstpl  %3           # 27                                     \n"  \
"fstpl  %1           # 28                                     \n"  \
"fstpl  %2           # 29                                     \n"  \
"fstpl  %0           # 30                                     \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):        \
 "st","memory");

#endif

#ifndef USE_ASM

#define cplx_load_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1)        \
	{                                                          \
		y_limb_t _ar,_ai,_br,_bi;                          \
		_br= *( _pd1 + _index );                           \
		_bi= *( _pd1 + _index + 1);                        \
                prefetch_p( _pd1 + _index);                        \
		_ar = (_br) * (_f##r) - (_bi) * (_f##i);           \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r);           \
		_t0##r = *( _pd0 + _index) + _ar;                  \
		_t1##r = *( _pd0 + _index) - _ar;                  \
		_t0##i = *( _pd0 + _index + 1) + _ai;              \
		_t1##i = *( _pd0 + _index + 1) - _ai;              \
                prefetch_p( _pd0 + _index);                        \
	}

#else

#define cplx_load_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1)        \
__asm__ volatile("                                            \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache        \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5 t1i,t1r                     \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                                \n"  \
"fldl %0     		     # 6 fr,t1i,t1r                   \n"  \
"fmul %%st(1)                 # 7-11 IR,t1i,t1r               \n"  \
"fldl %1                      # 8    fi,IR,t1i,t1r            \n"  \
"fmul %%st(3)                 # 9-13 RI,IR,t1i,t1r            \n"  \
"fldl %0                      # 10   fr,RI,IR,t1i,t1r         \n"  \
"fmulp %%st,%%st(4)           # 11-15   RI,IR,t1i,RR          \n"  \
"fldl %1			     # 12   fi,RI,ir,t1i,RR   \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR            \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr             \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I             \n"  \
"fldl (%%edi,%%ebx,8)         # 18     nt0r,rr,ii,NT1I        \n"  \
"fxch %%st(1)                 # 18     rr,nt0r,ii,NT1I        \n"  \
"fsubp %%st,%%st(2)           # 19-23  nt0r,NT1R,NT1I         \n"  \
"fldl 8(%%edi,%%ebx,8)        # 20     nt0i,nt0r,NT1R,NT1I    \n"  \
"fsub %%st(3)                 # 21-25  T1I,nt0r,NT1R,nt1i     \n"  \
"fxch %%st(3)                 # 21     nt1i,nt0r,NT1R,T1I     \n"  \
"faddl 8(%%edi,%%ebx,8)       # 22-26  T0I,nt0r,nt1r,T1I      \n"  \
"fxch %%st(1)                 # 22     nt0r,T0I,nt1r,T1I      \n"  \
"prefetchnta 64(%%edi,%%ebx,8)                                \n"  \
"fsub %%st(2)                 # 23-27  T1R,T0I,nt1r,T1I       \n"  \
"fxch %%st(2)                 # 23     nt1r,T0I,T1R,T1I       \n"  \
"faddl (%%edi,%%ebx,8)        # 24-28  T0R,T0I,T1R,t1i        \n"  \
"fxch %%st(3)                 # 25     t1i,T0I,T1R,T0R        \n"  \
: : "m"(_f##r),"m"(_f##i),"b"(_index),"S"(_pd1),"D"(_pd0) );       \
     __asm__ volatile("                                       \n"  \
"fstpl  %3           # 27                                     \n"  \
"fstpl  %1           # 28                                     \n"  \
"fstpl  %2           # 29                                     \n"  \
"fstpl  %0           # 30                                     \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):        \
"st","memory");

#endif


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
#ifndef USE_ASM

#define cplx_load_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1) \
	{                                                      \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci;              \
		_br= _array[addr(_k0<<1)];                     \
		_bi= _array[addr(_k0<<1)+1];                   \
                prefetch_p(&_array[addr(_k0<<1)]);             \
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i);     \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r);     \
		_br= _array[addr(_k1<<1)];                     \
		_bi= _array[addr(_k1<<1)+1];                   \
                prefetch_p(&_array[addr(_k1<<1)]);             \
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i);     \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r);     \
		_t0##r = _ar + _cr;                            \
		_t1##r = _ar - _cr;                            \
		_t0##i = _ai + _ci;                            \
		_t1##i = _ai - _ci;                            \
	}
#else

#define cplx_load_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1)        \
 __asm__ volatile("                                              \n"  \
"movl %5,%%ebx                                                   \n"  \
"movl %6,%%edx                                                   \n"  \
"movl %%ebx,%%eax             # 1  U                             \n"  \
"andl $-512,%%ebx             # 1  V                             \n"  \
"shrl $7,%%ebx                # 2  U                             \n"  \
"addl %%eax,%%eax	     # 2  V                              \n"  \
"addl %%eax,%%ebx	     # 3  U                              \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache           \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r                   \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                                   \n"  \
"fldl %0     		     # 6      fr,t1i,t1r                 \n"  \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r                \n"  \
"fldl %1                      # 8      fi,IR,t1i,t1r             \n"  \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r             \n"  \
"fldl %0                      # 10     fr,RI,IR,t1i,t1r          \n"  \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR              \n"  \
"fldl %1			     # 12     fi,RI,ir,t1i,RR    \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR               \n"  \
"movl %%edx,%%eax             # 14  U                            \n"  \
"andl $-512,%%edx             # 14  V                            \n"  \
"shrl $7,%%edx                # 15  U                            \n"  \
"addl %%eax,%%eax	     # 15  V                             \n"  \
"addl %%eax,%%edx	     # 16  U                             \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr                \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I                \n"  \
"fldl (%%esi,%%edx,8)         # 18     nt0r,rr,ii,NT1I           \n"  \
"fmull %2                     # 19-23  RR0,rr,ii,NT1I            \n"  \
"fxch %%st(1)                 # 19     rr,RR0,ii,NT1I            \n"  \
"fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I             \n"  \
"fldl 8(%%esi,%%edx,8)        # 21     t0i,RR,NT1R,NT1I          \n"  \
"fmull %3                     # 22-26  II,RR,NT1R,nt1i           \n"  \
"fldl (%%esi,%%edx,8)         # 23     t0r,II,RR,NT1R,nt1i       \n"  \
"fmull %3                     # 24-28  RI,II,rr,NT1R,nt1i        \n"  \
"fldl 8(%%esi,%%edx,8)        # 25     t0i,RI,II,rr,nt1r,nt1i    \n"  \
"prefetchnta 64(%%esi,%%edx,8)                                   \n"  \
"fmull %2                     # 26-30  IR,RI,II,rr,nt1r,nt1i     \n"  \
"fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i     \n"  \
"fsubp %%st,%%st(2)           # 27-31  RI,NTOR,IR,nt1r,nt1i      \n"  \
"fxch %%st(2)                 # 27     IR,NTOR,RI,nt1r,nt1i      \n"  \
"fld  %%st(4)                 # 28     nt1i,IR,NTOR,RI,nt1r,nt1i \n"  \
"fld  %%st(4)                 # 29     nt1r,nt1i,IR,NTOR,RI,nt1r,nt1i\n"\
"fadd %%st(3)                 # 30-34  T0R,nt1i,IR,NTOR,ri,nt1r,nt1i\n" \
"fxch %%st(4)                 # 30     ri,nt1i,IR,NTOR,TOR,nt1r,nt1i\n" \
"faddp %%st,%%st(2)           # 31-35  nt1i,NTOI,ntor,TOR,nt1r,nt1i\n"\
"fxch %%st(2)                 # 31     nt0r,NTOI,nt1i,TOR,nt1r,nt1i\n"\
"fsubp %%st,%%st(4)           # 32-36  NTOI,nt1i,TOR,T1R,nt1i    \n"  \
"fadd  %%st,%%st(1)           # 33-37  ntoi,TOI,TOR,T1R,nt1i     \n"  \
"fsubp %%st,%%st(4)           # 34-38  TOI,TOR,T1R,T1I           \n"  \
: : "m"(_f1##r),"m"(_f1##i),"m"(_f0##r),                              \
"m"(_f0##i),"S"(_array),"m"(_k1),"m"(_k0) :"ax","bx","dx","st");      \
      __asm__ volatile("                                         \n"  \
"fstpl  %1           # 37                                        \n"  \
"fstpl  %0           # 38                                        \n"  \
"fstpl  %2           # 39                                        \n"  \
"fstpl  %3           # 40                                        \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):           \
"memory", "st");

#endif

#ifndef USE_ASM

# define cplx_load_mulmuladdsub_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{                                                                   \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci;                           \
		_pd0 = _array + addr((_k0)<<1);                             \
		_br= *( _pd0 );                                             \
		_bi= *( _pd0 + 1);                                          \
                prefetch_p(_pd0);                                           \
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i);                  \
		_pd1 = _array + addr((_k1)<<1);                             \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r);                  \
		_br= *( _pd1 );                                             \
		_bi= *( _pd1 + 1);                                          \
                prefetch_p(_pd1);                                           \
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i);                  \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r);                  \
		_t0##r = _ar + _cr;                                         \
		_t1##r = _ar - _cr;                                         \
		_t0##i = _ai + _ci;                                         \
		_t1##i = _ai - _ci;                                         \
	}

#else

#define cplx_load_mulmuladdsub_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
 __asm__ volatile("                                                   \n"  \
"movl %%ebx,%%eax            # 1  U                                   \n"  \
"andl $-512,%%ebx            # 1  V                                   \n"  \
"shrl $7,%%ebx               # 2  U                                   \n"  \
"addl %%eax,%%eax	     # 2  V                                   \n"  \
"addl %%eax,%%ebx	     # 3  U                                   \n"  \
"fldl (%%esi,%%ebx,8)        # 4 (14) if not in cache                 \n"  \
"fldl 8(%%esi,%%ebx,8)       # 5      t1i,t1r                         \n"  \
"leal (%%esi,%%ebx,8),%%ebx                                           \n"  \
"prefetchnta 64(%%ebx)                                                \n"  \
"fldl %2     		     # 6      fr,t1i,t1r                      \n"  \
"fmul %%st(1)                # 7-11   IR,t1i,t1r                      \n"  \
"fldl %3                     # 8      fi,IR,t1i,t1r                   \n"  \
"fmul %%st(3)                # 9-13   RI,IR,t1i,t1r                   \n"  \
"fldl %2                     # 10     fr,RI,IR,t1i,t1r                \n"  \
"fmulp %%st,%%st(4)          # 11-15  RI,IR,t1i,RR                    \n"  \
"fldl %3		     # 12     fi,RI,ir,t1i,RR                 \n"  \
"fmulp %%st,%%st(3)          # 13-17  RI,ir,II,RR                     \n"  \
"movl %%edx,%%eax            # 14  U                                  \n"  \
"andl $-512,%%edx            # 14  V                                  \n"  \
"shrl $7,%%edx               # 15  U                                  \n"  \
"addl %%eax,%%eax	     # 15  V                                  \n"  \
"addl %%eax,%%edx	     # 16  U                                  \n"  \
"leal (%%esi,%%edx,8),%%edx                                           \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr                     \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I                     \n"  \
"fldl (%%edx)                 # 18     nt0r,rr,ii,NT1I                \n"  \
"fmull %4                     # 19-23  RR0,rr,ii,NT1I                 \n"  \
"fxch %%st(1)                 # 19     rr,RR0,ii,NT1I                 \n"  \
"fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I                  \n"  \
"fldl 8(%%edx)                # 21     t0i,RR,NT1R,NT1I               \n"  \
"fmull %5                     # 22-26  II,RR,NT1R,nt1i                \n"  \
"fldl (%%edx)                 # 23     t0r,II,RR,NT1R,nt1i            \n"  \
"fmull %5                     # 24-28  RI,II,rr,NT1R,nt1i             \n"  \
"fldl 8(%%edx)                # 25     t0i,RI,II,rr,nt1r,nt1i         \n"  \
"prefetchnta 64(%%edx)                                                \n"  \
"fmull %4                     # 26-30  IR,RI,II,rr,nt1r,nt1i          \n"  \
"fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i          \n"  \
"fsubp %%st,%%st(2)           # 27-31  RI,NTOR,IR,nt1r,nt1i           \n"  \
"fxch %%st(2)                 # 27     IR,NTOR,RI,nt1r,nt1i           \n"  \
"fld  %%st(4)                 # 28     nt1i,IR,NTOR,RI,nt1r,nt1i      \n"  \
"fld  %%st(4)                 # 29     nt1r,nt1i,IR,NTOR,RI,nt1r,nt1i \n"  \
"fadd %%st(3)                 # 30-34  T0R,nt1i,IR,NTOR,ri,nt1r,nt1i  \n"  \
"fxch %%st(4)                 # 30     ri,nt1i,IR,NTOR,TOR,nt1r,nt1i  \n"  \
"faddp %%st,%%st(2)           # 31-35  nt1i,NTOI,ntor,TOR,nt1r,nt1i   \n"  \
"fxch %%st(2)                 # 31     nt0r,NTOI,nt1i,TOR,nt1r,nt1i   \n"  \
"fsubp %%st,%%st(4)           # 32-36  NTOI,nt1i,TOR,T1R,nt1i         \n"  \
"fadd  %%st,%%st(1)           # 33-37  ntoi,TOI,TOR,T1R,nt1i          \n"  \
"fsubp %%st,%%st(4)           # 34-38  TOI,TOR,T1R,T1I                \n"  \
: "=&b"(_pd1) ,"=&d" (_pd0): "m"(_f1##r),"m"(_f1##i),"m"(_f0##r),          \
"m"(_f0##i),"S"(_array),"0"(_k1),"1"(_k0) :"ax","st");                     \
      __asm__ volatile("                                              \n"  \
"fstpl  %1           # 37                                             \n"  \
"fstpl  %0           # 38                                             \n"  \
"fstpl  %2           # 39                                             \n"  \
"fstpl  %3           # 40                                             \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):                \
"memory", "st");

#endif

#ifndef USE_ASM

# define cplx_load_mulmuladdsub_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
	{                                                            \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci;                    \
		_br= *( _pd0 + _index);                              \
		_bi= *( _pd0 + _index + 1);                          \
                prefetch_p( _pd0 + _index);                          \
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i);           \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r);           \
		_br= *( _pd1 + _index);                              \
		_bi= *( _pd1 + _index + 1);                          \
                prefetch_p( _pd1 + _index);                          \
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i);           \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r);           \
		_t0##r = _ar + _cr;                                  \
		_t1##r = _ar - _cr;                                  \
		_t0##i = _ai + _ci;                                  \
		_t1##i = _ai - _ci;                                  \
	}

#else

#define cplx_load_mulmuladdsub_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1)  \
 __asm__ volatile("                                           \n"    \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache        \n"    \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r                \n"    \
"prefetchnta 64(%%esi,%%ebx,8)                                \n"    \
"fldl %0     		     # 6      fr,t1i,t1r              \n"    \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r             \n"    \
"fldl %1                      # 8      fi,IR,t1i,t1r          \n"    \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r          \n"    \
"fldl %0                      # 10     fr,RI,IR,t1i,t1r       \n"    \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR           \n"    \
"fldl %1		      # 12     fi,RI,ir,t1i,RR        \n"    \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR            \n"    \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr             \n"    \
"fxch %%st(2)                 # 17     rr,II,NT1I             \n"    \
"fldl (%%edi,%%ebx,8)         # 18     nt0r,rr,ii,NT1I        \n"    \
"fmull %2                     # 19-23  RR0,rr,ii,NT1I         \n"    \
"fxch %%st(1)                 # 19     rr,RR0,ii,NT1I         \n"    \
"fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I          \n"    \
"fldl 8(%%edi,%%ebx,8)        # 21     t0i,RR,NT1R,NT1I       \n"    \
"fmull %3                     # 22-26  II,RR,NT1R,nt1i        \n"    \
"fldl (%%edi,%%ebx,8)         # 23     t0r,II,RR,NT1R,nt1i    \n"    \
"fmull %3                     # 24-28  RI,II,rr,NT1R,nt1i     \n"    \
"fldl 8(%%edi,%%ebx,8)        # 25     t0i,RI,II,rr,nt1r,nt1i \n"    \
"prefetchnta 64(%%edi,%%ebx,8)                                \n"    \
"fmull %2                     # 26-30  IR,RI,II,rr,nt1r,nt1i  \n"    \
"fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i  \n"    \
"fsubp %%st,%%st(2)           # 27-31  RI,NTOR,IR,nt1r,nt1i   \n"    \
"fxch %%st(2)                 # 27     IR,NTOR,RI,nt1r,nt1i   \n"    \
"fld  %%st(4)                 # 28     nt1i,IR,NTOR,RI,nt1r,nt1i \n" \
"fld  %%st(4)                 # 29     nt1r,nt1i,IR,NTOR,RI,nt1r,nt1i \n"\
"fadd %%st(3)                 # 30-34  T0R,nt1i,IR,NTOR,ri,nt1r,nt1i \n"\
"fxch %%st(4)                 # 30     ri,nt1i,IR,NTOR,TOR,nt1r,nt1i \n"\
"faddp %%st,%%st(2)           # 31-35  nt1i,NTOI,ntor,TOR,nt1r,nt1i \n"\
"fxch %%st(2)                 # 31     nt0r,NTOI,nt1i,TOR,nt1r,nt1i \n"\
"fsubp %%st,%%st(4)           # 32-36  NTOI,nt1i,TOR,T1R,nt1i \n"    \
"fadd  %%st,%%st(1)           # 33-37  ntoi,TOI,TOR,T1R,nt1i  \n"    \
"fsubp %%st,%%st(4)           # 34-38  TOI,TOR,T1R,T1I        \n"    \
	: : "m"(_f1##r),"m"(_f1##i),"m"(_f0##r),                     \
	"m"(_f0##i),"b"(_index),"S"(_pd1),"D"(_pd0) :"st");          \
      __asm__ volatile(" \n"                                         \
"fstpl  %1           # 37    \n"                                     \
"fstpl  %0           # 38    \n"                                     \
"fstpl  %2           # 39    \n"                                     \
"fstpl  %3           # 40    \n"                                     \
        : :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):  \
	      "memory", "st");

#endif

#ifndef USE_ASM

#define cplx_load_mul_p(_t,_pd,_f,_array,_k)    \
{                                               \
	y_limb_t _ar,_ai;                       \
	_pd = _array + addr((_k)<<1);           \
	_ar = *( _pd );                         \
	_ai = *( _pd + 1);                      \
         prefetch_p(_pd);                       \
	_t##r = _ar * _f##r - _ai * _f##i;      \
	_t##i = _ai * _f##r + _ar * _f##i;      \
}

#else

#define cplx_load_mul_p(_t,_pd,_f,_array,_k)              \
 __asm__ volatile("                                     \n"  \
"movl %%edx,%%eax             # 1  U                    \n"  \
"andl $-512,%%edx             # 1  V                    \n"  \
"shrl $7,%%edx                # 2  U                    \n"  \
"addl %%eax,%%eax	     # 2  V                     \n"  \
"addl %%eax,%%edx	     # 3  U                     \n"  \
"fldl (%%esi,%%edx,8)         # 4 (14) if not in cache  \n"  \
"fldl 8(%%esi,%%edx,8)        # 5      t1i,t1r          \n"  \
"leal (%%esi,%%edx,8),%%edx                             \n"  \
"fldl %1     		     # 6      fr,t1i,t1r        \n"  \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r       \n"  \
"fldl %2                      # 8      fi,IR,t1i,t1r    \n"  \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r    \n"  \
"fldl %1                      # 10     fr,RI,IR,t1i,t1r \n"  \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR     \n"  \
"fldl %2		     # 12     fi,RI,ir,t1i,RR   \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR      \n"  \
"prefetchnta 64(%%edx)                                  \n"  \
: "=&d" (_pd): "m"(_f##r),"m"(_f##i),"S"(_array),"0"(_k):    \
"ax","st");                                                  \
__asm__ volatile ("                                     \n"  \
"faddp %%st,%%st(1)           # 17-19  NT1I,ii,rr        \n" \
"fxch %%st(2)                 # 17     rr,II,NT1I        \n" \
"fsubp %%st,%%st(1)           # 18-20  NT1r,NT1I         \n" \
"fxch %%st(1)		     # 19     NT1I,NT1R         \n"  \
"fstpl %1                     # 20     NT1R              \n" \
"fstpl %0                     # 21     NT1I              \n" \
: :"m"( _t##r ),"m" ( _t##i ):"memory", "st");
#endif

#ifndef USE_ASM

#define cplx_load_mul_pp(_t,_f,_index,_pd)    \
{                                             \
	y_limb_t _ar,_ai;                     \
	_ar = *( _pd + _index );              \
	_ai = *( _pd + _index + 1);           \
         prefetch_p( _pd + _index);           \
	_t##r = _ar * _f##r - _ai * _f##i;    \
	_t##i = _ai * _f##r + _ar * _f##i;    \
}

#else

#define cplx_load_mul_pp(_t,_f,_index,_pd)                        \
 __asm__ volatile("                                          \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache       \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r               \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                               \n"  \
"fldl %0     		     # 6      fr,t1i,t1r             \n"  \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r            \n"  \
"fldl %1                      # 8      fi,IR,t1i,t1r         \n"  \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r         \n"  \
"fldl %0                      # 10     fr,RI,IR,t1i,t1r      \n"  \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR          \n"  \
"fldl %1			     # 12     fi,RI,ir,t1i,RR\n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR           \n"  \
: : "m"(_f##r),"m"(_f##i),                                        \
"S"(_pd),"b"(_index):"st");                                       \
__asm__ volatile ("                                          \n"  \
"faddp %%st,%%st(1)           # 17-19  NT1I,ii,rr            \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I            \n"  \
"fsubp %%st,%%st(1)           # 18-20  NT1r,NT1I             \n"  \
"fxch %%st(1)		     # 19     NT1I,NT1R              \n"  \
"fstpl %1                     # 20     NT1R                  \n"  \
"fstpl %0                     # 21     NT1I                  \n"  \
: :"m"( _t##r ),"m" ( _t##i ):"memory", "st");
#endif

#ifndef USE_ASM

#define cplx_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{                                                            \
		y_limb_t _ar,_ai,_b,_cr,_ci;                         \
		_pd0 = _array + addr((_k0)<<1);                      \
		_b = *( _pd0 );                                      \
		_ai= *( _pd0 + 1);                                   \
                prefetch_p(_pd0);                                    \
		_ar = (_b) * (_f0##r) - (_ai) * (_f0##i);            \
		_pd1 = _array + addr((_k1)<<1);                      \
		_ai = (_b) * (_f0##i) + (_ai) * (_f0##r);            \
		_b = *( _pd1 );                                      \
		_ci = *( _pd1 + 1);                                  \
                prefetch_p(_pd1);                                    \
		_cr = (_b) * (_f1##r) - (_ci) * (_f1##i);            \
		_ci = (_b) * (_f1##i) + (_ci) * (_f1##r);            \
		_t0##r = _ar;                                        \
		_t0##i = _ai;                                        \
		_t1##r = _cr;                                        \
		_t1##i = _ci;                                        \
	}

#else

#define cplx_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
 __asm__ volatile ("                                            \n"  \
"movl %%ebx,%%eax             # 1  U                            \n"  \
"andl $-512,%%ebx             # 1  V                            \n"  \
"shrl $7,%%ebx                # 2  U                            \n"  \
"addl %%eax,%%eax	     # 2  V                             \n"  \
"addl %%eax,%%ebx	     # 3  U                             \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache          \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r                  \n"  \
"leal (%%esi,%%ebx,8),%%ebx                                     \n"  \
"prefetchnta 64(%%ebx)                                          \n"  \
"fldl %2     		     # 6      fr,t1i,t1r                \n"  \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r               \n"  \
"fldl %3                      # 8      fi,IR,t1i,t1r            \n"  \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r            \n"  \
"fldl %2                      # 10     fr,RI,IR,t1i,t1r         \n"  \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR             \n"  \
"fldl %3			     # 12     fi,RI,ir,t1i,RR   \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR              \n"  \
"movl %%edx,%%eax             # 14  U                           \n"  \
"andl $-512,%%edx             # 14  V                           \n"  \
"shrl $7,%%edx                # 15  U                           \n"  \
"addl %%eax,%%eax	     # 15  V                            \n"  \
"addl %%eax,%%edx	     # 16  U                            \n"  \
"leal (%%esi,%%edx,8),%%edx                                     \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr               \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I               \n"  \
"fldl (%%edx)                 # 18     nt0r,rr,ii,NT1I          \n"  \
"fmull %4                     # 19-23  RR0,rr,ii,NT1I           \n"  \
"fxch %%st(1)                 # 19     rr,RR0,ii,NT1I           \n"  \
"fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I            \n"  \
"fldl 8(%%edx)                # 21     t0i,RR,NT1R,NT1I         \n"  \
"fmull %5                     # 22-26  II,RR,NT1R,nt1i          \n"  \
"fldl (%%edx)                 # 23     t0r,II,RR,NT1R,nt1i      \n"  \
"fmull %5                     # 24-28  RI,II,rr,NT1R,nt1i       \n"  \
"fldl 8(%%edx)                # 25     t0i,RI,II,rr,nt1r,nt1i   \n"  \
"fmull %4                     # 26-30  IR,RI,II,rr,nt1r,nt1i    \n"  \
"fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i    \n"  \
"fsubp %%st,%%st(2)           # 27-31  RI,NT0R,IR,nt1r,nt1i     \n"  \
"prefetchnta 64(%%edx)                                          \n"  \
: "=&b"(_pd1) ,"=&d" (_pd0): "m"(_f1##r),"m"(_f1##i),"m"(_f0##r),    \
"m"(_f0##i),"S"(_array),"0"(_k1),"1"(_k0) :"ax","st");               \
__asm__ volatile("                                              \n"  \
"faddp  %%st,%%st(2)  	# 30-34  NT0R,NT0I,nt1r,nt1i            \n"  \
"fxch   %%st(3)          # 31     nt1i,NT0I,nt1r,NT0R           \n"  \
"fstpl  %3           	# 32     NT0I,nt1r,nt0r                 \n"  \
"fxch   %%st(2)          # 33     nt0r,nt1r,NT0I                \n"  \
"fstpl  %0           	# 34     nt1r,NT0I                      \n"  \
"fstpl  %2           	# 35     nt0i                           \n"  \
"fstpl  %1           	# 36                                    \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):          \
"memory", "st");

#endif

#ifndef USE_ASM

#define cplx_load_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1)        \
	{                                                            \
		y_limb_t _ar,_ai,_b,_cr,_ci;                         \
		_b = *( _pd0 + _index );                             \
		_ai= *( _pd0 + _index + 1);                          \
                prefetch_p( _pd0 + _index);                          \
		_ar = (_b) * (_f0##r) - (_ai) * (_f0##i);            \
		_ai = (_b) * (_f0##i) + (_ai) * (_f0##r);            \
		_b = *( _pd1 + _index );                             \
		_ci = *( _pd1 + _index + 1);                         \
                prefetch_p( _pd1 + _index);                          \
		_cr = (_b) * (_f1##r) - (_ci) * (_f1##i);            \
		_ci = (_b) * (_f1##i) + (_ci) * (_f1##r);            \
		_t0##r = _ar;                                        \
		_t0##i = _ai;                                        \
		_t1##r = _cr;                                        \
		_t1##i = _ci;                                        \
	}

#else

#define cplx_load_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1)        \
 __asm__ volatile ("                                            \n"  \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache          \n"  \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r                  \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                                  \n"  \
"fldl %0     		     # 6      fr,t1i,t1r                \n"  \
"fmul %%st(1)                 # 7-11   IR,t1i,t1r               \n"  \
"fldl %1                      # 8      fi,IR,t1i,t1r            \n"  \
"fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r            \n"  \
"fldl %0                      # 10     fr,RI,IR,t1i,t1r         \n"  \
"fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR             \n"  \
"fldl %1		      # 12     fi,RI,ir,t1i,RR          \n"  \
"fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR              \n"  \
"faddp %%st,%%st(1)           # 17-21  NT1I,II,rr               \n"  \
"fxch %%st(2)                 # 17     rr,II,NT1I               \n"  \
"fldl (%%edi,%%ebx,8)         # 18     nt0r,rr,ii,NT1I          \n"  \
"fmull %2                     # 19-23  RR0,rr,ii,NT1I           \n"  \
"fxch %%st(1)                 # 19     rr,RR0,ii,NT1I           \n"  \
"fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I            \n"  \
"fldl 8(%%edi,%%ebx,8)        # 21     t0i,RR,NT1R,NT1I         \n"  \
"fmull %3                     # 22-26  II,RR,NT1R,nt1i          \n"  \
"fldl (%%edi,%%ebx,8)         # 23     t0r,II,RR,NT1R,nt1i      \n"  \
"fmull %3                     # 24-28  RI,II,rr,NT1R,nt1i       \n"  \
"fldl 8(%%edi,%%ebx,8)        # 25     t0i,RI,II,rr,nt1r,nt1i   \n"  \
"fmull %2                     # 26-30  IR,RI,II,rr,nt1r,nt1i    \n"  \
"fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i    \n"  \
"fsubp %%st,%%st(2)           # 27-31  RI,NT0R,IR,nt1r,nt1i     \n"  \
"prefetchnta 64(%%edi,%%ebx,8)                                  \n"  \
: : "m"(_f1##r),"m"(_f1##i),"m"(_f0##r),                             \
"m"(_f0##i),"b"(_index),"S"(_pd1),"D"(_pd0) :"st");                  \
__asm__ volatile("                                              \n"  \
"faddp  %%st,%%st(2)  	# 30-34  NT0R,NT0I,nt1r,nt1i            \n"  \
"fxch   %%st(3)          # 31     nt1i,NT0I,nt1r,NT0R           \n"  \
"fstpl  %3           	# 32     NT0I,nt1r,nt0r                 \n"  \
"fxch   %%st(2)          # 33     nt0r,nt1r,NT0I                \n"  \
"fstpl  %0           	# 34     nt1r,NT0I                      \n"  \
"fstpl  %2           	# 35     nt0i                           \n"  \
"fstpl  %1           	# 36                                    \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):          \
"memory", "st");

#endif

/**********************************************************************
   cplx_load_addsub macro:
		_t0 = _array( _k0) + _array (_k1);
		_t1 = _array( _k0) - _array (_k1);
 where all the variables are pseudo-complex .
**********************************************************************/
#ifndef USE_ASM
#define cplx_load_addsub(_t0,_t1,_array,_k0,_k1)      \
	{                                             \
		y_limb_t _ar,_ai;                     \
		_ar = _array[addr(_k1<<1)];           \
		_ai = _array[addr(_k1<<1)+1];         \
                prefetch_p(&_array[addr(_k1<<1)]);    \
		_t0##r = _array[addr(_k0<<1)] + _ar;  \
		_t1##r = _array[addr(_k0<<1)] - _ar;  \
		_t0##i = _array[addr(_k0<<1)+1] + _ai;\
		_t1##i = _array[addr(_k0<<1)+1] - _ai;\
                prefetch_p(&_array[addr(_k0<<1)]);    \
	}

#else

#define cplx_load_addsub(_t0,_t1,_array,_k0,_k1)           \
	__asm__ volatile( "                           \n"  \
"movl %1,%%ebx                                        \n"  \
"movl %2,%%edx                                        \n"  \
"movl %%ebx,%%eax             # 1  U                  \n"  \
"andl $-512,%%ebx             # 1  V                  \n"  \
"shrl $7,%%ebx                # 2  U                  \n"  \
"addl %%eax,%%eax	     # 2  V                   \n"  \
"addl %%eax,%%ebx	     # 3  U                   \n"  \
"fldl (%%esi,%%ebx,8)         # 4  t0r                \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                        \n"  \
"movl %%edx,%%eax             # 5  U                  \n"  \
"andl $-512,%%edx             # 5  V                  \n"  \
"shrl $7,%%edx                # 6  U                  \n"  \
"addl %%eax,%%eax	     # 6  V                   \n"  \
"addl %%eax,%%edx	     # 7  U                   \n"  \
"faddl (%%esi,%%edx,8)        # 8-12   NT0R           \n"  \
"fldl  (%%esi,%%ebx,8)        # 9      t0r,NT0R       \n"  \
"fsubl (%%esi,%%edx,8)        # 10-14  NT1R,NT0R      \n"  \
"fldl  8(%%esi,%%ebx,8)       # 11     t0i,NT1R,NT0R  \n"  \
"faddl 8(%%esi,%%edx,8)       # 12-16  NT0I,NT1R,ntor \n"  \
"fxch %%st(2)                 # 12     nt0r,NT1R,NT0I \n"  \
"fldl 8(%%esi,%%ebx,8)        # 13     t0i,nt0r,NT1R,NT0I\n"\
"fsubl 8(%%esi,%%edx,8)       # 14-18  NT1I,nt0r,NT1R,NTOI\n"\
"fxch %%st(3)                 # 15     NTOI,nt0r,nt1R,NT1I\n"\
"prefetchnta 64(%%esi,%%edx,8)                        \n"  \
: :"S"(_array),"m"(_k0),"m"(_k1):"ax","bx","dx","st");     \
      __asm__ volatile("                              \n"  \
"fstpl  %1                    # 18                    \n"  \
"fstpl  %0                    # 19                    \n"  \
"fstpl  %2                    # 20                    \n"  \
"fstpl  %3                    # 21                    \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):\
"memory", "st");
#endif

#ifndef USE_ASM

#define cplx_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1) \
	{                                                    \
		y_limb_t _ar,_ai;                            \
		_pd1 = _array + addr((_k1)<<1);              \
		_ar = *( _pd1 );                             \
		_ai = *( _pd1 + 1);                          \
                prefetch_p(_pd1);                            \
		_pd0 = _array + addr((_k0)<<1);              \
                prefetch_p(_pd0);                            \
		_t0##r = *(_pd0 ) + _ar;                     \
		_t1##r = *(_pd0 ) - _ar;                     \
		_t0##i = *(_pd0 + 1) + _ai;                  \
		_t1##i = *(_pd0 + 1) - _ai;                  \
	}

#else

#define cplx_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1) \
	__asm__ volatile( "                             \n"  \
"movl %%ebx,%%eax             # 1  U                    \n"  \
"andl $-512,%%ebx             # 1  V                    \n"  \
"shrl $7,%%ebx                # 2  U                    \n"  \
"addl %%eax,%%eax	     # 2  V                     \n"  \
"addl %%eax,%%ebx	     # 3  U                     \n"  \
"fldl (%%esi,%%ebx,8)         # 4  t0r                  \n"  \
"leal (%%esi,%%ebx,8),%%ebx                             \n"  \
"prefetchnta 64(%%ebx)                                  \n"  \
"movl %%edx,%%eax             # 5  U                    \n"  \
"andl $-512,%%edx             # 5  V                    \n"  \
"shrl $7,%%edx                # 6  U                    \n"  \
"addl %%eax,%%eax	     # 6  V                     \n"  \
"addl %%eax,%%edx	     # 7  U                     \n"  \
"leal (%%esi,%%edx,8),%%edx                             \n"  \
"faddl (%%edx)                # 8-12   NT0R             \n"  \
"fldl  (%%ebx)                # 9      t0r,NT0R         \n"  \
"fsubl (%%edx)                # 10-14  NT1R,NT0R        \n"  \
"fldl  8(%%ebx)               # 11     t0i,NT1R,NT0R    \n"  \
"faddl 8(%%edx)               # 12-16  NT0I,NT1R,ntor   \n"  \
"fxch %%st(2)                 # 12     nt0r,NT1R,NT0I   \n"  \
"fldl 8(%%ebx)                # 13     t0i,nt0r,NT1R,NT0I\n" \
"fsubl 8(%%edx)               # 14-18  NT1I,nt0r,NT1R,NTOI\n"\
"fxch %%st(3)                 # 15     NTOI,nt0r,nt1R,NT1I\n"\
"prefetchnta 64(%%edx)                                  \n"  \
:"=&b"(_pd0),"=&d"(_pd1):"S"(_array),"0"(_k0),"1"(_k1):      \
"ax","st");                                                  \
     __asm__ volatile("                                 \n"  \
"fstpl  %1                    # 18                      \n"  \
"fstpl  %0                    # 19                      \n"  \
"fstpl  %2                    # 20                      \n"  \
"fstpl  %3                    # 21                      \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):  \
"memory", "st");
#endif

#ifndef USE_ASM

#define cplx_load_addsub_pp(_t0,_t1,_index,_pd0,_pd1)        \
	{                                                    \
		y_limb_t _ar,_ai;                            \
		_ar = *( _pd1 + _index);                     \
		_ai = *( _pd1 + _index + 1);                 \
                prefetch_p( _pd1 + _index);                  \
		_t0##r = *( _pd0 + _index) + _ar;            \
		_t1##r = *( _pd0 + _index) - _ar;            \
		_t0##i = *( _pd0 + _index + 1) + _ai;        \
		_t1##i = *( _pd0 + _index + 1) - _ai;        \
                prefetch_p( _pd0 + _index);                  \
	}

#else

#define cplx_load_addsub_pp(_t0,_t1,_index,_pd0,_pd1)          \
	__asm__ volatile( "                               \n"  \
"fldl (%%esi,%%ebx,8)         # 4  t0r                    \n"  \
"faddl (%%edi,%%ebx,8)        # 8-12   NT0R               \n"  \
"fldl  (%%esi,%%ebx,8)        # 9      t0r,NT0R           \n"  \
"fsubl (%%edi,%%ebx,8)        # 10-14  NT1R,NT0R          \n"  \
"fldl  8(%%esi,%%ebx,8)       # 11     t0i,NT1R,NT0R      \n"  \
"faddl 8(%%edi,%%ebx,8)       # 12-16  NT0I,NT1R,ntor     \n"  \
"prefetchnta 64(%%esi,%%ebx,8)                            \n"  \
"fxch %%st(2)                 # 12     nt0r,NT1R,NT0I     \n"  \
"fldl 8(%%esi,%%ebx,8)        # 13     t0i,nt0r,NT1R,NT0I \n"  \
"fsubl 8(%%edi,%%ebx,8)       # 14-18  NT1I,nt0r,NT1R,NTOI \n" \
"fxch %%st(3)                 # 15     NTOI,nt0r,nt1R,NT1I \n" \
"prefetchnta 64(%%edi,%%ebx,8)                            \n"  \
::"b"(_index),"S"(_pd0),"D"(_pd1):"st")                       ;\
      __asm__ volatile("                                  \n"  \
"fstpl  %1                    # 18                        \n"  \
"fstpl  %0                    # 19                        \n"  \
"fstpl  %2                    # 20                        \n"  \
"fstpl  %3                    # 21                        \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):    \
"memory", "st");
#endif


/*$Id$*/










