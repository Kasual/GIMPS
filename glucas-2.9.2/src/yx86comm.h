/* $Id$ */
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

/* THIS FILE INCLUDES ALL THE ASSEMBLER MACROS COMMON TO X86 PLATFORMS
   NOT INCLUDED IN FF3X86.H FFT5X86.H FFT7X86.H                         */

#if defined(USE_ASM) && defined(__CYGWIN__)
/* Cygwin's GCC uses the variable exactly as given in asm, but exports with _ prefix in C */
#define MEMORY_VARIABLE_NAME(name) "_" #name
#else
#define MEMORY_VARIABLE_NAME(name) #name
#endif


#define cplx_addsub_store(_k0,_k1,_t0,_t1)    \
	_d[addr(_k0<<1)] = _t0##r + _t1##r;   \
	_d[addr(_k0<<1)+1] = _t0##i + _t1##i; \
	_d[addr(_k1<<1)] = _t0##r - _t1##r;   \
	_d[addr(_k1<<1)+1] = _t0##i - _t1##i;


#define cplx_addsub_store_p(_index,_pd0,_pd1,_t0,_t1)\
	*( _pd0 + _index) = _t0##r + _t1##r;         \
	*( _pd0 + _index + 1) = _t0##i + _t1##i;     \
	*( _pd1 + _index) = _t0##r - _t1##r;         \
	*( _pd1 + _index + 1) = _t0##i - _t1##i;


/* To store on read padded array from local data as complex */

#ifndef USE_ASM

#define cplx_local_to_data(_array,_k,_t) \
	_array[addr(( _k )<<1)]=_t##r ;  \
	_array[addr(( _k )<<1)+1]=_t##i;

#else

#define cplx_local_to_data(_array,_k,_t) \
	__asm__("                   \n"  \
"movl %%ebx,%%eax        #1  U      \n"  \
"andl $-512,%%ebx        #1  V      \n"  \
"shrl $7,%%ebx           #2  U      \n"  \
"addl %%eax,%%eax	#2  V       \n"  \
"addl %%eax,%%ebx	#3  U       \n"  \
"fldl %0                 #4         \n"  \
"fldl %1                 #5         \n"  \
"fxch %%st(1)            #6         \n"  \
"fstpl (%%esi,%%ebx,8)   #7         \n"  \
"fstpl 8(%%esi,%%ebx,8)  #8         \n"  \
: :"m"(_t##r),"m"(_t##i),"S"(_array),    \
"b"(_k):"ax","memory");

#endif

#ifndef USE_ASM

#define cplx_local_to_data_p(_pd,_t) \
	*( _pd ) = _t##r ;           \
	*( _pd + 1) = _t##i;

#else

#define cplx_local_to_data_p(_pd,_t) \
	__asm__("               \n"  \
"fldl %0                   #1   \n"  \
"fldl %1		   #2   \n"  \
"fxch %%st(1)              #3   \n"  \
"fstpl (%%esi)             #4   \n"  \
"fstpl 8(%%esi)            #5   \n"  \
: :"m"(_t##r),"m"(_t##i),"S"(_pd):   \
"memory");

#endif

#ifndef USE_ASM

/* To load from real non-padded array of data as a complex array */
#define cplx_mem_to_local(_t,_array,_k) \
	_t##r =_array[( _k )<<1];       \
	_t##i =_array[(( _k )<<1)+1];

#else

#define cplx_mem_to_local(_t,_array,_k) \
	__asm__("                \n"    \
"addl %%ebx,%%ebx           #1   \n"    \
"fldl (%%esi,%%ebx,8)       #2   \n"    \
"fldl 8(%%esi,%%ebx,8)      #2   \n"    \
"fxch %%st(1)               #3   \n"    \
"fstpl %0                   #4   \n"    \
"fstpl %1                   #5   \n"    \
: :"m"(_t##r),"m"(_t##i),"b"(_k),       \
"S"(_array): "memory");

#endif

#ifndef USE_ASM

/* To strore to a non-padded array from local */
#define cplx_local_to_mem(_array,_k,_t) \
	_array[( _k )<<1]=_t##r;        \
	_array[(( _k )<<1)+1]=_t##i;

#else

#define cplx_local_to_mem(_array,_k,_t) \
	__asm__("                \n"     \
"fldl %0                    #1   \n"     \
"fldl %1                    #2   \n"     \
"addl %ebx,%ebx             #3   \n"     \
"fstpl 8(%%esi,%%ebx,8)     #4   \n"     \
"fstpl (%%esi,%%ebx,8)      #5   \n"     \
: :"m"(_t##r),"m"(_t##i),"b"(_k),        \
"S"(_array):"memory");

#endif

/**********************************************************************
   cplx_muladdsub macro :
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex and _t is temporal aux var.
***********************************************************************/

#ifndef USE_ASM

# define cplx_muladdsub(_t0,_t1,_f)                                \
	{                                                          \
		y_limb_t _ar,_ai,_a=_t0##r;                        \
		_ar = (_t1##r ) * (_f##r ) - (_t1##i ) * (_f##i ); \
		_ai = (_t1##r ) * (_f##i ) + (_t1##i ) * (_f##r ); \
		_t0##r += _ar;                                     \
		_t1##r = _a - _ar;                                 \
		_a = _t0##i;                                       \
		_t0##i += _ai;                                     \
		_t1##i = _a - _ai;                                 \
	}

#else

# define cplx_muladdsub(_t0,_t1,_f)                             \
   __asm__ volatile( "                                    \n"   \
"fldl %2                      #1     t1r                  \n"   \
"fmull %4                     #2-6   RR                   \n"   \
"fldl %3                      #3     t1i,RR               \n"   \
"fmull %5                     #4-8   II,RR                \n"   \
"fldl %2                      #5     t1r,II,RR            \n"   \
"fmull %5                     #6-10  RI,II,RR             \n"   \
"fldl %3                      #7     t1i,RI,II,rr         \n"   \
"fmull %4                     #8-12  IR,RI,II,rr          \n"   \
"fxch %%st(3)                 #8     rr,RI,II,IR          \n"   \
"fsubp %%st,%%st(2)           #9-12  RI,R1,IR             \n"   \
"fldl %0                      #10    t0r,RI,R1,IR         \n"   \
"fldl %0                      #11    t0r,t0r,ri,R1,IR     \n"   \
"fxch %%st(2)                 #11    ri,t0r,t0r,R1,IR     \n"   \
"faddp %%st,%%st(4)           #12-16 t0r,t0r,R1,I1        \n"   \
"fadd %%st(2)                 #13-17 NT0R,t0r,r1,I1       \n"   \
"fxch %%st(1)                 #13    t0r,NT0R,r1,I1       \n"   \
"fsubp %%st,%%st(2)           #14-18 NT0R,NT1R,I1         \n"   \
"fldl %1                      #15    t0i,NT0R,NT1R,I1     \n"   \
"fadd %%st(3)                 #16-20 NT0I,NT0R,NT1R,I1    \n"   \
"fldl %1                      #17    t0i,NT0I,NT0R,NT1R,I1\n"   \
"fsubp %%st,%%st(4)           #18-22 NT0I,nt0r,NT1R,NT1I  \n"   \
"fxch %%st(2)                 #19    nt1r,nt0r,NT0I,NT1I  \n"   \
"fstpl %2                     #20    nt0r,NT0I,NT1I       \n"   \
"fstpl %0                     #21    nt0i,NT1I            \n"   \
"fstpl %1                     #22    nt1i                 \n"   \
"fstpl %3                     #23                         \n"   \
: :"m"( _t0##r ),"m" ( _t0##i ), "m" (_t1##r), "m" (_t1##i),    \
   "m"( _f##r ), "m" ( _f##i ): "memory");

#endif

#ifndef USE_ASM

#define cplx_load_muladdsub_pp_no_fetch(_t0,_t1,_f,_index,_pd0,_pd1) \
	{                                                            \
		y_limb_t _ar,_ai,_br,_bi;                            \
		_br= *( _pd1 + _index );                             \
		_bi= *( _pd1 + _index + 1);                          \
		_ar = (_br) * (_f##r) - (_bi) * (_f##i);             \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r);             \
		_t0##r = *( _pd0 + _index) + _ar;                    \
		_t1##r = *( _pd0 + _index) - _ar;                    \
		_t0##i = *( _pd0 + _index + 1) + _ai;                \
		_t1##i = *( _pd0 + _index + 1) - _ai;                \
	}

#else

#define cplx_load_muladdsub_pp_no_fetch(_t0,_t1,_f,_index,_pd0,_pd1) \
__asm__ volatile("                                            \n"    \
"fldl (%%esi,%%ebx,8)         # 1 (a lot if not in cache)     \n"    \
"fldl 8(%%esi,%%ebx,8)        # 2 t1i,t1r                     \n"    \
"fldl %0     		     # 4 fr,t1i,t1r                   \n"    \
"fmul %%st(1)                 # 5-9  IR,t1i,t1r               \n"    \
"fldl %1                      # 6    fi,IR,t1i,t1r            \n"    \
"fmul %%st(3)                 # 7-11 RI,IR,t1i,t1r            \n"    \
"fldl %0                      # 8    fr,RI,IR,t1i,t1r         \n"    \
"fmulp %%st,%%st(4)           # 9-13   RI,IR,t1i,RR           \n"    \
"fldl %1			     # 10     fi,RI,ir,t1i,RR \n"    \
"fmulp %%st,%%st(3)           # 11-15  RI,ir,II,RR            \n"    \
"faddp %%st,%%st(1)           # 12-16  NT1I,II,rr             \n"    \
"fxch %%st(2)                 # 12     rr,II,NT1I             \n"    \
"fldl (%%edi,%%ebx,8)         # 13     nt0r,rr,II,NT1I        \n"    \
"fxch %%st(1)                 # 13     rr,nt0r,II,NT1I        \n"    \
"fldl 8(%%edi,%%ebx,8)        # 14     ntoi,rr,ntor,II,NT1I   \n"    \
"fxch %%st(1)                 # 14     rr,nt0i,nt0r,II,NT1I   \n"    \
"fsubp %%st,%%st(3)           # 15-19  nt0i,nt0r,NT1R,NT1I    \n"    \
"fsub %%st(3)                 # 16-20  T1I,nt0r,NT1R,nt1i     \n"    \
"fxch %%st(3)                 # 16     nt1i,nt0r,NT1R,T1I     \n"    \
"faddl 8(%%edi,%%ebx,8)       # 17-21  T0I,nt0r,nt1r,T1I      \n"    \
"fxch %%st(1)                 # 17     nt0r,T0I,nt1r,T1I      \n"    \
"fsub %%st(2)                 # 19-23  T1R,T0I,nt1r,T1I       \n"    \
"fxch %%st(2)                 # 20     nt1r,T0I,T1R,T1I       \n"    \
"faddl (%%edi,%%ebx,8)        # 21-25  T0R,T0I,T1R,t1i        \n"    \
"fxch %%st(3)                 # 22     t1i,T0I,T1R,T0R        \n"    \
: : "m"(_f##r),"m"(_f##i),"b"(_index),"S"(_pd1),"D"(_pd0) );         \
      __asm__ volatile("                                      \n"    \
"fstpl  %3           # 23                                     \n"    \
"fstpl  %1           # 24                                     \n"    \
"fstpl  %2           # 25                                     \n"    \
"fstpl  %0           # 26                                     \n"    \
: :"m"( _t0##r ),"m" ( _t0##i ),"m" (_t1##r), "m" (_t1##i):          \
"st","memory");

#endif


/**********************************************************************
   cplx_muladdsub_store macro:
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
		_array( _k0) = _t0;
		_array( _k1) = _t1;
   where all the variables are pseudo-complex and _t is temporal aux var.
**********************************************************************/
#define cplx_muladdsub_store(_array,_k0,_k1,_t0,_t1,_f)        \
	{                                                      \
		y_limb_t _ar,_ai;                              \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i); \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r); \
		_array[addr(_k0<<1)]= _t0##r + _ar;            \
                _array[addr(_k1<<1)]= _t0##r - _ar;            \
		_array[addr(_k0<<1)+1]= _t0##i + _ai;          \
                _array[addr(_k1<<1)+1]= _t0##i - _ai;          \
        }

#ifndef USE_ASM

#define cplx_muladdsub_store_p(_index,_pd0,_pd1,_t0,_t1,_f)        \
	{                                                          \
		y_limb_t _ar,_ai;                                  \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i);     \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r);     \
		*( _pd0 + _index) = _t0##r + _ar;                  \
		*( _pd1 + _index) = _t0##r - _ar;                  \
		*( _pd0 + _index + 1) = _t0##i + _ai;              \
		*( _pd1 + _index + 1) = _t0##i - _ai;              \
	}

#else

#define cplx_muladdsub_store_p(_index,_pd0,_pd1,_t0,_t1,_f)        \
__asm__ volatile("                                          \n"    \
"fldl %7                      # 1      t1r                  \n"    \
"fldl %8                      # 2      t1i,t1r              \n"    \
"fldl %2     		     # 3      fr,t1i,t1r            \n"    \
"fmul %%st(1)                 # 4-8    IR,t1i,t1r           \n"    \
"fldl %3                      # 5      fi,IR,t1i,t1r        \n"    \
"fmul %%st(3)                 # 6-10   RI,IR,t1i,t1r        \n"    \
"fldl %2                      # 7      fr,RI,IR,t1i,t1r     \n"    \
"fmulp %%st,%%st(4)           # 8-12   RI,IR,t1i,RR         \n"    \
"fldl %3			     # 13     fi,RI,ir,t1i,RR\n"   \
"fmulp %%st,%%st(3)           # 14-18  RI,ir,II,RR          \n"    \
"faddp %%st,%%st(1)           # 15-19  NT1I,II,rr           \n"    \
"fxch %%st(2)                 # 15     rr,II,NT1I           \n"    \
"fldl %5                      # 16     nt0r,rr,ii,NT1I      \n"    \
"fxch %%st(1)                 # 16     rr,nt0r,ii,NT1I      \n"    \
"fsubp %%st,%%st(2)           # 17-21  nt0r,NT1R,NT1I       \n"    \
"fldl %6                      # 18     nt0i,nt0r,NT1R,NT1I  \n"    \
"fsub %%st(3)                 # 19-23  T1I,nt0r,NT1R,nt1i   \n"    \
"fxch %%st(3)                 # 19     nt1i,nt0r,NT1R,T1I   \n"    \
"faddl %6                     # 20-24  T0I,nt0r,nt1r,T1I    \n"    \
"fxch %%st(1)                 # 20     nt0r,T0I,nt1r,T1I    \n"    \
"fsub %%st(2)                 # 21-25  T1R,T0I,nt1r,T1I     \n"    \
"fxch %%st(2)                 # 21     nt1r,T0I,T1R,T1I     \n"    \
"faddl %5                     # 22-26  T0R,T0I,T1R,t1i      \n"    \
"fxch %%st(3)                 # 23     t1i,T0I,T1R,T0R      \n"    \
"fstpl 8(%%esi,%%ebx,8)       # 24                          \n"    \
"fstpl 8(%%edi,%%ebx,8)       # 25                          \n"    \
"fstpl (%%esi,%%ebx,8)        # 26                          \n"    \
"fstpl (%%edi,%%ebx,8)        # 27                          \n"    \
: :"S"(_pd1) ,"D" (_pd0), "m"(_f##r),"m"(_f##i),"b"(_index),       \
"m"(_t0##r),"m"(_t0##i),"m"(_t1##r),"m"(_t1##i) :"memory");

#endif


/**********************************************************************
   cplx_mulmuladdsub macro:
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/
#ifndef USE_ASM

# define cplx_mulmuladdsub(_t0,_t1,_f0,_f1)                      \
	{                                                        \
	        y_limb_t _ar,_ai,_br,_bi;                        \
		_ar = (_t1##r) * (_f1##r) - (_t1##i) * (_f1##i); \
		_br = (_t0##r) * (_f0##r) - (_t0##i) * (_f0##i); \
		_ai = (_t1##r) * (_f1##i) + (_t1##i) * (_f1##r); \
		_bi = (_t0##r) * (_f0##i) + (_t0##i) * (_f0##r); \
                _t0##r = _br + _ar;                              \
                _t1##r = _br - _ar;                              \
                _t0##i = _bi + _ai;                              \
		_t1##i = _bi - _ai;                              \
	}

#else

#define cplx_mulmuladdsub(_t0,_t1,_f0,_f1)                       \
   __asm__ volatile( "                                      \n"  \
"fldl %2                      #1-1 t1r                      \n"  \
"fmull %6                     #2-6 rr                       \n"  \
"fldl %3                      #3-3 t1i,rr                   \n"  \
"fmull %7                     #4-8 ii,rr                    \n"  \
"fldl %2                      #5-5 t1r,ii,rr                \n"  \
"fmull %7                     #6-10 ri,ii,rr                \n"  \
"fldl %3                      #7-7  t1i,ri,ii,rr            \n"  \
"fmull %6                     #8-12  ir,ri,ii,rr            \n"  \
"fxch %%st(3)                 #8-8   rr,ri,ii,ir            \n"  \
"fsubp %%st,%%st(2)           #9-12  ri,r1,ir               \n"  \
"fldl %0                      #10-10  t0r,ri,r1,ir          \n"  \
"fmull %4                     #11-15  rr0,ri,r1,ir          \n"  \
"fldl %1                      #12-12  t0i,rr0,ri,r1,ir      \n"  \
"fmull %5                     #13-17  ii0,rr0,ri,r1,ir      \n"  \
"fxch %%st(2)                 #13-13  ri,rr0,ii0,r1,ir      \n"  \
"faddp %%st,%%st(4)           #14-17  rr0,ii0,r1,i1         \n"  \
"fldl %0                      #15-15  t0r,rr0,ii0,r1,i1     \n"  \
"fmull %5                     #16-20  ri0,rr0,ii0,r1,i1     \n"  \
"fldl %1                      #17-17  t0i,ri0,rr0,ii0,r1,i1 \n"  \
"fmull %4                     #18-22  ir0,ri0,rr0,ii0,r1,i1 \n"  \
"fxch %%st(2)                 #18-18  rr0,ri0,ir0,ii0,r1,i1 \n"  \
"fsubp %%st,%%st(3)           #19-22  ri0,ir0,r0,r1,i1      \n"  \
"fld %%st(3)                  #20-20  r1,ri0,ir0,r0,r1,i1   \n"  \
"fld %%st(5)                  #21-21  i1,r1,ri0,ir0,r0,r1,i1\n"  \
"fxch %%st(4)                 #21-21  r0,r1,ri0,ir0,i1,r1,i1\n"  \
"fadd %%st,%%st(5)            #22-25  r0,r1,ri0,ir0,i1,nt0r,i1\n"\
"fxch %%st(2)                 #22-22  ri0,r1,r0,ir0,i1,nt0r,i1\n"\
"faddp %%st,%%st(3)           #23-25  r1,r0,i0,i1,nt0r,i1   \n"  \
"fxch %%st(1)                 #23-23  r0,r1,i0,i1,nt0r,i1   \n"  \
"fsubp %%st,%%st(1)           #24-27  nt1r,i0,i1,nt0r,i1    \n"  \
"fxch %%st(1)                 #24-24  i0,nt1r,i1,nt0r,i1    \n"  \
"fadd %%st,%%st(4)            #25-28  i0,nt1r,i1,nt0r,nt0i  \n"  \
"fsubp %%st,%%st(2)           #26-29  nt1r,nt1i,nt0r,nt0i   \n"  \
"fxch %%st(2)                 #27     nt0r,nt1i,nt1r,nt0i   \n"  \
"fstpl %0                     #28     nt1i,nt1r,nt0i        \n"  \
"fxch %%st(1)                 #29     nt1r,nt1i,nt0i        \n"  \
"fstpl %2                     #30     nt1i,nt0i             \n"  \
"fxch %%st(1)                 #31     nt0i,nt1i             \n"  \
"fstpl %1                     #32     nt1i                  \n"  \
"fstpl %3                     #33                           \n"  \
: :"m"( _t0##r ),"m" ( _t0##i ), "m" (_t1##r ), "m" (_t1##i ),   \
"m"( _f0##r ), "m" ( _f0##i ),"m" ( _f1##r ),"m" ( _f1##i) :     \
"memory", "st");

#endif

#ifndef USE_ASM

#define cplx_mulmul(_t0,_t1,_f0,_f1)                                \
	{                                                           \
		y_limb_t _ar,_br;                                   \
		_ar = (_t1##r) * (_f1##r) - (_t1##i) * (_f1##i);    \
		_br = (_t0##r) * (_f0##r) - (_t0##i) * (_f0##i);    \
		_t1##i = (_t1##r) * (_f1##i) + (_t1##i) * (_f1##r); \
		_t0##i = (_t0##r) * (_f0##i) + (_t0##i) * (_f0##r); \
		_t1##r = _ar;                                       \
		_t0##r = _br;                                       \
	}

/* The factors here are supossed trig factors */
#define cplx_divdiv(_t0,_t1,_f0,_f1)                                \
	{                                                           \
		y_limb_t _ar,_br;                                   \
		_ar = (_t1##r) * (_f1##r) + (_t1##i) * (_f1##i);    \
		_br = (_t0##r) * (_f0##r) + (_t0##i) * (_f0##i);    \
		_t1##i = (_t1##i) * (_f1##r) - (_t1##r) * (_f1##i); \
		_t0##i = (_t0##i) * (_f0##r) - (_t0##r) * (_f0##i); \
		_t1##r = _ar;\
		_t0##r = _br;\
	}


#else

#define cplx_mulmul(_t0,_t1,_f0,_f1)                           \
__asm__ volatile("                                        \n"  \
"fldl %0                #1     t0r                        \n"  \
"fmull %4               #2-5   RR0                        \n"  \
"fldl %1                #3     t0i,RR0                    \n"  \
"fmull %5               #4-7   II0,RR0                    \n"  \
"fldl %0                #5     t0r,II0,rr0                \n"  \
"fmull %5               #6-9   RI0,II0,rr0                \n"  \
"fldl %1                #7     t0i,RI0,II0,rr0            \n"  \
"fmull %4               #8-10  IR0,RI0,ii0,rr0            \n"  \
"fxch %%st(3)           #8     rr0,RI0,ii0,IR0            \n"  \
"fsubp %%st,%%st(2)     #9-11  RI0,NR0,ir0                \n"  \
"fldl %2                #10    t1r,ri0,NR0,ir0            \n"  \
"fmull %6               #11    RR1,ri0,NR0,ir0            \n"  \
"fxch %%st(1)           #12    ri0,RR1,nr0,ir0            \n"  \
"faddp %%st,%%st(3)     #13-15 RR1,nr0,NI0                \n"  \
"fldl %3                #14    t1i,RR1,nr0,IR0            \n"  \
"fmull %7               #15-17 II1,rr1,nr0,ir0            \n"  \
"fldl %2                #16    t1r,II1,rr1,nr0,ir0        \n"  \
"fmull %7               #17-19 RI1,II1,rr1,nr0,ir0        \n"  \
"fldl %3                #18    t1i,RI1,ii1,rr1,nr0,ir0    \n"  \
"fmull %6               #19-21 IR1,RI1,ii1,rr1,nr0,ir0    \n"  \
"fxch %%st(3)           #19    rr1,RI1,ii1,IR1,nr0,ir0    \n"  \
"fsubp %%st,%%st(2)     #20    RI1,NR1,IR1,nr0,ir0        \n"  \
"fxch %%st(4)           #21    ir0,NR1,IR1,nr0,RI1        \n"  \
"fstpl %1               #22    nr1,ir1,nr0,ri1            \n"  \
"fxch %%st(1)           #22    ir1,nr1,nr0,ri1            \n"  \
"faddp %%st,%%st(3)     #23    nr1,nr0,NI1                \n"  \
"fstpl %2               #24    nr0,NI1                    \n"  \
"fstpl %0               #25    NI1                        \n"  \
"fstpl %3               #26                               \n"  \
: :"m"(_t0##r),"m"(_t0##i),"m"(_t1##r),"m"(_t1##i),            \
   "m"(_f0##r),"m"(_f0##i),"m"(_f1##r),"m"(_f1##i): "memory");

/* The factors here are supossed trig factors */
#define cplx_divdiv(_t0,_t1,_f0,_f1)                           \
__asm__ volatile("                                        \n"  \
"fldl %0                #1     t0r                        \n"  \
"fmull %4               #2-5   RR0                        \n"  \
"fldl %1                #3     t0i,RR0                    \n"  \
"fmull %5               #4-7   II0,RR0                    \n"  \
"fldl %0                #5     t0r,II0,rr0                \n"  \
"fmull %5               #6-9   RI0,II0,rr0                \n"  \
"fldl %1                #7     t0i,RI0,II0,rr0            \n"  \
"fmull %4               #8-10  IR0,RI0,ii0,rr0            \n"  \
"fxch %%st(3)           #8     rr0,RI0,ii0,IR0            \n"  \
"faddp %%st,%%st(2)     #9-11  RI0,NR0,ir0                \n"  \
"fldl %2                #10    t1r,ri0,NR0,ir0            \n"  \
"fmull %6               #11    RR1,ri0,NR0,ir0            \n"  \
"fxch %%st(3)           #11    ir0,rio,NR0,RR1            \n"  \
"fsubp %%st,%%st(1)     #12    NIO,NR0,RR1                \n"  \
"fxch %%st(2)           #12    RR1,nr0,NI0                \n"  \
"fldl %3                #14    t1i,RR1,nr0,IR0            \n"  \
"fmull %7               #15-17 II1,rr1,nr0,ir0            \n"  \
"fldl %2                #16    t1r,II1,rr1,nr0,ir0        \n"  \
"fmull %7               #17-19 RI1,II1,rr1,nr0,ir0        \n"  \
"fldl %3                #18    t1i,RI1,ii1,rr1,nr0,ir0    \n"  \
"fmull %6               #19-21 IR1,RI1,ii1,rr1,nr0,ir0    \n"  \
"fxch %%st(3)           #19    rr1,RI1,ii1,IR1,nr0,ir0    \n"  \
"faddp %%st,%%st(2)     #20    RI1,NR1,IR1,nr0,ir0        \n"  \
"fxch %%st(4)           #21    ir0,NR1,IR1,nr0,RI1        \n"  \
"fstpl %1               #22    nr1,ir1,nr0,ri1            \n"  \
"fxch %%st(1)           #22    ir1,nr1,nr0,ri1            \n"  \
"fsubp %%st,%%st(3)     #23    nr1,nr0,NI1                \n"  \
"fstpl %2               #24    nr0,NI1                    \n"  \
"fstpl %0               #25    NI1                        \n"  \
"fstpl %3               #26                               \n"  \
: :"m"(_t0##r),"m"(_t0##i),"m"(_t1##r),"m"(_t1##i),            \
   "m"(_f0##r),"m"(_f0##i),"m"(_f1##r),"m"(_f1##i): "memory");

#endif

#ifndef USE_ASM

# define cplx_load_mulmuladdsub_pp_no_fetch(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
	{                                                                     \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci;                             \
		_br= *( _pd0 + _index);                                       \
		_bi= *( _pd0 + _index + 1);                                   \
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i);                    \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r);                    \
		_br= *( _pd1 + _index);                                       \
		_bi= *( _pd1 + _index + 1);                                   \
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i);                    \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r);                    \
		_t0##r = _ar + _cr;                                           \
		_t1##r = _ar - _cr;                                           \
		_t0##i = _ai + _ci;                                           \
		_t1##i = _ai - _ci;                                           \
	}

#else

#define cplx_load_mulmuladdsub_pp_no_fetch(_t0,_t1,_f0,_f1,_index,_pd0,_pd1)  \
 __asm__ volatile("                                           \n"    \
"fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache        \n"    \
"fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r                \n"    \
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

/**********************************************************************
   cplx_divmul macro:
		_t0= _f0 / _f1;
		_t1= _f0 * _f1;
   where all the variables  are pseudo-complex. _f0 and _f1 are module
   1 pseudo-complex (like trig. factros)
**********************************************************************/
#ifndef USE_ASM

#define cplx_divmul(_t0,_t1,_f0,_f1)      \
	{                                 \
		y_limb_t _w0,_w1,_w2,_w3; \
		_w0 = _f0##r * _f1##r;    \
		_w1 = _f0##i * _f1##i;    \
		_w2 = _f0##i * _f1##r;    \
		_w3 = _f0##r * _f1##i;    \
		_t0##r = _w0 + _w1;       \
		_t1##r = _w0 - _w1;       \
		_t0##i = _w2 - _w3;       \
		_t1##i = _w2 + _w3;       \
	}

#else

#define cplx_divmul(_t0,_t1,_f0,_f1)                   \
	__asm__ volatile("                         \n" \
"fldl %4	       	 # 1     f0r               \n" \
"fmull %6                # 2-6   RR                \n" \
"fldl %5                 # 3     f0i,RR            \n" \
"fmull %7                # 4-8   II,RR             \n" \
"fldl %5                 # 5     f0i,II,RR         \n" \
"fmull %6                # 6-10  IR,II,RR          \n" \
"fldl %4                 # 7     f0r,IR,II,rr      \n" \
"fmull %7                # 8-12  RI,IR,II,rr       \n" \
"fld %%st(3)             # 9     rr,RI,IR,ii,rr    \n" \
"fsub %%st(3)            # 10-14 T1R,RI,IR,ii,rr   \n" \
"fxch %%st(3)            # 10    ii,RI,IR,T1R,rr   \n" \
"faddp %%st,%%st(4)      # 11-15 RI,ir,T1R,T0R     \n" \
"fld %%st(1)             # 12    ir,RI,ir,T1R,T0R  \n" \
"fsub %%st(1)            # 13-17 T0I,ri,ir,T1R,T0R \n" \
"fxch %%st(2)            # 13    ir,ri,T0I,T1R,T0R \n" \
"faddp %%st,%%st(1)      # 14-18 T1I,T0I,T1R,T0R   \n" \
"fxch %%st(2)            # 15    t1r,T0I,t1I,T0R   \n" \
"fstpl %2                # 16    T0I,t1I,t0r       \n" \
"fstpl %1                # 17    T1I,t0r           \n" \
"fstpl %3                # 18    t0r               \n" \
"fstpl %0                # 19                      \n" \
: :"m"(_t0##r),"m"(_t0##i),"m"(_t1##r),"m"(_t1##i),    \
"m"(_f0##r),"m"(_f0##i),"m"(_f1##r),"m"(_f1##i):"memory");

#endif


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/
#ifndef USE_ASM

# define cplx_addsub(_t0,_t1)            \
	{                                \
		y_limb_t _a=_t0##r;      \
		_t0##r += (_t1##r);      \
		_t1##r = _a - (_t1##r);  \
		_a = (_t0##i);           \
		_t0##i += (_t1##i);      \
		_t1##i = _a - (_t1##i);  \
		}
#else
# define cplx_addsub(_t0,_t1)                         \
__asm__ volatile( "                               \n" \
"fldl %0             #1-1 t0r                     \n" \
"faddl %2            #2-5 NT0R                    \n" \
"fldl %0             #3-3 t0r,NT0R                \n" \
"fsubl %2            #4-7 NT1R,NT0R               \n" \
"fldl %1             #5-5 t0i,NT1R,NT0R           \n" \
"faddl %3            #6-9 NT0I,NT1R,ntor          \n" \
"fxch %%st(2)        #6   nt0r,NT1R,NT0I          \n" \
"fldl %1             #7   t0i,nt0r,NT1R,NT0I      \n" \
"fsubl %3            #8-11 NT1I,nt0r,NT1R,NTOI    \n" \
"fxch %%st(3)        #9    NTOI,nt0r,nt1R,NT1I    \n" \
"fstpl %1            #10   nt0r,NT1R,NT1I         \n" \
"fstpl %0            #11   NT1R,NT1I              \n" \
"fstpl %2            #12   nt1I                   \n" \
"fstpl %3            #13   void                   \n" \
: :"m"( _t0##r ),"m" ( _t0##i ), "m" (_t1##r),        \
"m" (_t1##i): "memory", "st");


#endif


#ifndef USE_ASM

#define cplx_load_addsub_pp_no_fetch(_t0,_t1,_index,_pd0,_pd1) \
	{                                                      \
		y_limb_t _ar,_ai;                              \
		_ar = *( _pd1 + _index);                       \
		_ai = *( _pd1 + _index + 1);                   \
		_t0##r = *( _pd0 + _index) + _ar;              \
		_t1##r = *( _pd0 + _index) - _ar;              \
		_t0##i = *( _pd0 + _index + 1) + _ai;          \
		_t1##i = *( _pd0 + _index + 1) - _ai;          \
	}

#else

#define cplx_load_addsub_pp_no_fetch(_t0,_t1,_index,_pd0,_pd1)    \
	__asm__ volatile( "                                   \n" \
"fldl (%%esi,%%ebx,8)         # 4  t0r                        \n" \
"faddl (%%edi,%%ebx,8)        # 8-12   NT0R                   \n" \
"fldl  (%%esi,%%ebx,8)        # 9      t0r,NT0R               \n" \
"fsubl (%%edi,%%ebx,8)        # 10-14  NT1R,NT0R              \n" \
"fldl  8(%%esi,%%ebx,8)       # 11     t0i,NT1R,NT0R          \n" \
"faddl 8(%%edi,%%ebx,8)       # 12-16  NT0I,NT1R,ntor         \n" \
"fxch %%st(2)                 # 12     nt0r,NT1R,NT0I         \n" \
"fldl 8(%%esi,%%ebx,8)        # 13     t0i,nt0r,NT1R,NT0I     \n" \
"fsubl 8(%%edi,%%ebx,8)       # 14-18  NT1I,nt0r,NT1R,NTOI    \n" \
"fxch %%st(3)                 # 15     NTOI,nt0r,nt1R,NT1I    \n" \
"fstpl  %4                    # 18                            \n" \
"fstpl  %3                    # 19                            \n" \
"fstpl  %5                    # 20                            \n" \
"fstpl  %6                    # 21                            \n" \
::"b"(_index),"S"(_pd0),"D"(_pd1),"m"( _t0##r ),"m" ( _t0##i ),   \
"m" (_t1##r), "m" (_t1##i):"memory", "st");

#endif


/**********************************************************************
   cplx_mul macro:
		_t= _f0*f1;
   where all the variables are pseudo complex. 
**********************************************************************/
#define cplx_mul(_t,_f0,_f1)                                          \
	{                                                             \
		y_limb_t _r;                                          \
		_r= ((_f0##r) * (_f1##r)) - ((_f0##i) * (_f1##i));    \
		_t##i= ((_f0##r) * (_f1##i)) + ((_f0##i) * (_f1##r)); \
		_t##r = _r;                                           \
	}

/**********************************************************************
   cplx_squar macro:
		_t= _f0 * _f0;
   where all the variables are pseudo-complex. 
***********************************************************************/
#define cplx_squar(_t,_f0)                                            \
	{                                                             \
		y_limb_t _r;                                          \
		_r= (_f0##r +_f0##i ) * (_f0##r - _f0##i);            \
		_t##i= 2.0 * (_f0##r) * (_f0##i) ;                    \
		_t##r = _r;                                           \
	}

/***********************************************************************
   cplx_add macro:
		_t= _s0 + _s1;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_add(_t,_s0,_s1)        \
	_t##r= (_s0##r) + (_s1##r); \
	_t##i= (_s0##i) + (_s1##i);

/***********************************************************************
   cplx_sub macro:
		_t= _m - _s;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_sub(_t,_m,_s)            \
	_t##r= (_m##r) - (_s##r);     \
	_t##i= (_m##i) - (_s##i);


/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM

# define cplx_mul_1_4_F_addsub(_t0,_t1)                          \
	{                                                        \
		y_limb_t _ar= _t0##r, _ai = _t0##i ,_a = _t1##r; \
		_t0##r += (_t1##i);                              \
		_t0##i -= _a;                                    \
		_t1##r = _ar - (_t1##i);                         \
		_t1##i = _ai + _a;                               \
	}

#else

# define cplx_mul_1_4_F_addsub(_t0,_t1)              \
	__asm__ volatile("                      \n"  \
"fldl %0            #1-1 t0r                    \n"  \
"faddl %3           #2-5 nt0r                   \n"  \
"fldl %1            #3-3 t0i,nt0r               \n"  \
"fsubl %2           #4-7 nt0i,nt0r              \n"  \
"fldl %0            #5-5 t0r,nt0i,nt0r          \n"  \
"fsubl %3           #6-9 nt1r,nt0i,nt0r         \n"  \
"fldl %1            #7-7 t0i,nt1r,nt0i,nt0r     \n"  \
"faddl %2           #8-11 nt1i,nt1r,nt0i,nt0r   \n"  \
"fxch %%st(3)       #8-8  nt0r,nt1r,nt0i,nt1i   \n"  \
"fstpl %0           #9-9  nt1r,nt0i,nt1i        \n"  \
"fxch %%st(1)       #9-9  nt0i,nt1r,nt1i        \n"  \
"fstpl %1           #10-10  nt1r,nt1i           \n"  \
"fstpl %2           #11-11  nt1i                \n"  \
"fstpl %3           #12-12  void                \n"  \
  : :"m"( _t0##r ),"m" ( _t0##i ), "m" ( _t1##r ),   \
  "m" (_t1##i ): "memory");

#endif

/***********************************************************************
   cplx_mul_1_4_F_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
	_d[addr( _k0<<1 )]= _t0##r + _t1##i;            \
	_d[addr( _k1<<1 )]= _t0##r - _t1##i;            \
	_d[addr( _k0<<1 )+ 1]= _t0##i - _t1##r;         \
	_d[addr( _k1<<1 )+ 1]= _t0##i + _t1##r;

#ifndef USE_ASM

#define cplx_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	*(_pd0 + _index) = _t0##r + _t1##i;                     \
	*(_pd1 + _index) = _t0##r - _t1##i;                     \
	*( _pd0 + _index + 1) = _t0##i - _t1##r;                \
	*( _pd1 + _index + 1) = _t0##i + _t1##r;

#else

#define cplx_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	__asm__ volatile("                                 \n"  \
"fldl %0             	#1-1 t0r                           \n"  \
"faddl %3            	#2-5 nt0r                          \n"  \
"fldl %1             	#3-3 t0i,nt0r                      \n"  \
"fsubl %2            	#4-7 nt0i,nt0r                     \n"  \
"fldl %0             	#5-5 t0r,nt0i,nt0r                 \n"  \
"fsubl %3            	#6-9 nt1r,nt0i,nt0r                \n"  \
"fldl %1             	#7-7 t0i,nt1r,nt0i,nt0r            \n"  \
"faddl %2            	#8-11 nt1i,nt1r,nt0i,nt0r          \n"  \
"fxch %%st(3)        	#8-8  nt0r,nt1r,nt0i,nt1i          \n"  \
"fstpl (%%esi,%%ebx,8)  #9-9  nt1r,nt0i,nt1i               \n"  \
"fxch %%st(1)        	#9-9  nt0i,nt1r,nt1i               \n"  \
"fstpl 8(%%esi,%%ebx,8) #10-10  nt1r,nt1i                  \n"  \
"fstpl (%%edi,%%ebx,8)  #11-11  nt1i                       \n"  \
"fstpl 8(%%edi,%%ebx,8) #12-12  void                       \n"  \
 : :"m"( _t0##r ),"m" ( _t0##i ), "m" ( _t1##r ), "m" (_t1##i ),\
 "b"(_index),"S"(_pd0),"D"(_pd1):"memory");

#endif


/***********************************************************************
   cplx_mul_1_4_B macro:
		_t = _t1 * G^(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM

#define cplx_mul_1_4_B_addsub(_t0,_t1)                       \
	{                                                    \
		y_limb_t _ar=_t0##r, _ai=_t0##i ,_a= _t1##r; \
		_t0##r -= (_t1##i);                          \
		_t0##i += _a;                                \
		_t1##r = _ar + (_t1##i);                     \
		_t1##i = _ai - _a;                           \
	}

#else
#define cplx_mul_1_4_B_addsub(_t0,_t1)                       \
	__asm__ volatile("                              \n"  \
"fldl %0             #1-1 t0r                           \n"  \
"fsubl %3            #2-5 nt0r                          \n"  \
"fldl %1             #3-3 t0i,nt0r                      \n"  \
"faddl %2            #4-7 nt0i,nt0r                     \n"  \
"fldl %0             #5-5 t0r,nt0i,nt0r                 \n"  \
"faddl %3            #6-9 nt1r,nt0i,nt0r                \n"  \
"fldl %1             #7-7 t0i,nt1r,nt0i,nt0r            \n"  \
"fsubl %2            #8-11 nt1i,nt1r,nt0i,nt0r          \n"  \
"fxch %%st(3)        #8-8  nt0r,nt1r,nt0i,nt1i          \n"  \
"fstpl %0            #9-9  nt1r,nt0i,nt1i               \n"  \
"fxch %%st(1)        #9-9  nt0i,nt1r,nt1i               \n"  \
"fstpl %1            #10-10  nt1r,nt1i                  \n"  \
"fstpl %2            #11-11  nt1i                       \n"  \
"fstpl %3            #12-12  void                       \n"  \
 : :"m"( _t0##r ),"m"( _t0##i ),"m"( _t1##r ), "m" (_t1##i ):\
	    "memory");

#endif


/***********************************************************************
   cplx_mul_1_4_B_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
	_d[addr( _k0<<1 )]= _t0##r - _t1##i;            \
	_d[addr( _k1<<1 )]= _t0##r + _t1##i;            \
	_d[addr( _k0<<1 )+ 1]= _t0##i + _t1##r;         \
	_d[addr( _k1<<1 )+ 1]= _t0##i - _t1##r;

#ifndef USE_ASM

#define cplx_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	*( _pd0 + _index) = _t0##r - _t1##i;                    \
	*( _pd1 + _index) = _t0##r + _t1##i;                    \
	*( _pd0 + _index + 1) = _t0##i + _t1##r;                \
	*( _pd1 + _index + 1) = _t0##i - _t1##r;

#else
#define cplx_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	__asm__ volatile("                                 \n"  \
"fldl %0              	#1-1 t0r                           \n"  \
"fsubl %3             	#2-5 nt0r                          \n"  \
"fldl %1              	#3-3 t0i,nt0r                      \n"  \
"faddl %2             	#4-7 nt0i,nt0r                     \n"  \
"fldl %0              	#5-5 t0r,nt0i,nt0r                 \n"  \
"faddl %3             	#6-9 nt1r,nt0i,nt0r                \n"  \
"fldl %1              	#7-7 t0i,nt1r,nt0i,nt0r            \n"  \
"fsubl %2             	#8-11 nt1i,nt1r,nt0i,nt0r          \n"  \
"fxch %%st(3)         	#8-8  nt0r,nt1r,nt0i,nt1i          \n"  \
"fstpl (%%esi,%%ebx,8) 	#9-9  nt1r,nt0i,nt1i               \n"  \
"fxch %%st(1)         	#9-9  nt0i,nt1r,nt1i               \n"  \
"fstpl 8(%%esi,%%ebx,8)	#10-10  nt1r,nt1i                  \n"  \
"fstpl (%%edi,%%ebx,8)	#11-11  nt1i                       \n"  \
"fstpl 8(%%edi,%%ebx,8) #12-12  void                       \n"  \
 : :"m"( _t0##r ),"m" ( _t0##i ), "m" ( _t1##r ), "m" (_t1##i ),\
   "b"(_index),"S"(_pd0),"D"(_pd1): "memory");

#endif


/***********************************************************************
   cplx_mul_1_8_F macro:
		_t = _t1 * G^-(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM
# define cplx_mul_1_8_F_addsub(_tt0,_tt1)          \
	{                                          \
		y_limb_t _ar,_ai;                  \
		_ar = (_tt1##r + _tt1##i )*F_1_8r; \
		_ai = (_tt1##i - _tt1##r )*F_1_8r; \
		_tt1##r = _tt0##r - _ar;           \
		_tt0##r += _ar;                    \
		_tt1##i = _tt0##i - _ai;           \
		_tt0##i += _ai;                    \
	}
#else
# define cplx_mul_1_8_F_addsub(_tt0,_tt1)                    \
        __asm__ volatile("                              \n"  \
"fldl %3                 #1-1 t1i                       \n"  \
"fld  %%st               #2-2 t1i,t1i                   \n"  \
"faddl %2                #3-6 t1i+t1r,t1i               \n"  \
"fxch %%st(1)            #3-3 t1i,t1i+t1r               \n"  \
"fsubl %2                #4-7 t1i-t1r,t1i+t1r           \n"  \
"fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r       \n"  \
"fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #6-10 nt1r,t1i-t1r,t0r         \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r      \n"  \
"fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r      \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #8-12 nt1i,nt1r,t0i,t0r        \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i         \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i     \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i    \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i\n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i\n"\
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i\n"\
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i\n"   \
"fstpl %0                #15-15 nnt1r,nnt0i,nnt1i       \n"  \
"fstpl %2                #16-16 nnt0i,nnt1i             \n"  \
"fstpl %1                #17-17 nnt1i                   \n"  \
"fstpl %3                #18-18 void                    \n"  \
: :"m"(_tt0##r ),"m" (_tt0##i ),"m" (_tt1##r),"m" (_tt1##i ) \
:   "memory");

#endif


/***********************************************************************
   cplx_mul_1_8_F_addsub_store macro:
		_t = _t1 * G^-(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_F_addsub_store(_d,_k0,_k1,_t0,_t1)    \
	{                                                  \
		y_limb_t _ar,_ai;                          \
		_ar = (_t1##r + _t1##i )*F_1_8r;           \
		_ai = (_t1##i - _t1##r )*F_1_8r;           \
		_d[addr( _k0<<1 )]= _t0##r + _ar;          \
		_d[addr( _k1<<1 )]= _t0##r - _ar;          \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai;       \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai;       \
	}


#ifndef USE_ASM

#define cplx_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{                                                       \
		y_limb_t _ar,_ai;                               \
		_ar = (_t1##r + _t1##i )*F_1_8r;                \
		_ai = (_t1##i - _t1##r )*F_1_8r;                \
		*( _pd0 + _index)= _t0##r + _ar;                \
		*( _pd1 + _index)= _t0##r - _ar;                \
		*( _pd0 + _index + 1)= _t0##i + _ai;            \
		*( _pd1 + _index + 1)= _t0##i - _ai;            \
	}
#else

#define cplx_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1)  \
        __asm__ volatile("                                  \n"  \
"fldl %3                 #1-1 t1i                           \n"  \
"fld  %%st               #2-2 t1i,t1i                       \n"  \
"faddl %2                #3-6 t1i+t1r,t1i                   \n"  \
"fxch %%st(1)            #3-3 t1i,t1i+t1r                   \n"  \
"fsubl %2                #4-7 t1i-t1r,t1i+t1r               \n"  \
"fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r           \n"  \
"fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r           \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #6-10 nt1r,t1i-t1r,t0r             \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r          \n"  \
"fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #8-12 nt1i,nt1r,t0i,t0r            \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i             \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i         \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i     \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i     \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i        \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i    \n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i  \n"  \
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i  \n"  \
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i     \n"  \
"fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i           \n"  \
"fstpl (%%edi,%%ebx,8)   #16-16 nnt0i,nnt1i                 \n"  \
"fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i                       \n"  \
"fstpl 8(%%edi,%%ebx,8)  #18-18 void                        \n"  \
  : :"m" (_t0##r ),"m" (_t0##i ),"m" (_t1##r),"m" (_t1##i ),     \
  "b"(_index),"S"(_pd0),"D"(_pd1): "memory");

#endif


/***********************************************************************
   cplx_mul_1_8_B macro:
		_t = _t1 * G^(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM

# define cplx_mul_1_8_B_addsub(_tt0,_tt1)          \
        {                                          \
		y_limb_t _ar,_ai;                  \
		_ar = (_tt1##r - _tt1##i )*F_1_8r; \
		_ai = (_tt1##i + _tt1##r )*F_1_8r; \
		_tt1##r = _tt0##r - _ar;           \
		_tt0##r += _ar;                    \
		_tt1##i = _tt0##i - _ai;           \
		_tt0##i += _ai;                    \
	}

#else

#define cplx_mul_1_8_B_addsub(_tt0,_tt1)                     \
        __asm__ volatile("                              \n"  \
"fldl %2                 #1-1 t1r                       \n"  \
"fld  %%st               #2-2 t1r,t1r                   \n"  \
"fsubl %3                #3-6 t1r-t1i,t1i               \n"  \
"fxch %%st(1)            #3-3 t1i,t1r-t1i               \n"  \
"faddl %3                #4-7 t1i+t1r,t1i-t1r           \n"  \
"fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r       \n"  \
"fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #6-10 nt1r,t1i+t1r,t0r         \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r      \n"  \
"fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r      \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #8-12 nt1i,nt1r,t0i,t0r        \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i         \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i     \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i    \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i\n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i\n"\
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i\n"\
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i \n"  \
"fstpl %0                #15-15 nnt1r,nnt0i,nnt1i       \n"  \
"fstpl %2                #16-16 nnt0i,nnt1i             \n"  \
"fstpl %1                #17-17 nnt1i                   \n"  \
"fstpl %3                #18-18 void                    \n"  \
: :"m" (_tt0##r ),"m"(_tt0##i ),"m" (_tt1##r),"m" (_tt1##i ) \
:   "memory");

#endif

/***********************************************************************
   cplx_mul_1_8_B_addsub_store macro:
		_t = _t1 * G^(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{                                               \
		y_limb_t _ar,_ai;                       \
		_ar = (_t1##r - _t1##i )*F_1_8r;        \
		_ai = (_t1##i + _t1##r )*F_1_8r;        \
		_d[addr( _k0<<1 )]= _t0##r + _ar;       \
		_d[addr( _k1<<1 )]= _t0##r - _ar;       \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai;    \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai;    \
	}

#ifndef USE_ASM

#define cplx_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{                                                       \
		y_limb_t _ar,_ai;                               \
		_ar = (_t1##r - _t1##i )*F_1_8r;                \
		_ai = (_t1##i + _t1##r )*F_1_8r;                \
		*( _pd0 + _index) = _t0##r + _ar;               \
		*( _pd1 + _index) = _t0##r - _ar;               \
		*( _pd0 + _index + 1) = _t0##i + _ai;           \
		*( _pd1 + _index + 1) = _t0##i - _ai;           \
	}
#else

#define cplx_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
        __asm__ volatile("                                 \n"  \
"fldl %2                 #1-1 t1r                          \n"  \
"fld  %%st               #2-2 t1r,t1r                      \n"  \
"fsubl %3                #3-6 t1r-t1i,t1i                  \n"  \
"fxch %%st(1)            #3-3 t1i,t1r-t1i                  \n"  \
"faddl %3                #4-7 t1i+t1r,t1i-t1r              \n"  \
"fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r          \n"  \
"fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #6-10 nt1r,t1i+t1r,t0r            \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r         \n"  \
"fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r         \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8r) "            #8-12 nt1i,nt1r,t0i,t0r           \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i            \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i        \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i    \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i    \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i       \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i   \n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i \n"  \
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i \n"  \
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i    \n"  \
"fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i          \n"  \
"fstpl (%%edi,%%ebx,8)   #16-16 nnt0i,nnt1i                \n"  \
"fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i                      \n"  \
"fstpl 8(%%edi,%%ebx,8)  #18-18 void                       \n"  \
      : :"m" (_t0##r ),"m" (_t0##i ),"m" (_t1##r),"m" (_t1##i ),\
	"b"(_index),"S"(_pd0),"D"(_pd1):  "memory");
#endif


/***********************************************************************
   cplx_mul_3_8_F macro:
		_t = _t1 * G^-(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/

#ifndef USE_ASM

# define cplx_mul_3_8_F_addsub(_tt0,_tt1)          \
	{                                          \
		y_limb_t _ar,_ai;                  \
		_ar = (_tt1##r - _tt1##i )*F_1_8i; \
		_ai = (_tt1##i + _tt1##r )*F_1_8i; \
		_tt1##r = _tt0##r - _ar;           \
		_tt0##r += _ar;                    \
		_tt1##i = _tt0##i - _ai;           \
		_tt0##i += _ai;                    \
	}

#else

#define cplx_mul_3_8_F_addsub(_tt0,_tt1)                     \
	__asm__ volatile("                              \n"  \
"fldl %2                 #1-1 t1r                       \n"  \
"fld  %%st               #2-2 t1r,t1r                   \n"  \
"fsubl %3                #3-6 t1r-t1i,t1i               \n"  \
"fxch %%st(1)            #3-3 t1i,t1r-t1i               \n"  \
"faddl %3                #4-7 t1i+t1r,t1i-t1r           \n"  \
"fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r       \n"  \
"fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #6-10 nt1r,t1i+t1r,t0r         \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r      \n"  \
"fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r      \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #8-12 nt1i,nt1r,t0i,t0r        \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i         \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i     \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i    \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i\n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i\n"\
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i\n"\
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i \n"  \
"fstpl %0                #15-15 nnt1r,nnt0i,nnt1i       \n"  \
"fstpl %2                #16-16 nnt0i,nnt1i             \n"  \
"fstpl %1                #17-17 nnt1i                   \n"  \
"fstpl %3                #18-18 void                    \n"  \
: :"m"(_tt0##r ),"m" (_tt0##i ),"m" (_tt1##r),"m" (_tt1##i ) \
:   "memory");

#endif

/***********************************************************************
   cplx_mul_3_8_F_addsub_store macro:
		_t = _t1 * G^-(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{                                               \
		y_limb_t _ar,_ai;                       \
		_ar = (_t1##r - _t1##i )*F_1_8i;        \
		_ai = (_t1##i + _t1##r )*F_1_8i;        \
		_d[addr( _k0<<1 )]= _t0##r + _ar;       \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai;    \
		_d[addr( _k1<<1 )]= _t0##r - _ar;       \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai;    \
	}

#ifndef USE_ASM

#define cplx_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{                                                       \
		y_limb_t _ar,_ai;                               \
		_ar = (_t1##r - _t1##i )*F_1_8i;                \
		_ai = (_t1##i + _t1##r )*F_1_8i;                \
		*( _pd0 + _index) = _t0##r + _ar;               \
		*( _pd0 + _index + 1) = _t0##i + _ai;           \
		*( _pd1 + _index) = _t0##r - _ar;               \
		*( _pd1 + _index + 1) = _t0##i - _ai;           \
	}

#else

#define cplx_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	__asm__ volatile("                                 \n"  \
"fldl %2                 #1-1 t1r                          \n"  \
"fld  %%st               #2-2 t1r,t1r                      \n"  \
"fsubl %3                #3-6 t1r-t1i,t1i                  \n"  \
"fxch %%st(1)            #3-3 t1i,t1r-t1i                  \n"  \
"faddl %3                #4-7 t1i+t1r,t1i-t1r              \n"  \
"fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r          \n"  \
"fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #6-10 nt1r,t1i+t1r,t0r            \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r         \n"  \
"fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r         \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #8-12 nt1i,nt1r,t0i,t0r           \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i            \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i        \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i    \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i    \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i       \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i   \n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i \n"  \
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i \n"  \
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i    \n"  \
"fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i          \n"  \
"fstpl (%%edi,%%ebx,8)   #16-16 nnt0i,nnt1i                \n"  \
"fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i                      \n"  \
"fstpl 8(%%edi,%%ebx,8)  #18-18 void                       \n"  \
      : :"m" (_t0##r ),"m" (_t0##i ),"m" (_t1##r),"m" (_t1##i ),\
	"b"(_index),"S"(_pd0),"D"(_pd1):   "memory");
#endif

/***********************************************************************
   cplx_mul_3_8_B macro:
		_t = _t1 * G^(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM

#define cplx_mul_3_8_B_addsub(_tt0,_tt1)           \
        {                                          \
		y_limb_t _ar,_ai;                  \
		_ar = (_tt1##r + _tt1##i )*F_1_8i; \
		_ai = (_tt1##i - _tt1##r )*F_1_8i; \
		_tt1##r = _tt0##r - _ar;           \
		_tt0##r += _ar;                    \
		_tt1##i = _tt0##i - _ai;           \
		_tt0##i += _ai;                    \
	}

#else

#define cplx_mul_3_8_B_addsub(_tt0,_tt1)                     \
	__asm__ volatile("                              \n"  \
"fldl %3                 #1-1 t1i                       \n"  \
"fld  %%st               #2-2 t1i,t1i                   \n"  \
"faddl %2                #3-6 t1i+t1r,t1i               \n"  \
"fxch %%st(1)            #3-3 t1i,t1i+t1r               \n"  \
"fsubl %2                #4-7 t1i-t1r,t1i+t1r           \n"  \
"fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r       \n"  \
"fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #6-10 nt1r,t1i-t1r,t0r         \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r      \n"  \
"fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r      \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #8-12 nt1i,nt1r,t0i,t0r        \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i         \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i     \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i    \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i\n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i\n"\
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i\n"\
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i \n"  \
"fstpl %0                #15-15 nnt1r,nnt0i,nnt1i       \n"  \
"fstpl %2                #16-16 nnt0i,nnt1i             \n"  \
"fstpl %1                #17-17 nnt1i                   \n"  \
"fstpl %3                #18-18 void                    \n"  \
: :"m"(_tt0##r ),"m" (_tt0##i ),"m" (_tt1##r),"m" (_tt1##i ) \
:   "memory");

#endif

/***********************************************************************
   cplx_mul_3_8_B_addsub_store macro:
		_t = _t1 * G^(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_B_addsub_store(_d,_k0,_k1,_t0,_t1)      \
	{                                                    \
		y_limb_t _ar,_ai;                            \
		_ar = (_t1##r + _t1##i )*F_1_8i;             \
		_ai = (_t1##i - _t1##r )*F_1_8i;             \
		_d[addr( _k0<<1 )]= _t0##r + _ar;            \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai;         \
		_d[addr( _k1<<1 )]= _t0##r - _ar;            \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai;         \
	}

#ifndef USE_ASM

#define cplx_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{                                                       \
		y_limb_t _ar,_ai;                               \
		_ar = (_t1##r + _t1##i )*F_1_8i;                \
		_ai = (_t1##i - _t1##r )*F_1_8i;                \
		*( _pd0 + _index) = _t0##r + _ar;               \
		*( _pd0 + _index + 1) = _t0##i + _ai;           \
		*( _pd1 + _index) = _t0##r - _ar;               \
		*( _pd1 + _index + 1) = _t0##i - _ai;           \
	}

#else

#define cplx_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	__asm__ volatile("                                \n"  \
"fldl %3                 #1-1 t1i                          \n"  \
"fld  %%st               #2-2 t1i,t1i                      \n"  \
"faddl %2                #3-6 t1i+t1r,t1i                  \n"  \
"fxch %%st(1)            #3-3 t1i,t1i+t1r                  \n"  \
"fsubl %2                #4-7 t1i-t1r,t1i+t1r              \n"  \
"fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r          \n"  \
"fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #6-10 nt1r,t1i-t1r,t0r            \n"  \
"fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r         \n"  \
"fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r         \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_8i) "            #8-12 nt1i,nt1r,t0i,t0r           \n"  \
"fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i            \n"  \
"fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i        \n"  \
"fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i    \n"  \
"fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i    \n"  \
"fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i       \n"  \
"fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i   \n"  \
"fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i \n"  \
"fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i \n"  \
"fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i    \n"  \
"fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i          \n"  \
"fstpl (%%edi,%%ebx,8)   #16-16 nnt0i,nnt1i                \n"  \
"fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i                      \n"  \
"fstpl 8(%%edi,%%ebx,8)  #18-18 void                       \n"  \
      : :"m" (_t0##r ),"m" (_t0##i ),"m" (_t1##r),"m" (_t1##i ),\
	"b"(_index),"S"(_pd0),"D"(_pd1): "memory");
#endif


/***********************************************************************
 This is the dyadic mul for nested complex representation 
 xk =(xk + xmk* ) * (yk + ymk* ) + 2*((xk * yk) - (xmk* * ymk* ) -
     G^(-k) * (xk - xmk* ) * (yk - ymk* );
 xkm =(xkm + xk* ) * (ymk + yk* ) + 2*((xmk * xk*) - (ymk * yk* ) -
     G^(k) * (xmk - xk* ) * (ymk - yk* );
 where all are pseudo-complex vars. 
***********************************************************************/
/* 17 Fadds 18 FMULS */
#define conv_nested(xk,xmk,yk,ymk,gkr,gki)          \
{                                                   \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;                      \
   _y0i=(xk##i + xmk##i)*0.25;                      \
   _y1r=(yk##r - ymk##r);                           \
   _y1i=(yk##i + ymk##i);                           \
   _y2r=_y0r*_y1r - _y0i*_y1i;                      \
   _y2i=_y0r*_y1i + _y0i*_y1r;                      \
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);          \
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);          \
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);      \
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);      \
   _y3r= ((1.0 +(gkr)) * _y2r - (gki) * _y2i);      \
   _y3i= ((1.0 +(gkr)) * _y2i + (gki) * _y2r);      \
   xk##r = _y0r - _y3r;                             \
   xk##i = _y0i - _y3i;                             \
   xmk##r = _y1r - _y3r;                            \
   xmk##i = _y1i + _y3i;                            \
}

/*********************************************************************
This is the dyadic square for nested-compls representation
 xk =(xk + xmk* ) * (xk + xmk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(-k) * (xk - xmk* ) * (xk - xmk* );
 xkm =(xkm + xk* ) * (xmk + xk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(k) * (xmk - xk* ) * (xmk - xk* );
 where all are pseudo-complex vars. 
**********************************************************************/
/* 18 fadds  12  Fmuls*/

#ifndef USE_ASM
#define square_nested(xk,xmk,gkr,gki)            \
{                                                \
   y_limb_t _r0,_r1,_r2,_r3;                     \
   _r0=(xk##r - xmk##r)*0.5;                     \
   _r1=(xk##i + xmk##i)*0.5;                     \
   _r2=(xk##r + xk##i) * (xk##r - xk##i);        \
   xk##i *= xk##r;                               \
   xk##r = _r2;                                  \
   _r2 =(_r0 + _r1) * (_r0 - _r1);               \
   _r1 *= _r0;                                   \
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);   \
   xmk##i *= xmk##r;                             \
   xmk##r = _r3;                                 \
   _r1 +=  _r1;                                  \
   xk##i += xk##i;                               \
   xmk##i += xmk##i;                             \
   _r0= (((gkr)+ 1.0) * _r2 - (gki) * _r1);      \
   _r1= (((gkr)+ 1.0) * _r1 + (gki) * _r2);      \
   xk##r -= _r0;                                 \
   xmk##r -= _r0;                                \
   xk##i -= _r1;                                 \
   xmk##i += _r1;                                \
}

#else

#define square_nested(xk,xmk,gkr,gki)                \
__asm__ volatile("                              \n"  \
"fldl %0 	      #1   xkr                  \n"  \
"fsubl %2             #2-6 XR0                  \n"  \
"fldl %1              #3   xki,XR0              \n"  \
"faddl %3             #4-8 XR1,XR0              \n"  \
"fldl %0              #5   xkr,XR1,XR0          \n"  \
"faddl %1             #6-10 XK+,XR1,XR0         \n"  \
"fxch %%st(2)         #6    XR0,XR1,XK+         \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #7-11 R0,XR1,XK+          \n"  \
"fldl %0              #8    xkr,R0,XR1,XK+      \n"  \
"fsubl %1             #9-12 XK-,R0,XR1,XK+      \n"  \
"fxch %%st(2)         #9    XR1,R0,XK-,XK+      \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #10-14 R1,R0,XK-,XK+      \n"  \
"fldl %0              #11    xkr,R1,R0,XK-,kx+  \n"  \
"fmull %1             #12-16 RI,R1,r0,kx-,xk+   \n"  \
"fxch %%st(4)         #12    xk+,R1,r0,xk-,RI   \n"  \
"fmulp %%st,%%st(3)   #14-18 r1,r0,RR,RI        \n"  \
"fld %%st(1)          #15    r0,r1,r0,RR,RI     \n"  \
"fsub %%st(1)         #16-20 r-,r1,r0,RR,ri     \n"  \
"fld %%st(2)          #17    r0,R-,r1,r0,RR,ri  \n"  \
"fadd %%st(2)         #18-22 R+,R-,r1,r0,RR,ri  \n"  \
"fxch %%st(3)         #18    r0,R-,r1,R+,RR,ri  \n"  \
"fmulp %%st,%%st(2)   #19-23 R-,R2,R+,rr,ri     \n"  \
"fxch %%st(4)         #19    ri,R2,R+,rr,R-     \n"  \
"fadd %%st            #20-24 RI,R2,R+,rr,r-     \n"  \
"fxch %%st(3)         #21    rr,R2,R+,RI,r-     \n"  \
"fstpl %0             #22    R2,r+,RI,r-        \n"  \
"fxch %%st(1)         #22    r+,R2,RI,r-        \n"  \
"fmulp %%st(3)        #23-27 r2,RI,R1           \n"  \
"fxch %%st(1)         #24    ri,r2,R1           \n"  \
"fstpl %1             #25    r2,R1              \n"  \
"fadd %%st            #26-30 R2,R1              \n"  \
"fldl %2              #27    xmr,R2,R1          \n"  \
"faddl %3             #28-32 XM+,R2,r1          \n"  \
"fldl %2              #29    xmi,XM+,R2,r1      \n"  \
"fsubl %3             #30    XM-,XM+,R2,r1      \n"  \
"fldl %2              #31    xmr,XM-,XM+,r2,r1  \n"  \
"fmull %3             #32-36 RI,XM-,xm+,r2,r1   \n"  \
"fxch %%st(2)         #32    xm+,XM-,RI,r2,r1   \n"  \
"fmulp %%st,%%st(1)   #34-38 RR,RI,r2,r1        \n"  \
"fld1                 #35    1.0,RR,RI,r2,r1 ****\n" \
"faddl %4             #36-40 GR,RR,RI,r2,r1  ****\n" \
"fxch %%st(2)         #37    RI,RR,GR,r2,r1     \n"  \
"fadd %%st            #38-42 RI,rr,GR,r2,r1     \n"  \
"fxch %%st(1)         #39    rr,RI,GR,r2,r1     \n"  \
"fstpl %2             #40    RI,GR,r2,r1        \n"  \
"fld %%st(3)          #41    r1,RI,gr,r2,r1     \n"  \
"fmul %%st(2)         #42-46 RR,RI,gr,r2,r1     \n"  \
"fxch %%st(1)         #43    RI,RR,gr,r2,r1     \n"  \
"fstpl %3             #44    RR,gr,r2,r1        \n"  \
"fld %%st(2)          #45    r2,RR,gr,r2,r1     \n"  \
"fmull %5             #46    II,RR,gr,r2,r1   ****\n"\
"fxch %%st(2)         #46    gr,RR,II,r2,r1     \n"  \
"fmulp %%st(3)        #47-51 rr,II,RI,r1        \n"  \
"fsubp %%st,%%st(1)   #48-52 R0,RI,r1    **     \n"  \
"fxch %%st(2)         #48    r1,RI,R0           \n"  \
"fmull %5             #49-53 IR,RI,R0    **     \n"  \
"fldl %0              #50    xr,IR,RI,R0        \n"  \
"fsub %%st(3)         #51-55 XR,IR,RI,r0        \n"  \
"fldl %2              #52    xmr,XR,IR,ri,r0    \n"  \
"fsubp %%st,%%st(4)   #53-57 XR,ir,ri,XMR       \n"  \
"fxch %%st(1)         #53    ir,XR,ri,XMR       \n"  \
"faddp %%st(2)        #54-58 XR,R1,XMR          \n"  \
"fstpl %0             #55    R1,XMR             \n"  \
"fldl %1              #56    xi,R1,XMR          \n"  \
"fsub %%st(1)         #58-62 XI,R1,xmr          \n"  \
"fxch %%st(1)         #58    r1,XI,xmr          \n"  \
"faddl %3             #59-63 XMI,XI,xmr         \n"  \
"fxch %%st(2)         #60    xmr,XI,XMI         \n"  \
"fstpl %2             #61    xi,XMI             \n"  \
"fstpl %1             #62    XMI                \n"  \
"fstpl %3             #63                       \n"  \
 : :"m" (xk##r),"m"(xk##i),"m"(xmk##r),"m"(xmk##i),  \
"m"(gkr), "m"(gki): "memory");

#endif


#define square_nested_eq(xk,gkr)                \
{                                               \
   y_limb_t _y0r=xk##r,_y0i=xk##i;              \
   xk##r =(_y0r * _y0r + (gkr) * _y0i * _y0i);  \
   xk##i = 2.0 * _y0r * _y0i;                   \
}

/* These are version to exploid better the symmetries */
#define conv_nested_1_4(xk,xmk,yk,ymk,gkr,gki)         \
{                                                      \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;   \
   _y0r=(xk##r - xmk##r)*0.25;                         \
   _y0i=(xk##i + xmk##i)*0.25;                         \
   _y1r=(yk##r - ymk##r);                              \
   _y1i=(yk##i + ymk##i);                              \
   _y2r=_y0r*_y1r - _y0i*_y1i;                         \
   _y2i=_y0r*_y1i + _y0i*_y1r;                         \
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);             \
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);             \
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);         \
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);         \
   _y3r= ((gki +1.0) * _y2r + gkr * _y2i);             \
   _y3i= ((gki +1.0) * _y2i - gkr * _y2r);             \
   xk##r = _y0r - _y3r;                                \
   xk##i = _y0i - _y3i;                                \
   xmk##r = _y1r - _y3r;                               \
   xmk##i = _y1i + _y3i;                               \
}

#define conv_nested_1_2(xk,xmk,yk,ymk,gkr,gki)         \
{                                                      \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;   \
   _y0r=(xk##r - xmk##r)*0.25;                         \
   _y0i=(xk##i + xmk##i)*0.25;                         \
   _y1r=(yk##r - ymk##r);                              \
   _y1i=(yk##i + ymk##i);                              \
   _y2r=_y0r*_y1r - _y0i*_y1i;                         \
   _y2i=_y0r*_y1i + _y0i*_y1r;                         \
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);             \
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);             \
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);         \
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);         \
   _y3r= ((gkr - 1.0) * _y2r - gki * _y2i);            \
   _y3i= ((gkr - 1.0) * _y2i + gki * _y2r);            \
   xk##r = _y0r + _y3r;                                \
   xk##i = _y0i + _y3i;                                \
   xmk##r = _y1r + _y3r;                               \
   xmk##i = _y1i - _y3i;                               \
}

#define conv_nested_3_4(xk,xmk,yk,ymk,gkr,gki)         \
{                                                      \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;   \
   _y0r=(xk##r - xmk##r)*0.25;                         \
   _y0i=(xk##i + xmk##i)*0.25;                         \
   _y1r=(yk##r - ymk##r);                              \
   _y1i=(yk##i + ymk##i);                              \
   _y2r=_y0r*_y1r - _y0i*_y1i;                         \
   _y2i=_y0r*_y1i + _y0i*_y1r;                         \
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);             \
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);             \
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);         \
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);         \
   _y3r= ((gki - 1.0) * _y2r + gkr * _y2i);            \
   _y3i= ((gki - 1.0) * _y2i - gkr * _y2r);            \
   xk##r = _y0r + _y3r;                                \
   xk##i = _y0i + _y3i;                                \
   xmk##r = _y1r + _y3r;                               \
   xmk##i = _y1i - _y3i;                               \
}

#ifndef USE_ASM
#define square_nested_1_4(xk,xmk,gkr,gki)              \
{                                                      \
   y_limb_t _r0,_r1,_r2,_r3;                           \
   _r0=(xk##r - xmk##r)*0.5;                           \
   _r1=(xk##i + xmk##i)*0.5;                           \
   _r2=(xk##r + xk##i) * (xk##r - xk##i);              \
   xk##i *= xk##r;                                     \
   xk##r = _r2;                                        \
   _r2 =(_r0 + _r1) * (_r0 - _r1);                     \
   _r1 *= _r0;                                         \
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);         \
   xmk##i *= xmk##r;                                   \
   xmk##r = _r3;                                       \
   _r1 +=  _r1;                                        \
   xk##i += xk##i;                                     \
   xmk##i += xmk##i;                                   \
   _r0= ((gki + 1.0) * _r2 + gkr * _r1);               \
   _r1= ((gki + 1.0) * _r1 - gkr * _r2);               \
   xk##r -= _r0;                                       \
   xmk##r -= _r0;                                      \
   xk##i -= _r1;                                       \
   xmk##i += _r1;                                      \
}

#else

#define square_nested_1_4(xk,xmk,gkr,gki)              \
__asm__ volatile("                                \n"  \
"fldl %0 	     #1   xkr                     \n"  \
"fsubl %2             #2-6 XR0                    \n"  \
"fldl %1              #3   xki,XR0                \n"  \
"faddl %3             #4-8 XR1,XR0                \n"  \
"fldl %0              #5   xkr,XR1,XR0            \n"  \
"faddl %1             #6-10 XK+,XR1,XR0           \n"  \
"fxch %%st(2)         #6    XR0,XR1,XK+           \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #7-11 R0,XR1,XK+            \n"  \
"fldl %0              #8    xkr,R0,XR1,XK+        \n"  \
"fsubl %1             #9-12 XK-,R0,XR1,XK+        \n"  \
"fxch %%st(2)         #9    XR1,R0,XK-,XK+        \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #10-14 R1,R0,XK-,XK+        \n"  \
"fldl %0              #11    xkr,R1,R0,XK-,kx+    \n"  \
"fmull %1             #12-16 RI,R1,r0,kx-,xk+     \n"  \
"fxch %%st(4)         #12    xk+,R1,r0,xk-,RI     \n"  \
"fmulp %%st,%%st(3)   #14-18 r1,r0,RR,RI          \n"  \
"fld %%st(1)          #15    r0,r1,r0,RR,RI       \n"  \
"fsub %%st(1)         #16-20 r-,r1,r0,RR,ri       \n"  \
"fld %%st(2)          #17    r0,R-,r1,r0,RR,ri    \n"  \
"fadd %%st(2)         #18-22 R+,R-,r1,r0,RR,ri    \n"  \
"fxch %%st(3)         #18    r0,R-,r1,R+,RR,ri    \n"  \
"fmulp %%st,%%st(2)   #19-23 R-,R2,R+,rr,ri       \n"  \
"fxch %%st(4)         #19    ri,R2,R+,rr,R-       \n"  \
"fadd %%st            #20-24 RI,R2,R+,rr,r-       \n"  \
"fxch %%st(3)         #21    rr,R2,R+,RI,r-       \n"  \
"fstpl %0             #22    R2,r+,RI,r-          \n"  \
"fxch %%st(1)         #22    r+,R2,RI,r-          \n"  \
"fmulp %%st(3)        #23-27 r2,RI,R1             \n"  \
"fxch %%st(1)         #24    ri,r2,R1             \n"  \
"fstpl %1             #25    r2,R1                \n"  \
"fadd %%st            #26-30 R2,R1                \n"  \
"fldl %2              #27    xmr,R2,R1            \n"  \
"faddl %3             #28-32 XM+,R2,r1            \n"  \
"fldl %2              #29    xmi,XM+,R2,r1        \n"  \
"fsubl %3             #30    XM-,XM+,R2,r1        \n"  \
"fldl %2              #31    xmr,XM-,XM+,r2,r1    \n"  \
"fmull %3             #32-36 RI,XM-,xm+,r2,r1     \n"  \
"fxch %%st(2)         #32    xm+,XM-,RI,r2,r1     \n"  \
"fmulp %%st,%%st(1)   #34-38 RR,RI,r2,r1          \n"  \
"fld1                 #35    1.0,RR,RI,r2,r1   ****\n" \
"faddl %5             #36-40 GR,RR,RI,r2,r1    ****\n" \
"fxch %%st(2)         #37    RI,RR,GR,r2,r1       \n"  \
"fadd %%st            #38-42 RI,rr,GR,r2,r1       \n"  \
"fxch %%st(1)         #39    rr,RI,GR,r2,r1       \n"  \
"fstpl %2             #40    RI,GR,r2,r1          \n"  \
"fld %%st(3)          #41    r1,RI,gr,r2,r1       \n"  \
"fmul %%st(2)         #42-46 RR,RI,gr,r2,r1       \n"  \
"fxch %%st(1)         #43    RI,RR,gr,r2,r1       \n"  \
"fstpl %3             #44    RR,gr,r2,r1          \n"  \
"fld %%st(2)          #45    r2,RR,gr,r2,r1       \n"  \
"fmull %4             #46    II,RR,gr,r2,r1   ****\n"  \
"fxch %%st(2)         #46    gr,RR,II,r2,r1       \n"  \
"fmulp %%st(3)        #47-51 rr,II,RI,r1          \n"  \
"faddp %%st,%%st(1)   #48-52 R0,RI,r1    **       \n"  \
"fxch %%st(2)         #48    r1,RI,R0             \n"  \
"fmull %4             #49-53 IR,RI,R0    **       \n"  \
"fldl %0              #50    xr,IR,RI,R0          \n"  \
"fsub %%st(3)         #51-55 XR,IR,RI,r0          \n"  \
"fldl %2              #52    xmr,XR,IR,ri,r0      \n"  \
"fsubp %%st,%%st(4)   #53-57 XR,ir,ri,XMR         \n"  \
"fxch %%st(1)         #53    ir,XR,ri,XMR         \n"  \
"fsubrp %%st,%%st(2)  #54-58 XR,R1,XMR            \n"  \
"fstpl %0             #55    R1,XMR               \n"  \
"fldl %1              #56    xi,R1,XMR            \n"  \
"fsub %%st(1)         #58-62 XI,R1,xmr            \n"  \
"fxch %%st(1)         #58    r1,XI,xmr            \n"  \
"faddl %3             #59-63 XMI,XI,xmr           \n"  \
"fxch %%st(2)         #60    xmr,XI,XMI           \n"  \
"fstpl %2             #61    xi,XMI               \n"  \
"fstpl %1             #62    XMI                  \n"  \
"fstpl %3             #63                         \n"  \
: :"m" (xk##r),"m"(xk##i),"m"(xmk##r),"m"(xmk##i),     \
"m"(gkr),"m"(gki): "memory");

#endif

#ifndef USE_ASM

#define square_nested_1_2(xk,xmk,gkr,gki)              \
{                                                      \
   y_limb_t _r0,_r1,_r2,_r3;                           \
   _r0=(xk##r - xmk##r)*0.5;                           \
   _r1=(xk##i + xmk##i)*0.5;                           \
   _r2=(xk##r + xk##i) * (xk##r - xk##i);              \
   xk##i *= xk##r;                                     \
   xk##r = _r2;                                        \
   _r2 =(_r0 + _r1) * (_r0 - _r1);                     \
   _r1 *= _r0;                                         \
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);         \
   xmk##i *= xmk##r;                                   \
   xmk##r = _r3;                                       \
   _r1 +=  _r1;                                        \
   xk##i += xk##i;                                     \
   xmk##i += xmk##i;                                   \
   _r0= ((gkr - 1.0) * _r2 - gki * _r1);               \
   _r1= ((gkr - 1.0) * _r1 + gki * _r2);               \
   xk##r += _r0;                                       \
   xmk##r += _r0;                                      \
   xk##i += _r1;                                       \
   xmk##i -= _r1;                                      \
}

#else
#define square_nested_1_2(xk,xmk,gkr,gki)              \
__asm__ volatile("                                \n"  \
"fldl %0 	     #1   xkr                     \n"  \
"fsubl %2             #2-6 XR0                    \n"  \
"fldl %1              #3   xki,XR0                \n"  \
"faddl %3             #4-8 XR1,XR0                \n"  \
"fldl %0              #5   xkr,XR1,XR0            \n"  \
"faddl %1             #6-10 XK+,XR1,XR0           \n"  \
"fxch %%st(2)         #6    XR0,XR1,XK+           \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #7-11 R0,XR1,XK+            \n"  \
"fldl %0              #8    xkr,R0,XR1,XK+        \n"  \
"fsubl %1             #9-12 XK-,R0,XR1,XK+        \n"  \
"fxch %%st(2)         #9    XR1,R0,XK-,XK+        \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #10-14 R1,R0,XK-,XK+        \n"  \
"fldl %0              #11    xkr,R1,R0,XK-,kx+    \n"  \
"fmull %1             #12-16 RI,R1,r0,kx-,xk+     \n"  \
"fxch %%st(4)         #12    xk+,R1,r0,xk-,RI     \n"  \
"fmulp %%st,%%st(3)   #14-18 r1,r0,RR,RI          \n"  \
"fld %%st(1)          #15    r0,r1,r0,RR,RI       \n"  \
"fsub %%st(1)         #16-20 r-,r1,r0,RR,ri       \n"  \
"fld %%st(2)          #17    r0,R-,r1,r0,RR,ri    \n"  \
"fadd %%st(2)         #18-22 R+,R-,r1,r0,RR,ri    \n"  \
"fxch %%st(3)         #18    r0,R-,r1,R+,RR,ri    \n"  \
"fmulp %%st,%%st(2)   #19-23 R-,R2,R+,rr,ri       \n"  \
"fxch %%st(4)         #19    ri,R2,R+,rr,R-       \n"  \
"fadd %%st            #20-24 RI,R2,R+,rr,r-       \n"  \
"fxch %%st(3)         #21    rr,R2,R+,RI,r-       \n"  \
"fstpl %0             #22    R2,r+,RI,r-          \n"  \
"fxch %%st(1)         #22    r+,R2,RI,r-          \n"  \
"fmulp %%st(3)        #23-27 r2,RI,R1             \n"  \
"fxch %%st(1)         #24    ri,r2,R1             \n"  \
"fstpl %1             #25    r2,R1                \n"  \
"fadd %%st            #26-30 R2,R1                \n"  \
"fldl %2              #27    xmr,R2,R1            \n"  \
"faddl %3             #28-32 XM+,R2,r1            \n"  \
"fldl %2              #29    xmi,XM+,R2,r1        \n"  \
"fsubl %3             #30    XM-,XM+,R2,r1        \n"  \
"fldl %2              #31    xmr,XM-,XM+,r2,r1    \n"  \
"fmull %3             #32-36 RI,XM-,xm+,r2,r1     \n"  \
"fxch %%st(2)         #32    xm+,XM-,RI,r2,r1     \n"  \
"fmulp %%st,%%st(1)   #34-38 RR,RI,r2,r1          \n"  \
"fld1                 #35    1.0,RR,RI,r2,r1   ****\n" \
"fsubrl %4            #36-40 GR,RR,RI,r2,r1    ****\n" \
"fxch %%st(2)         #37    RI,RR,GR,r2,r1       \n"  \
"fadd %%st            #38-42 RI,rr,GR,r2,r1       \n"  \
"fxch %%st(1)         #39    rr,RI,GR,r2,r1       \n"  \
"fstpl %2             #40    RI,GR,r2,r1          \n"  \
"fld %%st(3)          #41    r1,RI,gr,r2,r1       \n"  \
"fmul %%st(2)         #42-46 RR,RI,gr,r2,r1       \n"  \
"fxch %%st(1)         #43    RI,RR,gr,r2,r1       \n"  \
"fstpl %3             #44    RR,gr,r2,r1          \n"  \
"fld %%st(2)          #45    r2,RR,gr,r2,r1       \n"  \
"fmull %5             #46    II,RR,gr,r2,r1   ****\n"  \
"fxch %%st(2)         #46    gr,RR,II,r2,r1       \n"  \
"fmulp %%st(3)        #47-51 rr,II,RI,r1          \n"  \
"fsubp %%st,%%st(1)   #48-52 R0,RI,r1    **       \n"  \
"fxch %%st(2)         #48    r1,RI,R0             \n"  \
"fmull %5             #49-53 IR,RI,R0    **       \n"  \
"fldl %0              #50    xr,IR,RI,R0          \n"  \
"fadd %%st(3)         #51-55 XR,IR,RI,r0          \n"  \
"fldl %2              #52    xmr,XR,IR,ri,r0      \n"  \
"faddp %%st,%%st(4)   #53-57 XR,ir,ri,XMR         \n"  \
"fxch %%st(1)         #53    ir,XR,ri,XMR         \n"  \
"faddp %%st(2)        #54-58 XR,R1,XMR            \n"  \
"fstpl %0             #55    R1,XMR               \n"  \
"fldl %1              #56    xi,R1,XMR            \n"  \
"fadd %%st(1)         #58-62 XI,R1,xmr            \n"  \
"fxch %%st(1)         #58    r1,XI,xmr            \n"  \
"fsubrl %3            #59-63 XMI,XI,xmr           \n"  \
"fxch %%st(2)         #60    xmr,XI,XMI           \n"  \
"fstpl %2             #61    xi,XMI               \n"  \
"fstpl %1             #62    XMI                  \n"  \
"fstpl %3             #63                         \n"  \
: :"m" (xk##r),"m"(xk##i),"m"(xmk##r),"m"(xmk##i),     \
"m"(gkr),"m"(gki): "memory");

#endif

#ifndef USE_ASM

#define square_nested_3_4(xk,xmk,gkr,gki)              \
{                                                      \
   y_limb_t _r0,_r1,_r2,_r3;                           \
   _r0=(xk##r - xmk##r)*0.5;                           \
   _r1=(xk##i + xmk##i)*0.5;                           \
   _r2=(xk##r + xk##i) * (xk##r - xk##i);              \
   xk##i *= xk##r;                                     \
   xk##r = _r2;                                        \
   _r2 =(_r0 + _r1) * (_r0 - _r1);                     \
   _r1 *= _r0;                                         \
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);         \
   xmk##i *= xmk##r;                                   \
   xmk##r = _r3;                                       \
   _r1 +=  _r1;                                        \
   xk##i += xk##i;                                     \
   xmk##i += xmk##i;                                   \
   _r0= ((gki - 1.0) * _r2 + gkr * _r1);               \
   _r1= ((gki - 1.0) * _r1 - gkr * _r2);               \
   xk##r += _r0;                                       \
   xmk##r += _r0;                                      \
   xk##i += _r1;                                       \
   xmk##i -= _r1;                                      \
}

#else

#define square_nested_3_4(xk,xmk,gkr,gki)              \
__asm__ volatile("                                \n"  \
"fldl %0 	     #1   xkr                     \n"  \
"fsubl %2             #2-6 XR0                    \n"  \
"fldl %1              #3   xki,XR0                \n"  \
"faddl %3             #4-8 XR1,XR0                \n"  \
"fldl %0              #5   xkr,XR1,XR0            \n"  \
"faddl %1             #6-10 XK+,XR1,XR0           \n"  \
"fxch %%st(2)         #6    XR0,XR1,XK+           \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #7-11 R0,XR1,XK+            \n"  \
"fldl %0              #8    xkr,R0,XR1,XK+        \n"  \
"fsubl %1             #9-12 XK-,R0,XR1,XK+        \n"  \
"fxch %%st(2)         #9    XR1,R0,XK-,XK+        \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_6r) "         #10-14 R1,R0,XK-,XK+        \n"  \
"fldl %0              #11    xkr,R1,R0,XK-,kx+    \n"  \
"fmull %1             #12-16 RI,R1,r0,kx-,xk+     \n"  \
"fxch %%st(4)         #12    xk+,R1,r0,xk-,RI     \n"  \
"fmulp %%st,%%st(3)   #14-18 r1,r0,RR,RI          \n"  \
"fld %%st(1)          #15    r0,r1,r0,RR,RI       \n"  \
"fsub %%st(1)         #16-20 r-,r1,r0,RR,ri       \n"  \
"fld %%st(2)          #17    r0,R-,r1,r0,RR,ri    \n"  \
"fadd %%st(2)         #18-22 R+,R-,r1,r0,RR,ri    \n"  \
"fxch %%st(3)         #18    r0,R-,r1,R+,RR,ri    \n"  \
"fmulp %%st,%%st(2)   #19-23 R-,R2,R+,rr,ri       \n"  \
"fxch %%st(4)         #19    ri,R2,R+,rr,R-       \n"  \
"fadd %%st            #20-24 RI,R2,R+,rr,r-       \n"  \
"fxch %%st(3)         #21    rr,R2,R+,RI,r-       \n"  \
"fstpl %0             #22    R2,r+,RI,r-          \n"  \
"fxch %%st(1)         #22    r+,R2,RI,r-          \n"  \
"fmulp %%st(3)        #23-27 r2,RI,R1             \n"  \
"fxch %%st(1)         #24    ri,r2,R1             \n"  \
"fstpl %1             #25    r2,R1                \n"  \
"fadd %%st            #26-30 R2,R1                \n"  \
"fldl %2              #27    xmr,R2,R1            \n"  \
"faddl %3             #28-32 XM+,R2,r1            \n"  \
"fldl %2              #29    xmi,XM+,R2,r1        \n"  \
"fsubl %3             #30    XM-,XM+,R2,r1        \n"  \
"fldl %2              #31    xmr,XM-,XM+,r2,r1    \n"  \
"fmull %3             #32-36 RI,XM-,xm+,r2,r1     \n"  \
"fxch %%st(2)         #32    xm+,XM-,RI,r2,r1     \n"  \
"fmulp %%st,%%st(1)   #34-38 RR,RI,r2,r1          \n"  \
"fld1                 #35    1.0,RR,RI,r2,r1   ****\n" \
"fsubrl %5            #36-40 GR,RR,RI,r2,r1    ****\n" \
"fxch %%st(2)         #37    RI,RR,GR,r2,r1       \n"  \
"fadd %%st            #38-42 RI,rr,GR,r2,r1       \n"  \
"fxch %%st(1)         #39    rr,RI,GR,r2,r1       \n"  \
"fstpl %2             #40    RI,GR,r2,r1          \n"  \
"fld %%st(3)          #41    r1,RI,gr,r2,r1       \n"  \
"fmul %%st(2)         #42-46 RR,RI,gr,r2,r1       \n"  \
"fxch %%st(1)         #43    RI,RR,gr,r2,r1       \n"  \
"fstpl %3             #44    RR,gr,r2,r1          \n"  \
"fld %%st(2)          #45    r2,RR,gr,r2,r1       \n"  \
"fmull %4             #46    II,RR,gr,r2,r1   ****\n"  \
"fxch %%st(2)         #46    gr,RR,II,r2,r1       \n"  \
"fmulp %%st(3)        #47-51 rr,II,RI,r1          \n"  \
"faddp %%st,%%st(1)   #48-52 R0,RI,r1    **       \n"  \
"fxch %%st(2)         #48    r1,RI,R0             \n"  \
"fmull %4             #49-53 IR,RI,R0    **       \n"  \
"fldl %0              #50    xr,IR,RI,R0          \n"  \
"fadd %%st(3)         #51-55 XR,IR,RI,r0          \n"  \
"fldl %2              #52    xmr,XR,IR,ri,r0      \n"  \
"faddp %%st,%%st(4)   #53-57 XR,ir,ri,XMR         \n"  \
"fxch %%st(1)         #53    ir,XR,ri,XMR         \n"  \
"fsubrp %%st,%%st(2)  #54-58 XR,R1,XMR            \n"  \
"fstpl %0             #55    R1,XMR               \n"  \
"fldl %1              #56    xi,R1,XMR            \n"  \
"fadd %%st(1)         #58-62 XI,R1,xmr            \n"  \
"fxch %%st(1)         #58    r1,XI,xmr            \n"  \
"fsubrl %3            #59-63 XMI,XI,xmr           \n"  \
"fxch %%st(2)         #60    xmr,XI,XMI           \n"  \
"fstpl %2             #61    xi,XMI               \n"  \
"fstpl %1             #62    XMI                  \n"  \
"fstpl %3             #63                         \n"  \
: :"m" (xk##r),"m"(xk##i),"m"(xmk##r),"m"(xmk##i),     \
"m"(gkr),"m"(gki): "memory");

#endif

/* $Id$ */
