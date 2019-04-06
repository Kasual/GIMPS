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
    gbvvalor@worldonline.es
*/
/* This are macros to be used with complex FFT                        */

#define Y_REG register y_limb_t

/***************  SINGLE LOAD AND STORE MACROS ***************/
/* To load from real padded array of data as a complex array */

#define cplx_data_to_local(_t,_array,_k) \
	_t##.re = _array[addr(( _k )<<1)]; \
	_t##.im = _array[addr(( _k )<<1)+1];

#define cplx_data_to_local_p(_t,_pd,_array,_k) \
	_pd = addr(( _k)<<1);\
	_t##.re = _array[ _pd ]; \
	_t##.im = _array[ _pd + 1];


/* To store on read padded array from local data as complex */

#define cplx_local_to_data(_array,_k,_t) \
	_array[addr(( _k )<<1)]=_t##.re ; \
	_array[addr(( _k )<<1)+1]=_t##.im;

#define cplx_local_to_data_p(_array,_pd,_t) \
	_array[ _pd ] = _t##.re ; \
	_array[ _pd + 1] = _t##.im;


/* To load from real non-padded array of data as a complex array */
#define cplx_mem_to_local(_t,_array,_k) \
	_t##.re =_array[( _k )<<1]; \
	_t##.im =_array[(( _k )<<1)+1];


/* To strore to a non-padded array from local */
#define cplx_local_to_mem(_array,_k,_t) \
	_array[( _k )<<1]=_t##.re; \
	_array[(( _k )<<1)+1]=_t##.im;

/**************** PSEUDO_COMPLEX GENERAL AUXILIAR MACROS **************/
/**********************************************************************
   cplx_muladdsub macro :
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex and _t is temporal aux var.
***********************************************************************/

#ifndef USE_ASM
# define cplx_muladdsub(_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai,_a=_t0##.re; \
		_ar = (_t1##.re ) * (_f##.re ) - (_t1##.im ) * (_f##.im ); \
		_ai = (_t1##.re ) * (_f##.im ) + (_t1##.im ) * (_f##.re ); \
		_t0##.re += _ar;\
		_t1##.re = _a - _ar;\
		_a = _t0##.im;\
		_t0##.im += _ai;\
		_t1##.im = _a - _ai;\
	}

#else

# define cplx_muladdsub(_t0,_t1,_f) \
   __asm__( "
fldl %2                      #1     t1r
fmull %4                     #2-6   RR
fldl %3                      #3     t1i,RR
fmull %5                     #4-8   II,RR
fldl %2                      #5     t1r,II,RR
fmull %5                     #6-10  RI,II,RR
fldl %3                      #7     t1i,RI,II,rr
fmull %4                     #8-12  IR,RI,II,rr
fxch %%st(3)                 #8     rr,RI,II,IR
fsubp %%st,%%st(2)           #9-12  RI,R1,IR
fldl %0                      #10    t0r,RI,R1,IR
fldl %0                      #11    t0r,t0r,ri,R1,IR
fxch %%st(2)                 #11    ri,t0r,t0r,R1,IR
faddp %%st,%%st(4)           #12-16 t0r,t0r,R1,I1
fadd %%st(2)                 #13-17 NT0R,t0r,r1,I1
fxch %%st(1)                 #13    t0r,NT0R,r1,I1
fsubp %%st,%%st(2)           #14-18 NT0R,NT1R,I1
fldl %1                      #15    t0i,NT0R,NT1R,I1
fadd %%st(3)                 #16-20 NT0I,NT0R,NT1R,I1
fldl %1                      #17    t0i,NT0I,NT0R,NT1R,I1
fsubp %%st,%%st(4)           #18-22 NT0I,nt0r,NT1R,NT1I
fxch %%st(2)                 #19    nt1r,nt0r,NT0I,NT1I
fstpl %2                     #20    nt0r,NT0I,NT1I
fstpl %0                     #21    nt0i,NT1I
fstpl %1                     #22    nt1i
fstpl %3                     #23
"\
: :"m"( _t0##.re ),"m"( _t0##.im ), "m"(_t1##.re), "m"(_t1##.im),\
"m"( _f##.re ), "m" ( _f##.im ): "memory");

#endif




/**********************************************************************
   cplx_load_muladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
 where all the variables are pseudo-complex and _t is temporal aux var.
**********************************************************************/
#define cplx_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
		_ar = (_br) * (_f##.re) - (_bi) * (_f##.im); \
		_ai = (_br) * (_f##.im) + (_bi) * (_f##.re); \
		_t0##.re = _array[addr(_k0<<1)] + _ar; \
		_t1##.re = _array[addr(_k0<<1)] - _ar; \
		_t0##.im = _array[addr(_k0<<1)+1] + _ai; \
		_t1##.im = _array[addr(_k0<<1)+1] - _ai; \
	}

#ifndef USE_ASM
#define cplx_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_pd1= addr((_k1)<<1);\
		_br= _array[ _pd1 ];\
		_bi= _array[ _pd1 + 1];\
		_ar = (_br) * (_f##.re) - (_bi) * (_f##.im); \
		_pd0 = addr((_k0)<<1);\
		_ai = (_br) * (_f##.im) + (_bi) * (_f##.re); \
		_t0##.re = _array[ _pd0 ] + _ar; \
		_t1##.re = _array[ _pd0 ] - _ar; \
		_t0##.im = _array[ _pd0 + 1] + _ai; \
		_t1##.im = _array[ _pd0 + 1] - _ai; \
	}
#else

#define cplx_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
__asm__("
movl %%ebx,%%eax             # 1  U
andl $-256,%%ebx             # 1  V
shrl $7,%%ebx                # 2  U
addl %%eax,%%eax	     # 2  V
addl %%eax,%%ebx	     # 3  U
fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache
fldl 8(%%esi,%%ebx,8)        # 5 t1i,t1r
  fldl %2     		     # 6 fr,t1i,t1r
  fmul %%st(1)                 # 7-11 IR,t1i,t1r
  fldl %3                      # 8    fi,IR,t1i,t1r
  fmul %%st(3)                 # 9-13 RI,IR,t1i,t1r
  fldl %2                      # 10   fr,RI,IR,t1i,t1r
  fmulp %%st,%%st(4)           # 11-15   RI,IR,t1i,RR
  fldl %3			     # 12   fi,RI,ir,t1i,RR
  fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR
  movl %%edx,%%eax             # 14  U
  andl $-256,%%edx             # 14  V
  shrl $7,%%edx                # 15  U
  addl %%eax,%%eax	     # 15  V
  addl %%eax,%%edx	     # 16  U
  faddp %%st,%%st(1)           # 17-21  NT1I,II,rr
  fxch %%st(2)                 # 17     rr,II,NT1I
  fldl (%%esi,%%edx,8)         # 18     nt0r,rr,ii,NT1I
  fxch %%st(1)                 # 18     rr,nt0r,ii,NT1I
  fsubp %%st,%%st(2)           # 19-23  nt0r,NT1R,NT1I
  fldl 8(%%esi,%%edx,8)        # 20     nt0i,nt0r,NT1R,NT1I
  fsub %%st(3)                 # 21-25  T1I,nt0r,NT1R,nt1i
  fxch %%st(3)                 # 21     nt1i,nt0r,NT1R,T1I
  faddl 8(%%esi,%%edx,8)       # 22-26  T0I,nt0r,nt1r,T1I
  fxch %%st(1)                 # 22     nt0r,T0I,nt1r,T1I
  fsub %%st(2)                 # 23-27  T1R,T0I,nt1r,T1I
  fxch %%st(2)                 # 23     nt1r,T0I,T1R,T1I
  faddl (%%esi,%%edx,8)        # 24-28  T0R,T0I,T1R,t1i
  fxch %%st(3)                 # 25     t1i,T0I,T1R,T0R
  "\
  : "=&b"(_pd1) ,"=&d" (_pd0): "m"(_f##.re),"m"(_f##.im),"S"(_array),\
                   "0"(_k1),"1"(_k0) :"ax","st");
                   \
                   \
                   __asm__("
                   fstpl  %3           # 27
                   fstpl  %1           # 28
                   fstpl  %2           # 29
                   fstpl  %0           # 30
                   "\
                   : :"m"( _t0##.re ),"m"( _t0##.im ),"m"(_t1##.re),"m"(_t1##.im):\
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
#define cplx_muladdsub_store(_array,_k0,_k1,_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##.re) * (_f##.re) - (_t1##.im) * (_f##.im); \
		_ai = (_t1##.re) * (_f##.im) + (_t1##.im) * (_f##.re); \
		_array[addr(_k0<<1)]= _t0##.re + _ar; \
		_array[addr(_k1<<1)]= _t0##.re - _ar; \
		_array[addr(_k0<<1)+1]= _t0##.im + _ai; \
		_array[addr(_k1<<1)+1]= _t0##.im - _ai; \
        }

#ifndef USE_ASM
# define cplx_muladdsub_store_p(_array,_pd0,_pd1,_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##.re) * (_f##.re) - (_t1##.im) * (_f##.im); \
		_ai = (_t1##.re) * (_f##.im) + (_t1##.im) * (_f##.re); \
		_array[ _pd0 ] = _t0##.re + _ar; \
		_array[ _pd1 ] = _t0##.re - _ar; \
		_array[ _pd0 + 1] = _t0##.im + _ai; \
		_array[ _pd1 + 1] = _t0##.im - _ai; \
	}
#else
# define cplx_muladdsub_store_p(_array,_pd0,_pd1,_t0,_t1,_f) \
__asm__("
fldl %7                      # 1      t1r
fldl %8                      # 2      t1i,t1r
fldl %2     		     # 3      fr,t1i,t1r
fmul %%st(1)                 # 4-8    IR,t1i,t1r
fldl %3                      # 5      fi,IR,t1i,t1r
fmul %%st(3)                 # 6-10   RI,IR,t1i,t1r
fldl %2                      # 7      fr,RI,IR,t1i,t1r
fmulp %%st,%%st(4)           # 8-12   RI,IR,t1i,RR
fldl %3			     # 13     fi,RI,ir,t1i,RR
fmulp %%st,%%st(3)           # 14-18  RI,ir,II,RR
faddp %%st,%%st(1)           # 15-19  NT1I,II,rr
fxch %%st(2)                 # 15     rr,II,NT1I
fldl %5                      # 16     nt0r,rr,ii,NT1I
fxch %%st(1)                 # 16     rr,nt0r,ii,NT1I
fsubp %%st,%%st(2)           # 17-21  nt0r,NT1R,NT1I
fldl %6                      # 18     nt0i,nt0r,NT1R,NT1I
fsub %%st(3)                 # 19-23  T1I,nt0r,NT1R,nt1i
fxch %%st(3)                 # 19     nt1i,nt0r,NT1R,T1I
faddl %6                     # 20-24  T0I,nt0r,nt1r,T1I
fxch %%st(1)                 # 20     nt0r,T0I,nt1r,T1I
fsub %%st(2)                 # 21-25  T1R,T0I,nt1r,T1I
fxch %%st(2)                 # 21     nt1r,T0I,T1R,T1I
faddl %5                     # 22-26  T0R,T0I,T1R,t1i
fxch %%st(3)                 # 23     t1i,T0I,T1R,T0R
fstpl 8(%%esi,%%ebx,8)       # 24
fstpl 8(%%esi,%%edx,8)       # 25
fstpl (%%esi,%%ebx,8)        # 26
fstpl (%%esi,%%edx,8)        # 27
"\
: :"b"(_pd1) ,"d" (_pd0), "m"(_f##.re),"m"(_f##.im),"S"(_array),\
"m"(_t0##.re),"m"(_t0##.im),"m"(_t1##.re),"m"(_t1##.im) :"memory");

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
# define cplx_mulmuladdsub(_t0,_t1,_f0,_f1) \
	{ \
	        y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##.re) * (_f1##.re) - (_t1##.im) * (_f1##.im); \
		_br = (_t0##.re) * (_f0##.re) - (_t0##.im) * (_f0##.im); \
		_ai = (_t1##.re) * (_f1##.im) + (_t1##.im) * (_f1##.re); \
		_bi = (_t0##.re) * (_f0##.im) + (_t0##.im) * (_f0##.re); \
		_t0##.re = _br + _ar;\
		_t1##.re = _br - _ar;\
		_t0##.im = _bi + _ai;\
		_t1##.im = _bi - _ai;\
	}
#else
#define cplx_mulmuladdsub(_t0,_t1,_f0,_f1) \
   __asm__( "
fldl %2                      #1-1 t1r
fmull %6                     #2-6 rr
fldl %3                      #3-3 t1i,rr
fmull %7                     #4-8 ii,rr
fldl %2                      #5-5 t1r,ii,rr
fmull %7                     #6-10 ri,ii,rr
fldl %3                      #7-7  t1i,ri,ii,rr
fmull %6                     #8-12  ir,ri,ii,rr
fxch %%st(3)                 #8-8   rr,ri,ii,ir
fsubp %%st,%%st(2)           #9-12  ri,r1,ir
fldl %0                      #10-10  t0r,ri,r1,ir
fmull %4                     #11-15  rr0,ri,r1,ir
fldl %1                      #12-12  t0i,rr0,ri,r1,ir
fmull %5                     #13-17  ii0,rr0,ri,r1,ir
fxch %%st(2)                 #13-13  ri,rr0,ii0,r1,ir
faddp %%st,%%st(4)           #14-17  rr0,ii0,r1,i1
fldl %0                      #15-15  t0r,rr0,ii0,r1,i1
fmull %5                     #16-20  ri0,rr0,ii0,r1,i1
fldl %1                      #17-17  t0i,ri0,rr0,ii0,r1,i1
fmull %4                     #18-22  ir0,ri0,rr0,ii0,r1,i1
fxch %%st(2)                 #18-18  rr0,ri0,ir0,ii0,r1,i1
fsubp %%st,%%st(3)           #19-22  ri0,ir0,r0,r1,i1
fld %%st(3)                  #20-20  r1,ri0,ir0,r0,r1,i1
fld %%st(5)                  #21-21  i1,r1,ri0,ir0,r0,r1,i1
fxch %%st(4)                 #21-21  r0,r1,ri0,ir0,i1,r1,i1
fadd %%st,%%st(5)            #22-25  r0,r1,ri0,ir0,i1,nt0r,i1
fxch %%st(2)                 #22-22  ri0,r1,r0,ir0,i1,nt0r,i1
faddp %%st,%%st(3)           #23-25  r1,r0,i0,i1,nt0r,i1
fxch %%st(1)                 #23-23  r0,r1,i0,i1,nt0r,i1
fsubp %%st,%%st(1)           #24-27  nt1r,i0,i1,nt0r,i1
fxch %%st(1)                 #24-24  i0,nt1r,i1,nt0r,i1
fadd %%st,%%st(4)            #25-28  i0,nt1r,i1,nt0r,nt0i
fsubp %%st,%%st(2)           #26-29  nt1r,nt1i,nt0r,nt0i
fxch %%st(2)                 #27     nt0r,nt1i,nt1r,nt0i
fstpl %0                     #28     nt1i,nt1r,nt0i
fxch %%st(1)                 #29     nt1r,nt1i,nt0i
fstpl %2                     #30     nt1i,nt0i
fxch %%st(1)                 #31     nt0i,nt1i
fstpl %1                     #32     nt1i
fstpl %3                     #33
"\
: :"m"( _t0##.re ),"m"( _t0##.im ),"m"(_t1##.re ),"m"(_t1##.im ),\
"m"( _f0##.re ), "m" ( _f0##.im ),"m" ( _f1##.re ),"m" ( _f1##.im):\
"memory");

#endif




#define cplx_mulmul(_t0,_t1,_f0,_f1) \
	{ \
		y_limb_t _ar,_br; \
		_ar = (_t1##.re) * (_f1##.re) - (_t1##.im) * (_f1##.im); \
		_br = (_t0##.re) * (_f0##.re) - (_t0##.im) * (_f0##.im); \
		_t1##.im = (_t1##.re) * (_f1##.im) + (_t1##.im) * (_f1##.re); \
		_t0##.im = (_t0##.re) * (_f0##.im) + (_t0##.im) * (_f0##.re); \
		_t1##.re = _ar;\
		_t0##.re = _br;\
	}


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
#define cplx_load_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= _array[addr(_k0<<1)];\
		_bi= _array[addr(_k0<<1)+1];\
		_ar = (_br) * (_f0##.re) - (_bi) * (_f0##.im); \
		_ai = (_br) * (_f0##.im) + (_bi) * (_f0##.re); \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
		_cr = (_br) * (_f1##.re) - (_bi) * (_f1##.im); \
		_ci = (_br) * (_f1##.im) + (_bi) * (_f1##.re); \
		_t0##.re = _ar + _cr; \
		_t1##.re = _ar - _cr; \
		_t0##.im = _ai + _ci; \
		_t1##.im = _ai - _ci; \
	}

#ifndef USE_ASM
# define cplx_load_mulmuladdsub_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_pd0 = addr((_k0)<<1);\
		_br= _array[ _pd0 ];\
		_bi= _array[ _pd0 + 1];\
		_ar = (_br) * (_f0##.re) - (_bi) * (_f0##.im); \
		_pd1 = addr((_k1)<<1);\
		_ai = (_br) * (_f0##.im) + (_bi) * (_f0##.re); \
		_br= _array[ _pd1 ];\
		_bi= _array[ _pd1 + 1];\
		_cr = (_br) * (_f1##.re) - (_bi) * (_f1##.im); \
		_ci = (_br) * (_f1##.im) + (_bi) * (_f1##.re); \
		_t0##.re = _ar + _cr; \
		_t1##.re = _ar - _cr; \
		_t0##.im = _ai + _ci; \
		_t1##.im = _ai - _ci; \
	}

#else
# define cplx_load_mulmuladdsub_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
                __asm__("
movl %%ebx,%%eax             # 1  U
andl $-256,%%ebx             # 1  V
shrl $7,%%ebx                # 2  U
addl %%eax,%%eax	     # 2  V
addl %%eax,%%ebx	     # 3  U
fldl (%%esi,%%ebx,8)         # 4 (14) if not in cache
fldl 8(%%esi,%%ebx,8)        # 5      t1i,t1r
  fldl %2     		     # 6      fr,t1i,t1r
  fmul %%st(1)                 # 7-11   IR,t1i,t1r
  fldl %3                      # 8      fi,IR,t1i,t1r
  fmul %%st(3)                 # 9-13   RI,IR,t1i,t1r
  fldl %2                      # 10     fr,RI,IR,t1i,t1r
  fmulp %%st,%%st(4)           # 11-15  RI,IR,t1i,RR
  fldl %3			     # 12     fi,RI,ir,t1i,RR
  fmulp %%st,%%st(3)           # 13-17  RI,ir,II,RR
  movl %%edx,%%eax             # 14  U
  andl $-256,%%edx             # 14  V
  shrl $7,%%edx                # 15  U
  addl %%eax,%%eax	     # 15  V
  addl %%eax,%%edx	     # 16  U
  faddp %%st,%%st(1)           # 17-21  NT1I,II,rr
  fxch %%st(2)                 # 17     rr,II,NT1I
  fldl (%%esi,%%edx,8)         # 18     nt0r,rr,ii,NT1I
  fmull %4                     # 19-23  RR0,rr,ii,NT1I
  fxch %%st(1)                 # 19     rr,RR0,ii,NT1I
  fsubp %%st,%%st(2)           # 20-24  RR0,NT1R,NT1I
  fldl 8(%%esi,%%edx,8)        # 21     t0i,RR,NT1R,NT1I
  fmull %5                     # 22-26  II,RR,NT1R,nt1i
  fldl (%%esi,%%edx,8)         # 23     t0r,II,RR,NT1R,nt1i
  fmull %5                     # 24-28  RI,II,rr,NT1R,nt1i
  fldl 8(%%esi,%%edx,8)        # 25     t0i,RI,II,rr,nt1r,nt1i
  fmull %4                     # 26-30  IR,RI,II,rr,nt1r,nt1i
  fxch %%st(3)                 # 26     rr,RI,II,IR,nt1r,nt1i
  fsubp %%st,%%st(2)           # 27-31  RI,NTOR,IR,nt1r,nt1i
  fxch %%st(2)                 # 27     IR,NTOR,RI,nt1r,nt1i
  fld  %%st(4)                 # 28     nt1i,IR,NTOR,RI,nt1r,nt1i
  fld  %%st(4)                 # 29     nt1r,nt1i,IR,NTOR,RI,nt1r,nt1i
  fadd %%st(3)                 # 30-34  T0R,nt1i,IR,NTOR,ri,nt1r,nt1i
  fxch %%st(4)                 # 30     ri,nt1i,IR,NTOR,TOR,nt1r,nt1i
  faddp %%st,%%st(2)           # 31-35  nt1i,NTOI,ntor,TOR,nt1r,nt1i
  fxch %%st(2)                 # 31     nt0r,NTOI,nt1i,TOR,nt1r,nt1i
  fsubp %%st,%%st(4)           # 32-36  NTOI,nt1i,TOR,T1R,nt1i
  fadd  %%st,%%st(1)           # 33-37  ntoi,TOI,TOR,T1R,nt1i
  fsubp %%st,%%st(4)           # 34-38  TOI,TOR,T1R,T1I
  "\
  : "=&b"(_pd1) ,"=&d" (_pd0): "m"(_f1##.re),"m"(_f1##.im),"m"(_f0##.re),\
                   "m"(_f0##.im),"S"(_array),"0"(_k1),"1"(_k0) :"ax","st");\
                   \
                   __asm__("
                   fstpl  %1           # 37
                   fstpl  %0           # 38
                   fstpl  %2           # 39
                   fstpl  %3           # 40
                   "\
                   : :"m"( _t0##.re ),"m" ( _t0##.im ),"m"(_t1##.re),"m"(_t1##.im):\
                   "memory", "st");

#endif


#define cplx_load_mul_p(_t,_pd,_f,_array,_k)\
{\
	y_limb_t _ar,_ai;\
	_pd = addr((_k)<<1);\
	_ar = _array[ _pd ];\
	_ai = _array[ _pd + 1];\
	_t##.re = _ar * _f##.re - _ai * _f##.im;\
	_t##.im = _ai * _f##.re + _ar * _f##.im;\
}

#define cplx_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_b,_cr,_ci; \
		_pd0 = addr((_k0)<<1);\
		_b = _array[ _pd0 ];\
		_ai= _array[ _pd0 + 1];\
		_ar = (_b) * (_f0##.re) - (_ai) * (_f0##.im); \
		_pd1 = addr((_k1)<<1);\
		_ai = (_b) * (_f0##.im) + (_ai) * (_f0##.re); \
		_b = _array[ _pd1 ];\
		_ci = _array[ _pd1 + 1];\
		_cr = (_b) * (_f1##.re) - (_ci) * (_f1##.im); \
		_ci = (_b) * (_f1##.im) + (_ci) * (_f1##.re); \
		_t0##.re = _ar; \
		_t0##.im = _ai; \
		_t1##.re = _cr; \
		_t1##.im = _ci; \
	}

/**********************************************************************
   cplx_divmul macro:
		_t0= _f0 / _f1;
		_t1= _f0 * _f1;
   where all the variables  are pseudo-complex. _f0 and _f1 are module
   1 pseudo-complex (like trig. factros)
**********************************************************************/
#ifndef USE_ASM
#define cplx_divmul(_t0,_t1,_f0,_f1)\
	{\
		y_limb_t _w0,_w1,_w2,_w3;\
		_w0 = _f0##.re * _f1##.re;\
		_w1 = _f0##.im * _f1##.im;\
		_w2 = _f0##.im * _f1##.re;\
		_w3 = _f0##.re * _f1##.im;\
		_t0##.re = _w0 + _w1;\
		_t1##.re = _w0 - _w1;\
		_t0##.im = _w2 - _w3;\
		_t1##.im = _w2 + _w3;\
	}
#else
#define cplx_divmul(_t0,_t1,_f0,_f1)\
	__asm__("
fldl %4	       		# 1     f0r
fmull %6                # 2-6   RR
fldl %5                 # 3     f0i,RR
fmull %7                # 4-8   II,RR
fldl %5                 # 5     f0i,II,RR
fmull %6                # 6-10  IR,II,RR
fldl %4                 # 7     f0r,IR,II,rr
fmull %7                # 8-12  RI,IR,II,rr
fld %%st(3)             # 9     rr,RI,IR,ii,rr
fsub %%st(3)            # 10-14 T1R,RI,IR,ii,rr
fxch %%st(3)            # 10    ii,RI,IR,T1R,rr
faddp %%st,%%st(4)      # 11-15 RI,ir,T1R,T0R
fld %%st(1)             # 12    ir,RI,ir,T1R,T0R
fsub %%st(1)            # 13-17 T0I,ri,ir,T1R,T0R
fxch %%st(2)            # 13    ir,ri,T0I,T1R,T0R
faddp %%st,%%st(1)      # 14-18 T1I,T0I,T1R,T0R
fxch %%st(2)            # 15    t1r,T0I,t1I,T0R
fstpl %2                # 16    T0I,t1I,t0r
fstpl %1                # 17    T1I,t0r
fstpl %3                # 18    t0r
fstpl %0                # 19
"\
: :"m"(_t0##.re),"m"(_t0##.im),"m"(_t1##.re),"m"(_t1##.im),\
"m"(_f0##.re),"m"(_f0##.im),"m"(_f1##.re),"m"(_f1##.im):"memory");
#endif


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/
#ifndef USE_ASM
# define cplx_addsub(_t0,_t1) \
	{ \
		y_limb_t _a=_t0##.re; \
		_t0##.re += (_t1##.re); \
		_t1##.re = _a - (_t1##.re); \
		_a = (_t0##.im); \
		_t0##.im += (_t1##.im); \
		_t1##.im = _a - (_t1##.im); \
	}
#else
# define cplx_addsub(_t0,_t1) \
__asm__( "
fldl %0             #1-1 t0r
faddl %2            #2-5 NT0R
fldl %0             #3-3 t0r,NT0R
fsubl %2            #4-7 NT1R,NT0R
fldl %1             #5-5 t0i,NT1R,NT0R
faddl %3            #6-9 NT0I,NT1R,ntor
fxch %%st(2)        #6   nt0r,NT1R,NT0I
fldl %1             #7   t0i,nt0r,NT1R,NT0I
fsubl %3            #8-11 NT1I,nt0r,NT1R,NTOI
fxch %%st(3)        #9    NTOI,nt0r,nt1R,NT1I
fstpl %1            #10   nt0r,NT1R,NT1I
fstpl %0            #11   NT1R,NT1I
fstpl %2            #12   nt1I
fstpl %3            #13   void
"\
: :"m"(_t0##.re),"m"(_t0##.im ),"m"(_t1##.re),"m"(_t1##.im):\
"memory");
#endif



/**********************************************************************
   cplx_load_addsub macro:
		_t0 = _array( _k0) + _array (_k1);
		_t1 = _array( _k0) - _array (_k1);
 where all the variables are pseudo-complex .
**********************************************************************/

#define cplx_load_addsub(_t0,_t1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = _array[addr(_k1<<1)];\
		_ai = _array[addr(_k1<<1)+1];\
		_t0##.re = _array[addr(_k0<<1)] + _ar; \
		_t1##.re = _array[addr(_k0<<1)] - _ar; \
		_t0##.im = _array[addr(_k0<<1)+1] + _ai; \
		_t1##.im = _array[addr(_k0<<1)+1] - _ai; \
	}

#ifndef USE_ASM
#define cplx_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_pd1 = addr((_k1)<<1);\
		_ar = _array[ _pd1 ];\
		_ai = _array[ _pd1 + 1];\
		_pd0 = addr((_k0)<<1);\
		_t0##.re = _array[ _pd0 ] + _ar; \
		_t1##.re = _array[ _pd0 ] - _ar; \
		_t0##.im = _array[ _pd0 + 1] + _ai; \
		_t1##.im = _array[ _pd0 + 1] - _ai; \
	}
#else
#define cplx_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1) \
	__asm__( "
movl %%ebx,%%eax             # 1  U
andl $-256,%%ebx             # 1  V
shrl $7,%%ebx                # 2  U
addl %%eax,%%eax	     # 2  V
addl %%eax,%%ebx	     # 3  U
fldl (%%esi,%%ebx,8)         # 4  t0r
movl %%edx,%%eax             # 5  U
andl $-256,%%edx             # 5  V
shrl $7,%%edx                # 6  U
addl %%eax,%%eax	     # 6  V
addl %%eax,%%edx	     # 7  U
faddl (%%esi,%%edx,8)        # 8-12   NT0R
fldl  (%%esi,%%ebx,8)        # 9      t0r,NT0R
fsubl (%%esi,%%edx,8)        # 10-14  NT1R,NT0R
fldl  8(%%esi,%%ebx,8)       # 11     t0i,NT1R,NT0R
faddl 8(%%esi,%%edx,8)       # 12-16  NT0I,NT1R,ntor
fxch %%st(2)                 # 12     nt0r,NT1R,NT0I
fldl 8(%%esi,%%ebx,8)        # 13     t0i,nt0r,NT1R,NT0I
fsubl 8(%%esi,%%edx,8)       # 14-18  NT1I,nt0r,NT1R,NTOI
fxch %%st(3)                 # 15     NTOI,nt0r,nt1R,NT1I
"\
:"=&b"(_pd0),"=&d"(_pd1):"S"(_array),"0"(_k0),"1"(_k1):"ax","st");\
               \
               __asm__("
               fstpl  %1                    # 18
               fstpl  %0                    # 19
               fstpl  %2                    # 20
               fstpl  %3                    # 21
               "\
               : :"m"( _t0##.re ),"m"( _t0##.im ),"m"(_t1##.re),"m"(_t1##.im):\
               "memory", "st");
#endif


#define cplx_addsub_store(_d,_k0,_k1,_t0,_t1)\
	_d[addr(_k0<<1)] = _t0##.re + _t1##.re;\
	_d[addr(_k0<<1)+1] = _t0##.im + _t1##.im;\
	_d[addr(_k1<<1)] = _t0##.re - _t1##.re;\
	_d[addr(_k1<<1)+1] = _t0##.im - _t1##.im;


#ifndef USE_ASM
#define cplx_addsub_store_p(_d,_pd0,_pd1,_t0,_t1)\
	_d[ _pd0 ] = _t0##.re + _t1##.re;\
	_d[ _pd0 + 1] = _t0##.im + _t1##.im;\
	_d[ _pd1 ] = _t0##.re - _t1##.re;\
	_d[ _pd1 + 1] = _t0##.im - _t1##.im;
#else

#define cplx_addsub_store_p(_d,_pd0,_pd1,_t0,_t1)\
__asm__( "
fldl %0                    #1   t0r
faddl %2                   #2-5 NT0R
fldl %0                    #3-3 t0r,NT0R
fsubl %2                   #4-7 NT1R,NT0R
fldl %1                    #5-5 t0i,NT1R,NT0R
faddl %3                   #6-9 NT0I,NT1R,ntor
fxch %%st(2)               #6   nt0r,NT1R,NT0I
fldl %1                    #7   t0i,nt0r,NT1R,NT0I
fsubl %3                   #8-11 NT1I,nt0r,NT1R,NTOI
fxch %%st(3)               #9    NTOI,nt0r,nt1R,NT1I
fstpl 8(%%esi,%%ebx,8)     #10   nt0r,NT1R,NT1I
fstpl (%%esi,%%ebx,8)      #11   NT1R,NT1I
fstpl (%%esi,%%edx,8)      #12   nt1I
fstpl 8(%%esi,%%edx,8)     #13   void
"\
: :"m"(_t0##.re),"m"(_t0##.im),"m"(_t1##.re),"m"(_t1##.im),\
"S"(_d),"b"(_pd0),"d"(_pd1):"memory");

#endif



/**********************************************************************
   cplx_mul macro:
		_t= _f0*f1;
   where all the variables are pseudo complex. 
**********************************************************************/
#define cplx_mul(_t,_f0,_f1) \
	{ \
		y_limb_t _r; \
		_r= ((_f0##.re) * (_f1##.re)) - ((_f0##.im) * (_f1##.im)); \
		_t##.im= ((_f0##.re) * (_f1##.im)) + ((_f0##.im) * (_f1##.re)); \
		_t##.re = _r; \
	}

/**********************************************************************
   cplx_squar macro:
		_t= _f0 * _f0;
   where all the variables are pseudo-complex. 
***********************************************************************/
#define cplx_squar(_t,_f0) \
	{ \
		y_limb_t _r; \
		_r= (_f0##.re +_f0##.im ) * (_f0##.re - _f0##.im);\
		_t##.im= 2.0 * (_f0##.re) * (_f0##.im) ;\
		_t##.re = _r; \
	}

/***********************************************************************
   cplx_add macro:
		_t= _s0 + _s1;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_add(_t,_s0,_s1) \
	_t##.re= (_s0##.re) + (_s1##.re); \
	_t##.im= (_s0##.im) + (_s1##.im);

/***********************************************************************
   cplx_sub macro:
		_t= _m - _s;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_sub(_t,_m,_s) \
	_t##.re= (_m##.re) - (_s##.re); \
	_t##.im= (_m##.im) - (_s##.im);


/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#ifndef USE_ASM
# define cplx_mul_1_4_F_addsub(_t0,_t1) \
	{ \
		y_limb_t _ar= _t0##.re, _ai = _t0##.im ,_a = _t1##.re; \
		_t0##.re += (_t1##.im); \
		_t0##.im -= _a; \
		_t1##.re = _ar - (_t1##.im); \
		_t1##.im = _ai + _a; \
	}
#else
# define cplx_mul_1_4_F_addsub(_t0,_t1) \
	__asm__("
fldl %0             #1-1 t0r
faddl %3            #2-5 nt0r
fldl %1             #3-3 t0i,nt0r
fsubl %2            #4-7 nt0i,nt0r
fldl %0             #5-5 t0r,nt0i,nt0r
fsubl %3            #6-9 nt1r,nt0i,nt0r
fldl %1             #7-7 t0i,nt1r,nt0i,nt0r
faddl %2            #8-11 nt1i,nt1r,nt0i,nt0r
fxch %%st(3)        #8-8  nt0r,nt1r,nt0i,nt1i
fstpl %0            #9-9  nt1r,nt0i,nt1i
fxch %%st(1)        #9-9  nt0i,nt1r,nt1i
fstpl %1            #10-10  nt1r,nt1i
fstpl %2            #11-11  nt1i
fstpl %3            #12-12  void
"\
: :"m"( _t0##.re),"m"( _t0##.im), "m"( _t1##.re),"m"(_t1##.im):\
"memory");

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
	_d[addr( _k0<<1 )]= _t0##.re + _t1##.im; \
	_d[addr( _k1<<1 )]= _t0##.re - _t1##.im; \
	_d[addr( _k0<<1 )+ 1]= _t0##.im - _t1##.re; \
	_d[addr( _k1<<1 )+ 1]= _t0##.im + _t1##.re;

#ifndef USE_ASM
#define cplx_mul_1_4_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
	_d[ _pd0 ] = _t0##.re + _t1##.im; \
	_d[ _pd1 ] = _t0##.re - _t1##.im; \
	_d[ _pd0 + 1] = _t0##.im - _t1##.re; \
	_d[ _pd1 + 1] = _t0##.im + _t1##.re;
#else
#define cplx_mul_1_4_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
	__asm__("
fldl %0             	#1-1 t0r
faddl %3            	#2-5 nt0r
fldl %1             	#3-3 t0i,nt0r
fsubl %2            	#4-7 nt0i,nt0r
fldl %0             	#5-5 t0r,nt0i,nt0r
fsubl %3            	#6-9 nt1r,nt0i,nt0r
fldl %1             	#7-7 t0i,nt1r,nt0i,nt0r
faddl %2            	#8-11 nt1i,nt1r,nt0i,nt0r
fxch %%st(3)        	#8-8  nt0r,nt1r,nt0i,nt1i
fstpl (%%esi,%%ebx,8)   #9-9  nt1r,nt0i,nt1i
fxch %%st(1)        	#9-9  nt0i,nt1r,nt1i
fstpl 8(%%esi,%%ebx,8)  #10-10  nt1r,nt1i
fstpl (%%esi,%%edx,8)   #11-11  nt1i
fstpl 8(%%esi,%%edx,8)  #12-12  void
"\
: :"m"( _t0##.re),"m" ( _t0##.im), "m" ( _t1##.re), "m" (_t1##.im),\
"S"(_d),"b"(_pd0),"d"(_pd1):"memory");
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
#define cplx_mul_1_4_B_addsub(_t0,_t1) \
	{ \
		y_limb_t _ar=_t0##.re, _ai=_t0##.im ,_a= _t1##.re; \
		_t0##.re -= (_t1##.im); \
		_t0##.im += _a; \
		_t1##.re = _ar + (_t1##.im); \
		_t1##.im = _ai - _a; \
	}
#else
#define cplx_mul_1_4_B_addsub(_t0,_t1) \
	__asm__("
fldl %0              #1-1 t0r
fsubl %3             #2-5 nt0r
fldl %1              #3-3 t0i,nt0r
faddl %2             #4-7 nt0i,nt0r
fldl %0              #5-5 t0r,nt0i,nt0r
faddl %3             #6-9 nt1r,nt0i,nt0r
fldl %1              #7-7 t0i,nt1r,nt0i,nt0r
fsubl %2             #8-11 nt1i,nt1r,nt0i,nt0r
fxch %%st(3)         #8-8  nt0r,nt1r,nt0i,nt1i
fstpl %0             #9-9  nt1r,nt0i,nt1i
fxch %%st(1)         #9-9  nt0i,nt1r,nt1i
fstpl %1             #10-10  nt1r,nt1i
fstpl %2             #11-11  nt1i
fstpl %3             #12-12  void
"\
: :"m"(_t0##.re),"m"( _t0##.im ),"m"( _t1##.re),"m"(_t1##.im):\
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
	_d[addr( _k0<<1 )]= _t0##.re - _t1##.im; \
	_d[addr( _k1<<1 )]= _t0##.re + _t1##.im; \
	_d[addr( _k0<<1 )+ 1]= _t0##.im + _t1##.re; \
	_d[addr( _k1<<1 )+ 1]= _t0##.im - _t1##.re;

#ifndef USE_ASM
#define cplx_mul_1_4_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
	_d[ _pd0 ] = _t0##.re - _t1##.im; \
	_d[ _pd1 ] = _t0##.re + _t1##.im; \
	_d[ _pd0 + 1] = _t0##.im + _t1##.re; \
	_d[ _pd1 + 1] = _t0##.im - _t1##.re;

#else
#define cplx_mul_1_4_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
	__asm__("
fldl %0              	#1-1 t0r
fsubl %3             	#2-5 nt0r
fldl %1              	#3-3 t0i,nt0r
faddl %2             	#4-7 nt0i,nt0r
fldl %0              	#5-5 t0r,nt0i,nt0r
faddl %3             	#6-9 nt1r,nt0i,nt0r
fldl %1              	#7-7 t0i,nt1r,nt0i,nt0r
fsubl %2             	#8-11 nt1i,nt1r,nt0i,nt0r
fxch %%st(3)         	#8-8  nt0r,nt1r,nt0i,nt1i
fstpl (%%esi,%%ebx,8) 	#9-9  nt1r,nt0i,nt1i
fxch %%st(1)         	#9-9  nt0i,nt1r,nt1i
fstpl 8(%%esi,%%ebx,8)	#10-10  nt1r,nt1i
fstpl (%%esi,%%edx,8)	#11-11  nt1i
fstpl 8(%%esi,%%edx,8)  #12-12  void
"\
: :"m"( _t0##.re),"m"( _t0##.im), "m" (_t1##.re),"m"(_t1##.im),\
"S"(_d),"b"(_pd0),"d"(_pd1): "memory");

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
#define cplx_mul_1_8_F_addsub(_tt0,_tt1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##.re + _tt1##.im )*F_1_8r; \
		_ai = (_tt1##.im - _tt1##.re )*F_1_8r; \
		_tt1##.re = _tt0##.re - _ar;\
		_tt0##.re += _ar;\
		_tt1##.im = _tt0##.im - _ai;\
		_tt0##.im += _ai;\
	}
#else

# define cplx_mul_1_8_F_addsub(_tt0,_tt1) \
        __asm__("
fldl %3                 #1-1 t1i
fld  %%st               #2-2 t1i,t1i
faddl %2                #3-6 t1i+t1r,t1i
fxch %%st(1)            #3-3 t1i,t1i+t1r
fsubl %2                #4-7 t1i-t1r,t1i+t1r
fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r
fmull F_1_8r            #6-10 nt1r,t1i-t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r
fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r
fmull F_1_8r            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl %0                #15-15 nnt1r,nnt0i,nnt1i
fstpl %2                #16-16 nnt0i,nnt1i
fstpl %1                #17-17 nnt1i
fstpl %3                #18-18 void
"\
: :"m" (_tt0##.re ),"m" (_tt0##.im ),"m" (_tt1##.re),"m" (_tt1##.im )\
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
#define cplx_mul_1_8_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re + _t1##.im )*F_1_8r; \
_ai = (_t1##.im - _t1##.re )*F_1_8r; \
_d[addr( _k0<<1 )]= _t0##.re + _ar; \
_d[addr( _k1<<1 )]= _t0##.re - _ar; \
_d[addr( _k0<<1 )+ 1]= _t0##.im + _ai; \
_d[addr( _k1<<1 )+ 1]= _t0##.im - _ai; \
}

#ifndef USE_ASM
#define cplx_mul_1_8_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re + _t1##.im )*F_1_8r; \
_ai = (_t1##.im - _t1##.re )*F_1_8r; \
_d[ _pd0 ]= _t0##.re + _ar; \
_d[ _pd1 ]= _t0##.re - _ar; \
_d[ _pd0 + 1]= _t0##.im + _ai; \
_d[ _pd1 + 1]= _t0##.im - _ai; \
}
#else
#define cplx_mul_1_8_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
__asm__("
fldl %3                 #1-1 t1i
fld  %%st               #2-2 t1i,t1i
faddl %2                #3-6 t1i+t1r,t1i
fxch %%st(1)            #3-3 t1i,t1i+t1r
fsubl %2                #4-7 t1i-t1r,t1i+t1r
fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r
fmull F_1_8r            #6-10 nt1r,t1i-t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r
fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r
fmull F_1_8r            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%edx,8)   #16-16 nnt0i,nnt1i
fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i
fstpl 8(%%esi,%%edx,8)  #18-18 void
"\
: :"m" (_t0##.re ),"m" (_t0##.im ),"m" (_t1##.re),"m" (_t1##.im ),\
"S"(_d),"b"(_pd0),"d"(_pd1): "memory");

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
# define cplx_mul_1_8_B_addsub(_tt0,_tt1) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##.re - _tt1##.im )*F_1_8r; \
		_ai = (_tt1##.im + _tt1##.re )*F_1_8r; \
		_tt1##.re = _tt0##.re - _ar;\
		_tt0##.re += _ar;\
		_tt1##.im = _tt0##.im - _ai;\
		_tt0##.im += _ai;\
	}
#else
#define cplx_mul_1_8_B_addsub(_tt0,_tt1) \
        __asm__("
fldl %2                 #1-1 t1r
fld  %%st               #2-2 t1r,t1r
fsubl %3                #3-6 t1r-t1i,t1i
fxch %%st(1)            #3-3 t1i,t1r-t1i
faddl %3                #4-7 t1i+t1r,t1i-t1r
fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r
fmull F_1_8r            #6-10 nt1r,t1i+t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r
fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r
fmull F_1_8r            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl %0                #15-15 nnt1r,nnt0i,nnt1i
fstpl %2                #16-16 nnt0i,nnt1i
fstpl %1                #17-17 nnt1i
fstpl %3                #18-18 void
"\
: :"m" (_tt0##.re ),"m" (_tt0##.im ),"m" (_tt1##.re),"m" (_tt1##.im )\
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
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re - _t1##.im )*F_1_8r; \
_ai = (_t1##.im + _t1##.re )*F_1_8r; \
_d[addr( _k0<<1 )]= _t0##.re + _ar; \
_d[addr( _k1<<1 )]= _t0##.re - _ar; \
_d[addr( _k0<<1 )+ 1]= _t0##.im + _ai; \
_d[addr( _k1<<1 )+ 1]= _t0##.im - _ai; \
}

#ifndef USE_ASM
#define cplx_mul_1_8_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re - _t1##.im )*F_1_8r; \
_ai = (_t1##.im + _t1##.re )*F_1_8r; \
_d[ _pd0 ] = _t0##.re + _ar; \
_d[ _pd1 ] = _t0##.re - _ar; \
_d[ _pd0 + 1] = _t0##.im + _ai; \
_d[ _pd1 + 1] = _t0##.im - _ai; \
}
#else
#define cplx_mul_1_8_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
__asm__("
fldl %2                 #1-1 t1r
fld  %%st               #2-2 t1r,t1r
fsubl %3                #3-6 t1r-t1i,t1i
fxch %%st(1)            #3-3 t1i,t1r-t1i
faddl %3                #4-7 t1i+t1r,t1i-t1r
fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r
fmull F_1_8r            #6-10 nt1r,t1i+t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r
fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r
fmull F_1_8r            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%edx,8)   #16-16 nnt0i,nnt1i
fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i
fstpl 8(%%esi,%%edx,8)  #18-18 void
"\
: :"m" (_t0##.re ),"m" (_t0##.im ),"m" (_t1##.re),"m" (_t1##.im ),\
"S"(_d),"b"(_pd0),"d"(_pd1):  "memory");
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
# define cplx_mul_3_8_F_addsub(_tt0,_tt1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##.re - _tt1##.im )*F_1_8i; \
		_ai = (_tt1##.im + _tt1##.re )*F_1_8i; \
		_tt1##.re = _tt0##.re - _ar;\
		_tt0##.re += _ar;\
		_tt1##.im = _tt0##.im - _ai;\
		_tt0##.im += _ai;\
	}
#else
#define cplx_mul_3_8_F_addsub(_tt0,_tt1) \
	__asm__("
fldl %2                 #1-1 t1r
fld  %%st               #2-2 t1r,t1r
fsubl %3                #3-6 t1r-t1i,t1i
fxch %%st(1)            #3-3 t1i,t1r-t1i
faddl %3                #4-7 t1i+t1r,t1i-t1r
fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r
fmull F_1_8i            #6-10 nt1r,t1i+t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r
fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r
fmull F_1_8i            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl %0                #15-15 nnt1r,nnt0i,nnt1i
fstpl %2                #16-16 nnt0i,nnt1i
fstpl %1                #17-17 nnt1i
fstpl %3                #18-18 void
"\
: :"m" (_tt0##.re ),"m" (_tt0##.im ),"m" (_tt1##.re),"m" (_tt1##.im )\
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
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re - _t1##.im )*F_1_8i; \
_ai = (_t1##.im + _t1##.re )*F_1_8i; \
_d[addr( _k0<<1 )]= _t0##.re + _ar; \
_d[addr( _k0<<1 )+ 1]= _t0##.im + _ai; \
_d[addr( _k1<<1 )]= _t0##.re - _ar; \
_d[addr( _k1<<1 )+ 1]= _t0##.im - _ai; \
}

#ifndef USE_ASM
#define cplx_mul_3_8_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re - _t1##.im )*F_1_8i; \
_ai = (_t1##.im + _t1##.re )*F_1_8i; \
_d[ _pd0 ] = _t0##.re + _ar; \
_d[ _pd0 + 1] = _t0##.im + _ai; \
_d[ _pd1 ] = _t0##.re - _ar; \
_d[ _pd1 + 1] = _t0##.im - _ai; \
}
#else
#define cplx_mul_3_8_F_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
__asm__("
fldl %2                 #1-1 t1r
fld  %%st               #2-2 t1r,t1r
fsubl %3                #3-6 t1r-t1i,t1i
fxch %%st(1)            #3-3 t1i,t1r-t1i
faddl %3                #4-7 t1i+t1r,t1i-t1r
fldl  %0                #5-5 t0r,t1i+t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i-t1r,t1i+t1r,t0r
fmull F_1_8i            #6-10 nt1r,t1i+t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i+t1r,t0r
fxch %%st(2)            #7-7 t1i+t1r,nt1r,t0i,t0r
fmull F_1_8i            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%edx,8)   #16-16 nnt0i,nnt1i
fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i
fstpl 8(%%esi,%%edx,8)  #18-18 void
"\
: :"m" (_t0##.re ),"m" (_t0##.im ),"m" (_t1##.re),"m" (_t1##.im ),\
"S"(_d),"b"(_pd0),"d"(_pd1):   "memory");
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
#define cplx_mul_3_8_B_addsub(_tt0,_tt1) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##.re + _tt1##.im )*F_1_8i; \
		_ai = (_tt1##.im - _tt1##.re )*F_1_8i; \
		_tt1##.re = _tt0##.re - _ar;\
		_tt0##.re += _ar;\
		_tt1##.im = _tt0##.im - _ai;\
		_tt0##.im += _ai;\
	}
#else
#define cplx_mul_3_8_B_addsub(_tt0,_tt1) \
	__asm__("
fldl %3                 #1-1 t1i
fld  %%st               #2-2 t1i,t1i
faddl %2                #3-6 t1i+t1r,t1i
fxch %%st(1)            #3-3 t1i,t1i+t1r
fsubl %2                #4-7 t1i-t1r,t1i+t1r
fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r
fmull F_1_8i            #6-10 nt1r,t1i-t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r
fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r
fmull F_1_8i            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl %0                #15-15 nnt1r,nnt0i,nnt1i
fstpl %2                #16-16 nnt0i,nnt1i
fstpl %1                #17-17 nnt1i
fstpl %3                #18-18 void
"\
: :"m" (_tt0##.re ),"m" (_tt0##.im ),"m" (_tt1##.re),"m" (_tt1##.im )\
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
#define cplx_mul_3_8_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re + _t1##.im )*F_1_8i; \
_ai = (_t1##.im - _t1##.re )*F_1_8i; \
_d[addr( _k0<<1 )]= _t0##.re + _ar; \
_d[addr( _k0<<1 )+ 1]= _t0##.im + _ai; \
_d[addr( _k1<<1 )]= _t0##.re - _ar; \
_d[addr( _k1<<1 )+ 1]= _t0##.im - _ai; \
}

#ifndef USE_ASM
#define cplx_mul_3_8_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
{ \
y_limb_t _ar,_ai; \
_ar = (_t1##.re + _t1##.im )*F_1_8i; \
_ai = (_t1##.im - _t1##.re )*F_1_8i; \
_d[ _pd0 ] = _t0##.re + _ar; \
_d[ _pd0 + 1] = _t0##.im + _ai; \
_d[ _pd1 ] = _t0##.re - _ar; \
_d[ _pd1 + 1] = _t0##.im - _ai; \
}
#else
#define cplx_mul_3_8_B_addsub_store_p(_d,_pd0,_pd1,_t0,_t1) \
__asm__("
fldl %3                 #1-1 t1i
fld  %%st               #2-2 t1i,t1i
faddl %2                #3-6 t1i+t1r,t1i
fxch %%st(1)            #3-3 t1i,t1i+t1r
fsubl %2                #4-7 t1i-t1r,t1i+t1r
fldl  %0                #5-5 t0r,t1i-t1r,t1i+t1r
fxch %%st(2)            #5-5 t1i+t1r,t1i-t1r,t0r
fmull F_1_8i            #6-10 nt1r,t1i-t1r,t0r
fldl %1                 #7-7 t0i,nt1r,t1i-t1r,t0r
fxch %%st(2)            #7-7 t1i-t1r,nt1r,t0i,t0r
fmull F_1_8i            #8-12 nt1i,nt1r,t0i,t0r
fxch %%st(3)            #8-8 t0r,nt1r,t0i,nt1i
fld  %%st               #9-9 t0r,t0r,nt1r,t0i,nt1i
fadd %%st(2)            #10-13 nnt0r,t0r,nt1r,t0i,nt1i
fxch %%st(1)            #10-10 t0r,nnt0r,nt1r,t0i,nt1i
fsubp %%st,%%st(2)      #11-14 nnt0r,nnt1r,t0i,nt1i
fld %%st(2)             #12-12 t0i,nnt0r,nnt1r,t0i,nt1i
fadd %%st(4)            #13-15 nnt0i,nnt0r,nnt1r,t0i,nt1i
fxch %%st(3)            #13-13 t0i,nnt0r,nnt1r,nnt0i,nt1i
fsubp %%st,%%st(4)      #14-17 nnt0r,nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%ebx,8)   #15-15 nnt1r,nnt0i,nnt1i
fstpl (%%esi,%%edx,8)   #16-16 nnt0i,nnt1i
fstpl 8(%%esi,%%ebx,8)  #17-17 nnt1i
fstpl 8(%%esi,%%edx,8)  #18-18 void
"\
: :"m" (_t0##.re ),"m" (_t0##.im ),"m" (_t1##.re),"m" (_t1##.im ),\
"S"(_d),"b"(_pd0),"d"(_pd1): "memory");
#endif


/***********************************************************************
 This is the dyadic mul for nested complex representation 
 xk =(xk + xmk* ) * (yk + ymk* ) + 2*((xk * yk) - (xmk* * ymk* ) -
     G^(-k) * (xk - xmk* ) * (yk - ymk* );
 xkm =(xkm + xk* ) * (ymk + yk* ) + 2*((xmk * xk*) - (ymk * yk* ) -
     G^(k) * (xmk - xk* ) * (ymk - yk* );
 where all are complex vars. 
***********************************************************************/
/* This is the dyadic mul for nested complex representation */
#define yconv_nested(_xk,_xmk,_yk,_ymk,_i,_j) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i,_gkr,_gki;\
   _gkr=(tw[1].re * Y_DYADIC[ _i << 1]) - (tw[1].im * Y_DYADIC[(_i << 1)+1]) + 1.0;\
   _gki=(tw[1].re * Y_DYADIC[ (_i << 1)+1 ]) + (tw[1].im * Y_DYADIC[_i << 1]);\
   _y0r=(_xk[_i]##.re - _xmk[_j]##.re);\
   _y0i=(_xk[_i]##.im + _xmk[_j]##.im);\
   _y1r=(_yk[_i]##.re - _ymk[_j]##.re);\
   _y1i=(_yk[_i]##.im + _ymk[_j]##.im);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(_xk[_i]##.re * _yk[_i]##.re) - (_xk[_i]##.im * _yk[_i]##.im);\
   _y0i=(_xk[_i]##.re * _yk[_i]##.im) + (_xk[_i]##.im * _yk[_i]##.re);\
   _y1r=(_xmk[_j]##.re * _ymk[_j]##.re) - (_xmk[_j]##.im * _ymk[_j]##.im);\
   _y1i=(_xmk[_j]##.re * _ymk[_j]##.im) + (_xmk[_j]##.im * _ymk[_j]##.re);\
   _y3r= (_gkr * _y2r - _gki * _y2i)*0.25;\
   _y3i= (_gkr * _y2i + _gki * _y2r)*0.25;\
   _xk[_i]##.re = _y0r - _y3r;\
   _xk[_i]##.im = _y0i - _y3i;\
   _xmk[_j]##.re = _y1r - _y3r;\
   _xmk[_j]##.im = _y1i + _y3i;\
}

/* This is the dyadic square for nested-compls representation */
#define ysquare_nested(_xk,_xmk,_i,_j) \
{\
   y_limb_t _r0,_r1,_r2,_r3,_r4,_r5,_r6,_r7,_gkr,_gki;\
   _gkr=(tw[1].re * Y_DYADIC[ _i << 1]) - (tw[1].im * Y_DYADIC[(_i << 1)+1]) +1.0;\
   _gki=(tw[1].re * Y_DYADIC[ (_i << 1)+1 ]) + (tw[1].im * Y_DYADIC[_i << 1]);\
   _r0=(_xk[_i]##.re - _xmk[_j]##.re)*0.5;\
   _r1=(_xk[_i]##.im + _xmk[_j]##.im)*0.5;\
   _r2=(_xk[_i]##.re + _xk[_i]##.im) * (_xk[_i]##.re - _xk[_i]##.im);\
   _r3=(_xk[_i]##.im * _xk[_i]##.re);\
   _r4 =(_r0 + _r1) * (_r0 - _r1);\
   _r5 = _r0 * _r1;\
   _r6 =(_xmk[_j]##.re + _xmk[_j]##.im) * (_xmk[_j]##.re - _xmk[_j]##.im);\
   _r7 = _xmk[_j]##.im * _xmk[_j]##.re;\
   _r3 += _r3;\
   _r5 += _r5;\
   _r7 += _r7;\
   _r0= (_gkr * _r4 - _gki * _r5);\
   _r1= (_gkr * _r5 + _gki * _r4);\
   _xk[_i]##.re = _r2 - _r0;\
   _xmk[_j]##.re = _r6 - _r0;\
   _xk[_i]##.im = _r3 - _r1;\
   _xmk[_j]##.im = _r7 + _r1;\
}

#define ysquare_nested_eq(_xk,_i) \
{\
   y_limb_t _y0r=_xk[_i]##.re, _y0i=_xk[_i]##.im, _gkr;\
   _gkr=(tw[1].re * Y_DYADIC[ _i << 1]) - (tw[1].im * Y_DYADIC[(_i << 1)+1]);\
   _xk[_i]##.re =(_y0r * _y0r + _gkr * (_y0i * _y0i));\
   _xk[_i]##.im = 2.0 * _y0r * _y0i;\
}


/*$Id$*/











