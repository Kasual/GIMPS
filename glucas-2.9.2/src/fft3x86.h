/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
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

#if defined(USE_ASM) && defined(__CYGWIN__)
/* Cygwin's GCC uses the variable exactly as given in asm, but exports with _ prefix in C */
#define MEMORY_VARIABLE_NAME(name) "_" #name
#else
#define MEMORY_VARIABLE_NAME(name) #name
#endif

/* macros to do a FFT tarnsform length 3*/
#ifdef USE_ASM

#define asm_cplx_fft3_F(_r0, _r1, _r2) \
__asm__ volatile ("                                    \n"         \
"fldl %2		#1	r1r                    \n"         \
"faddl %4		#2-4    R1R                    \n"         \
"fldl %2		#3      r1r                    \n"         \
"fsubl %4               #4-6    R2R,R1R                \n"         \
"fldl %3		#5      r1i,R2R,R1R            \n"         \
"faddl %5		#6-8    R1I,R2R,r1r            \n"         \
"fldl %3		#7      r1i,R1I,R2R,r1r        \n"         \
"fsubl %5               #8-10   R2I,R1I,r2r,r1r        \n"         \
"fxch %%st(2)		#8	r2r,R1I,R2I,r1r        \n"         \
"fmull " MEMORY_VARIABLE_NAME(F_1_3i) "         #9-11   R2R,R1I,R2I,r1r        \n"         \
"fldl  %0		#10	r0r,R2R,r1i,R2I,r1r    \n"         \
"fadd %%st(4)		#11-13  R0R,R2R,r1i,r2i,r1r    \n"         \
"fxch %%st(3)		#11     r2i,R2R,r1i,R0R,r1r    \n"         \
"fmull " MEMORY_VARIABLE_NAME(F_1_3i) "         #12-14  R2I,r2r,r1i,R0R,r1r    \n"         \
"fldl %1		#13     r0i,R2I,r2r,r1i,r0r,r1r \n"        \
"fadd %%st(3)		#14-16  R0I,R2I,r2r,r1i,r0r,r1r \n"        \
"fxch %%st(5)		#14 	r1r,R2I,r2r,r1i,r0r,R0I \n"        \
"fmull " MEMORY_VARIABLE_NAME(FM150) "                  #15-17  R1R,r2i,r2r,r1i,r0r,R0I \n"        \
"fxch %%st(3)		#16     r1i,r2i,r2r,R1R,r0r,r0i \n"        \
"fmull " MEMORY_VARIABLE_NAME(FM150) "          #17-19  R1I,r2i,r2r,R1R,r0r,r0i \n"        \
"fxch %%st(4)		#18	r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"fstl %0	       	#19     r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"faddp %%st(0),%%st(3)  #20-22  r2i,r2r,R1R,r1i,r0i     \n"        \
"fxch %%st(4)	        #21     r0i,r2r,R1R,r1i,r2i    \n"         \
"fstl %1		#22     r0i,r2r,R1R,r1i,r2i    \n"         \
"faddp %%st(0),%%st(3)  #23-25  r2r,r1r,R1I,r2i        \n"         \
"fld %%st(1)		#24	r1r,r2r,r1r,R1I,r2i    \n"         \
"fadd %%st(4)		#25-27  NR2R,r2r,r1r,r1i,r2i   \n"         \
"fld %%st(3)		#26     r1i,NR2R,r2r,r1r,r1i,r2i \n"       \
"fsub %%st(2)		#27-29	NR2I,NR2R,r2r,r1r,r1i,r2i \n"      \
"fxch %%st(2)		#28     r2r,NR2R,NR2I,r1r,r1i,r2i \n"      \
"faddp %%st(0),%%st(4)  #29-31	NR2R,NR2I,r1r,NR1I,r2i  \n"        \
"fxch %%st(4)		#29	r2i,NR2I,r1r,NR1I,NR2R  \n"        \
"fsubrp %%st(0),%%st(2) #30-32  nr2i,NR1R,NR1I,nr2r    \n"         \
"fstpl %5		#31     NR1R,NR1I,nr2r         \n"         \
"fxch %%st(2)		#32     nr2r,nr1i,nr1r         \n"         \
"fstpl %4		#33     nr1i,nr1r              \n"         \
"fstpl %3               #34     nr1r                   \n"         \
"fstpl %2		#35     void                   \n"         \
   : :"m"(_r0##r),"m"(_r0##i),"m"(_r1##r),"m"(_r1##i),"m"(_r2##r), \
   "m"(_r2##i) : "memory");

#define asm_cplx_fft3_B(_r0, _r1, _r2) \
__asm__ volatile ("                                    \n"         \
"fldl %2		#1	r1r                    \n"         \
"faddl %4		#2-4    R1R                    \n"         \
"fldl %2		#3      r1r                    \n"         \
"fsubl %4               #4-6    R2R,R1R                \n"         \
"fldl %3		#5      r1i,R2R,R1R            \n"         \
"faddl %5		#6-8    R1I,R2R,r1r            \n"         \
"fldl %3		#7      r1i,R1I,R2R,r1r        \n"         \
"fsubl %5               #8-10   R2I,R1I,r2r,r1r        \n"         \
"fxch %%st(2)		#8	r2r,R1I,R2I,r1r        \n"         \
"fmull " MEMORY_VARIABLE_NAME(B_1_3i) "           #9-11   R2R,R1I,R2I,r1r        \n"         \
"fldl  %0		#10	r0r,R2R,r1i,R2I,r1r    \n"         \
"fadd %%st(4)		#11-13  R0R,R2R,r1i,r2i,r1r    \n"         \
"fxch %%st(3)		#11     r2i,R2R,r1i,R0R,r1r    \n"         \
"fmull " MEMORY_VARIABLE_NAME(B_1_3i) "           #12-14  R2I,r2r,r1i,R0R,r1r    \n"         \
"fldl %1		#13     r0i,R2I,r2r,r1i,r0r,r1r \n"        \
"fadd %%st(3)		#14-16  R0I,R2I,r2r,r1i,r0r,r1r \n"        \
"fxch %%st(5)		#14 	r1r,R2I,r2r,r1i,r0r,R0I \n"        \
"fmull " MEMORY_VARIABLE_NAME(BM150) "            #15-17  R1R,r2i,r2r,r1i,r0r,R0I \n"        \
"fxch %%st(3)		#16     r1i,r2i,r2r,R1R,r0r,r0i \n"        \
"fmull " MEMORY_VARIABLE_NAME(BM150) "            #17-19  R1I,r2i,r2r,R1R,r0r,r0i \n"        \
"fxch %%st(4)		#18	r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"fstl %0		#19     r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"faddp %%st(0),%%st(3)  #20-22  r2i,r2r,R1R,r1i,r0i    \n"         \
"fxch %%st(4)		#21     r0i,r2r,R1R,r1i,r2i    \n"         \
"fstl %1		#22     r0i,r2r,R1R,r1i,r2i    \n"         \
"faddp %%st(0),%%st(3)  #23-25  r2r,r1r,R1I,r2i        \n"         \
"fld %%st(1)	        #24	r1r,r2r,r1r,R1I,r2i    \n"         \
"fadd %%st(4)		#25-27  NR2R,r2r,r1r,r1i,r2i   \n"         \
"fld %%st(3)		#26     r1i,NR2R,r2r,r1r,r1i,r2i \n"       \
"fsub %%st(2)		#27-29	NR2I,NR2R,r2r,r1r,r1i,r2i \n"      \
"fxch %%st(2)		#28     r2r,NR2R,NR2I,r1r,r1i,r2i \n"      \
"faddp %%st(0),%%st(4)  #29-31	NR2R,NR2I,r1r,NR1I,r2i \n"         \
"fxch %%st(4)		#29	r2i,NR2I,r1r,NR1I,NR2R \n"         \
"fsubrp %%st(0),%%st(2) #30-32  nr2i,NR1R,NR1I,nr2r    \n"         \
"fstpl %5		#31     NR1R,NR1I,nr2r         \n"         \
"fxch %%st(2)		#32     nr2r,nr1i,nr1r         \n"         \
"fstpl %4		#33     nr1i,nr1r              \n"         \
"fstpl %3               #34     nr1r                   \n"         \
"fstpl %2		#35     void                   \n"         \
   : :"m"(_r0##r),"m"(_r0##i),"m"(_r1##r),"m"(_r1##i),"m"(_r2##r), \
   "m"(_r2##i) : "memory");

#define asm_cplx_fft3_store_p_F(_d,_pd0,_pd1,_pd2,_r0, _r1, _r2)   \
__asm__ volatile ("                                    \n"         \
"fldl %2		#1	r1r                    \n"         \
"faddl %4		#2-4    R1R                    \n"         \
"fldl %2		#3      r1r                    \n"         \
"fsubl %4               #4-6    R2R,R1R                \n"         \
"fldl %3		#5      r1i,R2R,R1R            \n"         \
"faddl %5	        #6-8    R1I,R2R,r1r            \n"         \
"fldl %3		#7      r1i,R1I,R2R,r1r        \n"         \
"fsubl %5               #8-10   R2I,R1I,r2r,r1r        \n"         \
"fxch %%st(2)		#8	r2r,R1I,R2I,r1r        \n"         \
"fmull " MEMORY_VARIABLE_NAME(F_1_3i) "         #9-11   R2R,R1I,R2I,r1r        \n"         \
"fldl  %0		#10	r0r,R2R,r1i,R2I,r1r    \n"         \
"fadd %%st(4)		#11-13  R0R,R2R,r1i,r2i,r1r    \n"         \
"fxch %%st(3)		#11     r2i,R2R,r1i,R0R,r1r    \n"         \
"fmull " MEMORY_VARIABLE_NAME(F_1_3i) "         #12-14  R2I,r2r,r1i,R0R,r1r    \n"         \
"fldl %1		#13     r0i,R2I,r2r,r1i,r0r,r1r \n"        \
"fadd %%st(3)		#14-16  R0I,R2I,r2r,r1i,r0r,r1r \n"        \
"fxch %%st(5)		#14 	r1r,R2I,r2r,r1i,r0r,R0I \n"        \
"fmull " MEMORY_VARIABLE_NAME(FM150) "                  #15-17  R1R,r2i,r2r,r1i,r0r,R0I \n"        \
"fxch %%st(3)		#16     r1i,r2i,r2r,R1R,r0r,r0i \n"        \
"fmull " MEMORY_VARIABLE_NAME(FM150) "          #17-19  R1I,r2i,r2r,R1R,r0r,r0i \n"        \
"fxch %%st(4)		#18	r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"fstl (%%eax,%%edi,8)	#19     r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"faddp %%st(0),%%st(3)  #20-22  r2i,r2r,R1R,r1i,r0i    \n"         \
"fxch %%st(4)		#21     r0i,r2r,R1R,r1i,r2i    \n"         \
"fstl 8(%%eax,%%edi,8)	#22     r0i,r2r,R1R,r1i,r2i    \n"         \
"faddp %%st(0),%%st(3)  #23-25  r2r,r1r,R1I,r2i        \n"         \
"fld %%st(1)		#24	r1r,r2r,r1r,R1I,r2i    \n"         \
"fadd %%st(4)		#25-27  NR2R,r2r,r1r,r1i,r2i   \n"         \
"fld %%st(3)		#26     r1i,NR2R,r2r,r1r,r1i,r2i \n"       \
"fsub %%st(2)		#27-29	NR2I,NR2R,r2r,r1r,r1i,r2i \n"      \
"fxch %%st(2)		#28     r2r,NR2R,NR2I,r1r,r1i,r2i \n"      \
"faddp %%st(0),%%st(4)  #29-31	NR2R,NR2I,r1r,NR1I,r2i \n"         \
"fxch %%st(4)		#29	r2i,NR2I,r1r,NR1I,NR2R \n"         \
"fsubrp %%st(0),%%st(2) #30-32  nr2i,NR1R,NR1I,nr2r    \n"         \
"fstpl 8(%%edx,%%edi,8)	#31     NR1R,NR1I,nr2r         \n"         \
"fxch %%st(2)		#32     nr2r,nr1i,nr1r         \n"         \
"fstpl (%%edx,%%edi,8)	#33     nr1i,nr1r              \n"         \
"fstpl 8(%%ecx,%%edi,8) #34     nr1r                   \n"         \
"fstpl (%%ecx,%%edi,8)	#35     void                   \n"         \
   : :"m"(_r0##r),"m"(_r0##i),"m"(_r1##r),"m"(_r1##i),"m"(_r2##r), \
   "m"(_r2##i),"D"(_d),"a"(_pd0),"c"(_pd1),"d"(_pd2): "memory");

#define asm_cplx_fft3_store_p_B(_d,_pd0,_pd1,_pd2,_r0, _r1, _r2) \
__asm__ volatile ("                                    \n"         \
"fldl %2		#1	r1r                    \n"         \
"faddl %4		#2-4    R1R                    \n"         \
"fldl %2		#3      r1r                    \n"         \
"fsubl %4               #4-6    R2R,R1R                \n"         \
"fldl %3		#5      r1i,R2R,R1R            \n"         \
"faddl %5		#6-8    R1I,R2R,r1r            \n"         \
"fldl %3		#7      r1i,R1I,R2R,r1r        \n"         \
"fsubl %5               #8-10   R2I,R1I,r2r,r1r        \n"         \
"fxch %%st(2)		#8	r2r,R1I,R2I,r1r        \n"         \
"fmull " MEMORY_VARIABLE_NAME(B_1_3i) "           #9-11   R2R,R1I,R2I,r1r        \n"         \
"fldl  %0		#10	r0r,R2R,r1i,R2I,r1r    \n"         \
"fadd %%st(4)		#11-13  R0R,R2R,r1i,r2i,r1r    \n"         \
"fxch %%st(3)		#11     r2i,R2R,r1i,R0R,r1r    \n"         \
"fmull " MEMORY_VARIABLE_NAME(B_1_3i) "           #12-14  R2I,r2r,r1i,R0R,r1r    \n"         \
"fldl %1		#13     r0i,R2I,r2r,r1i,r0r,r1r \n"        \
"fadd %%st(3)		#14-16  R0I,R2I,r2r,r1i,r0r,r1r \n"        \
"fxch %%st(5)		#14 	r1r,R2I,r2r,r1i,r0r,R0I \n"        \
"fmull " MEMORY_VARIABLE_NAME(BM150) "            #15-17  R1R,r2i,r2r,r1i,r0r,R0I \n"        \
"fxch %%st(3)		#16     r1i,r2i,r2r,R1R,r0r,r0i \n"        \
"fmull " MEMORY_VARIABLE_NAME(BM150) "            #17-19  R1I,r2i,r2r,R1R,r0r,r0i \n"        \
"fxch %%st(4)		#18	r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"fstl (%%eax,%%edi,8)	#19     r0r,r2i,r2r,r1r,R1I,r0i \n"        \
"faddp %%st(0),%%st(3)  #20-22  r2i,r2r,R1R,r1i,r0i    \n"         \
"fxch %%st(4)		#21     r0i,r2r,R1R,r1i,r2i    \n"         \
"fstl 8(%%eax,%%edi,8)	#22     r0i,r2r,R1R,r1i,r2i    \n"         \
"faddp %%st(0),%%st(3)  #23-25  r2r,r1r,R1I,r2i        \n"         \
"fld %%st(1)		#24	r1r,r2r,r1r,R1I,r2i    \n"         \
"fadd %%st(4)		#25-27  NR2R,r2r,r1r,r1i,r2i   \n"         \
"fld %%st(3)		#26     r1i,NR2R,r2r,r1r,r1i,r2i \n"       \
"fsub %%st(2)		#27-29	NR2I,NR2R,r2r,r1r,r1i,r2i \n"      \
"fxch %%st(2)		#28     r2r,NR2R,NR2I,r1r,r1i,r2i \n"      \
"faddp %%st(0),%%st(4)  #29-31	NR2R,NR2I,r1r,NR1I,r2i  \n"        \
"fxch %%st(4)		#29	r2i,NR2I,r1r,NR1I,NR2R \n"         \
"fsubrp %%st(0),%%st(2) #30-32  nr2i,NR1R,NR1I,nr2r    \n"         \
"fstpl 8(%%edx,%%edi,8)	#31     NR1R,NR1I,nr2r         \n"         \
"fxch %%st(2)		#32     nr2r,nr1i,nr1r         \n"         \
"fstpl (%%edx,%%edi,8)	#33     nr1i,nr1r              \n"         \
"fstpl 8(%%ecx,%%edi,8) #34     nr1r                   \n"         \
"fstpl (%%ecx,%%edi,8)	#35     void                   \n"         \
   : :"m"(_r0##r),"m"(_r0##i),"m"(_r1##r),"m"(_r1##i),"m"(_r2##r), \
   "m"(_r2##i),"D"(_d),"a"(_pd0),"c"(_pd1),"d"(_pd2): "memory");

#endif


/* These are the Nussbaumer version. Taken from E. Mayer (Mlucas code) */

#ifndef USE_ASM

#define cplx_fft3(_S, _r0, _r1, _r2) \
	{ \
	   y_limb_t _ar,_ai;\
	   cplx_addsub( _r1, _r2);\
	   _r0##r += _r1##r;\
	   _r0##i += _r1##i;\
	   _r2##r *= _S##_1_3i;\
	   _r2##i *= _S##_1_3i;\
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   _ar = _r2##r;\
	   _ai = _r2##i;\
	   _r2##r = _r1##r + _ai;\
	   _r2##i = _r1##i - _ar;\
           _r1##r = _r1##r - _ai;\
           _r1##i = _r1##i + _ar;\
        }

#else

#define cplx_fft3(_S, _r0, _r1, _r2) asm_cplx_fft3_##_S( _r0,_r1,_r2)

#endif

#ifndef USE_ASM

#define cplx_fft3_store_p( _index, _pd0, _pd1, _pd2, _S, _r0, _r1, _r2) \
	{ \
	   cplx_addsub( _r1, _r2);\
	   _r0##r += _r1##r;\
	   _r0##i += _r1##i;\
	   _r2##r *= _S##_1_3i;\
	   _r2##i *= _S##_1_3i;\
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   *( _pd0 + _index) = _r0##r;\
	   *( _pd0 + _index + 1) = _r0##i;\
	   *( _pd2 + _index) = _r1##r + _r2##i;\
	   *( _pd2 + _index + 1) = _r1##i - _r2##r;\
	   *( _pd1 + _index) = _r1##r - _r2##i;\
	   *( _pd1 + _index+ 1) = _r1##i + _r2##r;\
	}

#else

#define cplx_fft3_store_p( _d, _pd0, _pd1, _pd2, _S, _r0, _r1, _r2) \
	asm_cplx_fft3_store_p_##_S(_d,_pd0,_pd1,_pd2,_r0,_r1,_r2);

#endif

/*$Id$*/

