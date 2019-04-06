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

/* Nussbaumer version for FFT-5. Code adapted from E. Mayer's (Mlucas)     */

#ifdef USE_ASM

#define asm_cplx_fft5_F(_x0, _x1, _x2, _x3, _x4)                       \
__asm__ volatile("                                                \n"  \
"fldl %2                  #1      x1r                             \n"  \
"faddl %8                 #2-5    NX1R                            \n"  \
"fldl %2                  #3      x1r,NX1R                        \n"  \
"fsubl %8                 #4-7    NX4R,NX1R                       \n"  \
"fldl %4                  #5      x2r,NX4R,NX1R                   \n"  \
"fsubl %6                 #6-9    NX3R,NX4R,nx1r                  \n"  \
"fldl %4                  #7      x2r,NX3R,NX4R,nx1r              \n"  \
"faddl %6                 #8-11   NX2R,NX3R,n4xr,nx1r             \n"  \
"fldl %3                  #9      x1i,NX2R,NX3R,n4xr,nx1r         \n"  \
"fsubl %9                 #10-13  NX4I,NX2R,nx3r,n4xr,nx1r        \n"  \
"fxch %%st(3)             #10     n4xr,nX2R,nx3r,NX4I,nx1r        \n"  \
"fstpl %8                 #11     nx2r,nx3r,NX4I,nx1r             \n"  \
"fxch %%st(2)             #12     NX4I,nx3r,nx2r,nx1r             \n"  \
"fldl %3                  #13     x1i,NX4I,nx3r,nx2r,nx1r         \n"  \
"faddl %9                 #14-17  NX1I,nx4i,nx3r,nx2r,nx1r        \n"  \
"fxch %%st(2)             #14     nx3r,nx4i,NX1I,nx2r,nx1r        \n"  \
"fldl %5                  #15     x2i,nx3r,nx4i,NX1I,nx2r,nx1r    \n"  \
"fsubl %7                 #16-19  NX3I,nx3r,nx4i,NX1I,nx2r,nx1r   \n"  \
"fxch %%st(2)             #16     nx4i,nx3r,NX3I,nx1i,nx2r,nx1r   \n"  \
"fldl %5                  #17     x2i,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n" \
"faddl %7                 #18-21  NX2I,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n"\
"fxch %%st(2)             #19     nx3r,nx4i,NX2I,NX3I,nx1i,nx2r,nx1r\n"\
"fstpl %6                 #20     nx4i,NX2I,NX3I,nx1i,nx2r,nx1r   \n"  \
"fstpl %9                 #21     NX2I,NX3I,nx1i,nx2r,nx1r        \n"  \
"fxch %%st(1)             #22     nx3i,nx2i,nx1i,nx2r,nx1r        \n"  \
"fstpl %7                 #23     x2i,x1i,x2r,x1r                 \n"  \
"fld %%st(3)              #24     x1r,x2i,x1i,x2r,x1r             \n"  \
"fadd %%st(3)             #25-28  AR,x2i,x1i,x2r,x1r              \n"  \
"fld %%st(2)              #26     x1i,AR,x2i,x1i,x2r,x1r          \n"  \
"fadd %%st(2)             #27-30  AI,AR,x2i,x1i,x2r,x1r           \n"  \
"fldl %0                  #28     x0r,AI,AR,x2i,x1i,x2r,x1r       \n"  \
"fadd %%st(2)             #29-32  NX0R,AI,AR,x2i,x1i,x2r,x1r      \n"  \
"fxch %%st(6)             #29     x1r,AI,ar,x2i,x1i,x2r,NX0R      \n"  \
"fsubp %%st,%%st(5)       #30-33  AI,ar,x2i,x1i,NX2R,NX0R         \n"  \
"fxch %%st(3)             #30     x1i,ar,x2i,ai,NX2R,NX0R         \n"  \
"fsubp %%st,%%st(2)       #31-34  ar,NX2I,ai,NX2R,NX0R            \n"  \
"fldl %1                  #32     x0i,ar,NX2I,ai,NX2R,nx0r        \n"  \
"fadd %%st(3)             #33     NX0I,ar,NX2I,ai,nx2r,nx0r       \n"  \
"fxch %%st(4)             #34     nx2r,ar,NX2I,ai,NX0i,nx0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5r) "             #35-38  X2R,ar,nx2i,ai,NX0I,nx0r        \n"  \
"fxch %%st(5)             #36     nx0r,ar,nx2i,ai,NX0I,X2R        \n"  \
"fstpl %0                 #37     ar,nx2i,ai,nx0i,X2R             \n"  \
"fxch %%st(1)             #37     nx2i,ar,ai,nx0i,X2R             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5r) "             #38-41  X2I,ar,ai,nx0i,X2R              \n"  \
"fxch %%st(3)             #39     nx0i,ar,ai,X2I,x2r              \n"  \
"fstpl %1                 #40     ar,ai,X2I,x2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(FNM125) "             #41-44  AR,ai,X2I,x2r                   \n"  \
"fxch %%st(1)             #42     ai,AR,X2i,X2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(FNM125) "             #43-45  AI,AR,x2i,x2r                   \n"  \
"fxch %%st(1)             #43     AR,AI,x2i,x2r                   \n"  \
"faddl %0                 #44-47  NAR,AI,x2i,x2r                  \n"  \
"fxch %%st(1)             #45     AI,NAR,x2i,x2r                  \n"  \
"faddl %1                 #46     NAI,NAR,x2i,x2r                 \n"  \
"fld %%st(1)              #47     nar,NAI,nar,x2i,x2r             \n"  \
"fadd %%st(4)             #48     X1R,NAI,nar,x2i,x2r             \n"  \
"fxch %%st(2)             #46     nar,NAI,X1R,x2i,x2r             \n"  \
"fsubp %%st,%%st(4)       #47-50  nai,X1R,x2i,NX2R                \n"  \
"fld %%st                 #48     nai,nai,X1R,x2i,NX2R,           \n"  \
"fadd %%st(3)             #49-52  X1I,nai,X1R,x2i,NX2R            \n"  \
"fxch %%st(3)             #49     x2i,nai,x1r,X1I,NX2R            \n"  \
"fsubrp %%st,%%st(1)      #50     X2I,x1r,X1I,x2r                 \n"  \
"fxch %%st(3)             #51     x2r,x1r,x1I,X2I                 \n"  \
"fstpl %4                 #52     x1r,x1i,x2i                     \n"  \
"fstpl %2                 #53     x1i,x2i                         \n"  \
"fstpl %3                 #54     x2i                             \n"  \
"fstpl %5                 #55                                     \n"  \
"fldl %8		  #56     x4r                             \n"  \
"fsubl %6		  #57-60  AR                              \n"  \
"fldl %9 		  #58     x4i,AR                          \n"  \
"fsubl %7                 #59-63  AI,AR                           \n"  \
"fldl %6                  #60     x3r,AI,AR                       \n"  \
"fxch %%st(2)             #60	  ar,AI,x3r                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_5i) "             #61-64  NAR,Ai,x3r                      \n"  \
"fldl %7                  #62     x3i,NAR,Ai,x3r                  \n"  \
"fxch %%st(2)             #63     ai,NAR,x3i,x3r                  \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_5i) "             #64-67  NAI,NAR,x3i,x3r                 \n"  \
"fldl %8		  #65     x4r,NAI,NAR,x3i,x3r             \n"  \
"fxch %%st(4)             #65     x3r,NAI,nar,x3i,x4r             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_5i) "             #66-69  NX3R,NAI,nar,x3i,x4r            \n"  \
"fldl %9                  #67     x4i,NX3R,NAI,nar,x3i,x4r        \n"  \
"fxch %%st(4)             #67     x3i,NX3R,nai,nar,x4i,x4r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_5i) "             #68-71  NX3I,NX3R,nai,nar,x4i,x4r       \n"  \
"fxch %%st(1)             #68     NX3R,NX3I,nai,nar,x4i,x4r       \n"  \
"fadd %%st(3)             #69-72  X3R,NX3I,nai,nar,x4i,x4r        \n"  \
"fxch %%st(5)             #69     x4r,NX3I,nai,nar,x4i,X3R        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5i) "             #70-73  NX4R,NX3I,nai,nar,x4i,X3R       \n"  \
"fxch %%st(1)             #70     nx3i,NX4R,nai,nar,x4i,X3r       \n"  \
"fadd %%st(2)             #71-74  X3I,NX4R,nai,nar,x4i,x3r        \n"  \
"fxch %%st(4)             #71     x4i,NX4r,nai,nar,X3I,x3r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5i) "             #72-75  NX4I,NX4R,nai,nar,X3I,x3r       \n"  \
"fxch %%st(1)             #72     NX4r,NX4I,nai,nar,X3I,x3r       \n"  \
"fsubrp %%st,%%st(3)      #73-76  NX4I,nai,AR,x3i,x3r             \n"  \
"fldl %2                  #74     x1r,NX4I,nai,AR,x3i,x3r         \n"  \
"fldl %3                  #75     x1i,x1r,NX4I,nai,AR,x3i,X3r     \n"  \
"fxch %%st(2)             #75     nx4i,x1r,x1i,nai,AR,x3i,x3r     \n"  \
"fsubrp %%st,%%st(3)      #76-79  x1r,x1i,AI,AR,x3i,x3r           \n"  \
"fadd %%st(4)		  #77-80  X4R,x1i,AI,AR,x3i,x3r           \n"  \
"fxch %%st(1)             #77     x1i,X4R,AI,AR,x3i,x3r           \n"  \
"fsub %%st(5)             #78-81  X4I,X4R,Ai,ar,x3i,x3r           \n"  \
"fldl %2                  #79     x1r,X4I,X4R,Ai,ar,x3i,x3r       \n"  \
"fsubp %%st,%%st(5)       #80-83  x4i,x4r,ai,ar,X1R,x3r           \n"  \
"fldl %3                  #81     x1i,x4i,x4r,ai,ar,X1R,x3r       \n"  \
"faddp %%st,%%st(6)       #82     x4i,x4r,ai,ar,X1R,X1I           \n"  \
"fstpl %9                 #83     x4r,ai,ar,X1R,X1I               \n"  \
"fstpl %8                 #84     ai,ar,X1R,X1I                   \n"  \
"fldl %4                  #85     x2r,ai,ar,x1r,x1i               \n"  \
"fadd %%st(1)             #86-89  X3R,ai,ar,x1r,x1i               \n"  \
"fxch %%st(3)             #86     x1r,ai,ar,X3R,x1i               \n"  \
"fldl %5                  #87     x2i,x1r,ai,ar,X3R,x1i           \n"  \
"fsub %%st(3)             #88-91  X3I,x1r,ai,ar,X3R,x1i           \n"  \
"fxch %%st(5)             #89     x1i,x1r,ai,ar,X3R,X3I           \n"  \
"fldl %4                  #90     x2r,x1i,x1r,ai,ar,x3r,x3i       \n"  \
"fsubp %%st,%%st(3)       #91-94  x1i,x1r,X2R,ar,x3r,x3i          \n"  \
"fldl %5                  #92     x2i,x1i,x1r,X2R,ar,x3r,x3i      \n"  \
"faddp %%st,%%st(4)       #93-96  x1i,x1r,X2R,X2I,x3r,x3i         \n"  \
"fstpl %3                 #94     x1r,x2r,X2I,x3r,x3i             \n"  \
"fstpl %2                 #95     x2r,X2I,x3r,x3i                 \n"  \
"fstpl %4                 #96     X2I,x3r,x3i                     \n"  \
"fstpl %5                 #97     x3r,x3i                         \n"  \
"fstpl %6                 #98     x3i                             \n"  \
"fstpl %7                 #99                                     \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),     \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory");

#define asm_cplx_fft5_B(_x0, _x1, _x2, _x3, _x4)                       \
__asm__ volatile("                                                \n"  \
"fldl %2                  #1      x1r                             \n"  \
"faddl %8                 #2-5    NX1R                            \n"  \
"fldl %2                  #3      x1r,NX1R                        \n"  \
"fsubl %8                 #4-7    NX4R,NX1R                       \n"  \
"fldl %4                  #5      x2r,NX4R,NX1R                   \n"  \
"fsubl %6                 #6-9    NX3R,NX4R,nx1r                  \n"  \
"fldl %4                  #7      x2r,NX3R,NX4R,nx1r              \n"  \
"faddl %6                 #8-11   NX2R,NX3R,n4xr,nx1r             \n"  \
"fldl %3                  #9      x1i,NX2R,NX3R,n4xr,nx1r         \n"  \
"fsubl %9                 #10-13  NX4I,NX2R,nx3r,n4xr,nx1r        \n"  \
"fxch %%st(3)             #10     n4xr,nX2R,nx3r,NX4I,nx1r        \n"  \
"fstpl %8                 #11     nx2r,nx3r,NX4I,nx1r             \n"  \
"fxch %%st(2)             #12     NX4I,nx3r,nx2r,nx1r             \n"  \
"fldl %3                  #13     x1i,NX4I,nx3r,nx2r,nx1r         \n"  \
"faddl %9                 #14-17  NX1I,nx4i,nx3r,nx2r,nx1r        \n"  \
"fxch %%st(2)             #14     nx3r,nx4i,NX1I,nx2r,nx1r        \n"  \
"fldl %5                  #15     x2i,nx3r,nx4i,NX1I,nx2r,nx1r    \n"  \
"fsubl %7                 #16-19  NX3I,nx3r,nx4i,NX1I,nx2r,nx1r   \n"  \
"fxch %%st(2)             #16     nx4i,nx3r,NX3I,nx1i,nx2r,nx1r   \n"  \
"fldl %5                  #17     x2i,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r \n"\
"faddl %7                 #18-21  NX2I,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n"\
"fxch %%st(2)             #19     nx3r,nx4i,NX2I,NX3I,nx1i,nx2r,nx1r\n"\
"fstpl %6                 #20     nx4i,NX2I,NX3I,nx1i,nx2r,nx1r   \n"  \
"fstpl %9                 #21     NX2I,NX3I,nx1i,nx2r,nx1r        \n"  \
"fxch %%st(1)             #22     nx3i,nx2i,nx1i,nx2r,nx1r        \n"  \
"fstpl %7                 #23     x2i,x1i,x2r,x1r                 \n"  \
"fld %%st(3)              #24     x1r,x2i,x1i,x2r,x1r             \n"  \
"fadd %%st(3)             #25-28  AR,x2i,x1i,x2r,x1r              \n"  \
"fld %%st(2)              #26     x1i,AR,x2i,x1i,x2r,x1r          \n"  \
"fadd %%st(2)             #27-30  AI,AR,x2i,x1i,x2r,x1r           \n"  \
"fldl %0                  #28     x0r,AI,AR,x2i,x1i,x2r,x1r       \n"  \
"fadd %%st(2)             #29-32  NX0R,AI,AR,x2i,x1i,x2r,x1r      \n"  \
"fxch %%st(6)             #29     x1r,AI,ar,x2i,x1i,x2r,NX0R      \n"  \
"fsubp %%st,%%st(5)       #30-33  AI,ar,x2i,x1i,NX2R,NX0R         \n"  \
"fxch %%st(3)             #30     x1i,ar,x2i,ai,NX2R,NX0R         \n"  \
"fsubp %%st,%%st(2)       #31-34  ar,NX2I,ai,NX2R,NX0R            \n"  \
"fldl %1                  #32     x0i,ar,NX2I,ai,NX2R,nx0r        \n"  \
"fadd %%st(3)             #33     NX0I,ar,NX2I,ai,nx2r,nx0r       \n"  \
"fxch %%st(4)             #34     nx2r,ar,NX2I,ai,NX0i,nx0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5r) "             #35-38  X2R,ar,nx2i,ai,NX0I,nx0r        \n"  \
"fxch %%st(5)             #36     nx0r,ar,nx2i,ai,NX0I,X2R        \n"  \
"fstpl %0                 #37     ar,nx2i,ai,nx0i,X2R             \n"  \
"fxch %%st(1)             #37     nx2i,ar,ai,nx0i,X2R             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5r) "             #38-41  X2I,ar,ai,nx0i,X2R              \n"  \
"fxch %%st(3)             #39     nx0i,ar,ai,X2I,x2r              \n"  \
"fstpl %1                 #40     ar,ai,X2I,x2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(BNM125) "             #41-44  AR,ai,X2I,x2r                   \n"  \
"fxch %%st(1)             #42     ai,AR,X2i,X2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(BNM125) "             #43-45  AI,AR,x2i,x2r                   \n"  \
"fxch %%st(1)             #43     AR,AI,x2i,x2r                   \n"  \
"faddl %0                 #44-47  NAR,AI,x2i,x2r                  \n"  \
"fxch %%st(1)             #45     AI,NAR,x2i,x2r                  \n"  \
"faddl %1                 #46     NAI,NAR,x2i,x2r                 \n"  \
"fld %%st(1)              #47     nar,NAI,nar,x2i,x2r             \n"  \
"fadd %%st(4)             #48     X1R,NAI,nar,x2i,x2r             \n"  \
"fxch %%st(2)             #46     nar,NAI,X1R,x2i,x2r             \n"  \
"fsubp %%st,%%st(4)       #47-50  nai,X1R,x2i,NX2R                \n"  \
"fld  %%st                #48     nai,nai,X1R,x2i,NX2R,           \n"  \
"fadd %%st(3)             #49-52  X1I,nai,X1R,x2i,NX2R            \n"  \
"fxch %%st(3)             #49     x2i,nai,x1r,X1I,NX2R            \n"  \
"fsubrp %%st,%%st(1)      #50     X2I,x1r,X1I,x2r                 \n"  \
"fxch %%st(3)             #51     x2r,x1r,x1I,X2I                 \n"  \
"fstpl %4                 #52     x1r,x1i,x2i                     \n"  \
"fstpl %2                 #53     x1i,x2i                         \n"  \
"fstpl %3                 #54     x2i                             \n"  \
"fstpl %5                 #55                                     \n"  \
"fldl %8		  #56     x4r                             \n"  \
"fsubl %6		  #57-60 AR                              \n"  \
"fldl %9 		  #58     x4i,AR                          \n"  \
"fsubl %7                 #59-63  AI,AR                           \n"  \
"fldl %6                  #60     x3r,AI,AR                       \n"  \
"fxch %%st(2)             #60	  ar,AI,x3r                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(B_1_5i) "             #61-64  NAR,Ai,x3r                      \n"  \
"fldl %7                  #62     x3i,NAR,Ai,x3r                  \n"  \
"fxch %%st(2)             #63     ai,NAR,x3i,x3r                  \n"  \
"fmull " MEMORY_VARIABLE_NAME(B_1_5i) "             #64-67  NAI,NAR,x3i,x3r                 \n"  \
"fldl %8		  #65     x4r,NAI,NAR,x3i,x3r             \n"  \
"fxch %%st(4)             #65     x3r,NAI,nar,x3i,x4r             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_5i) "             #66-69  NX3R,NAI,nar,x3i,x4r            \n"  \
"fldl %9                  #67     x4i,NX3R,NAI,nar,x3i,x4r        \n"  \
"fxch %%st(4)             #67     x3i,NX3R,nai,nar,x4i,x4r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_5i) "             #68-71  NX3I,NX3R,nai,nar,x4i,x4r       \n"  \
"fxch %%st(1)             #68     NX3R,NX3I,nai,nar,x4i,x4r       \n"  \
"fadd %%st(3)             #69-72  X3R,NX3I,nai,nar,x4i,x4r        \n"  \
"fxch %%st(5)             #69     x4r,NX3I,nai,nar,x4i,X3R        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5i) "             #70-73  NX4R,NX3I,nai,nar,x4i,X3R       \n"  \
"fxch %%st(1)             #70     nx3i,NX4R,nai,nar,x4i,X3r       \n"  \
"fadd %%st(2)             #71-74  X3I,NX4R,nai,nar,x4i,x3r        \n"  \
"fxch %%st(4)             #71     x4i,NX4r,nai,nar,X3I,x3r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5i) "             #72-75  NX4I,NX4R,nai,nar,X3I,x3r       \n"  \
"fxch %%st(1)             #72     NX4r,NX4I,nai,nar,X3I,x3r       \n"  \
"fsubrp %%st,%%st(3)      #73-76  NX4I,nai,AR,x3i,x3r             \n"  \
"fldl %2                  #74     x1r,NX4I,nai,AR,x3i,x3r         \n"  \
"fldl %3                  #75     x1i,x1r,NX4I,nai,AR,x3i,X3r     \n"  \
"fxch %%st(2)             #75     nx4i,x1r,x1i,nai,AR,x3i,x3r     \n"  \
"fsubrp %%st,%%st(3)      #76-79  x1r,x1i,AI,AR,x3i,x3r           \n"  \
"fadd %%st(4)		  #77-80  X4R,x1i,AI,AR,x3i,x3r           \n"  \
"fxch %%st(1)             #77     x1i,X4R,AI,AR,x3i,x3r           \n"  \
"fsub %%st(5)             #78-81  X4I,X4R,Ai,ar,x3i,x3r           \n"  \
"fldl %2                  #79     x1r,X4I,X4R,Ai,ar,x3i,x3r       \n"  \
"fsubp %%st,%%st(5)       #80-83  x4i,x4r,ai,ar,X1R,x3r           \n"  \
"fldl %3                  #81     x1i,x4i,x4r,ai,ar,X1R,x3r       \n"  \
"faddp %%st,%%st(6)       #82     x4i,x4r,ai,ar,X1R,X1I           \n"  \
"fstpl %9                 #83     x4r,ai,ar,X1R,X1I               \n"  \
"fstpl %8                 #84     ai,ar,X1R,X1I                   \n"  \
"fldl %4                  #85     x2r,ai,ar,x1r,x1i               \n"  \
"fadd %%st(1)             #86-89  X3R,ai,ar,x1r,x1i               \n"  \
"fxch %%st(3)             #86     x1r,ai,ar,X3R,x1i               \n"  \
"fldl %5                  #87     x2i,x1r,ai,ar,X3R,x1i           \n"  \
"fsub %%st(3)             #88-91  X3I,x1r,ai,ar,X3R,x1i           \n"  \
"fxch %%st(5)             #89     x1i,x1r,ai,ar,X3R,X3I           \n"  \
"fldl %4                  #90     x2r,x1i,x1r,ai,ar,x3r,x3i       \n"  \
"fsubp %%st,%%st(3)       #91-94  x1i,x1r,X2R,ar,x3r,x3i          \n"  \
"fldl %5                  #92     x2i,x1i,x1r,X2R,ar,x3r,x3i      \n"  \
"faddp %%st,%%st(4)       #93-96  x1i,x1r,X2R,X2I,x3r,x3i         \n"  \
"fstpl %3                 #94     x1r,x2r,X2I,x3r,x3i             \n"  \
"fstpl %2                 #95     x2r,X2I,x3r,x3i                 \n"  \
"fstpl %4                 #96     X2I,x3r,x3i                     \n"  \
"fstpl %5                 #97     x3r,x3i                         \n"  \
"fstpl %6                 #98     x3i                             \n"  \
"fstpl %7                 #99                                     \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),     \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory");

#define asm_cplx_fft5_store_p_F(_j,_pd0,_pd1,_pd2,_pd3,_pd4,_x0, _x1, _x2, _x3, _x4) \
__asm__ volatile("                                                \n"  \
"fldl %2                  #1      x1r                             \n"  \
"faddl %8                 #2-5    NX1R                            \n"  \
"fldl %2                  #3      x1r,NX1R                        \n"  \
"fsubl %8                 #4-7    NX4R,NX1R                       \n"  \
"fldl %4                  #5      x2r,NX4R,NX1R                   \n"  \
"fsubl %6                 #6-9    NX3R,NX4R,nx1r                  \n"  \
"fldl %4                  #7      x2r,NX3R,NX4R,nx1r              \n"  \
"faddl %6                 #8-11   NX2R,NX3R,n4xr,nx1r             \n"  \
"fldl %3                  #9      x1i,NX2R,NX3R,n4xr,nx1r         \n"  \
"fsubl %9                 #10-13  NX4I,NX2R,nx3r,n4xr,nx1r        \n"  \
"fxch %%st(3)             #10     n4xr,nX2R,nx3r,NX4I,nx1r        \n"  \
"fstpl %8                 #11     nx2r,nx3r,NX4I,nx1r             \n"  \
"fxch %%st(2)             #12     NX4I,nx3r,nx2r,nx1r             \n"  \
"fldl %3                  #13     x1i,NX4I,nx3r,nx2r,nx1r         \n"  \
"faddl %9                 #14-17  NX1I,nx4i,nx3r,nx2r,nx1r        \n"  \
"fxch %%st(2)             #14     nx3r,nx4i,NX1I,nx2r,nx1r        \n"  \
"fldl %5                  #15     x2i,nx3r,nx4i,NX1I,nx2r,nx1r    \n"  \
"fsubl %7                 #16-19  NX3I,nx3r,nx4i,NX1I,nx2r,nx1r   \n"  \
"fxch %%st(2)             #16     nx4i,nx3r,NX3I,nx1i,nx2r,nx1r   \n"  \
"fldl %5                  #17     x2i,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r \n"\
"faddl %7                 #18-21  NX2I,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n"\
"fxch %%st(2)             #19     nx3r,nx4i,NX2I,NX3I,nx1i,nx2r,nx1r\n"\
"fstpl %6                 #20     nx4i,NX2I,NX3I,nx1i,nx2r,nx1r   \n"  \
"fstpl %9                 #21     NX2I,NX3I,nx1i,nx2r,nx1r        \n"  \
"fxch %%st(1)             #22     nx3i,nx2i,nx1i,nx2r,nx1r        \n"  \
"fstpl %7                 #23     x2i,x1i,x2r,x1r                 \n"  \
"fld %%st(3)              #24     x1r,x2i,x1i,x2r,x1r             \n"  \
"fadd %%st(3)             #25-28  AR,x2i,x1i,x2r,x1r              \n"  \
"fld %%st(2)              #26     x1i,AR,x2i,x1i,x2r,x1r          \n"  \
"fadd %%st(2)             #27-30  AI,AR,x2i,x1i,x2r,x1r           \n"  \
"fldl %0                  #28     x0r,AI,AR,x2i,x1i,x2r,x1r       \n"  \
"fadd %%st(2)             #29-32  NX0R,AI,AR,x2i,x1i,x2r,x1r      \n"  \
"fxch %%st(6)             #29     x1r,AI,ar,x2i,x1i,x2r,NX0R      \n"  \
"fsubp %%st,%%st(5)       #30-33  AI,ar,x2i,x1i,NX2R,NX0R         \n"  \
"fxch %%st(3)             #30     x1i,ar,x2i,ai,NX2R,NX0R         \n"  \
"fsubp %%st,%%st(2)       #31-34  ar,NX2I,ai,NX2R,NX0R            \n"  \
"fldl %1                  #32     x0i,ar,NX2I,ai,NX2R,nx0r        \n"  \
"fadd %%st(3)             #33     NX0I,ar,NX2I,ai,nx2r,nx0r       \n"  \
"fxch %%st(4)             #34     nx2r,ar,NX2I,ai,NX0i,nx0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5r) "             #35-38  X2R,ar,nx2i,ai,NX0I,nx0r        \n"  \
"fxch %%st(5)             #36     nx0r,ar,nx2i,ai,NX0I,X2R        \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),     \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory");                                                          \
__asm__ volatile("                                                \n"  \
"fstl (%%eax,%%edi,8)     #37                                     \n"  \
"fstpl %0                 #       ar,nx2i,ai,nx0i,X2R             \n"  \
"fxch %%st(1)             #38     nx2i,ar,ai,nx0i,X2R             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5r) "             #39-42  X2I,ar,ai,nx0i,X2R              \n"  \
"fxch %%st(3)             #40     nx0i,ar,ai,X2I,x2r              \n"  \
"fstl 8(%%eax,%%edi,8)    #41                                     \n"  \
"fstpl %1                 #42     ar,ai,X2I,x2r                   \n"  \
	: :"m"(_x0##r),"m"(_x0##i),"D"(_j),"a"(_pd0):"memory","st");   \
__asm__ volatile("                                                \n"  \
"fmull " MEMORY_VARIABLE_NAME(FNM125) "             #41-44  AR,ai,X2I,x2r                   \n"  \
"fxch %%st(1)             #42     ai,AR,X2i,X2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(FNM125) "             #43-45  AI,AR,x2i,x2r                   \n"  \
"fxch %%st(1)             #43     AR,AI,x2i,x2r                   \n"  \
"faddl %0                 #44-47  NAR,AI,x2i,x2r                  \n"  \
"fxch %%st(1)             #45     AI,NAR,x2i,x2r                  \n"  \
"faddl %1                 #46     NAI,NAR,x2i,x2r                 \n"  \
"fld %%st(1)              #47     nar,NAI,nar,x2i,x2r             \n"  \
"fadd %%st(4)             #48     X1R,NAI,nar,x2i,x2r             \n"  \
"fxch %%st(2)             #46     nar,NAI,X1R,x2i,x2r             \n"  \
"fsubp %%st,%%st(4)       #47-50  nai,X1R,x2i,NX2R                \n"  \
"fld  %%st                #48     nai,nai,X1R,x2i,NX2R,           \n"  \
"fadd %%st(3)             #49-52  X1I,nai,X1R,x2i,NX2R            \n"  \
"fxch %%st(3)             #49     x2i,nai,x1r,X1I,NX2R            \n"  \
"fsubrp %%st,%%st(1)      #50     X2I,x1r,X1I,x2r                 \n"  \
"fxch %%st(3)             #51     x2r,x1r,x1I,X2I                 \n"  \
"fstpl %4                 #52     x1r,x1i,x2i                     \n"  \
"fstpl %2                 #53     x1i,x2i                         \n"  \
"fstpl %3                 #54     x2i                             \n"  \
"fstpl %5                 #55                                     \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),                 \
   "m"(_x2##r),"m"(_x2##i):"memory","st");                             \
__asm__ volatile("                                                \n"  \
"fldl %8		  #56     x4r                             \n"  \
"fsubl %6		  #57-60  AR                              \n"  \
"fldl %9 		  #58     x4i,AR                          \n"  \
"fsubl %7                 #59-63  AI,AR                           \n"  \
"fldl %6                  #60     x3r,AI,AR                       \n"  \
"fxch %%st(2)             #60	  ar,AI,x3r                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_5i) "             #61-64  NAR,Ai,x3r                      \n"  \
"fldl %7                  #62     x3i,NAR,Ai,x3r                  \n"  \
"fxch %%st(2)             #63     ai,NAR,x3i,x3r                  \n"  \
"fmull " MEMORY_VARIABLE_NAME(F_1_5i) "             #64-67  NAI,NAR,x3i,x3r                 \n"  \
"fldl %8		  #65     x4r,NAI,NAR,x3i,x3r             \n"  \
"fxch %%st(4)             #65     x3r,NAI,nar,x3i,x4r             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_5i) "             #66-69  NX3R,NAI,nar,x3i,x4r            \n"  \
"fldl %9                  #67     x4i,NX3R,NAI,nar,x3i,x4r        \n"  \
"fxch %%st(4)             #67     x3i,NX3R,nai,nar,x4i,x4r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_5i) "             #68-71  NX3I,NX3R,nai,nar,x4i,x4r       \n"  \
"fxch %%st(1)             #68     NX3R,NX3I,nai,nar,x4i,x4r       \n"  \
"fadd %%st(3)             #69-72  X3R,NX3I,nai,nar,x4i,x4r        \n"  \
"fxch %%st(5)             #69     x4r,NX3I,nai,nar,x4i,X3R        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5i) "             #70-73  NX4R,NX3I,nai,nar,x4i,X3R       \n"  \
"fxch %%st(1)             #70     nx3i,NX4R,nai,nar,x4i,X3r       \n"  \
"fadd %%st(2)             #71-74  X3I,NX4R,nai,nar,x4i,x3r        \n"  \
"fxch %%st(4)             #71     x4i,NX4r,nai,nar,X3I,x3r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_5i) "             #72-75  NX4I,NX4R,nai,nar,X3I,x3r       \n"  \
"fxch %%st(1)             #72     NX4r,NX4I,nai,nar,X3I,x3r       \n"  \
"fsubrp %%st,%%st(3)      #73-76  NX4I,nai,AR,x3i,x3r             \n"  \
"fldl %2                  #74     x1r,NX4I,nai,AR,x3i,x3r         \n"  \
"fldl %3                  #75     x1i,x1r,NX4I,nai,AR,x3i,X3r     \n"  \
"fxch %%st(2)             #75     nx4i,x1r,x1i,nai,AR,x3i,x3r     \n"  \
"fsubrp %%st,%%st(3)      #76-79  x1r,x1i,AI,AR,x3i,x3r           \n"  \
"fadd %%st(4)		  #77-80  X4R,x1i,AI,AR,x3i,x3r           \n"  \
"fxch %%st(1)             #77     x1i,X4R,AI,AR,x3i,x3r           \n"  \
"fsub %%st(5)             #78-81  X4I,X4R,Ai,ar,x3i,x3r           \n"  \
"fldl %2                  #79     x1r,X4I,X4R,Ai,ar,x3i,x3r       \n"  \
"fsubp %%st,%%st(5)       #80-83  x4i,x4r,ai,ar,X1R,x3r           \n"  \
"fldl %3                  #81     x1i,x4i,x4r,ai,ar,X1R,x3r       \n"  \
"faddp %%st,%%st(6)       #82     x4i,x4r,ai,ar,X1R,X1I           \n"  \
"fstpl 8(%%eax,%%edi,8)   #83     x4r,ai,ar,X1R,X1I               \n"  \
"fstpl (%%eax,%%edi,8)    #84     ai,ar,X1R,X1I                   \n"  \
"fldl %4                  #85     x2r,ai,ar,x1r,x1i               \n"  \
"fadd %%st(1)             #86-89  X3R,ai,ar,x1r,x1i               \n"  \
"fxch %%st(3)             #86     x1r,ai,ar,X3R,x1i               \n"  \
"fldl %5                  #87     x2i,x1r,ai,ar,X3R,x1i           \n"  \
"fsub %%st(3)             #88-91  X3I,x1r,ai,ar,X3R,x1i           \n"  \
"fxch %%st(5)             #89     x1i,x1r,ai,ar,X3R,X3I           \n"  \
"fldl %4                  #90     x2r,x1i,x1r,ai,ar,x3r,x3i       \n"  \
"fsubp %%st,%%st(3)       #91-94  x1i,x1r,X2R,ar,x3r,x3i          \n"  \
"fldl %5                  #92     x2i,x1i,x1r,X2R,ar,x3r,x3i      \n"  \
"faddp %%st,%%st(4)       #93-96  x1i,x1r,X2R,X2I,x3r,x3i         \n"  \
   : :"D"(_j),"a"(_pd4),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),           \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory","st");                                                     \
__asm__ volatile("                                                \n"  \
"fstpl 8(%%eax,%%edi,8)   #94     x1r,x2r,X2I,x3r,x3i             \n"  \
"fstpl (%%eax,%%edi,8)    #95     x2r,X2I,x3r,x3i                 \n"  \
"fstpl (%%ebx,%%edi,8)    #96     X2I,x3r,x3i                     \n"  \
"fstpl 8(%%ebx,%%edi,8)   #97     x3r,x3i                         \n"  \
"fstpl (%%edx,%%edi,8)    #98     x3i                             \n"  \
"fstpl 8(%%edx,%%edi,8)   #99                                     \n"  \
   : :"D"(_j),"a"(_pd1),"b"(_pd2),"d"(_pd3):"memory","st");

#define asm_cplx_fft5_store_p_B(_j,_pd0,_pd1,_pd2,_pd3,_pd4,_x0, _x1, _x2, _x3, _x4) \
__asm__ volatile("                                                \n"  \
"fldl %2                  #1      x1r                             \n"  \
"faddl %8                 #2-5    NX1R                            \n"  \
"fldl %2                  #3      x1r,NX1R                        \n"  \
"fsubl %8                 #4-7    NX4R,NX1R                       \n"  \
"fldl %4                  #5      x2r,NX4R,NX1R                   \n"  \
"fsubl %6                 #6-9    NX3R,NX4R,nx1r                  \n"  \
"fldl %4                  #7      x2r,NX3R,NX4R,nx1r              \n"  \
"faddl %6                 #8-11   NX2R,NX3R,n4xr,nx1r             \n"  \
"fldl %3                  #9      x1i,NX2R,NX3R,n4xr,nx1r         \n"  \
"fsubl %9                 #10-13  NX4I,NX2R,nx3r,n4xr,nx1r        \n"  \
"fxch %%st(3)             #10     n4xr,nX2R,nx3r,NX4I,nx1r        \n"  \
"fstpl %8                 #11     nx2r,nx3r,NX4I,nx1r             \n"  \
"fxch %%st(2)             #12     NX4I,nx3r,nx2r,nx1r             \n"  \
"fldl %3                  #13     x1i,NX4I,nx3r,nx2r,nx1r         \n"  \
"faddl %9                 #14-17  NX1I,nx4i,nx3r,nx2r,nx1r        \n"  \
"fxch %%st(2)             #14     nx3r,nx4i,NX1I,nx2r,nx1r        \n"  \
"fldl %5                  #15     x2i,nx3r,nx4i,NX1I,nx2r,nx1r    \n"  \
"fsubl %7                 #16-19  NX3I,nx3r,nx4i,NX1I,nx2r,nx1r   \n"  \
"fxch %%st(2)             #16     nx4i,nx3r,NX3I,nx1i,nx2r,nx1r   \n"  \
"fldl %5                  #17     x2i,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n" \
"faddl %7                 #18-21  NX2I,nx4i,nx3r,NX3I,nx1i,nx2r,nx1r\n"\
"fxch %%st(2)             #19     nx3r,nx4i,NX2I,NX3I,nx1i,nx2r,nx1r\n"\
"fstpl %6                 #20     nx4i,NX2I,NX3I,nx1i,nx2r,nx1r   \n"  \
"fstpl %9                 #21     NX2I,NX3I,nx1i,nx2r,nx1r        \n"  \
"fxch %%st(1)             #22     nx3i,nx2i,nx1i,nx2r,nx1r        \n"  \
"fstpl %7                 #23     x2i,x1i,x2r,x1r                 \n"  \
"fld %%st(3)              #24     x1r,x2i,x1i,x2r,x1r             \n"  \
"fadd %%st(3)             #25-28  AR,x2i,x1i,x2r,x1r              \n"  \
"fld %%st(2)              #26     x1i,AR,x2i,x1i,x2r,x1r          \n"  \
"fadd %%st(2)             #27-30  AI,AR,x2i,x1i,x2r,x1r           \n"  \
"fldl %0                  #28     x0r,AI,AR,x2i,x1i,x2r,x1r       \n"  \
"fadd %%st(2)             #29-32  NX0R,AI,AR,x2i,x1i,x2r,x1r      \n"  \
"fxch %%st(6)             #29     x1r,AI,ar,x2i,x1i,x2r,NX0R      \n"  \
"fsubp %%st,%%st(5)       #30-33  AI,ar,x2i,x1i,NX2R,NX0R         \n"  \
"fxch %%st(3)             #30     x1i,ar,x2i,ai,NX2R,NX0R         \n"  \
"fsubp %%st,%%st(2)       #31-34  ar,NX2I,ai,NX2R,NX0R            \n"  \
"fldl %1                  #32     x0i,ar,NX2I,ai,NX2R,nx0r        \n"  \
"fadd %%st(3)             #33     NX0I,ar,NX2I,ai,nx2r,nx0r       \n"  \
"fxch %%st(4)             #34     nx2r,ar,NX2I,ai,NX0i,nx0r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5r) "             #35-38  X2R,ar,nx2i,ai,NX0I,nx0r        \n"  \
"fxch %%st(5)             #36     nx0r,ar,nx2i,ai,NX0I,X2R        \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),     \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory","st");                                                     \
__asm__ volatile("                                                \n"  \
"fstl (%%eax,%%edi,8)     #37                                     \n"  \
"fstpl %0                 #       ar,nx2i,ai,nx0i,X2R             \n"  \
"fxch %%st(1)             #38     nx2i,ar,ai,nx0i,X2R             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5r) "             #39-42  X2I,ar,ai,nx0i,X2R              \n"  \
"fxch %%st(3)             #40     nx0i,ar,ai,X2I,x2r              \n"  \
"fstl 8(%%eax,%%edi,8)    #41                                     \n"  \
"fstpl %1                 #42     ar,ai,X2I,x2r                   \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"D"(_j),"a"(_pd0):"memory","st");        \
__asm__ volatile("                                                \n"  \
"fmull " MEMORY_VARIABLE_NAME(BNM125) "             #41-44  AR,ai,X2I,x2r                   \n"  \
"fxch %%st(1)             #42     ai,AR,X2i,X2r                   \n"  \
"fmull " MEMORY_VARIABLE_NAME(BNM125) "             #43-45  AI,AR,x2i,x2r                   \n"  \
"fxch %%st(1)             #43     AR,AI,x2i,x2r                   \n"  \
"faddl %0                 #44-47  NAR,AI,x2i,x2r                  \n"  \
"fxch %%st(1)             #45     AI,NAR,x2i,x2r                  \n"  \
"faddl %1                 #46     NAI,NAR,x2i,x2r                 \n"  \
"fld %%st(1)              #47     nar,NAI,nar,x2i,x2r             \n"  \
"fadd %%st(4)             #48     X1R,NAI,nar,x2i,x2r             \n"  \
"fxch %%st(2)             #46     nar,NAI,X1R,x2i,x2r             \n"  \
"fsubp %%st,%%st(4)       #47-50  nai,X1R,x2i,NX2R                \n"  \
"fld  %%st                #48     nai,nai,X1R,x2i,NX2R,           \n"  \
"fadd %%st(3)             #49-52  X1I,nai,X1R,x2i,NX2R            \n"  \
"fxch %%st(3)             #49     x2i,nai,x1r,X1I,NX2R            \n"  \
"fsubrp %%st,%%st(1)      #50     X2I,x1r,X1I,x2r                 \n"  \
"fxch %%st(3)             #51     x2r,x1r,x1I,X2I                 \n"  \
"fstpl %4                 #52     x1r,x1i,x2i                     \n"  \
"fstpl %2                 #53     x1i,x2i                         \n"  \
"fstpl %3                 #54     x2i                             \n"  \
"fstpl %5                 #55                                     \n"  \
   : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),                 \
   "m"(_x2##r),"m"(_x2##i):"memory","st");                             \
__asm__ volatile("                                                \n"  \
"fldl %8		  #56     x4r                             \n"  \
"fsubl %6		  #57-60  AR                              \n"  \
"fldl %9 		  #58     x4i,AR                          \n"  \
"fsubl %7                 #59-63  AI,AR                           \n"  \
"fldl %6                  #60     x3r,AI,AR                       \n"  \
"fxch %%st(2)             #60	  ar,AI,x3r                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(B_1_5i) "             #61-64  NAR,Ai,x3r                      \n"  \
"fldl %7                  #62     x3i,NAR,Ai,x3r                  \n"  \
"fxch %%st(2)             #63     ai,NAR,x3i,x3r                  \n"  \
"fmull " MEMORY_VARIABLE_NAME(B_1_5i) "             #64-67  NAI,NAR,x3i,x3r                 \n"  \
"fldl %8		  #65     x4r,NAI,NAR,x3i,x3r             \n"  \
"fxch %%st(4)             #65     x3r,NAI,nar,x3i,x4r             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_5i) "             #66-69  NX3R,NAI,nar,x3i,x4r            \n"  \
"fldl %9                  #67     x4i,NX3R,NAI,nar,x3i,x4r        \n"  \
"fxch %%st(4)             #67     x3i,NX3R,nai,nar,x4i,x4r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_5i) "             #68-71  NX3I,NX3R,nai,nar,x4i,x4r       \n"  \
"fxch %%st(1)             #68     NX3R,NX3I,nai,nar,x4i,x4r       \n"  \
"fadd %%st(3)             #69-72  X3R,NX3I,nai,nar,x4i,x4r        \n"  \
"fxch %%st(5)             #69     x4r,NX3I,nai,nar,x4i,X3R        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5i) "             #70-73  NX4R,NX3I,nai,nar,x4i,X3R       \n"  \
"fxch %%st(1)             #70     nx3i,NX4R,nai,nar,x4i,X3r       \n"  \
"fadd %%st(2)             #71-74  X3I,NX4R,nai,nar,x4i,x3r        \n"  \
"fxch %%st(4)             #71     x4i,NX4r,nai,nar,X3I,x3r        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_5i) "             #72-75  NX4I,NX4R,nai,nar,X3I,x3r       \n"  \
"fxch %%st(1)             #72     NX4r,NX4I,nai,nar,X3I,x3r       \n"  \
"fsubrp %%st,%%st(3)      #73-76  NX4I,nai,AR,x3i,x3r             \n"  \
"fldl %2                  #74     x1r,NX4I,nai,AR,x3i,x3r         \n"  \
"fldl %3                  #75     x1i,x1r,NX4I,nai,AR,x3i,X3r     \n"  \
"fxch %%st(2)             #75     nx4i,x1r,x1i,nai,AR,x3i,x3r     \n"  \
"fsubrp %%st,%%st(3)      #76-79  x1r,x1i,AI,AR,x3i,x3r           \n"  \
"fadd %%st(4)		  #77-80  X4R,x1i,AI,AR,x3i,x3r           \n"  \
"fxch %%st(1)             #77     x1i,X4R,AI,AR,x3i,x3r           \n"  \
"fsub %%st(5)             #78-81  X4I,X4R,Ai,ar,x3i,x3r           \n"  \
"fldl %2                  #79     x1r,X4I,X4R,Ai,ar,x3i,x3r       \n"  \
"fsubp %%st,%%st(5)       #80-83  x4i,x4r,ai,ar,X1R,x3r           \n"  \
"fldl %3                  #81     x1i,x4i,x4r,ai,ar,X1R,x3r       \n"  \
"faddp %%st,%%st(6)       #82     x4i,x4r,ai,ar,X1R,X1I           \n"  \
"fstpl 8(%%eax,%%edi,8)   #83     x4r,ai,ar,X1R,X1I               \n"  \
"fstpl (%%eax,%%edi,8)    #84     ai,ar,X1R,X1I                   \n"  \
"fldl %4                  #85     x2r,ai,ar,x1r,x1i               \n"  \
"fadd %%st(1)             #86-89  X3R,ai,ar,x1r,x1i               \n"  \
"fxch %%st(3)             #86     x1r,ai,ar,X3R,x1i               \n"  \
"fldl %5                  #87     x2i,x1r,ai,ar,X3R,x1i           \n"  \
"fsub %%st(3)             #88-91  X3I,x1r,ai,ar,X3R,x1i           \n"  \
"fxch %%st(5)             #89     x1i,x1r,ai,ar,X3R,X3I           \n"  \
"fldl %4                  #90     x2r,x1i,x1r,ai,ar,x3r,x3i       \n"  \
"fsubp %%st,%%st(3)       #91-94  x1i,x1r,X2R,ar,x3r,x3i          \n"  \
"fldl %5                  #92     x2i,x1i,x1r,X2R,ar,x3r,x3i      \n"  \
"faddp %%st,%%st(4)       #93-96  x1i,x1r,X2R,X2I,x3r,x3i         \n"  \
   : :"D"(_j),"a"(_pd4),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),           \
   "m"(_x2##i),"m"(_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):        \
   "memory","st");                                                     \
__asm__ volatile("                                                \n"  \
"fstpl 8(%%eax,%%edi,8)   #94     x1r,x2r,X2I,x3r,x3i             \n"  \
"fstpl (%%eax,%%edi,8)    #95     x2r,X2I,x3r,x3i                 \n"  \
"fstpl (%%ebx,%%edi,8)    #96     X2I,x3r,x3i                     \n"  \
"fstpl 8(%%ebx,%%edi,8)   #97     x3r,x3i                         \n"  \
"fstpl (%%edx,%%edi,8)    #98     x3i                             \n"  \
"fstpl 8(%%edx,%%edi,8)   #99                                     \n"  \
	: :"D"(_j),"a"(_pd1),"b"(_pd2),"d"(_pd3):"memory","st");

#endif

#ifndef USE_ASM

# define cplx_fft5(_S, _x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = _S##N2_5r * ( _x1##r - _x2##r);\
     _x2##i = _S##N2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = _S##_1_5i * (_x4##r - _x3##r);\
     _ai = _S##_1_5i * (_x4##i - _x3##i);\
     _x3##r *= _S##N1_5i;\
     _x3##i *= _S##N1_5i;\
     _x4##r *= _S##N2_5i;\
     _x4##i *= _S##N2_5i;\
     _x3##r += _ar;\
     _x3##i += _ai;\
     _x4##r = _ar - _x4##r;\
     _x4##i = _ai - _x4##i;\
     \
     _ar = _x4##r;\
     _ai = _x4##i;\
     _x4##r = _x1##r + _x3##i;\
     _x4##i = _x1##i - _x3##r;\
     _x1##r -= _x3##i;\
     _x1##i += _x3##r;\
     _x3##r = _x2##r + _ai;\
     _x3##i = _x2##i - _ar;\
     _x2##r -= _ai;\
     _x2##i += _ar;\
}
#else
# define cplx_fft5(_S, _x0, _x1, _x2, _x3, _x4) \
	asm_cplx_fft5_##_S(_x0, _x1, _x2, _x3, _x4)

#endif

#ifndef USE_ASM

# define cplx_fft5_store_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_S, _x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = _S##N2_5r * ( _x1##r - _x2##r);\
     _x2##i = _S##N2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     *(_pd0 + _index) = _x0##r;\
     *(_pd0 + _index + 1) = _x0##i;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = _S##_1_5i * (_x4##r - _x3##r);\
     _ai = _S##_1_5i * (_x4##i - _x3##i);\
     _x3##r *= _S##N1_5i;\
     _x3##i *= _S##N1_5i;\
     _x4##r *= _S##N2_5i;\
     _x4##i *= _S##N2_5i;\
     _x3##r += _ar;\
     _x3##i += _ai;\
     _x4##r = _ar - _x4##r;\
     _x4##i = _ai - _x4##i;\
     *(_pd4 + _index ) = _x1##r + _x3##i;\
     *(_pd4 + _index +1) = _x1##i - _x3##r;\
     *(_pd1 + _index ) = _x1##r - _x3##i;\
     *(_pd1 + _index + 1) = _x1##i + _x3##r;\
     *(_pd3 + _index ) = _x2##r + _x4##i;\
     *(_pd3 + _index + 1) = _x2##i - _x4##r;\
     *(_pd2 + _index ) = _x2##r - _x4##i;\
     *(_pd2 + _index + 1) = _x2##i + _x4##r;\
}

#else

# define cplx_fft5_store_p(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_S, _x0, _x1, _x2, _x3, _x4) \
   asm_cplx_fft5_store_p_##_S(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_x0,_x1,_x2,_x3,_x4);

#endif

/*$Id$*/




