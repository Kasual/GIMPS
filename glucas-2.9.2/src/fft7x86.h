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
/* This is a c-macro style of Mlucas (Ernst Mayer) version for Nusbaumer
   FFT length 7 */

#ifdef USE_ASM
# define asm_cplx_fft7_F(_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
{                                                                \
    BIG_DOUBLE _a0r, _a0i, _a1r, _a1i, _a2r, _a2i,               \
      _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i;            \
    __asm__ volatile("                                      \n"  \
"fldl %4	#1                                          \n"  \
"faddl %6       #2     %0                                   \n"  \
"fldl %4        #3                                          \n"  \
"fsubl %6       #4     %2,%0                                \n"  \
"fldl %5	#5                                          \n"  \
"faddl %7       #6     %1,%2,%0                             \n"  \
"fldl %5        #7                                          \n"  \
"fsubl %7       #8     %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #9     %0,%1,%2,%3                          \n"  \
"fstpl %0       #10                                         \n"  \
"fstpl %1       #11                                         \n"  \
"fstpl %2       #12                                         \n"  \
"fstpl %3       #13                                         \n"  \
 : : "m" (_a1r),"m"(_a1i),"m"(_a6r),"m"(_a6i),                   \
 "m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");     \
    __asm__ volatile("                                      \n"  \
"fldl %4	#14                                         \n"  \
"faddl %6       #15     %0                                  \n"  \
"fldl %4        #16                                         \n"  \
"fsubl %6       #17     %2,%0                               \n"  \
"fldl %5	#18                                         \n"  \
"faddl %7       #19     %1,%2,%0                            \n"  \
"fldl %5        #20                                         \n"  \
"fsubl %7       #21    %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #22     %0,%1,%2,%3                         \n"  \
"fstpl %0       #23                                         \n"  \
"fstpl %1       #24                                         \n"  \
"fstpl %2       #25                                         \n"  \
"fstpl %3       #26                                         \n"  \
 : : "m" (_a2r),"m"(_a2i),"m"(_a5r),"m"(_a5i),                   \
 "m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i):"memory");     \
    __asm__ volatile("                                      \n"  \
"fldl %4	#27                                         \n"  \
"faddl %6       #28    %0                                   \n"  \
"fldl %4        #29                                         \n"  \
"fsubl %6       #30     %2,%0                               \n"  \
"fldl %5	#31                                         \n"  \
"faddl %7       #32     %1,%2,%0                            \n"  \
"fldl %5        #33                                         \n"  \
"fsubl %7       #34     %3,%1,%2,%0                         \n"  \
"fxch %%st(3)   #35     %0,%1,%2,%3                         \n"  \
"fstpl %0       #36                                         \n"  \
"fstpl %1       #37                                         \n"  \
"fstpl %2       #38                                         \n"  \
"fstpl %3       #39                                         \n"  \
 : : "m" (_a3r),"m"(_a3i),"m"(_a4r),"m"(_a4i),                   \
 "m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):"memory");     \
__asm__ volatile("                                          \n"  \
"fldl %2	#40     x1r                                 \n"  \
"faddl %4       #41-44  X12R                                \n"  \
"fldl %3        #42     x1i,X12R                            \n"  \
"faddl %5       #43-46  X12I,X12R                           \n"  \
"fldl %0	#44     x0r,X12I,X12R                       \n"  \
"fldl %1        #45     x0i,x0r,X12I,x12r                   \n"  \
"fxch %%st(3)   #45     x12r,x0r,X12I,x0i                   \n"  \
"faddl %6       #46-49  AR,x0r,x12i,x0i                     \n"  \
"fxch %%st(2)   #46     x12i,x0r,AR,x0i                     \n"  \
"faddl %7       #47-50  AI,x0r,AR,x0i                       \n"  \
"fxch %%st(1)   #47     x0r,AI,AR,x0i                       \n"  \
"fstl %8        #48     x0r,AI,AR,x0i                       \n"  \
"fadd %%st(2)   #49-52  X0R,AI,AR,x0i                       \n"  \
"fxch %%st(3)   #50     x0i,ai,ar,X0R                       \n"  \
"fstl %9        #51     x0i,ai,ar,X0R                       \n"  \
"fadd %%st(1)   #52-55  X0I,ai,ar,x0r                       \n"  \
:  :"m"(_x0##r),"m"(_x0##i),"m"(_a1r),"m"(_a1i),"m"(_a2r),"m"(_a2i),\
    "m"(_a3r),"m"(_a3i),"m"(_a0r),"m"(_a0i):"memory");           \
__asm__ volatile ("                                         \n"  \
"fxch %%st(2)   #52     ar,ai,X0I,X0R                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7r) " #53-56  AR,ai,X0I,x0r                       \n"  \
"fxch %%st(3)   #54     x0r,ai,X0I,AR                       \n"  \
"fstpl %0       #55     ai,X0I,AR                           \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7r) "   #56-59  AI,x0i,AR                           \n"  \
"fxch %%st(1)   #57     x0i,AI,ar                           \n"  \
"fstpl %1       #58     AI,ar                               \n"  \
   : :"m"(_x0##r),"m"(_x0##i):"memory");                         \
__asm__ volatile ("                                         \n"  \
"fxch %%st(1)         #59     ar,AI                         \n"  \
"faddl %0             #60-63  AR,ai                         \n"  \
"fxch %%st(1)         #60     ai,AR                         \n"  \
"faddl %1             #61-64  AI,AR                         \n"  \
"fldl %4	      #62     a2r,AI,AR                     \n"  \
"fsubl %2             #63-66  A0R,AI,AR                     \n"  \
"fldl %5              #64     a2i,A0R,AI,AR                 \n"  \
"fsubl %3             #65-68  A0I,A0R,AI,AR                 \n"  \
"fldl %2              #66     a1r,A0I,A0R,ai,ar             \n"  \
"fsubl %6             #67-70  A1R,A0I,A0R,ai,ar             \n"  \
"fxch %%st(2)         #68     a0r,A0I,A1R,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7r) "         #69-72  A0R,A0I,A1r,ai,ar             \n"  \
"fxch %%st(1)         #70     a0i,A0R,A1r,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7r) "         #71-74  A0I,A0R,A1r,ai,ar             \n"  \
"fldl %3              #72     a1i,A0I,A0R,A1r,ai,ar         \n"  \
"fsubl %7             #73-76  A1i,A0I,A0R,A1r,ai,ar         \n"  \
"fxch %%st(2)         #74     a0r,A0I,A1I,a1r,ai,ar         \n"  \
"fstpl %0             #75     a0i,A1I,a1r,ai,ar             \n"  \
"fstpl %1             #76     A1I,a1r,ai,ar                 \n"  \
"fxch %%st(1)         #77     a1r,a1I,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7r) "         #78-81  A1R,a1I,ai,ar                 \n"  \
"fxch %%st(1)         #79     a1i,A1R,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7r) "         #80-83  A1I,A1R,ai,ar                 \n"  \
"fldl %6              #81     a3r,A1I,A1R,ai,ar             \n"  \
"fsubl %4             #82-85  A2R,A1I,A1R,ai,ar             \n"  \
"fxch %%st(2)         #83     a1r,A1I,A2r,ai,ar             \n"  \
"fstpl %2             #84     a1i,A2r,ai,ar                 \n"  \
"fstpl %3             #85     A2r,ai,ar                     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7r) "         #86-89  A2R,ai,ar                     \n"  \
"fldl %7              #87     a3i,A2R,ai,ar                 \n"  \
"fsubl %5             #88-91  A2I,A2R,ai,ar                 \n"  \
"fldl %2              #89     a1r,A2I,a2r,ai,ar             \n"  \
"fadd %%st(4)         #90     NA3R,A2I,a2r,ai,ar            \n"  \
"fxch %%st(1)         #90     A2I,NA3R,a2r,ai,ar            \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7r) "         #91-94  A2i,NA3R,a2r,ai,ar            \n"  \
"fxch %%st(1)         #91     na3R,A2I,a2r,ai,ar            \n"  \
"fldl %3              #92     a1i,na3r,a2i,a2r,ai,ar        \n"  \
"fadd %%st(4)         #93-96  NA3I,na3r,a2i,a2r,ai,ar       \n"  \
"fxch %%st(1)         #93     na3r,NA3I,a2i,a2r,ai,ar       \n"  \
"fadd %%st(3)         #94     A3R,NA3I,a2i,a2r,ai,ar        \n"  \
"fldl %2              #95     a1r,A3R,NA3I,a2i,a2r,ai,ar    \n"  \
"faddl %0             #96-99  NA1R,A3R,NA3I,a2i,a2r,ai,ar   \n"  \
"fxch %%st(2)         #97     na3i,A3R,NA1R,a2i,a2r,ai,ar   \n"  \
"fadd %%st(3)         #98-101 A3I,a3r,NA1R,a2i,a2r,ai,ar    \n"  \
"fxch %%st(1)         #99     a3r,A3I,NA1R,a2i,a2r,ai,ar    \n"  \
"fstpl %6             #100    A3I,NA1R,a2i,a2r,ai,ar        \n"  \
"fstpl %7             #102    NA1R,a2i,a2r,ai,ar     STALL !\n"  \
"fldl %3              #103    a1i,na1r,a2i,a2r,ai,ar        \n"  \
"faddl %1             #104-107 Na1i,na1r,a2i,a2r,ai,ar      \n"  \
"fxch %%st(1)         #104     na1r,Na1i,a2i,a2r,ai,ar      \n"  \
"fsubr %%st(5)        #105-108 A1R,na1i,a2i,a2r,ai,ar       \n"  \
"fxch %%st(3)         #105     a2r,na1i,a2i,A1R,ai,ar       \n"  \
"fsubl %0             #106-109 NA2R,na1i,a2i,A1R,ai,ar      \n"  \
"fxch %%st(1)         #106     na1i,NA2R,a2i,A1R,ai,ar      \n"  \
"fsubr %%st(4)        #107-110 A1I,NA2R,a2i,a1r,ai,ar       \n"  \
"fxch %%st(2)         #107     a2i,NA2R,A1I,a1r,ai,ar       \n"  \
"fsubl %1             #108     NA2I,na2r,A1I,a1r,ai,ar      \n"  \
"fxch %%st(1)         #109     na2r,NA2I,a1i,a1r,ai,ar      \n"  \
"fsubrp %%st,%%st(5)  #110-113 NA2I,a1i,a1r,ai,A2R          \n"  \
"fsubrp %%st,%%st(3)  #111-114 a1i,a1r,A2I,A2R              \n"  \
"fstpl %3             #112  a1r,A2I,A2r                     \n"  \
"fstpl %2             #113  A2I,A2r                         \n"  \
"fstpl %5             #114  a2r                             \n"  \
"fstpl %4             #115                                  \n"  \
  : :"m"(_a0r),"m"(_a0i),"m"(_a1r),"m"(_a1i),"m"(_a2r),"m"(_a2i),\
	"m"(_a3r),"m"(_a3i):"memory");                           \
__asm__ volatile ("                                         \n"  \
"fldl %0            #116      a4r                           \n"  \
"fsubl %2           #117-120  A45r                          \n"  \
"fldl %1            #118      a4i                           \n"  \
"fsubl %3           #119-122  A45i,A45r                     \n"  \
"fldl %4            #120      a6r,A45i,A45r                 \n"  \
"fsubl %2           #121-124  A65r,A45i,A45r                \n"  \
"fldl %5            #122      a6i,A65r,A45i,a45r            \n"  \
"fsubl %3           #123-126  a65i,A65r,a45i,a45r           \n"  \
"fxch %%st(3)       #123      a45r,A65r,a45i,A65i           \n"  \
"fsubl %4           #124-127  Nar,A65r,a45i,A65i            \n"  \
"fxch %%st(2)       #124      a45i,a65r,Nar,A65i            \n"  \
"fsubl %5           #125-128  Nai,a65r,Nar,A65i             \n"  \
"fxch %%st(1)       #125      a65r,Nai,Nar,A65i             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7i) "       #126-129  A0r,Nai,Nar,A65i              \n"  \
"fldl %4            #127      a6r,A0r,Nai,Nar,A65i          \n"  \
"fxch %%st(4)       #127      a65i,A0r,Nai,Nar,a6r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7i) "       #128-131  A0i,A0r,Nai,Nar,a6r           \n"  \
"fldl %5            #129      a6i,A0i,A0r,Nai,Nar,a6r       \n"  \
"fxch %%st(4)       #130      nar,A0i,A0r,Nai,a6i,a6r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7i) "       #131-134  Ar,A0i,a0r,nai,a6i,a6r        \n"  \
"fxch %%st(5)       #131      a6r,A0i,a0r,nai,a6i,Ar        \n"  \
"faddl %0           #132-135  A6r,a0i,a0r,nai,a6i,Ar        \n"  \
"fxch %%st(3)       #132      nai,a0i,a0r,A6r,a6i,Ar        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7i) "       #133-136  Ai,a0i,a0r,a6r,a6i,Ar         \n"  \
"fxch %%st(4)       #133      a6i,a0i,a0r,a6r,Ai,ar         \n"  \
"faddl %1           #134-137  A6i,a0i,a0r,a6r,Ai,ar         \n"  \
"fldl %2            #135      a5r,A6i,a0i,a0r,a6r,Ai,ar     \n"  \
"faddl %0           #136-139  A5r,A6i,a0i,a0r,a6r,ai,ar     \n"  \
"fxch %%st(4)       #136      a6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7i) "       #137-140  A6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fxch %%st(1)       #138      A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7i) "       #139-142  A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fxch %%st(4)       #140      a5r,A6r,a0i,a0r,A6i,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7i) "       #141-144  A5r,a6r,a0i,a0r,A6i,ai,ar     \n"  \
"fxch %%st(1)       #142      a6r,A5r,a0i,a0r,A6i,ai,ar     \n"  \
"fstpl %4           #143      A5r,a0i,a0r,A6i,ai,ar         \n"  \
"fldl %3            #144      a5i,A5r,a0i,a0r,a6i,ai,ar     \n"  \
"faddl %1           #145-148  A5i,a5r,a0i,a0r,a6i,ai,ar     \n"  \
"fxch %%st(4)       #146      a6i,a5r,a0i,a0r,A5i,ai,ar     \n"  \
"fstpl %5           #147      a5r,a0i,a0r,A5i,ai,ar         \n"  \
"fstpl %2           #148      a0i,a0r,A5i,ai,ar             \n"  \
"fxch %%st(2)       #148      a5i,a0r,a0i,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7i) "       #149-152  A5i,a0r,a0i,ai,ar             \n"  \
"fld %%st(1)        #150      a0r,A5i,a0r,a0i,ai,ar         \n"  \
"fsubl %2           #151-154  Na4,A5i,a0r,a0i,ai,ar         \n"  \
"fld %%st(3)        #152      a0i,Na4r,a5i,a0r,a0i,ai,ar    \n"  \
"fsub %%st(2)       #153-156  Na4i,Na4r,a5i,a0r,a0i,ai,ar   \n"  \
"fxch %%st(1)       #153      Na4r,Na4i,a5i,a0r,a0i,ai,ar   \n"  \
"fadd %%st(6)       #154-157  A4r,Na4i,a5i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(2)       #154      a5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #155-158  A5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #155      Na4i,A5i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #156-159  A4i,A5i,A4r,a0r,a0i,ai,ar     \n"  \
"fxch %%st(2)       #157      A4r,A5i,A4i,a0r,a0i,ai,ar     \n"  \
"fstpl %0           #158      A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fsubl %5           #159-162  A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fld %%st(5)        #160      ar,A5i,A4i,a0r,a0i,ai,ar      \n"  \
"fsubl %4           #161      Na5r,a5i,a4i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #162      a5i,Na5r,a4i,a0r,a0i,ai,ar    \n"  \
"fstpl %3           #163      Na5r,a4i,a0r,a0i,ai,ar        \n"  \
"fxch %%st(1)       #164      a4i,Na5r,a0r,a0i,ai,ar        \n"  \
"fstpl %1           #165      Na5r,a0r,a0i,ai,ar            \n"  \
"faddl %2           #166      A5r,a0r,a0i,ai,ar             \n"  \
"fxch %%st(4)       #166      ar,a0r,a0i,ai,A5r             \n"  \
"fsubp %%st,%%st(1) #167-170  Na6r,a0i,ai,A5r               \n"  \
"fxch %%st(2)       #167      ai,a0i,Na6r,A5r               \n"  \
"fsubp %%st,%%st(1) #168-171  Na6i,Na6r,A5r                 \n"  \
"fxch %%st(2)       #169      A5r,Na6r,Na6i                 \n"  \
"fstpl %2           #170      Na6r,Na6i                     \n"  \
"faddl %4           #171-174  A6r,Na6i                      \n"  \
"fxch %%st(1)       #171      Na6i,A6r                      \n"  \
"faddl %5           #172      A6i,A6r                       \n"  \
  : :"m"(_a4r),"m"(_a4i),"m"(_a5r),"m"(_a5i),"m"(_a6r),"m"(_a6i):\
  "memory");                                                     \
__asm__ volatile ("                                         \n"  \
"fldl %0            #173      a3r,a6i,a6r                   \n"  \
"faddl %3           #174-177  X1r,a6i,a6r                   \n"  \
"fldl %0            #175      a3r,X1r,a6i,a6r               \n"  \
"fsubl %3           #176-179  X6r,X1r,a6i,a6r               \n"  \
"fldl %1            #177      a3i,X6r,X1r,a6i,a6r           \n"  \
"fsubl %2           #178-181  X1i,X6r,x1r,a6i,a6r           \n"  \
"fldl %1            #179      a3i,X1i,x6r,x1r,a6i,a6r       \n"  \
"faddl %2           #180-183  X6i,X1i,x6r,x1r,a6i,a6r       \n"  \
"fxch %%st(3)       #181      x1r,x1i,x6r,x6i,a6i,a6r       \n"  \
"fstpl %4           #182                                    \n"  \
"fstpl %5           #183                                    \n"  \
"fstpl %6           #184                                    \n"  \
"fstpl %7           #185      a6i,a6r                       \n"  \
  : :"m"(_a3r),"m"(_a3i),"m"(_a5r),"m"(_a5i),                    \
 "m"(_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");      \
__asm__ volatile ("                                         \n"  \
"fldl %0            #186      a1r,a6i,a6r                   \n"  \
"fadd %%st(1)       #187-190  X2r,a6i,a6r                   \n"  \
"fldl %0            #188      a1r,X2r,a6i,a6r               \n"  \
"fsubp %%st,%%st(2) #189-192  X2r,x5r,a6r                   \n"  \
"fldl %1            #190      a1i,X2r,x5r,a6r               \n"  \
"fsub %%st(3)       #191      X2i,X2r,x5r,a6r               \n"  \
"fldl %1            #192      a1i,X2i,x2r,x5r,a6r           \n"  \
"faddp %%st,%%st(4) #191-194  x2i,x2r,X5r,x5i               \n"  \
"fstpl %3           #192                                    \n"  \
"fstpl %2           #193                                    \n"  \
"fstpl %4           #195                                    \n"  \
"fstpl %5           #196                                    \n"  \
  : :"m"(_a1r),"m"(_a1i),"m"(_x2##r),"m"(_x2##i),"m"(_x5##r),    \
  "m"(_x5##i): "memory");                                        \
__asm__ volatile ("                                         \n"  \
"fldl %0            #197      a2r                           \n"  \
"faddl %3           #198-201  X4r                           \n"  \
"fldl %0            #199      a2r,X4r                       \n"  \
"fsubl %3           #200-203  X3r,X4r                       \n"  \
"fldl %1            #201      a2i,X3r,X4r                   \n"  \
"fsubl %2           #202-205  X4i,X3r,x4r                   \n"  \
"fldl %1            #206      a2i,X4i,x3r,x4r               \n"  \
"faddl %2           #207-210  X3i,X4i,x3r,x4r               \n"  \
"fxch %%st(3)       #208      x4r,x4i,x3r,x3i               \n"  \
"fstpl %4           #209                                    \n"  \
"fstpl %5           #210                                    \n"  \
"fstpl %6           #211                                    \n"  \
"fstpl %7           #212                                    \n"  \
  : :"m"(_a2r),"m"(_a2i),"m"(_a4r),"m"(_a4i),                    \
  "m"(_x4##r),"m"(_x4##i),"m"(_x3##r),"m"(_x3##i):"memory");     \
}




# define asm_cplx_fft7_B(_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
{                                                                \
    BIG_DOUBLE _a0r, _a0i, _a1r, _a1i, _a2r, _a2i,               \
      _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i;            \
    __asm__ volatile("                                      \n"  \
"fldl %4	#1                                          \n"  \
"faddl %6       #2     %0                                   \n"  \
"fldl %4        #3                                          \n"  \
"fsubl %6       #4     %2,%0                                \n"  \
"fldl %5	#5                                          \n"  \
"faddl %7       #6     %1,%2,%0                             \n"  \
"fldl %5        #7                                          \n"  \
"fsubl %7       #8     %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #9     %0,%1,%2,%3                          \n"  \
"fstpl %0       #10                                         \n"  \
"fstpl %1       #11                                         \n"  \
"fstpl %2       #12                                         \n"  \
"fstpl %3       #13                                         \n"  \
  : : "m" (_a1r),"m"(_a1i),"m"(_a6r),"m"(_a6i),                  \
  "m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");    \
    __asm__ volatile("                                      \n"  \
"fldl %4	#14                                         \n"  \
"faddl %6       #15     %0                                  \n"  \
"fldl %4        #16                                         \n"  \
"fsubl %6       #17     %2,%0                               \n"  \
"fldl %5	#18                                         \n"  \
"faddl %7       #19     %1,%2,%0                            \n"  \
"fldl %5        #20                                         \n"  \
"fsubl %7       #21    %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #22     %0,%1,%2,%3                         \n"  \
"fstpl %0       #23                                         \n"  \
"fstpl %1       #24                                         \n"  \
"fstpl %2       #25                                         \n"  \
"fstpl %3       #26                                         \n"  \
 : : "m" (_a2r),"m"(_a2i),"m"(_a5r),"m"(_a5i),                   \
 "m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i):"memory");     \
    __asm__ volatile("                                      \n"  \
"fldl %4	#27                                         \n"  \
"faddl %6       #28    %0                                   \n"  \
"fldl %4        #29                                         \n"  \
"fsubl %6       #30     %2,%0                               \n"  \
"fldl %5	#31                                         \n"  \
"faddl %7       #32     %1,%2,%0                            \n"  \
"fldl %5        #33                                         \n"  \
"fsubl %7       #34     %3,%1,%2,%0                         \n"  \
"fxch %%st(3)   #35     %0,%1,%2,%3                         \n"  \
"fstpl %0       #36                                         \n"  \
"fstpl %1       #37                                         \n"  \
"fstpl %2       #38                                         \n"  \
"fstpl %3       #39                                         \n"  \
 : : "m" (_a3r),"m"(_a3i),"m"(_a4r),"m"(_a4i),                   \
 "m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):"memory");     \
__asm__ volatile("                                          \n"  \
"fldl %2	#40     x1r                                 \n"  \
"faddl %4       #41-44  X12R                                \n"  \
"fldl %3        #42     x1i,X12R                            \n"  \
"faddl %5       #43-46  X12I,X12R                           \n"  \
"fldl %0	#44     x0r,X12I,X12R                       \n"  \
"fldl %1        #45     x0i,x0r,X12I,x12r                   \n"  \
"fxch %%st(3)   #45     x12r,x0r,X12I,x0i                   \n"  \
"faddl %6       #46-49  AR,x0r,x12i,x0i                     \n"  \
"fxch %%st(2)   #46     x12i,x0r,AR,x0i                     \n"  \
"faddl %7       #47-50  AI,x0r,AR,x0i                       \n"  \
"fxch %%st(1)   #47     x0r,AI,AR,x0i                       \n"  \
"fstl %8        #48     x0r,AI,AR,x0i                       \n"  \
"fadd %%st(2)   #49-52  X0R,AI,AR,x0i                       \n"  \
"fxch %%st(3)   #50     x0i,ai,ar,X0R                       \n"  \
"fstl %9        #51     x0i,ai,ar,X0R                       \n"  \
"fadd %%st(1)   #52-55  X0I,ai,ar,x0r                       \n"  \
 : :"m"(_x0##r),"m"(_x0##i),"m"(_a1r),"m"(_a1i),"m"(_a2r),"m"(_a2i),\
    "m"(_a3r),"m"(_a3i),"m"(_a0r),"m"(_a0i):"memory");           \
__asm__ volatile ("                                         \n"  \
"fxch %%st(2)   #52     ar,ai,X0I,X0R                       \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7r) "   #53-56  AR,ai,X0I,x0r                       \n"  \
"fxch %%st(3)   #54     x0r,ai,X0I,AR                       \n"  \
"fstpl %0       #55     ai,X0I,AR                           \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7r) "   #56-59  AI,x0i,AR                           \n"  \
"fxch %%st(1)   #57     x0i,AI,ar                           \n"  \
"fstpl %1       #58     AI,ar                               \n"  \
  : :"m"(_x0##r),"m"(_x0##i):"memory");                          \
__asm__ volatile ("                                         \n"  \
"fxch %%st(1)         #59     ar,AI                         \n"  \
"faddl %0             #60-63  AR,ai                         \n"  \
"fxch %%st(1)         #60     ai,AR                         \n"  \
"faddl %1             #61-64  AI,AR                         \n"  \
"fldl %4	      #62     a2r,AI,AR                     \n"  \
"fsubl %2             #63-66  A0R,AI,AR                     \n"  \
"fldl %5              #64     a2i,A0R,AI,AR                 \n"  \
"fsubl %3             #65-68  A0I,A0R,AI,AR                 \n"  \
"fldl %2              #66     a1r,A0I,A0R,ai,ar             \n"  \
"fsubl %6             #67-70  A1R,A0I,A0R,ai,ar             \n"  \
"fxch %%st(2)         #68     a0r,A0I,A1R,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7r) "         #69-72  A0R,A0I,A1r,ai,ar             \n"  \
"fxch %%st(1)         #70     a0i,A0R,A1r,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7r) "         #71-74  A0I,A0R,A1r,ai,ar             \n"  \
"fldl %3              #72     a1i,A0I,A0R,A1r,ai,ar         \n"  \
"fsubl %7             #73-76  A1i,A0I,A0R,A1r,ai,ar         \n"  \
"fxch %%st(2)         #74     a0r,A0I,A1I,a1r,ai,ar         \n"  \
"fstpl %0             #75     a0i,A1I,a1r,ai,ar             \n"  \
"fstpl %1             #76     A1I,a1r,ai,ar                 \n"  \
"fxch %%st(1)         #77     a1r,a1I,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7r) "         #78-81  A1R,a1I,ai,ar                 \n"  \
"fxch %%st(1)         #79     a1i,A1R,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7r) "         #80-83  A1I,A1R,ai,ar                 \n"  \
"fldl %6              #81     a3r,A1I,A1R,ai,ar             \n"  \
"fsubl %4             #82-85  A2R,A1I,A1R,ai,ar             \n"  \
"fxch %%st(2)         #83     a1r,A1I,A2r,ai,ar             \n"  \
"fstpl %2             #84     a1i,A2r,ai,ar                 \n"  \
"fstpl %3             #85     A2r,ai,ar                     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7r) "         #86-89  A2R,ai,ar                     \n"  \
"fldl %7              #87     a3i,A2R,ai,ar                 \n"  \
"fsubl %5             #88-91  A2I,A2R,ai,ar                 \n"  \
"fldl %2              #89     a1r,A2I,a2r,ai,ar             \n"  \
"fadd %%st(4)         #90     NA3R,A2I,a2r,ai,ar            \n"  \
"fxch %%st(1)         #90     A2I,NA3R,a2r,ai,ar            \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7r) "         #91-94  A2i,NA3R,a2r,ai,ar            \n"  \
"fxch %%st(1)         #91     na3R,A2I,a2r,ai,ar            \n"  \
"fldl %3              #92     a1i,na3r,a2i,a2r,ai,ar        \n"  \
"fadd %%st(4)         #93-96  NA3I,na3r,a2i,a2r,ai,ar       \n"  \
"fxch %%st(1)         #93     na3r,NA3I,a2i,a2r,ai,ar       \n"  \
"fadd %%st(3)         #94     A3R,NA3I,a2i,a2r,ai,ar        \n"  \
"fldl %2              #95     a1r,A3R,NA3I,a2i,a2r,ai,ar    \n"  \
"faddl %0             #96-99  NA1R,A3R,NA3I,a2i,a2r,ai,ar   \n"  \
"fxch %%st(2)         #97     na3i,A3R,NA1R,a2i,a2r,ai,ar   \n"  \
"fadd %%st(3)         #98-101 A3I,a3r,NA1R,a2i,a2r,ai,ar    \n"  \
"fxch %%st(1)         #99     a3r,A3I,NA1R,a2i,a2r,ai,ar    \n"  \
"fstpl %6             #100    A3I,NA1R,a2i,a2r,ai,ar        \n"  \
"fstpl %7             #102    NA1R,a2i,a2r,ai,ar     STALL !\n"  \
"fldl %3              #103    a1i,na1r,a2i,a2r,ai,ar        \n"  \
"faddl %1             #104-107 Na1i,na1r,a2i,a2r,ai,ar      \n"  \
"fxch %%st(1)         #104     na1r,Na1i,a2i,a2r,ai,ar      \n"  \
"fsubr %%st(5)        #105-108 A1R,na1i,a2i,a2r,ai,ar       \n"  \
"fxch %%st(3)         #105     a2r,na1i,a2i,A1R,ai,ar       \n"  \
"fsubl %0             #106-109 NA2R,na1i,a2i,A1R,ai,ar      \n"  \
"fxch %%st(1)         #106     na1i,NA2R,a2i,A1R,ai,ar      \n"  \
"fsubr %%st(4)        #107-110 A1I,NA2R,a2i,a1r,ai,ar       \n"  \
"fxch %%st(2)         #107     a2i,NA2R,A1I,a1r,ai,ar       \n"  \
"fsubl %1             #108     NA2I,na2r,A1I,a1r,ai,ar      \n"  \
"fxch %%st(1)         #109     na2r,NA2I,a1i,a1r,ai,ar      \n"  \
"fsubrp %%st,%%st(5)  #110-113 NA2I,a1i,a1r,ai,A2R          \n"  \
"fsubrp %%st,%%st(3)  #111-114 a1i,a1r,A2I,A2R              \n"  \
"fstpl %3             #112  a1r,A2I,A2r                     \n"  \
"fstpl %2             #113  A2I,A2r                         \n"  \
"fstpl %5             #114  a2r                             \n"  \
"fstpl %4             #115                                  \n"  \
  : :"m"(_a0r),"m"(_a0i),"m"(_a1r),"m"(_a1i),"m"(_a2r),"m"(_a2i),\
	"m"(_a3r),"m"(_a3i):"memory");                           \
__asm__ volatile ("                                         \n"  \
"fldl %0            #116      a4r                           \n"  \
"fsubl %2           #117-120  A45r                          \n"  \
"fldl %1            #118      a4i                           \n"  \
"fsubl %3           #119-122  A45i,A45r                     \n"  \
"fldl %4            #120      a6r,A45i,A45r                 \n"  \
"fsubl %2           #121-124  A65r,A45i,A45r                \n"  \
"fldl %5            #122      a6i,A65r,A45i,a45r            \n"  \
"fsubl %3           #123-126  a65i,A65r,a45i,a45r           \n"  \
"fxch %%st(3)       #123      a45r,A65r,a45i,A65i           \n"  \
"fsubl %4           #124-127  Nar,A65r,a45i,A65i            \n"  \
"fxch %%st(2)       #124      a45i,a65r,Nar,A65i            \n"  \
"fsubl %5           #125-128  Nai,a65r,Nar,A65i             \n"  \
"fxch %%st(1)       #125      a65r,Nai,Nar,A65i             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7i) "       #126-129  A0r,Nai,Nar,A65i              \n"  \
"fldl %4            #127      a6r,A0r,Nai,Nar,A65i          \n"  \
"fxch %%st(4)       #127      a65i,A0r,Nai,Nar,a6r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7i) "       #128-131  A0i,A0r,Nai,Nar,a6r           \n"  \
"fldl %5            #129      a6i,A0i,A0r,Nai,Nar,a6r       \n"  \
"fxch %%st(4)       #130      nar,A0i,A0r,Nai,a6i,a6r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7i) "       #131-134  Ar,A0i,a0r,nai,a6i,a6r        \n"  \
"fxch %%st(5)       #131      a6r,A0i,a0r,nai,a6i,Ar        \n"  \
"faddl %0           #132-135  A6r,a0i,a0r,nai,a6i,Ar        \n"  \
"fxch %%st(3)       #132      nai,a0i,a0r,A6r,a6i,Ar        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7i) "       #133-136  Ai,a0i,a0r,a6r,a6i,Ar         \n"  \
"fxch %%st(4)       #133      a6i,a0i,a0r,a6r,Ai,ar         \n"  \
"faddl %1           #134-137  A6i,a0i,a0r,a6r,Ai,ar         \n"  \
"fldl %2            #135      a5r,A6i,a0i,a0r,a6r,Ai,ar     \n"  \
"faddl %0           #136-139  A5r,A6i,a0i,a0r,a6r,ai,ar     \n"  \
"fxch %%st(4)       #136      a6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7i) "       #137-140  A6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fxch %%st(1)       #138      A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7i) "       #139-142  A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fxch %%st(4)       #140      a5r,A6r,a0i,a0r,A6i,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7i) "       #141-144  A5r,a6r,a0i,a0r,A6i,ai,ar     \n"  \
"fxch %%st(1)       #142      a6r,A5r,a0i,a0r,A6i,ai,ar     \n"  \
"fstpl %4           #143      A5r,a0i,a0r,A6i,ai,ar         \n"  \
"fldl %3            #144      a5i,A5r,a0i,a0r,a6i,ai,ar     \n"  \
"faddl %1           #145-148  A5i,a5r,a0i,a0r,a6i,ai,ar     \n"  \
"fxch %%st(4)       #146      a6i,a5r,a0i,a0r,A5i,ai,ar     \n"  \
"fstpl %5           #147      a5r,a0i,a0r,A5i,ai,ar         \n"  \
"fstpl %2           #148      a0i,a0r,A5i,ai,ar             \n"  \
"fxch %%st(2)       #148      a5i,a0r,a0i,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7i) "       #149-152  A5i,a0r,a0i,ai,ar             \n"  \
"fld %%st(1)        #150      a0r,A5i,a0r,a0i,ai,ar         \n"  \
"fsubl %2           #151-154  Na4,A5i,a0r,a0i,ai,ar         \n"  \
"fld %%st(3)        #152      a0i,Na4r,a5i,a0r,a0i,ai,ar    \n"  \
"fsub %%st(2)       #153-156  Na4i,Na4r,a5i,a0r,a0i,ai,ar   \n"  \
"fxch %%st(1)       #153      Na4r,Na4i,a5i,a0r,a0i,ai,ar   \n"  \
"fadd %%st(6)       #154-157  A4r,Na4i,a5i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(2)       #154      a5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #155-158  A5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #155      Na4i,A5i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #156-159  A4i,A5i,A4r,a0r,a0i,ai,ar     \n"  \
"fxch %%st(2)       #157      A4r,A5i,A4i,a0r,a0i,ai,ar     \n"  \
"fstpl %0           #158      A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fsubl %5           #159-162  A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fld %%st(5)        #160      ar,A5i,A4i,a0r,a0i,ai,ar      \n"  \
"fsubl %4           #161      Na5r,a5i,a4i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #162      a5i,Na5r,a4i,a0r,a0i,ai,ar    \n"  \
"fstpl %3           #163      Na5r,a4i,a0r,a0i,ai,ar        \n"  \
"fxch %%st(1)       #164      a4i,Na5r,a0r,a0i,ai,ar        \n"  \
"fstpl %1           #165      Na5r,a0r,a0i,ai,ar            \n"  \
"faddl %2           #166      A5r,a0r,a0i,ai,ar             \n"  \
"fxch %%st(4)       #166      ar,a0r,a0i,ai,A5r             \n"  \
"fsubp %%st,%%st(1) #167-170  Na6r,a0i,ai,A5r               \n"  \
"fxch %%st(2)       #167      ai,a0i,Na6r,A5r               \n"  \
"fsubp %%st,%%st(1) #168-171  Na6i,Na6r,A5r                 \n"  \
"fxch %%st(2)       #169      A5r,Na6r,Na6i                 \n"  \
"fstpl %2           #170      Na6r,Na6i                     \n"  \
"faddl %4           #171-174  A6r,Na6i                      \n"  \
"fxch %%st(1)       #171      Na6i,A6r                      \n"  \
"faddl %5           #172      A6i,A6r                       \n"  \
  : :"m"(_a4r),"m"(_a4i),"m"(_a5r),"m"(_a5i),"m"(_a6r),"m"(_a6i):"memory");\
__asm__ volatile ("                                         \n"  \
"fldl %0            #173      a3r,a6i,a6r                   \n"  \
"faddl %3           #174-177  X1r,a6i,a6r                   \n"  \
"fldl %0            #175      a3r,X1r,a6i,a6r               \n"  \
"fsubl %3           #176-179  X6r,X1r,a6i,a6r               \n"  \
"fldl %1            #177      a3i,X6r,X1r,a6i,a6r           \n"  \
"fsubl %2           #178-181  X1i,X6r,x1r,a6i,a6r           \n"  \
"fldl %1            #179      a3i,X1i,x6r,x1r,a6i,a6r       \n"  \
"faddl %2           #180-183  X6i,X1i,x6r,x1r,a6i,a6r       \n"  \
"fxch %%st(3)       #181      x1r,x1i,x6r,x6i,a6i,a6r       \n"  \
"fstpl %4           #182                                    \n"  \
"fstpl %5           #183                                    \n"  \
"fstpl %6           #184                                    \n"  \
"fstpl %7           #185      a6r,a6i                       \n"  \
  : :"m"(_a3r),"m"(_a3i),"m"(_a5r),"m"(_a5i),                    \
  "m"(_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");     \
__asm__ volatile ("                                         \n"  \
"fldl %0            #186      a1r,a6i,a6r                   \n"  \
"fadd %%st(1)       #187-190  X2r,a6i,a6r                   \n"  \
"fldl %0            #188      a1r,X2r,a6i,a6r               \n"  \
"fsubp %%st,%%st(2) #189-192  X2r,x5r,a6r                   \n"  \
"fldl %1            #190      a1i,X2r,x5r,a6r               \n"  \
"fsub %%st(3)       #191      X2i,X2r,x5r,a6r               \n"  \
"fldl %1            #192      a1i,X2i,x2r,x5r,a6r           \n"  \
"faddp %%st,%%st(4) #191-194  x2i,x2r,X5r,x5i               \n"  \
"fstpl %3           #192                                    \n"  \
"fstpl %2           #193                                    \n"  \
"fstpl %4           #195                                    \n"  \
"fstpl %5           #196                                    \n"  \
  : :"m"(_a1r),"m"(_a1i),"m"(_x2##r),"m"(_x2##i),"m"(_x5##r),    \
   "m"(_x5##i): "memory");                                       \
__asm__ volatile ("                                         \n"  \
"fldl %0            #197      a2r                           \n"  \
"faddl %3           #198-201  X4r                           \n"  \
"fldl %0            #199      a2r,X4r                       \n"  \
"fsubl %3           #200-203  X3r,X4r                       \n"  \
"fldl %1            #201      a2i,X3r,X4r                   \n"  \
"fsubl %2           #202-205  X4i,X3r,x4r                   \n"  \
"fldl %1            #206      a2i,X4i,x3r,x4r               \n"  \
"faddl %2           #207-210  X3i,X4i,x3r,x4r               \n"  \
"fxch %%st(3)       #208      x4r,x4i,x3r,x3i               \n"  \
"fstpl %4           #209                                    \n"  \
"fstpl %5           #210                                    \n"  \
"fstpl %6           #211                                    \n"  \
"fstpl %7           #212                                    \n"  \
  : :"m"(_a2r),"m"(_a2i),"m"(_a4r),"m"(_a4i),                    \
  "m"(_x4##r),"m"(_x4##i),"m"(_x3##r),"m"(_x3##i):"memory");     \
}



# define asm_cplx_fft7_store_p_B(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3, _x4,_x5,_x6)\
{                                                                \
    __asm__ volatile("                                      \n"  \
"fldl %4	#1                                          \n"  \
"faddl %6       #2     %0                                   \n"  \
"fldl %4        #3                                          \n"  \
"fsubl %6       #4     %2,%0                                \n"  \
"fldl %5	#5                                          \n"  \
"faddl %7       #6     %1,%2,%0                             \n"  \
"fldl %5        #7                                          \n"  \
"fsubl %7       #8     %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #9     %0,%1,%2,%3                          \n"  \
"fstpl %0       #10                                         \n"  \
"fstpl %1       #11                                         \n"  \
"fstpl %2       #12                                         \n"  \
"fstpl %3       #13                                         \n"  \
  : :"m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i),           \
 "m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");     \
    __asm__ volatile("                                      \n"  \
"fldl %4	#14                                         \n"  \
"faddl %6       #15     %0                                  \n"  \
"fldl %4        #16                                         \n"  \
"fsubl %6       #17     %2,%0                               \n"  \
"fldl %5	#18                                         \n"  \
"faddl %7       #19     %1,%2,%0                            \n"  \
"fldl %5        #20                                         \n"  \
"fsubl %7       #21    %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #22     %0,%1,%2,%3                         \n"  \
"fstpl %0       #23                                         \n"  \
"fstpl %1       #24                                         \n"  \
"fstpl %2       #25                                         \n"  \
"fstpl %3       #26                                         \n"  \
  : :"m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i),           \
  "m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i):"memory");    \
    __asm__ volatile("                                      \n"  \
"fldl %4	#27                                         \n"  \
"faddl %6       #28    %0                                   \n"  \
"fldl %4        #29                                         \n"  \
"fsubl %6       #30     %2,%0                               \n"  \
"fldl %5	#31                                         \n"  \
"faddl %7       #32     %1,%2,%0                            \n"  \
"fldl %5        #33                                         \n"  \
"fsubl %7       #34     %3,%1,%2,%0                         \n"  \
"fxch %%st(3)   #35     %0,%1,%2,%3                         \n"  \
"fstpl %0       #36                                         \n"  \
"fstpl %1       #37                                         \n"  \
"fstpl %2       #38                                         \n"  \
"fstpl %3       #39                                         \n"  \
 : :"m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i),            \
 "m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):"memory");     \
__asm__ volatile("                                          \n"  \
"fldl %2	#40     x1r                                 \n"  \
"faddl %4       #41-44  X12R                                \n"  \
"fldl %3        #42     x1i,X12R                            \n"  \
"faddl %5       #43-46  X12I,X12R                           \n"  \
"fldl %0	#44     x0r,X12I,X12R                       \n"  \
"fldl %1        #45     x0i,x0r,X12I,x12r                   \n"  \
"fxch %%st(3)   #45     x12r,x0r,X12I,x0i                   \n"  \
"faddl %6       #46-49  AR,x0r,x12i,x0i                     \n"  \
"fxch %%st(2)   #46     x12i,x0r,AR,x0i                     \n"  \
"faddl %7       #47-50  AI,x0r,AR,x0i                       \n"  \
"fxch %%st(1)   #47     x0r,AI,AR,x0i                       \n"  \
"fstl %8        #48     x0r,AI,AR,x0i                       \n"  \
"fadd %%st(2)   #49-52  X0R,AI,AR,x0i                       \n"  \
"fxch %%st(3)   #50     x0i,ai,ar,X0R                       \n"  \
"fstl %9        #51     x0i,ai,ar,X0R                       \n"  \
"fadd %%st(1)   #52-55  X0I,ai,ar,x0r                       \n"  \
:  :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),"m"(_x2##i),\
    "m"(_x3##r),"m"(_x3##i),"m"(_x0##r),"m"(_x0##i):"memory");   \
__asm__ volatile ("                                         \n"  \
"fxch %%st(2)          #52     ar,ai,X0I,X0R                \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7r) "          #53-56  AR,ai,X0I,x0r                \n"  \
"fxch %%st(3)          #54     x0r,ai,X0I,AR                \n"  \
"fstpl (%%eax,%%edi,8) #55     ai,X0I,AR                    \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7r) "          #56-59  AI,x0i,AR                    \n"  \
"fxch %%st(1)          #57     x0i,AI,ar                    \n"  \
"fstpl 8(%%eax,%%edi,8)#58     AI,ar                        \n"  \
  : :"D"(_d),"a"(_pd0):"memory");                                \
__asm__ volatile ("                                         \n"  \
"fxch %%st(1)         #59     ar,AI                         \n"  \
"faddl %0             #60-63  AR,ai                         \n"  \
"fxch %%st(1)         #60     ai,AR                         \n"  \
"faddl %1             #61-64  AI,AR                         \n"  \
"fldl %4	      #62     a2r,AI,AR                     \n"  \
"fsubl %2             #63-66  A0R,AI,AR                     \n"  \
"fldl %5              #64     a2i,A0R,AI,AR                 \n"  \
"fsubl %3             #65-68  A0I,A0R,AI,AR                 \n"  \
"fldl %2              #66     a1r,A0I,A0R,ai,ar             \n"  \
"fsubl %6             #67-70  A1R,A0I,A0R,ai,ar             \n"  \
"fxch %%st(2)         #68     a0r,A0I,A1R,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7r) "         #69-72  A0R,A0I,A1r,ai,ar             \n"  \
"fxch %%st(1)         #70     a0i,A0R,A1r,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7r) "         #71-74  A0I,A0R,A1r,ai,ar             \n"  \
"fldl %3              #72     a1i,A0I,A0R,A1r,ai,ar         \n"  \
"fsubl %7             #73-76  A1i,A0I,A0R,A1r,ai,ar         \n"  \
"fxch %%st(2)         #74     a0r,A0I,A1I,a1r,ai,ar         \n"  \
"fstpl %0             #75     a0i,A1I,a1r,ai,ar             \n"  \
"fstpl %1             #76     A1I,a1r,ai,ar                 \n"  \
"fxch %%st(1)         #77     a1r,a1I,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7r) "         #78-81  A1R,a1I,ai,ar                 \n"  \
"fxch %%st(1)         #79     a1i,A1R,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7r) "         #80-83  A1I,A1R,ai,ar                 \n"  \
"fldl %6              #81     a3r,A1I,A1R,ai,ar             \n"  \
"fsubl %4             #82-85  A2R,A1I,A1R,ai,ar             \n"  \
"fxch %%st(2)         #83     a1r,A1I,A2r,ai,ar             \n"  \
"fstpl %2             #84     a1i,A2r,ai,ar                 \n"  \
"fstpl %3             #85     A2r,ai,ar                     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7r) "         #86-89  A2R,ai,ar                     \n"  \
"fldl %7              #87     a3i,A2R,ai,ar                 \n"  \
"fsubl %5             #88-91  A2I,A2R,ai,ar                 \n"  \
"fldl %2              #89     a1r,A2I,a2r,ai,ar             \n"  \
"fadd %%st(4)         #90     NA3R,A2I,a2r,ai,ar            \n"  \
"fxch %%st(1)         #90     A2I,NA3R,a2r,ai,ar            \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7r) "         #91-94  A2i,NA3R,a2r,ai,ar            \n"  \
"fxch %%st(1)         #91     na3R,A2I,a2r,ai,ar            \n"  \
"fldl %3              #92     a1i,na3r,a2i,a2r,ai,ar        \n"  \
"fadd %%st(4)         #93-96  NA3I,na3r,a2i,a2r,ai,ar       \n"  \
"fxch %%st(1)         #93     na3r,NA3I,a2i,a2r,ai,ar       \n"  \
"fadd %%st(3)         #94     A3R,NA3I,a2i,a2r,ai,ar        \n"  \
"fldl %2              #95     a1r,A3R,NA3I,a2i,a2r,ai,ar    \n"  \
"faddl %0             #96-99  NA1R,A3R,NA3I,a2i,a2r,ai,ar   \n"  \
"fxch %%st(2)         #97     na3i,A3R,NA1R,a2i,a2r,ai,ar   \n"  \
"fadd %%st(3)         #98-101 A3I,a3r,NA1R,a2i,a2r,ai,ar    \n"  \
"fxch %%st(1)         #99     a3r,A3I,NA1R,a2i,a2r,ai,ar    \n"  \
"fstpl %6             #100    A3I,NA1R,a2i,a2r,ai,ar        \n"  \
"fstpl %7             #102    NA1R,a2i,a2r,ai,ar     STALL !\n"  \
"fldl %3              #103    a1i,na1r,a2i,a2r,ai,ar        \n"  \
"faddl %1             #104-107 Na1i,na1r,a2i,a2r,ai,ar      \n"  \
"fxch %%st(1)         #104     na1r,Na1i,a2i,a2r,ai,ar      \n"  \
"fsubr %%st(5)        #105-108 A1R,na1i,a2i,a2r,ai,ar       \n"  \
"fxch %%st(3)         #105     a2r,na1i,a2i,A1R,ai,ar       \n"  \
"fsubl %0             #106-109 NA2R,na1i,a2i,A1R,ai,ar      \n"  \
"fxch %%st(1)         #106     na1i,NA2R,a2i,A1R,ai,ar      \n"  \
"fsubr %%st(4)        #107-110 A1I,NA2R,a2i,a1r,ai,ar       \n"  \
"fxch %%st(2)         #107     a2i,NA2R,A1I,a1r,ai,ar       \n"  \
"fsubl %1             #108     NA2I,na2r,A1I,a1r,ai,ar      \n"  \
"fxch %%st(1)         #109     na2r,NA2I,a1i,a1r,ai,ar      \n"  \
"fsubrp %%st,%%st(5)  #110-113 NA2I,a1i,a1r,ai,A2R          \n"  \
"fsubrp %%st,%%st(3)  #111-114 a1i,a1r,A2I,A2R              \n"  \
"fstpl %3             #112  a1r,A2I,A2r                     \n"  \
"fstpl %2             #113  A2I,A2r                         \n"  \
"fstpl %5             #114  a2r                             \n"  \
"fstpl %4             #115                                  \n"  \
: :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),"m"(_x2##i),\
	"m"(_x3##r),"m"(_x3##i):"memory");                       \
__asm__ volatile ("                                         \n"  \
"fldl %0            #116      a4r                           \n"  \
"fsubl %2           #117-120  A45r                          \n"  \
"fldl %1            #118      a4i                           \n"  \
"fsubl %3           #119-122  A45i,A45r                     \n"  \
"fldl %4            #120      a6r,A45i,A45r                 \n"  \
"fsubl %2           #121-124  A65r,A45i,A45r                \n"  \
"fldl %5            #122      a6i,A65r,A45i,a45r            \n"  \
"fsubl %3           #123-126  a65i,A65r,a45i,a45r           \n"  \
"fxch %%st(3)       #123      a45r,A65r,a45i,A65i           \n"  \
"fsubl %4           #124-127  Nar,A65r,a45i,A65i            \n"  \
"fxch %%st(2)       #124      a45i,a65r,Nar,A65i            \n"  \
"fsubl %5           #125-128  Nai,a65r,Nar,A65i             \n"  \
"fxch %%st(1)       #125      a65r,Nai,Nar,A65i             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7i) "       #126-129  A0r,Nai,Nar,A65i              \n"  \
"fldl %4            #127      a6r,A0r,Nai,Nar,A65i          \n"  \
"fxch %%st(4)       #127      a65i,A0r,Nai,Nar,a6r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN4_7i) "       #128-131  A0i,A0r,Nai,Nar,a6r           \n"  \
"fldl %5            #129      a6i,A0i,A0r,Nai,Nar,a6r       \n"  \
"fxch %%st(4)       #130      nar,A0i,A0r,Nai,a6i,a6r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7i) "       #131-134  Ar,A0i,a0r,nai,a6i,a6r        \n"  \
"fxch %%st(5)       #131      a6r,A0i,a0r,nai,a6i,Ar        \n"  \
"faddl %0           #132-135  A6r,a0i,a0r,nai,a6i,Ar        \n"  \
"fxch %%st(3)       #132      nai,a0i,a0r,A6r,a6i,Ar        \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN1_7i) "       #133-136  Ai,a0i,a0r,a6r,a6i,Ar         \n"  \
"fxch %%st(4)       #133      a6i,a0i,a0r,a6r,Ai,ar         \n"  \
"faddl %1           #134-137  A6i,a0i,a0r,a6r,Ai,ar         \n"  \
"fldl %2            #135      a5r,A6i,a0i,a0r,a6r,Ai,ar     \n"  \
"faddl %0           #136-139  A5r,A6i,a0i,a0r,a6r,ai,ar     \n"  \
"fxch %%st(4)       #136      a6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7i) "       #137-140  A6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fxch %%st(1)       #138      A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN2_7i) "       #139-142  A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fxch %%st(4)       #140      a5r,A6r,a0i,a0r,A6i,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7i) "       #141-144  A5r,a6r,a0i,a0r,A6i,ai,ar     \n"  \
"fxch %%st(1)       #142      a6r,A5r,a0i,a0r,A6i,ai,ar     \n"  \
"fstpl %4           #143      A5r,a0i,a0r,A6i,ai,ar         \n"  \
"fldl %3            #144      a5i,A5r,a0i,a0r,a6i,ai,ar     \n"  \
"faddl %1           #145-148  A5i,a5r,a0i,a0r,a6i,ai,ar     \n"  \
"fxch %%st(4)       #146      a6i,a5r,a0i,a0r,A5i,ai,ar     \n"  \
"fstpl %5           #147      a5r,a0i,a0r,A5i,ai,ar         \n"  \
"fstpl %2           #148      a0i,a0r,A5i,ai,ar             \n"  \
"fxch %%st(2)       #148      a5i,a0r,a0i,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(BN3_7i) "       #149-152  A5i,a0r,a0i,ai,ar             \n"  \
"fld %%st(1)        #150      a0r,A5i,a0r,a0i,ai,ar         \n"  \
"fsubl %2           #151-154  Na4,A5i,a0r,a0i,ai,ar         \n"  \
"fld %%st(3)        #152      a0i,Na4r,a5i,a0r,a0i,ai,ar    \n"  \
"fsub %%st(2)       #153-156  Na4i,Na4r,a5i,a0r,a0i,ai,ar   \n"  \
"fxch %%st(1)       #153      Na4r,Na4i,a5i,a0r,a0i,ai,ar   \n"  \
"fadd %%st(6)       #154-157  A4r,Na4i,a5i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(2)       #154      a5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #155-158  A5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #155      Na4i,A5i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #156-159  A4i,A5i,A4r,a0r,a0i,ai,ar     \n"  \
"fxch %%st(2)       #157      A4r,A5i,A4i,a0r,a0i,ai,ar     \n"  \
"fstpl %0           #158      A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fsubl %5           #159-162  A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fld %%st(5)        #160      ar,A5i,A4i,a0r,a0i,ai,ar      \n"  \
"fsubl %4           #161      Na5r,a5i,a4i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #162      a5i,Na5r,a4i,a0r,a0i,ai,ar    \n"  \
"fstpl %3           #163      Na5r,a4i,a0r,a0i,ai,ar        \n"  \
"fxch %%st(1)       #164      a4i,Na5r,a0r,a0i,ai,ar        \n"  \
"fstpl %1           #165      Na5r,a0r,a0i,ai,ar            \n"  \
"faddl %2           #166      A5r,a0r,a0i,ai,ar             \n"  \
"fxch %%st(4)       #166      ar,a0r,a0i,ai,A5r             \n"  \
"fsubp %%st,%%st(1) #167-170  Na6r,a0i,ai,A5r               \n"  \
"fxch %%st(2)       #167      ai,a0i,Na6r,A5r               \n"  \
"fsubp %%st,%%st(1) #168-171  Na6i,Na6r,A5r                 \n"  \
"fxch %%st(2)       #169      A5r,Na6r,Na6i                 \n"  \
"fstpl %2           #170      Na6r,Na6i                     \n"  \
"faddl %4           #171-174  A6r,Na6i                      \n"  \
"fxch %%st(1)       #171      Na6i,A6r                      \n"  \
"faddl %5           #172      A6i,A6r                       \n"  \
  : :"m"(_x4##r),"m"(_x4##i),"m"(_x5##r),"m"(_x5##i),"m"(_x6##r),"m"(_x6##i):"memory");\
__asm__ volatile ("                                         \n"  \
"fldl %0            #173      a3r,a6i,a6r                   \n"  \
"faddl %3           #174-177  X1r,a6i,a6r                   \n"  \
"fldl %0            #175      a3r,X1r,a6i,a6r               \n"  \
"fsubl %3           #176-179  X6r,X1r,a6i,a6r               \n"  \
"fldl %1            #177      a3i,X6r,X1r,a6i,a6r           \n"  \
"fsubl %2           #178-181  X1i,X6r,x1r,a6i,a6r           \n"  \
"fldl %1            #179      a3i,X1i,x6r,x1r,a6i,a6r       \n"  \
"faddl %2           #180-183  X6i,X1i,x6r,x1r,a6i,a6r       \n"  \
"fxch %%st(3)       #181      x1r,x1i,x6r,x6i,a6i,a6r       \n"  \
"fstpl (%%eax,%%edi,8)     #182                             \n"  \
"fstpl 8(%%eax,%%edi,8)    #183                             \n"  \
"fstpl (%%edx,%%edi,8)     #184                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #185      a6r,a6i                \n"  \
  : :"m"(_x3##r),"m"(_x3##i),"m"(_x5##r),"m"(_x5##i),            \
  "D"(_d),"a"(_pd1),"d"(_pd6):"memory");                         \
__asm__ volatile ("                                         \n"  \
"fldl %0            #186      a1r,a6i,a6r                   \n"  \
"fadd %%st(1)       #187-190  X2r,a6i,a6r                   \n"  \
"fldl %0            #188      a1r,X2r,a6i,a6r               \n"  \
"fsubp %%st,%%st(2) #189-192  X2r,x5r,a6r                   \n"  \
"fldl %1            #190      a1i,X2r,x5r,a6r               \n"  \
"fsub %%st(3)       #191      X2i,X2r,x5r,a6r               \n"  \
"fldl %1            #192      a1i,X2i,x2r,x5r,a6r           \n"  \
"faddp %%st,%%st(4) #191-194  x2i,x2r,X5r,x5i               \n"  \
"fstpl 8(%%eax,%%edi,8)    #192                             \n"  \
"fstpl (%%eax,%%edi,8)     #193                             \n"  \
"fstpl (%%edx,%%edi,8)     #195                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #196                             \n"  \
  : :"m"(_x1##r),"m"(_x1##i),"D"(_d),"a"(_pd2),"d"(_pd5):        \
		  "memory");                                     \
__asm__ volatile ("                                         \n"  \
"fldl %0            #197      a2r                           \n"  \
"faddl %3           #198-201  X4r                           \n"  \
"fldl %0            #199      a2r,X4r                       \n"  \
"fsubl %3           #200-203  X3r,X4r                       \n"  \
"fldl %1            #201      a2i,X3r,X4r                   \n"  \
"fsubl %2           #202-205  X4i,X3r,x4r                   \n"  \
"fldl %1            #206      a2i,X4i,x3r,x4r               \n"  \
"faddl %2           #207-210  X3i,X4i,x3r,x4r               \n"  \
"fxch %%st(3)       #208      x4r,x4i,x3r,x3i               \n"  \
"fstpl (%%eax,%%edi,8)     #209                             \n"  \
"fstpl 8(%%eax,%%edi,8)    #210                             \n"  \
"fstpl (%%edx,%%edi,8)     #211                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #212                             \n"  \
  : :"m"(_x2##r),"m"(_x2##i),"m"(_x4##r),"m"(_x4##i),            \
  "D"(_d),"a"(_pd4),"d"(_pd3):"memory");                         \
}


# define asm_cplx_fft7_store_p_F(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3, _x4,_x5,_x6) \
{                                                                \
    __asm__ volatile("                                      \n"  \
"fldl %4	#1                                          \n"  \
"faddl %6       #2     %0                                   \n"  \
"fldl %4        #3                                          \n"  \
"fsubl %6       #4     %2,%0                                \n"  \
"fldl %5	#5                                          \n"  \
"faddl %7       #6     %1,%2,%0                             \n"  \
"fldl %5        #7                                          \n"  \
"fsubl %7       #8     %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #9     %0,%1,%2,%3                          \n"  \
"fstpl %0       #10                                         \n"  \
"fstpl %1       #11                                         \n"  \
"fstpl %2       #12                                         \n"  \
"fstpl %3       #13                                         \n"  \
  : : "m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i),          \
  "m" (_x1##r),"m"(_x1##i),"m"(_x6##r),"m"(_x6##i):"memory");    \
   __asm__ volatile("                                       \n"  \
"fldl %4	#14                                         \n"  \
"faddl %6       #15     %0                                  \n"  \
"fldl %4        #16                                         \n"  \
"fsubl %6       #17     %2,%0                               \n"  \
"fldl %5	#18                                         \n"  \
"faddl %7       #19     %1,%2,%0                            \n"  \
"fldl %5        #20                                         \n"  \
"fsubl %7       #21    %3,%1,%2,%0                          \n"  \
"fxch %%st(3)   #22     %0,%1,%2,%3                         \n"  \
"fstpl %0       #23                                         \n"  \
"fstpl %1       #24                                         \n"  \
"fstpl %2       #25                                         \n"  \
"fstpl %3       #26                                         \n"  \
: : "m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i),            \
 "m" (_x2##r),"m"(_x2##i),"m"(_x5##r),"m"(_x5##i):"memory");     \
    __asm__ volatile("                                      \n"  \
"fldl %4	#27                                         \n"  \
"faddl %6       #28    %0                                   \n"  \
"fldl %4        #29                                         \n"  \
"fsubl %6       #30     %2,%0                               \n"  \
"fldl %5	#31                                         \n"  \
"faddl %7       #32     %1,%2,%0                            \n"  \
"fldl %5        #33                                         \n"  \
"fsubl %7       #34     %3,%1,%2,%0                         \n"  \
"fxch %%st(3)   #35     %0,%1,%2,%3                         \n"  \
"fstpl %0       #36                                         \n"  \
"fstpl %1       #37                                         \n"  \
"fstpl %2       #38                                         \n"  \
"fstpl %3       #39                                         \n"  \
  : : "m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i),          \
  "m" (_x3##r),"m"(_x3##i),"m"(_x4##r),"m"(_x4##i):"memory");    \
__asm__ volatile("                                          \n"  \
"fldl %2	#40     x1r                                 \n"  \
"faddl %4       #41-44  X12R                                \n"  \
"fldl %3        #42     x1i,X12R                            \n"  \
"faddl %5       #43-46  X12I,X12R                           \n"  \
"fldl %0	#44     x0r,X12I,X12R                       \n"  \
"fldl %1        #45     x0i,x0r,X12I,x12r                   \n"  \
"fxch %%st(3)   #45     x12r,x0r,X12I,x0i                   \n"  \
"faddl %6       #46-49  AR,x0r,x12i,x0i                     \n"  \
"fxch %%st(2)   #46     x12i,x0r,AR,x0i                     \n"  \
"faddl %7       #47-50  AI,x0r,AR,x0i                       \n"  \
"fxch %%st(1)   #47     x0r,AI,AR,x0i                       \n"  \
"fstl %8        #48     x0r,AI,AR,x0i                       \n"  \
"fadd %%st(2)   #49-52  X0R,AI,AR,x0i                       \n"  \
"fxch %%st(3)   #50     x0i,ai,ar,X0R                       \n"  \
"fstl %9        #51     x0i,ai,ar,X0R                       \n"  \
"fadd %%st(1)   #52-55  X0I,ai,ar,x0r                       \n"  \
  :  :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),"m"(_x2##i),\
    "m"(_x3##r),"m"(_x3##i),"m"(_x0##r),"m"(_x0##i):"memory");   \
__asm__ volatile ("                                         \n"  \
"fxch %%st(2)          #52     ar,ai,X0I,X0R                \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7r) "        #53-56  AR,ai,X0I,x0r                \n"  \
"fxch %%st(3)          #54     x0r,ai,X0I,AR                \n"  \
"fstpl (%%eax,%%edi,8) #55     ai,X0I,AR                    \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7r) "          #56-59  AI,x0i,AR                    \n"  \
"fxch %%st(1)          #57     x0i,AI,ar                    \n"  \
"fstpl 8(%%eax,%%edi,8)#58     AI,ar                        \n"  \
  : :"D"(_d),"a"(_pd0):"memory");                                \
__asm__ volatile ("                                         \n"  \
"fxch %%st(1)         #59     ar,AI                         \n"  \
"faddl %0             #60-63  AR,ai                         \n"  \
"fxch %%st(1)         #60     ai,AR                         \n"  \
"faddl %1             #61-64  AI,AR                         \n"  \
"fldl %4	      #62     a2r,AI,AR                     \n"  \
"fsubl %2             #63-66  A0R,AI,AR                     \n"  \
"fldl %5              #64     a2i,A0R,AI,AR                 \n"  \
"fsubl %3             #65-68  A0I,A0R,AI,AR                 \n"  \
"fldl %2              #66     a1r,A0I,A0R,ai,ar             \n"  \
"fsubl %6             #67-70  A1R,A0I,A0R,ai,ar             \n"  \
"fxch %%st(2)         #68     a0r,A0I,A1R,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7r) "         #69-72  A0R,A0I,A1r,ai,ar             \n"  \
"fxch %%st(1)         #70     a0i,A0R,A1r,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7r) "         #71-74  A0I,A0R,A1r,ai,ar             \n"  \
"fldl %3              #72     a1i,A0I,A0R,A1r,ai,ar         \n"  \
"fsubl %7             #73-76  A1i,A0I,A0R,A1r,ai,ar         \n"  \
"fxch %%st(2)         #74     a0r,A0I,A1I,a1r,ai,ar         \n"  \
"fstpl %0             #75     a0i,A1I,a1r,ai,ar             \n"  \
"fstpl %1             #76     A1I,a1r,ai,ar                 \n"  \
"fxch %%st(1)         #77     a1r,a1I,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7r) "         #78-81  A1R,a1I,ai,ar                 \n"  \
"fxch %%st(1)         #79     a1i,A1R,ai,ar                 \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7r) "         #80-83  A1I,A1R,ai,ar                 \n"  \
"fldl %6              #81     a3r,A1I,A1R,ai,ar             \n"  \
"fsubl %4             #82-85  A2R,A1I,A1R,ai,ar             \n"  \
"fxch %%st(2)         #83     a1r,A1I,A2r,ai,ar             \n"  \
"fstpl %2             #84     a1i,A2r,ai,ar                 \n"  \
"fstpl %3             #85     A2r,ai,ar                     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7r) "         #86-89  A2R,ai,ar                     \n"  \
"fldl %7              #87     a3i,A2R,ai,ar                 \n"  \
"fsubl %5             #88-91  A2I,A2R,ai,ar                 \n"  \
"fldl %2              #89     a1r,A2I,a2r,ai,ar             \n"  \
"fadd %%st(4)         #90     NA3R,A2I,a2r,ai,ar            \n"  \
"fxch %%st(1)         #90     A2I,NA3R,a2r,ai,ar            \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7r) "         #91-94  A2i,NA3R,a2r,ai,ar            \n"  \
"fxch %%st(1)         #91     na3R,A2I,a2r,ai,ar            \n"  \
"fldl %3              #92     a1i,na3r,a2i,a2r,ai,ar        \n"  \
"fadd %%st(4)         #93-96  NA3I,na3r,a2i,a2r,ai,ar       \n"  \
"fxch %%st(1)         #93     na3r,NA3I,a2i,a2r,ai,ar       \n"  \
"fadd %%st(3)         #94     A3R,NA3I,a2i,a2r,ai,ar        \n"  \
"fldl %2              #95     a1r,A3R,NA3I,a2i,a2r,ai,ar    \n"  \
"faddl %0             #96-99  NA1R,A3R,NA3I,a2i,a2r,ai,ar   \n"  \
"fxch %%st(2)         #97     na3i,A3R,NA1R,a2i,a2r,ai,ar   \n"  \
"fadd %%st(3)         #98-101 A3I,a3r,NA1R,a2i,a2r,ai,ar    \n"  \
"fxch %%st(1)         #99     a3r,A3I,NA1R,a2i,a2r,ai,ar    \n"  \
"fstpl %6             #100    A3I,NA1R,a2i,a2r,ai,ar        \n"  \
"fstpl %7             #102    NA1R,a2i,a2r,ai,ar     STALL !\n"  \
"fldl %3              #103    a1i,na1r,a2i,a2r,ai,ar        \n"  \
"faddl %1             #104-107 Na1i,na1r,a2i,a2r,ai,ar      \n"  \
"fxch %%st(1)         #104     na1r,Na1i,a2i,a2r,ai,ar      \n"  \
"fsubr %%st(5)        #105-108 A1R,na1i,a2i,a2r,ai,ar       \n"  \
"fxch %%st(3)         #105     a2r,na1i,a2i,A1R,ai,ar       \n"  \
"fsubl %0             #106-109 NA2R,na1i,a2i,A1R,ai,ar      \n"  \
"fxch %%st(1)         #106     na1i,NA2R,a2i,A1R,ai,ar      \n"  \
"fsubr %%st(4)        #107-110 A1I,NA2R,a2i,a1r,ai,ar       \n"  \
"fxch %%st(2)         #107     a2i,NA2R,A1I,a1r,ai,ar       \n"  \
"fsubl %1             #108     NA2I,na2r,A1I,a1r,ai,ar      \n"  \
"fxch %%st(1)         #109     na2r,NA2I,a1i,a1r,ai,ar      \n"  \
"fsubrp %%st,%%st(5)  #110-113 NA2I,a1i,a1r,ai,A2R          \n"  \
"fsubrp %%st,%%st(3)  #111-114 a1i,a1r,A2I,A2R              \n"  \
"fstpl %3             #112  a1r,A2I,A2r                     \n"  \
"fstpl %2             #113  A2I,A2r                         \n"  \
"fstpl %5             #114  a2r                             \n"  \
"fstpl %4             #115                                  \n"  \
 : :"m"(_x0##r),"m"(_x0##i),"m"(_x1##r),"m"(_x1##i),"m"(_x2##r),"m"(_x2##i),\
	"m"(_x3##r),"m"(_x3##i):"memory");                       \
__asm__ volatile ("                                         \n"  \
"fldl %0            #116      a4r                           \n"  \
"fsubl %2           #117-120  A45r                          \n"  \
"fldl %1            #118      a4i                           \n"  \
"fsubl %3           #119-122  A45i,A45r                     \n"  \
"fldl %4            #120      a6r,A45i,A45r                 \n"  \
"fsubl %2           #121-124  A65r,A45i,A45r                \n"  \
"fldl %5            #122      a6i,A65r,A45i,a45r            \n"  \
"fsubl %3           #123-126  a65i,A65r,a45i,a45r           \n"  \
"fxch %%st(3)       #123      a45r,A65r,a45i,A65i           \n"  \
"fsubl %4           #124-127  Nar,A65r,a45i,A65i            \n"  \
"fxch %%st(2)       #124      a45i,a65r,Nar,A65i            \n"  \
"fsubl %5           #125-128  Nai,a65r,Nar,A65i             \n"  \
"fxch %%st(1)       #125      a65r,Nai,Nar,A65i             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7i) "       #126-129  A0r,Nai,Nar,A65i              \n"  \
"fldl %4            #127      a6r,A0r,Nai,Nar,A65i          \n"  \
"fxch %%st(4)       #127      a65i,A0r,Nai,Nar,a6r          \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN4_7i) "       #128-131  A0i,A0r,Nai,Nar,a6r           \n"  \
"fldl %5            #129      a6i,A0i,A0r,Nai,Nar,a6r       \n"  \
"fxch %%st(4)       #130      nar,A0i,A0r,Nai,a6i,a6r       \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7i) "       #131-134  Ar,A0i,a0r,nai,a6i,a6r        \n"  \
"fxch %%st(5)       #131      a6r,A0i,a0r,nai,a6i,Ar        \n"  \
"faddl %0           #132-135  A6r,a0i,a0r,nai,a6i,Ar        \n"  \
"fxch %%st(3)       #132      nai,a0i,a0r,A6r,a6i,Ar        \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN1_7i) "       #133-136  Ai,a0i,a0r,a6r,a6i,Ar         \n"  \
"fxch %%st(4)       #133      a6i,a0i,a0r,a6r,Ai,ar         \n"  \
"faddl %1           #134-137  A6i,a0i,a0r,a6r,Ai,ar         \n"  \
"fldl %2            #135      a5r,A6i,a0i,a0r,a6r,Ai,ar     \n"  \
"faddl %0           #136-139  A5r,A6i,a0i,a0r,a6r,ai,ar     \n"  \
"fxch %%st(4)       #136      a6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7i) "       #137-140  A6r,a6i,a0i,a0r,A5r,ai,ar     \n"  \
"fxch %%st(1)       #138      A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN2_7i) "       #139-142  A6i,A6r,a0i,a0r,a5r,ai,ar     \n"  \
"fxch %%st(4)       #140      a5r,A6r,a0i,a0r,A6i,ai,ar     \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7i) "       #141-144  A5r,a6r,a0i,a0r,A6i,ai,ar     \n"  \
"fxch %%st(1)       #142      a6r,A5r,a0i,a0r,A6i,ai,ar     \n"  \
"fstpl %4           #143      A5r,a0i,a0r,A6i,ai,ar         \n"  \
"fldl %3            #144      a5i,A5r,a0i,a0r,a6i,ai,ar     \n"  \
"faddl %1           #145-148  A5i,a5r,a0i,a0r,a6i,ai,ar     \n"  \
"fxch %%st(4)       #146      a6i,a5r,a0i,a0r,A5i,ai,ar     \n"  \
"fstpl %5           #147      a5r,a0i,a0r,A5i,ai,ar         \n"  \
"fstpl %2           #148      a0i,a0r,A5i,ai,ar             \n"  \
"fxch %%st(2)       #148      a5i,a0r,a0i,ai,ar             \n"  \
"fmull " MEMORY_VARIABLE_NAME(FN3_7i) "       #149-152  A5i,a0r,a0i,ai,ar             \n"  \
"fld %%st(1)        #150      a0r,A5i,a0r,a0i,ai,ar         \n"  \
"fsubl %2           #151-154  Na4,A5i,a0r,a0i,ai,ar         \n"  \
"fld %%st(3)        #152      a0i,Na4r,a5i,a0r,a0i,ai,ar    \n"  \
"fsub %%st(2)       #153-156  Na4i,Na4r,a5i,a0r,a0i,ai,ar   \n"  \
"fxch %%st(1)       #153      Na4r,Na4i,a5i,a0r,a0i,ai,ar   \n"  \
"fadd %%st(6)       #154-157  A4r,Na4i,a5i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(2)       #154      a5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #155-158  A5i,Na4i,A4r,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #155      Na4i,A5i,A4r,a0r,a0i,ai,ar    \n"  \
"fadd %%st(5)       #156-159  A4i,A5i,A4r,a0r,a0i,ai,ar     \n"  \
"fxch %%st(2)       #157      A4r,A5i,A4i,a0r,a0i,ai,ar     \n"  \
"fstpl %0           #158      A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fsubl %5           #159-162  A5i,A4i,a0r,a0i,ai,ar         \n"  \
"fld %%st(5)        #160      ar,A5i,A4i,a0r,a0i,ai,ar      \n"  \
"fsubl %4           #161      Na5r,a5i,a4i,a0r,a0i,ai,ar    \n"  \
"fxch %%st(1)       #162      a5i,Na5r,a4i,a0r,a0i,ai,ar    \n"  \
"fstpl %3           #163      Na5r,a4i,a0r,a0i,ai,ar        \n"  \
"fxch %%st(1)       #164      a4i,Na5r,a0r,a0i,ai,ar        \n"  \
"fstpl %1           #165      Na5r,a0r,a0i,ai,ar            \n"  \
"faddl %2           #166      A5r,a0r,a0i,ai,ar             \n"  \
"fxch %%st(4)       #166      ar,a0r,a0i,ai,A5r             \n"  \
"fsubp %%st,%%st(1) #167-170  Na6r,a0i,ai,A5r               \n"  \
"fxch %%st(2)       #167      ai,a0i,Na6r,A5r               \n"  \
"fsubp %%st,%%st(1) #168-171  Na6i,Na6r,A5r                 \n"  \
"fxch %%st(2)       #169      A5r,Na6r,Na6i                 \n"  \
"fstpl %2           #170      Na6r,Na6i                     \n"  \
"faddl %4           #171-174  A6r,Na6i                      \n"  \
"fxch %%st(1)       #171      Na6i,A6r                      \n"  \
"faddl %5           #172      A6i,A6r                       \n"  \
  : :"m"(_x4##r),"m"(_x4##i),"m"(_x5##r),"m"(_x5##i),"m"(_x6##r),"m"(_x6##i):"memory");\
__asm__ volatile ("                                         \n"  \
"fldl %0            #173      a3r,a6i,a6r                   \n"  \
"faddl %3           #174-177  X1r,a6i,a6r                   \n"  \
"fldl %0            #175      a3r,X1r,a6i,a6r               \n"  \
"fsubl %3           #176-179  X6r,X1r,a6i,a6r               \n"  \
"fldl %1            #177      a3i,X6r,X1r,a6i,a6r           \n"  \
"fsubl %2           #178-181  X1i,X6r,x1r,a6i,a6r           \n"  \
"fldl %1            #179      a3i,X1i,x6r,x1r,a6i,a6r       \n"  \
"faddl %2           #180-183  X6i,X1i,x6r,x1r,a6i,a6r       \n"  \
"fxch %%st(3)       #181      x1r,x1i,x6r,x6i,a6i,a6r       \n"  \
"fstpl (%%eax,%%edi,8)     #182                             \n"  \
"fstpl 8(%%eax,%%edi,8)    #183                             \n"  \
"fstpl (%%edx,%%edi,8)     #184                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #185      a6r,a6i                \n"  \
  : :"m"(_x3##r),"m"(_x3##i),"m"(_x5##r),"m"(_x5##i),            \
  "D"(_d),"a"(_pd1),"d"(_pd6):"memory");                         \
__asm__ volatile ("                                         \n"  \
"fldl %0            #186      a1r,a6i,a6r                   \n"  \
"fadd %%st(1)       #187-190  X2r,a6i,a6r                   \n"  \
"fldl %0            #188      a1r,X2r,a6i,a6r               \n"  \
"fsubp %%st,%%st(2) #189-192  X2r,x5r,a6r                   \n"  \
"fldl %1            #190      a1i,X2r,x5r,a6r               \n"  \
"fsub %%st(3)       #191      X2i,X2r,x5r,a6r               \n"  \
"fldl %1            #192      a1i,X2i,x2r,x5r,a6r           \n"  \
"faddp %%st,%%st(4) #191-194  x2i,x2r,X5r,x5i               \n"  \
"fstpl 8(%%eax,%%edi,8)    #192                             \n"  \
"fstpl (%%eax,%%edi,8)     #193                             \n"  \
"fstpl (%%edx,%%edi,8)     #195                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #196                             \n"  \
  : :"m"(_x1##r),"m"(_x1##i),"D"(_d),"a"(_pd2),"d"(_pd5):        \
  "memory");                                                     \
__asm__ volatile ("                                         \n"  \
"fldl %0            #197      a2r                           \n"  \
"faddl %3           #198-201  X4r                           \n"  \
"fldl %0            #199      a2r,X4r                       \n"  \
"fsubl %3           #200-203  X3r,X4r                       \n"  \
"fldl %1            #201      a2i,X3r,X4r                   \n"  \
"fsubl %2           #202-205  X4i,X3r,x4r                   \n"  \
"fldl %1            #206      a2i,X4i,x3r,x4r               \n"  \
"faddl %2           #207-210  X3i,X4i,x3r,x4r               \n"  \
"fxch %%st(3)       #208      x4r,x4i,x3r,x3i               \n"  \
"fstpl (%%eax,%%edi,8)     #209                             \n"  \
"fstpl 8(%%eax,%%edi,8)    #210                             \n"  \
"fstpl (%%edx,%%edi,8)     #211                             \n"  \
"fstpl 8(%%edx,%%edi,8)    #212                             \n"  \
  : :"m"(_x2##r),"m"(_x2##i),"m"(_x4##r),"m"(_x4##i),            \
  "D"(_d),"a"(_pd4),"d"(_pd3):"memory");                         \
}

#endif


#ifndef USE_ASM

# define cplx_fft7(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{ \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);\
            cplx_sub( _a6, _x1, _x6);\
            cplx_add( _a2, _x2, _x5);\
            cplx_sub( _a5, _x2, _x5);\
            cplx_add( _a3, _x3, _x4);\
            cplx_sub( _a4, _x3, _x4);\
            _ar = (_a1r + _a2r + _a3r);\
            _ai = (_a1i + _a2i + _a3i);\
            _x0##r += _ar;\
            _x0##i += _ai;\
	    _ar = _a0r + _ar * _S##N1_7r;\
	    _ai = _a0i + _ai * _S##N1_7r;\
	    _a0r = _S##N4_7r * ( _a2r - _a1r );\
	    _a0i = _S##N4_7r * ( _a2i - _a1i );\
	    _a1r = _S##N2_7r * ( _a1r - _a3r );\
	    _a1i = _S##N2_7r * ( _a1i - _a3i );\
	    _a2r = _S##N3_7r * ( _a3r - _a2r );\
	    _a2i = _S##N3_7r * ( _a3i - _a2i );\
	    _a3r = _ar + _a1r + _a2r;\
	    _a3i = _ai + _a1i + _a2i;\
	    _a1r = _ar - _a1r - _a0r;\
	    _a1i = _ai - _a1i - _a0i;\
	    _a2r = _ar - _a2r + _a0r;\
	    _a2i = _ai - _a2i + _a0i;\
	    \
	    _ar = _S##N1_7i * ( _a4r - _a5r - _a6r);\
	    _ai = _S##N1_7i * ( _a4i - _a5i - _a6i);\
	    _a0r = _S##N4_7i * ( _a6r - _a5r );\
	    _a0i = _S##N4_7i * ( _a6i - _a5i );\
	    _a6r = _S##N2_7i * ( _a6r + _a4r );\
	    _a6i = _S##N2_7i * ( _a6i + _a4i );\
	    _a5r = _S##N3_7i * ( _a5r + _a4r );\
	    _a5i = _S##N3_7i * ( _a5i + _a4i );\
	    _a4r = _ar - _a5r + _a0r;\
	    _a4i = _ai - _a5i + _a0i;\
	    _a5r += _ar - _a6r;\
	    _a5i += _ai - _a6i;\
	    _a6r += _ar - _a0r;\
	    _a6i += _ai - _a0i;\
	    \
	    _x1##r = _a3r + _a5i;\
            _x6##r = _a3r - _a5i;\
            _x1##i = _a3i - _a5r;\
	    _x6##i = _a3i + _a5r;\
            _x2##r = _a1r + _a6i;\
	    _x5##r = _a1r - _a6i;\
            _x2##i = _a1i - _a6r;\
            _x5##i = _a1i + _a6r;\
            _x3##r = _a2r - _a4i;\
            _x4##r = _a2r + _a4i;\
	    _x3##i = _a2i + _a4r;\
            _x4##i = _a2i - _a4r;\
	}

#else
# define cplx_fft7(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
  asm_cplx_fft7_##_S(_x0,_x1,_x2,_x3,_x4,_x5, _x6)

#endif

#ifndef USE_ASM

#define cplx_fft7_store_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{ \
	    BIG_DOUBLE  _ar, _ai;\
	    cplx_addsub( _x1, _x6);\
	    cplx_addsub( _x2, _x5);\
	    cplx_addsub( _x3, _x4);\
	    _ar = (_x1##r + _x2##r + _x3##r);\
	    _ai = (_x1##i + _x2##i + _x3##i);\
	    *(_pd0 + _index) = _x0##r +  _ar;\
	    *(_pd0 + _index + 1) = _x0##i + _ai;\
	    _ar = _x0##r + _ar * _S##N1_7r;\
	    _ai = _x0##i + _ai * _S##N1_7r;\
	    _x0##r = _S##N4_7r * ( _x2##r - _x1##r );\
	    _x0##i = _S##N4_7r * ( _x2##i - _x1##i );\
	    _x1##r = _S##N2_7r * ( _x1##r - _x3##r );\
	    _x1##i = _S##N2_7r * ( _x1##i - _x3##i );\
	    _x2##r = _S##N3_7r * ( _x3##r - _x2##r );\
	    _x2##i = _S##N3_7r * ( _x3##i - _x2##i );\
	    _x3##r = _ar + _x1##r + _x2##r;\
	    _x3##i = _ai + _x1##i + _x2##i;\
	    _x1##r = _ar - _x1##r - _x0##r;\
	    _x1##i = _ai - _x1##i - _x0##i;\
	    _x2##r = _ar - _x2##r + _x0##r;\
	    _x2##i = _ai - _x2##i + _x0##i;\
	    \
	    _ar = _S##N1_7i * ( _x4##r - _x5##r - _x6##r);\
	    _ai = _S##N1_7i * ( _x4##i - _x5##i - _x6##i);\
	    _x0##r = _S##N4_7i * ( _x6##r - _x5##r );\
	    _x0##i = _S##N4_7i * ( _x6##i - _x5##i );\
	    _x6##r = _S##N2_7i * ( _x6##r + _x4##r );\
	    _x6##i = _S##N2_7i * ( _x6##i + _x4##i );\
	    _x5##r = _S##N3_7i * ( _x5##r + _x4##r );\
	    _x5##i = _S##N3_7i * ( _x5##i + _x4##i );\
	    _x4##r = _ar - _x5##r + _x0##r;\
	    _x4##i = _ai - _x5##i + _x0##i;\
	    _x5##r += _ar - _x6##r;\
	    _x5##i += _ai - _x6##i;\
	    _x6##r += _ar - _x0##r;\
	    _x6##i += _ai - _x0##i;\
	    \
	    *( _pd1 + _index ) = _x3##r + _x5##i;\
	    *( _pd6 + _index ) = _x3##r - _x5##i;\
	    *( _pd1 + _index + 1) = _x3##i - _x5##r;\
	    *( _pd6 + _index + 1) = _x3##i + _x5##r;\
	    *( _pd2 + _index) = _x1##r + _x6##i;\
	    *( _pd5 + _index) = _x1##r - _x6##i;\
	    *( _pd2 + _index + 1 ) = _x1##i - _x6##r;\
	    *( _pd5 + _index + 1 ) = _x1##i + _x6##r;\
	    *( _pd3 + _index) = _x2##r - _x4##i;\
	    *( _pd4 + _index) = _x2##r + _x4##i;\
	    *( _pd3 + _index + 1 ) = _x2##i + _x4##r;\
	    *( _pd4 + _index + 1 ) = _x2##i - _x4##r;\
	}
#else
#define cplx_fft7_store_p(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
 asm_cplx_fft7_store_p_##_S(_d,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

#endif


/*$Id$*/






