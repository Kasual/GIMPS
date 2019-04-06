/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006 Guillermo Ballester Valor
 
    This is the SSE2 version of fft3.h file
 
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

/* macros to do a FFT tarnsform length 3*/

#define sse2_fft3(_S, _r0, _r1, _r2) sse2_fft3_##_S(_r0, _r1, _r2)

#define sse2_fft3_store_p( _index, _pd0, _pd1, _pd2, _S, _r0, _r1, _r2) \
   sse2_fft3_store_p_##_S(_index, _pd0, _pd1, _pd2, _r0, _r1, _r2)

#define sse2_fft3_store_inter_p( _index, _pd0, _pd1, _pd2, _S, _r0, _r1, _r2) \
   sse2_fft3_store_inter_p_##_S(_index, _pd0, _pd1, _pd2, _r0, _r1, _r2)

/* These are the Nussbaumer version. Taken from E. Mayer (Mlucas code) */


#define sse2_fft3_F(_r0, _r1, _r2)                \
	{                                         \
	   Y__M128D _ar,_ai;                      \
	   sse2_addsub( _r1, _r2);                \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);  \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);\
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);  \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);\
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);  \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);  \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);  \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);  \
           Y_MM_ADD_PD( _r2##r, _r1##r, _ai);     \
           Y_MM_SUB_PD( _r2##i, _r1##i, _ar);     \
           Y_MM_SUB_PD( _r1##r, _r1##r, _ai);     \
           Y_MM_ADD_PD( _r1##i, _r1##i, _ar);     \
        }

#define sse2_fft3_B(_r0, _r1, _r2)                \
	{                                         \
	   Y__M128D _ar,_ai;                      \
	   sse2_addsub( _r1, _r2);                \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);  \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);\
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);  \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);\
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);  \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);  \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);  \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);  \
           Y_MM_SUB_PD( _r2##r, _r1##r, _ai);     \
           Y_MM_ADD_PD( _r2##i, _r1##i, _ar);     \
           Y_MM_ADD_PD( _r1##r, _r1##r, _ai);     \
           Y_MM_SUB_PD( _r1##i, _r1##i, _ar);     \
        }



#define sse2_fft3_store_p_F( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2)      \
	{                                                             \
	   Y__M128D _ar,_ai;                                          \
	   sse2_addsub( _r1, _r2);                                    \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);                      \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);                    \
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);                      \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);                    \
           Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _r0##r, _r0##i);      \
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);                      \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);                      \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);                      \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);                      \
           Y_MM_ADD_PD( _r2##r, _r1##r, _ai);                         \
           Y_MM_SUB_PD( _r2##i, _r1##i, _ar);                         \
           Y_MM_SUB_PD( _r1##r, _r1##r, _ai);                         \
           Y_MM_ADD_PD( _r1##i, _r1##i, _ar);                         \
           Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _r2##r, _r2##i);      \
           Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _r1##r, _r1##i);      \
        }

#define sse2_fft3_store_p_B( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2) \
	{                                                             \
	   Y__M128D _ar,_ai;                                          \
	   sse2_addsub( _r1, _r2);                                    \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);                      \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);                    \
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);                      \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);                    \
           Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _r0##r, _r0##i);      \
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);                      \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);                      \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);                      \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);                      \
           Y_MM_SUB_PD( _r2##r, _r1##r, _ai);                         \
           Y_MM_ADD_PD( _r2##i, _r1##i, _ar);                         \
           Y_MM_ADD_PD( _r1##r, _r1##r, _ai);                         \
           Y_MM_SUB_PD( _r1##i, _r1##i, _ar);                         \
           Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _r2##r, _r2##i);      \
           Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _r1##r, _r1##i);      \
        }

#define sse2_fft3_store_inter_p_F( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2)      \
	{                                                                   \
	   Y__M128D _ar,_ai;                                                \
	   sse2_addsub( _r1, _r2);                                          \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);                            \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);                          \
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);                            \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);                          \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _r0##r, _r0##i);      \
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);                            \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);                            \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);                            \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);                            \
           Y_MM_ADD_PD( _r2##r, _r1##r, _ai);                               \
           Y_MM_SUB_PD( _r2##i, _r1##i, _ar);                               \
           Y_MM_SUB_PD( _r1##r, _r1##r, _ai);                               \
           Y_MM_ADD_PD( _r1##i, _r1##i, _ar);                               \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _r2##r, _r2##i);      \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _r1##r, _r1##i);      \
        }

#define sse2_fft3_store_inter_p_B( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2) \
	{                                                                   \
	   Y__M128D _ar,_ai;                                                \
	   sse2_addsub( _r1, _r2);                                          \
           Y_MM_ADD_PD( _r0##r, _r0##r, _r1##r);                            \
           Y_MM_MUL_PD( _r1##r, _r1##r, MM_FM150);                          \
           Y_MM_ADD_PD( _r0##i, _r0##i, _r1##i);                            \
           Y_MM_MUL_PD( _r1##i, _r1##i, MM_FM150);                          \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _r0##r, _r0##i);      \
           Y_MM_MUL_PD( _ar, _r2##r, MM_F_1_3i);                            \
           Y_MM_ADD_PD( _r1##r, _r1##r, _r0##r);                            \
           Y_MM_MUL_PD( _ai, _r2##i, MM_F_1_3i);                            \
           Y_MM_ADD_PD( _r1##i, _r1##i, _r0##i);                            \
           Y_MM_SUB_PD( _r2##r, _r1##r, _ai);                               \
           Y_MM_ADD_PD( _r2##i, _r1##i, _ar);                               \
           Y_MM_ADD_PD( _r1##r, _r1##r, _ai);                               \
           Y_MM_SUB_PD( _r1##i, _r1##i, _ar);                               \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _r2##r, _r2##i);      \
           Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _r1##r, _r1##i);      \
        }

/*$Id$*/
