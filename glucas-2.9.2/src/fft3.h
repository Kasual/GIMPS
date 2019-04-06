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

/* macros to do a FFT tarnsform length 3*/
#if (Y_TARGET == 1 || Y_TARGET == 11 || Y_TARGET == 12 || Y_TARGET == 13)
# include "fft3x86.h"
#else

# define cplx_fft3(_S, _r0, _r1, _r2) cplx_fft3_##_S(_r0,_r1,_r2)

# define cplx_fft3_store_p( _index, _pd0, _pd1, _pd2, _S, _r0, _r1, _r2) \
   cplx_fft3_store_p_##_S(_index,_pd0,_pd1,_pd2,_r0,_r1,_r2)

/* These are the Nussbaumer version. Taken from E. Mayer (Mlucas code) */


# define cplx_fft3_F(_r0, _r1, _r2)       \
	{                                 \
	   y_limb_t _ar,_ai;              \
	   cplx_addsub( _r1, _r2);        \
	   _r0##r += _r1##r;              \
	   _r0##i += _r1##i;              \
	   _ar = _r2##r * F_1_3i;         \
           _ai = _r2##i * F_1_3i;         \
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   _r2##r = _r1##r + _ai;         \
	   _r2##i = _r1##i - _ar;         \
           _r1##r = _r1##r - _ai;         \
           _r1##i = _r1##i + _ar;         \
        }


# define cplx_fft3_B(_r0, _r1, _r2)       \
	{                                 \
	   y_limb_t _ar,_ai;              \
	   cplx_addsub( _r1, _r2);        \
	   _r0##r += _r1##r;              \
	   _r0##i += _r1##i;              \
	   _ar = _r2##r * F_1_3i;         \
           _ai = _r2##i * F_1_3i;         \
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   _r2##r = _r1##r - _ai;         \
	   _r2##i = _r1##i + _ar;         \
           _r1##r = _r1##r + _ai;         \
           _r1##i = _r1##i - _ar;         \
        }


# define cplx_fft3_store_p_F( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2) \
	{ \
	   cplx_addsub( _r1, _r2);\
	   _r0##r += _r1##r;\
	   _r0##i += _r1##i;\
	   _r2##r *= F_1_3i;\
	   _r2##i *= F_1_3i;\
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   *( _pd0 + _index) = _r0##r;\
	   *( _pd0 + _index + 1) = _r0##i;\
	   *( _pd2 + _index) = _r1##r + _r2##i;\
	   *( _pd2 + _index + 1) = _r1##i - _r2##r;\
	   *( _pd1 + _index) = _r1##r - _r2##i;\
	   *( _pd1 + _index+ 1) = _r1##i + _r2##r;\
	}



# define cplx_fft3_store_p_B( _index, _pd0, _pd1, _pd2, _r0, _r1, _r2) \
	{ \
	   cplx_addsub( _r1, _r2);\
	   _r0##r += _r1##r;\
	   _r0##i += _r1##i;\
	   _r2##r *= F_1_3i;\
	   _r2##i *= F_1_3i;\
	   _r1##r = _r0##r - 1.5 * _r1##r;\
	   _r1##i = _r0##i - 1.5 * _r1##i;\
	   *( _pd0 + _index) = _r0##r;\
	   *( _pd0 + _index + 1) = _r0##i;\
	   *( _pd2 + _index) = _r1##r - _r2##i;\
	   *( _pd2 + _index + 1) = _r1##i + _r2##r;\
	   *( _pd1 + _index) = _r1##r + _r2##i;\
	   *( _pd1 + _index+ 1) = _r1##i - _r2##r;\
	}



#endif

# define cplx_two_fft3(_S,_r0,_r1,_r2,_s0,_s1,_s2)    \
 cplx_two_fft3_##_S(_r0,_r1,_r2,_s0,_s1,_s2)


# define cplx_three_fft3(_S,_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2) \
 cplx_three_fft3_##_S(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2)

# define cplx_two_fft3_F(_r0,_r1,_r2,_s0,_s1,_s2)     \
	{                                             \
	   y_limb_t _ar,_ai;                          \
	   y_limb_t _br,_bi;                          \
	   cplx_addsub( _r1, _r2);                    \
	   cplx_addsub( _s1, _s2);                    \
	   _r0##r += _r1##r;                          \
	   _r0##i += _r1##i;                          \
	   _s0##r += _s1##r;                          \
	   _s0##i += _s1##i;                          \
	   _ar = _r2##r * F_1_3i;                     \
	   _ai = _r2##i * F_1_3i;                     \
	   _br = _s2##r * F_1_3i;                     \
	   _bi = _s2##i * F_1_3i;                     \
	   _r1##r = _r0##r - 1.5 * _r1##r;            \
	   _r1##i = _r0##i - 1.5 * _r1##i;            \
	   _s1##r = _s0##r - 1.5 * _s1##r;            \
	   _s1##i = _s0##i - 1.5 * _s1##i;            \
	   _r2##r = _r1##r + _ai;                     \
	   _r2##i = _r1##i - _ar;                     \
	   _s2##r = _s1##r + _bi;                     \
	   _s2##i = _s1##i - _br;                     \
           _r1##r = _r1##r - _ai;                     \
           _r1##i = _r1##i + _ar;                     \
           _s1##r = _s1##r - _bi;                     \
           _s1##i = _s1##i + _br;                     \
        }

# define cplx_two_fft3_B(_r0,_r1,_r2,_s0,_s1,_s2)     \
	{                                             \
	   y_limb_t _ar,_ai;                          \
	   y_limb_t _br,_bi;                          \
	   cplx_addsub( _r1, _r2);                    \
	   cplx_addsub( _s1, _s2);                    \
	   _r0##r += _r1##r;                          \
	   _r0##i += _r1##i;                          \
	   _s0##r += _s1##r;                          \
	   _s0##i += _s1##i;                          \
	   _ar = _r2##r * F_1_3i;                     \
	   _ai = _r2##i * F_1_3i;                     \
	   _br = _s2##r * F_1_3i;                     \
	   _bi = _s2##i * F_1_3i;                     \
	   _r1##r = _r0##r - 1.5 * _r1##r;            \
	   _r1##i = _r0##i - 1.5 * _r1##i;            \
	   _s1##r = _s0##r - 1.5 * _s1##r;            \
	   _s1##i = _s0##i - 1.5 * _s1##i;            \
	   _r2##r = _r1##r - _ai;                     \
	   _r2##i = _r1##i + _ar;                     \
	   _s2##r = _s1##r - _bi;                     \
	   _s2##i = _s1##i + _br;                     \
           _r1##r = _r1##r + _ai;                     \
           _r1##i = _r1##i - _ar;                     \
           _s1##r = _s1##r + _bi;                     \
           _s1##i = _s1##i - _br;                     \
        }


# define cplx_three_fft3_F(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2) \
	{                                             \
	   y_limb_t _ar,_ai;                          \
	   y_limb_t _br,_bi;                          \
	   y_limb_t _cr,_ci;                          \
	   cplx_addsub( _r1, _r2);                    \
	   cplx_addsub( _s1, _s2);                    \
	   cplx_addsub( _t1, _t2);                    \
	   _r0##r += _r1##r;                          \
	   _r0##i += _r1##i;                          \
	   _s0##r += _s1##r;                          \
	   _s0##i += _s1##i;                          \
	   _t0##r += _t1##r;                          \
	   _t0##i += _t1##i;                          \
	   _ar = _r2##r * F_1_3i;                     \
	   _ai = _r2##i * F_1_3i;                     \
	   _br = _s2##r * F_1_3i;                     \
	   _bi = _s2##i * F_1_3i;                     \
	   _cr = _t2##r * F_1_3i;                     \
	   _ci = _t2##i * F_1_3i;                     \
	   _r1##r = _r0##r - 1.5 * _r1##r;            \
	   _r1##i = _r0##i - 1.5 * _r1##i;            \
	   _s1##r = _s0##r - 1.5 * _s1##r;            \
	   _s1##i = _s0##i - 1.5 * _s1##i;            \
	   _t1##r = _t0##r - 1.5 * _t1##r;            \
	   _t1##i = _t0##i - 1.5 * _t1##i;            \
	   _r2##r = _r1##r + _ai;                     \
	   _r2##i = _r1##i - _ar;                     \
	   _s2##r = _s1##r + _bi;                     \
	   _s2##i = _s1##i - _br;                     \
	   _t2##r = _t1##r + _ci;                     \
	   _t2##i = _t1##i - _cr;                     \
           _r1##r = _r1##r - _ai;                     \
           _r1##i = _r1##i + _ar;                     \
           _s1##r = _s1##r - _bi;                     \
           _s1##i = _s1##i + _br;                     \
           _t1##r = _t1##r - _ci;                     \
           _t1##i = _t1##i + _cr;                     \
        }


# define cplx_three_fft3_B(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2) \
	{                                             \
	   y_limb_t _ar,_ai;                          \
	   y_limb_t _br,_bi;                          \
	   y_limb_t _cr,_ci;                          \
	   cplx_addsub( _r1, _r2);                    \
	   cplx_addsub( _s1, _s2);                    \
	   cplx_addsub( _t1, _t2);                    \
	   _r0##r += _r1##r;                          \
	   _r0##i += _r1##i;                          \
	   _s0##r += _s1##r;                          \
	   _s0##i += _s1##i;                          \
	   _t0##r += _t1##r;                          \
	   _t0##i += _t1##i;                          \
	   _ar = _r2##r * F_1_3i;                     \
	   _ai = _r2##i * F_1_3i;                     \
	   _br = _s2##r * F_1_3i;                     \
	   _bi = _s2##i * F_1_3i;                     \
	   _cr = _t2##r * F_1_3i;                     \
	   _ci = _t2##i * F_1_3i;                     \
	   _r1##r = _r0##r - 1.5 * _r1##r;            \
	   _r1##i = _r0##i - 1.5 * _r1##i;            \
	   _s1##r = _s0##r - 1.5 * _s1##r;            \
	   _s1##i = _s0##i - 1.5 * _s1##i;            \
	   _t1##r = _t0##r - 1.5 * _t1##r;            \
	   _t1##i = _t0##i - 1.5 * _t1##i;            \
	   _r2##r = _r1##r - _ai;                     \
	   _r2##i = _r1##i + _ar;                     \
	   _s2##r = _s1##r - _bi;                     \
	   _s2##i = _s1##i + _br;                     \
	   _t2##r = _t1##r - _ci;                     \
	   _t2##i = _t1##i + _cr;                     \
           _r1##r = _r1##r + _ai;                     \
           _r1##i = _r1##i - _ar;                     \
           _s1##r = _s1##r + _bi;                     \
           _s1##i = _s1##i - _br;                     \
           _t1##r = _t1##r + _ci;                     \
           _t1##i = _t1##i - _cr;                     \
        }

# define cplx_two_fft3_store_p_F(_index,_pd0,_pd1,_pd2,_qd0,_qd1,_qd2,_r0,_r1,_r2,_s0,_s1,_s2)     \
	{                                             \
	   y_limb_t _ar,_ai;                          \
	   y_limb_t _br,_bi;                          \
	   cplx_addsub( _r1, _r2);                    \
	   cplx_addsub( _s1, _s2);                    \
           _r0##r += _r1##r;                          \
	   _r0##i += _r1##i;                          \
	   _s0##r += _s1##r;                          \
	   _s0##i += _s1##i;                          \
	   _r1##r = _r0##r - 1.5 * _r1##r;            \
	   _r1##i = _r0##i - 1.5 * _r1##i;            \
	   _s1##r = _s0##r - 1.5 * _s1##r;            \
	   _s1##i = _s0##i - 1.5 * _s1##i;            \
	   _ar = _r2##r * F_1_3i;                     \
	   _ai = _r2##i * F_1_3i;                     \
	   _br = _s2##r * F_1_3i;                     \
	   _bi = _s2##i * F_1_3i;                     \
           *(_pd0 + _index) = _r0##r;                 \
           *(_pd0 + _index + 1) = _r0##i;             \
           *(_qd0 + _index) = _s0##r;                 \
           *(_qd0 + _index + 1) = _s0##i;             \
	   *(_pd2 + _index) = _r1##r + _ai;           \
           *(_pd2 + _index + 1) = _r1##i - _ar;       \
	   *(_qd2 + _index) = _s1##r + _bi;           \
	   *(_qd2 + _index + 1) = _s1##i - _br;       \
           *(_pd1 + _index) = _r1##r - _ai;           \
           *(_pd1 + _index + 1) = _r1##i + _ar;       \
           *(_qd1 + _index) = _s1##r - _bi;           \
           *(_qd1 + _index + 1) = _s1##i + _br;       \
        }


/*$Id$*/







