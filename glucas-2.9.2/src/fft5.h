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

/* Nussbaumer version for FFT-5. Code adapted from E. Mayer's (Mlucas) */
#if (Y_TARGET == 1 || Y_TARGET == 11 || Y_TARGET == 12 || Y_TARGET == 13)
# include "fft5x86.h"
#else

# if defined(__GNUC__)
#  define  cplx_fft5_F(_x0, _x1, _x2, _x3, _x4) \
{                                            \
     y_limb_t _ar, _ai, _br, _bi;            \
     cplx_add(_a,_x1,_x4);                   \
     cplx_add(_b,_x2,_x3);                   \
     cplx_sub(_x4,_x1,_x4);                  \
     cplx_sub(_x3,_x2,_x3);                  \
     cplx_sub(_x2,_a,_b);                    \
     cplx_add(_a,_a,_b);                     \
     cplx_sub(_b,_x4,_x3);                   \
     _x2##r *= FN2_5r;                       \
     _x2##i *= FN2_5r;                       \
     cplx_add(_x0,_x0,_a);                   \
     _br *= F_1_5i;                          \
     _bi *= F_1_5i;                          \
     _ar = _x0##r - 1.25 * _ar;              \
     _ai = _x0##i - 1.25 * _ai;              \
     _x3##r *= FN1_5i;                       \
     _x3##i *= FN1_5i;                       \
     _x4##r *= FN2_5i;                       \
     _x4##i *= FN2_5i;                       \
     _x1##r = _ar + _x2##r;                  \
     _x1##i = _ai + _x2##i;                  \
     _x2##r = _ar - _x2##r;                  \
     _x2##i = _ai - _x2##i;                  \
     _x3##r += _br;                          \
     _x3##i += _bi;                          \
     _x4##r = _br - _x4##r;                  \
     _x4##i = _bi - _x4##i;                  \
                                             \
     _ar = _x4##r;                           \
     _ai = _x4##i;                           \
     _x4##r = _x1##r + _x3##i;               \
     _x4##i = _x1##i - _x3##r;               \
     _x1##r -= _x3##i;                       \
     _x1##i += _x3##r;                       \
     _x3##r = _x2##r + _ai;                  \
     _x3##i = _x2##i - _ar;                  \
     _x2##r -= _ai;                          \
     _x2##i += _ar;                          \
}

#  define cplx_fft5_B(_x0, _x1, _x2, _x3, _x4) \
{                                            \
     y_limb_t _ar, _ai, _br, _bi;            \
     cplx_add(_a,_x1,_x4);                   \
     cplx_add(_b,_x2,_x3);                   \
     cplx_sub(_x4,_x1,_x4);                  \
     cplx_sub(_x3,_x2,_x3);                  \
     cplx_sub(_x2,_a,_b);                    \
     cplx_add(_a,_a,_b);                     \
     cplx_sub(_b,_x4,_x3);                   \
     _x2##r *= FN2_5r;                       \
     _x2##i *= FN2_5r;                       \
     cplx_add(_x0,_x0,_a);                   \
     _br *= F_1_5i;                          \
     _bi *= F_1_5i;                          \
     _ar = _x0##r - 1.25 * _ar;              \
     _ai = _x0##i - 1.25 * _ai;              \
     _x3##r *= FN1_5i;                       \
     _x3##i *= FN1_5i;                       \
     _x4##r *= FN2_5i;                       \
     _x4##i *= FN2_5i;                       \
     _x1##r = _ar + _x2##r;                  \
     _x1##i = _ai + _x2##i;                  \
     _x2##r = _ar - _x2##r;                  \
     _x2##i = _ai - _x2##i;                  \
     _x3##r += _br;                          \
     _x3##i += _bi;                          \
     _x4##r = _br - _x4##r;                  \
     _x4##i = _bi - _x4##i;                  \
                                             \
     _ar = _x4##r;                           \
     _ai = _x4##i;                           \
     _x4##r = _x1##r - _x3##i;               \
     _x4##i = _x1##i + _x3##r;               \
     _x1##r += _x3##i;                       \
     _x1##i -= _x3##r;                       \
     _x3##r = _x2##r - _ai;                  \
     _x3##i = _x2##i + _ar;                  \
     _x2##r += _ai;                          \
     _x2##i -= _ar;                          \
}

#  define cplx_fft5_store_p_F(_index,_pd0,_pd1,_pd2,_pd3,_pd4, _x0, _x1, _x2, _x3, _x4) \
{                                            \
     y_limb_t _ar, _ai, _br, _bi;            \
     cplx_add(_a,_x1,_x4);                   \
     cplx_add(_b,_x2,_x3);                   \
     cplx_sub(_x4,_x1,_x4);                  \
     cplx_sub(_x3,_x2,_x3);                  \
     cplx_sub(_x2,_a,_b);                    \
     cplx_add(_a,_a,_b);                     \
     cplx_sub(_b,_x4,_x3);                   \
     _x2##r *= FN2_5r;                       \
     _x2##i *= FN2_5r;                       \
     cplx_add(_x0,_x0,_a);                   \
     _br *= F_1_5i;                          \
     _bi *= F_1_5i;                          \
     _ar = _x0##r - 1.25 * _ar;              \
     _ai = _x0##i - 1.25 * _ai;              \
     _x3##r *= FN1_5i;                       \
     _x3##i *= FN1_5i;                       \
     _x4##r *= FN2_5i;                       \
     _x4##i *= FN2_5i;                       \
     _x1##r = _ar + _x2##r;                  \
     _x1##i = _ai + _x2##i;                  \
     _x2##r = _ar - _x2##r;                  \
     _x2##i = _ai - _x2##i;                  \
     _x3##r += _br;                          \
     _x3##i += _bi;                          \
     _x4##r = _br - _x4##r;                  \
     _x4##i = _bi - _x4##i;                  \
     *( _pd0 + _index) = _x0##r;             \
     *( _pd0 + _index + 1) = _x0##i;         \
     *( _pd4 + _index) = _x1##r + _x3##i;    \
     *( _pd4 + _index + 1) = _x1##i - _x3##r;\
     *( _pd1 + _index ) = _x1##r - _x3##i;   \
     *( _pd1 + _index + 1 ) = _x1##i + _x3##r;\
     *( _pd3 + _index ) = _x2##r + _x4##i;   \
     *( _pd3 + _index + 1) = _x2##i - _x4##r;\
     *( _pd2 + _index) = _x2##r - _x4##i;    \
     *( _pd2 + _index + 1) = _x2##i + _x4##r;\
}

#  define cplx_fft5_store_p_B(_index,_pd0,_pd1,_pd2,_pd3,_pd4, _x0, _x1, _x2, _x3, _x4) \
{                                            \
     y_limb_t _ar, _ai, _br, _bi;            \
     cplx_add(_a,_x1,_x4);                   \
     cplx_add(_b,_x2,_x3);                   \
     cplx_sub(_x4,_x1,_x4);                  \
     cplx_sub(_x3,_x2,_x3);                  \
     cplx_sub(_x2,_a,_b);                    \
     cplx_add(_a,_a,_b);                     \
     cplx_sub(_b,_x4,_x3);                   \
     _x2##r *= FN2_5r;                       \
     _x2##i *= FN2_5r;                       \
     cplx_add(_x0,_x0,_a);                   \
     _br *= F_1_5i;                          \
     _bi *= F_1_5i;                          \
     _ar = _x0##r - 1.25 * _ar;              \
     _ai = _x0##i - 1.25 * _ai;              \
     _x3##r *= FN1_5i;                       \
     _x3##i *= FN1_5i;                       \
     _x4##r *= FN2_5i;                       \
     _x4##i *= FN2_5i;                       \
     _x1##r = _ar + _x2##r;                  \
     _x1##i = _ai + _x2##i;                  \
     _x2##r = _ar - _x2##r;                  \
     _x2##i = _ai - _x2##i;                  \
     _x3##r += _br;                          \
     _x3##i += _bi;                          \
     _x4##r = _br - _x4##r;                  \
     _x4##i = _bi - _x4##i;                  \
     *( _pd0 + _index) = _x0##r;             \
     *( _pd0 + _index + 1) = _x0##i;         \
     *( _pd4 + _index) = _x1##r - _x3##i;    \
     *( _pd4 + _index + 1) = _x1##i + _x3##r;\
     *( _pd1 + _index ) = _x1##r + _x3##i;   \
     *( _pd1 + _index + 1 ) = _x1##i - _x3##r;\
     *( _pd3 + _index ) = _x2##r - _x4##i;   \
     *( _pd3 + _index + 1) = _x2##i + _x4##r;\
     *( _pd2 + _index) = _x2##r + _x4##i;    \
     *( _pd2 + _index + 1) = _x2##i - _x4##r;\
}
# else

/* for non GNU-Gcc compiler */

#  define cplx_fft5_F(_x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = FN2_5r * ( _x1##r - _x2##r);\
     _x2##i = FN2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = F_1_5i * (_x4##r - _x3##r);\
     _ai = F_1_5i * (_x4##i - _x3##i);\
     _x3##r *= FN1_5i;\
     _x3##i *= FN1_5i;\
     _x4##r *= FN2_5i;\
     _x4##i *= FN2_5i;\
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

#  define cplx_fft5_B(_x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = FN2_5r * ( _x1##r - _x2##r);\
     _x2##i = FN2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = F_1_5i * (_x4##r - _x3##r);\
     _ai = F_1_5i * (_x4##i - _x3##i);\
     _x3##r *= FN1_5i;\
     _x3##i *= FN1_5i;\
     _x4##r *= FN2_5i;\
     _x4##i *= FN2_5i;\
     _x3##r += _ar;\
     _x3##i += _ai;\
     _x4##r = _ar - _x4##r;\
     _x4##i = _ai - _x4##i;\
     \
     _ar = _x4##r;\
     _ai = _x4##i;\
     _x4##r = _x1##r - _x3##i;\
     _x4##i = _x1##i + _x3##r;\
     _x1##r += _x3##i;\
     _x1##i -= _x3##r;\
     _x3##r = _x2##r - _ai;\
     _x3##i = _x2##i + _ar;\
     _x2##r += _ai;\
     _x2##i -= _ar;\
}

#  define cplx_fft5_store_p_F(_index,_pd0,_pd1,_pd2,_pd3,_pd4, _x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = FN2_5r * ( _x1##r - _x2##r);\
     _x2##i = FN2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     *( _pd0 + _index) = _x0##r;\
     *( _pd0 + _index + 1) = _x0##i;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = F_1_5i * (_x4##r - _x3##r);\
     _ai = F_1_5i * (_x4##i - _x3##i);\
     _x3##r *= FN1_5i;\
     _x3##i *= FN1_5i;\
     _x4##r *= FN2_5i;\
     _x4##i *= FN2_5i;\
     _x3##r += _ar;\
     _x3##i += _ai;\
     _x4##r = _ar - _x4##r;\
     _x4##i = _ai - _x4##i;\
     *( _pd4 + _index) = _x1##r + _x3##i;\
     *( _pd4 + _index + 1) = _x1##i - _x3##r;\
     *( _pd1 + _index ) = _x1##r - _x3##i;\
     *( _pd1 + _index + 1 ) = _x1##i + _x3##r;\
     *( _pd3 + _index ) = _x2##r + _x4##i;\
     *( _pd3 + _index + 1) = _x2##i - _x4##r;\
     *( _pd2 + _index) = _x2##r - _x4##i;\
     *( _pd2 + _index + 1) = _x2##i + _x4##r;\
}

#  define cplx_fft5_store_p_B(_index,_pd0,_pd1,_pd2,_pd3,_pd4, _x0, _x1, _x2, _x3, _x4) \
{ \
     y_limb_t _ar, _ai;\
     cplx_addsub( _x1, _x4);\
     cplx_addsub( _x2, _x3);\
     _ar = _x1##r + _x2##r;\
     _ai = _x1##i + _x2##i;\
     _x0##r += _ar;\
     _x0##i += _ai;\
     _x2##r = FN2_5r * ( _x1##r - _x2##r);\
     _x2##i = FN2_5r * ( _x1##i - _x2##i);\
     _ar = _x0##r - 1.25 * _ar;\
     _ai = _x0##i - 1.25 * _ai;\
     *( _pd0 + _index) = _x0##r;\
     *( _pd0 + _index + 1) = _x0##i;\
     _x1##r = _ar + _x2##r;\
     _x1##i = _ai + _x2##i;\
     _x2##r = _ar - _x2##r;\
     _x2##i = _ai - _x2##i;\
     _ar = F_1_5i * (_x4##r - _x3##r);\
     _ai = F_1_5i * (_x4##i - _x3##i);\
     _x3##r *= FN1_5i;\
     _x3##i *= FN1_5i;\
     _x4##r *= FN2_5i;\
     _x4##i *= FN2_5i;\
     _x3##r += _ar;\
     _x3##i += _ai;\
     _x4##r = _ar - _x4##r;\
     _x4##i = _ai - _x4##i;\
     *( _pd4 + _index) = _x1##r - _x3##i;\
     *( _pd4 + _index + 1) = _x1##i + _x3##r;\
     *( _pd1 + _index ) = _x1##r + _x3##i;\
     *( _pd1 + _index + 1 ) = _x1##i - _x3##r;\
     *( _pd3 + _index ) = _x2##r - _x4##i;\
     *( _pd3 + _index + 1) = _x2##i + _x4##r;\
     *( _pd2 + _index) = _x2##r + _x4##i;\
     *( _pd2 + _index + 1) = _x2##i - _x4##r;\
}

# endif


# define cplx_fft5(_S, _x0, _x1, _x2, _x3, _x4) \
    cplx_fft5_##_S(_x0,_x1,_x2,_x3,_x4);

# define cplx_fft5_store_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_S, _x0, _x1, _x2, _x3, _x4) \
  cplx_fft5_store_p_##_S(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_x0,_x1,_x2,_x3,_x4);

#endif

/*$Id$*/



