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
#if (Y_TARGET == 1 || Y_TARGET == 11 || Y_TARGET == 12 || Y_TARGET == 13)

# include "fft7x86.h"

#else

# if defined(__GNUC__) && ( __GNUC__ < 3 ) /*|| defined(__INTEL_COMPILER)*/

#  define cplx_fft7_F( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);         \
            cplx_add( _a2, _x2, _x5);         \
            cplx_add( _a3, _x3, _x4);         \
            cplx_sub( _a6, _x1, _x6);         \
            cplx_add( _a, _a1, _a2);          \
            cplx_sub( _a5, _x2, _x5);         \
            cplx_sub( _x2, _a2, _a1);         \
            cplx_sub( _a4, _x3, _x4);         \
            cplx_add( _a, _a, _a3);           \
            cplx_sub( _x1, _a1, _a3);         \
            cplx_sub( _x3, _a3, _a2);         \
            cplx_add( _x0, _x0, _a);          \
	    _ar = _a0r + _ar * FN1_7r;        \
	    _ai = _a0i + _ai * FN1_7r;        \
            cplx_add(_x4,_a5,_a6);            \
	    _a0r = FN4_7r * _x2##r;           \
	    _a0i = FN4_7r * _x2##i;           \
            cplx_sub(_x6,_a6,_a5);            \
	    _a1r = FN2_7r * _x1##r;           \
	    _a1i = FN2_7r * _x1##i;           \
            cplx_sub(_x4,_a4,_x4);            \
	    _a2r = FN3_7r * _x3##r;           \
	    _a2i = FN3_7r * _x3##i;           \
            cplx_add(_x5,_a5,_a4);            \
            cplx_add(_a6,_a6,_a4);            \
            cplx_add(_x1,_a1,_a0);            \
            _x4##r *= FN1_7i;                 \
            _x4##i *= FN1_7i;                 \
            cplx_add(_x3,_a1,_a2);            \
            cplx_sub(_x2,_a0,_a2);            \
            _a0r = FN4_7i * _x6##r;           \
            _a0i = FN4_7i * _x6##i;           \
            cplx_sub(_a1,_a,_x1);             \
            _a5r = FN3_7i * _x5##r;           \
            _a5i = FN3_7i * _x5##i;           \
            cplx_add(_a3,_a,_x3);             \
            _a6r *= FN2_7i ;                  \
            _a6i *= FN2_7i ;                  \
            cplx_add(_a2,_a,_x2);             \
            cplx_sub(_x1,_a0,_a5);            \
            cplx_sub(_a5,_a5,_a6);            \
            cplx_sub(_a6,_a6,_a0);            \
            cplx_add(_a4,_x4,_x1);            \
            cplx_add(_a5,_a5,_x4);            \
            cplx_add(_a6,_a6,_x4);            \
	                                      \
            _x1##i = _a3i - _a5r;             \
	    _x6##i = _a3i + _a5r;             \
	    _x1##r = _a3r + _a5i;             \
            _x6##r = _a3r - _a5i;             \
            _x2##r = _a1r + _a6i;             \
	    _x5##r = _a1r - _a6i;             \
            _x2##i = _a1i - _a6r;             \
            _x5##i = _a1i + _a6r;             \
            _x3##r = _a2r - _a4i;             \
            _x4##r = _a2r + _a4i;             \
	    _x3##i = _a2i + _a4r;             \
            _x4##i = _a2i - _a4r;             \
	}

#  define cplx_fft7_B( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);         \
            cplx_add( _a2, _x2, _x5);         \
            cplx_add( _a3, _x3, _x4);         \
            cplx_sub( _a6, _x1, _x6);         \
            cplx_add( _a, _a1, _a2);          \
            cplx_sub( _a5, _x2, _x5);         \
            cplx_sub( _x2, _a2, _a1);         \
            cplx_sub( _a4, _x3, _x4);         \
            cplx_add( _a, _a, _a3);           \
            cplx_sub( _x1, _a1, _a3);         \
            cplx_sub( _x3, _a3, _a2);         \
            cplx_add( _x0, _x0, _a);          \
	    _ar = _a0r + _ar * FN1_7r;        \
	    _ai = _a0i + _ai * FN1_7r;        \
            cplx_add(_x4,_a5,_a6);            \
	    _a0r = FN4_7r * _x2##r;           \
	    _a0i = FN4_7r * _x2##i;           \
            cplx_sub(_x6,_a6,_a5);            \
	    _a1r = FN2_7r * _x1##r;           \
	    _a1i = FN2_7r * _x1##i;           \
            cplx_sub(_x4,_a4,_x4);            \
	    _a2r = FN3_7r * _x3##r;           \
	    _a2i = FN3_7r * _x3##i;           \
            cplx_add(_x5,_a5,_a4);            \
            cplx_add(_a6,_a6,_a4);            \
            cplx_add(_x1,_a1,_a0);            \
            _x4##r *= FN1_7i;                 \
            _x4##i *= FN1_7i;                 \
            cplx_add(_x3,_a1,_a2);            \
            cplx_sub(_x2,_a0,_a2);            \
            _a0r = FN4_7i * _x6##r;           \
            _a0i = FN4_7i * _x6##i;           \
            cplx_sub(_a1,_a,_x1);             \
            _a5r = FN3_7i * _x5##r;           \
            _a5i = FN3_7i * _x5##i;           \
            cplx_add(_a3,_a,_x3);             \
            _a6r *= FN2_7i ;                  \
            _a6i *= FN2_7i ;                  \
            cplx_add(_a2,_a,_x2);             \
            cplx_sub(_x1,_a0,_a5);            \
            cplx_sub(_a5,_a5,_a6);            \
            cplx_sub(_a6,_a6,_a0);            \
            cplx_add(_a4,_x4,_x1);            \
            cplx_add(_a5,_a5,_x4);            \
            cplx_add(_a6,_a6,_x4);            \
                                              \
            _x1##i = _a3i + _a5r;             \
	    _x6##i = _a3i - _a5r;             \
	    _x1##r = _a3r - _a5i;             \
            _x6##r = _a3r + _a5i;             \
            _x2##r = _a1r - _a6i;             \
	    _x5##r = _a1r + _a6i;             \
            _x2##i = _a1i + _a6r;             \
            _x5##i = _a1i - _a6r;             \
            _x3##r = _a2r + _a4i;             \
            _x4##r = _a2r - _a4i;             \
	    _x3##i = _a2i - _a4r;             \
            _x4##i = _a2i + _a4r;             \
	}

/*
#  define cplx_fft7(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{ \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _ar, _ai;\
            cplx_addsub(_x1,_x6);\
            cplx_addsub(_x2,_x5);\
            cplx_addsub(_x3,_x4);\
            _ar = (_x1##r + _x2##r + _x3##r);\
            _ai = (_x1##i + _x2##i + _x3##i);\
            _x0##r += _ar;\
            _x0##i += _ai;\
	    _ar = _a0r + _ar * _S##N1_7r;\
	    _ai = _a0i + _ai * _S##N1_7r;\
	    _a0r = _S##N4_7r * ( _x2##r - _x1##r );\
	    _a0i = _S##N4_7r * ( _x2##i - _x1##i );\
	    _x1##r = _S##N2_7r * ( _x1##r - _x3##r );\
	    _x1##i = _S##N2_7r * ( _x1##i - _x3##i );\
	    _x2##r = _S##N3_7r * ( _x3##r - _x2##r );\
	    _x2##i = _S##N3_7r * ( _x3##i - _x2##i );\
	    _x3##r = _ar + _x1##r + _x2##r;\
	    _x3##i = _ai + _x1##i + _x2##i;\
	    _x1##r = _ar - _x1##r - _a0r;\
	    _x1##i = _ai - _x1##i - _a0i;\
	    _x2##r = _ar - _x2##r + _a0r;\
	    _x2##i = _ai - _x2##i + _a0i;\
	    \
	    _ar = _S##N1_7i * ( _x4##r - _x5##r - _x6##r);\
	    _ai = _S##N1_7i * ( _x4##i - _x5##i - _x6##i);\
	    _a0r = _S##N4_7i * ( _x6##r - _x5##r );\
	    _a0i = _S##N4_7i * ( _x6##i - _x5##i );\
	    _x6##r = _S##N2_7i * ( _x6##r + _x4##r );\
	    _x6##i = _S##N2_7i * ( _x6##i + _x4##i );\
	    _x5##r = _S##N3_7i * ( _x5##r + _x4##r );\
	    _x5##i = _S##N3_7i * ( _x5##i + _x4##i );\
	    _x4##r = _ar - _x5##r + _a0r;\
	    _x4##i = _ai - _x5##i + _a0i;\
	    _x5##r += _ar - _x6##r;\
	    _x5##i += _ai - _x6##i;\
	    _ar += _x6##r - _a0r;\
	    _ai += _x6##i - _a0i;\
            _a0r = _x1##r;\
            _a0i = _x1##i;\
	    \
	    _x1##r = _x3##r + _x5##i;\
            _x6##r = _x3##r - _x5##i;\
            _x1##i = _x3##i - _x5##r;\
	    _x6##i = _x3##i + _x5##r;\
            _a1r = _x2##r;\
            _a1i = _x2##i;\
            _x2##r = _a0r + _ai;\
	    _x5##r = _a0r - _ai;\
            _x2##i = _a0i - _ar;\
            _x5##i = _a0i + _ar;\
            _a0r=_x4##r;\
            _a0i=_x4##i;\
            _x3##r = _a1r - _a0i;\
            _x4##r = _a1r + _a0i;\
	    _x3##i = _a1i + _a0r;\
            _x4##i = _a1i - _a0r;\
	}

*/

# define cplx_fft7_store_p_F(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3,_x4,_x5,_x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);         \
            cplx_add( _a2, _x2, _x5);         \
            cplx_add( _a3, _x3, _x4);         \
            cplx_sub( _a6, _x1, _x6);         \
            cplx_add( _a, _a1, _a2);          \
            cplx_sub( _a5, _x2, _x5);         \
            cplx_sub( _x2, _a2, _a1);         \
            cplx_sub( _a4, _x3, _x4);         \
            cplx_add( _a, _a, _a3);           \
            cplx_sub( _x1, _a1, _a3);         \
            cplx_sub( _x3, _a3, _a2);         \
            cplx_add( _x0, _x0, _a);          \
	    _ar = _a0r + _ar * FN1_7r;        \
	    _ai = _a0i + _ai * FN1_7r;        \
            cplx_add(_x4,_a5,_a6);            \
	    _a0r = FN4_7r * _x2##r;           \
	    _a0i = FN4_7r * _x2##i;           \
            cplx_sub(_x6,_a6,_a5);            \
	    _a1r = FN2_7r * _x1##r;           \
	    _a1i = FN2_7r * _x1##i;           \
            *(_pd0 + _index) = _x0##r;        \
            *(_pd0 + _index + 1) = _x0##i;    \
            cplx_sub(_x4,_a4,_x4);            \
	    _a2r = FN3_7r * _x3##r;           \
	    _a2i = FN3_7r * _x3##i;           \
            cplx_add(_x5,_a5,_a4);            \
            cplx_add(_a6,_a6,_a4);            \
            cplx_add(_x1,_a1,_a0);            \
            _x4##r *= FN1_7i;                 \
            _x4##i *= FN1_7i;                 \
            cplx_add(_x3,_a1,_a2);            \
            cplx_sub(_x2,_a0,_a2);            \
            _a0r = FN4_7i * _x6##r;           \
            _a0i = FN4_7i * _x6##i;           \
            cplx_sub(_a1,_a,_x1);             \
            _a5r = FN3_7i * _x5##r;           \
            _a5i = FN3_7i * _x5##i;           \
            cplx_add(_a3,_a,_x3);             \
            _a6r *= FN2_7i ;                  \
            _a6i *= FN2_7i ;                  \
            cplx_add(_a2,_a,_x2);             \
            cplx_sub(_x1,_a0,_a5);            \
            cplx_sub(_a5,_a5,_a6);            \
            cplx_sub(_a6,_a6,_a0);            \
            cplx_add(_a4,_x4,_x1);            \
            cplx_add(_a5,_a5,_x4);            \
            cplx_add(_a6,_a6,_x4);            \
	                                      \
	    *(_pd1 + _index)= _a3r + _a5i;    \
            *(_pd6 + _index)= _a3r - _a5i;    \
            *(_pd1 + _index + 1)= _a3i - _a5r;\
	    *(_pd6 + _index + 1)= _a3i + _a5r;\
            *(_pd2 + _index)= _a1r + _a6i;    \
	    *(_pd5 + _index)= _a1r - _a6i;    \
            *(_pd2 + _index + 1)= _a1i - _a6r;\
            *(_pd5 + _index + 1)= _a1i + _a6r;\
            *(_pd3 + _index)= _a2r - _a4i;    \
            *(_pd4 + _index)= _a2r + _a4i;    \
	    *(_pd3 + _index + 1)= _a2i + _a4r;\
            *(_pd4 + _index + 1)= _a2i - _a4r;\
	}

# define cplx_fft7_store_p_B(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3,_x4,_x5,_x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);         \
            cplx_add( _a2, _x2, _x5);         \
            cplx_add( _a3, _x3, _x4);         \
            cplx_sub( _a6, _x1, _x6);         \
            cplx_add( _a, _a1, _a2);          \
            cplx_sub( _a5, _x2, _x5);         \
            cplx_sub( _x2, _a2, _a1);         \
            cplx_sub( _a4, _x3, _x4);         \
            cplx_add( _a, _a, _a3);           \
            cplx_sub( _x1, _a1, _a3);         \
            cplx_sub( _x3, _a3, _a2);         \
            cplx_add( _x0, _x0, _a);          \
	    _ar = _a0r + _ar * FN1_7r;        \
	    _ai = _a0i + _ai * FN1_7r;        \
            cplx_add(_x4,_a5,_a6);            \
	    _a0r = FN4_7r * _x2##r;           \
	    _a0i = FN4_7r * _x2##i;           \
            cplx_sub(_x6,_a6,_a5);            \
	    _a1r = FN2_7r * _x1##r;           \
	    _a1i = FN2_7r * _x1##i;           \
            *(_pd0 + _index) = _x0##r;        \
            *(_pd0 + _index + 1) = _x0##i;    \
            cplx_sub(_x4,_a4,_x4);            \
	    _a2r = FN3_7r * _x3##r;           \
	    _a2i = FN3_7r * _x3##i;           \
            cplx_add(_x5,_a5,_a4);            \
            cplx_add(_a6,_a6,_a4);            \
            cplx_add(_x1,_a1,_a0);            \
            _x4##r *= FN1_7i;                 \
            _x4##i *= FN1_7i;                 \
            cplx_add(_x3,_a1,_a2);            \
            cplx_sub(_x2,_a0,_a2);            \
            _a0r = FN4_7i * _x6##r;           \
            _a0i = FN4_7i * _x6##i;           \
            cplx_sub(_a1,_a,_x1);             \
            _a5r = FN3_7i * _x5##r;           \
            _a5i = FN3_7i * _x5##i;           \
            cplx_add(_a3,_a,_x3);             \
            _a6r *= FN2_7i ;                  \
            _a6i *= FN2_7i ;                  \
            cplx_add(_a2,_a,_x2);             \
            cplx_sub(_x1,_a0,_a5);            \
            cplx_sub(_a5,_a5,_a6);            \
            cplx_sub(_a6,_a6,_a0);            \
            cplx_add(_a4,_x4,_x1);            \
            cplx_add(_a5,_a5,_x4);            \
            cplx_add(_a6,_a6,_x4);            \
	                                      \
	    *(_pd1 + _index)= _a3r - _a5i;    \
            *(_pd6 + _index)= _a3r + _a5i;    \
            *(_pd1 + _index + 1)= _a3i + _a5r;\
	    *(_pd6 + _index + 1)= _a3i - _a5r;\
            *(_pd2 + _index)= _a1r - _a6i;    \
	    *(_pd5 + _index)= _a1r + _a6i;    \
            *(_pd2 + _index + 1)= _a1i + _a6r;\
            *(_pd5 + _index + 1)= _a1i - _a6r;\
            *(_pd3 + _index)= _a2r + _a4i;    \
            *(_pd4 + _index)= _a2r - _a4i;    \
	    *(_pd3 + _index + 1)= _a2i - _a4r;\
            *(_pd4 + _index + 1)= _a2i + _a4r;\
	}

# else

/* For other compilers */


#  define cplx_fft7_F( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);\
            cplx_add( _a2, _x2, _x5);\
            cplx_add( _a3, _x3, _x4);\
            cplx_sub( _a6, _x1, _x6);\
            cplx_sub( _a5, _x2, _x5);\
            cplx_sub( _a4, _x3, _x4);\
            _ar = (_a1r + _a2r + _a3r);\
            _ai = (_a1i + _a2i + _a3i);\
            _x0##r += _ar;\
            _x0##i += _ai;\
	    _ar = _a0r + _ar * FN1_7r;\
	    _ai = _a0i + _ai * FN1_7r;\
	    _a0r = FN4_7r * ( _a2r - _a1r );\
	    _a0i = FN4_7r * ( _a2i - _a1i );\
	    _a1r = FN2_7r * ( _a1r - _a3r );\
	    _a1i = FN2_7r * ( _a1i - _a3i );\
	    _a2r = FN3_7r * ( _a3r - _a2r );\
	    _a2i = FN3_7r * ( _a3i - _a2i );\
	    _a3r = _ar + _a1r + _a2r;\
	    _a3i = _ai + _a1i + _a2i;\
	    _a1r = _ar - _a1r - _a0r;\
	    _a1i = _ai - _a1i - _a0i;\
	    _a2r = _ar - _a2r + _a0r;\
	    _a2i = _ai - _a2i + _a0i;\
	    \
	    _ar = FN1_7i * ( _a4r - _a5r - _a6r);\
	    _ai = FN1_7i * ( _a4i - _a5i - _a6i);\
	    _a0r = FN4_7i * ( _a6r - _a5r );\
	    _a0i = FN4_7i * ( _a6i - _a5i );\
	    _a6r = FN2_7i * ( _a6r + _a4r );\
	    _a6i = FN2_7i * ( _a6i + _a4i );\
	    _a5r = FN3_7i * ( _a5r + _a4r );\
	    _a5i = FN3_7i * ( _a5i + _a4i );\
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

#  define cplx_fft7_B( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{                                                     \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;\
            cplx_add( _a1, _x1, _x6);\
            cplx_add( _a2, _x2, _x5);\
            cplx_add( _a3, _x3, _x4);\
            cplx_sub( _a6, _x1, _x6);\
            cplx_sub( _a5, _x2, _x5);\
            cplx_sub( _a4, _x3, _x4);\
            _ar = (_a1r + _a2r + _a3r);\
            _ai = (_a1i + _a2i + _a3i);\
            _x0##r += _ar;\
            _x0##i += _ai;\
	    _ar = _a0r + _ar * FN1_7r;\
	    _ai = _a0i + _ai * FN1_7r;\
	    _a0r = FN4_7r * ( _a2r - _a1r );\
	    _a0i = FN4_7r * ( _a2i - _a1i );\
	    _a1r = FN2_7r * ( _a1r - _a3r );\
	    _a1i = FN2_7r * ( _a1i - _a3i );\
	    _a2r = FN3_7r * ( _a3r - _a2r );\
	    _a2i = FN3_7r * ( _a3i - _a2i );\
	    _a3r = _ar + _a1r + _a2r;\
	    _a3i = _ai + _a1i + _a2i;\
	    _a1r = _ar - _a1r - _a0r;\
	    _a1i = _ai - _a1i - _a0i;\
	    _a2r = _ar - _a2r + _a0r;\
	    _a2i = _ai - _a2i + _a0i;\
	    \
	    _ar = FN1_7i * ( _a4r - _a5r - _a6r);\
	    _ai = FN1_7i * ( _a4i - _a5i - _a6i);\
	    _a0r = FN4_7i * ( _a6r - _a5r );\
	    _a0i = FN4_7i * ( _a6i - _a5i );\
	    _a6r = FN2_7i * ( _a6r + _a4r );\
	    _a6i = FN2_7i * ( _a6i + _a4i );\
	    _a5r = FN3_7i * ( _a5r + _a4r );\
	    _a5i = FN3_7i * ( _a5i + _a4i );\
	    _a4r = _ar - _a5r + _a0r;\
	    _a4i = _ai - _a5i + _a0i;\
	    _a5r += _ar - _a6r;\
	    _a5i += _ai - _a6i;\
	    _a6r += _ar - _a0r;\
	    _a6i += _ai - _a0i;\
	    \
	    _x1##r = _a3r - _a5i;\
            _x6##r = _a3r + _a5i;\
            _x1##i = _a3i + _a5r;\
	    _x6##i = _a3i - _a5r;\
            _x2##r = _a1r - _a6i;\
	    _x5##r = _a1r + _a6i;\
            _x2##i = _a1i + _a6r;\
            _x5##i = _a1i - _a6r;\
            _x3##r = _a2r + _a4i;\
            _x4##r = _a2r - _a4i;\
	    _x3##i = _a2i - _a4r;\
            _x4##i = _a2i + _a4r;\
	}

#  define cplx_fft7_store_p_F(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{ \
	    BIG_DOUBLE  _ar, _ai;\
	    cplx_addsub( _x1, _x6);\
	    cplx_addsub( _x2, _x5);\
	    cplx_addsub( _x3, _x4);\
	    _ar = (_x1##r + _x2##r + _x3##r);\
	    _ai = (_x1##i + _x2##i + _x3##i);\
	    *(_pd0 + _index) = _x0##r +  _ar;\
	    *(_pd0 + _index + 1) = _x0##i + _ai;\
	    _ar = _x0##r + _ar * FN1_7r;\
	    _ai = _x0##i + _ai * FN1_7r;\
	    _x0##r = FN4_7r * ( _x2##r - _x1##r );\
	    _x0##i = FN4_7r * ( _x2##i - _x1##i );\
	    _x1##r = FN2_7r * ( _x1##r - _x3##r );\
	    _x1##i = FN2_7r * ( _x1##i - _x3##i );\
	    _x2##r = FN3_7r * ( _x3##r - _x2##r );\
	    _x2##i = FN3_7r * ( _x3##i - _x2##i );\
	    _x3##r = _ar + _x1##r + _x2##r;\
	    _x3##i = _ai + _x1##i + _x2##i;\
	    _x1##r = _ar - _x1##r - _x0##r;\
	    _x1##i = _ai - _x1##i - _x0##i;\
	    _x2##r = _ar - _x2##r + _x0##r;\
	    _x2##i = _ai - _x2##i + _x0##i;\
	    \
	    _ar = FN1_7i * ( _x4##r - _x5##r - _x6##r);\
	    _ai = FN1_7i * ( _x4##i - _x5##i - _x6##i);\
	    _x0##r = FN4_7i * ( _x6##r - _x5##r );\
	    _x0##i = FN4_7i * ( _x6##i - _x5##i );\
	    _x6##r = FN2_7i * ( _x6##r + _x4##r );\
	    _x6##i = FN2_7i * ( _x6##i + _x4##i );\
	    _x5##r = FN3_7i * ( _x5##r + _x4##r );\
	    _x5##i = FN3_7i * ( _x5##i + _x4##i );\
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
            _x0##r = *( _pd0 + _index);\
            _x0##i = *( _pd0 + _index + 1);\
	}


#  define cplx_fft7_store_p_B(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
	{ \
	    BIG_DOUBLE  _ar, _ai;\
	    cplx_addsub( _x1, _x6);\
	    cplx_addsub( _x2, _x5);\
	    cplx_addsub( _x3, _x4);\
	    _ar = (_x1##r + _x2##r + _x3##r);\
	    _ai = (_x1##i + _x2##i + _x3##i);\
	    *(_pd0 + _index) = _x0##r +  _ar;\
	    *(_pd0 + _index + 1) = _x0##i + _ai;\
	    _ar = _x0##r + _ar * FN1_7r;\
	    _ai = _x0##i + _ai * FN1_7r;\
	    _x0##r = FN4_7r * ( _x2##r - _x1##r );\
	    _x0##i = FN4_7r * ( _x2##i - _x1##i );\
	    _x1##r = FN2_7r * ( _x1##r - _x3##r );\
	    _x1##i = FN2_7r * ( _x1##i - _x3##i );\
	    _x2##r = FN3_7r * ( _x3##r - _x2##r );\
	    _x2##i = FN3_7r * ( _x3##i - _x2##i );\
	    _x3##r = _ar + _x1##r + _x2##r;\
	    _x3##i = _ai + _x1##i + _x2##i;\
	    _x1##r = _ar - _x1##r - _x0##r;\
	    _x1##i = _ai - _x1##i - _x0##i;\
	    _x2##r = _ar - _x2##r + _x0##r;\
	    _x2##i = _ai - _x2##i + _x0##i;\
	    \
	    _ar = FN1_7i * ( _x4##r - _x5##r - _x6##r);\
	    _ai = FN1_7i * ( _x4##i - _x5##i - _x6##i);\
	    _x0##r = FN4_7i * ( _x6##r - _x5##r );\
	    _x0##i = FN4_7i * ( _x6##i - _x5##i );\
	    _x6##r = FN2_7i * ( _x6##r + _x4##r );\
	    _x6##i = FN2_7i * ( _x6##i + _x4##i );\
	    _x5##r = FN3_7i * ( _x5##r + _x4##r );\
	    _x5##i = FN3_7i * ( _x5##i + _x4##i );\
	    _x4##r = _ar - _x5##r + _x0##r;\
	    _x4##i = _ai - _x5##i + _x0##i;\
	    _x5##r += _ar - _x6##r;\
	    _x5##i += _ai - _x6##i;\
	    _x6##r += _ar - _x0##r;\
	    _x6##i += _ai - _x0##i;\
	    \
	    *( _pd1 + _index ) = _x3##r - _x5##i;\
	    *( _pd6 + _index ) = _x3##r + _x5##i;\
	    *( _pd1 + _index + 1) = _x3##i + _x5##r;\
	    *( _pd6 + _index + 1) = _x3##i - _x5##r;\
	    *( _pd2 + _index) = _x1##r - _x6##i;\
	    *( _pd5 + _index) = _x1##r + _x6##i;\
	    *( _pd2 + _index + 1 ) = _x1##i + _x6##r;\
	    *( _pd5 + _index + 1 ) = _x1##i - _x6##r;\
	    *( _pd3 + _index) = _x2##r + _x4##i;\
	    *( _pd4 + _index) = _x2##r - _x4##i;\
	    *( _pd3 + _index + 1 ) = _x2##i - _x4##r;\
	    *( _pd4 + _index + 1 ) = _x2##i + _x4##r;\
	}

# endif

# define cplx_fft7(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
  cplx_fft7_##_S( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

# define cplx_fft7_store_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
 cplx_fft7_store_p_##_S(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

#endif

/*$Id$*/






