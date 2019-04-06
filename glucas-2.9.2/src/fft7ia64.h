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


/* This code is written with Intel ia64 in mind. The cycles in the comments
are, of course, only aproximated. These shedule pattern is what we would like
the compiler to do. But, who knows, the compiler can find better solutions */


/* This is a c-macro style of Mlucas (Ernst Mayer) version for Nusbaumer
   FFT length 7 */

#  define cplx_fft7_ia64(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)    \
 cplx_fft7_ia64_##_S(_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)


#  define cplx_fft7_ia64_F(_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)    \
	{                                                               \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai,\
               _b0r,_b0i,_b1r,_b1i,_b2r,_b2i,_b3r,_b3i;                 \
            /* Cycle  1 */                                              \
            cplx_add( _a1, _x1, _x6);                                   \
            /* Cycle  2 */                                              \
            cplx_add( _a2, _x2, _x5);                                   \
            /* Cycle  3 */                                              \
            cplx_add( _a3, _x3, _x4);                                   \
            /* Cycle  4 */                                              \
            cplx_sub( _a6, _x1, _x6);                                   \
            /* Cycle  5 */                                              \
            cplx_sub( _a5, _x2, _x5);                                   \
            /* Cycle  6 */                                              \
            cplx_sub( _a4, _x3, _x4);                                   \
            /* Cycle  7 */                                              \
            cplx_sub( _b0, _a2, _a1);                                   \
            /* Cycle  8 */                                              \
            cplx_add( _b3, _a2, _a3);                                   \
            /* Cycle  9 */                                              \
            cplx_sub( _b1, _a1, _a3);                                   \
            /* Cycle  10 */                                             \
            cplx_sub( _b2, _a3, _a2);                                   \
            /* Cycle  11 */                                             \
	    _b0r *= FN4_7r;                                             \
	    _b0i *= FN4_7r;                                             \
            /* Cycle  12 */                                             \
	    _b1r *= FN2_7r;                                             \
	    _b1i *= FN2_7r;                                             \
            /* Cycle  13 */                                             \
            cplx_add(  _a, _a1, _b3);                                   \
            /* Cycle  14 */                                             \
	    _b2r *= FN3_7r;                                             \
	    _b2i *= FN3_7r;                                             \
            /* Cycle  15 */                                             \
            cplx_add( _a3, _a5, _a6);                                   \
            /* Cycle  16 */                                             \
            /* Cycle  17 */                                             \
            /* Cycle  18 */                                             \
            cplx_add( _x0, _x0, _a);                                    \
            /* Cycle  19 */                                             \
	    _ar = _ar * FN1_7r + _a0r;                                  \
	    _ai = _ai * FN1_7r + _a0i;                                  \
            /* Cycle  20 */                                             \
            cplx_sub( _a0, _a6, _a5);                                   \
            /* Cycle  21 */                                             \
            cplx_add( _a6, _a6, _a4);                                   \
            /* Cycle  22 */                                             \
            cplx_add( _b3, _b1, _b2);                                   \
            /* Cycle  23 */                                             \
            cplx_add( _a5, _a5, _a4);                                   \
            /* Cycle  24 */                                             \
            cplx_sub( _a3, _a4, _a3);                                   \
            /* Cycle  25 */                                             \
            cplx_add( _b1, _b1, _b0);                                   \
            /* Cycle  26 */                                             \
            cplx_sub( _b2, _b0, _b2);                                   \
            /* Cycle  27 */                                             \
            _a0r *= FN4_7i;                                             \
            _a0i *= FN4_7i;                                             \
            /* Cycle  28 */                                             \
	    _a6r *= FN2_7i;                                             \
	    _a6i *= FN2_7i;                                             \
            /* Cycle  29 */                                             \
	    _a5r *= FN3_7i;                                             \
	    _a5i *= FN3_7i;                                             \
            /* Cycle  30 */                                             \
            cplx_add( _b3, _a, _b3);                                    \
            /* Cycle  31 */                                             \
            cplx_sub( _b1, _a, _b1);                                    \
            /* Cycle  32 */                                             \
            cplx_add( _b2, _a, _b2);                                    \
            /* Cycle  33 */                                             \
	    _ar = FN1_7i * _a3r;                                        \
	    _ai = FN1_7i * _a3i;                                        \
            /* Cycle  34 */                                             \
            cplx_sub( _a4, _a0, _a5);                                   \
            /* Cycle  35 */                                             \
            cplx_sub( _a5, _a5, _a6);                                   \
            /* Cycle  36 */                                             \
            cplx_sub( _a6, _a6, _a0);                                   \
            /* Cycle  37 */                                             \
            /* Cycle  38 */                                             \
            /* Cycle  39 */                                             \
            cplx_add( _a4, _a, _a4);                                    \
            /* Cycle  40 */                                             \
            cplx_add( _a5, _a, _a5);                                    \
            /* Cycle  41 */                                             \
            cplx_add( _a6, _a, _a6);                                    \
            /* Cycle  42 */                                             \
            /* Cycle  43 */                                             \
            /* Cycle  44 */                                             \
            _x3##r = _b2r - _a4i;                                       \
            _x4##r = _b2r + _a4i;                                       \
            /* Cycle  45 */                                             \
	    _x3##i = _b2i + _a4r;                                       \
            _x4##i = _b2i - _a4r;                                       \
            /* Cycle  46 */                                             \
	    _x1##r = _b3r + _a5i;                                       \
            _x6##r = _b3r - _a5i;                                       \
            /* Cycle  47 */                                             \
            _x1##i = _b3i - _a5r;                                       \
	    _x6##i = _b3i + _a5r;                                       \
            /* Cycle  48 */                                             \
            _x2##r = _b1r + _a6i;                                       \
	    _x5##r = _b1r - _a6i;                                       \
            /* Cycle  49 */                                             \
            _x2##i = _b1i - _a6r;                                       \
            _x5##i = _b1i + _a6r;                                       \
}

#  define cplx_fft7_ia64_B(_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)    \
	{                                                               \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai,\
               _b0r,_b0i,_b1r,_b1i,_b2r,_b2i,_b3r,_b3i;                 \
            /* Cycle  1 */                                              \
            cplx_add( _a1, _x1, _x6);                                   \
            /* Cycle  2 */                                              \
            cplx_add( _a2, _x2, _x5);                                   \
            /* Cycle  3 */                                              \
            cplx_add( _a3, _x3, _x4);                                   \
            /* Cycle  4 */                                              \
            cplx_sub( _a6, _x1, _x6);                                   \
            /* Cycle  5 */                                              \
            cplx_sub( _a5, _x2, _x5);                                   \
            /* Cycle  6 */                                              \
            cplx_sub( _a4, _x3, _x4);                                   \
            /* Cycle  7 */                                              \
            cplx_sub( _b0, _a2, _a1);                                   \
            /* Cycle  8 */                                              \
            cplx_add( _b3, _a2, _a3);                                   \
            /* Cycle  9 */                                              \
            cplx_sub( _b1, _a1, _a3);                                   \
            /* Cycle  10 */                                             \
            cplx_sub( _b2, _a3, _a2);                                   \
            /* Cycle  11 */                                             \
	    _b0r *= FN4_7r;                                             \
	    _b0i *= FN4_7r;                                             \
            /* Cycle  12 */                                             \
	    _b1r *= FN2_7r;                                             \
	    _b1i *= FN2_7r;                                             \
            /* Cycle  13 */                                             \
            cplx_add(  _a, _a1, _b3);                                   \
            /* Cycle  14 */                                             \
	    _b2r *= FN3_7r;                                             \
	    _b2i *= FN3_7r;                                             \
            /* Cycle  15 */                                             \
            cplx_add( _a3, _a5, _a6);                                   \
            /* Cycle  16 */                                             \
            /* Cycle  17 */                                             \
            /* Cycle  18 */                                             \
            cplx_add( _x0, _x0, _a);                                    \
            /* Cycle  19 */                                             \
	    _ar = _ar * FN1_7r + _a0r;                                  \
	    _ai = _ai * FN1_7r + _a0i;                                  \
            /* Cycle  20 */                                             \
            cplx_sub( _a0, _a6, _a5);                                   \
            /* Cycle  21 */                                             \
            cplx_add( _a6, _a6, _a4);                                   \
            /* Cycle  22 */                                             \
            cplx_add( _b3, _b1, _b2);                                   \
            /* Cycle  23 */                                             \
            cplx_add( _a5, _a5, _a4);                                   \
            /* Cycle  24 */                                             \
            cplx_sub( _a3, _a4, _a3);                                   \
            /* Cycle  25 */                                             \
            cplx_add( _b1, _b1, _b0);                                   \
            /* Cycle  26 */                                             \
            cplx_sub( _b2, _b0, _b2);                                   \
            /* Cycle  27 */                                             \
            _a0r *= FN4_7i;                                             \
            _a0i *= FN4_7i;                                             \
            /* Cycle  28 */                                             \
	    _a6r *= FN2_7i;                                             \
	    _a6i *= FN2_7i;                                             \
            /* Cycle  29 */                                             \
	    _a5r *= FN3_7i;                                             \
	    _a5i *= FN3_7i;                                             \
            /* Cycle  30 */                                             \
            cplx_add( _b3, _a, _b3);                                    \
            /* Cycle  31 */                                             \
            cplx_sub( _b1, _a, _b1);                                    \
            /* Cycle  32 */                                             \
            cplx_add( _b2, _a, _b2);                                    \
            /* Cycle  33 */                                             \
	    _ar = FN1_7i * _a3r;                                        \
	    _ai = FN1_7i * _a3i;                                        \
            /* Cycle  34 */                                             \
            cplx_sub( _a4, _a0, _a5);                                   \
            /* Cycle  35 */                                             \
            cplx_sub( _a5, _a5, _a6);                                   \
            /* Cycle  36 */                                             \
            cplx_sub( _a6, _a6, _a0);                                   \
            /* Cycle  37 */                                             \
            /* Cycle  38 */                                             \
            /* Cycle  39 */                                             \
            cplx_add( _a4, _a, _a4);                                    \
            /* Cycle  40 */                                             \
            cplx_add( _a5, _a, _a5);                                    \
            /* Cycle  41 */                                             \
            cplx_add( _a6, _a, _a6);                                    \
            /* Cycle  42 */                                             \
            /* Cycle  43 */                                             \
            /* Cycle  44 */                                             \
            _x3##r = _b2r + _a4i;                                       \
            _x4##r = _b2r - _a4i;                                       \
            /* Cycle  45 */                                             \
	    _x3##i = _b2i - _a4r;                                       \
            _x4##i = _b2i + _a4r;                                       \
            /* Cycle  46 */                                             \
	    _x1##r = _b3r - _a5i;                                       \
            _x6##r = _b3r + _a5i;                                       \
            /* Cycle  47 */                                             \
            _x1##i = _b3i + _a5r;                                       \
	    _x6##i = _b3i - _a5r;                                       \
            /* Cycle  48 */                                             \
            _x2##r = _b1r - _a6i;                                       \
	    _x5##r = _b1r + _a6i;                                       \
            /* Cycle  49 */                                             \
            _x2##i = _b1i + _a6r;                                       \
            _x5##i = _b1i - _a6r;                                       \
}



# define cplx_fft7_store_ia64(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S,_x0,_x1,_x2,_x3,_x4,_x5,_x6) \
  cplx_fft7_store_ia64_##_S(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3,_x4,_x5,_x6)


# define cplx_fft7_store_ia64_F(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3,_x4,_x5,_x6) \
	{                                                               \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai,\
               _b0r,_b0i,_b1r,_b1i,_b2r,_b2i,_b3r,_b3i;                 \
            /* Cycle  1 */                                              \
            cplx_add( _a1, _x1, _x6);                                   \
            /* Cycle  2 */                                              \
            cplx_add( _a2, _x2, _x5);                                   \
            /* Cycle  3 */                                              \
            cplx_add( _a3, _x3, _x4);                                   \
            /* Cycle  4 */                                              \
            cplx_sub( _a6, _x1, _x6);                                   \
            /* Cycle  5 */                                              \
            cplx_sub( _a5, _x2, _x5);                                   \
            /* Cycle  6 */                                              \
            cplx_sub( _a4, _x3, _x4);                                   \
            /* Cycle  7 */                                              \
            cplx_sub( _b0, _a2, _a1);                                   \
            /* Cycle  8 */                                              \
            cplx_add( _b3, _a2, _a3);                                   \
            /* Cycle  9 */                                              \
            cplx_sub( _b1, _a1, _a3);                                   \
            /* Cycle  10 */                                             \
            cplx_sub( _b2, _a3, _a2);                                   \
            /* Cycle  11 */                                             \
	    _b0r *= FN4_7r;                                             \
	    _b0i *= FN4_7r;                                             \
            /* Cycle  12 */                                             \
	    _b1r *= FN2_7r;                                             \
	    _b1i *= FN2_7r;                                             \
            /* Cycle  13 */                                             \
            cplx_add(  _a, _a1, _b3);                                   \
            /* Cycle  14 */                                             \
	    _b2r *= FN3_7r;                                             \
	    _b2i *= FN3_7r;                                             \
            /* Cycle  15 */                                             \
            cplx_add( _a3, _a5, _a6);                                   \
            /* Cycle  16 */                                             \
            /* Cycle  17 */                                             \
            /* Cycle  18 */                                             \
            cplx_add( _x0, _x0, _a);                                    \
            /* Cycle  19 */                                             \
	    _ar = _ar * FN1_7r + _a0r;                                  \
	    _ai = _ai * FN1_7r + _a0i;                                  \
            /* Cycle  20 */                                             \
            cplx_sub( _a0, _a6, _a5);                                   \
            /* Cycle  21 */                                             \
            cplx_add( _a6, _a6, _a4);                                   \
            /* Cycle  22 */                                             \
            cplx_add( _b3, _b1, _b2);                                   \
            /* Cycle  23 */                                             \
            cplx_add( _a5, _a5, _a4);                                   \
            /* Cycle  24 */                                             \
            cplx_sub( _a3, _a4, _a3);                                   \
            /* Cycle  25 */                                             \
            cplx_add( _b1, _b1, _b0);                                   \
            /* Cycle  26 */                                             \
            cplx_sub( _b2, _b0, _b2);                                   \
            /* Cycle  27 */                                             \
            _a0r *= FN4_7i;                                             \
            _a0i *= FN4_7i;                                             \
            /* Cycle  28 */                                             \
	    _a6r *= FN2_7i;                                             \
	    _a6i *= FN2_7i;                                             \
            /* Cycle  29 */                                             \
	    _a5r *= FN3_7i;                                             \
	    _a5i *= FN3_7i;                                             \
            /* Cycle  30 */                                             \
            cplx_add( _b3, _a, _b3);                                    \
            /* Cycle  31 */                                             \
            cplx_sub( _b1, _a, _b1);                                    \
            /* Cycle  32 */                                             \
            cplx_add( _b2, _a, _b2);                                    \
            /* Cycle  33 */                                             \
	    _ar = FN1_7i * _a3r;                                        \
	    _ai = FN1_7i * _a3i;                                        \
            /* Cycle  34 */                                             \
            cplx_sub( _a4, _a0, _a5);                                   \
            /* Cycle  35 */                                             \
            cplx_sub( _a5, _a5, _a6);                                   \
            /* Cycle  36 */                                             \
            cplx_sub( _a6, _a6, _a0);                                   \
            /* Cycle  37 */                                             \
            *(_pd0##r) = _x0##r;                                        \
            *(_pd0##i) = _x0##i;                                        \
            /* Cycle  38 */                                             \
            /* Cycle  39 */                                             \
            cplx_add( _a4, _a, _a4);                                    \
            /* Cycle  40 */                                             \
            cplx_add( _a5, _a, _a5);                                    \
            /* Cycle  41 */                                             \
            cplx_add( _a6, _a, _a6);                                    \
            /* Cycle  42 */                                             \
            /* Cycle  43 */                                             \
            /* Cycle  44 */                                             \
            _x3##r = _b2r - _a4i;                                       \
            _x4##r = _b2r + _a4i;                                       \
            /* Cycle  45 */                                             \
	    _x3##i = _b2i + _a4r;                                       \
            _x4##i = _b2i - _a4r;                                       \
            /* Cycle  46 */                                             \
	    _x1##r = _b3r + _a5i;                                       \
            _x6##r = _b3r - _a5i;                                       \
            /* Cycle  47 */                                             \
            _x1##i = _b3i - _a5r;                                       \
	    _x6##i = _b3i + _a5r;                                       \
            /* Cycle  48 */                                             \
            _x2##r = _b1r + _a6i;                                       \
	    _x5##r = _b1r - _a6i;                                       \
            /* Cycle  49 */                                             \
            *(_pd3##r) = _x3##r;                                        \
            _x2##i = _b1i - _a6r;                                       \
            *(_pd4##r) = _x4##r;                                        \
            _x5##i = _b1i + _a6r;                                       \
            /* Cycle  50 */                                             \
            *(_pd3##i) = _x3##i;                                        \
            *(_pd4##i) = _x4##i;                                        \
            /* Cycle  51 */                                             \
            *(_pd1##r) = _x1##r;                                        \
            *(_pd6##r) = _x6##r;                                        \
            /* Cycle  52 */                                             \
            *(_pd1##i) = _x1##i;                                        \
            *(_pd6##i) = _x6##i;                                        \
            /* Cycle  53 */                                             \
            *(_pd2##r) = _x2##r;                                        \
            *(_pd5##r) = _x5##r;                                        \
            /* Cycle  54 */                                             \
            *(_pd2##i) = _x2##i;                                        \
            *(_pd5##i) = _x5##i;                                        \
}

# define cplx_fft7_store_ia64_B(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0,_x1,_x2,_x3,_x4,_x5,_x6) \
	{                                                               \
            BIG_DOUBLE _a0r=_x0##r, _a0i=_x0##i, _a1r, _a1i, _a2r, _a2i,\
               _a3r, _a3i, _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai,\
               _b0r,_b0i,_b1r,_b1i,_b2r,_b2i,_b3r,_b3i;                 \
            /* Cycle  1 */                                              \
            cplx_add( _a1, _x1, _x6);                                   \
            /* Cycle  2 */                                              \
            cplx_add( _a2, _x2, _x5);                                   \
            /* Cycle  3 */                                              \
            cplx_add( _a3, _x3, _x4);                                   \
            /* Cycle  4 */                                              \
            cplx_sub( _a6, _x1, _x6);                                   \
            /* Cycle  5 */                                              \
            cplx_sub( _a5, _x2, _x5);                                   \
            /* Cycle  6 */                                              \
            cplx_sub( _a4, _x3, _x4);                                   \
            /* Cycle  7 */                                              \
            cplx_sub( _b0, _a2, _a1);                                   \
            /* Cycle  8 */                                              \
            cplx_add( _b3, _a2, _a3);                                   \
            /* Cycle  9 */                                              \
            cplx_sub( _b1, _a1, _a3);                                   \
            /* Cycle  10 */                                             \
            cplx_sub( _b2, _a3, _a2);                                   \
            /* Cycle  11 */                                             \
	    _b0r *= FN4_7r;                                             \
	    _b0i *= FN4_7r;                                             \
            /* Cycle  12 */                                             \
	    _b1r *= FN2_7r;                                             \
	    _b1i *= FN2_7r;                                             \
            /* Cycle  13 */                                             \
            cplx_add(  _a, _a1, _b3);                                   \
            /* Cycle  14 */                                             \
	    _b2r *= FN3_7r;                                             \
	    _b2i *= FN3_7r;                                             \
            /* Cycle  15 */                                             \
            cplx_add( _a3, _a5, _a6);                                   \
            /* Cycle  16 */                                             \
            /* Cycle  17 */                                             \
            /* Cycle  18 */                                             \
            cplx_add( _x0, _x0, _a);                                    \
            /* Cycle  19 */                                             \
	    _ar = _ar * FN1_7r + _a0r;                                  \
	    _ai = _ai * FN1_7r + _a0i;                                  \
            /* Cycle  20 */                                             \
            cplx_sub( _a0, _a6, _a5);                                   \
            /* Cycle  21 */                                             \
            cplx_add( _a6, _a6, _a4);                                   \
            /* Cycle  22 */                                             \
            cplx_add( _b3, _b1, _b2);                                   \
            /* Cycle  23 */                                             \
            cplx_add( _a5, _a5, _a4);                                   \
            /* Cycle  24 */                                             \
            cplx_sub( _a3, _a4, _a3);                                   \
            /* Cycle  25 */                                             \
            cplx_add( _b1, _b1, _b0);                                   \
            /* Cycle  26 */                                             \
            cplx_sub( _b2, _b0, _b2);                                   \
            /* Cycle  27 */                                             \
            _a0r *= FN4_7i;                                             \
            _a0i *= FN4_7i;                                             \
            /* Cycle  28 */                                             \
	    _a6r *= FN2_7i;                                             \
	    _a6i *= FN2_7i;                                             \
            /* Cycle  29 */                                             \
	    _a5r *= FN3_7i;                                             \
	    _a5i *= FN3_7i;                                             \
            /* Cycle  30 */                                             \
            cplx_add( _b3, _a, _b3);                                    \
            /* Cycle  31 */                                             \
            cplx_sub( _b1, _a, _b1);                                    \
            /* Cycle  32 */                                             \
            cplx_add( _b2, _a, _b2);                                    \
            /* Cycle  33 */                                             \
	    _ar = FN1_7i * _a3r;                                        \
	    _ai = FN1_7i * _a3i;                                        \
            /* Cycle  34 */                                             \
            cplx_sub( _a4, _a0, _a5);                                   \
            /* Cycle  35 */                                             \
            cplx_sub( _a5, _a5, _a6);                                   \
            /* Cycle  36 */                                             \
            cplx_sub( _a6, _a6, _a0);                                   \
            /* Cycle  37 */                                             \
            *(_pd0##r) = _x0##r;                                        \
            *(_pd0##i) = _x0##i;                                        \
            /* Cycle  38 */                                             \
            /* Cycle  39 */                                             \
            cplx_add( _a4, _a, _a4);                                    \
            /* Cycle  40 */                                             \
            cplx_add( _a5, _a, _a5);                                    \
            /* Cycle  41 */                                             \
            cplx_add( _a6, _a, _a6);                                    \
            /* Cycle  42 */                                             \
            /* Cycle  43 */                                             \
            /* Cycle  44 */                                             \
            _x3##r = _b2r + _a4i;                                       \
            _x4##r = _b2r - _a4i;                                       \
            /* Cycle  45 */                                             \
	    _x3##i = _b2i - _a4r;                                       \
            _x4##i = _b2i + _a4r;                                       \
            /* Cycle  46 */                                             \
	    _x1##r = _b3r - _a5i;                                       \
            _x6##r = _b3r + _a5i;                                       \
            /* Cycle  47 */                                             \
            _x1##i = _b3i + _a5r;                                       \
	    _x6##i = _b3i - _a5r;                                       \
            /* Cycle  48 */                                             \
            _x2##r = _b1r - _a6i;                                       \
	    _x5##r = _b1r + _a6i;                                       \
            /* Cycle  49 */                                             \
            *(_pd3##r) = _x3##r;                                        \
            _x2##i = _b1i + _a6r;                                       \
            *(_pd4##r) = _x4##r;                                        \
            _x5##i = _b1i - _a6r;                                       \
            /* Cycle  50 */                                             \
            *(_pd3##i) = _x3##i;                                        \
            *(_pd4##i) = _x4##i;                                        \
            /* Cycle  51 */                                             \
            *(_pd1##r) = _x1##r;                                        \
            *(_pd6##r) = _x6##r;                                        \
            /* Cycle  52 */                                             \
            *(_pd1##i) = _x1##i;                                        \
            *(_pd6##i) = _x6##i;                                        \
            /* Cycle  53 */                                             \
            *(_pd2##r) = _x2##r;                                        \
            *(_pd5##r) = _x5##r;                                        \
            /* Cycle  54 */                                             \
            *(_pd2##i) = _x2##i;                                        \
            *(_pd5##i) = _x5##i;                                        \
}



/*$Id$*/






