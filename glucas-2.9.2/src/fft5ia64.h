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

/* Nussbaumer version for FFT-5. Code adapted from E. Mayer's (Mlucas) */


#define cplx_fft5_ia64(_S, _x0, _x1, _x2, _x3, _x4)         \
{                                                           \
     y_limb_t _ar, _ai, _b0r, _b0i, _b1r, _b1i, _b2r, _b2i; \
     y_limb_t _b3r, _b3i, _b4r, _b4i;                       \
     /* Cycle  1 */                                         \
     cplx_sub( _b4, _x1, _x4);                              \
     /* Cycle  2 */                                         \
     cplx_sub( _b3, _x2, _x3);                              \
     /* Cycle  3 */                                         \
     cplx_add( _b1, _x1, _x4);                              \
     /* Cycle  4 */                                         \
     cplx_add( _b2, _x2, _x3);                              \
     /* Cycle  5 */                                         \
     /* Cycle  6 */                                         \
     /* Cycle  7 */                                         \
     cplx_sub( _b0, _b4, _b3);                              \
     /* Cycle  8 */                                         \
     _x3##r = _b3r * _S##N1_5i;                             \
     _x3##i = _b3i * _S##N1_5i;                             \
     /* Cycle  9 */                                         \
     cplx_add( _a, _b1, _b2);                               \
     /* Cycle  10 */                                        \
     _x4##r = _b4r * _S##N2_5i;                             \
     _x4##i = _b4i * _S##N2_5i;                             \
     /* Cycle  11 */                                        \
     cplx_sub( _x2, _b1, _b2);                              \
     /* Cycle  12 */                                        \
     _b0r *= _S##_1_5i;                                     \
     _b0i *= _S##_1_5i;                                     \
     /* Cycle  13 */                                        \
     _b3r = _ar * (-0.25) + _x0##r;                         \
     _b3i = _ai * (-0.25) + _x0##i;                         \
     /* Cycle  14 */                                        \
     cplx_add( _x0, _a, _x0);                               \
     /* Cycle  15 */                                        \
     /* Cycle  16 */                                        \
     _x2##r *= _S##N2_5r;                                   \
     _x2##i *= _S##N2_5r;                                   \
     /* Cycle  17 */                                        \
     cplx_add( _x3, _x3, _b0);                              \
     /* Cycle  18 */                                        \
     cplx_sub( _b4, _b0, _x4);                              \
     /* Cycle  20 */                                        \
     /* Cycle  21 */                                        \
     cplx_add( _x1, _b3, _x2);                              \
     /* Cycle  22 */                                        \
     cplx_sub( _x2, _b3, _x2);                              \
     /* Cycle  23 */                                        \
     /* Cycle  24 */                                        \
     /* Cycle  25 */                                        \
     _x4##r = _x1##r + _x3##i;                              \
     _x4##i = _x1##i - _x3##r;                              \
     /* Cycle  26 */                                        \
     _x1##r -= _x3##i;                                      \
     _x1##i += _x3##r;                                      \
     /* Cycle  27 */                                        \
     _x3##r = _x2##r + _b4i;                                \
     _x3##i = _x2##i - _b4r;                                \
     /* Cycle  28 */                                        \
     _x2##r -= _b4i;                                        \
     _x2##i += _b4r;                                        \
}

#define cplx_fft5_store_ia64(_pd0,_pd1,_pd2,_pd3,_pd4,_S, _x0, _x1, _x2, _x3, _x4) \
{                                                           \
     y_limb_t _ar, _ai, _b0r, _b0i, _b1r, _b1i, _b2r, _b2i; \
     y_limb_t _b3r, _b3i, _b4r, _b4i;                      \
     /* Cycle  1 */                                         \
     cplx_sub( _b4, _x1, _x4);                              \
     /* Cycle  2 */                                         \
     cplx_sub( _b3, _x2, _x3);                              \
     /* Cycle  3 */                                         \
     cplx_add( _b1, _x1, _x4);                              \
     /* Cycle  4 */                                         \
     cplx_add( _b2, _x2, _x3);                              \
     /* Cycle  5 */                                         \
     /* Cycle  6 */                                         \
     /* Cycle  7 */                                         \
     cplx_sub( _b0, _b4, _b3);                              \
     /* Cycle  8 */                                         \
     _x3##r = _b3r * _S##N1_5i;                             \
     _x3##i = _b3i * _S##N1_5i;                             \
     /* Cycle  9 */                                         \
     cplx_add( _a, _b1, _b2);                               \
     /* Cycle  10 */                                        \
     _x4##r = _b4r * _S##N2_5i;                             \
     _x4##i = _b4i * _S##N2_5i;                             \
     /* Cycle  11 */                                        \
     cplx_sub( _x2, _b1, _b2);                              \
     /* Cycle  12 */                                        \
     _b0r *= _S##_1_5i;                                     \
     _b0i *= _S##_1_5i;                                     \
     /* Cycle  13 */                                        \
     _b3r = _ar * (-0.25) + _x0##r;                         \
     _b3i = _ai * (-0.25) + _x0##i;                         \
     /* Cycle  14 */                                        \
     cplx_add( _x0, _a, _x0);                               \
     /* Cycle  15 */                                        \
     /* Cycle  16 */                                        \
     _x2##r *= _S##N2_5r;                                   \
     _x2##i *= _S##N2_5r;                                   \
     /* Cycle  17 */                                        \
     cplx_add( _x3, _x3, _b0);                              \
     /* Cycle  18 */                                        \
     cplx_sub( _b4, _b0, _x4);                              \
     /* Cycle  19 */                                        \
     /* Cycle  20 */                                        \
     cplx_add( _x1, _b3, _x2);                              \
     /* Cycle  21 */                                        \
     cplx_sub( _x2, _b3, _x2);                              \
     /* Cycle  22 */                                        \
     *( _pd0##r) = _x0##r;                                  \
     *( _pd0##i) = _x0##i;                                  \
     /* Cycle  23 */                                        \
     /* Cycle  24 */                                        \
     _x4##r = _x1##r + _x3##i;                              \
     _x4##i = _x1##i - _x3##r;                              \
     /* Cycle  25 */                                        \
     _x1##r -= _x3##i;                                      \
     _x1##i += _x3##r;                                      \
     /* Cycle  26 */                                        \
     _x3##r = _x2##r + _b4i;                                \
     _x3##i = _x2##i - _b4r;                                \
     /* Cycle  27 */                                        \
     _x2##r -= _b4i;                                        \
     _x2##i += _b4r;                                        \
     /* Cycle  28 */                                        \
     *( _pd4##r) = _x4##r;                                  \
     *( _pd4##i) = _x4##i;                                  \
     /* Cycle  29 */                                        \
     *( _pd1##r) = _x1##r;                                  \
     *( _pd1##i) = _x1##i;                                  \
     /* Cycle  30 */                                        \
     *( _pd3##r) = _x3##r;                                  \
     *( _pd3##i) = _x3##i;                                  \
     /* Cycle  31 */                                        \
     *( _pd2##r) = _x2##r;                                  \
     *( _pd2##i) = _x2##i;                                  \
}

/*$Id$*/



