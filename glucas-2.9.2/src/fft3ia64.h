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

/*
This code is written with Intel ia64 in mind. The cycles in the comments
are, of course, only aproximated. These schedule pattern is what we would like
the compiler to do. But, sure, the compiler can find better solutions 
 
Every estimated cycle is, by default, using two .mfi boundles
*/


#define cplx_two_fft3_ia64(_S,_r0,_r1,_r2,_s0,_s1,_s2) \
{                                                 \
    y_limb_t a0,a1,b0,b1;                         \
/*  Cycle  1 */                                   \
    a0 = _r1##r;                                  \
    b0 = _r1##i;                                  \
/*  Cycle  2 */                                   \
    a1 = _s1##r;                                  \
    b1 = _s1##i;                                  \
/*  Cycle  3 */                                   \
    _r1##r += _r2##r;                             \
    _r1##i += _r2##i;                             \
/*  Cycle  4 */                                   \
    _s1##r += _s2##r;                             \
    _s1##i += _s2##i;                             \
/*  Cycle  5 */                                   \
    _r2##r = a0 - _r2##r;                         \
    _r2##i = b0 - _r2##i;                         \
/*  Cycle  6 */                                   \
    _s2##r = a1 - _s2##r;                         \
    _s2##i = b1 - _s2##i;                         \
/*  Cycle  7 */                                   \
/*  Cycle  8 */                                   \
   _r0##r += _r1##r;                              \
   _r0##i += _r1##i;                              \
/*  Cycle  9 */                                   \
   _s0##r += _s1##r;                              \
   _s0##i += _s1##i;                              \
/*  Cycle 10 */                                   \
   a0 = _r2##r * _S##_1_3i;                       \
   b0 = _r2##i * _S##_1_3i;                       \
/*  Cycle 11 */                                   \
   a1 = _s2##r * _S##_1_3i;                       \
   b1 = _s2##i * _S##_1_3i;                       \
/*  Cycle 12 */                                   \
/*  Cycle 13 */                                   \
   _r1##r = minus_3_2 * _r1##r + _r0##r;          \
   _r1##i = minus_3_2 * _r1##i + _r0##i;          \
/*  Cycle 14 */                                   \
   _s1##r = minus_3_2 * _s1##r + _s0##r;          \
   _s1##i = minus_3_2 * _s1##i + _s0##i;          \
/*  Cycle 15 */                                   \
/*  Cycle 16 */                                   \
/*  Cycle 17 */                                   \
/*  Cycle 18 */                                   \
   _r2##r = _r1##r + b0;                          \
   _r2##i = _r1##i - a0;                          \
/*  Cycle 19 */                                   \
   _s2##r = _s1##r + b1;                          \
   _s2##i = _s1##i - a1;                          \
/*  Cycle 20 */                                   \
   _r1##r = _r1##r - b0;                          \
   _r1##i = _r1##i + a0;                          \
/*  Cycle 21 */                                   \
   _s1##r = _s1##r - b1;                          \
   _s1##i = _s1##i + a1;                          \
}

#define cplx_three_fft3_ia64_F(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2) \
{                                                               \
    y_limb_t a0,a1,a2,b0,b1,b2;                                 \
/*  Cycle  1 */                                                 \
    a0 = _r1##r;                                                \
    b0 = _r1##i;                                                \
/*  Cycle  2 */                                                 \
    a1 = _s1##r;                                                \
    b1 = _s1##i;                                                \
/*  Cycle  3 */                                                 \
    a2 = _t1##r;                                                \
    b2 = _t1##i;                                                \
/*  Cycle  4 */                                                 \
    _r1##r += _r2##r;                                           \
    _r1##i += _r2##i;                                           \
/*  Cycle  5 */                                                 \
    _s1##r += _s2##r;                                           \
    _s1##i += _s2##i;                                           \
/*  Cycle  5 */                                                 \
    _t1##r += _t2##r;                                           \
    _t1##i += _t2##i;                                           \
/*  Cycle  6 */                                                 \
    _r2##r = a0 - _r2##r;                                       \
    _r2##i = b0 - _r2##i;                                       \
/*  Cycle  7 */                                                 \
    _s2##r = a1 - _s2##r;                                       \
    _s2##i = b1 - _s2##i;                                       \
/*  Cycle  8 */                                                 \
    _t2##r = a2 - _t2##r;                                       \
    _t2##i = b2 - _t2##i;                                       \
/*  Cycle  9 */                                                 \
   _r0##r += _r1##r;                                            \
   _r0##i += _r1##i;                                            \
/*  Cycle 10 */                                                 \
   _s0##r += _s1##r;                                            \
   _s0##i += _s1##i;                                            \
/*  Cycle 11 */                                                 \
   _t0##r += _t1##r;                                            \
   _t0##i += _t1##i;                                            \
/*  Cycle 12 */                                                 \
   a0 = _r2##r * F_1_3i;                                        \
   b0 = _r2##i * F_1_3i;                                        \
/*  Cycle 13 */                                                 \
   a1 = _s2##r * F_1_3i;                                        \
   b1 = _s2##i * F_1_3i;                                        \
/*  Cycle 14 */                                                 \
   a2 = _t2##r * F_1_3i;                                        \
   b2 = _t2##i * F_1_3i;                                        \
/*  Cycle 15 */                                                 \
   _r1##r = minus_3_2 * _r1##r + _r0##r;                        \
   _r1##i = minus_3_2 * _r1##i + _r0##i;                        \
/*  Cycle 16 */                                                 \
   _s1##r = minus_3_2 * _s1##r + _s0##r;                        \
   _s1##i = minus_3_2 * _s1##i + _s0##i;                        \
/*  Cycle 17 */                                                 \
   _t1##r = minus_3_2 * _t1##r + _t0##r;                        \
   _t1##i = minus_3_2 * _t1##i + _t0##i;                        \
/*  Cycle 18 */                                                 \
/*  Cycle 19 */                                                 \
/*  Cycle 20 */                                                 \
   _r2##r = _r1##r + b0;                                        \
   _r2##i = _r1##i - a0;                                        \
/*  Cycle 21 */                                                 \
   _s2##r = _s1##r + b1;                                        \
   _s2##i = _s1##i - a1;                                        \
/*  Cycle 22 */                                                 \
   _t2##r = _t1##r + b2;                                        \
   _t2##i = _t1##i - a2;                                        \
/*  Cycle 23 */                                                 \
   _r1##r = _r1##r - b0;                                        \
   _r1##i = _r1##i + a0;                                        \
/*  Cycle 24 */                                                 \
   _s1##r = _s1##r - b1;                                        \
   _s1##i = _s1##i + a1;                                        \
/*  Cycle 25 */                                                 \
   _t1##r = _t1##r - b2;                                        \
   _t1##i = _t1##i + a2;                                        \
}

#define cplx_three_fft3_ia64_B(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2) \
{                                                               \
    y_limb_t a0,a1,a2,b0,b1,b2;                                 \
/*  Cycle  1 */                                                 \
    a0 = _r1##r;                                                \
    b0 = _r1##i;                                                \
/*  Cycle  2 */                                                 \
    a1 = _s1##r;                                                \
    b1 = _s1##i;                                                \
/*  Cycle  3 */                                                 \
    a2 = _t1##r;                                                \
    b2 = _t1##i;                                                \
/*  Cycle  4 */                                                 \
    _r1##r += _r2##r;                                           \
    _r1##i += _r2##i;                                           \
/*  Cycle  5 */                                                 \
    _s1##r += _s2##r;                                           \
    _s1##i += _s2##i;                                           \
/*  Cycle  5 */                                                 \
    _t1##r += _t2##r;                                           \
    _t1##i += _t2##i;                                           \
/*  Cycle  6 */                                                 \
    _r2##r = a0 - _r2##r;                                       \
    _r2##i = b0 - _r2##i;                                       \
/*  Cycle  7 */                                                 \
    _s2##r = a1 - _s2##r;                                       \
    _s2##i = b1 - _s2##i;                                       \
/*  Cycle  8 */                                                 \
    _t2##r = a2 - _t2##r;                                       \
    _t2##i = b2 - _t2##i;                                       \
/*  Cycle  9 */                                                 \
   _r0##r += _r1##r;                                            \
   _r0##i += _r1##i;                                            \
/*  Cycle 10 */                                                 \
   _s0##r += _s1##r;                                            \
   _s0##i += _s1##i;                                            \
/*  Cycle 11 */                                                 \
   _t0##r += _t1##r;                                            \
   _t0##i += _t1##i;                                            \
/*  Cycle 12 */                                                 \
   a0 = _r2##r * F_1_3i;                                        \
   b0 = _r2##i * F_1_3i;                                        \
/*  Cycle 13 */                                                 \
   a1 = _s2##r * F_1_3i;                                        \
   b1 = _s2##i * F_1_3i;                                        \
/*  Cycle 14 */                                                 \
   a2 = _t2##r * F_1_3i;                                        \
   b2 = _t2##i * F_1_3i;                                        \
/*  Cycle 15 */                                                 \
   _r1##r = minus_3_2 * _r1##r + _r0##r;                        \
   _r1##i = minus_3_2 * _r1##i + _r0##i;                        \
/*  Cycle 16 */                                                 \
   _s1##r = minus_3_2 * _s1##r + _s0##r;                        \
   _s1##i = minus_3_2 * _s1##i + _s0##i;                        \
/*  Cycle 17 */                                                 \
   _t1##r = minus_3_2 * _t1##r + _t0##r;                        \
   _t1##i = minus_3_2 * _t1##i + _t0##i;                        \
/*  Cycle 18 */                                                 \
/*  Cycle 19 */                                                 \
/*  Cycle 20 */                                                 \
   _r2##r = _r1##r - b0;                                        \
   _r2##i = _r1##i + a0;                                        \
/*  Cycle 21 */                                                 \
   _s2##r = _s1##r - b1;                                        \
   _s2##i = _s1##i + a1;                                        \
/*  Cycle 22 */                                                 \
   _t2##r = _t1##r - b2;                                        \
   _t2##i = _t1##i + a2;                                        \
/*  Cycle 23 */                                                 \
   _r1##r = _r1##r + b0;                                        \
   _r1##i = _r1##i - a0;                                        \
/*  Cycle 24 */                                                 \
   _s1##r = _s1##r + b1;                                        \
   _s1##i = _s1##i - a1;                                        \
/*  Cycle 25 */                                                 \
   _t1##r = _t1##r + b2;                                        \
   _t1##i = _t1##i - a2;                                        \
}



#define cplx_three_fft3_first_dif_9_store_F(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2,_pr0,_pr1,_pr2,_ps0,_ps1,_ps2,_pt0,_pt1,_pt2) \
{                                                               \
    y_limb_t a0,a1,a2,b0,b1,b2;                                 \
/*  Cycle  1 */                                                 \
    a0 = _r1##r;                                                \
    b0 = _r1##i;                                                \
/*  Cycle  2 */                                                 \
    a1 = _s1##r;                                                \
    b1 = _s1##i;                                                \
/*  Cycle  3 */                                                 \
    a2 = _t1##r;                                                \
    b2 = _t1##i;                                                \
/*  Cycle  4 */                                                 \
    _r1##r += _r2##r;                                           \
    _r1##i += _r2##i;                                           \
/*  Cycle  5 */                                                 \
    _s1##r += _s2##r;                                           \
    _s1##i += _s2##i;                                           \
/*  Cycle  5 */                                                 \
    _t1##r += _t2##r;                                           \
    _t1##i += _t2##i;                                           \
/*  Cycle  6 */                                                 \
    _r2##r = a0 - _r2##r;                                       \
    _r2##i = b0 - _r2##i;                                       \
/*  Cycle  7 */                                                 \
    _s2##r = a1 - _s2##r;                                       \
    _s2##i = b1 - _s2##i;                                       \
/*  Cycle  8 */                                                 \
    _t2##r = a2 - _t2##r;                                       \
    _t2##i = b2 - _t2##i;                                       \
/*  Cycle  9 */                                                 \
   _r0##r += _r1##r;                                            \
   _r0##i += _r1##i;                                            \
/*  Cycle 10 */                                                 \
   _s0##r += _s1##r;                                            \
   _s0##i += _s1##i;                                            \
/*  Cycle 11 */                                                 \
   _t0##r += _t1##r;                                            \
   _t0##i += _t1##i;                                            \
/*  Cycle 12 */                                                 \
   a0 = _r2##r * F_1_3i;                                        \
   b0 = _r2##i * F_1_3i;                                        \
/*  Cycle 13 */                                                 \
   a1 = _s2##r * F_1_3i;                                        \
   b1 = _s2##i * F_1_3i;                                        \
/*  Cycle 14 */                                                 \
   a2 = _t2##r * F_1_3i;                                        \
   b2 = _t2##i * F_1_3i;                                        \
/*  Cycle 15 */                                                 \
   _r1##r = minus_3_2 * _r1##r + _r0##r;                        \
   _r1##i = minus_3_2 * _r1##i + _r0##i;                        \
/*  Cycle 16 */                                                 \
   _s1##r = minus_3_2 * _s1##r + _s0##r;                        \
   _s1##i = minus_3_2 * _s1##i + _s0##i;                        \
/*  Cycle 17 */                                                 \
   _t1##r = minus_3_2 * _t1##r + _t0##r;                        \
   _t1##i = minus_3_2 * _t1##i + _t0##i;                        \
/*  Cycle 18 */                                                 \
   * (_pr0##r) = _r0##r;                                        \
   * (_pr0##i) = _r0##i;                                        \
/*  Cycle 19 */                                                 \
   * (_ps0##r) = _s0##r;                                        \
   * (_ps0##i) = _s0##i;                                        \
/*  Cycle 20 */                                                 \
   * (_pt0##r) = _t0##r;                                        \
   _r2##r = _r1##r + b0;                                        \
   * (_pt0##i) = _t0##i;                                        \
   _r2##i = _r1##i - a0;                                        \
/*  Cycle 21 */                                                 \
   _s2##r = _s1##r + b1;                                        \
   _s2##i = _s1##i - a1;                                        \
/*  Cycle 22 */                                                 \
   _t2##r = _t1##r + b2;                                        \
   _t2##i = _t1##i - a2;                                        \
/*  Cycle 23 */                                                 \
   _r1##r = _r1##r - b0;                                        \
   _r1##i = _r1##i + a0;                                        \
/*  Cycle 24 */                                                 \
   * (_pr2##r) = _r2##r;                                        \
   _s1##r = _s1##r - b1;                                        \
   * (_pr2##i) = _r2##i;                                        \
   _s1##i = _s1##i + a1;                                        \
/*  Cycle 25 */                                                 \
   * (_ps2##r) = _s2##r;                                        \
   _t1##r = _t1##r - b2;                                        \
   * (_ps2##i) = _s2##i;                                        \
   _t1##i = _t1##i + a2;                                        \
/*  Cycle 26 */                                                 \
   * (_pt2##r) = _t2##r;                                        \
   * (_pt2##i) = _t2##i;                                        \
/*  Cycle 27 */                                                 \
   * (_pr1##r) = _r1##r;                                        \
   * (_pr1##i) = _r1##i;                                        \
/*  Cycle 28 */                                                 \
   * (_ps1##r) = _s1##r;                                        \
   * (_ps1##i) = _s1##i;                                        \
/*  Cycle 29 */                                                 \
   * (_pt1##r) = _t1##r;                                        \
   * (_pt1##i) = _t1##i;                                        \
}


#define cplx_three_fft3_first_dif_9_store_B(_r0,_r1,_r2,_s0,_s1,_s2,_t0,_t1,_t2,_pr0,_pr1,_pr2,_ps0,_ps1,_ps2,_pt0,_pt1,_pt2) \
{                                                               \
    y_limb_t a0,a1,a2,b0,b1,b2;                                 \
/*  Cycle  1 */                                                 \
    a0 = _r1##r;                                                \
    b0 = _r1##i;                                                \
/*  Cycle  2 */                                                 \
    a1 = _s1##r;                                                \
    b1 = _s1##i;                                                \
/*  Cycle  3 */                                                 \
    a2 = _t1##r;                                                \
    b2 = _t1##i;                                                \
/*  Cycle  4 */                                                 \
    _r1##r += _r2##r;                                           \
    _r1##i += _r2##i;                                           \
/*  Cycle  5 */                                                 \
    _s1##r += _s2##r;                                           \
    _s1##i += _s2##i;                                           \
/*  Cycle  5 */                                                 \
    _t1##r += _t2##r;                                           \
    _t1##i += _t2##i;                                           \
/*  Cycle  6 */                                                 \
    _r2##r = a0 - _r2##r;                                       \
    _r2##i = b0 - _r2##i;                                       \
/*  Cycle  7 */                                                 \
    _s2##r = a1 - _s2##r;                                       \
    _s2##i = b1 - _s2##i;                                       \
/*  Cycle  8 */                                                 \
    _t2##r = a2 - _t2##r;                                       \
    _t2##i = b2 - _t2##i;                                       \
/*  Cycle  9 */                                                 \
   _r0##r += _r1##r;                                            \
   _r0##i += _r1##i;                                            \
/*  Cycle 10 */                                                 \
   _s0##r += _s1##r;                                            \
   _s0##i += _s1##i;                                            \
/*  Cycle 11 */                                                 \
   _t0##r += _t1##r;                                            \
   _t0##i += _t1##i;                                            \
/*  Cycle 12 */                                                 \
   a0 = _r2##r * F_1_3i;                                        \
   b0 = _r2##i * F_1_3i;                                        \
/*  Cycle 13 */                                                 \
   a1 = _s2##r * F_1_3i;                                        \
   b1 = _s2##i * F_1_3i;                                        \
/*  Cycle 14 */                                                 \
   a2 = _t2##r * F_1_3i;                                        \
   b2 = _t2##i * F_1_3i;                                        \
/*  Cycle 15 */                                                 \
   _r1##r = minus_3_2 * _r1##r + _r0##r;                        \
   _r1##i = minus_3_2 * _r1##i + _r0##i;                        \
/*  Cycle 16 */                                                 \
   _s1##r = minus_3_2 * _s1##r + _s0##r;                        \
   _s1##i = minus_3_2 * _s1##i + _s0##i;                        \
/*  Cycle 17 */                                                 \
   _t1##r = minus_3_2 * _t1##r + _t0##r;                        \
   _t1##i = minus_3_2 * _t1##i + _t0##i;                        \
/*  Cycle 18 */                                                 \
   * (_pr0##r) = _r0##r;                                        \
   * (_pr0##i) = _r0##i;                                        \
/*  Cycle 19 */                                                 \
   * (_ps0##r) = _s0##r;                                        \
   * (_ps0##i) = _s0##i;                                        \
/*  Cycle 20 */                                                 \
   * (_pt0##r) = _t0##r;                                        \
   _r2##r = _r1##r - b0;                                        \
   * (_pt0##i) = _t0##i;                                        \
   _r2##i = _r1##i + a0;                                        \
/*  Cycle 21 */                                                 \
   _s2##r = _s1##r - b1;                                        \
   _s2##i = _s1##i + a1;                                        \
/*  Cycle 22 */                                                 \
   _t2##r = _t1##r - b2;                                        \
   _t2##i = _t1##i + a2;                                        \
/*  Cycle 23 */                                                 \
   _r1##r = _r1##r + b0;                                        \
   _r1##i = _r1##i - a0;                                        \
/*  Cycle 24 */                                                 \
   * (_pr2##r) = _r2##r;                                        \
   _s1##r = _s1##r + b1;                                        \
   * (_pr2##i) = _r2##i;                                        \
   _s1##i = _s1##i - a1;                                        \
/*  Cycle 25 */                                                 \
   * (_ps2##r) = _s2##r;                                        \
   _t1##r = _t1##r + b2;                                        \
   * (_ps2##i) = _s2##i;                                        \
   _t1##i = _t1##i - a2;                                        \
/*  Cycle 26 */                                                 \
   * (_pt2##r) = _t2##r;                                        \
   * (_pt2##i) = _t2##i;                                        \
/*  Cycle 27 */                                                 \
   * (_pr1##r) = _r1##r;                                        \
   * (_pr1##i) = _r1##i;                                        \
/*  Cycle 28 */                                                 \
   * (_ps1##r) = _s1##r;                                        \
   * (_ps1##i) = _s1##i;                                        \
/*  Cycle 29 */                                                 \
   * (_pt1##r) = _t1##r;                                        \
   * (_pt1##i) = _t1##i;                                        \
}



#define cplx_fft6_ia64(_S,_r0,_r1,_r2,_r3,_r4,_r5) \
{                                                 \
    y_limb_t a0,a1,b0,b1,a2,b2;                   \
/*  Cycle  1 */                                   \
    a0 = _r1##r;                                  \
    b0 = _r1##i;                                  \
/*  Cycle  2 */                                   \
    a1 = _r4##r;                                  \
    b1 = _r4##i;                                  \
/*  Cycle  3 */                                   \
    _r1##r += _r2##r;                             \
    _r1##i += _r2##i;                             \
/*  Cycle  4 */                                   \
    _r4##r += _r5##r;                             \
    _r4##i += _r5##i;                             \
/*  Cycle  5 */                                   \
    _r2##r = a0 - _r2##r;                         \
    _r2##i = b0 - _r2##i;                         \
/*  Cycle  6 */                                   \
    _r5##r = a1 - _r5##r;                         \
    _r5##i = b1 - _r5##i;                         \
/*  Cycle  7 */                                   \
/*  Cycle  8 */                                   \
   _r0##r += _r1##r;                              \
   _r0##i += _r1##i;                              \
/*  Cycle  9 */                                   \
   _r3##r += _r4##r;                              \
   _r3##i += _r4##i;                              \
/*  Cycle 10 */                                   \
   a0 = _r2##r * _S##_1_3i;                       \
   b0 = _r2##i * _S##_1_3i;                       \
/*  Cycle 11 */                                   \
   a1 = _r5##r * _S##_1_3i;                       \
   b1 = _r5##i * _S##_1_3i;                       \
/*  Cycle 12 */                                   \
/*  Cycle 13 */                                   \
   _r1##r = minus_3_2 * _r1##r + _r0##r;          \
   _r1##i = minus_3_2 * _r1##i + _r0##i;          \
/*  Cycle 14 */                                   \
   _r4##r = minus_3_2 * _r4##r + _r3##r;          \
   _r4##i = minus_3_2 * _r4##i + _r3##i;          \
/*  Cycle 15 */                                   \
   a2 = _r0##r;                                   \
   b2 = _r0##i;                                   \
/*  Cycle 16 */                                   \
   _r0##r += _r3##r;                              \
   _r0##i += _r3##i;                              \
/*  Cycle 17 */                                   \
/*  Cycle 18 */                                   \
   _r2##r = _r1##r + b0;                          \
   _r2##i = _r1##i - a0;                          \
/*  Cycle 19 */                                   \
   _r5##r = _r4##r + b1;                          \
   _r5##i = _r4##i - a1;                          \
/*  Cycle 20 */                                   \
   _r1##r = _r1##r - b0;                          \
   _r1##i = _r1##i + a0;                          \
/*  Cycle 21 */                                   \
   _r4##r = _r4##r - b1;                          \
   _r4##i = _r4##i + a1;                          \
/*  Cycle 22 */                                   \
   _r3##r = a2 - _r3##r;                          \
   _r3##i = b2 - _r3##i;                          \
/*  Cycle 23 */                                   \
/*  Cycle 24 */                                   \
   a0 = _r5##i * _S##_1_3i;                       \
   b0 = _r5##r * _S##_1_3i;                       \
/*  Cycle 25 */                                   \
   a1 = _r4##i * _S##_1_3i;                       \
   b1 = _r4##r * _S##_1_3i;                       \
/*  Cycle 26 */                                   \
/*  Cycle 27 */                                   \
/*  Cycle 28 */                                   \
   a0 = _r5##r * F_1_3r - a0;                     \
   b0 = _r5##i * F_1_3r + b0;                     \
/*  Cycle 29 */                                   \
   a1 = _r4##r * F_1_3r + a1;                     \
   b1 = _r4##i * F_1_3r - b1;                     \
/*  Cycle 30 */                                   \
/*  Cycle 31 */                                   \
/*  Cycle 32 */                                   \
   _r5##r = _r2##r - a0;                          \
   _r5##i = _r2##i - b0;                          \
/*  Cycle 33 */                                   \
   _r4##r = _r1##r + a1;                          \
   _r4##i = _r1##i + b1;                          \
/*  Cycle 34 */                                   \
   _r2##r += a0;                                  \
   _r2##i += b0;                                  \
/*  Cycle 35 */                                   \
   _r1##r -= a1;                                  \
   _r1##i -= b1;                                  \
}

#define cplx_fft6_ia64_store(_S,_r0,_r1,_r2,_r3,_r4,_r5,_p0,_p1,_p2,_p3,_p4,_p5) \
{                                                 \
    y_limb_t a0,a1,b0,b1,a2,b2;                   \
/*  Cycle  1 */                                   \
    a0 = _r1##r;                                  \
    b0 = _r1##i;                                  \
/*  Cycle  2 */                                   \
    a1 = _r4##r;                                  \
    b1 = _r4##i;                                  \
/*  Cycle  3 */                                   \
    _r1##r += _r2##r;                             \
    _r1##i += _r2##i;                             \
/*  Cycle  4 */                                   \
    _r4##r += _r5##r;                             \
    _r4##i += _r5##i;                             \
/*  Cycle  5 */                                   \
    _r2##r = a0 - _r2##r;                         \
    _r2##i = b0 - _r2##i;                         \
/*  Cycle  6 */                                   \
    _r5##r = a1 - _r5##r;                         \
    _r5##i = b1 - _r5##i;                         \
/*  Cycle  7 */                                   \
/*  Cycle  8 */                                   \
   _r0##r += _r1##r;                              \
   _r0##i += _r1##i;                              \
/*  Cycle  9 */                                   \
   _r3##r += _r4##r;                              \
   _r3##i += _r4##i;                              \
/*  Cycle 10 */                                   \
   a0 = _r2##r * _S##_1_3i;                       \
   b0 = _r2##i * _S##_1_3i;                       \
/*  Cycle 11 */                                   \
   a1 = _r5##r * _S##_1_3i;                       \
   b1 = _r5##i * _S##_1_3i;                       \
/*  Cycle 12 */                                   \
/*  Cycle 13 */                                   \
   _r1##r = minus_3_2 * _r1##r + _r0##r;          \
   _r1##i = minus_3_2 * _r1##i + _r0##i;          \
/*  Cycle 14 */                                   \
   _r4##r = minus_3_2 * _r4##r + _r3##r;          \
   _r4##i = minus_3_2 * _r4##i + _r3##i;          \
/*  Cycle 15 */                                   \
   a2 = _r0##r;                                   \
   b2 = _r0##i;                                   \
/*  Cycle 16 */                                   \
   _r0##r += _r3##r;                              \
   _r0##i += _r3##i;                              \
/*  Cycle 17 */                                   \
/*  Cycle 18 */                                   \
   _r2##r = _r1##r + b0;                          \
   _r2##i = _r1##i - a0;                          \
/*  Cycle 19 */                                   \
   _r5##r = _r4##r + b1;                          \
   _r5##i = _r4##i - a1;                          \
/*  Cycle 20 */                                   \
   _r1##r = _r1##r - b0;                          \
   _r1##i = _r1##i + a0;                          \
/*  Cycle 21 */                                   \
   _r4##r = _r4##r - b1;                          \
   _r4##i = _r4##i + a1;                          \
/*  Cycle 22 */                                   \
   _r3##r = a2 - _r3##r;                          \
   _r3##i = b2 - _r3##i;                          \
/*  Cycle 23 */                                   \
/*  Cycle 24 */                                   \
   a0 = _r5##i * _S##_1_3i;                       \
   b0 = _r5##r * _S##_1_3i;                       \
/*  Cycle 25 */                                   \
   a1 = _r4##i * _S##_1_3i;                       \
   b1 = _r4##r * _S##_1_3i;                       \
/*  Cycle 26 */                                   \
   *(_p0##r) = _r0##r;                            \
   *(_p0##i) = _r0##i;                            \
/*  Cycle 27 */                                   \
   *(_p3##r) = _r3##r;                            \
   *(_p3##i) = _r3##i;                            \
/*  Cycle 28 */                                   \
   a0 = _r5##r * F_1_3r - a0;                     \
   b0 = _r5##i * F_1_3r + b0;                     \
/*  Cycle 29 */                                   \
   a1 = _r4##r * F_1_3r + a1;                     \
   b1 = _r4##i * F_1_3r - b1;                     \
/*  Cycle 30 */                                   \
/*  Cycle 31 */                                   \
/*  Cycle 32 */                                   \
   _r5##r = _r2##r - a0;                          \
   _r5##i = _r2##i - b0;                          \
/*  Cycle 33 */                                   \
   _r4##r = _r1##r + a1;                          \
   _r4##i = _r1##i + b1;                          \
/*  Cycle 34 */                                   \
   _r2##r += a0;                                  \
   _r2##i += b0;                                  \
/*  Cycle 35 */                                   \
   _r1##r -= a1;                                  \
   _r1##i -= b1;                                  \
/*  Cycle 36 */                                   \
   *(_p5##r) = _r5##r;                            \
   *(_p5##i) = _r5##i;                            \
/*  Cycle 37 */                                   \
   *(_p4##r) = _r4##r;                            \
   *(_p4##i) = _r4##i;                            \
/*  Cycle 38 */                                   \
   *(_p2##r) = _r2##r;                            \
   *(_p2##i) = _r2##i;                            \
/*  Cycle 39 */                                   \
   *(_p1##r) = _r1##r;                            \
   *(_p1##i) = _r1##i;                            \
}

/*$Id$*/





