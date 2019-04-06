/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2003-2006  Guillermo Ballester Valor
 
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

#define sse2_fft7(_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
  sse2_fft7_##_S( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

#define sse2_fft7_store_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
 sse2_fft7_store_p_##_S(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

#define sse2_fft7_store_inter_p(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_S, _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
 sse2_fft7_store_inter_p_##_S(_index,_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6)

#define sse2_fft7_F( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_SUB_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_SUB_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_ADD_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_ADD_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_SUB_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_SUB_PD (_x4##i, _a2i, _a4r);                          \
}

#define sse2_fft7_B( _x0 , _x1 , _x2 , _x3,  _x4 , _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_ADD_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_ADD_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_SUB_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_SUB_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_ADD_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_ADD_PD (_x4##i, _a2i, _a4r);                          \
}

#define sse2_fft7_store_p_F(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _pd5, _pd6, _x0, _x1, _x2, _x3, _x4, _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);      \
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_SUB_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_SUB_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd6 + _index, _x6##r, _x6##i);      \
  Y_MM_ADD_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_ADD_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd5 + _index, _x5##r, _x5##i);      \
  Y_MM_SUB_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_SUB_PD (_x4##i, _a2i, _a4r);                          \
  Y_MM_STORE_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);      \
}

#define sse2_fft7_store_p_B(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _pd5, _pd6, _x0, _x1, _x2, _x3, _x4, _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);      \
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_ADD_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_ADD_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd6 + _index, _x6##r, _x6##i);      \
  Y_MM_SUB_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_SUB_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd5 + _index, _x5##r, _x5##i);      \
  Y_MM_ADD_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_ADD_PD (_x4##i, _a2i, _a4r);                          \
  Y_MM_STORE_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);      \
  Y_MM_STORE_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);      \
}

#define sse2_fft7_store_inter_p_F(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _pd5, _pd6, _x0, _x1, _x2, _x3, _x4, _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_SUB_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_ADD_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_SUB_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd6 + _index, _x6##r, _x6##i);\
  Y_MM_ADD_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_SUB_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_ADD_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd5 + _index, _x5##r, _x5##i);\
  Y_MM_SUB_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_ADD_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_SUB_PD (_x4##i, _a2i, _a4r);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
}

#define sse2_fft7_store_inter_p_B(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _pd5, _pd6, _x0, _x1, _x2, _x3, _x4, _x5, _x6) \
{                                                            \
  Y__M128D _a0r, _a0i, _a1r, _a1i, _a2r, _a2i, _a3r, _a3i,   \
           _a4r, _a4i, _a5r, _a5i, _a6r, _a6i, _ar, _ai;     \
  Y_MM_MOV_PD (_a0r, _x0##r);                                \
  Y_MM_MOV_PD (_a0i, _x0##i);                                \
  sse2_add (_a1, _x1, _x6);                                  \
  sse2_sub (_a6, _x1, _x6);                                  \
  sse2_add (_a2, _x2, _x5);                                  \
  sse2_sub (_a5, _x2, _x5);                                  \
  sse2_add (_a3, _x3, _x4);                                  \
  sse2_add (_a, _a1, _a2);                                   \
  sse2_sub (_x2, _a2, _a1);                                  \
  sse2_sub (_a4, _x3, _x4);                                  \
  sse2_add ( _a, _a, _a3);                                   \
  sse2_sub ( _x1, _a1, _a3);                                 \
  sse2_sub ( _x3, _a3, _a2);                                 \
  sse2_add ( _x0, _x0, _a);                                  \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
  Y_MM_MUL_PD (_ar, _ar, MM_FN1_7r);                         \
  Y_MM_MUL_PD (_ai, _ai, MM_FN1_7r);                         \
  Y_MM_ADD_PD (_ar, _ar, _a0r);                              \
  Y_MM_ADD_PD (_ai, _ai, _a0i);                              \
  sse2_add (_x4, _a5, _a6);                                  \
  Y_MM_MUL_PD (_a0r, _x2##r, MM_FN4_7r);                     \
  Y_MM_MUL_PD (_a0i, _x2##i, MM_FN4_7r);                     \
  sse2_sub (_x6, _a6, _a5);                                  \
  Y_MM_MUL_PD (_a1r, _x1##r, MM_FN2_7r);                     \
  Y_MM_MUL_PD (_a1i, _x1##i, MM_FN2_7r);                     \
  sse2_sub (_x4, _a4, _x4);                                  \
  Y_MM_MUL_PD (_a2r, _x3##r, MM_FN3_7r);                     \
  Y_MM_MUL_PD (_a2i, _x3##i, MM_FN3_7r);                     \
  sse2_add (_x5, _a5, _a4);                                  \
  sse2_add (_a6, _a6, _a4);                                  \
  sse2_add (_x1, _a1, _a0);                                  \
  Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN1_7i);                   \
  Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN1_7i);                   \
  sse2_add (_x3, _a1, _a2);                                  \
  sse2_sub (_x2, _a0, _a2);                                  \
  Y_MM_MUL_PD (_a0r, _x6##r, MM_FN4_7i);                     \
  Y_MM_MUL_PD (_a0i, _x6##i, MM_FN4_7i);                     \
  sse2_sub(_a1, _a, _x1);                                    \
  Y_MM_MUL_PD (_a5r, _x5##r, MM_FN3_7i);                     \
  Y_MM_MUL_PD (_a5i, _x5##i, MM_FN3_7i);                     \
  sse2_add (_a3, _a, _x3);                                   \
  Y_MM_MUL_PD (_a6r, _a6r, MM_FN2_7i);                       \
  Y_MM_MUL_PD (_a6i, _a6i, MM_FN2_7i);                       \
  sse2_add (_a2, _a, _x2);                                   \
  sse2_sub (_x1, _a0, _a5);                                  \
  sse2_sub (_a5, _a5, _a6);                                  \
  sse2_sub (_a6, _a6, _a0);                                  \
  sse2_add (_a4, _x4, _x1);                                  \
  sse2_add (_a5, _a5, _x4);                                  \
  sse2_add (_a6, _a6, _x4);                                  \
	                                                     \
  Y_MM_ADD_PD (_x1##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x6##i, _a3i, _a5r);                          \
  Y_MM_SUB_PD (_x1##r, _a3r, _a5i);                          \
  Y_MM_ADD_PD (_x6##r, _a3r, _a5i);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd6 + _index, _x6##r, _x6##i);\
  Y_MM_SUB_PD (_x2##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x5##r, _a1r, _a6i);                          \
  Y_MM_ADD_PD (_x2##i, _a1i, _a6r);                          \
  Y_MM_SUB_PD (_x5##i, _a1i, _a6r);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd5 + _index, _x5##r, _x5##i);\
  Y_MM_ADD_PD (_x3##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x4##r, _a2r, _a4i);                          \
  Y_MM_SUB_PD (_x3##i, _a2i, _a4r);                          \
  Y_MM_ADD_PD (_x4##i, _a2i, _a4r);                          \
  Y_MM_STORE_INTER_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
  Y_MM_STORE_INTER_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
}

/*$Id$*/






