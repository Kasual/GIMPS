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

/*
   Nussbaumer version for FFT-5. Code adapted from E. Mayer's (Mlucas) 
   to SSE2
*/


#define sse2_fft5(_S, _x0, _x1, _x2, _x3, _x4) \
  sse2_fft5_##_S(_x0, _x1, _x2, _x3, _x4);

#define sse2_fft5_store_p(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _S, _x0, _x1, _x2, _x3, _x4) \
  sse2_fft5_store_p_##_S(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4);

#define sse2_fft5_store_inter_p(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _S, _x0, _x1, _x2, _x3, _x4) \
  sse2_fft5_store_inter_p_##_S(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4);


#define  sse2_fft5_F(_x0, _x1, _x2, _x3, _x4)   \
{                                               \
     Y__M128D _ar, _ai, _br, _bi;               \
     sse2_add (_a, _x1, _x4);                   \
     sse2_add (_b, _x2, _x3);                   \
     sse2_sub (_x4, _x1, _x4);                  \
     sse2_sub (_x3, _x2, _x3);                  \
     sse2_sub (_x2, _a, _b);                    \
     sse2_add (_a, _a, _b);                     \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);   \
     sse2_sub (_b, _x4, _x3);                   \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);   \
     sse2_add (_x0, _x0, _a);                   \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);         \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);         \
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);         \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);         \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);            \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);            \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);   \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);   \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);   \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);   \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);         \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);         \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);         \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);         \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);         \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);         \
     Y_MM_SUB_PD (_x4##r, _br, _x4##r);         \
     Y_MM_SUB_PD (_x4##i, _bi, _x4##i);         \
     Y_MM_MOV_PD (_ar, _x4##r);                 \
     Y_MM_MOV_PD (_ai, _x4##i);                 \
                                                \
     Y_MM_ADD_PD (_x4##r, _x1##r, _x3##i);      \
     Y_MM_SUB_PD (_x4##i, _x1##i, _x3##r);      \
     Y_MM_SUB_PD (_x1##r, _x1##r, _x3##i);      \
     Y_MM_ADD_PD (_x1##i, _x1##i, _x3##r);      \
     Y_MM_ADD_PD (_x3##r, _x2##r, _ai);         \
     Y_MM_SUB_PD (_x3##i, _x2##i, _ar);         \
     Y_MM_SUB_PD (_x2##r, _x2##r, _ai);         \
     Y_MM_ADD_PD (_x2##i, _x2##i, _ar);         \
}

#define  sse2_fft5_B(_x0, _x1, _x2, _x3, _x4)   \
{                                               \
     Y__M128D _ar, _ai, _br, _bi;               \
     sse2_add (_a, _x1, _x4);                   \
     sse2_add (_b, _x2, _x3);                   \
     sse2_sub (_x4, _x1, _x4);                  \
     sse2_sub (_x3, _x2, _x3);                  \
     sse2_sub (_x2, _a, _b);                    \
     sse2_add (_a, _a, _b);                     \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);   \
     sse2_sub (_b, _x4, _x3);                   \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);   \
     sse2_add (_x0, _x0, _a);                   \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);         \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);         \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);            \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);            \
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);         \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);         \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);   \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);   \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);   \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);   \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);         \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);         \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);         \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);         \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);         \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);         \
     Y_MM_SUB_PD (_ar, _br, _x4##r);            \
     Y_MM_SUB_PD (_ai, _bi, _x4##i);            \
                                                \
     Y_MM_SUB_PD (_x4##r, _x1##r, _x3##i);      \
     Y_MM_ADD_PD (_x4##i, _x1##i, _x3##r);      \
     Y_MM_ADD_PD (_x1##r, _x1##r, _x3##i);      \
     Y_MM_SUB_PD (_x1##i, _x1##i, _x3##r);      \
     Y_MM_SUB_PD (_x3##r, _x2##r, _ai);         \
     Y_MM_ADD_PD (_x3##i, _x2##i, _ar);         \
     Y_MM_ADD_PD (_x2##r, _x2##r, _ai);         \
     Y_MM_SUB_PD (_x2##i, _x2##i, _ar);         \
}

#define sse2_fft5_store_p_F(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4) \
{                                                         \
     Y__M128D _ar, _ai, _br, _bi;                         \
     sse2_add (_a, _x1, _x4);                             \
     sse2_add (_b, _x2, _x3);                             \
     sse2_sub (_x4, _x1, _x4);                            \
     sse2_sub (_x3, _x2, _x3);                            \
     sse2_sub (_x2, _a, _b);                              \
     sse2_add (_a, _a, _b);                               \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);             \
     sse2_sub (_b, _x4, _x3);                             \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);             \
     sse2_add (_x0, _x0, _a);                             \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);                   \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);                   \
     Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);                   \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);                   \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);                      \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);                      \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);             \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);             \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);             \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);             \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);                   \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);                   \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);                   \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);                   \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);                   \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);                   \
     Y_MM_SUB_PD (_x4##r, _br, _x4##r);                   \
     Y_MM_SUB_PD (_x4##i, _bi, _x4##i);                   \
     Y_MM_MOV_PD (_ar, _x4##r);                           \
     Y_MM_MOV_PD (_ai, _x4##i);                           \
                                                          \
     Y_MM_ADD_PD (_x4##r, _x1##r, _x3##i);                \
     Y_MM_SUB_PD (_x4##i, _x1##i, _x3##r);                \
     Y_MM_SUB_PD (_x1##r, _x1##r, _x3##i);                \
     Y_MM_ADD_PD (_x1##i, _x1##i, _x3##r);                \
     Y_MM_STORE_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
     Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
     Y_MM_ADD_PD (_x3##r, _x2##r, _ai);                   \
     Y_MM_SUB_PD (_x3##i, _x2##i, _ar);                   \
     Y_MM_SUB_PD (_x2##r, _x2##r, _ai);                   \
     Y_MM_ADD_PD (_x2##i, _x2##i, _ar);                   \
     Y_MM_STORE_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
     Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
}

#define sse2_fft5_store_p_B(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4) \
{                                                         \
     Y__M128D _ar, _ai, _br, _bi;                         \
     sse2_add (_a, _x1, _x4);                             \
     sse2_add (_b, _x2, _x3);                             \
     sse2_sub (_x4, _x1, _x4);                            \
     sse2_sub (_x3, _x2, _x3);                            \
     sse2_sub (_x2, _a, _b);                              \
     sse2_add (_a, _a, _b);                               \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);             \
     sse2_sub (_b, _x4, _x3);                             \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);             \
     sse2_add (_x0, _x0, _a);                             \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);                   \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);                   \
     Y_MM_STORE_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);                   \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);                   \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);                      \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);                      \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);             \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);             \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);             \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);             \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);                   \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);                   \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);                   \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);                   \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);                   \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);                   \
     Y_MM_SUB_PD (_x4##r, _br, _x4##r);                   \
     Y_MM_SUB_PD (_x4##i, _bi, _x4##i);                   \
     Y_MM_MOV_PD (_ar, _x4##r);                           \
     Y_MM_MOV_PD (_ai, _x4##i);                           \
                                                          \
     Y_MM_SUB_PD (_x4##r, _x1##r, _x3##i);                \
     Y_MM_ADD_PD (_x4##i, _x1##i, _x3##r);                \
     Y_MM_ADD_PD (_x1##r, _x1##r, _x3##i);                \
     Y_MM_SUB_PD (_x1##i, _x1##i, _x3##r);                \
     Y_MM_STORE_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
     Y_MM_STORE_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
     Y_MM_SUB_PD (_x3##r, _x2##r, _ai);                   \
     Y_MM_ADD_PD (_x3##i, _x2##i, _ar);                   \
     Y_MM_ADD_PD (_x2##r, _x2##r, _ai);                   \
     Y_MM_SUB_PD (_x2##i, _x2##i, _ar);                   \
     Y_MM_STORE_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
     Y_MM_STORE_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
}

#define sse2_fft5_store_inter_p_F(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4) \
{                                                               \
     Y__M128D _ar, _ai, _br, _bi;                               \
     sse2_add (_a, _x1, _x4);                                   \
     sse2_add (_b, _x2, _x3);                                   \
     sse2_sub (_x4, _x1, _x4);                                  \
     sse2_sub (_x3, _x2, _x3);                                  \
     sse2_sub (_x2, _a, _b);                                    \
     sse2_add (_a, _a, _b);                                     \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);                   \
     sse2_sub (_b, _x4, _x3);                                   \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);                   \
     sse2_add (_x0, _x0, _a);                                   \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);                         \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);                         \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);                         \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);                         \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);                            \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);                            \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);                   \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);                   \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);                   \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);                   \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);                         \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);                         \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);                         \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);                         \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);                         \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);                         \
     Y_MM_SUB_PD (_x4##r, _br, _x4##r);                         \
     Y_MM_SUB_PD (_x4##i, _bi, _x4##i);                         \
     Y_MM_MOV_PD (_ar, _x4##r);                                 \
     Y_MM_MOV_PD (_ai, _x4##i);                                 \
                                                                \
     Y_MM_ADD_PD (_x4##r, _x1##r, _x3##i);                      \
     Y_MM_SUB_PD (_x4##i, _x1##i, _x3##r);                      \
     Y_MM_SUB_PD (_x1##r, _x1##r, _x3##i);                      \
     Y_MM_ADD_PD (_x1##i, _x1##i, _x3##r);                      \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
     Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
     Y_MM_ADD_PD (_x3##r, _x2##r, _ai);                         \
     Y_MM_SUB_PD (_x3##i, _x2##i, _ar);                         \
     Y_MM_SUB_PD (_x2##r, _x2##r, _ai);                         \
     Y_MM_ADD_PD (_x2##i, _x2##i, _ar);                         \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
     Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
}

#define sse2_fft5_store_inter_p_B(_index, _pd0, _pd1, _pd2, _pd3, _pd4, _x0, _x1, _x2, _x3, _x4) \
{                                                               \
     Y__M128D _ar, _ai, _br, _bi;                               \
     sse2_add (_a, _x1, _x4);                                   \
     sse2_add (_b, _x2, _x3);                                   \
     sse2_sub (_x4, _x1, _x4);                                  \
     sse2_sub (_x3, _x2, _x3);                                  \
     sse2_sub (_x2, _a, _b);                                    \
     sse2_add (_a, _a, _b);                                     \
     Y_MM_MUL_PD( _x2##r, _x2##r, MM_FN2_5r);                   \
     sse2_sub (_b, _x4, _x3);                                   \
     Y_MM_MUL_PD( _x2##i, _x2##i, MM_FN2_5r);                   \
     sse2_add (_x0, _x0, _a);                                   \
     Y_MM_MUL_PD (_ar, _ar, MM_FNM125);                         \
     Y_MM_MUL_PD (_ai, _ai, MM_FNM125);                         \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd0 + _index, _x0##r, _x0##i);\
     Y_MM_MUL_PD (_br, _br, MM_F_1_5i);                         \
     Y_MM_MUL_PD (_bi, _bi, MM_F_1_5i);                         \
     Y_MM_ADD_PD (_ar, _ar, _x0##r);                            \
     Y_MM_ADD_PD (_ai, _ai, _x0##i);                            \
     Y_MM_MUL_PD (_x3##r, _x3##r, MM_FN1_5i);                   \
     Y_MM_MUL_PD (_x3##i, _x3##i, MM_FN1_5i);                   \
     Y_MM_MUL_PD (_x4##r, _x4##r, MM_FN2_5i);                   \
     Y_MM_MUL_PD (_x4##i, _x4##i, MM_FN2_5i);                   \
     Y_MM_ADD_PD (_x1##r, _ar, _x2##r);                         \
     Y_MM_ADD_PD (_x1##i, _ai, _x2##i);                         \
     Y_MM_SUB_PD (_x2##r, _ar, _x2##r);                         \
     Y_MM_SUB_PD (_x2##i, _ai, _x2##i);                         \
     Y_MM_ADD_PD (_x3##r, _x3##r, _br);                         \
     Y_MM_ADD_PD (_x3##i, _x3##i, _bi);                         \
     Y_MM_SUB_PD (_x4##r, _br, _x4##r);                         \
     Y_MM_SUB_PD (_x4##i, _bi, _x4##i);                         \
     Y_MM_MOV_PD (_ar, _x4##r);                                 \
     Y_MM_MOV_PD (_ai, _x4##i);                                 \
                                                                \
     Y_MM_SUB_PD (_x4##r, _x1##r, _x3##i);                      \
     Y_MM_ADD_PD (_x4##i, _x1##i, _x3##r);                      \
     Y_MM_ADD_PD (_x1##r, _x1##r, _x3##i);                      \
     Y_MM_SUB_PD (_x1##i, _x1##i, _x3##r);                      \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd4 + _index, _x4##r, _x4##i);\
     Y_MM_STORE_INTER_PAIR_NEXT( _pd1 + _index, _x1##r, _x1##i);\
     Y_MM_SUB_PD (_x3##r, _x2##r, _ai);                         \
     Y_MM_ADD_PD (_x3##i, _x2##i, _ar);                         \
     Y_MM_ADD_PD (_x2##r, _x2##r, _ai);                         \
     Y_MM_SUB_PD (_x2##i, _x2##i, _ar);                         \
     Y_MM_STORE_INTER_PAIR_NEXT( _pd3 + _index, _x3##r, _x3##i);\
     Y_MM_STORE_INTER_PAIR_NEXT( _pd2 + _index, _x2##r, _x2##i);\
}

/*$Id$*/
