/*$Id$*/
/*  This file is a part of
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
   This the subroutines are related with dyadic mul and square 
   The main difference among difdit_*.c are the use of subroutines 
   instead of macros
   Thus the size is reduced
*/
#include <stdio.h>
#include <stdlib.h>
#include "yeafft.h"
#include "generic.h"
#include "ydebug.h"
#include "prefetch.h"
#define NDEBUG1


/* asign the correct values to every pad */

#if (Y_AVAL > 4)

# if defined(Y_MINIMUM)
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 2



# define cplx_trig_min_load(_t0r,_t0i,_p) \
    _t0r = *(_p);\
    _t0i = *(_p + 1);\
    prefetch_data_trig( _p , 2);

# define cplx_trig_min_square(_tr,_ti,_sr,_si)\
    _tr = (_sr - _si)*(_sr + _si);\
    _ti = 2.0 * _sr * _si ;

# define cplx_trig_min_mul(_tr,_ti,_s0r,_s0i,_s1r,_s1i)\
    _tr = (_s0r * _s1r) - (_s0i * _s1i);\
    _ti = (_s0r * _s1i) + (_s0i * _s1r);

# define cplx_trig_min_div(_tr,_ti,_s0r,_s0i,_s1r,_s1i)\
    _tr = (_s0r * _s1r) + (_s0i * _s1i);\
    _ti = (_s0i * _s1r) - (_s0r * _s1i);

# define cplx_trig_min_muldiv(_tmr,_tmi,_tdr,_tdi,_s0r,_s0i,_s1r,_s1i)\
    _tmr = (_s0r * _s1r) - (_s0i * _s1i);\
    _tdr = (_s0r * _s1r) + (_s0i * _s1i);\
    _tmi = (_s0r * _s1i) + (_s0i * _s1r);\
    _tdi = (_s0i * _s1r) - (_s0r * _s1i);

void fget_twiddle_factors_32_last(struct y_complex * tw, y_size_t _j,
                                  y_size_t inc)
{
  y_ptr px;
  px = &(Y_TWDF[Y_NRADICES-2][( _j * inc) << 1 ]);
  cplx_trig_min_load(tw[1].r,tw[1].i,px);
  cplx_trig_min_square(tw[2].r,tw[2].i,tw[1].r,tw[1].i);
  cplx_trig_min_square(tw[4].r,tw[4].i,tw[2].r,tw[2].i);
  cplx_trig_min_square(tw[8].r,tw[8].i,tw[4].r,tw[4].i);
  cplx_divmul(tw[3].r,tw[3].i,tw[5].r,tw[5].i,tw[4].r,tw[4].i,tw[1].r,tw[1].i);
  cplx_divmul(tw[7].r,tw[7].i,tw[9].r,tw[9].i,tw[8].r,tw[8].i,tw[1].r,tw[1].i);
  cplx_divmul(tw[6].r,tw[6].i,tw[10].r,tw[10].i,tw[8].r,tw[8].i,tw[2].r,tw[2].i);
  cplx_trig_min_square(tw[16].r,tw[16].i,tw[8].r,tw[8].i);
  cplx_divmul(tw[15].r,tw[15].i,tw[17].r,tw[17].i,tw[16].r,tw[16].i,tw[1].r,tw[1].i);
  cplx_divmul(tw[14].r,tw[14].i,tw[18].r,tw[18].i,tw[16].r,tw[16].i,tw[2].r,tw[2].i);
  cplx_divmul(tw[13].r,tw[13].i,tw[19].r,tw[19].i,tw[16].r,tw[16].i,tw[3].r,tw[3].i);
  cplx_divmul(tw[12].r,tw[12].i,tw[20].r,tw[20].i,tw[16].r,tw[16].i,tw[4].r,tw[4].i);
  cplx_divmul(tw[11].r,tw[11].i,tw[21].r,tw[21].i,tw[16].r,tw[16].i,tw[5].r,tw[5].i);
  cplx_trig_min_square(tw[26].r,tw[26].i,tw[13].r,tw[13].i);
  cplx_divmul(tw[25].r,tw[25].i,tw[27].r,tw[27].i,tw[26].r,tw[26].i,tw[1].r,tw[1].i);
  cplx_divmul(tw[24].r,tw[24].i,tw[28].r,tw[28].i,tw[26].r,tw[26].i,tw[2].r,tw[2].i);
  cplx_divmul(tw[23].r,tw[23].i,tw[29].r,tw[29].i,tw[26].r,tw[26].i,tw[3].r,tw[3].i);
  cplx_divmul(tw[22].r,tw[22].i,tw[30].r,tw[30].i,tw[26].r,tw[26].i,tw[4].r,tw[4].i);
  cplx_trig_min_mul(tw[31].r,tw[31].i,tw[16].r,tw[16].i,tw[15].r,tw[15].i);
}

# elif defined(Y_MAXIMUM)
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 62

void fget_twiddle_factors_32_last(struct y_complex * tw, y_size_t j,
                                  y_size_t inc)
{
  y_ptr px;
  ASSERT_ALIGNED_DOUBLE();
  px = &(Y_TWDF[Y_NRADICES-2][ (j * inc)<<1 ] );
  tw[1].r= px[0];
  tw[1].i= px[1];
  tw[2].r= px[2];
  tw[2].i= px[3];
  tw[3].r= px[4];
  tw[3].i= px[5];
  tw[4].r= px[6];
  tw[4].i= px[7];
  tw[5].r= px[8];
  tw[5].i= px[9];
  tw[6].r= px[10];
  tw[6].i= px[11];
  tw[7].r= px[12];
  tw[7].i= px[13];
  tw[8].r= px[14];
  tw[8].i= px[15];
  tw[9].r= px[16];
  tw[9].i= px[17];
  tw[10].r= px[18];
  tw[10].i= px[19];
  tw[11].r= px[20];
  tw[11].i= px[21];
  tw[12].r= px[22];
  tw[12].i= px[23];
  tw[13].r= px[24];
  tw[13].i= px[25];
  tw[14].r= px[26];
  tw[14].i= px[27];
  tw[15].r= px[28];
  tw[15].i= px[29];
  prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 2*Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 3*Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 4*Y_CACHE_LINE);
  tw[16].r= px[30];
  tw[16].i= px[31];
  tw[17].r= px[32];
  tw[17].i= px[33];
  tw[18].r= px[34];
  tw[18].i= px[35];
  tw[19].r= px[36];
  tw[19].i= px[37];
  tw[20].r= px[38];
  tw[20].i= px[39];
  tw[21].r= px[40];
  tw[21].i= px[41];
  tw[22].r= px[42];
  tw[22].i= px[43];
  tw[23].r= px[44];
  tw[23].i= px[45];
  tw[24].r= px[46];
  tw[24].i= px[47];
  tw[25].r= px[48];
  tw[25].i= px[49];
  tw[26].r= px[50];
  tw[26].i= px[51];
  tw[27].r= px[52];
  tw[27].i= px[53];
  tw[28].r= px[54];
  tw[28].i= px[55];
  tw[29].r= px[56];
  tw[29].i= px[57];
  tw[30].r= px[58];
  tw[30].i= px[59];
  tw[31].r= px[60];
  tw[31].i= px[61];
  prefetch_data_trig(px, Y_STEP + 5*Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 6*Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 7*Y_CACHE_LINE);
  prefetch_data_trig(px, Y_STEP + 8*Y_CACHE_LINE);
}

# else
#  ifdef Y_STEP
#   undef Y_STEP
#  endif
#  define Y_STEP 18

void fget_twiddle_factors_32_last(struct y_complex * tw, y_size_t j,
                                  y_size_t inc)
{
  y_ptr px;
  ASSERT_ALIGNED_DOUBLE();
  px = &(Y_TWDF[Y_NRADICES-2][ (j * inc)<<1 ] );

  tw[1].r= px[0];
  tw[1].i= px[1];
  tw[2].r= px[2];
  tw[2].i= px[3];
  tw[3].r= px[4];
  tw[3].i= px[5];
  tw[4].r= px[6];
  tw[4].i= px[7];
  tw[5].r= px[8];
  tw[5].i= px[9];
  tw[6].r= px[10];
  tw[6].i= px[11];
  tw[7].r= px[12];
  tw[7].i= px[13];
  tw[15].r= px[14];
  tw[15].i= px[15];
  cplx_divmul(tw[8].r,tw[8].i,tw[22].r,tw[22].i,tw[15].r,tw[15].i,tw[7].r,tw[7].i);
  cplx_divmul(tw[9].r,tw[9].i,tw[21].r,tw[21].i,tw[15].r,tw[15].i,tw[6].r,tw[6].i);
  cplx_divmul(tw[10].r,tw[10].i,tw[20].r,tw[20].i,tw[15].r,tw[15].i,tw[5].r,tw[5].i);
  cplx_divmul(tw[11].r,tw[11].i,tw[19].r,tw[19].i,tw[15].r,tw[15].i,tw[4].r,tw[4].i);
  cplx_divmul(tw[12].r,tw[12].i,tw[18].r,tw[18].i,tw[15].r,tw[15].i,tw[3].r,tw[3].i);
  cplx_divmul(tw[13].r,tw[13].i,tw[17].r,tw[17].i,tw[15].r,tw[15].i,tw[2].r,tw[2].i);
  cplx_divmul(tw[14].r,tw[14].i,tw[16].r,tw[16].i,tw[15].r,tw[15].i,tw[1].r,tw[1].i);
  tw[27].r= px[16];
  tw[27].i= px[17];
  cplx_divmul(tw[23].r,tw[23].i,tw[31].r,tw[31].i,tw[27].r,tw[27].i,tw[4].r,tw[4].i);
  cplx_divmul(tw[24].r,tw[24].i,tw[30].r,tw[30].i,tw[27].r,tw[27].i,tw[3].r,tw[3].i);
  cplx_divmul(tw[25].r,tw[25].i,tw[29].r,tw[29].i,tw[27].r,tw[27].i,tw[2].r,tw[2].i);
  cplx_divmul(tw[26].r,tw[26].i,tw[28].r,tw[28].i,tw[27].r,tw[27].i,tw[1].r,tw[1].i);
}

# endif
#endif

#define ysquare_nested_c(_xk,_xmk,_i,_j) \
{\
   y_limb_t _r0,_r1,_r2,_r3,_r4,_r5,_r6,_r7,_gkr,_gki;\
   _gkr=(tw[1].r * Y_DYADIC[ _i << 1]) - (tw[1].i * Y_DYADIC[(_i << 1)+1]) +1.0;\
   _gki=(tw[1].r * Y_DYADIC[ (_i << 1)+1 ]) + (tw[1].i * Y_DYADIC[_i << 1]);\
   _r0=(_xk[_i].r - _xmk[_j].r)*0.5;\
   _r1=(_xk[_i].i + _xmk[_j].i)*0.5;\
   _r2=(_xk[_i].r + _xk[_i].i) * (_xk[_i].r - _xk[_i].i);\
   _r3=(_xk[_i].i * _xk[_i].r);\
   _r4 =(_r0 + _r1) * (_r0 - _r1);\
   _r5 = _r0 * _r1;\
   _r6 =(_xmk[_j].r + _xmk[_j].i) * (_xmk[_j].r - _xmk[_j].i);\
   _r7 = _xmk[_j].i * _xmk[_j].r;\
   _r3 += _r3;\
   _r5 += _r5;\
   _r7 += _r7;\
   _r0= (_gkr * _r4 - _gki * _r5);\
   _r1= (_gkr * _r5 + _gki * _r4);\
   _xk[_i].r = _r2 - _r0;\
   _xmk[_j].r = _r6 - _r0;\
   _xk[_i].i = _r3 - _r1;\
   _xmk[_j].i = _r7 + _r1;\
}


#define ysquare_nested_eq_c(_xk,_i) \
{\
   y_limb_t _y0r=_xk[_i].r, _y0i=_xk[_i].i, _gkr;\
   _gkr=(tw[1].r * Y_DYADIC[ _i << 1]) - (tw[1].i * Y_DYADIC[(_i << 1)+1]);\
   _xk[_i].r =(_y0r * _y0r + _gkr * (_y0i * _y0i));\
   _xk[_i].i = 2.0 * _y0r * _y0i;\
}

/***********************************************************************
 This is the dyadic mul for nested complex representation 
 xk =(xk + xmk* ) * (yk + ymk* ) + 2*((xk * yk) - (xmk* * ymk* ) -
     G^(-k) * (xk - xmk* ) * (yk - ymk* );
 xkm =(xkm + xk* ) * (ymk + yk* ) + 2*((xmk * xk*) - (ymk * yk* ) -
     G^(k) * (xmk - xk* ) * (ymk - yk* );
 where all are complex vars. 
***********************************************************************/
/* This is the dyadic mul for nested complex representation */
#define yconv_nested_c(_xk,_xmk,_yk,_ymk,_i,_j) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i,_gkr,_gki;\
   _gkr=(tw[1].r * Y_DYADIC[ _i << 1]) - (tw[1].i * Y_DYADIC[(_i << 1)+1]) + 1.0;\
   _gki=(tw[1].r * Y_DYADIC[ (_i << 1)+1 ]) + (tw[1].i * Y_DYADIC[_i << 1]);\
   _y0r=(_xk[_i].r - _xmk[_j].r);\
   _y0i=(_xk[_i].i + _xmk[_j].i);\
   _y1r=(_yk[_i].r - _ymk[_j].r);\
   _y1i=(_yk[_i].i + _ymk[_j].i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(_xk[_i].r * _yk[_i].r) - (_xk[_i].i * _yk[_i].i);\
   _y0i=(_xk[_i].r * _yk[_i].i) + (_xk[_i].i * _yk[_i].r);\
   _y1r=(_xmk[_j].r * _ymk[_j].r) - (_xmk[_j].i * _ymk[_j].i);\
   _y1i=(_xmk[_j].r * _ymk[_j].i) + (_xmk[_j].i * _ymk[_j].r);\
   _y3r= (_gkr * _y2r - _gki * _y2i)*0.25;\
   _y3i= (_gkr * _y2i + _gki * _y2r)*0.25;\
   _xk[_i].r = _y0r - _y3r;\
   _xk[_i].i = _y0i - _y3i;\
   _xmk[_j].r = _y1r - _y3r;\
   _xmk[_j].i = _y1i + _y3i;\
}

#if (Y_AVAL > 4)

y_ptr fradix_32_twd_last(struct y_complex * _z, y_ptr  _d,
                         y_size_t _j, struct y_complex  * tw )
{
  y_ptr pd;
  ASSERT_ALIGNED_DOUBLE();
# ifndef NDEBUG1

  printf(" radix_32_twd_last \n");
# endif

  pd = _d + addr((_j)<<1);
  cplx_load_muladdsub_pp(_z[0].r,_z[0].i,_z[1].r,_z[1].i,tw[16].r,tw[16].i, 0 , pd , pd + 32);
  cplx_load_mulmuladdsub_pp(_z[2].r,_z[2].i,_z[3].r,_z[3].i,tw[8].r,tw[8].i,tw[24].r,tw[24].i,0,pd + 16,pd + 48);
  cplx_addsub(_z[0].r,_z[0].i,_z[2].r,_z[2].i);
  cplx_mul_1_4_F_addsub(_z[1].r,_z[1].i,_z[3].r,_z[3].i);
  cplx_load_mulmuladdsub_pp(_z[4].r,_z[4].i,_z[5].r,_z[5].i,tw[4].r,tw[4].i,tw[20].r,tw[20].i,0,pd + 8, pd + 40);
  cplx_load_mulmuladdsub_pp(_z[6].r,_z[6].i,_z[7].r,_z[7].i,tw[12].r,tw[12].i,tw[28].r,tw[28].i,0,pd +24,pd + 56);
  cplx_addsub(_z[4].r,_z[4].i,_z[6].r,_z[6].i);
  cplx_mul_1_4_F_addsub(_z[5].r,_z[5].i,_z[7].r,_z[7].i);
  cplx_addsub(_z[0].r,_z[0].i,_z[4].r,_z[4].i);
  cplx_mul_1_8_F_addsub(_z[1].r,_z[1].i,_z[5].r,_z[5].i);
  cplx_mul_1_4_F_addsub(_z[2].r,_z[2].i,_z[6].r,_z[6].i);
  cplx_mul_3_8_F_addsub(_z[3].r,_z[3].i,_z[7].r,_z[7].i);
  cplx_load_mulmuladdsub_pp(_z[8].r,_z[8].i,_z[9].r,_z[9].i,tw[2].r,tw[2].i,tw[18].r,tw[18].i,0,pd + 4, pd + 36);
  cplx_load_mulmuladdsub_pp(_z[10].r,_z[10].i,_z[11].r,_z[11].i,tw[10].r,tw[10].i,tw[26].r,tw[26].i,0,pd+20,pd+52);
  cplx_addsub(_z[8].r,_z[8].i,_z[10].r,_z[10].i);
  cplx_mul_1_4_F_addsub(_z[9].r,_z[9].i,_z[11].r,_z[11].i);
  cplx_load_mulmuladdsub_pp(_z[12].r,_z[12].i,_z[13].r,_z[13].i,tw[6].r,tw[6].i,tw[22].r,tw[22].i,0,pd+12,pd+44);
  cplx_load_mulmuladdsub_pp(_z[14].r,_z[14].i,_z[15].r,_z[15].i,tw[14].r,tw[14].i,tw[30].r,tw[30].i,0,pd+28,pd+60);
  cplx_addsub(_z[12].r,_z[12].i,_z[14].r,_z[14].i);
  cplx_mul_1_4_F_addsub(_z[13].r,_z[13].i,_z[15].r,_z[15].i);
  cplx_addsub(_z[8].r,_z[8].i,_z[12].r,_z[12].i);
  cplx_mul_1_8_F_addsub(_z[9].r,_z[9].i,_z[13].r,_z[13].i);
  cplx_mul_1_4_F_addsub(_z[10].r,_z[10].i,_z[14].r,_z[14].i);
  cplx_mul_3_8_F_addsub(_z[11].r,_z[11].i,_z[15].r,_z[15].i);
  cplx_addsub(_z[0].r,_z[0].i,_z[8].r,_z[8].i);
  cplx_muladdsub(_z[1].r,_z[1].i,_z[9].r,_z[9].i,F_1_16r,F_1_16i);
  cplx_mul_1_8_F_addsub(_z[2].r,_z[2].i,_z[10].r,_z[10].i);
  cplx_muladdsub(_z[3].r,_z[3].i,_z[11].r,_z[11].i,F_3_16r,F_3_16i);
  cplx_mul_1_4_F_addsub(_z[4].r,_z[4].i,_z[12].r,_z[12].i);
  cplx_muladdsub(_z[5].r,_z[5].i,_z[13].r,_z[13].i,F_5_16r,F_5_16i);
  cplx_mul_3_8_F_addsub(_z[6].r,_z[6].i,_z[14].r,_z[14].i);
  cplx_muladdsub(_z[7].r,_z[7].i,_z[15].r,_z[15].i,F_7_16r,F_7_16i);

  cplx_load_mulmuladdsub_pp(_z[16].r,_z[16].i,_z[17].r,_z[17].i,tw[1].r,tw[1].i,tw[17].r,tw[17].i,0,pd+2,pd+34);
  cplx_load_mulmuladdsub_pp(_z[18].r,_z[18].i,_z[19].r,_z[19].i,tw[9].r,tw[9].i,tw[25].r,tw[25].i,0,pd+18,pd+50);
  cplx_addsub(_z[16].r,_z[16].i,_z[18].r,_z[18].i);
  cplx_mul_1_4_F_addsub(_z[17].r,_z[17].i,_z[19].r,_z[19].i);
  cplx_load_mulmuladdsub_pp(_z[20].r,_z[20].i,_z[21].r,_z[21].i,tw[5].r,tw[5].i,tw[21].r,tw[21].i,0,pd+10,pd+42);
  cplx_load_mulmuladdsub_pp(_z[22].r,_z[22].i,_z[23].r,_z[23].i,tw[13].r,tw[13].i,tw[29].r,tw[29].i,0,pd+26,pd+58);
  cplx_addsub(_z[20].r,_z[20].i,_z[22].r,_z[22].i);
  cplx_mul_1_4_F_addsub(_z[21].r,_z[21].i,_z[23].r,_z[23].i);
  cplx_addsub(_z[16].r,_z[16].i,_z[20].r,_z[20].i);
  cplx_mul_1_8_F_addsub(_z[17].r,_z[17].i,_z[21].r,_z[21].i);
  cplx_mul_1_4_F_addsub(_z[18].r,_z[18].i,_z[22].r,_z[22].i);
  cplx_mul_3_8_F_addsub(_z[19].r,_z[19].i,_z[23].r,_z[23].i);
  cplx_load_mulmuladdsub_pp(_z[24].r,_z[24].i,_z[25].r,_z[25].i,tw[3].r,tw[3].i,tw[19].r,tw[19].i,0,pd+6,pd+38);
  cplx_load_mulmuladdsub_pp(_z[26].r,_z[26].i,_z[27].r,_z[27].i,tw[11].r,tw[11].i,tw[27].r,tw[27].i,0,pd+22,pd+54);
  cplx_addsub(_z[24].r,_z[24].i,_z[26].r,_z[26].i);
  cplx_mul_1_4_F_addsub(_z[25].r,_z[25].i,_z[27].r,_z[27].i);
  cplx_load_mulmuladdsub_pp(_z[28].r,_z[28].i,_z[29].r,_z[29].i,tw[7].r,tw[7].i,tw[23].r,tw[23].i,0,pd+14,pd+46);
  cplx_load_mulmuladdsub_pp(_z[30].r,_z[30].i,_z[31].r,_z[31].i,tw[15].r,tw[15].i,tw[31].r,tw[31].i,0,pd+30,pd+62);
  cplx_addsub(_z[28].r,_z[28].i,_z[30].r,_z[30].i);
  cplx_mul_1_4_F_addsub(_z[29].r,_z[29].i,_z[31].r,_z[31].i);
  cplx_addsub(_z[24].r,_z[24].i,_z[28].r,_z[28].i);
  cplx_mul_1_8_F_addsub(_z[25].r,_z[25].i,_z[29].r,_z[29].i);
  cplx_mul_1_4_F_addsub(_z[26].r,_z[26].i,_z[30].r,_z[30].i);
  cplx_mul_3_8_F_addsub(_z[27].r,_z[27].i,_z[31].r,_z[31].i);
  cplx_addsub(_z[16].r,_z[16].i,_z[24].r,_z[24].i);
  cplx_muladdsub(_z[17].r,_z[17].i,_z[25].r,_z[25].i,F_1_16r,F_1_16i);
  cplx_mul_1_8_F_addsub(_z[18].r,_z[18].i,_z[26].r,_z[26].i);
  cplx_muladdsub(_z[19].r,_z[19].i,_z[27].r,_z[27].i,F_3_16r,F_3_16i);
  cplx_mul_1_4_F_addsub(_z[20].r,_z[20].i,_z[28].r,_z[28].i);
  cplx_muladdsub(_z[21].r,_z[21].i,_z[29].r,_z[29].i,F_5_16r,F_5_16i);
  cplx_mul_3_8_F_addsub(_z[22].r,_z[22].i,_z[30].r,_z[30].i);
  cplx_muladdsub(_z[23].r,_z[23].i,_z[31].r,_z[31].i,F_7_16r,F_7_16i);

  cplx_addsub(_z[0].r,_z[0].i,_z[16].r,_z[16].i);
  cplx_muladdsub(_z[1].r,_z[1].i,_z[17].r,_z[17].i,F_1_32r,F_1_32i);
  cplx_muladdsub(_z[2].r,_z[2].i,_z[18].r,_z[18].i,F_1_16r,F_1_16i);
  cplx_muladdsub(_z[3].r,_z[3].i,_z[19].r,_z[19].i,F_3_32r,F_3_32i);
  cplx_mul_1_8_F_addsub(_z[4].r,_z[4].i,_z[20].r,_z[20].i);
  cplx_muladdsub(_z[5].r,_z[5].i,_z[21].r,_z[21].i,F_5_32r,F_5_32i);
  cplx_muladdsub(_z[6].r,_z[6].i,_z[22].r,_z[22].i,F_3_16r,F_3_16i);
  cplx_muladdsub(_z[7].r,_z[7].i,_z[23].r,_z[23].i,F_7_32r,F_7_32i);
  cplx_mul_1_4_F_addsub(_z[8].r,_z[8].i,_z[24].r,_z[24].i);
  cplx_muladdsub(_z[9].r,_z[9].i,_z[25].r,_z[25].i,F_9_32r,F_9_32i);
  cplx_muladdsub(_z[10].r,_z[10].i,_z[26].r,_z[26].i,F_5_16r,F_5_16i);
  cplx_muladdsub(_z[11].r,_z[11].i,_z[27].r,_z[27].i,F_11_32r,F_11_32i);
  cplx_mul_3_8_F_addsub(_z[12].r,_z[12].i,_z[28].r,_z[28].i);
  cplx_muladdsub(_z[13].r,_z[13].i,_z[29].r,_z[29].i,F_13_32r,F_13_32i);
  cplx_muladdsub(_z[14].r,_z[14].i,_z[30].r,_z[30].i,F_7_16r,F_7_16i);
  cplx_muladdsub(_z[15].r,_z[15].i,_z[31].r,_z[31].i,F_15_32r,F_15_32i);
  return pd;
}



void fradix_32_notwd_first(struct y_complex * _z, y_ptr pd)
{
  ASSERT_ALIGNED_DOUBLE();
# ifndef NDEBUG1

  printf(" radix_32_notwd_first \n");
# endif

  cplx_addsub(_z[0].r,_z[0].i,_z[16].r,_z[16].i);
  cplx_addsub(_z[8].r,_z[8].i,_z[24].r,_z[24].i);
  cplx_addsub(_z[0].r,_z[0].i,_z[8].r,_z[8].i);
  cplx_mul_1_4_B_addsub(_z[16].r,_z[16].i,_z[24].r,_z[24].i);
  cplx_addsub(_z[4].r,_z[4].i,_z[20].r,_z[20].i);
  cplx_addsub(_z[12].r,_z[12].i,_z[28].r,_z[28].i);
  cplx_addsub(_z[4].r,_z[4].i,_z[12].r,_z[12].i);
  cplx_mul_1_4_B_addsub(_z[20].r,_z[20].i,_z[28].r,_z[28].i);
  cplx_addsub(_z[0].r,_z[0].i,_z[4].r,_z[4].i);
  cplx_mul_1_8_B_addsub(_z[16].r,_z[16].i,_z[20].r,_z[20].i);
  cplx_mul_1_4_B_addsub(_z[8].r,_z[8].i,_z[12].r,_z[12].i);
  cplx_mul_3_8_B_addsub(_z[24].r,_z[24].i,_z[28].r,_z[28].i);
  cplx_addsub(_z[2].r,_z[2].i,_z[18].r,_z[18].i);
  cplx_addsub(_z[10].r,_z[10].i,_z[26].r,_z[26].i);
  cplx_addsub(_z[2].r,_z[2].i,_z[10].r,_z[10].i);
  cplx_mul_1_4_B_addsub(_z[18].r,_z[18].i,_z[26].r,_z[26].i);
  cplx_addsub(_z[6].r,_z[6].i,_z[22].r,_z[22].i);
  cplx_addsub(_z[14].r,_z[14].i,_z[30].r,_z[30].i);
  cplx_addsub(_z[6].r,_z[6].i,_z[14].r,_z[14].i);
  cplx_mul_1_4_B_addsub(_z[22].r,_z[22].i,_z[30].r,_z[30].i);
  cplx_addsub(_z[2].r,_z[2].i,_z[6].r,_z[6].i);
  cplx_mul_1_8_B_addsub(_z[18].r,_z[18].i,_z[22].r,_z[22].i);
  cplx_mul_1_4_B_addsub(_z[10].r,_z[10].i,_z[14].r,_z[14].i);
  cplx_mul_3_8_B_addsub(_z[26].r,_z[26].i,_z[30].r,_z[30].i);
  cplx_addsub(_z[0].r,_z[0].i,_z[2].r,_z[2].i);
  cplx_muladdsub(_z[16].r,_z[16].i,_z[18].r,_z[18].i,B_1_16r,B_1_16i);
  cplx_mul_1_8_B_addsub(_z[8].r,_z[8].i,_z[10].r,_z[10].i);
  cplx_muladdsub(_z[24].r,_z[24].i,_z[26].r,_z[26].i,B_3_16r,B_3_16i);
  cplx_mul_1_4_B_addsub(_z[4].r,_z[4].i,_z[6].r,_z[6].i);
  cplx_muladdsub(_z[20].r,_z[20].i,_z[22].r,_z[22].i,B_5_16r,B_5_16i);
  cplx_mul_3_8_B_addsub(_z[12].r,_z[12].i,_z[14].r,_z[14].i);
  cplx_muladdsub(_z[28].r,_z[28].i,_z[30].r,_z[30].i,B_7_16r,B_7_16i);

  cplx_addsub(_z[1].r,_z[1].i,_z[17].r,_z[17].i);
  cplx_addsub(_z[9].r,_z[9].i,_z[25].r,_z[25].i);
  cplx_addsub(_z[1].r,_z[1].i,_z[9].r,_z[9].i);
  cplx_mul_1_4_B_addsub(_z[17].r,_z[17].i,_z[25].r,_z[25].i);
  cplx_addsub(_z[5].r,_z[5].i,_z[21].r,_z[21].i);
  cplx_addsub(_z[13].r,_z[13].i,_z[29].r,_z[29].i);
  cplx_addsub(_z[5].r,_z[5].i,_z[13].r,_z[13].i);
  cplx_mul_1_4_B_addsub(_z[21].r,_z[21].i,_z[29].r,_z[29].i);
  cplx_addsub(_z[1].r,_z[1].i,_z[5].r,_z[5].i);
  cplx_mul_1_8_B_addsub(_z[17].r,_z[17].i,_z[21].r,_z[21].i);
  cplx_mul_1_4_B_addsub(_z[9].r,_z[9].i,_z[13].r,_z[13].i);
  cplx_mul_3_8_B_addsub(_z[25].r,_z[25].i,_z[29].r,_z[29].i);
  cplx_addsub(_z[3].r,_z[3].i,_z[19].r,_z[19].i);
  cplx_addsub(_z[11].r,_z[11].i,_z[27].r,_z[27].i);
  cplx_addsub(_z[3].r,_z[3].i,_z[11].r,_z[11].i);
  cplx_mul_1_4_B_addsub(_z[19].r,_z[19].i,_z[27].r,_z[27].i);
  cplx_addsub(_z[7].r,_z[7].i,_z[23].r,_z[23].i);
  cplx_addsub(_z[15].r,_z[15].i,_z[31].r,_z[31].i);
  cplx_addsub(_z[7].r,_z[7].i,_z[15].r,_z[15].i);
  cplx_mul_1_4_B_addsub(_z[23].r,_z[23].i,_z[31].r,_z[31].i);
  cplx_addsub(_z[3].r,_z[3].i,_z[7].r,_z[7].i);
  cplx_mul_1_8_B_addsub(_z[19].r,_z[19].i,_z[23].r,_z[23].i);
  cplx_mul_1_4_B_addsub(_z[11].r,_z[11].i,_z[15].r,_z[15].i);
  cplx_mul_3_8_B_addsub(_z[27].r,_z[27].i,_z[31].r,_z[31].i);
  cplx_addsub(_z[1].r,_z[1].i,_z[3].r,_z[3].i);
  cplx_muladdsub(_z[17].r,_z[17].i,_z[19].r,_z[19].i,B_1_16r,B_1_16i);
  cplx_mul_1_8_B_addsub(_z[9].r,_z[9].i,_z[11].r,_z[11].i);
  cplx_muladdsub(_z[25].r,_z[25].i,_z[27].r,_z[27].i,B_3_16r,B_3_16i);
  cplx_mul_1_4_B_addsub(_z[5].r,_z[5].i,_z[7].r,_z[7].i);
  cplx_muladdsub(_z[21].r,_z[21].i,_z[23].r,_z[23].i,B_5_16r,B_5_16i);
  cplx_mul_3_8_B_addsub(_z[13].r,_z[13].i,_z[15].r,_z[15].i);
  cplx_muladdsub(_z[29].r,_z[29].i,_z[31].r,_z[31].i,B_7_16r,B_7_16i);

  cplx_addsub_store_p(0,pd,pd+32,_z[0].r,_z[0].i,_z[1].r,_z[1].i);
  cplx_muladdsub_store_p(0,pd+2,pd+34,_z[16].r,_z[16].i,_z[17].r,_z[17].i,B_1_32r,B_1_32i);
  cplx_muladdsub_store_p(0,pd+4,pd+36,_z[8].r,_z[8].i,_z[9].r,_z[9].i,B_1_16r,B_1_16i);
  cplx_muladdsub_store_p(0,pd+6,pd+38,_z[24].r,_z[24].i,_z[25].r,_z[25].i,B_3_32r,B_3_32i);
  cplx_mul_1_8_B_addsub_store_p(0,pd+8,pd+40,_z[4].r,_z[4].i,_z[5].r,_z[5].i);
  cplx_muladdsub_store_p(0,pd+10,pd+42,_z[20].r,_z[20].i,_z[21].r,_z[21].i,B_5_32r,B_5_32i);
  cplx_muladdsub_store_p(0,pd+12,pd+44,_z[12].r,_z[12].i,_z[13].r,_z[13].i,B_3_16r,B_3_16i);
  cplx_muladdsub_store_p(0,pd+14,pd+46,_z[28].r,_z[28].i,_z[29].r,_z[29].i,B_7_32r,B_7_32i);
  cplx_mul_1_4_B_addsub_store_p(0,pd+16,pd+48,_z[2].r,_z[2].i,_z[3].r,_z[3].i);
  cplx_muladdsub_store_p(0,pd+18,pd+50,_z[18].r,_z[18].i,_z[19].r,_z[19].i,B_9_32r,B_9_32i);
  cplx_muladdsub_store_p(0,pd+20,pd+52,_z[10].r,_z[10].i,_z[11].r,_z[11].i,B_5_16r,B_5_16i);
  cplx_muladdsub_store_p(0,pd+22,pd+54,_z[26].r,_z[26].i,_z[27].r,_z[27].i,B_11_32r,B_11_32i);
  cplx_mul_3_8_B_addsub_store_p(0,pd+24,pd+56,_z[6].r,_z[6].i,_z[7].r,_z[7].i);
  cplx_muladdsub_store_p(0,pd+26,pd+58,_z[22].r,_z[22].i,_z[23].r,_z[23].i,B_13_32r,B_13_32i);
  cplx_muladdsub_store_p(0,pd+28,pd+60,_z[14].r,_z[14].i,_z[15].r,_z[15].i,B_7_16r,B_7_16i);
  cplx_muladdsub_store_p(0,pd+30,pd+62,_z[30].r,_z[30].i,_z[31].r,_z[31].i,B_15_32r,B_15_32i);
}


void radix_32_dif_square_dit(y_ptr d, y_size_t nc)
{
  struct y_complex t[32],u[32],tw[32];
  y_ptr pt,pu;
  y_size_t ki,kj,j,nr,bi,bj,jj,ii,l,inc=Y_POWERS[31];

  /*HACK_ALIGN_STACK_ODD();*/

  /*DEBUG_VARS;*/
#ifndef NDEBUG1

  printf(" radix_32_dif_square_dit \n");
#endif

  ASSERT_ALIGNED_DOUBLE();

  /*if(Y_PLAN[Y_NRADICES-1] != 32) return;*/
  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 8, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 16 
  */

  /* Cycle 0 block 0 */
  j=0;
  fget_twiddle_factors_32_last( tw , j,  inc);

  pt=fradix_32_twd_last(t,d,j,tw);
  /* dyadic mul for k=0 */
  ysquare_nested_eq_c( t , 0 );
  for (ki=1,kj=31;ki<=15;ki++,kj--)
    {
      ysquare_nested_c( t, t, ki, kj);
    }
  ysquare_nested_eq_c( t , 16 );
  fradix_32_notwd_first(t,pt);

  if(nc==0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-2-(int)(nc));nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>5;
      bj=(Y_LRIGHT[nr]>>5)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          fget_twiddle_factors_32_last(tw, jj, inc);
          l=jj<<5;
          pu=fradix_32_twd_last(u,d,l,tw);
          fget_twiddle_factors_32_last(tw, ii, inc);
          l=ii<<5;
          pt=fradix_32_twd_last(t,d,l,tw);
          for(ki=0,kj=31;ki<=31;ki++,kj--)
            {
              ysquare_nested_c( t, u, ki, kj);
            }
          fradix_32_notwd_first(t,pt);
          fradix_32_notwd_first(u,pu);
        }
      if(ii==jj)
        {
          fget_twiddle_factors_32_last(tw ,ii ,inc);
          l=ii<<5;
          pt=fradix_32_twd_last(t,d,l,tw);
          for(ki=0,kj=31;ki<=15;ki++,kj--)
            {
              ysquare_nested_c( t, t, ki, kj);
            }
          fradix_32_notwd_first(t,pt);
        }
    }
}


void radix_32_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  struct y_complex t[32],u[32],s[32],v[32],tw[32];
  y_ptr pu,ps,pp;
  y_size_t ki,kj,j,bi,bj,jj,ii,l,inc=Y_POWERS[31];
  int nr;


  ASSERT_ALIGNED_DOUBLE();


  /*
  Now, for every radix reduction in the Y_PLAN we make:
    1) The last radix dif reduction
    2) The dyadic square.
    3) The first radix dit-notwd reduction  
    The dyadic mul is between two blocks of length 8, so we have
    to compute which blocks we need
    
    The first cycle, asociated with block 0, is special because of
    k=0 frecuency.
    
    The pad is 1 and bigpad 32 
  */

  /* Cycle 0 block 0 */
  j=0;
  fget_twiddle_factors_32_last(tw,j,inc);
  ps=fradix_32_twd_last(s,d1,j,tw);

  pp=fradix_32_twd_last(t,d2,j,tw);
  /* dyadic mul for k=0 */
  yconv_nested_c( s , s , t, t, 0, 0);
  for (ki=1,kj=31;ki<=16;ki++,kj--)
    {
      yconv_nested_c( s, s, t, t, ki, kj);
    }
  fradix_32_notwd_first(s,ps);
  if(nc == 0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>5;
      bj=(Y_LRIGHT[nr]>>5)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          fget_twiddle_factors_32_last(tw,jj,inc);
          l=jj<<5;
          pu=fradix_32_twd_last(u,d1,l,tw);
          pp=fradix_32_twd_last(v,d2,l,tw);
          fget_twiddle_factors_32_last(tw,ii,inc);
          l=ii<<5;
          ps=fradix_32_twd_last(s,d1,l,tw);
          pp=fradix_32_twd_last(t,d2,l,tw);
          for (ki=0,kj=31;ki<=31;ki++,kj--)
            {
              yconv_nested_c( s, u, t, v, ki, kj);
            }
          fradix_32_notwd_first(s,ps);
          fradix_32_notwd_first(u,pu);
        }
      if(ii==jj)
        {
          fget_twiddle_factors_32_last(tw,ii,inc);
          l=ii<<5;
          ps=fradix_32_twd_last(s,d1,l,tw);
          pp=fradix_32_twd_last(t,d2,l,tw);
          for (ki=0,kj=31;ki<=15;ki++,kj--)
            {
              yconv_nested_c( s, s, t, t, ki, kj);
            }
          fradix_32_notwd_first(s,ps);
        }
    }
}

void radix_32_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                                y_size_t  jj)
{
  struct y_complex t[32],u[32],s[32],v[32],tw[32];
  y_ptr pu,ps,pp;
  y_size_t i,j,ki,kj,nb,bi,bj,ni,l,inc=Y_POWERS[31];

  ASSERT_ALIGNED_DOUBLE();


  nb=bs>>5;
  bi=(bs*ii)>>5;
  bj=((bs*jj)>>5)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      fget_twiddle_factors_32_last(tw,j,inc);
      l=j<<5;
      pu=fradix_32_twd_last(u,d1,l,tw);
      pp=fradix_32_twd_last(v,d2,l,tw);
      fget_twiddle_factors_32_last(tw,i,inc);
      l=i<<5;
      ps=fradix_32_twd_last(s,d1,l,tw);
      pp=fradix_32_twd_last(t,d2,l,tw);
      for (ki=0,kj=31;ki<=31;ki++,kj--)
        {
          yconv_nested_c( s, u, t, v, ki, kj);
        }
      fradix_32_notwd_first(s,ps);
      fradix_32_notwd_first(u,pu);
    }
  if((i==j)&&(ni<nb))
    {
      fget_twiddle_factors_32_last(tw,i,inc);
      l=i<<5;
      ps=fradix_32_twd_last(s,d1,l,tw);
      pp=fradix_32_twd_last(t,d2,l,tw);
      for (ki=0,kj=31;ki<=15;ki++,kj--)
        {
          yconv_nested_c( s, s, t, t, ki, kj);
        }
      fradix_32_notwd_first(s,ps);
    }
}

void radix_32_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                   y_size_t jj)
{
  struct y_complex t[32],u[32],tw[32];
  y_ptr pu,pt;
  y_size_t i,j,ki,kj,nb,bi,bj,ni,l,inc=Y_POWERS[31];

  /*DEBUG_VARS;*/
#ifndef NDEBUG1

  printf(" radix_32_dif_square_dit_block \n");
#endif

  ASSERT_ALIGNED_DOUBLE();

  /*HACK_ALIGN_STACK_ODD();*/

  nb=bs>>5;
  bi=(bs*ii)>>5;
  bj=((bs*jj)>>5)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      fget_twiddle_factors_32_last(tw,j,inc);
      l=j<<5;
      pu=fradix_32_twd_last(u,d,l,tw);
      fget_twiddle_factors_32_last(tw,i,inc);
      l=i<<5;
      pt=fradix_32_twd_last(t,d,l,tw);
      for(ki=0,kj=31;ki<=31;ki++,kj--)
        {
          ysquare_nested_c( t, u, ki, kj);
        }
      fradix_32_notwd_first(t,pt);
      fradix_32_notwd_first(u,pu);
    }
  if((i==j)&&(ni<nb))
    {
      fget_twiddle_factors_32_last(tw,i,inc);
      l=i<<5;
      pt=fradix_32_twd_last(t,d,l,tw);
      for(ki=0,kj=31;ki<=15;ki++,kj--)
        {
          ysquare_nested_c( t, t, ki, kj);
        }
      fradix_32_notwd_first(t,pt);
    }
}

#else

void radix_32_dif_square_dit(y_ptr d, y_size_t nc)
{
  printf("Called radix_32_dif_square_dit , d=%p, nc=%d\n",d,nc);
  printf("Please compile again the file dyadic32.c with \n -DY_AVAL=5\n");
  exit(-1);
}

void radix_32_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  printf("Called radix_32_dif_mul_dit , d1=%p, d2=%p, nc=%d\n",d1,d2,nc);
  printf("Please compile again the file dyadic32.c with \n -DY_AVAL=5\n");
  exit(-1);
}

void radix_32_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                                y_size_t  jj)
{
  printf("Called radix_32_dif_mul_dit_block , d1=%p, d2=%p, bs=%d, ii=%d, jj=%d\n",d1,d2,bs,ii,jj);
  printf("Please compile again the file dyadic32.c with \n -DY_AVAL=5\n");
  exit(-1);
}

void radix_32_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                   y_size_t jj)
{
  printf("Called radix_32_dif_square_dit_block , d=%p, bs=%d, ii=%d, jj=%d\n",d,bs,ii,jj);
  printf("Please compile again the file dyadic32.c with \n -DY_AVAL=5\n");
  exit(-1);
}

#endif
/*$Id$*/








