/*$Id$*/
/*
   (c) 2000-2006 Guillermo Ballester Valor 
   
   It includes ideas from many Authors.
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
#define YNORM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#include "glucas.h"
#include "fft3.h"
#include "gsetup.h"
#include "ynorm.h"


#define get_pads_12 \
  pad2=(pad<<1);\
  pad3=pad+(pad<<1);\
  pad4=(pad<<2);\
  pad5=pad2+pad3;\
  pad6=pad3+pad3;\
  pad7=pad4+pad3;\
  pad8=(pad<<3);\
  pad9=pad8+pad;\
  pad10=pad8+pad2;\
  pad11=pad8+pad3;\
  pd0= x;\
  pd1= x + addr(pad<<1);\
  pd2= x + addr(pad2<<1);\
  pd3= x + addr(pad3<<1);\
  pd4= x + addr(pad4<<1);\
  pd5= x + addr(pad5<<1);\
  pd6= x + addr(pad6<<1);\
  pd7= x + addr(pad7<<1);\
  pd8= x + addr(pad8<<1);\
  pd9= x + addr(pad9<<1);\
  pd10= x + addr(pad10<<1);\
  pd11= x + addr(pad11<<1);



#define init_bjs_12 \
  bj0=N;\
  bj1=(b%12)<<(Y_K+1); \
  bj2=((2*b)%12)<<(Y_K+1); \
  bj3=((3*b)%12)<<(Y_K+1); \
  bj4=((4*b)%12)<<(Y_K+1); \
  bj5=((5*b)%12)<<(Y_K+1); \
  bj6=((6*b)%12)<<(Y_K+1); \
  bj7=((7*b)%12)<<(Y_K+1); \
  bj8=((8*b)%12)<<(Y_K+1); \
  bj9=((9*b)%12)<<(Y_K+1); \
  bj10=((10*b)%12)<<(Y_K+1); \
  bj11=((11*b)%12)<<(Y_K+1);


#define init_twiddle_factors_12\
  tw1r=1.0; tw1i=0.0;\
  tw2r=tw1r; tw2i=tw1i;\
  tw3r=tw1r; tw3i=tw1i;\
  tw4r=tw1r; tw4i=tw1i;\
  tw5r=tw1r; tw5i=tw1i;\
  tw6r=tw1r; tw6i=tw1i;\
  tw7r=tw1r; tw7i=tw1i;\
  tw8r=tw1r; tw8i=tw1i;\
  tw9r=tw1r; tw9i=tw1i;\
  tw10r=tw1r; tw10i=tw1i;\
  tw11r=tw1r; tw11i=tw1i;

#ifdef Y_KILL_BRANCHES

# define load_dwf_12(_i)                                                      \
	tt0[0] = two_to_phi[addr(_i)];                                        \
	tt0[1] = two_to_minusphi[addr(_i)];                                   \
	tt1[0] = two_to_phi[addr(((pad2 >> SHIFT_UPDATE) + _i))];             \
	tt1[1] = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];        \
	tt2[0] = two_to_phi[addr(((pad4 >> SHIFT_UPDATE) + _i))];             \
	tt2[1] = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];        \
	tt3[0] = two_to_phi[addr(((pad6 >> SHIFT_UPDATE) + _i))];             \
	tt3[1] = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];        \
	tt4[0] = two_to_phi[addr(((pad8 >> SHIFT_UPDATE) + _i))];             \
	tt4[1] = two_to_minusphi[addr(((pad8 >> SHIFT_UPDATE) + _i))];        \
	tt5[0] = two_to_phi[addr(((pad10 >> SHIFT_UPDATE) + _i))];            \
	tt5[1] = two_to_minusphi[addr(((pad10 >> SHIFT_UPDATE) + _i))];       \
	tt6[0] = two_to_phi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];       \
	tt6[1] = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];  \
	tt7[0] = two_to_phi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];       \
	tt7[1] = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];  \
	tt8[0] = two_to_phi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];       \
	tt8[1] = two_to_minusphi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];  \
	tt9[0] = two_to_phi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];       \
	tt9[1] = two_to_minusphi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];  \
	tt10[0] = two_to_phi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt10[1] = two_to_minusphi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt11[0] = two_to_phi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt11[1] = two_to_minusphi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];


# define  load_two_to_minusphi_12(_i) \
	tt0[1] = two_to_minusphi[addr(_i)]; \
	tt1[1] = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];\
	tt2[1] = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];\
	tt3[1] = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];\
	tt4[1] = two_to_minusphi[addr(((pad8 >> SHIFT_UPDATE) + _i))];\
	tt5[1] = two_to_minusphi[addr(((pad10 >> SHIFT_UPDATE) + _i))];\
	tt6[1] = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt7[1] = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt8[1] = two_to_minusphi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt9[1] = two_to_minusphi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt10[1] = two_to_minusphi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt11[1] = two_to_minusphi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];


# define  load_two_to_phi_12(_i) \
	tt0[0] = two_to_phi[addr(_i)]; \
	tt1[0] = two_to_phi[addr(((pad2 >> SHIFT_UPDATE) + _i))];\
	tt2[0] = two_to_phi[addr(((pad4 >> SHIFT_UPDATE) + _i))];\
	tt3[0] = two_to_phi[addr(((pad6 >> SHIFT_UPDATE) + _i))];\
	tt4[0] = two_to_phi[addr(((pad8 >> SHIFT_UPDATE) + _i))];\
	tt5[0] = two_to_phi[addr(((pad10 >> SHIFT_UPDATE) + _i))];\
	tt6[0] = two_to_phi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt7[0] = two_to_phi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt8[0] = two_to_phi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt9[0] = two_to_phi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt10[0] = two_to_phi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt11[0] = two_to_phi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];

#else

#define load_dwf_12(_i)                                                      \
	ttp0 = two_to_phi[addr(_i)];                                         \
	ttmp0 = two_to_minusphi[addr(_i)];                                   \
	ttp1 = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];                \
	ttmp1 = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];        \
	ttp2 = two_to_phi[addr(((pad4>>SHIFT_UPDATE) + _i))];                \
	ttmp2 = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];        \
	ttp3 = two_to_phi[addr(((pad6>>SHIFT_UPDATE) + _i))];                \
	ttmp3 = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];        \
	ttp4 = two_to_phi[addr(((pad8>>SHIFT_UPDATE) + _i))];                \
	ttmp4 = two_to_minusphi[addr(((pad8 >> SHIFT_UPDATE) + _i))];        \
	ttp5 = two_to_phi[addr(((pad10>>SHIFT_UPDATE) + _i))];               \
	ttmp5 = two_to_minusphi[addr(((pad10 >> SHIFT_UPDATE) + _i))];       \
	ttp6 = two_to_phi[addr(((pad6>>(SHIFT_UPDATE - 1)) + _i))];          \
	ttmp6 = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];  \
	ttp7 = two_to_phi[addr(((pad7>>(SHIFT_UPDATE - 1)) + _i))];          \
	ttmp7 = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];  \
	ttp8 = two_to_phi[addr(((pad8>>(SHIFT_UPDATE - 1)) + _i))];          \
	ttmp8 = two_to_minusphi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];  \
	ttp9 = two_to_phi[addr(((pad9>>(SHIFT_UPDATE - 1)) + _i))];          \
	ttmp9 = two_to_minusphi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];  \
	ttp10 = two_to_phi[addr(((pad10>>(SHIFT_UPDATE - 1)) + _i))];        \
	ttmp10 = two_to_minusphi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttp11 = two_to_phi[addr(((pad11>>(SHIFT_UPDATE - 1)) + _i))];        \
	ttmp11 = two_to_minusphi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];

# define  load_two_to_minusphi_12(_i) \
	ttmp0 = two_to_minusphi[addr(_i)]; \
	ttmp1 = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];\
	ttmp2 = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];\
	ttmp3 = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];\
	ttmp4 = two_to_minusphi[addr(((pad8 >> SHIFT_UPDATE) + _i))];\
	ttmp5 = two_to_minusphi[addr(((pad10 >> SHIFT_UPDATE) + _i))];\
	ttmp6 = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttmp7 = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttmp8 = two_to_minusphi[addr(((pad8 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttmp9 = two_to_minusphi[addr(((pad9 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttmp10 = two_to_minusphi[addr(((pad10 >> (SHIFT_UPDATE - 1)) + _i))];\
	ttmp11 = two_to_minusphi[addr(((pad11 >> (SHIFT_UPDATE - 1)) + _i))];


# define  load_two_to_phi_12(_i) \
	ttp0 = two_to_phi[addr(_i)]; \
	ttp1 = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	ttp2 = two_to_phi[addr(((pad4>>SHIFT_UPDATE) + _i))];\
	ttp3 = two_to_phi[addr(((pad6>>SHIFT_UPDATE) + _i))];\
	ttp4 = two_to_phi[addr(((pad8>>SHIFT_UPDATE) + _i))];\
	ttp5 = two_to_phi[addr(((pad10>>SHIFT_UPDATE) + _i))];\
	ttp6 = two_to_phi[addr(((pad6>>(SHIFT_UPDATE-1)) + _i))];\
	ttp7 = two_to_phi[addr(((pad7>>(SHIFT_UPDATE-1)) + _i))];\
	ttp8 = two_to_phi[addr(((pad8>>(SHIFT_UPDATE-1)) + _i))];\
	ttp9 = two_to_phi[addr(((pad9>>(SHIFT_UPDATE-1)) + _i))];\
	ttp10 = two_to_phi[addr(((pad10>>(SHIFT_UPDATE-1)) + _i))];\
	ttp11 = two_to_phi[addr(((pad11>>(SHIFT_UPDATE-1)) + _i))];

#endif

#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

/* get the basic twiddle from memory and computes its powers */
# include "yminimum.h"
# define get_twiddle_factors_12 \
          cplx_trig_min_load(tw1,px);\
          cplx_trig_min_square(tw2,tw1);\
          cplx_trig_min_square(tw4,tw2);\
          cplx_trig_min_square(tw8,tw4);\
          cplx_divmul(tw3,tw5,tw4,tw1);\
          cplx_divmul(tw7,tw9,tw8,tw1);\
          cplx_divmul(tw6,tw10,tw8,tw2);\
          cplx_trig_min_mul(tw11,tw8,tw3);\
          prefetch_p(px);

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 22

#  define get_twiddle_factors_12       \
        tw1r= px[0]; tw1i= px[1];      \
        tw2r= px[2]; tw2i= px[3];      \
        tw3r= px[4]; tw3i= px[5];      \
        tw4r= px[6]; tw4i= px[7];      \
        tw5r= px[8]; tw5i= px[9];      \
        tw6r= px[10]; tw6i= px[11];    \
        tw7r= px[12]; tw7i= px[13];    \
        tw8r= px[14]; tw8i= px[15];    \
        tw9r= px[16]; tw9i= px[17];    \
        tw10r= px[18]; tw10i= px[19];  \
        tw11r= px[20]; tw11i= px[21];  \
        prefetch_data_trig(px, Y_STEP); \
        prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);\
        prefetch_data_trig(px, Y_STEP + 2*Y_CACHE_LINE);


#else

# define get_twiddle_factors_12 \
          tw1r= px[0];\
          tw1i= px[1];\
	  tw2r= px[2];\
	  tw2i= px[3];\
	  tw4r= (tw2r + tw2i)*(tw2r - tw2i);\
	  tw4i= 2.0*tw2r*tw2i;\
	  tw3r= px[4];\
	  tw3i= px[5];\
	  tw8r= px[6];\
	  tw8i= px[7];\
	  cplx_divmul(tw5,tw11,tw8,tw3);\
	  cplx_divmul(tw6,tw10,tw8,tw2);\
	  cplx_divmul(tw7,tw9,tw8,tw1);\
          prefetch_data(px,8);\
          prefetch_data(px,8 + Y_CACHE_LINE);


# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 8
#endif

#define radix_12_twd_last_dit(_j) \
	  cplx_data_to_local_pp(t0,_j,pd0);\
	  cplx_load_mulmul_pp(t1,t2,tw4,tw8,_j,pd4,pd8);\
	  cplx_fft3(B,t0,t1,t2);\
	  cplx_load_mul_pp(t3,tw2,_j,pd2);\
	  cplx_load_mulmul_pp(t4,t5,tw6,tw10,_j,pd6,pd10);\
	  cplx_fft3(B,t3,t4,t5);\
	  cplx_addsub(t0,t3);\
	  cplx_muladdsub(t1,t4, B_1_6);\
	  cplx_muladdsub(t2,t5, B_1_3);\
	  cplx_load_mul_pp(t6,tw1,_j,pd1);\
	  cplx_load_mulmul_pp(t7,t8,tw5,tw9,_j,pd5,pd9);\
	  cplx_fft3(B,t6,t7,t8);\
	  cplx_load_mul_pp(t9,tw3,_j,pd3);\
	  cplx_load_mulmul_pp(t10,t11,tw7,tw11,_j,pd7,pd11);\
	  cplx_fft3(B,t9,t10,t11);\
	  cplx_addsub(t6,t9);\
	  cplx_muladdsub(t7,t10, B_1_6);\
	  cplx_muladdsub(t8,t11, B_1_3);\
	  \
	  cplx_addsub(t0,t6);\
	  cplx_muladdsub(t1,t7, B_1_12);\
	  cplx_muladdsub(t2,t8, B_1_6);\
	  cplx_mul_1_4_B_addsub(t3,t9);\
	  cplx_muladdsub(t4,t10, B_1_3);\
	  cplx_muladdsub(t5,t11, B_5_12);


#define radix_12_first_dif(_j) \
	  cplx_fft3(F,t0,t4,t8);\
	  cplx_fft3(F,t2,t6,t10);\
	  cplx_addsub(t0,t2);\
	  cplx_muladdsub(t4,t6, F_1_6);\
	  cplx_muladdsub(t8,t10, F_1_3);\
	  cplx_fft3(F,t1,t5,t9);\
	  cplx_fft3(F,t3,t7,t11);\
	  cplx_addsub(t1,t3);\
	  cplx_muladdsub(t5,t7, F_1_6);\
	  cplx_muladdsub(t9,t11, F_1_3);\
	  \
	  cplx_addsub_store_p(_j,pd0,pd6,t0,t1);\
	  cplx_muladdsub_store_p(_j,pd1,pd7,t4,t5, F_1_12);\
	  cplx_muladdsub_store_p(_j,pd2,pd8,t8,t9, F_1_6);\
	  cplx_mul_1_4_F_addsub_store_p(_j,pd3,pd9,t2,t3);\
	  cplx_muladdsub_store_p(_j,pd4,pd10,t6,t7, F_1_3);\
	  cplx_muladdsub_store_p(_j,pd5,pd11,t10,t11,F_5_12);


#define radix_12_first_dif_add(_d) \
	  cplx_fft3(F,t0,t4,t8);\
	  cplx_fft3(F,t2,t6,t10);\
	  cplx_addsub(t0,t2);\
	  cplx_muladdsub(t4,t6, F_1_6);\
	  cplx_muladdsub(t8,t10, F_1_3);\
	  cplx_fft3(F,t1,t5,t9);\
	  cplx_fft3(F,t3,t7,t11);\
	  cplx_addsub(t1,t3);\
	  cplx_muladdsub(t5,t7, F_1_6);\
	  cplx_muladdsub(t9,t11, F_1_3);\
	  \
	  cplx_addsub(t0,t1);\
	  cplx_local_to_data_add(_d, 0,t0);\
	  cplx_local_to_data_add(_d,pad6,t1);\
	  cplx_muladdsub(t4,t5, F_1_12);\
	  cplx_local_to_data_add(_d,pad,t4);\
	  cplx_local_to_data_add(_d,pad7,t5);\
	  cplx_muladdsub(t8,t9, F_1_6);\
	  cplx_local_to_data_add(_d,pad2,t8);\
	  cplx_local_to_data_add(_d,pad8,t9);\
	  cplx_mul_1_4_F_addsub(t2,t3);\
	  cplx_local_to_data_add(_d,pad3,t2);\
	  cplx_local_to_data_add(_d,pad9,t3);\
	  cplx_muladdsub(t6,t7, F_1_3);\
	  cplx_local_to_data_add(_d,pad4,t6);\
	  cplx_local_to_data_add(_d,pad10,t7);\
	  cplx_muladdsub(t10,t11,F_5_12);\
	  cplx_local_to_data_add(_d,pad5,t10);\
	  cplx_local_to_data_add(_d,pad11,t11);

#define radix_12_first_dif_add_shift(_d,_ofs) \
	  cplx_fft3(F,t0,t4,t8);              \
	  cplx_fft3(F,t2,t6,t10);             \
	  cplx_addsub(t0,t2);                 \
	  cplx_muladdsub(t4,t6, F_1_6);       \
	  cplx_muladdsub(t8,t10, F_1_3);      \
	  cplx_fft3(F,t1,t5,t9);              \
	  cplx_fft3(F,t3,t7,t11);             \
	  cplx_addsub(t1,t3);                 \
	  cplx_muladdsub(t5,t7, F_1_6);       \
	  cplx_muladdsub(t9,t11, F_1_3);      \
	                                      \
	  cplx_addsub(t0,t1);                 \
	  cplx_local_to_data_add(_d,_ofs,t0); \
	  cplx_local_to_data_add(_d,pad6,t1); \
	  cplx_muladdsub(t4,t5, F_1_12);      \
	  cplx_local_to_data_add(_d,pad,t4);  \
	  cplx_local_to_data_add(_d,pad7,t5); \
	  cplx_muladdsub(t8,t9, F_1_6);       \
	  cplx_local_to_data_add(_d,pad2,t8); \
	  cplx_local_to_data_add(_d,pad8,t9); \
	  cplx_mul_1_4_F_addsub(t2,t3);       \
	  cplx_local_to_data_add(_d,pad3,t2); \
	  cplx_local_to_data_add(_d,pad9,t3); \
	  cplx_muladdsub(t6,t7, F_1_3);       \
	  cplx_local_to_data_add(_d,pad4,t6); \
	  cplx_local_to_data_add(_d,pad10,t7);\
	  cplx_muladdsub(t10,t11,F_5_12);     \
	  cplx_local_to_data_add(_d,pad5,t10);\
	  cplx_local_to_data_add(_d,pad11,t11);


#ifdef SUM_CHECK
#define sum_12(_Sum) \
     _Sum += ((t0r + t0i) + (t1r + t1i)) + ((t2r + t2i) + (t3r + t3i)) + \
             ((t4r + t4i) + (t5r + t5i)) + ((t6r + t6i) + (t7r + t7i)) + \
             ((t8r + t8i) + (t9r + t9i)) + ((t10r + t10i) + (t11r + t11i));
#endif

#if (Y_AVAL > 3) && defined(Y_MANY_REGISTERS)

void substract_two_12(BIG_DOUBLE *x, UL N)
{
  double t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,
  t5r,t5i,t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,xtwo;
  UL pad=Y_LRIGHT[1],pad2,pad3,pad4,pad5,pad6,pad7,pad8,pad9,pad10,pad11,
         ofs,ni,who,ib,a;

  /* the bit to substract */
  ib=Y_SBIT +1;
  if(ib == Y_EXPONENT)
    ib=0;
  /* The element ni where bit 'ib' is */
  ni=floor((BIG_DOUBLE)ib/Y_XBITS);
  /* now ib is the bit within x[ni] */
  ib -= (UL)ceil(Y_XBITS*ni);
  /* The offset, in size of complex, of group */
  ofs = (ni % (N/Y_PLAN[0]))>>1;
  /* What element of group is x[ni] */
  who=(ni*Y_PLAN[0])/N;
  /* The shifted value of minustwo */
  xtwo = -(BIG_DOUBLE)(1<<ib);
  /*computes (N - ni*b % N). We suposse up to 256 million exponents */
  if(ni)
    {
#if ULONG_MAX > 0xFFFFFFFF
      /* 64 bits integer machine */
      a = N -(ni*b)%N;
#else
      /* 32 bits integer machine */
      modmul_32(a,b,ni,N);
      a = N - a;
#endif

    }
  else
    a=0;
  /* Multiply by two_to_phi() */
  xtwo *= exp(a * M_LN2 / N);
  /* Do some easy work while is computing long latency functions */
  t0r=0.0;
  t0i=0.0;
  t1r=0.0;
  t1i=0.0;
  t2r=0.0;
  t2i=0.0;
  t3r=0.0;
  t3i=0.0;
  t4r=0.0;
  t4i=0.0;
  t5r=0.0;
  t5i=0.0;
  t6r=0.0;
  t6i=0.0;
  t7r=0.0;
  t7i=0.0;
  t8r=0.0;
  t8i=0.0;
  t9r=0.0;
  t9i=0.0;
  t10r=0.0;
  t10i=0.0;
  t11r=0.0;
  t11i=0.0;
  switch (who)
    {
    case 0:
      if(ni & 1U)
        t0i = xtwo;
      else
        t0r = xtwo;
      break;
    case 1:
      if(ni & 1U)
        t1i = xtwo;
      else
        t1r = xtwo;
      break;
    case 2:
      if(ni & 1U)
        t2i = xtwo;
      else
        t2r = xtwo;
      break;
    case 3:
      if(ni & 1U)
        t3i = xtwo;
      else
        t3r = xtwo;
      break;
    case 4:
      if(ni & 1U)
        t4i = xtwo;
      else
        t4r = xtwo;
      break;
    case 5:
      if(ni & 1U)
        t5i = xtwo;
      else
        t5r = xtwo;
      break;
    case 6:
      if(ni & 1U)
        t6i = xtwo;
      else
        t6r = xtwo;
      break;
    case 7:
      if(ni & 1U)
        t7i = xtwo;
      else
        t7r = xtwo;
      break;
    case 8:
      if(ni & 1U)
        t8i = xtwo;
      else
        t8r = xtwo;
      break;
    case 9:
      if(ni & 1U)
        t9i = xtwo;
      else
        t9r = xtwo;
      break;
    case 10:
      if(ni & 1U)
        t10i = xtwo;
      else
        t10r = xtwo;
      break;
    case 11:
      if(ni & 1U)
        t11i = xtwo;
      else
        t11r = xtwo;
      break;
    }

  /* Sumcheck */
#if defined(SUM_CHECK)
  SumIn += xtwo;
#endif

  /* prepare pointers */
  pad2=(pad<<1);
  pad4=(pad<<2);
  pad8=(pad<<3);
  pad3=pad+pad2;
  pad5=pad+pad4;
  pad6=pad3+pad3;
  pad7=pad4+pad3;
  pad9=pad5+pad4;
  pad10=pad6+pad4;
  pad11=pad7+pad4;
  pad += ofs;
  pad2 += ofs;
  pad3 += ofs;
  pad4 += ofs;
  pad5 += ofs;
  pad6 += ofs;
  pad7 += ofs;
  pad8 += ofs;
  pad9 += ofs;
  pad10 += ofs;
  pad11 += ofs;

  /* make first DIF and add to array*/
  radix_12_first_dif_add_shift(x,ofs);

  /* New Y_SBIT and Y_SBJ*/
  Y_SBIT += Y_SBIT;
  if(Y_SBIT >= Y_EXPONENT)
    Y_SBIT -= Y_EXPONENT;
}


#if !defined(_OPENMP) && !defined(_SUNMP)  && !defined(_PTHREADS)

void dit_carry_norm_dif_12(BIG_DOUBLE *x,UL N, UL err_flag )
{
#  ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6],maxerr=0.0;
  double Y_ALIGNED(16) tt0[2],tt1[2],tt2[2],tt3[2],tt4[2],tt5[2],tt6[2],tt7[2],tt8[2],
  tt9[2],tt10[2],tt11[2];
#  else

  BIG_DOUBLE hiinv=highinv, loinv=lowinv;
  BIG_DOUBLE maxerr=0.0,ttmpSmall=Gsmall,ttmpBig=Gbig,
                                  ttpSmall=Hsmall;
  BIG_DOUBLE ttmp0,ttmp1,ttmp2,ttmp3,ttmp4,ttmp5,ttmp6,ttmp7,ttmp8,
  ttmp9,ttmp10,ttmp11;
  BIG_DOUBLE ttp0,ttp1,ttp2,ttp3,ttp4,ttp5,ttp6,ttp7,ttp8,ttp9,ttp10,ttp11;
#  endif
#  if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A=bigA,B=bigB;
#  endif

  BIG_DOUBLE carry0=0.0,carry1=0.0,carry2=0.0,carry3=0.0,
                                          carry4=0.0,carry5=0.0,carry6=0.0,carry7=0.0,carry8=0.0,carry9=0.0,
                                                                                  carry10=0.0,carry11=0.0;
  BIG_DOUBLE tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i;
  BIG_DOUBLE t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,
  t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i;
  BIG_DOUBLE *px;
  UL bj0,bj1,bj2,bj3,bj4,bj5,bj6,bj7,bj8,bj9,bj10,bj11;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  UL pad=Y_LRIGHT[1],pad2,pad3,pad4,pad5,pad6,pad7,pad8,pad9,pad10,pad11;
  UL i,j,k,l,ll;
#  ifdef Y_KILL_BRANCHES

  long int nolast;
  UL ib[2];
#   ifndef USE_ASM

  UL issmall;
#   ifdef Y_VECTORIZE

  UL issmall1;
#    endif
#   endif

  aux[0]=highinv;
  aux[1]=lowinv;
  aux[2]=Gbig;
  aux[3]=Hsmall;
  aux[4]=bigA;
  aux[5]=0.5;
  ib[0]=(UL)(-c);
  ib[1]=(UL)(b);
#  endif


  ASSERT_ALIGNED_DOUBLE();

  if(Y_SBIT==0)
    carry0=-2.0;
  Err=0.0;
  get_pads_12;
  init_bjs_12;
  init_twiddle_factors_12;
  px=Y_TWDB[Y_NRADICES-2];

  for (i=0,j=0; j<pad2; j+=UPDATE,i++)
    {
      load_dwf_12 (i);
      for (k=0; k<UPDATE; k+=2, px+=Y_STEP)
        {
          l=(j+k)>>1;
#  ifdef Y_KILL_BRANCHES

          nolast = l - (pad-1);
          nolast >>=(BITS_PER_UL - 1);
#  endif

          ll= addr(k+j);
          radix_12_twd_last_dit( ll );
#  ifdef SUM_CHECK

          sum_12(SumOut);
#  endif

          if(err_flag)
            {
#  if defined(Y_VECTORIZE)
              cplx_carry_norm_check_vector(t,0,1);
              cplx_carry_norm_check_vector(t,2,3);
              cplx_carry_norm_check_vector(t,4,5);
              cplx_carry_norm_check_vector(t,6,7);
              cplx_carry_norm_check_vector(t,8,9);
#  else

              cplx_carry_norm_check(t,0);
              cplx_carry_norm_check(t,1);
              cplx_carry_norm_check(t,2);
              cplx_carry_norm_check(t,3);
              cplx_carry_norm_check(t,4);
              cplx_carry_norm_check(t,5);
              cplx_carry_norm_check(t,6);
              cplx_carry_norm_check(t,7);
              cplx_carry_norm_check(t,8);
              cplx_carry_norm_check(t,9);
#  endif

              cplx_carry_norm_check(t,10);
              cplx_carry_norm_last_check(t,11);
            }
          else
            {
#  if defined(Y_VECTORIZE)
              cplx_carry_norm_vector(t,0,1);
              cplx_carry_norm_vector(t,2,3);
              cplx_carry_norm_vector(t,4,5);
              cplx_carry_norm_vector(t,6,7);
              cplx_carry_norm_vector(t,8,9);
#  else

              cplx_carry_norm(t,0);
              cplx_carry_norm(t,1);
              cplx_carry_norm(t,2);
              cplx_carry_norm(t,3);
              cplx_carry_norm(t,4);
              cplx_carry_norm(t,5);
              cplx_carry_norm(t,6);
              cplx_carry_norm(t,7);
              cplx_carry_norm(t,8);
              cplx_carry_norm(t,9);
#  endif

              cplx_carry_norm(t,10);
              cplx_carry_norm_last(t,11);
            }
          radix_12_first_dif( ll );

# ifdef SUM_CHECK

          SumIn += (t0r + t0i);
# endif

          get_twiddle_factors_12;
        }
    }
  /* adjust the last carries */
  t0r=carry11;
  t1r=carry0;
  t2r=carry1;
  t3r=carry2;
  t4r=carry3;
  t5r=carry4;
  t6r=carry5;
  t7r=carry6;
  t8r=carry7;
  t9r=carry8;
  t10r=carry9;
  t11r=carry10;

  init_bjs_12;

  j=0;
  load_two_to_phi_12(j);
  cplx_carry(t,0);
  cplx_carry(t,1);
  cplx_carry(t,2);
  cplx_carry(t,3);
  cplx_carry(t,4);
  cplx_carry(t,5);
  cplx_carry(t,6);
  cplx_carry(t,7);
  cplx_carry(t,8);
  cplx_carry(t,9);
  cplx_carry(t,10);
  cplx_carry(t,11);

  /* add carry signal */
  radix_12_first_dif_add( x );

# ifdef SUM_CHECK

  SumIn += (t0r + t0i);
# endif

  if(err_flag)
    Err = maxerr;
  if(Y_SBIT)
    substract_two_12( x, N);
}

# else

/***********************************************************************/
/* These are the Multithread adds needed                               */
/***********************************************************************/
/*
_ofs is the offset, in complex size,  of every thread nthr
_ofs = _n0 * nthr;  0 <= nthr < max_num_of_threads
where _n0 is the basic length of a thread, and MUST be a
power of two
*/
#  define get_pads_12_mp(_ofs)     \
  pad2=(pad<<1);                   \
  pad4=(pad<<2);                   \
  pad3=pad+pad2;                   \
  pad5=pad+pad4;                   \
  pad6=pad3+pad3;                  \
  pad7=pad4+pad3;                  \
  pad8=(pad<<3);                   \
  pad9=pad8+pad;                   \
  pad10=pad8+pad2;                 \
  pad11=pad8+pad3;                 \
  pd0= x + addr((_ofs<<1));        \
  pd1= x + addr((pad  + _ofs)<<1); \
  pd2= x + addr((pad2 + _ofs)<<1); \
  pd3= x + addr((pad3 + _ofs)<<1); \
  pd4= x + addr((pad4 + _ofs)<<1); \
  pd5= x + addr((pad5 + _ofs)<<1); \
  pd6= x + addr((pad6 + _ofs)<<1); \
  pd7= x + addr((pad7 + _ofs)<<1); \
  pd8= x + addr((pad8 + _ofs)<<1); \
  pd9= x + addr((pad9 + _ofs)<<1); \
  pd10= x + addr((pad10 + _ofs)<<1);\
  pd11= x + addr((pad11 + _ofs)<<1);

/*
To compute the initial bj of every thread we use:

b(i,nthr)= (((nthr + i*(pad/n0))*b)% (N/(2*n0)) * 2 * n0;
where
pad/n0 must be integer and 0< (pad/n0) < max_num_of_threads
n0 must be a power of two.
It will be computed at first time. so we avoid to compute it every
iteration
*/

#  define init_bjs_12_mp \
  bj0= bp[0];\
  bj1= bp[1];\
  bj2= bp[2];\
  bj3= bp[3];\
  bj4= bp[4];\
  bj5= bp[5];\
  bj6= bp[6];\
  bj7= bp[7];\
  bj8= bp[8];\
  bj9= bp[9];\
  bj10= bp[10];\
  bj11= bp[11];


/*
The twidle factors now are not (1.0,0.0) at the begining
so it should be computed

inc=(Y_POWERS[Y_PLAN[0]-1])<<1;
px = Y_TWDB[Y_NRADICES - 2] + _ofs * inc - 2;

It could be computed at a first iteration, so we save
some work
*/

# if defined(Y_MINIMUM)
#  define init_twiddle_factors_12_mp(_px)  \
          tw1r=_px[0]; tw1i=_px[1];        \
          cplx_trig_min_square(tw2,tw1);   \
          cplx_trig_min_square(tw4,tw2);   \
          cplx_trig_min_square(tw8,tw4);   \
          cplx_divmul(tw3,tw5,tw4,tw1);    \
          cplx_divmul(tw7,tw9,tw8,tw1);    \
          cplx_divmul(tw6,tw10,tw8,tw2);   \
          cplx_trig_min_mul(tw11,tw8,tw3); \
          prefetch_p(_px);

#elif defined(Y_MAXIMUM)
#  define init_twiddle_factors_12_mp(_px)  \
         tw1r= _px[0];	 tw1i= _px[1];    \
         tw2r= _px[2];	 tw2i= _px[3];    \
	 tw3r= _px[4];   tw3i= _px[5];    \
         tw4r= _px[6];	 tw4i= _px[7];    \
         tw5r= _px[8];	 tw5i= _px[9];    \
         tw6r= _px[10];	 tw6i= _px[11];   \
         tw7r= _px[12];	 tw7i= _px[13];   \
         tw8r= _px[14];	 tw8i= _px[15];   \
         tw9r= _px[16];	 tw9i= _px[17];   \
         tw10r= _px[18]; tw10i= _px[19];   \
         tw11r= _px[20]; tw11i= _px[21];   \
         prefetch_data_trig(_px, Y_STEP); \
         prefetch_data_trig(_px, Y_STEP + Y_CACHE_LINE); \
         prefetch_data_trig(_px, Y_STEP + 2*Y_CACHE_LINE);\


# else


#  define init_twiddle_factors_12_mp(_px)     \
          tw1r= _px[0];  tw1i= _px[1];        \
	  tw2r= _px[2];  tw2i= _px[3];        \
	  tw4r= (tw2r + tw2i)*(tw2r - tw2i);  \
	  tw4i= 2.0*tw2r*tw2i;                \
	  tw3r= _px[4];  tw3i= _px[5];        \
	  tw8r= _px[6];  tw8i= _px[7];        \
	  cplx_divmul(tw5,tw11,tw8,tw3);      \
	  cplx_divmul(tw6,tw10,tw8,tw2);      \
	  cplx_divmul(tw7,tw9,tw8,tw1);       \
          prefetch_data(_px,8);               \
          prefetch_data(_px,8 + Y_CACHE_LINE);

# endif

# define save_carries_12_mp(_pcarr)  \
        _pcarr[0]=carry0;            \
        _pcarr[1]=carry1;            \
        _pcarr[2]=carry2;            \
        _pcarr[3]=carry3;            \
        _pcarr[4]=carry4;            \
        _pcarr[5]=carry5;            \
        _pcarr[6]=carry6;            \
        _pcarr[7]=carry7;            \
        _pcarr[8]=carry8;            \
        _pcarr[9]=carry9;            \
        _pcarr[10]=carry10;          \
        _pcarr[11]=carry11;


# define radix_12_first_dif_add_mp \
	  cplx_fft3(F,t0,t4,t8);\
	  cplx_fft3(F,t2,t6,t10);\
	  cplx_addsub(t0,t2);\
	  cplx_muladdsub(t4,t6, F_1_6);\
	  cplx_muladdsub(t8,t10, F_1_3);\
	  cplx_fft3(F,t1,t5,t9);\
	  cplx_fft3(F,t3,t7,t11);\
	  cplx_addsub(t1,t3);\
	  cplx_muladdsub(t5,t7, F_1_6);\
	  cplx_muladdsub(t9,t11, F_1_3);\
	  \
	  cplx_addsub(t0,t1);\
	  cplx_local_to_data_add(pd0, 0,t0);\
	  cplx_local_to_data_add(pd6, 0,t1);\
	  cplx_muladdsub(t4,t5, F_1_12);\
	  cplx_local_to_data_add(pd1, 0,t4);\
	  cplx_local_to_data_add(pd7, 0,t5);\
	  cplx_muladdsub(t8,t9, F_1_6);\
	  cplx_local_to_data_add(pd2, 0,t8);\
	  cplx_local_to_data_add(pd8, 0,t9);\
	  cplx_mul_1_4_F_addsub(t2,t3);\
	  cplx_local_to_data_add(pd3, 0,t2);\
	  cplx_local_to_data_add(pd9, 0,t3);\
	  cplx_muladdsub(t6,t7, F_1_3);\
	  cplx_local_to_data_add(pd4, 0,t6);\
	  cplx_local_to_data_add(pd10, 0,t7);\
	  cplx_muladdsub(t10,t11,F_5_12);\
	  cplx_local_to_data_add(pd5, 0,t10);\
	  cplx_local_to_data_add(pd11, 0,t11);

BIG_DOUBLE dit_carry_norm_dif_12_mp(BIG_DOUBLE *x,BIG_DOUBLE *ptx, BIG_DOUBLE *pcarr, UL *bp, UL N ,UL n0, UL nthr,UL err_flag)
{
#  ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6],maxerr=0.0;
  double Y_ALIGNED(16) tt0[2],tt1[2],tt2[2],tt3[2],tt4[2],tt5[2],tt6[2],tt7[2],tt8[2],
  tt9[2],tt10[2],tt11[2];
#  else

BIG_DOUBLE hiinv=highinv, loinv=lowinv;
BIG_DOUBLE maxerr=0.0,ttmpSmall=Gsmall,ttmpBig=Gbig,
                                ttpSmall=Hsmall;
BIG_DOUBLE ttmp0,ttmp1,ttmp2,ttmp3,ttmp4,ttmp5,ttmp6,ttmp7,ttmp8,
ttmp9,ttmp10,ttmp11;
BIG_DOUBLE ttp0,ttp1,ttp2,ttp3,ttp4,ttp5,ttp6,ttp7,ttp8,ttp9,ttp10,ttp11;
#  endif
#  if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A=bigA,B=bigB;
#  endif

  BIG_DOUBLE carry0=0.0,carry1=0.0,carry2=0.0,carry3=0.0,
                                          carry4=0.0,carry5=0.0,carry6=0.0,carry7=0.0,carry8=0.0,carry9=0.0,
                                                                                  carry10=0.0,carry11=0.0;
  BIG_DOUBLE tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i;
  BIG_DOUBLE t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,
  t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i;
  BIG_DOUBLE *px;
  UL bj0,bj1,bj2,bj3,bj4,bj5,bj6,bj7,bj8,bj9,bj10,bj11;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  UL pad=Y_LRIGHT[1],pad2,pad3,pad4,pad5,pad6,pad7,pad8,pad9,pad10,pad11;
  UL i, j, k, l, ll, ofs, lim;
#  ifdef Y_KILL_BRANCHES

  long int nolast;
  UL ib[2];
#   ifndef USE_ASM

  UL issmall;
#   ifdef Y_VECTORIZE

  UL issmall1;
#    endif
#   endif

  aux[0]=highinv;
  aux[1]=lowinv;
  aux[2]=Gbig;
  aux[3]=Hsmall;
  aux[4]=bigA;
  aux[5]=0.5;
  ib[0]=(UL)(-c);
  ib[1]=(UL)(b);
#  endif


  ASSERT_ALIGNED_DOUBLE();

  ofs = n0 * nthr;
  if(nthr < (Y_NTHREADS - 1))
    lim = ((ofs + n0) << 1);
  else
    lim = (Y_LRIGHT[1] << 1);

# if defined(SUM_CHECK)

  Y_SUMOUT[nthr] = 0.0;
  Y_SUMIN[nthr] = 0.0;
# endif

  if(nthr)
    {
      init_bjs_12_mp;
      get_pads_12;
      init_twiddle_factors_12_mp(ptx);
      px = ptx + Y_STEP;
    }
  else
    {
      if(Y_SBIT==0)
        carry0=-2.0;
      get_pads_12;
      init_bjs_12;
      init_twiddle_factors_12;
      px=Y_TWDB[Y_NRADICES-2];
    }

  for (i = (ofs << 1) / UPDATE, j = (ofs << 1); j < lim; j += UPDATE, i++)
    {
      load_dwf_12 (i);
      for (k=0; k<UPDATE; k+=2, px+=Y_STEP)
        {
          l=(j+k)>>1;

# ifdef Y_KILL_BRANCHES

          nolast = l - (pad-1);
          nolast >>=(BITS_PER_UL - 1);
# endif

          ll=addr(j+k);
          radix_12_twd_last_dit( ll );
# ifdef SUM_CHECK

          sum_12 (Y_SUMOUT[nthr]);
# endif

          if(err_flag)
            {
#  if defined(Y_VECTORIZE)
              cplx_carry_norm_check_vector(t,0,1);
              cplx_carry_norm_check_vector(t,2,3);
              cplx_carry_norm_check_vector(t,4,5);
              cplx_carry_norm_check_vector(t,6,7);
              cplx_carry_norm_check_vector(t,8,9);
#  else

cplx_carry_norm_check(t,0);
cplx_carry_norm_check(t,1);
cplx_carry_norm_check(t,2);
cplx_carry_norm_check(t,3);
cplx_carry_norm_check(t,4);
cplx_carry_norm_check(t,5);
cplx_carry_norm_check(t,6);
cplx_carry_norm_check(t,7);
cplx_carry_norm_check(t,8);
cplx_carry_norm_check(t,9);
#  endif

              cplx_carry_norm_check(t,10);
              cplx_carry_norm_last_check(t,11);
            }
          else
            {
#  if defined(Y_VECTORIZE)
              cplx_carry_norm_vector(t,0,1);
              cplx_carry_norm_vector(t,2,3);
              cplx_carry_norm_vector(t,4,5);
              cplx_carry_norm_vector(t,6,7);
              cplx_carry_norm_vector(t,8,9);
#  else

cplx_carry_norm(t,0);
cplx_carry_norm(t,1);
cplx_carry_norm(t,2);
cplx_carry_norm(t,3);
cplx_carry_norm(t,4);
cplx_carry_norm(t,5);
cplx_carry_norm(t,6);
cplx_carry_norm(t,7);
cplx_carry_norm(t,8);
cplx_carry_norm(t,9);
#  endif

              cplx_carry_norm(t,10);
              cplx_carry_norm_last(t,11);
            }
          radix_12_first_dif( ll );

# if defined(SUM_CHECK)

          Y_SUMIN[nthr] += (t0r + t0i);
# endif

          get_twiddle_factors_12;
        }
    }
  /* save the carries */
  save_carries_12_mp(pcarr);

  /* return the error */
  return maxerr;
}

void dit_carry_norm_dif_12_last_carries_mp(BIG_DOUBLE *x, UL *bp, UL N ,UL n0, int nthr)
{
#  ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6];
  double Y_ALIGNED(16) tt0[2],tt1[2],tt2[2],tt3[2],tt4[2],tt5[2],tt6[2],tt7[2],tt8[2],
  tt9[2],tt10[2],tt11[2];
#  else

BIG_DOUBLE hiinv=highinv, loinv=lowinv;
BIG_DOUBLE ttmpSmall=Gsmall,ttmpBig=Gbig,
                                    ttpSmall=Hsmall;
BIG_DOUBLE ttmp0,ttmp1,ttmp2,ttmp3,ttmp4,ttmp5,ttmp6,ttmp7,ttmp8,ttmp9,
ttmp10,ttmp11;
BIG_DOUBLE ttp0,ttp1,ttp2,ttp3,ttp4,ttp5,ttp6,ttp7,ttp8,ttp9,ttp10,ttp11;
#  endif
#  if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A=bigA,B=bigB;
#  endif

  BIG_DOUBLE carry0,carry1,carry2,carry3,carry4,carry5,carry6,carry7,carry8,
  carry9,carry10,carry11;
  BIG_DOUBLE t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,
  t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i;
  UL bj0,bj1,bj2,bj3,bj4,bj5,bj6,bj7,bj8,bj9,bj10,bj11;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11;
  UL pad=Y_LRIGHT[1],pad2,pad3,pad4,pad5,pad6,pad7,pad8,pad9,pad10,pad11;
  UL j,ofs;
#  ifdef Y_KILL_BRANCHES

  UL ib[2];
#   ifndef USE_ASM

  UL issmall;
#   endif

  aux[0]=highinv;
  aux[1]=lowinv;
  aux[2]=Gbig;
  aux[3]=Hsmall;
  aux[4]=bigA;
  aux[5]=0.5;
  ib[0]=(UL)(-c);
  ib[1]=(UL)(b);
#  endif

  ASSERT_ALIGNED_DOUBLE();

  /* Init the data */
  if(nthr)
    {
      t0r = Y_CARRIES[nthr-1][0];
      t1r = Y_CARRIES[nthr-1][1];
      t2r = Y_CARRIES[nthr-1][2];
      t3r = Y_CARRIES[nthr-1][3];
      t4r = Y_CARRIES[nthr-1][4];
      t5r = Y_CARRIES[nthr-1][5];
      t6r = Y_CARRIES[nthr-1][6];
      t7r = Y_CARRIES[nthr-1][7];
      t8r = Y_CARRIES[nthr-1][8];
      t9r = Y_CARRIES[nthr-1][9];
      t10r = Y_CARRIES[nthr-1][10];
      t11r = Y_CARRIES[nthr-1][11];
    }
  else
    {
      t0r = Y_CARRIES[Y_NTHREADS-1][11];
      t1r = Y_CARRIES[Y_NTHREADS-1][0];
      t2r = Y_CARRIES[Y_NTHREADS-1][1];
      t3r = Y_CARRIES[Y_NTHREADS-1][2];
      t4r = Y_CARRIES[Y_NTHREADS-1][3];
      t5r = Y_CARRIES[Y_NTHREADS-1][4];
      t6r = Y_CARRIES[Y_NTHREADS-1][5];
      t7r = Y_CARRIES[Y_NTHREADS-1][6];
      t8r = Y_CARRIES[Y_NTHREADS-1][7];
      t9r = Y_CARRIES[Y_NTHREADS-1][8];
      t10r = Y_CARRIES[Y_NTHREADS-1][9];
      t11r = Y_CARRIES[Y_NTHREADS-1][10];
    }

  init_bjs_12_mp;
  ofs=(n0*nthr);

  get_pads_12;
  j=(ofs<<1)/UPDATE;

  load_two_to_phi_12(j);
  load_two_to_minusphi_12(j);

  cplx_carry(t,0);
  cplx_carry(t,1);
  cplx_carry(t,2);
  cplx_carry(t,3);
  cplx_carry(t,4);
  cplx_carry(t,5);
  cplx_carry(t,6);
  cplx_carry(t,7);
  cplx_carry(t,8);
  cplx_carry(t,9);
  cplx_carry(t,10);
  cplx_carry(t,11);

  /* add carry signal */

  get_pads_12_mp(ofs)

  radix_12_first_dif_add_mp;

# if defined(SUM_CHECK)

  Y_SUMIN[nthr] += (t0r + t0i);
# endif

}

void dit_carry_norm_dif_12(BIG_DOUBLE *x, UL N ,UL err_flag)
{
  int i,n0;
  n0=Y_LRIGHT[1]/Y_NTHREADS;

  /* make n0 divisble by Y_UPDATE */
  n0 = ((n0 >> (Y_SHIFT_UPDATE - 1)) << (Y_SHIFT_UPDATE - 1));

  /* The multithread loop */
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(i)
#  pragma MP taskloop shared(x,Y_ERRS,Y_TWX,Y_CARRIES,Y_BJS,N,n0,err_flag)
# endif

  for(i=0;i<Y_NTHREADS;i++)
    {
      Y_ERRS[i]=dit_carry_norm_dif_12_mp( x, Y_TWX[i], Y_CARRIES[i], Y_BJS[i],
                                          N, n0, i, err_flag);
    }

  /* maxerr */
  Err=0;
  if(err_flag)
    {
      for(i=0;i<Y_NTHREADS;i++)
        if(Y_ERRS[i] > Err)
          Err = Y_ERRS[i];
    }

  /* And finally, the last carry propagation */
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(x,Y_BJS,N,n0)
# endif
  for(i=0;i<Y_NTHREADS;i++)
    {
      dit_carry_norm_dif_12_last_carries_mp( x, Y_BJS[i],
                                             N, n0, i);
    }
  /* Substract two if shifts*/
  if(Y_SBIT)
    substract_two_12( x, N);

#if defined(SUM_CHECK)

  for (i = 0; i < Y_NTHREADS; i++)
    {
      SumOut += Y_SUMOUT[i];
      SumIn += Y_SUMIN[i];
    }
#endif
}


# endif

#else

void dit_carry_norm_dif_12(BIG_DOUBLE *x,UL N, UL err_flag )
{
  printf("Called dit_carry_norm_dif_12 , x=%p, N=%ld, Err_flag=%ld\n",x,N,err_flag);
  printf("Please compile again the file ynorm_12.c with \n -DY_AVAL=4 -DY_MANY_REGISTERS flags\n");
  exit(-1);
}

#endif
/*$Id$*/







