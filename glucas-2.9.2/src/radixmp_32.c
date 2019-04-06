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
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1


/**********************************************************************/
/*    common blocks code and macros                                   */
/**********************************************************************/

/* asign the correct values to every pad */

#define get_pads_32       \
  jj = pad;               \
  pd0= d;                 \
  pd1= d + addr((jj<<1)); \
  jj+=pad;                \
  pd2= d + addr((jj<<1)); \
  jj+=pad;                \
  pd3= d + addr((jj<<1)); \
  jj+=pad;                \
  pd4= d + addr((jj<<1)); \
  jj+=pad;                \
  pd5= d + addr((jj<<1)); \
  jj+=pad;                \
  pd6= d + addr((jj<<1)); \
  jj+=pad;                \
  pd7= d + addr((jj<<1)); \
  jj+=pad;                \
  pd8= d + addr((jj<<1)); \
  jj+=pad;                \
  pd9= d + addr((jj<<1)); \
  jj+=pad;                \
  pd10= d + addr((jj<<1));\
  jj+=pad;                \
  pd11= d + addr((jj<<1));\
  jj+=pad;                \
  pd12= d + addr((jj<<1));\
  jj+=pad;                \
  pd13= d + addr((jj<<1));\
  jj+=pad;                \
  pd14= d + addr((jj<<1));\
  jj+=pad;                \
  pd15= d + addr((jj<<1));\
  jj+=pad;                \
  pd16= d + addr((jj<<1));\
  jj+=pad;                \
  pd17= d + addr((jj<<1));\
  jj+=pad;                \
  pd18= d + addr((jj<<1));\
  jj+=pad;                \
  pd19= d + addr((jj<<1));\
  jj+=pad;                \
  pd20= d + addr((jj<<1));\
  jj+=pad;                \
  pd21= d + addr((jj<<1));\
  jj+=pad;                \
  pd22= d + addr((jj<<1));\
  jj+=pad;                \
  pd23= d + addr((jj<<1));\
  jj+=pad;                \
  pd24= d + addr((jj<<1));\
  jj+=pad;                \
  pd25= d + addr((jj<<1));\
  jj+=pad;                \
  pd26= d + addr((jj<<1));\
  jj+=pad;                \
  pd27= d + addr((jj<<1));\
  jj+=pad;                \
  pd28= d + addr((jj<<1));\
  jj+=pad;                \
  pd29= d + addr((jj<<1));\
  jj+=pad;                \
  pd30= d + addr((jj<<1));\
  jj+=pad;                \
  pd31= d + addr((jj<<1));

/* get the basic twiddle from memory and computes its powers */
#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# include "yminimum.h"
# define get_twiddle_factors_32            \
          cplx_trig_min_load(tw1,px);      \
          cplx_trig_min_square(tw2,tw1);   \
          cplx_trig_min_square(tw4,tw2);   \
          cplx_trig_min_square(tw8,tw4);   \
          cplx_divmul(tw3,tw5,tw4,tw1);    \
          cplx_divmul(tw7,tw9,tw8,tw1);    \
          cplx_divmul(tw6,tw10,tw8,tw2);   \
          cplx_trig_min_square(tw16,tw8);  \
          cplx_divmul(tw15,tw17,tw16,tw1); \
          cplx_divmul(tw14,tw18,tw16,tw2); \
          cplx_divmul(tw13,tw19,tw16,tw3); \
          cplx_divmul(tw12,tw20,tw16,tw4); \
          cplx_divmul(tw11,tw21,tw16,tw5); \
          cplx_trig_min_square(tw26,tw13); \
          cplx_divmul(tw25,tw27,tw26,tw1); \
          cplx_divmul(tw24,tw28,tw26,tw2); \
          cplx_divmul(tw23,tw29,tw26,tw3); \
          cplx_divmul(tw22,tw30,tw26,tw4); \
          cplx_trig_min_mul(tw31,tw16,tw15);

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 62
# define get_twiddle_factors_32            \
     tw1r = px[ 0]; tw1i = px[ 1];         \
     tw2r = px[ 2]; tw2i = px[ 3];         \
     tw3r = px[ 4]; tw3i = px[ 5];         \
     tw4r = px[ 6]; tw4i = px[ 7];         \
     tw5r = px[ 8]; tw5i = px[ 9];         \
     tw6r = px[10]; tw6i = px[11];         \
     tw7r = px[12]; tw7i = px[13];         \
     tw8r = px[14]; tw8i = px[15];         \
     tw9r = px[16]; tw9i = px[17];         \
     tw10r= px[18]; tw10i= px[19];         \
     tw11r= px[20]; tw11i= px[21];         \
     tw12r= px[22]; tw12i= px[23];         \
     tw13r= px[24]; tw13i= px[25];         \
     tw14r= px[26]; tw14i= px[27];         \
     tw15r= px[28]; tw15i= px[29];         \
     prefetch_data_trig(px, Y_STEP + Y_CACHE_LINE);  \
     prefetch_data_trig(px, Y_STEP + 2*Y_CACHE_LINE);\
     prefetch_data_trig(px, Y_STEP + 3*Y_CACHE_LINE);\
     prefetch_data_trig(px, Y_STEP + 4*Y_CACHE_LINE);\
     tw16r= px[30]; tw16i= px[31];         \
     tw17r= px[32]; tw17i= px[33];         \
     tw18r= px[34]; tw18i= px[35];         \
     tw19r= px[36]; tw19i= px[37];         \
     tw20r= px[38]; tw20i= px[39];         \
     tw21r= px[40]; tw21i= px[41];         \
     tw22r= px[42]; tw22i= px[43];         \
     tw23r= px[44]; tw23i= px[45];         \
     tw24r= px[46]; tw24i= px[47];         \
     tw25r= px[48]; tw25i= px[49];         \
     tw26r= px[50]; tw26i= px[51];         \
     tw27r= px[52]; tw27i= px[53];         \
     tw28r= px[54]; tw28i= px[55];         \
     tw29r= px[56]; tw29i= px[57];         \
     tw30r= px[58]; tw30i= px[59];         \
     tw31r= px[60]; tw31i= px[61];         \
     prefetch_data_trig(px, Y_STEP + 5*Y_CACHE_LINE);\
     prefetch_data_trig(px, Y_STEP + 6*Y_CACHE_LINE);\
     prefetch_data_trig(px, Y_STEP + 7*Y_CACHE_LINE);\
     prefetch_data_trig(px, Y_STEP + 8*Y_CACHE_LINE);\



#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 18

# define get_twiddle_factors_32            \
	  tw1r= px[0];                     \
	  tw1i= px[1];                     \
	  tw2r= px[2];                     \
	  tw2i= px[3];                     \
	  tw3r= px[4];                     \
	  tw3i= px[5];                     \
	  tw4r= px[6];                     \
	  tw4i= px[7];                     \
	  tw5r= px[8];                     \
	  tw5i= px[9];                     \
	  tw6r= px[10];                    \
	  tw6i= px[11];                    \
	  tw7r= px[12];                    \
	  tw7i= px[13];                    \
	  tw15r= px[14];                   \
	  tw15i= px[15];                   \
	  cplx_divmul(tw8,tw22,tw15,tw7);  \
	  cplx_divmul(tw9,tw21,tw15,tw6);  \
	  cplx_divmul(tw10,tw20,tw15,tw5); \
	  cplx_divmul(tw11,tw19,tw15,tw4); \
	  cplx_divmul(tw12,tw18,tw15,tw3); \
	  cplx_divmul(tw13,tw17,tw15,tw2); \
	  cplx_divmul(tw14,tw16,tw15,tw1); \
	  tw27r= px[16];                   \
	  tw27i= px[17];                   \
	  cplx_divmul(tw23,tw31,tw27,tw4); \
	  cplx_divmul(tw24,tw30,tw27,tw3); \
	  cplx_divmul(tw25,tw29,tw27,tw2); \
	  cplx_divmul(tw26,tw28,tw27,tw1);


#endif

/*
	Get data from memory, mul by twiddle trig factors,
        make a complex 16-length FFT and store the data to memory 
*/

# define radix_32_twd_pp(_S,_j)                                     \
	  cplx_load_muladdsub_pp(t0,t1,tw16,_j,pd0,pd16);           \
	  cplx_load_mulmuladdsub_pp(t2,t3,tw8,tw24,_j,pd8,pd24);    \
	  cplx_addsub(t0,t2);                                       \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                        \
	  cplx_load_mulmuladdsub_pp(t4,t5,tw4,tw20,_j,pd4,pd20);    \
	  cplx_load_mulmuladdsub_pp(t6,t7,tw12,tw28,_j,pd12,pd28);  \
	  cplx_addsub(t4,t6);                                       \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                        \
	  cplx_addsub(t0,t4);                                       \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                        \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                        \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                        \
	  cplx_load_mulmuladdsub_pp(t8,t9,tw2,tw18,_j,pd2,pd18);    \
	  cplx_load_mulmuladdsub_pp(t10,t11,tw10,tw26,_j,pd10,pd26);\
	  cplx_addsub(t8,t10);                                      \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                       \
	  cplx_load_mulmuladdsub_pp(t12,t13,tw6,tw22,_j,pd6,pd22);  \
	  cplx_load_mulmuladdsub_pp(t14,t15,tw14,tw30,_j,pd14,pd30);\
	  cplx_addsub(t12,t14);                                     \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                      \
	  cplx_addsub(t8,t12);                                      \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                       \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                      \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                      \
	  cplx_addsub(t0,t8);                                       \
	  cplx_muladdsub(t1,t9,_S##_1_16);                          \
	  cplx_mul_1_8_##_S##_addsub(t2,t10);                       \
	  cplx_muladdsub(t3,t11,_S##_3_16);                         \
	  cplx_mul_1_4_##_S##_addsub(t4,t12);                       \
	  cplx_muladdsub(t5,t13,_S##_5_16);                         \
	  cplx_mul_3_8_##_S##_addsub(t6,t14);                       \
	  cplx_muladdsub(t7,t15,_S##_7_16);                         \
	                                                            \
	  cplx_load_mulmuladdsub_pp(t16,t17,tw1,tw17,_j,pd1,pd17);  \
	  cplx_load_mulmuladdsub_pp(t18,t19,tw9,tw25,_j,pd9,pd25);  \
	  cplx_addsub(t16,t18);                                     \
	  cplx_mul_1_4_##_S##_addsub(t17,t19);                      \
	  cplx_load_mulmuladdsub_pp(t20,t21,tw5,tw21,_j,pd5,pd21);  \
	  cplx_load_mulmuladdsub_pp(t22,t23,tw13,tw29,_j,pd13,pd29);\
	  cplx_addsub(t20,t22);                                     \
	  cplx_mul_1_4_##_S##_addsub(t21,t23);                      \
	  cplx_addsub(t16,t20);                                     \
	  cplx_mul_1_8_##_S##_addsub(t17,t21);                      \
	  cplx_mul_1_4_##_S##_addsub(t18,t22);                      \
	  cplx_mul_3_8_##_S##_addsub(t19,t23);                      \
	  cplx_load_mulmuladdsub_pp(t24,t25,tw3,tw19,_j,pd3,pd19);  \
	  cplx_load_mulmuladdsub_pp(t26,t27,tw11,tw27,_j,pd11,pd27);\
	  cplx_addsub(t24,t26);                                     \
	  cplx_mul_1_4_##_S##_addsub(t25,t27);                      \
	  cplx_load_mulmuladdsub_pp(t28,t29,tw7,tw23,_j,pd7,pd23);  \
	  cplx_load_mulmuladdsub_pp(t30,t31,tw15,tw31,_j,pd15,pd31);\
	  cplx_addsub(t28,t30);                                     \
	  cplx_mul_1_4_##_S##_addsub(t29,t31);                      \
	  cplx_addsub(t24,t28);                                     \
	  cplx_mul_1_8_##_S##_addsub(t25,t29);                      \
	  cplx_mul_1_4_##_S##_addsub(t26,t30);                      \
	  cplx_mul_3_8_##_S##_addsub(t27,t31);                      \
	  cplx_addsub(t16,t24);                                     \
	  cplx_muladdsub(t17,t25,_S##_1_16);                        \
	  cplx_mul_1_8_##_S##_addsub(t18,t26);                      \
	  cplx_muladdsub(t19,t27,_S##_3_16);                        \
	  cplx_mul_1_4_##_S##_addsub(t20,t28);                      \
	  cplx_muladdsub(t21,t29,_S##_5_16);                        \
	  cplx_mul_3_8_##_S##_addsub(t22,t30);                      \
	  cplx_muladdsub(t23,t31,_S##_7_16);                        \
	                                                            \
	  cplx_addsub_store_p(_j,pd0,pd16,t0,t16);                  \
	  cplx_muladdsub_store_p(_j,pd1,pd17,t1,t17,_S##_1_32);     \
	  cplx_muladdsub_store_p(_j,pd2,pd18,t2,t18,_S##_1_16);     \
	  cplx_muladdsub_store_p(_j,pd3,pd19,t3,t19,_S##_3_32);     \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd4,pd20,t4,t20);   \
	  cplx_muladdsub_store_p(_j,pd5,pd21,t5,t21,_S##_5_32);     \
	  cplx_muladdsub_store_p(_j,pd6,pd22,t6,t22,_S##_3_16);     \
	  cplx_muladdsub_store_p(_j,pd7,pd23,t7,t23,_S##_7_32);     \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd8,pd24,t8,t24);   \
	  cplx_muladdsub_store_p(_j,pd9,pd25,t9,t25,_S##_9_32);     \
	  cplx_muladdsub_store_p(_j,pd10,pd26,t10,t26,_S##_5_16);   \
	  cplx_muladdsub_store_p(_j,pd11,pd27,t11,t27,_S##_11_32);  \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd12,pd28,t12,t28); \
	  cplx_muladdsub_store_p(_j,pd13,pd29,t13,t29,_S##_13_32);  \
	  cplx_muladdsub_store_p(_j,pd14,pd30,t14,t30,_S##_7_16);   \
	  cplx_muladdsub_store_p(_j,pd15,pd31,t15,t31,_S##_15_32);



/*
	Get data from memory, make a complex 8-length FFT and
	store the data to memory 
  		  
*/
# define radix_32_notwd_pp(_S,_j)                                  \
	  cplx_load_addsub_pp(t0,t1,_j,pd0,pd16);                  \
	  cplx_load_addsub_pp(t2,t3,_j,pd8,pd24);                  \
	  cplx_addsub(t0,t2);                                      \
	  cplx_mul_1_4_##_S##_addsub(t1,t3);                       \
	  cplx_load_addsub_pp(t4,t5,_j,pd4,pd20);                  \
	  cplx_load_addsub_pp(t6,t7,_j,pd12,pd28);                 \
	  cplx_addsub(t4,t6);                                      \
	  cplx_mul_1_4_##_S##_addsub(t5,t7);                       \
	  cplx_addsub(t0,t4);                                      \
	  cplx_mul_1_8_##_S##_addsub(t1,t5);                       \
	  cplx_mul_1_4_##_S##_addsub(t2,t6);                       \
	  cplx_mul_3_8_##_S##_addsub(t3,t7);                       \
	  cplx_load_addsub_pp(t8,t9,_j,pd2,pd18);                  \
	  cplx_load_addsub_pp(t10,t11,_j,pd10,pd26);               \
	  cplx_addsub(t8,t10);                                     \
	  cplx_mul_1_4_##_S##_addsub(t9,t11);                      \
	  cplx_load_addsub_pp(t12,t13,_j,pd6,pd22);                \
	  cplx_load_addsub_pp(t14,t15,_j,pd14,pd30);               \
	  cplx_addsub(t12,t14);                                    \
	  cplx_mul_1_4_##_S##_addsub(t13,t15);                     \
	  cplx_addsub(t8,t12);                                     \
	  cplx_mul_1_8_##_S##_addsub(t9,t13);                      \
	  cplx_mul_1_4_##_S##_addsub(t10,t14);                     \
	  cplx_mul_3_8_##_S##_addsub(t11,t15);                     \
	  cplx_addsub(t0,t8);                                      \
	  cplx_muladdsub(t1,t9,_S##_1_16);                         \
	  cplx_mul_1_8_##_S##_addsub(t2,t10);                      \
	  cplx_muladdsub(t3,t11,_S##_3_16);                        \
	  cplx_mul_1_4_##_S##_addsub(t4,t12);                      \
	  cplx_muladdsub(t5,t13,_S##_5_16);                        \
	  cplx_mul_3_8_##_S##_addsub(t6,t14);                      \
	  cplx_muladdsub(t7,t15,_S##_7_16);                        \
	                                                           \
	  cplx_load_addsub_pp(t16,t17,_j,pd1,pd17);                \
	  cplx_load_addsub_pp(t18,t19,_j,pd9,pd25);                \
	  cplx_addsub(t16,t18);                                    \
	  cplx_mul_1_4_##_S##_addsub(t17,t19);                     \
	  cplx_load_addsub_pp(t20,t21,_j,pd5,pd21);                \
	  cplx_load_addsub_pp(t22,t23,_j,pd13,pd29);               \
	  cplx_addsub(t20,t22);                                    \
	  cplx_mul_1_4_##_S##_addsub(t21,t23);                     \
	  cplx_addsub(t16,t20);                                    \
	  cplx_mul_1_8_##_S##_addsub(t17,t21);                     \
	  cplx_mul_1_4_##_S##_addsub(t18,t22);                     \
	  cplx_mul_3_8_##_S##_addsub(t19,t23);                     \
	  cplx_load_addsub_pp(t24,t25,_j,pd3,pd19);                \
	  cplx_load_addsub_pp(t26,t27,_j,pd11,pd27);               \
	  cplx_addsub(t24,t26);                                    \
	  cplx_mul_1_4_##_S##_addsub(t25,t27);                     \
	  cplx_load_addsub_pp(t28,t29,_j,pd7,pd23);                \
	  cplx_load_addsub_pp(t30,t31,_j,pd15,pd31);               \
	  cplx_addsub(t28,t30);                                    \
	  cplx_mul_1_4_##_S##_addsub(t29,t31);                     \
	  cplx_addsub(t24,t28);                                    \
	  cplx_mul_1_8_##_S##_addsub(t25,t29);                     \
	  cplx_mul_1_4_##_S##_addsub(t26,t30);                     \
	  cplx_mul_3_8_##_S##_addsub(t27,t31);                     \
	  cplx_addsub(t16,t24);                                    \
	  cplx_muladdsub(t17,t25,_S##_1_16);                       \
	  cplx_mul_1_8_##_S##_addsub(t18,t26);                     \
	  cplx_muladdsub(t19,t27,_S##_3_16);                       \
	  cplx_mul_1_4_##_S##_addsub(t20,t28);                     \
	  cplx_muladdsub(t21,t29,_S##_5_16);                       \
	  cplx_mul_3_8_##_S##_addsub(t22,t30);                     \
	  cplx_muladdsub(t23,t31,_S##_7_16);                       \
	                                                           \
	  cplx_addsub_store_p(_j,pd0,pd16,t0,t16);                 \
	  cplx_muladdsub_store_p(_j,pd1,pd17,t1,t17,_S##_1_32);    \
	  cplx_muladdsub_store_p(_j,pd2,pd18,t2,t18,_S##_1_16);    \
	  cplx_muladdsub_store_p(_j,pd3,pd19,t3,t19,_S##_3_32);    \
	  cplx_mul_1_8_##_S##_addsub_store_p(_j,pd4,pd20,t4,t20);  \
	  cplx_muladdsub_store_p(_j,pd5,pd21,t5,t21,_S##_5_32);    \
	  cplx_muladdsub_store_p(_j,pd6,pd22,t6,t22,_S##_3_16);    \
	  cplx_muladdsub_store_p(_j,pd7,pd23,t7,t23,_S##_7_32);    \
	  cplx_mul_1_4_##_S##_addsub_store_p(_j,pd8,pd24,t8,t24);  \
	  cplx_muladdsub_store_p(_j,pd9,pd25,t9,t25,_S##_9_32);    \
	  cplx_muladdsub_store_p(_j,pd10,pd26,t10,t26,_S##_5_16);  \
	  cplx_muladdsub_store_p(_j,pd11,pd27,t11,t27,_S##_11_32); \
	  cplx_mul_3_8_##_S##_addsub_store_p(_j,pd12,pd28,t12,t28);\
	  cplx_muladdsub_store_p(_j,pd13,pd29,t13,t29,_S##_13_32); \
	  cplx_muladdsub_store_p(_j,pd14,pd30,t14,t30,_S##_7_16);  \
	  cplx_muladdsub_store_p(_j,pd15,pd31,t15,t31,_S##_15_32);

#if (Y_AVAL > 4)

# define RADIX 32

/*************************************************************************
   This is a prototype of an Inplace  Forward Decimation in Frecuency 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = Scrambled Twidle factors. Not padded.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radixmp_32_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,t21r,t21i,t22r,
  t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,
  t29i,t30r,t30i,t31r,t31i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i,tw16r,tw16i,tw17r,tw17i,tw18r,
  tw18i,tw19r,tw19i,tw20r,tw20i,tw21r,tw21i,tw22r,tw22i,tw23r,tw23i,
  tw24r,tw24i,tw25r,tw25i,tw26r,tw26i,tw27r,tw27i,tw28r,tw28i,tw29r,
  tw29i,tw30r,tw30i,tw31r,tw31i;
  y_ptr px;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,
  pd14,pd15,pd16,pd17,pd18,pd19,pd20,pd21,pd22,pd23,pd24,pd25,
  pd26,pd27,pd28,pd29,pd30,pd31;
  y_size_t i,j,jj,nc;
  y_size_t done, mask, ioff, ibase, n00;


  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();
# ifndef NDEBUG1

  printf(" Radixmp_32_dif n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  bigpad=(pad<<5);
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad);
  ioff = (n0 & mask) >> 5;

  get_pads_32;
  /* If ibase == 0 first elements has twiddle factors equal to one */
  if(ibase == 0)
    {
      for(j = ioff;j < pad && done < n; j++)
        {
          jj = addr(j<<1);
          radix_32_notwd_pp(F,jj);
          done += RADIX;
        }
      if(done == n)
        return;
      px = tw + Y_STEP;
      for(i = bigpad; done < n; i += bigpad, px += Y_STEP)
        {
          get_twiddle_factors_32;
          for(j = i; j < (pad + i) && done < n; j++)
            {
              /*BEGIN_TIME;*/
              jj = addr(j<<1);
              radix_32_twd_pp(F,jj);
              done += RADIX;
              /*END_TIME;*/
            }
        }
    }
  else
    {
      px = tw + Y_STEP * ibase;
      for(i = ibase * bigpad; done < n; i += bigpad, px += Y_STEP)
        {
          get_twiddle_factors_32;
          for(j = i + ioff; j < (pad + i) && done < n; j++)
            {
              /*BEGIN_TIME;*/
              jj = addr(j<<1);
              radix_32_twd_pp(F,jj);
              done += RADIX;
              /*END_TIME;*/
            }
          ioff = 0;
        }
    }
}



/*************************************************************************
   This is a prototype of an Inplace Backward Decimation in Time 
   Fast Numeric Transform Radix - r 
   INPUTS:
   d[] =all the data. Because of padding, d[0] must be the first
           element of the whole data array. 
   tw[] = Twidle Backward factors. Not scrambled.
   n = number of data to transform in this call. 
   n0 = index of first data.
   pad = the pad (logical) between data.
*************************************************************************/

void radixmp_32_dit(y_ptr d, y_ptr tw, y_size_t n , y_size_t n0, y_size_t pad)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,
  t8r,t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,
  t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,t21r,t21i,t22r,
  t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,
  t29i,t30r,t30i,t31r,t31i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i,tw16r,tw16i,tw17r,tw17i,tw18r,
  tw18i,tw19r,tw19i,tw20r,tw20i,tw21r,tw21i,tw22r,tw22i,tw23r,tw23i,
  tw24r,tw24i,tw25r,tw25i,tw26r,tw26i,tw27r,tw27i,tw28r,tw28i,tw29r,
  tw29i,tw30r,tw30i,tw31r,tw31i;
  y_size_t bigpad;
  y_ptr pd0,pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,pd13,
  pd14,pd15,pd16,pd17,pd18,pd19,pd20,pd21,pd22,pd23,pd24,pd25,
  pd26,pd27,pd28,pd29,pd30,pd31;

  y_ptr px;
  y_size_t i, j, jj;
  y_size_t done, mask, ioff, ibase;

  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

# ifndef NDEBUG1

  printf(" Radixmp_32_dit n=%i n0=%i pad=%i \n",n,n0,pad);
# endif

  bigpad = (pad << 5);
  done = 0;
  mask = bigpad - 1;
  ibase = (n0 / bigpad) * bigpad;
  ioff = (n0 & mask) >> 5;

  get_pads_32;

  for(i = ibase; done < n; i += bigpad)
    {
      if (ioff == 0)
        {
          jj = addr(i<<1);
          radix_32_notwd_pp(B,jj);
          done += RADIX;
          j = i + 1;
          px = tw;
        }
      else
        {
          j = i + ioff;
          px = tw + (ioff - 1) * Y_STEP;
        }
      for( ;j < (pad + i) && done < n; j++, px += Y_STEP)
        {
          get_twiddle_factors_32;
          jj = addr(j<<1);
          radix_32_twd_pp(B,jj);
          done += RADIX;
        }
      ioff = 0;
    }
}


#else

void radixmp_32_dif(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radixmp_32_dif , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_32.c with -DY_AVAL=5 flag\n");
  exit(-1);
}

void radixmp_32_dit(y_ptr d, y_ptr tw, y_size_t n, y_size_t n0, y_size_t pad)
{
  printf("Called radixmp_32_dit , d=%p, tw=%p, n=%d, n0=%d, pad=%d\n",d,tw,n,n0,pad);
  printf("Please compile again the file radix_32.c with -DY_AVAL=5 flag\n");
  exit(-1);
}

#endif


/*$Id$*/

















