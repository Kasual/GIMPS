/*$Id$*/
/*
   (c) 2000-2006 Guillermo Ballester Valor 
   
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
   
   It includes ideas from many Authors.
 
*/
/* this is to mark the code as norm-carry pass */
#define YNORM
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#include "glucas.h"
#include "gsetup.h"
#include "ynorm.h"
#ifdef Y_ITANIUM
# include "yia64.h"
# include "ynormia64.h"
# ifndef Y_PRELOAD_TRIGS
#  define Y_PRELOAD_TRIGS
# endif
# ifndef Y_PRELOAD_DATA
#  define Y_PRELOAD_DATA
# endif
# ifdef Y_LONG_MACROS
#  undef Y_LONG_MACROS
# endif
# ifdef Y_KILL_BRANCHES
#  undef Y_KILL_BRANCHES
# endif
# ifndef Y_VECTORIZE
#  define Y_VECTORIZE
# endif
#else
# ifdef Y_LONG_MACROS
#  include "ygeneric8.h"
# endif
#endif


#define get_pads_8       \
  pad2=(pad<<1);         \
  pad4=(pad<<2);         \
  pad3=pad+pad2;         \
  pad5=pad+pad4;         \
  pad6=pad3+pad3;        \
  pad7=pad4+pad3;        \
  pd0= x;                \
  prefetch_data( pd0, 0);\
  pd1= x + addr(pad<<1); \
  prefetch_data( pd1, 0);\
  pd2= x + addr(pad2<<1);\
  prefetch_data( pd2, 0);\
  pd3= x + addr(pad3<<1);\
  prefetch_data( pd3, 0);\
  pd4= x + addr(pad4<<1);\
  prefetch_data( pd4, 0);\
  pd5= x + addr(pad5<<1);\
  prefetch_data( pd5, 0);\
  pd6= x + addr(pad6<<1);\
  prefetch_data( pd6, 0);\
  pd7= x + addr(pad7<<1);\
  prefetch_data( pd7, 0);



#define init_bjs_8 \
  bj0=N;\
  bj1=(b & (UL)7 ) <<(Y_K+1); \
  bj2=((2*b) & (UL)7 ) <<(Y_K+1); \
  bj3=((3*b) & (UL)7 ) <<(Y_K+1); \
  bj4=((4*b) & (UL)7 ) <<(Y_K+1); \
  bj5=((5*b) & (UL)7 ) <<(Y_K+1); \
  bj6=((6*b) & (UL)7 ) <<(Y_K+1); \
  bj7=((7*b) & (UL)7 ) <<(Y_K+1);

#define init_twiddle_factors_8   \
  tw1r=1.0; tw1i=0.0;            \
  tw2r=1.0; tw2i=0.0;            \
  tw3r=1.0; tw3i=0.0;            \
  tw4r=1.0; tw4i=0.0;            \
  tw5r=1.0; tw5i=0.0;            \
  tw6r=1.0; tw6i=0.0;            \
  tw7r=1.0; tw7i=0.0;


#define init_twiddle_factors_8_t \
  t1r=1.0; t1i=0.0;              \
  t2r=1.0; t2i=0.0;              \
  t3r=1.0; t3i=0.0;              \
  t4r=1.0; t4i=0.0;              \
  t5r=1.0; t5i=0.0;              \
  t6r=1.0; t6i=0.0;              \
  t7r=1.0; t7i=0.0;

#ifdef Y_KILL_BRANCHES

# define  load_dwf_8(_i)                                                    \
	tt0[0] = two_to_phi[addr(_i)];                                      \
	tt0[1] = two_to_minusphi[addr(_i)];                                 \
	tt1[0] = two_to_phi[addr(((pad2 >> SHIFT_UPDATE) + _i))];           \
	tt1[1] = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];      \
	tt2[0] = two_to_phi[addr(((pad4 >> SHIFT_UPDATE) + _i))];           \
	tt2[1] = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];      \
	tt3[0] = two_to_phi[addr(((pad6 >> SHIFT_UPDATE) + _i))];           \
	tt3[1] = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];      \
	tt4[0] = two_to_phi[addr(((pad4 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt4[1] = two_to_minusphi[addr(((pad4>>(SHIFT_UPDATE - 1)) + _i))];  \
	tt5[0] = two_to_phi[addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt5[1] = two_to_minusphi[addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt6[0] = two_to_phi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt6[1] = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];\
	tt7[0] = two_to_phi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];     \
	tt7[1] = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];

# define  load_two_to_minusphi_8(_i) \
	tt0[1] = two_to_minusphi[addr(_i)]; \
	tt1[1] = two_to_minusphi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	tt2[1] = two_to_minusphi[addr(((pad4>>SHIFT_UPDATE) + _i))];\
	tt3[1] = two_to_minusphi[addr(((pad6>>SHIFT_UPDATE) + _i))];\
	tt4[1] = two_to_minusphi[addr(((pad4>>(SHIFT_UPDATE-1)) + _i))];\
	tt5[1] = two_to_minusphi[addr(((pad5>>(SHIFT_UPDATE-1)) + _i))];\
	tt6[1] = two_to_minusphi[addr(((pad6>>(SHIFT_UPDATE-1)) + _i))];\
	tt7[1] = two_to_minusphi[addr(((pad7>>(SHIFT_UPDATE-1)) + _i))];

# define  load_two_to_phi_8(_i) \
	tt0[0] = two_to_phi[addr(_i)]; \
	tt1[0] = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	tt2[0] = two_to_phi[addr(((pad4>>SHIFT_UPDATE) + _i))];\
	tt3[0] = two_to_phi[addr(((pad6>>SHIFT_UPDATE) + _i))];\
	tt4[0] = two_to_phi[addr(((pad4>>(SHIFT_UPDATE-1)) + _i))];\
	tt5[0] = two_to_phi[addr(((pad5>>(SHIFT_UPDATE-1)) + _i))];\
	tt6[0] = two_to_phi[addr(((pad6>>(SHIFT_UPDATE-1)) + _i))];\
	tt7[0] = two_to_phi[addr(((pad7>>(SHIFT_UPDATE-1)) + _i))];


#else

# define  load_dwf_8(_i)                                                   \
	ttp0 = two_to_phi[addr(_i)];                                       \
	ttmp0 = two_to_minusphi[addr(_i)];                                 \
        prefetch_data_trig(two_to_minusphi, addr(_i) + Y_CACHE_LINE);      \
        prefetch_data_trig(two_to_phi, addr(_i) + Y_CACHE_LINE);           \
	ttp1 = two_to_phi[addr(((pad2 >> SHIFT_UPDATE) + _i))];            \
	ttmp1 = two_to_minusphi[addr(((pad2 >> SHIFT_UPDATE) + _i))];      \
        prefetch_data_trig(two_to_minusphi, addr(((pad2 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad2 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
	ttp2 = two_to_phi[addr(((pad4 >> SHIFT_UPDATE) + _i))];            \
	ttmp2 = two_to_minusphi[addr(((pad4 >> SHIFT_UPDATE) + _i))];      \
        prefetch_data_trig(two_to_minusphi, addr(((pad4 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad4 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
	ttp3 = two_to_phi[addr(((pad6 >> SHIFT_UPDATE) + _i))];            \
	ttmp3 = two_to_minusphi[addr(((pad6 >> SHIFT_UPDATE) + _i))];      \
        prefetch_data_trig(two_to_minusphi, addr(((pad6 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad6 >> SHIFT_UPDATE) + _i)) + Y_CACHE_LINE);\
	ttp4 = two_to_phi[addr(((pad4 >> (SHIFT_UPDATE - 1)) + _i))];      \
	ttmp4 = two_to_minusphi[addr(((pad4 >> (SHIFT_UPDATE - 1)) + _i))];\
        prefetch_data_trig(two_to_minusphi, addr(((pad4 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad4 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
	ttp5 = two_to_phi[addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i))];      \
	ttmp5 = two_to_minusphi[addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i))];\
        prefetch_data_trig(two_to_minusphi, addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad5 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
	ttp6 = two_to_phi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];      \
	ttmp6 = two_to_minusphi[addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i))];\
        prefetch_data_trig(two_to_minusphi, addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad6 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
	ttp7 = two_to_phi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))];      \
	ttmp7 = two_to_minusphi[addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i))]; \
        prefetch_data_trig(two_to_minusphi, addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);\
        prefetch_data_trig(two_to_phi, addr(((pad7 >> (SHIFT_UPDATE - 1)) + _i)) + Y_CACHE_LINE);


# define  load_two_to_minusphi_8(_i) \
	ttmp0 = two_to_minusphi[addr(_i)]; \
	ttmp1 = two_to_minusphi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	ttmp2 = two_to_minusphi[addr(((pad4>>SHIFT_UPDATE) + _i))];\
	ttmp3 = two_to_minusphi[addr(((pad6>>SHIFT_UPDATE) + _i))];\
	ttmp4 = two_to_minusphi[addr(((pad4>>(SHIFT_UPDATE-1)) + _i))];\
	ttmp5 = two_to_minusphi[addr(((pad5>>(SHIFT_UPDATE-1)) + _i))];\
	ttmp6 = two_to_minusphi[addr(((pad6>>(SHIFT_UPDATE-1)) + _i))];\
	ttmp7 = two_to_minusphi[addr(((pad7>>(SHIFT_UPDATE-1)) + _i))];

# define  load_two_to_phi_8(_i) \
	ttp0 = two_to_phi[addr(_i)]; \
	ttp1 = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	ttp2 = two_to_phi[addr(((pad4>>SHIFT_UPDATE) + _i))];\
	ttp3 = two_to_phi[addr(((pad6>>SHIFT_UPDATE) + _i))];\
	ttp4 = two_to_phi[addr(((pad4>>(SHIFT_UPDATE-1)) + _i))];\
	ttp5 = two_to_phi[addr(((pad5>>(SHIFT_UPDATE-1)) + _i))];\
	ttp6 = two_to_phi[addr(((pad6>>(SHIFT_UPDATE-1)) + _i))];\
	ttp7 = two_to_phi[addr(((pad7>>(SHIFT_UPDATE-1)) + _i))];

#endif

/* get the basic twiddle from memory and computes its powers */
#ifdef Y_MINIMUM
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2

# include "yminimum.h"

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)              \
          tpw1r = *(_px );                   \
          tpw1i = *(_px + 1);

#  define trig_8_preload(_px)                \
          tpw1r = *(_px + Y_STEP);           \
          tpw1i = *(_px + Y_STEP + 1);

#  define get_twiddle_factors_8_preload      \
          tw1r = tpw1r; tw1i = tpw1i;        \
          trig_8_preload (px);               \
          cplx_trig_min_square (tw2, tw1);   \
          cplx_trig_min_square (tw4, tw2);   \
          cplx_divmul (tw3, tw5, tw4, tw1);  \
          cplx_trig_min_square (tw6, tw3);   \
          cplx_trig_min_mul (tw7, tw4, tw3);

#  define get_twiddle_factors_8_preload_consumed   \
          t1r = tpw1r; t1i = tpw1i;                \
          trig_8_preload (px);                     \
          cplx_trig_min_square (t2, t1);           \
          cplx_trig_min_square (t4, t2);           \
          cplx_divmul (t3, t5, t4, t1);            \
          cplx_trig_min_square (t6, t3);           \
          cplx_trig_min_mul (t7, t4, t3);

# else

#  define get_twiddle_factors_8             \
          cplx_trig_min_load (tw1, px);     \
          cplx_trig_min_square (tw2, tw1);  \
          cplx_trig_min_square (tw4, tw2);  \
          cplx_divmul (tw3, tw5, tw4, tw1); \
          cplx_trig_min_square (tw6, tw3);  \
          cplx_trig_min_mul (tw7, tw4, tw3);\
          prefetch_p_trig (px);

# endif

#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 14

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)                          \
          tpw1r = *(_px);         tpw1i = *(_px + 1);    \
          tpw2r = *(_px + 2);     tpw2i = *(_px + 3);    \
          tpw3r = *(_px + 4);     tpw3i = *(_px + 5);    \
          tpw4r = *(_px + 6);     tpw4i = *(_px + 7);    \
          tpw5r = *(_px + 8);     tpw5i = *(_px + 9);    \
          tpw6r = *(_px + 10);    tpw6i = *(_px + 11);   \
          tpw7r = *(_px + 12);    tpw7i = *(_px + 13);

#  define trig_8_preload(_px)                                              \
          tpw1r = *(_px + Y_STEP);         tpw1i = *(_px + Y_STEP + 1);    \
          tpw2r = *(_px + Y_STEP + 2);     tpw2i = *(_px + Y_STEP + 3);    \
          tpw3r = *(_px + Y_STEP + 4);     tpw3i = *(_px + Y_STEP + 5);    \
          tpw4r = *(_px + Y_STEP + 6);     tpw4i = *(_px + Y_STEP + 7);    \
          tpw5r = *(_px + Y_STEP + 8);     tpw5i = *(_px + Y_STEP + 9);    \
          tpw6r = *(_px + Y_STEP + 10);    tpw6i = *(_px + Y_STEP + 11);   \
          tpw7r = *(_px + Y_STEP + 12);    tpw7i = *(_px + Y_STEP + 13);

#  define get_twiddle_factors_8_preload       \
   cplx_complete_twiddle_8_ia64(tw,tpw);      \
   trig_8_preload(px);                        \

#  define get_twiddle_factors_8_preload_consumed   \
   cplx_complete_twiddle_8_ia64(t,tpw);            \
   trig_8_preload(px);                             \

# else

#  define get_twiddle_factors_8                                \
	  tw1r = px[0];  tw1i = px[1];                         \
	  tw2r = px[2];  tw2i = px[3];                         \
	  tw3r = px[4];  tw3i = px[5];                         \
	  tw4r = px[6];  tw4i = px[7];                         \
	  tw5r = px[8];  tw5i = px[9];                         \
	  tw6r = px[10];  tw6i = px[11];                       \
	  tw7r = px[12];  tw7i = px[13];                       \
          prefetch_data_trig (px, 2 * Y_STEP);                 \
          prefetch_data_trig (px, 2 * Y_STEP + Y_CACHE_LINE);  \
          prefetch_data_trig (px, 2 * Y_STEP + 2 * Y_CACHE_LINE);

# endif

#else
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 6

# ifdef Y_PRELOAD_TRIGS

#  define trig_8_load_init(_px)      \
          tpw1r = *(_px);            \
          tpw1i = *(_px + 1);        \
          tpw2r = *(_px + 2);        \
          tpw2i = *(_px + 3);        \
          tpw5r = *(_px + 4);        \
          tpw5i = *(_px + 5);        \


#  define trig_8_preload(_px)                 \
          tpw1r = *(_px + Y_STEP);            \
          tpw1i = *(_px + Y_STEP + 1);        \
          tpw2r = *(_px + Y_STEP + 2);        \
          tpw2i = *(_px + Y_STEP + 3);        \
          tpw5r = *(_px + Y_STEP + 4);        \
          tpw5i = *(_px + Y_STEP + 5);        \

#  define get_twiddle_factors_8_preload       \
   cplx_complete_twiddle_8_ia64 (tw, tpw);    \
   trig_8_preload (px);                       \

#  define get_twiddle_factors_8_preload_consumed   \
   cplx_complete_twiddle_8_ia64 (t, tpw);          \
   trig_8_preload(px);                             \

# else

#  define get_twiddle_factors_8                        \
	  tw1r = px[0];                                \
	  tw1i = px[1];                                \
          tw2r = px[2];                                \
          tw2i = px[3];                                \
          tw5r = px[4];                                \
	  tw5i = px[5];                                \
          prefetch_data_trig (px, Y_STEP);             \
	  cplx_divmul (tw3, tw7, tw5, tw2);            \
	  cplx_divmul (tw4, tw6, tw5, tw1);            \
          prefetch_data_trig (px, Y_STEP + Y_CACHE_LINE);

# endif

#endif

#ifdef Y_PRELOAD_DATA
# define radix_8_twd_last_dit(_k)                                          \
    t0r = s0r; t0i = s0i; ppd0r = ppdl0r; ppd0i = ppdl0i;                  \
    ppdl0r += _k; ppdl0i += _k; s0r= *(ppdl0r); s0i=*(ppdl0i);              \
    cplx_three_muls_fftn_preload_consumed_trigs (t2, t5, t7, 2, 5, 7, s, t, ppd, ppdl, _k);\
    cplx_four_muls_fftn_preload_consumed_trigs_paired(t1, t4, t3, t6, 1, 4, 3, 6, s, t, ppd, ppdl, _k);\
    cplx_fft8_ia64_B(t, 0, 1, 2, 3, 4, 5, 6, 7);


# define radix_8_first_dif(_j)                                     \
    cplx_fft8_ia64_store_F (ppd0, ppd1, ppd2, ppd3, ppd4, ppd5, ppd6, ppd7, t, 0, 4, 2, 6, 1, 5, 3, 7);

#elif defined(Y_LONG_MACROS)

# define radix_8_twd_last_dit(_j)                                  \
        cplx_radix_8_twd_pp_last_dit_pass1 (t, tw, _j, pd);        \
        cplx_radix_8_pp_pass2_B (t);                               \
        cplx_radix_8_pp_last_dit_pass3 (t);

# define radix_8_first_dif(_j)                                     \
        cplx_radix_8_pp_first_dif_pass1 (t);                       \
        cplx_radix_8_pp_first_dif_pass2 (t);                       \
        cplx_radix_8_pp_first_dif_pass3 (t, pd, _j);

#else

# define radix_8_twd_last_dit(_j)                                  \
        cplx_load_muladdsub_pp (t0, t1, tw4, _j, pd0, pd4);        \
	cplx_load_mulmuladdsub_pp (t2, t3, tw2, tw6, _j, pd2, pd6);\
	cplx_load_mulmuladdsub_pp (t4, t5, tw1, tw5, _j, pd1, pd5);\
	cplx_load_mulmuladdsub_pp (t6, t7, tw3, tw7, _j, pd3, pd7);\
	                                                           \
	cplx_addsub (t0, t2);                                      \
	cplx_mul_1_4_B_addsub (t1, t3);                            \
	cplx_addsub (t4, t6);                                      \
	cplx_mul_1_4_B_addsub (t5, t7);                            \
	                                                           \
	cplx_addsub (t0, t4);                                      \
	cplx_mul_1_8_B_addsub (t1, t5);                            \
	cplx_mul_1_4_B_addsub (t2, t6);                            \
	cplx_mul_3_8_B_addsub (t3, t7);

# define radix_8_first_dif(_j)                               \
	cplx_addsub (t0, t4);                                \
	cplx_addsub (t2, t6);                                \
	cplx_addsub (t0, t2);                                \
	cplx_mul_1_4_F_addsub (t4, t6);                      \
	cplx_addsub (t1, t5);                                \
	cplx_addsub (t3, t7);                                \
	cplx_addsub (t1, t3);                                \
	cplx_mul_1_4_F_addsub (t5, t7);                      \
	cplx_addsub_store_p (_j, pd0, pd4, t0, t1);          \
        cplx_mul_1_8_F_addsub_store_p (_j, pd1, pd5, t4, t5);\
        cplx_mul_1_4_F_addsub_store_p (_j, pd2, pd6, t2, t3);\
        cplx_mul_3_8_F_addsub_store_p (_j, pd3, pd7, t6, t7);\


#endif



#define radix_8_first_dif_add(_d)             \
	cplx_addsub (t0, t4);                 \
	cplx_addsub (t2, t6);                 \
	cplx_addsub (t0, t2);                 \
	cplx_mul_1_4_F_addsub (t4, t6);       \
	cplx_addsub (t1, t5);                 \
	cplx_addsub (t3, t7);                 \
	cplx_addsub (t1, t3);                 \
	cplx_mul_1_4_F_addsub (t5, t7);       \
	cplx_addsub (t0, t1);                 \
	cplx_local_to_data_add (_d, 0, t0);   \
	cplx_local_to_data_add (_d, pad4, t1);\
	cplx_mul_1_8_F_addsub (t4, t5);       \
	cplx_local_to_data_add (_d, pad, t4); \
	cplx_local_to_data_add (_d, pad5, t5);\
	cplx_mul_1_4_F_addsub (t2, t3);       \
	cplx_local_to_data_add (_d, pad2, t2);\
	cplx_local_to_data_add (_d, pad6, t3);\
	cplx_mul_3_8_F_addsub (t6, t7);       \
	cplx_local_to_data_add (_d, pad3, t6);\
	cplx_local_to_data_add (_d, pad7, t7);

#define radix_8_first_dif_add_shift(_d, _ofs) \
	cplx_addsub (t0, t4);                 \
	cplx_addsub (t2, t6);                 \
	cplx_addsub (t0, t2);                 \
	cplx_mul_1_4_F_addsub (t4, t6);       \
	cplx_addsub (t1, t5);                 \
	cplx_addsub (t3, t7);                 \
	cplx_addsub (t1, t3);                 \
	cplx_mul_1_4_F_addsub (t5, t7);       \
	cplx_addsub (t0, t1);                 \
	cplx_local_to_data_add (_d, _ofs, t0);\
	cplx_local_to_data_add (_d, pad4, t1);\
	cplx_mul_1_8_F_addsub (t4, t5);       \
	cplx_local_to_data_add (_d, pad, t4); \
	cplx_local_to_data_add (_d, pad5, t5);\
	cplx_mul_1_4_F_addsub (t2,t3);        \
	cplx_local_to_data_add (_d, pad2, t2);\
	cplx_local_to_data_add (_d, pad6, t3);\
	cplx_mul_3_8_F_addsub (t6,t7);        \
	cplx_local_to_data_add (_d, pad3, t6);\
	cplx_local_to_data_add (_d, pad7, t7);



#ifdef SUM_CHECK
#define sum_8(_Sum)                                                       \
     _Sum += (((t0r + t0i) + (t1r + t1i)) + ((t2r + t2i) + (t3r + t3i)) + \
	      ((t4r + t4i) + (t5r + t5i)) + ((t6r + t6i) + (t7r + t7i)));
#endif


void substract_two_8(BIG_DOUBLE *x, UL N)
{
  double t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i,
  t5r, t5i, t6r, t6i, t7r, t7i, xtwo;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4, pad5, pad6, pad7, ofs, ni, who,
           ib, a;

  /* the bit to substract */
  ib = Y_SBIT +1;
  if(ib == Y_EXPONENT)
    ib = 0;
  /* The element ni where bit 'ib' is */
  ni=floor((BIG_DOUBLE)ib / Y_XBITS);
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
    }

  /* Sumcheck */
#if defined(SUM_CHECK)
  SumIn += xtwo;
#endif

  /* prepare pointers */
  pad2 = (pad << 1);
  pad4 = (pad << 2);
  pad3 = pad + pad2;
  pad5 = pad + pad4;
  pad6 = pad3 + pad3;
  pad7 = pad4 + pad3;
  pad += ofs;
  pad2 += ofs;
  pad3 += ofs;
  pad4 += ofs;
  pad5 += ofs;
  pad6 += ofs;
  pad7 += ofs;

  /* make first DIF and add to array*/
  radix_8_first_dif_add_shift (x, ofs);

  /* New Y_SBIT and Y_SBJ*/
  Y_SBIT += Y_SBIT;
  if(Y_SBIT >= Y_EXPONENT)
    Y_SBIT -= Y_EXPONENT;
}




#if !defined(_OPENMP) && !defined(_SUNMP) && !defined(_PTHREADS)

void dit_carry_norm_dif_8(BIG_DOUBLE *x, UL N ,UL err_flag)
{
# ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6], maxerr = 0.0;
  double Y_ALIGNED(16) tt0[2], tt1[2], tt2[2], tt3[2], tt4[2], tt5[2], tt6[2], tt7[2];
# else

  BIG_DOUBLE hiinv = highinv, loinv = lowinv;
  BIG_DOUBLE maxerr = 0.0, ttmpSmall = Gsmall, ttmpBig = Gbig,
                                       ttpSmall = Hsmall;
  BIG_DOUBLE ttmp0, ttmp1, ttmp2, ttmp3, ttmp4, ttmp5, ttmp6, ttmp7;
  BIG_DOUBLE ttp0, ttp1, ttp2, ttp3, ttp4, ttp5, ttp6, ttp7;
# endif
# if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A = bigA, B = bigB;
# endif

  BIG_DOUBLE carry0 = 0.0, carry1 = 0.0, carry2 = 0.0, carry3 = 0.0,
                                    carry4 = 0.0, carry5 = 0.0, carry6 = 0.0, carry7 = 0.0;
# ifndef Y_PRELOAD_DATA

  BIG_DOUBLE tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i, tw5r, tw5i,
  tw6r, tw6i, tw7r, tw7i;
# endif

  BIG_DOUBLE t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i,
  t6r, t6i, t7r, t7i;
# ifdef Y_ITANIUM

  y_limb_t s0r, s0i, s1r, s1i, s2r, s2i, s3r, s3i, s4r, s4i, s5r, s5i, s6r,
  s6i, s7r, s7i;
  y_limb_t fr = F_1_8r, fi = F_1_8i;
# endif
# ifdef Y_PRELOAD_TRIGS
# if defined(Y_MAXIMUM)

  BIG_DOUBLE tpw1r, tpw1i, tpw2r, tpw2i, tpw3r, tpw3i, tpw4r, tpw4i, tpw5r,
  tpw5i, tpw6r, tpw6i, tpw7r, tpw7i;
# elif defined(Y_MINIMUM)

  BIG_DOUBLE tpw1r, tpw1i;
# else

  BIG_DOUBLE tpw1r, tpw1i, tpw2r, tpw2i, tpw5r, tpw5i;
#  endif
# endif

  BIG_DOUBLE *px;
  UL bj0, bj1, bj2, bj3, bj4, bj5, bj6, bj7;
  y_ptr pd0, pd1, pd2, pd3, pd4, pd5, pd6, pd7;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4, pad5, pad6, pad7;
  UL i, j, k, l, ll;
# ifdef Y_PRELOAD_DATA

  y_ptr ppd0r, ppd0i, ppd1r, ppd1i, ppd2r, ppd2i, ppd3r, ppd3i, ppd4r, ppd4i,
  ppd5r, ppd5i, ppd6r, ppd6i, ppd7r, ppd7i;
  y_ptr ppdl0r, ppdl0i, ppdl1r, ppdl1i, ppdl2r, ppdl2i, ppdl3r, ppdl3i,
  ppdl4r, ppdl4i, ppdl5r, ppdl5i, ppdl6r, ppdl6i, ppdl7r, ppdl7i;
# endif
# ifdef Y_KILL_BRANCHES

  long int nolast;
  UL ib[2];
#  ifndef USE_ASM

  UL issmall;
#   ifdef Y_VECTORIZE

  UL issmall1;
#   endif
#  endif

  aux[0] = highinv;
  aux[1] = lowinv;
  aux[2] = Gbig;
  aux[3] = Hsmall;
  aux[4] = bigA;
  aux[5] = 0.5;
  ib[0] = (UL)(-c);
  ib[1] = (UL)(b);
# endif
  /*DEBUG_VARS;*/

  ASSERT_ALIGNED_DOUBLE();

  if (Y_SBIT == 0)
    carry0 = -2.0;

  Err = 0.0;
  get_pads_8;
  init_bjs_8;

# ifdef Y_PRELOAD_DATA

  init_twiddle_factors_8_t;
# else

  init_twiddle_factors_8;
# endif

  px = Y_TWDB[Y_NRADICES-2];
# ifdef Y_PRELOAD_TRIGS

  trig_8_load_init (px);
# endif
# ifdef Y_PRELOAD_DATA

  cplx_radix_8_dif_dit_init (s, pd, ppdl);
# endif

  for (i = 0,j = 0; j < pad2; j += UPDATE, i++)
    {
      load_dwf_8 (i);
      for (k = 0; k < UPDATE; k += 2, px += Y_STEP)
        {
          l = (j + k) >> 1;
# ifdef Y_KILL_BRANCHES

          nolast = l - (pad-1);
          nolast >>= (BITS_PER_UL - 1);
# endif
# ifdef Y_PRELOAD_DATA

          ll = (addr(j + k + 2)) - (addr(j + k));
# else

          ll = addr(j + k);
# endif
          /*BEGIN_TIME;*/
          radix_8_twd_last_dit ( ll );
          /*END_TIME;*/
# ifdef SUM_CHECK

          sum_8 (SumOut);
# endif

          if(err_flag)
            {
              /*BEGIN_TIME;*/
# if defined(Y_ITANIUM)
              cplx_carry_norm_check_vector_4_ia64 (t, r, 0, 1, 2, 3);
              cplx_carry_norm_check_vector_4_ia64 (t, i, 0, 1, 2, 3);
              cplx_carry_norm_check_vector_4_ia64 (t, r, 4, 5, 6, 7);
              if(l == (pad - 1))
                bj7 = 0;
              cplx_carry_norm_check_vector_4_ia64 (t, i, 4, 5, 6, 7);
# else
#  if defined(Y_VECTORIZE)

              cplx_carry_norm_check_vector (t, 0, 1);
              cplx_carry_norm_check_vector (t, 2, 3);
              cplx_carry_norm_check_vector (t, 4, 5);
#  else

              cplx_carry_norm_check (t, 0);
              cplx_carry_norm_check (t, 1);
              cplx_carry_norm_check (t, 2);
              cplx_carry_norm_check (t, 3);
              cplx_carry_norm_check (t, 4);
              cplx_carry_norm_check (t, 5);
#  endif

              cplx_carry_norm_check (t, 6);
              cplx_carry_norm_last_check (t, 7);
# endif
              /*END_TIME;*/
            }
          else
            {
              /*BEGIN_TIME;*/
# ifdef Y_ITANIUM
              cplx_carry_norm_vector_8_ia64 (t, r, 0, 1, 2, 3, 4, 5, 6, 7);
              if(l == (pad - 1))
                bj7 = 0;
              cplx_carry_norm_vector_8_ia64 (t, i, 0, 1, 2, 3, 4, 5, 6, 7);
# else
#  ifdef Y_VECTORIZE

              cplx_carry_norm_vector (t, 0, 1);
              cplx_carry_norm_vector (t, 2, 3);
              cplx_carry_norm_vector (t, 4, 5);
#  else

              cplx_carry_norm (t, 0);
              cplx_carry_norm (t, 1);
              cplx_carry_norm (t, 2);
              cplx_carry_norm (t, 3);
              cplx_carry_norm (t, 4);
              cplx_carry_norm (t, 5);
#  endif

              cplx_carry_norm (t, 6);
              cplx_carry_norm_last (t, 7);
# endif
              /*END_TIME;*/
            }

          radix_8_first_dif (ll);


# ifdef SUM_CHECK

          SumIn += (t0r + t0i);
# endif

# ifdef Y_PRELOAD_TRIGS
#  ifdef Y_PRELOAD_DATA

          get_twiddle_factors_8_preload_consumed;
#  else

          get_twiddle_factors_8_preload;
#  endif
# else

          get_twiddle_factors_8;
# endif

        }
    }
  /* adjust the last carries */
  t0r = carry7;
  t1r = carry0;
  t2r = carry1;
  t3r = carry2;
  t4r = carry3;
  t5r = carry4;
  t6r = carry5;
  t7r = carry6;

  init_bjs_8;

  j = 0;
  load_two_to_phi_8 (j);
  cplx_carry (t, 0);
  cplx_carry (t, 1);
  cplx_carry (t, 2);
  cplx_carry (t, 3);
  cplx_carry (t, 4);
  cplx_carry (t, 5);
  cplx_carry (t, 6);
  cplx_carry (t, 7);

  /* add carry signal */
  radix_8_first_dif_add (x);

# ifdef SUM_CHECK

  SumIn += (t0r + t0i);
# endif

  if(err_flag)
    Err = maxerr;

  if(Y_SBIT)
    substract_two_8 (x, N);

#if defined(SUM_CHECK)
  /*printf("%lf\n", (2 * SumOut / N) - SumErr * SumErr);*/
#endif
}

#else

/***********************************************************************/
/* These are the Multithread adds needed                               */
/***********************************************************************/
/*
To work with mutithread environment, data is cut into independent slims.
As example, assume we are working with N=256 complex data size and plan
is [4,....]. Also assume we want four threads.

The access patern is
c[0], c[64], c[128], c[192] work with it and store new values in same
place
c[1], c[65], c[129], c[193] in carry and normalization we work with them
as r[]
....................
c[63], c[127], c[191], c[255]


Then there are 64 groups of 8 complex data. Divide them in four slims and
assign a slim per thread

Thread 0 :
c[0],c[64]...
c[1],c[65]...
..........
c[15],c[79]...

Thread 1 :
c[16],c[80]...
c[17],c[81]...
...........
c[31],c[47]

Thread 2 :
c[32],c[96]...
..........

Thread 3:
c[48],c[112]...
..........
c[63],c[127]....

It will require, of course, a last carry adjust among threads.

_ofs is the offset, in complex size,  of every thread nthr
_ofs = _n0 * nthr;  0 <= nthr < max_num_of_threads
where _n0 is the basic length of a thread, and MUST be a
power of two
*/
# define get_pads_8_mp(_ofs)         \
  pad2 = (pad << 1);                 \
  pad4 = (pad << 2);                 \
  pad3 = pad + pad2;                 \
  pad5 = pad + pad4;                 \
  pad6 = pad3 + pad3;                \
  pad7 = pad4 + pad3;                \
  pd0 = x + addr((_ofs << 1));       \
  prefetch_data( pd0, 0);            \
  pd1 = x + addr((pad  + _ofs) << 1);\
  prefetch_data( pd1, 0);            \
  pd2 = x + addr((pad2 + _ofs) << 1);\
  prefetch_data( pd2, 0);            \
  pd3 = x + addr((pad3 + _ofs) << 1);\
  prefetch_data( pd3, 0);            \
  pd4 = x + addr((pad4 + _ofs) << 1);\
  prefetch_data( pd4, 0);            \
  pd5 = x + addr((pad5 + _ofs) << 1);\
  prefetch_data( pd5, 0);            \
  pd6 = x + addr((pad6 + _ofs) << 1);\
  prefetch_data( pd6, 0);            \
  pd7 = x + addr((pad7 + _ofs) << 1);\
  prefetch_data( pd7, 0);

/*
To compute the initial bj of every thread
It will be computed at first time. so we avoid to compute it every
iteration
*/

# define init_bjs_8_mp \
  bj0 = bp[0];\
  bj1 = bp[1];\
  bj2 = bp[2];\
  bj3 = bp[3];\
  bj4 = bp[4];\
  bj5 = bp[5];\
  bj6 = bp[6];\
  bj7 = bp[7];

/*
The twidle factors now are not (1.0,0.0) at the begining
so it should be computed

inc=(Y_POWERS[Y_PLAN[0]-1])<<1;
px = Y_TWDB[Y_NRADICES - 2] + _ofs * inc - 2;

It could be computed at a first iteration, so we save
some work
*/

# if defined(Y_MINIMUM)
#  define init_twiddle_factors_8_mp(_px)  \
  tw1r = _px[0]; tw1i = _px[1];           \
  cplx_trig_min_square (tw2, tw1);        \
  cplx_trig_min_square (tw4, tw2);        \
  cplx_divmul (tw3, tw5, tw4, tw1);       \
  cplx_trig_min_square (tw6, tw3);        \
  cplx_trig_min_mul (tw7, tw4, tw3);


#  define init_twiddle_factors_8_t_mp(_px) \
  t1r = _px[0]; t1i = _px[1];              \
  cplx_trig_min_square (t2, t1);           \
  cplx_trig_min_square (t4, t2);           \
  cplx_divmul (t3, t5, t4, t1);            \
  cplx_trig_min_square (t6, t3);           \
  cplx_trig_min_mul (t7, t4, t3);

# elif defined(Y_MAXIMUM)

#  define init_twiddle_factors_8_mp(_px)                   \
         tw1r = *(_px);      tw1i = *(_px + 1);            \
         tw2r = *(_px + 2);  tw2i = *(_px + 3);            \
         tw3r = *(_px + 4);  tw3i = *(_px + 5);            \
         tw4r = *(_px + 6);  tw4i = *(_px + 7);            \
         tw5r = *(_px + 8);  tw5i = *(_px + 9);            \
         tw6r = *(_px + 10); tw6i = *(_px + 11);           \
         tw7r = *(_px + 12); tw7i = *(_px + 13);           \
         prefetch_data_trig (_px, Y_STEP);                 \
         prefetch_data_trig (_px, Y_STEP + Y_CACHE_LINE);  \
         prefetch_data_trig (_px, Y_STEP + 2 *Y_CACHE_LINE);

#  define init_twiddle_factors_8_t_mp(_px)                  \
         t1r = *(_px);      t1i = *(_px + 1);               \
         t2r = *(_px + 2);  t2i = *(_px + 3);               \
         t3r = *(_px + 4);  t3i = *(_px + 5);               \
         t4r = *(_px + 6);  t4i = *(_px + 7);               \
         t5r = *(_px + 8);  t5i = *(_px + 9);               \
         t6r = *(_px + 10); t6i = *(_px + 11);              \
         t7r = *(_px + 12); t7i = *(_px + 13);              \
         prefetch_data_trig (_px, Y_STEP);                  \
         prefetch_data_trig (_px, Y_STEP + Y_CACHE_LINE);   \
         prefetch_data_trig (_px, Y_STEP + 2 * Y_CACHE_LINE);


# else
#  define init_twiddle_factors_8_mp(_px)                 \
         tw1r = *(_px);                                  \
         tw1i = *(_px + 1);                              \
         tw2r = *(_px + 2);                              \
         tw2i = *(_px + 3);                              \
         tw5r = *(_px + 4);                              \
         tw5i = *(_px + 5);                              \
         cplx_divmul (tw3, tw7, tw5, tw2);               \
	 cplx_divmul (tw4, tw6, tw5, tw1);               \
         prefetch_data_trig (_px, Y_STEP);               \
         prefetch_data_trig (_px, Y_STEP + Y_CACHE_LINE);

#  define init_twiddle_factors_8_t_mp(_px)         \
         t1r = *(_px);                             \
         t1i = *(_px + 1);                         \
         t2r = *(_px + 2);                         \
         t2i = *(_px + 3);                         \
         t5r = *(_px + 4);                         \
         t5i = *(_px + 5);                         \
         cplx_divmul (t3, t7, t5, t2);             \
	 cplx_divmul (t4, t6, t5, t1);             \
         prefetch_data_trig (_px, 6);              \
         prefetch_data_trig (_px, 6 + Y_CACHE_LINE);

# endif

# define save_carries_8_mp(_pcarr)   \
        _pcarr[0] = carry0;          \
        _pcarr[1] = carry1;          \
        _pcarr[2] = carry2;          \
        _pcarr[3] = carry3;          \
        _pcarr[4] = carry4;          \
        _pcarr[5] = carry5;          \
        _pcarr[6] = carry6;          \
        _pcarr[7] = carry7;


# define radix_8_first_dif_add_mp          \
	cplx_addsub (t0, t4);               \
	cplx_addsub (t2, t6);               \
	cplx_addsub (t0, t2);               \
	cplx_mul_1_4_F_addsub (t4, t6);     \
	cplx_addsub (t1, t5);               \
	cplx_addsub (t3, t7);               \
	cplx_addsub (t1, t3);               \
	cplx_mul_1_4_F_addsub (t5, t7);     \
	cplx_addsub (t0, t1);               \
	cplx_local_to_data_add (pd0, 0, t0);\
	cplx_local_to_data_add (pd4, 0, t1);\
	cplx_mul_1_8_F_addsub (t4, t5);     \
	cplx_local_to_data_add (pd1, 0, t4);\
	cplx_local_to_data_add (pd5, 0, t5);\
	cplx_mul_1_4_F_addsub (t2, t3);     \
	cplx_local_to_data_add (pd2, 0, t2);\
	cplx_local_to_data_add (pd6, 0, t3);\
	cplx_mul_3_8_F_addsub (t6, t7);     \
	cplx_local_to_data_add (pd3, 0, t6);\
	cplx_local_to_data_add (pd7, 0, t7);



BIG_DOUBLE dit_carry_norm_dif_8_mp(BIG_DOUBLE *x,BIG_DOUBLE *ptx, BIG_DOUBLE *pcarr, UL *bp, UL N ,UL n0, UL nthr,UL err_flag)
{
# ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6], maxerr = 0.0;
  double Y_ALIGNED(16) tt0[2], tt1[2], tt2[2], tt3[2], tt4[2], tt5[2], tt6[2], tt7[2];
# else

BIG_DOUBLE hiinv = highinv, loinv = lowinv;
BIG_DOUBLE maxerr = 0.0, ttmpSmall = Gsmall, ttmpBig = Gbig,
                                     ttpSmall = Hsmall;
BIG_DOUBLE ttmp0, ttmp1, ttmp2, ttmp3, ttmp4, ttmp5, ttmp6, ttmp7;
BIG_DOUBLE ttp0, ttp1, ttp2, ttp3, ttp4, ttp5, ttp6, ttp7;
# endif
# if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A = bigA, B = bigB;
# endif

  BIG_DOUBLE carry0 = 0.0, carry1 = 0.0, carry2 = 0.0, carry3 = 0.0,
                                    carry4 = 0.0, carry5 = 0.0, carry6 = 0.0, carry7 = 0.0;
# ifndef Y_PRELOAD_DATA

  BIG_DOUBLE tw1r, tw1i, tw2r, tw2i, tw3r, tw3i, tw4r, tw4i, tw5r, tw5i,
  tw6r, tw6i, tw7r, tw7i;
# endif

  BIG_DOUBLE t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i,
  t6r, t6i, t7r, t7i;
# ifdef Y_ITANIUM

  y_limb_t s0r, s0i, s1r, s1i, s2r, s2i, s3r, s3i, s4r, s4i, s5r, s5i,
  s6r, s6i, s7r, s7i;
  y_limb_t fr = F_1_8r,fi = F_1_8i;
# endif
# ifdef Y_PRELOAD_TRIGS
# if defined(Y_MAXIMUM)

  BIG_DOUBLE tpw1r, tpw1i, tpw2r, tpw2i, tpw3r, tpw3i, tpw4r, tpw4i,
  tpw5r, tpw5i, tpw6r, tpw6i, tpw7r, tpw7i;
# elif defined(Y_MINIMUM)

BIG_DOUBLE tpw1r, tpw1i;
# else

BIG_DOUBLE tpw1r, tpw1i, tpw2r, tpw2i, tpw5r, tpw5i;
#  endif
# endif

  BIG_DOUBLE *px;
  UL bj0, bj1, bj2, bj3, bj4, bj5, bj6, bj7;
  y_ptr pd0, pd1, pd2, pd3, pd4, pd5, pd6, pd7;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4, pad5, pad6, pad7, lim;
  UL i, j, k, l, ll, ofs;
# ifdef Y_PRELOAD_DATA

  y_ptr ppd0r, ppd0i, ppd1r, ppd1i, ppd2r, ppd2i, ppd3r, ppd3i, ppd4r, ppd4i,
  ppd5r, ppd5i, ppd6r, ppd6i, ppd7r, ppd7i;
  y_ptr ppdl0r, ppdl0i, ppdl1r, ppdl1i, ppdl2r, ppdl2i, ppdl3r, ppdl3i,
  ppdl4r, ppdl4i, ppdl5r, ppdl5i, ppdl6r, ppdl6i, ppdl7r, ppdl7i;
# endif
# ifdef Y_KILL_BRANCHES

  long int nolast;
  UL ib[2];
#  ifndef USE_ASM

  UL issmall;
#   ifdef Y_VECTORIZE

  UL issmall1;
#   endif
#  endif

  aux[0] = highinv;
  aux[1] = lowinv;
  aux[2] = Gbig;
  aux[3] = Hsmall;
  aux[4] = bigA;
  aux[5] = 0.5;
  ib[0] = (UL)(-c);
  ib[1] = (UL)(b);
# endif

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
      init_bjs_8_mp;
# ifdef Y_PRELOAD_DATA

      get_pads_8_mp (ofs);
      init_twiddle_factors_8_t_mp (ptx);
# else

get_pads_8;
init_twiddle_factors_8_mp (ptx);
# endif

      px = ptx + Y_STEP;
    }
  else
    {
      if(Y_SBIT == 0)
        carry0 = -2.0;
      get_pads_8;
      init_bjs_8;
# ifdef Y_PRELOAD_DATA

      init_twiddle_factors_8_t;
# else

init_twiddle_factors_8;
# endif

      px = Y_TWDB[Y_NRADICES - 2];
    }
# ifdef Y_PRELOAD_TRIGS
  trig_8_load_init (px);
# endif
# ifdef Y_PRELOAD_DATA

  cplx_radix_8_dif_dit_init (s, pd, ppdl);
# endif

  for (i = (ofs << 1) / UPDATE, j = (ofs << 1); j < lim; j += UPDATE, i++)
    {
      load_dwf_8 (i);

      for (k = 0; k < UPDATE; k += 2, px += Y_STEP)
        {
          l = (j + k) >> 1;

# ifdef Y_KILL_BRANCHES

          nolast = l - (pad-1);
          nolast >>= (BITS_PER_UL - 1);
# endif
# ifdef Y_PRELOAD_DATA

          ll = (addr(j + k + 2)) - (addr(j + k));
# else

ll = addr(j + k);
# endif

          radix_8_twd_last_dit (ll);

# ifdef SUM_CHECK

          sum_8 (Y_SUMOUT[nthr]);
# endif

          if(err_flag)
            {
# if defined(Y_ITANIUM)
              cplx_carry_norm_check_vector_4_ia64 (t, r, 0, 1, 2, 3);
              cplx_carry_norm_check_vector_4_ia64 (t, i, 0, 1, 2, 3);
              cplx_carry_norm_check_vector_4_ia64 (t, r, 4, 5, 6, 7);
              if(l == (pad - 1))
                bj7 = 0;
              cplx_carry_norm_check_vector_4_ia64 (t, i, 4, 5, 6 ,7);
# else
#  if defined(Y_VECTORIZE)

cplx_carry_norm_check_vector (t, 0, 1);
cplx_carry_norm_check_vector (t, 2, 3);
cplx_carry_norm_check_vector (t, 4, 5);
#  else

cplx_carry_norm_check (t, 0);
cplx_carry_norm_check (t, 1);
cplx_carry_norm_check (t, 2);
cplx_carry_norm_check (t, 3);
cplx_carry_norm_check (t, 4);
cplx_carry_norm_check (t, 5);
#  endif

cplx_carry_norm_check (t, 6);
cplx_carry_norm_last_check (t, 7);
# endif

            }
          else
            {
# ifdef Y_ITANIUM
              cplx_carry_norm_vector_8_ia64 (t, r, 0, 1, 2, 3, 4, 5, 6, 7);
              if(l == (pad - 1))
                bj7 = 0;
              cplx_carry_norm_vector_8_ia64 (t, i, 0, 1, 2, 3, 4, 5, 6, 7);
# else
#  ifdef Y_VECTORIZE

cplx_carry_norm_vector (t, 0, 1);
cplx_carry_norm_vector (t, 2, 3);
cplx_carry_norm_vector (t, 4, 5);
#  else

cplx_carry_norm (t, 0);
cplx_carry_norm (t, 1);
cplx_carry_norm (t, 2);
cplx_carry_norm (t, 3);
cplx_carry_norm (t, 4);
cplx_carry_norm (t, 5);
#  endif

cplx_carry_norm (t, 6);
cplx_carry_norm_last (t, 7);
# endif

            }
          radix_8_first_dif (ll);

# if defined(SUM_CHECK)

          Y_SUMIN[nthr] += (t0r + t0i);
# endif

# ifdef Y_PRELOAD_TRIGS
#  ifdef Y_PRELOAD_DATA

          get_twiddle_factors_8_preload_consumed;
#  else

get_twiddle_factors_8_preload;
#  endif
# else

get_twiddle_factors_8;
# endif

        }
    }
  /* save the carries */
  save_carries_8_mp (pcarr);

  /* return the error */
  return maxerr;
}


void dit_carry_norm_dif_8_last_carries_mp(BIG_DOUBLE *x, UL *bp, UL N ,UL n0, int nthr)
{
# ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6];
  double Y_ALIGNED(16) tt0[2], tt1[2], tt2[2], tt3[2], tt4[2], tt5[2], tt6[2], tt7[2];
# else

BIG_DOUBLE hiinv = highinv, loinv = lowinv;
BIG_DOUBLE ttmpSmall = Gsmall, ttmpBig = Gbig,
                       ttpSmall = Hsmall;
BIG_DOUBLE ttmp0, ttmp1, ttmp2, ttmp3, ttmp4, ttmp5, ttmp6, ttmp7;
BIG_DOUBLE ttp0, ttp1, ttp2, ttp3, ttp4, ttp5, ttp6, ttp7;
# endif
# if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A = bigA, B = bigB;
# endif

  BIG_DOUBLE carry0, carry1, carry2, carry3, carry4, carry5, carry6, carry7;
  BIG_DOUBLE t0r, t0i, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i, t5r, t5i,
  t6r, t6i, t7r, t7i;
  UL bj0, bj1, bj2, bj3, bj4, bj5, bj6, bj7;
  y_ptr pd0, pd1, pd2, pd3, pd4, pd5, pd6, pd7;
  UL pad = Y_LRIGHT[1], pad2, pad3, pad4, pad5, pad6, pad7;
  UL j, ofs;
# ifdef Y_KILL_BRANCHES

  UL ib[2];
#  ifndef USE_ASM

  UL issmall;
#  endif

  aux[0] = highinv;
  aux[1] = lowinv;
  aux[2] = Gbig;
  aux[3] = Hsmall;
  aux[4] = bigA;
  aux[5] = 0.5;
  ib[0] = (UL)(-c);
  ib[1] = (UL)(b);
# endif

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
    }
  else
    {
      t0r = Y_CARRIES[Y_NTHREADS-1][7];
      t1r = Y_CARRIES[Y_NTHREADS-1][0];
      t2r = Y_CARRIES[Y_NTHREADS-1][1];
      t3r = Y_CARRIES[Y_NTHREADS-1][2];
      t4r = Y_CARRIES[Y_NTHREADS-1][3];
      t5r = Y_CARRIES[Y_NTHREADS-1][4];
      t6r = Y_CARRIES[Y_NTHREADS-1][5];
      t7r = Y_CARRIES[Y_NTHREADS-1][6];
    }

  init_bjs_8_mp;
  ofs = (n0 * nthr);

  get_pads_8;
  j = (ofs << 1) / UPDATE;

  load_two_to_phi_8 (j);
  load_two_to_minusphi_8 (j);

  cplx_carry (t, 0);
  cplx_carry (t, 1);
  cplx_carry (t, 2);
  cplx_carry (t, 3);
  cplx_carry (t, 4);
  cplx_carry (t, 5);
  cplx_carry (t, 6);
  cplx_carry (t, 7);

  /* add carry signal */

  get_pads_8_mp (ofs)

  radix_8_first_dif_add_mp;

# if defined(SUM_CHECK)

  Y_SUMIN[nthr] += (t0r + t0i);
# endif

}


void dit_carry_norm_dif_8(BIG_DOUBLE *x, UL N ,UL err_flag)
{
  int i, n0;
  n0 = Y_LRIGHT[1] / Y_NTHREADS;

  /* make n0 divisble by Y_UPDATE */
  n0 = ((n0 >> (Y_SHIFT_UPDATE - 1)) << (Y_SHIFT_UPDATE - 1));

  /* The multithread loop */
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop storeback(i)
#  pragma MP taskloop shared(x, Y_ERRS, Y_TWX, Y_CARRIES, Y_BJS, N, n0, err_flag)
# endif

  for(i = 0; i < Y_NTHREADS; i++)
    {
      Y_ERRS[i] = dit_carry_norm_dif_8_mp (x, Y_TWX[i], Y_CARRIES[i], Y_BJS[i],
                                           N, n0, i, err_flag);
    }

  /* maxerr */
  Err = 0.0;
  if(err_flag)
    {
      for(i = 0; i < Y_NTHREADS; i++)
        if(Y_ERRS[i] > Err)
          Err = Y_ERRS[i];
    }

  /* And finally, the last carry propagation */
# if defined(_OPENMP)
#  pragma omp parallel for default(shared)
# elif defined(_SUNMP)
#  pragma MP taskloop maxcpus(Y_NUM_THREADS)
#  pragma MP taskloop shared(x, Y_BJS, N, n0)
# endif
  for(i = 0; i < Y_NTHREADS; i++)
    {
      dit_carry_norm_dif_8_last_carries_mp (x, Y_BJS[i],
                                            N, n0, i);
    }

  /* Substract two if shifts*/
  if(Y_SBIT)
    substract_two_8 (x, N);

#if defined(SUM_CHECK)

  for (i = 0; i < Y_NTHREADS; i++)
    {
      SumOut += Y_SUMOUT[i];
      SumIn += Y_SUMIN[i];
    }
#endif

}

#endif
/*$Id$*/










