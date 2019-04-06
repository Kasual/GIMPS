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
   It includes routines for a 16_last_radix reduction
   It is a huge routine
*/
#include <stdio.h>
#include <stdlib.h>

#include "yeafft.h"
#include "mccomp.h"
#include "ydebug.h"
#define NDEBUG1

#if defined(Y_MINIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 2
# include "yminimum.h"

# define get_twiddle_factors_16_last(_j) \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          cplx_trig_min_load(tw1,px);\
          cplx_trig_min_square(tw2,tw1);\
          cplx_trig_min_square(tw4,tw2);\
          cplx_trig_min_square(tw8,tw4);\
          cplx_divmul(tw3,tw5,tw4,tw1);\
          cplx_divmul(tw7,tw9,tw8,tw1);\
          cplx_divmul(tw6,tw10,tw8,tw2);\
          cplx_trig_min_square(tw12,tw6);\
          cplx_divmul(tw11,tw13,tw12,tw1);\
          cplx_trig_min_square(tw14,tw7);\
          cplx_trig_min_mul(tw15,tw7,tw8);


#elif defined(Y_MAXIMUM)
# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 30

# define get_twiddle_factors_16_last(_j) \
     px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
     tw1r = px[ 0]; tw1i = px[ 1];             \
     tw2r = px[ 2]; tw2i = px[ 3];             \
     tw3r = px[ 4]; tw3i = px[ 5];             \
     tw4r = px[ 6]; tw4i = px[ 7];             \
     tw5r = px[ 8]; tw5i = px[ 9];             \
     tw6r = px[10]; tw6i = px[11];             \
     tw7r = px[12]; tw7i = px[13];             \
     tw8r = px[14]; tw8i = px[15];             \
     tw9r = px[16]; tw9i = px[17];             \
     tw10r= px[18]; tw10i= px[19];             \
     tw11r= px[20]; tw11i= px[21];             \
     tw12r= px[22]; tw12i= px[23];             \
     tw13r= px[24]; tw13i= px[25];             \
     tw14r= px[26]; tw14i= px[27];             \
     tw15r= px[28]; tw15i= px[29];             \
     prefetch_data_trig_nta(px, Y_STEP + Y_CACHE_LINE);\
     prefetch_data_trig_nta(px, Y_STEP + 2*Y_CACHE_LINE);\
     prefetch_data_trig_nta(px, Y_STEP + 3*Y_CACHE_LINE);\
     prefetch_data_trig_nta(px, Y_STEP + 4*Y_CACHE_LINE);\

#else

# define get_twiddle_factors_16_last(_j) \
          px=&(Y_TWDF[Y_NRADICES-2][(_j * inc)<<1]);\
          tw1r= px[0];\
          tw1i= px[1];\
          tw4r= px[4];\
          tw4i= px[5];\
	  cplx_divmul(tw3,tw5,tw4,tw1);\
          tw2r= px[2];\
          tw2i= px[3];\
	  tw8r= px[6];\
          tw8i= px[7];\
	  cplx_divmul(tw6,tw10,tw8,tw2);\
	  cplx_divmul(tw7,tw9,tw8,tw1);\
	  tw13r=px[8];\
          tw13i=px[9];\
	  cplx_divmul(tw11,tw15,tw13,tw2);\
	  cplx_divmul(tw12,tw14,tw13,tw1);



# ifdef Y_STEP
#  undef Y_STEP
# endif
# define Y_STEP 10
#endif


#define radix_16_twd_last(_z,_pd,_d,_j)\
          _pd = _d + addr((_j)<<1);\
	  cplx_load_muladdsub_pp(_z##0 , _z##1 , tw8, 0, _pd, _pd+16);\
	  cplx_load_mulmuladdsub_pp( _z##2 , _z##3 ,tw4,tw12,0,_pd+8,_pd+24);\
	  cplx_addsub( _z##0 , _z##2 );\
	  cplx_mul_1_4_F_addsub( _z##1 , _z##3 );\
	  cplx_load_mulmuladdsub_pp( _z##4 , _z##5 ,tw2,tw10,0,_pd+4,_pd+20);\
	  cplx_load_mulmuladdsub_pp( _z##6 , _z##7 ,tw6,tw14,0,_pd+12,_pd+28);\
	  cplx_addsub( _z##4 , _z##6 );\
	  cplx_mul_1_4_F_addsub( _z##5 , _z##7 );\
	  cplx_addsub( _z##0 , _z##4);\
	  cplx_mul_1_8_F_addsub( _z##1 , _z##5 );\
	  cplx_mul_1_4_F_addsub( _z##2 , _z##6 );\
	  cplx_mul_3_8_F_addsub( _z##3 , _z##7 );\
	  cplx_load_mulmuladdsub_pp( _z##8 , _z##9 ,tw1,tw9,0,_pd+2,_pd+18);\
	  cplx_load_mulmuladdsub_pp( _z##10 , _z##11 ,tw5,tw13,0,_pd+10,_pd+26);\
	  cplx_addsub( _z##8 , _z##10 );\
	  cplx_mul_1_4_F_addsub( _z##9 , _z##11 );\
	  cplx_load_mulmuladdsub_pp( _z##12 , _z##13 ,tw3,tw11,0,_pd+6,_pd+22);\
	  cplx_load_mulmuladdsub_pp( _z##14 , _z##15 ,tw7,tw15,0,_pd+14,_pd+30);\
          prefetch_p( _pd + 32 - Y_CACHE_LINE);\
          prefetch_p( _pd + 32);\
          prefetch_p( _pd + 32 + Y_CACHE_LINE);\
          prefetch_p( _pd + 32 + 2*Y_CACHE_LINE);\
	  cplx_addsub( _z##12 , _z##14 );\
	  cplx_mul_1_4_F_addsub( _z##13 , _z##15 );\
	  cplx_addsub( _z##8 , _z##12 );\
	  cplx_mul_1_8_F_addsub( _z##9 , _z##13 );\
	  cplx_mul_1_4_F_addsub( _z##10 , _z##14 );\
	  cplx_mul_3_8_F_addsub( _z##11 , _z##15 );\
	  \
	  cplx_addsub( _z##0 , _z##8 );\
	  cplx_muladdsub( _z##1 , _z##9 , F_1_16);\
	  cplx_mul_1_8_F_addsub( _z##2 , _z##10 );\
	  cplx_muladdsub( _z##3 , _z##11 , F_3_16);\
	  cplx_mul_1_4_F_addsub( _z##4 , _z##12 );\
	  cplx_muladdsub( _z##5 , _z##13 , F_5_16);\
	  cplx_mul_3_8_F_addsub( _z##6 , _z##14);\
	  cplx_muladdsub( _z##7 , _z##15 , F_7_16);\



#define radix_16_notwd_first(_pd,_z,_j)\
	  cplx_addsub( _z##0 , _z##8 );\
	  cplx_addsub( _z##4 , _z##12 );\
	  cplx_addsub( _z##0 , _z##4 );\
          cplx_mul_1_4_B_addsub( _z##8 , _z##12 );\
	  cplx_addsub( _z##2 , _z##10 );\
	  cplx_addsub( _z##6 , _z##14 );\
	  cplx_addsub( _z##2 , _z##6 );\
	  cplx_mul_1_4_B_addsub( _z##10 , _z##14 );\
	  cplx_addsub( _z##0 , _z##2 );\
	  cplx_mul_1_8_B_addsub( _z##8 , _z##10 );\
	  cplx_mul_1_4_B_addsub( _z##4 , _z##6 );\
	  cplx_mul_3_8_B_addsub( _z##12 , _z##14 );\
	  cplx_addsub( _z##1 , _z##9 );\
	  cplx_addsub( _z##5 , _z##13 );\
	  cplx_addsub( _z##1 , _z##5 );\
	  cplx_mul_1_4_B_addsub( _z##9 , _z##13 );\
	  cplx_addsub( _z##3 , _z##11 );\
	  cplx_addsub( _z##7 , _z##15 );\
	  cplx_addsub( _z##3 , _z##7 );\
	  cplx_mul_1_4_B_addsub( _z##11 , _z##15 );\
	  cplx_addsub( _z##1 , _z##3 );\
	  cplx_mul_1_8_B_addsub( _z##9 , _z##11 );\
	  cplx_mul_1_4_B_addsub( _z##5 , _z##7 );\
	  cplx_mul_3_8_B_addsub( _z##13 , _z##15 );\
	  \
	  cplx_addsub_store_p(0, _pd, _pd+16, _z##0 , _z##1 );\
	  cplx_muladdsub_store_p(0,_pd+2,_pd+18,_z##8 , _z##9 , B_1_16 );\
	  cplx_mul_1_8_B_addsub_store_p(0,_pd+4,_pd+20, _z##4, _z##5 );\
	  cplx_muladdsub_store_p(0,_pd+6,_pd+22,_z##12 , _z##13 , B_3_16 );\
	  cplx_mul_1_4_B_addsub_store_p(0,_pd+8,_pd+24, _z##2 , _z##3 );\
	  cplx_muladdsub_store_p(0,_pd+10,_pd+26,_z##10 , _z##11 , B_5_16 );\
	  cplx_mul_3_8_B_addsub_store_p(0,_pd+12,_pd+28, _z##6 , _z##7 );\
	  cplx_muladdsub_store_p(0,_pd+14,_pd+30, _z##14 , _z##15 , B_7_16);

#if (Y_AVAL > 4) || ((Y_AVAL == 4) && defined(Y_MANY_REGISTERS))

void radix_16_dif_square_dit(y_ptr d, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,t8r,
  t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i,u8r,
  u8i,u9r,u9i,u10r,u10i,u11r,u11i,u12r,u12i,u13r,u13i,u14r,u14i,u15r,u15i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_limb_t wkr,wki;
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[15];
  int nr;
  y_ptr px,pt,pu;
  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();


  /*if(Y_PLAN[Y_NRADICES-1] != 16) return;*/
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
  get_twiddle_factors_16_last(j);
  radix_16_twd_last(t,pt,d,j);
  /* dyadic mul for k=0 */
  square_nested_eq( t0 , 1.0);
  square_nested_1_4( t4 , t12 , 1.0, 0.0 );
  square_nested_eq( t8 , -1.0);

  square_nested( t1 , t15, F_1_16r, F_1_16i );
  square_nested_1_4( t5 , t11 , F_1_16r, F_1_16i );

  square_nested( t2 , t14 , F_1_8r, F_1_8i );
  square_nested_1_4( t6 , t10 , F_1_8r, F_1_8i );

  square_nested( t3 , t13 , F_3_16r, F_3_16i);
  square_nested_1_4( t7 , t9 ,  F_3_16r, F_3_16i);

  radix_16_notwd_first(pt,t,j);

  if(nc==0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>4;
      bj=(Y_LRIGHT[nr]>>4)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_16_last(jj);
          l=jj<<4;
          radix_16_twd_last(u,pu,d,l);
          get_twiddle_factors_16_last(ii);
          l=ii<<4;
          radix_16_twd_last(t,pt,d,l);

          square_nested(t0, u15, tw1r, tw1i);
          square_nested_1_4(t4, u11, tw1r, tw1i);
          square_nested_1_2(t8, u7, tw1r, tw1i);
          square_nested_3_4(t12, u3, tw1r, tw1i);

          cplx_mul(wk,tw1,F_1_16);
          square_nested(t1, u14, wkr, wki);
          square_nested_1_4(t5, u10, wkr, wki);
          square_nested_1_2(t9, u6, wkr, wki);
          square_nested_3_4(t13, u2, wkr, wki);

          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
          square_nested(t2, u13, wkr, wki);
          square_nested_1_4(t6, u9, wkr, wki);
          square_nested_1_2(t10, u5, wkr, wki);
          square_nested_3_4(t14, u1, wkr, wki);

          cplx_mul(wk,tw1,F_3_16);
          square_nested(t3, u12, wkr, wki);
          square_nested_1_4(t7, u8, wkr, wki);
          square_nested_1_2(t11, u4, wkr, wki);
          square_nested_3_4(t15, u0, wkr, wki);

          radix_16_notwd_first(pt,t,l);
          l=jj<<4;
          radix_16_notwd_first(pu,u,l);
        }
      if(ii==jj)
        {
          get_twiddle_factors_16_last(ii);
          l=ii<<4;
          radix_16_twd_last(t,pt,d,l);

          square_nested(t0, t15, tw1r,tw1i);
          square_nested_1_4(t4, t11, tw1r, tw1i);

          cplx_mul(wk,tw1,F_1_16);
          square_nested(t1, t14, wkr,wki);
          square_nested_1_4(t5, t10, wkr, wki);

          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
          square_nested(t2, t13, wkr, wki);
          square_nested_1_4(t6, t9, wkr, wki);

          cplx_mul(wk,tw1,F_3_16);
          square_nested(t3, t12, wkr,wki);
          square_nested_1_4(t7, t8, wkr, wki);

          radix_16_notwd_first(pt,t,l);
        }
    }
}


void radix_16_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,t8r,
  t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i,u8r,
  u8i,u9r,u9i,u10r,u10i,u11r,u11i,u12r,u12i,u13r,u13i,u14r,u14i,u15r,u15i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i,s8r,
  s8i,s9r,s9i,s10r,s10i,s11r,s11i,s12r,s12i,s13r,s13i,s14r,s14i,s15r,s15i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i,v8r,
  v8i,v9r,v9i,v10r,v10i,v11r,v11i,v12r,v12i,v13r,v13i,v14r,v14i,v15r,v15i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_limb_t wkr,wki;
  y_size_t j,bi,bj,jj,ii,l,inc=Y_POWERS[15];
  int nr;
  y_ptr px,pt,pu,ps,pv;

  ASSERT_ALIGNED_DOUBLE();

  /*if(Y_PLAN[Y_NRADICES-1] != 16) return;*/

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
  get_twiddle_factors_16_last(j);
  radix_16_twd_last(s,ps,d1,j);
  radix_16_twd_last(t,pt,d2,j);
  /* dyadic mul for k=0 */
  conv_nested( s0, s0, t0, t0, 1.0, 0.0 );
  conv_nested_1_4( s4, s12, t4, t12, 1.0, 0.0 );
  conv_nested_1_2( s8, s8, t8, t8, 1.0, 0.0 );

  conv_nested( s1, s15, t1, t15, F_1_16r, F_1_16i );
  conv_nested_1_4( s5, s11, t5, t11, F_1_16r, F_1_16i );

  conv_nested( s2, s14, t2, t14, F_1_8r, F_1_8i );
  conv_nested_1_4( s6, s10, t6, t10, F_1_8r, F_1_8i );

  conv_nested( s3, s13, t3, t13, F_3_16r, F_3_16i );
  conv_nested_1_4( s7, s9, t7, t9, F_3_16r, F_3_16i );

  radix_16_notwd_first(ps,s,j);
  if(nc == 0)
    return;

  for (nr=Y_NRADICES-2;nr>((int)(Y_NRADICES)-(int)(nc)-2);nr--)
    {
      bi=(Y_LRIGHT[nr+1])>>4;
      bj=(Y_LRIGHT[nr]>>4)-1;
      for(ii=bi,jj=bj;ii<jj;ii++,jj--)
        {
          get_twiddle_factors_16_last(jj);
          l=jj<<4;
          radix_16_twd_last(u,pu,d1,l);
          radix_16_twd_last(v,pv,d2,l);
          get_twiddle_factors_16_last(ii);
          l=ii<<4;
          radix_16_twd_last(s,ps,d1,l);
          radix_16_twd_last(t,pt,d2,l);

          conv_nested( s0, u15, t0, v15, tw1r, tw1i);
          conv_nested_1_4( s4, u11, t4, v11, tw1r, tw1i);
          conv_nested_1_2( s8, u7, t8, v7, tw1r, tw1i);
          conv_nested_3_4( s12, u3, t12, v3, tw1r, tw1i);

          cplx_mul(wk,tw1,F_1_16);
          conv_nested(s1, u14, t1, v14, wkr, wki);
          conv_nested_1_4(s5, u10, t5, v10, wkr, wki);
          conv_nested_1_2(s9, u6, t9, v6, wkr, wki);
          conv_nested_3_4(s13, u2, t13, v2, wkr, wki);

          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
          conv_nested(s2, u13, t2, v13, wkr, wki);
          conv_nested_1_4(s6, u9, t6, v9, wkr, wki);
          conv_nested_1_2(s10, u5, t10, v5, wkr, wki);
          conv_nested_3_4(s14, u1, t14, v1, wkr, wki);

          cplx_mul(wk,tw1,F_3_16);
          conv_nested(s3, u12, t3, v12, wkr, wki);
          conv_nested_1_4(s7, u8, t7, v8, wkr, wki);
          conv_nested_1_2(s11, u4, t11, v4, wkr, wki);
          conv_nested_3_4(s15, u0, t15, v0, wkr, wki);

          radix_16_notwd_first(ps,s,l);
          l=jj<<4;
          radix_16_notwd_first(pu,u,l);
        }
      if(ii==jj)
        {
          get_twiddle_factors_16_last(ii);
          l=ii<<4;
          radix_16_twd_last(s,ps,d1,l);
          radix_16_twd_last(t,pt,d2,l);

          conv_nested(s0, s15, t0, t15, tw1r, tw1i);
          conv_nested_1_4(s4, s11, t4, t11, tw1r, tw1i);

          cplx_mul(wk,tw1,F_1_16);
          conv_nested(s1, s14, t1, t14, wkr, wki);
          conv_nested_1_4(s5, s10, t5, t10, wkr, wki);

          wkr=(tw1r+tw1i)*F_1_8r;
          wki=(tw1i-tw1r)*F_1_8r;
          conv_nested(s2, s13, t2, t13, wkr, wki);
          conv_nested_1_4(s6, s9, t6, t9, wkr, wki);

          cplx_mul(wk,tw1,F_3_16);
          conv_nested(s3, s12, t3, t12, wkr, wki);
          conv_nested_1_4(s7, s8, t7, t8, wkr, wki);

          radix_16_notwd_first(ps,s,l);
        }
    }
}

void radix_16_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                                y_size_t  jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,t8r,
  t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i,u8r,
  u8i,u9r,u9i,u10r,u10i,u11r,u11i,u12r,u12i,u13r,u13i,u14r,u14i,u15r,u15i;
  y_limb_t s0r,s0i,s1r,s1i,s2r,s2i,s3r,s3i,s4r,s4i,s5r,s5i,s6r,s6i,s7r,s7i,s8r,
  s8i,s9r,s9i,s10r,s10i,s11r,s11i,s12r,s12i,s13r,s13i,s14r,s14i,s15r,s15i;
  y_limb_t v0r,v0i,v1r,v1i,v2r,v2i,v3r,v3i,v4r,v4i,v5r,v5i,v6r,v6i,v7r,v7i,v8r,
  v8i,v9r,v9i,v10r,v10i,v11r,v11i,v12r,v12i,v13r,v13i,v14r,v14i,v15r,v15i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_limb_t wkr,wki;
  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[15];
  y_ptr px,pt,pu,ps,pv;

  ASSERT_ALIGNED_DOUBLE();


  /*  if(Y_PLAN[Y_NRADICES-1] != 8) return;*/

  nb=bs>>4;
  bi=(bs*ii)>>4;
  bj=((bs*jj)>>4)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      get_twiddle_factors_16_last(j);
      l=j<<4;
      radix_16_twd_last(u,pu,d1,l);
      radix_16_twd_last(v,pv,d2,l);
      get_twiddle_factors_16_last(i);
      l=i<<4;
      radix_16_twd_last(s,ps,d1,l);
      radix_16_twd_last(t,pt,d2,l);

      conv_nested( s0, u15, t0, v15, tw1r, tw1i);
      conv_nested_1_4( s4, u11, t4, v11, tw1r, tw1i);
      conv_nested_1_2( s8, u7, t8, v7, tw1r, tw1i);
      conv_nested_3_4( s12, u3, t12, v3, tw1r, tw1i);

      cplx_mul(wk,tw1,F_1_16);
      conv_nested(s1, u14, t1, v14, wkr, wki);
      conv_nested_1_4(s5, u10, t5, v10, wkr, wki);
      conv_nested_1_2(s9, u6, t9, v6, wkr, wki);
      conv_nested_3_4(s13, u2, t13, v2, wkr, wki);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
      conv_nested(s2, u13, t2, v13, wkr, wki);
      conv_nested_1_4(s6, u9, t6, v9, wkr, wki);
      conv_nested_1_2(s10, u5, t10, v5, wkr, wki);
      conv_nested_3_4(s14, u1, t14, v1, wkr, wki);

      cplx_mul(wk,tw1,F_3_16);
      conv_nested(s3, u12, t3, v12, wkr, wki);
      conv_nested_1_4(s7, u8, t7, v8, wkr, wki);
      conv_nested_1_2(s11, u4, t11, v4, wkr, wki);
      conv_nested_3_4(s15, u0, t15, v0, wkr, wki);

      radix_16_notwd_first(ps,s,l);
      l=j<<4;
      radix_16_notwd_first(pu,u,l);
    }
  if((i==j)&&(ni<nb))
    {
      get_twiddle_factors_16_last(i);
      l=i<<4;
      radix_16_twd_last(s,ps,d1,l);
      radix_16_twd_last(t,pt,d2,l);
      conv_nested(s0, s15, t0, t15, tw1r, tw1i);
      conv_nested_1_4(s4, s11, t4, t11, tw1r, tw1i);

      cplx_mul(wk,tw1,F_1_16);
      conv_nested(s1, s14, t1, t14, wkr, wki);
      conv_nested_1_4(s5, s10, t5, t10, wkr, wki);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
      conv_nested(s2, s13, t2, t13, wkr, wki);
      conv_nested_1_4(s6, s9, t6, t9, wkr, wki);

      cplx_mul(wk,tw1,F_3_16);
      conv_nested(s3, s12, t3, t12, wkr, wki);
      conv_nested_1_4(s7, s8, t7, t8, wkr, wki);

      radix_16_notwd_first(ps,s,l);
    }
}

void radix_16_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                   y_size_t jj)
{
  y_limb_t t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i,t8r,
  t8i,t9r,t9i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i;
  y_limb_t u0r,u0i,u1r,u1i,u2r,u2i,u3r,u3i,u4r,u4i,u5r,u5i,u6r,u6i,u7r,u7i,u8r,
  u8i,u9r,u9i,u10r,u10i,u11r,u11i,u12r,u12i,u13r,u13i,u14r,u14i,u15r,u15i;
  y_limb_t tw1r,tw1i,tw2r,tw2i,tw3r,tw3i,tw4r,tw4i,tw5r,tw5i,tw6r,tw6i,
  tw7r,tw7i,tw8r,tw8i,tw9r,tw9i,tw10r,tw10i,tw11r,tw11i,tw12r,tw12i,
  tw13r,tw13i,tw14r,tw14i,tw15r,tw15i;
  y_limb_t wkr,wki;
  y_size_t i,j,nb,bi,bj,ni,l,inc=Y_POWERS[15];
  y_ptr px,pt,pu;
  /*DEBUG_VARS;*/
  ASSERT_ALIGNED_DOUBLE();

  /* if(Y_PLAN[Y_NRADICES-1] != 8) return; */

  nb=bs>>4;
  bi=(bs*ii)>>4;
  bj=((bs*jj)>>4)+nb-1;
  for(i=bi,j=bj,ni=0;(ni<nb)&&(i<j);i++,j--,ni++)
    {
      get_twiddle_factors_16_last(j);
      l=j<<4;
      radix_16_twd_last(u,pu,d,l);
      get_twiddle_factors_16_last(i);
      l=i<<4;
      radix_16_twd_last(t,pt,d,l);

      square_nested(t0, u15, tw1r, tw1i);
      square_nested_1_4(t4, u11, tw1r, tw1i);
      square_nested_1_2(t8, u7, tw1r, tw1i);
      square_nested_3_4(t12, u3, tw1r, tw1i);

      cplx_mul(wk,tw1,F_1_16);
      square_nested(t1, u14, wkr, wki);
      square_nested_1_4(t5, u10, wkr, wki);
      square_nested_1_2(t9, u6, wkr, wki);
      square_nested_3_4(t13, u2, wkr, wki);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
      square_nested(t2, u13, wkr, wki);
      square_nested_1_4(t6, u9, wkr, wki);
      square_nested_1_2(t10, u5, wkr, wki);
      square_nested_3_4(t14, u1, wkr, wki);

      cplx_mul(wk,tw1,F_3_16);
      square_nested(t3, u12, wkr, wki);
      square_nested_1_4(t7, u8, wkr, wki);
      square_nested_1_2(t11, u4, wkr, wki);
      square_nested_3_4(t15, u0, wkr, wki);

      radix_16_notwd_first(pt,t,l);
      l=j<<4;
      radix_16_notwd_first(pu,u,l);
    }
  if((i==j) && (ni<nb))
    {
      get_twiddle_factors_16_last(i);
      l=i<<4;
      radix_16_twd_last(t,pt,d,l);

      square_nested(t0, t15, tw1r,tw1i);
      square_nested_1_4(t4, t11, tw1r,tw1i);

      cplx_mul(wk,tw1,F_1_16);
      square_nested(t1, t14, wkr,wki);
      square_nested_1_4(t5, t10, wkr,wki);

      wkr=(tw1r+tw1i)*F_1_8r;
      wki=(tw1i-tw1r)*F_1_8r;
      square_nested(t2, t13, wkr,wki);
      square_nested_1_4(t6, t9, wkr,wki);

      cplx_mul(wk,tw1,F_3_16);
      square_nested(t3, t12, wkr,wki);
      square_nested_1_4(t7, t8, wkr,wki);

      radix_16_notwd_first(pt,t,l);
    }
}

#else

void radix_16_dif_square_dit(y_ptr d, y_size_t nc)
{
  printf("Called radix_16_dif_square_dit , d=%p, nc=%d\n",d,nc);
  printf("Please compile again the file difdit_16.c with \n -DY_AVAL=4\n");
  exit(-1);
}

void radix_16_dif_mul_dit(y_ptr d1, y_ptr d2, y_size_t nc)
{
  printf("Called radix_16_dif_mul_dit , d1=%p, d2=%p, nc=%d\n",d1,d2,nc);
  printf("Please compile again the file difdit_16.c with \n -DY_AVAL=4\n");
  exit(-1);
}

void radix_16_dif_mul_dit_block(y_ptr d1, y_ptr d2, y_size_t bs, y_size_t ii,
                                y_size_t  jj)
{
  printf("Called radix_16_dif_mul_dit_block , d1=%p, d2=%p, bs=%d, ii=%d, jj=%d\n",d1,d2,bs,ii,jj);
  printf("Please compile again the file difdit_16.c with \n -DY_AVAL=4\n");
  exit(-1);
}

void radix_16_dif_square_dit_block(y_ptr d, y_size_t bs, y_size_t ii,
                                   y_size_t jj)
{
  printf("Called radix_16_dif_square_dit_block , d=%p, bs=%d, ii=%d, jj=%d\n",d,bs,ii,jj);
  printf("Please compile again the file difdit_16.c with \n -DY_AVAL=4\n");
  exit(-1);
}

#endif
/*$Id$*/








