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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "yeafft.h"
#include "mccomp.h"
#include "glucas.h"
#include "gsetup.h"
#include "ynorm.h"

#define get_pads_4 \
  pad2=(pad<<1);\
  pad3=pad+pad2;\
  pd0= x;\
  pd1= x + addr(pad<<1);\
  pd2= x + addr(pad2<<1);\
  pd3= x + addr(pad3<<1);


#define init_bjs_4 \
  bj0=N;\
  bj1=(b & (UL)3 )<<(Y_K+1);\
  bj2=((2*b) & (UL)3 ) <<(Y_K+1); \
  bj3=((3*b) & (UL)3 ) <<(Y_K+1);


#define init_twiddle_factors_4 \
	tw1r=1.0; tw1i=0.0;\
	tw2r=tw1r; tw2i=tw1i;\
	tw3r=tw1r; tw3i=tw1i;

#ifdef Y_KILL_BRANCHES
# define  load_two_to_minusphi(_i) \
	tt0[1] = two_to_minusphi[addr(_i)]; \
	tt1[1] = two_to_minusphi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	tt2[1] = two_to_minusphi[addr(((pad2>>(SHIFT_UPDATE-1)) + _i))];\
	tt3[1] = two_to_minusphi[addr(((pad3>>(SHIFT_UPDATE-1)) + _i))];

# define  load_two_to_phi(_i) \
	tt0[0] = two_to_phi[_i]; \
	tt1[0] = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	tt2[0] = two_to_phi[addr(((pad2>>(SHIFT_UPDATE-1)) + _i))];\
	tt3[0] = two_to_phi[addr(((pad3>>(SHIFT_UPDATE-1)) + _i))];

#else

# define  load_two_to_minusphi(_i) \
	ttmp0 = two_to_minusphi[addr(_i)]; \
	ttmp1 = two_to_minusphi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	ttmp2 = two_to_minusphi[addr(((pad2>>(SHIFT_UPDATE-1)) + _i))];\
	ttmp3 = two_to_minusphi[addr(((pad3>>(SHIFT_UPDATE-1)) + _i))];

# define  load_two_to_phi(_i) \
	ttp0 = two_to_phi[addr(_i)]; \
	ttp1 = two_to_phi[addr(((pad2>>SHIFT_UPDATE) + _i))];\
	ttp2 = two_to_phi[addr(((pad2>>(SHIFT_UPDATE-1)) + _i))];\
	ttp3 = two_to_phi[addr(((pad3>>(SHIFT_UPDATE-1)) + _i))];

#endif

/* get the basic twiddle from memory and computes its powers */
#define get_twiddle_factors_4 \
	  tw2r=(px[0]+px[1])*(px[0]-px[1]);\
	  tw2i=2.0*px[0]*px[1];\
	  tw1r=px[0];\
	  tw1i=px[1];\
	  tw3r=px[0]*tw2r - px[1]*tw2i;\
	  tw3i=px[0]*tw2i + px[1]*tw2r;\
          prefetch_data(px, 0);

#define Y_STEP 2


# define radix_4_twd_last_dit(_j)\
	  cplx_load_muladdsub_pp(t0,t1,tw2,_j,pd0,pd2);\
	  cplx_load_mulmuladdsub_pp(t2,t3,tw1,tw3,_j,pd1,pd3);\
	  \
	  cplx_addsub(t0,t2);\
	  cplx_mul_1_4_B_addsub(t1,t3);


#define radix_4_first_dif(_j)\
	cplx_addsub(t0,t2);\
	cplx_addsub(t1,t3);\
	cplx_addsub_store_p(_j,pd0,pd2,t0,t1);\
	cplx_mul_1_4_F_addsub_store_p(_j,pd1,pd3,t2,t3);


#define radix_4_first_dif_add(_d)\
	cplx_addsub(t0,t2);\
	cplx_addsub(t1,t3);\
	cplx_addsub(t0,t1);\
	cplx_local_to_data(_d, 0 ,t0);\
	cplx_local_to_data(_d, pad2 ,t1);\
	cplx_mul_1_4_F_addsub(t2,t3);\
	cplx_local_to_data_add(_d, pad ,t2);\
	cplx_local_to_data_add(_d, pad3 ,t3);



#ifdef SUM_CHECK
#define sum_4(_Sum) \
	_Sum += (t0r + t0i) + (t1r + t1i) + (t2r + t2i) + (t3r + t3i);
#endif

void dit_carry_norm_dif_4(y_ptr x, UL N, UL err_flag)
{
#ifdef Y_KILL_BRANCHES
  double Y_ALIGNED(16) aux[6],maxerr=0.0;
  double Y_ALIGNED(16) tt0[2],tt1[2],tt2[2],tt3[2];
#else

  BIG_DOUBLE hiinv=highinv, loinv=lowinv;
  BIG_DOUBLE maxerr=0.0,ttmpSmall=Gsmall,ttmpBig=Gbig,
                                  ttpSmall=Hsmall;
  BIG_DOUBLE ttmp0,ttmp1,ttmp2,ttmp3;
  BIG_DOUBLE ttp0,ttp1,ttp2,ttp3;
#endif
#if defined(TRICKY_ROUND) && !defined(USE_ASM)

  BIG_DOUBLE A=bigA,B=bigB;
#endif

  BIG_DOUBLE carry0=-2.0,carry1=0.0,carry2=0.0,carry3=0.0;
  BIG_DOUBLE tw1r,tw1i,tw2r,tw2i,tw3r,tw3i;
  BIG_DOUBLE t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i;
  BIG_DOUBLE *px;
  UL bj0,bj1,bj2,bj3;
  y_ptr pd0,pd1,pd2,pd3;
  UL pad=Y_LRIGHT[1],pad2,pad3;
  UL i,j,k,l,ll;
#ifdef Y_KILL_BRANCHES

  long int nolast;
  UL ib[2];
# ifndef USE_ASM

  UL issmall;
# endif

  aux[0]=highinv;
  aux[1]=lowinv;
  aux[2]=Gbig;
  aux[3]=Hsmall;
  aux[4]=bigA;
  aux[5]=0.5;
  ib[0]=(UL)(-c);
  ib[1]=(UL)(b);
#endif

  Err=0.0;
  get_pads_4;

  init_bjs_4;

  init_twiddle_factors_4;
  px=Y_TWDB[Y_NRADICES-2];

  for (i=0,j=0; j<pad2; j+=UPDATE,i++)
    {
      load_two_to_minusphi(i);
      load_two_to_phi(i);
      for (k=0; k<UPDATE; k+=2, px+=Y_STEP)
        {
          l=(k+j)>>1;
#ifdef Y_KILL_BRANCHES

          nolast = l - (pad-1);
          nolast >>=(BITS_PER_UL-1);
#endif

          ll=addr(k+j);
          radix_4_twd_last_dit( ll );
#ifdef SUM_CHECK

          sum_4(SumOut);
#endif

          if(err_flag)
            {
              cplx_carry_norm_check(t,0);
              cplx_carry_norm_check(t,1);
              cplx_carry_norm_check(t,2);
              cplx_carry_norm_last_check(t,3);
            }
          else
            {
              cplx_carry_norm(t,0);
              cplx_carry_norm(t,1);
              cplx_carry_norm(t,2);
              cplx_carry_norm_last(t,3);
            }
#ifdef SUM_CHECK
          sum_4(SumIn);
#endif

          radix_4_first_dif( ll );
          get_twiddle_factors_4;
        }
    }
  /* adjust the last carries */
  t0r=carry3;
  t1r=carry0;
  t2r=carry1;
  t3r=carry2;

  init_bjs_4;
  j=0;
  load_two_to_phi(j);
  cplx_carry(t,0);
  cplx_carry(t,1);
  cplx_carry(t,2);
  cplx_carry(t,3);
#ifdef SUM_CHECK

  sum_4(SumIn);
#endif
  /* add carry signal */
  radix_4_first_dif_add( x );
  if (err_flag)
    Err=maxerr;
  printf("*\n");
}
/*$Id$*/


























