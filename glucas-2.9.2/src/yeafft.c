/*$Id$*/
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor, Klaus Kastens
 
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
/************************************************************************
	Guillermo Ballester Valor
	Granada, Spain.
 
These routines are designed specially to multiply two big integers using
	FAST FOURIER TRANSFORM
 
This file includes the interface to be used by other good multiprecision
packages (as GNU-MP multiprecision package).
 
This can 'emulate' the functions mpn_mul(), mpn_mul_fft_full() and 
mpn_mul_fft_full_sqr() GMP functions. See GMP info about the usage of this 
functions, and the following lines to the usage of the routines.
 
The routine y_mul() has the most general use. It can be used as mpn_mul in GMP
package. The routines y_mpn_mul_fft_full() and y_mpn_mul_fft_sqr() are the 
emulations of mpn_mul_fft_full() and mpn_mul_fft_full_sqr respectively. 
There are some important differences about how the routines runs.   
 
*************************************************************************
The usage of y_mul() is as follows
 
The types used by calling routines:
	mp_limb_t : the basic limb, or digit. It uses to be an UNSIGNED
		    integer with the same amount of bits like general 
		    registers in processor. 
		    Example, in an intel x86 32 architecture, an unsigned 
		    int with 32 bits (in most C compilers uint or ulong).
	y_size_t : A type with suficient bits to get the extension in
		    limbs of a big integer. It uses to be int type.
	mp_ptr    : A pointer to mp_limb_t.
	mp_srcptr : A pointer to a const mp_limb_t.
 
How the big integers are stored:
	A big integer is an array of mp_limb_t pointed by mp_ptr.
	The least significant limb is the first element mp_ptr[0].
	The limbs are considered unsigned.
 
The arguments:
	mp_ptr op1 pointer to factor 1, constant.
	mp_ptr op2 pointer to factor 2, constant. Both op1 and op2 can
		be the same (for square options). They will be the same
		number if op1==op2 and n1==n2
	n1 number of significant limbs of factor1 
		   op1 could have more than n1 limbs alocated *but*
		   we only want to use the first n1 of them.
	n2 idem. for operator 2.
	mp_ptr opt pointer to target. It contains the product. It can
		   overlap op1 and/or op2.
		   *IMPORTANT*
		   opt must have allocate sufficient space for the product.
		   I.E. opt *MUST* have at least n1+n2 limbs already
		   allocate when calling the routine.
The output:
	On exit, y_mul_fft returns the most significant limb opt[n1+n2-1]
	(which can be zero). This information is redundant because of
	all the product is stored on opt[0].....opt[n1+n2-1].
 
	IF AN ERROR HAS DETECTED, then global var Y_ERR=1, otherwise 
	Y_ERR is set to 0. 
 
	NOTE: PERHAPS IT SHOULD BE MORE C0NVENIENT TO USE THE SUBROUTINES 
	BELOW. BECAUSE OF NATURE OF FLOATING POINT FFT, IT IS VERY IMPROBABLE,
	BUT NOT IMPOSSIBLE,IT RETURNS A FALSE RESULT AND Y_ERR WAS SET TO 1.
	IF THE SOURCE OPERANDS HAS BEEN ALREADY DISTROYED, IT IS A 
	CATASTROPHIC ERROR (WE HAVE NO MANNER TO RECOVER). 
 
	*******************************************************************
 
	The usage of y_mpn_mul_fft_full()
	
	int y_mpn_mul_fft_full(mp_ptr opt, mp_ptr op1, size_t n1, 
	                          mp_ptr op2, size_t n2)
  
	y_mpn_mul_fft_full() TRIES TO MUL the big-integer pointed by 'op1' 
	with n1 limbs by the big-integer 'op2' with n2 limbs and assign the 
	result to 'opt' if success. The target result 'opt' have to be 
	allocated at least n1+n2 LIMBS by the calling routine. The interface 
	has four almost independent control check errors. Because of the 
	nature of Floating Point FFT, the 100% correctness is not guarranted,
	so there is built four control check errors (two independent 
	(the default option) and two semi-independent we can optionally add) 
	to see whether the reults are correct. 
 
	There is still other test we can select, the maximum round off error. 
 
	For a 32-bits integer machine there is a risk, with the default 
	options about 2^(-64) of a false good result, in a 64 bit machine
	about 2^(-128). In the worst case, the risk is very small compared 
	to a hardware error induced by cosmic rays or other software effects, 
	it should give a false good result error every some tenths of 
	million years or so 8-)). 
 
	You can see the document YEAFFT.htm included in this package for more
	details
	
	The routine returns 0 if failure, and left all the operands untouched, 
	so you must try to multiply again by 'secure' all-integer method, 
	otherwise it returns 1 if success , and the target opt is set 
	properly.
 
	The calling sequence from GMP could be something like 
 
	if(!y_mpn_mul_fft_full(opt,op1,n1,op2,n2)) 
	                mpn_mul_fft_full(opt,op1,n1,op2,n2);
 
	*******************************************************
 
	The usage of y_mpn_mul_fft_full_sqr()
 
	int y_mpn_mul_fft_full_sqr(mp_limb_t opt, mp_limb_t  op, n)
	
	It squares the big-integer pointed by op with n LIMBS and 
	assign the result, if success, to opt.  
 
	y_mpn_mul_fft_full() call this routine whether it detects the 
	operands are the same.
 
	The routine returns 0 if failure, and left all the operands untouched, 
	so you must try to multiply again by 'secure' all-integer method, 
	otherwise it returns 1 if success , and the target opt is set 
	properly.
 
	The calling sequence from GMP could be something like 
 
	if(!y_mpn_mul_fft_full_sqr(opt,op1,n)) 
	                mpn_mul_fft_full_sqr(opt,op1,n);
 
	***********************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#define USE_GMP*/
/******************************************************************
 yeafft is the include file for FFT. Read it to obtain more 
 information about it 
******************************************************************/
#ifdef USE_GMP
# define FFT_CHECK

/*
   Uncomment if you want to use this feature 
   Indeed it is not necessary with FFT_CHECK
*/
/*# define Y_RANDOMIZE_DIGITS*/
#endif

#include "yeafft.h"

/*#define Y_KILL_BRANCHES*/

/* Uncomment to activate the roundoff check error */
/*#define Y_ROUNDOFF_CHECK 0.4*/
/*#define Y_DEBUG*/

/* Uncomment to activate the additional check */
/*#define FFT_CHECK2*/

/* This is a variant of GNU longlong.h */
#ifdef FFT_CHECK
# include "longlong.h"
#endif

#ifdef Y_MOD_TEST2
# define FFT_CHECK2
#endif

/* This is to use the round trick
#if defined(TRICKY_ROUND)
BIG_DOUBLE bigA=0.0, bigB=0.0;
#endif
*/
unsigned long Y_BITS,Y_MASQ,Y_LIMB,Y_MAX;
int Y_FFTERR=0;
BIG_DOUBLE Y_INVL,Y_TINV;
y_size_t last_nsize=0;
#ifdef Y_RANDOMIZE_DIGITS
mp_ptr Y_RN=NULL;
#endif
#ifdef Y_ROUNDOFF_CHECK
y_limb_t Y_ERR;
#endif

/******************* CHECK BLOCK *************************************/
#ifdef FFT_CHECK

/* Some assembler tricks to go fast          */
/*to pass from real to integer in a fast way*/

# if ((Y_TARGET == 1)||(Y_TARGET == 11)||(Y_TARGET == 12))&& defined(__GNUC__)
/* x86 target*/

#  define Y_FAST_INT(_i,_fl)                   \
__asm__ volatile("fistl %1"                    \
        : :"t"(_fl),"m"(_i):"memory");

#  define Y_SUM_MOD_NINE(_r,_s)                \
__asm__ volatile("addl %2,%0 \n   adcl $0,%0"  \
:"=&r"(_r) :"0"(_r),"g"(_s));

# else

#  define Y_FAST_INT(i,fl) i=(long int)(RINT(fl))

#  define Y_SUM_MOD_NINE(_r,_s)                 \
    _r += (_s);                                 \
    if( _r < _s) _r++;

# endif

/* To compute mod (2^BITS_PER_LIMB + 1) */

# define Y_SUM_MOD_ELEVEN(_h,_l)                \
   sub_ddmmss(_h,_l,0,_l,0,_h);                 \
   if(_h>1) add_ssaaaa(_h,_l,_h,_l,one,one);

/* this performs a big integer mod NINE=((2^BITS_PER_LIMB)-1) to
   check the result. This is the scholar named NINE'S PROOF 
   Is very easy and effective 
      if  S1 is op1 mod NINE
      and  S2 is op2 mod NINE
   then (op1*op2) mod NINE = S1*S2 mod NINE 
 
   Well, compute the mod NINE is so easy as to sum the limbs, so is a very
   fast proof. The chance of a false good result is very very small 
   (2^(-BITS_PER_LIMB)
*/

/* This function returns the mod NINE of the argument w */
mp_limb_t y_ninemod(mp_ptr w, y_size_t n)
{
  mp_limb_t aux=0;
  y_size_t i;
  for (i=0;i<n;i++)
    {
      aux+=w[i];
      if(aux < w[i])
        aux++; /* carry*/
    }
  return aux;
}

/*
   This function returns whether the first nine's proof is passed
   (mod 2^(BITS_PER_LIMB) - 1)
*/
int y_fntcheck(mp_limb_t s1, mp_limb_t s2, mp_limb_t sp)
{
  mp_limb_t a,h,l;
  umul_ppmm(h,l,s1,s2);
  a=h;
  h+=l;
  if(h<a)
    h++; /* carry */
  /*
  if(sp  !=h )
    {
      printf("YEAFFT: nine proof1 failed! incorrect results\n");
      printf("%lX %lX %lX\n",sp,h,sp-h);
    }
  */
  if(sp == h)
    return 1;
  else
    return 0;
  /*assert (sp == h);*/

}

# ifdef FFT_CHECK2
/*
   This function returns whether the second nine's proof is passed
   (mod 2^(Y_BITS) - 1)
*/
int y_fntcheck2(mp_limb_t s1, mp_limb_t s2, mp_limb_t sp)
{
  mp_limb_t h,l,q,r;
  umul_ppmm(h,l,s1,s2);
  if(h != 0)
    udiv_qrnnd(q,r,h,l,Y_MASQ);
  else
    r = l%(Y_MASQ);
  /*
  if(sp  != r)
    {
      printf("YEAFFT: nine proof2 failed! incorrect results\n");
      printf("%lX %lX %lX\n",sp,r,sp-r);
    }
  */
  if(sp == r)
    return 1;
  else
    return 0;
}
# endif

/* test mod (2^bits_per_limb + 1) */
int y_fntcheck3(mod_test *s1, mod_test *s2, mod_test *sp)
{
  mp_limb_t h=0,l=0,one=1;
  if(((s1->_bph) == 0) && ((s2->_bph) == 0))
    {
      umul_ppmm(h,l,s1->_bpl,s2->_bpl);
      Y_SUM_MOD_ELEVEN(h,l);
    }
  else if(s1->_bph)
    sub_ddmmss(h,l,one,one,s2->_bph,s2->_bpl);
  else
    sub_ddmmss(h,l,one,one,s1->_bph,s1->_bpl);

  if( (l == sp->_bpl) && (h == sp->_bph) )
    return 1;
  else
    return 0;
}

# ifdef FFT_CHECK2
int y_fntcheck4(mp_limb_t s1, mp_limb_t s2, mp_limb_t sp)
{
  mp_limb_t h,l,q,r,yp=Y_LIMB+1;
  umul_ppmm(h,l,s1,s2);
  if(h != 0)
    udiv_qrnnd(q,r,h,l,yp);
  else
    r = l%(yp);
  if(sp == r)
    return 1;
  else
    return 0;
}
# endif

/* it makes all the modular test  */
int y_mod_tests(mod_test *s1, mod_test *s2, mod_test *sp)
{
  if(!y_fntcheck(s1->_bm,s2->_bm,sp->_bm))
    {
      /*printf(" Failed check mod (2^bits_per_limb - 1)\n");*/
      return 0;
    }
  if(!y_fntcheck3(s1 ,s2 ,sp))
    {
      /*printf(" Failed check mod (2^bits_per_limb + 1)\n");*/
      return 0;
    }

# ifdef FFT_CHECK2
  if(!y_fntcheck2(s1->_sm,s2->_sm,sp->_sm))
    {
      /*printf(" Failed check mod (2^bits_per_double - 1)\n");*/
      return 0;
    }
  if(!y_fntcheck4(s1->_sp,s2->_sp,sp->_sp))
    {
      /*printf(" Failed check mod (2^bits_per_double + 1)\n");*/
      return 0;
    }
# endif
  return 1;
}

/* This allocs memory for fft of n limbs. This work area can be padded
   for cache eficiency */
y_ptr get_fftworkspace(y_size_t n, y_ptr *P)
{
  (*P) = ALLOC_DOUBLES(addr(n));
  return (y_ptr) ALIGN_DOUBLES(*P);
}

# ifdef Y_RANDOMIZE_DIGITS
void random_bits(y_size_t needed)
{
  y_size_t n;
  if(Y_RN != NULL)
    free (Y_RN);
  n=needed/BITS_PER_LIMB + 1;
  Y_RN = (mp_limb_t *)malloc(2 * n * BYTES_PER_LIMB);
  mpn_random( Y_RN, 2*n);
}
# endif


/*
   This routine cuts the big integer into chunks and pass them balanced
   to the float array to be used by fft's routines. 
   Here n1 is the lengh of op in limbs, n is the lengh of FFT*2.
*/
void y_limb_to_workspace( y_ptr w1, mp_ptr op, y_size_t n1, y_size_t n)
{
  mp_limb_t carry=0,auxh,auxl,res;
  int ibit=0,jbit,j=0,i=1;
# ifdef Y_KILL_BRANCHES

  long int ib,aux;
# endif

  auxl=op[0];
  do
    {
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & Y_MASQ;
# ifdef Y_KILL_BRANCHES

      ibit &= (BITS_PER_LIMB - 1);
      ib=((long int)(ibit-jbit))>>(BITS_PER_LIMB-1);
      auxh = (op[i] & ib);
      res+=(auxh &((1U<<ibit)-1U))<<(BITS_PER_LIMB-jbit);
      auxl = (auxl & ~(ib)) + auxh;
      i +=(1U & ib);
# else

      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxh=op[i++];
          res+=(auxh &((1U<<ibit)-1U))<<(BITS_PER_LIMB-jbit);
          auxl=auxh;
        }
# endif
      res+=carry;
# ifdef Y_KILL_BRANCHES

      ib = ((long int)(Y_MAX-res))>>(BITS_PER_LIMB-1);
      carry = (1U & ib);
      aux = res - (Y_LIMB & ib);
      w1[addr(j)]=(y_limb_t)(aux);
# else

      carry=(res>Y_MAX);
      if (carry)
        w1[addr(j)]=(y_limb_t)(res)-(y_limb_t)Y_LIMB;
      else
        w1[addr(j)]=(y_limb_t) res;
# endif

      j++;
    }
  while(i<n1);
  /* last limb */
  do
    {
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & Y_MASQ;
      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxl=0;
        }
      res+=carry;
# ifdef Y_KILL_BRANCHES

      ib = ((long int)(Y_MAX-res))>>(BITS_PER_LIMB-1);
      carry = (1U & ib);
      aux = res - (Y_LIMB & ib);
      w1[addr(j)]=(y_limb_t)(aux);
# else

      carry=(res>Y_MAX);
      if (carry)
        w1[addr(j)]=(y_limb_t)(res)-(y_limb_t)Y_LIMB;
      else
        w1[addr(j)]=(y_limb_t) res;
# endif

      j++;
    }
  while (auxl !=0 || carry !=0);
  for(;j<n;j++)
    w1[addr(j)]=0.0;
}

/*
   This routine cuts the big integer into chunks and pass them balanced
   to the float array to be used by fft's routines. 
   Here n1 is the lengh of op in limbs, n is the lengh of FFT*2.
   It computes op mod (2^(2*BITS_PER_LIMB) - 1)
   op mod (2^(2*bits_per_double) - 1)  
   This is useful to check the results.
*/
void y_limb_to_workspace_four_tests( y_ptr w1, mp_ptr op, y_size_t n1, y_size_t n, mod_test *sum)
{
  mp_limb_t carry=0,auxh,auxl,res,ym=Y_MASQ,yl=Y_LIMB,yp;
  mp_limb_t bph,bpl,bmh,bml,zero=0,one=1,sh,sl;
# ifdef FFT_CHECK2

  mp_limb_t sph,spl,smh,sml;
# endif

  int ibit=0,jbit,j=0,i=1;

  /* inits the mod-test vars*/
  bph=0;
  bpl=op[0];
  bmh=0;
  bml=0;
# ifdef FFT_CHECK2

  sph=0;
  spl=0;
  smh=0;
  sml=0;
# endif

  yp = Y_LIMB+1;

  auxl=op[0];
  do
    {
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & ym;
      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxh=op[i];
          if( i & one )
            add_ssaaaa(bmh,bml,bmh,bml,zero,auxh);
          else
            add_ssaaaa(bph,bpl,bph,bpl,zero,auxh);
          i++;
          res+=(auxh &((1U<<ibit)-1U))<<(BITS_PER_LIMB-jbit);
          auxl=auxh;
        }
# ifdef FFT_CHECK2
      if( j & one)
        add_ssaaaa(smh,sml,smh,sml,zero,res);
      else
        add_ssaaaa(sph,spl,sph,spl,zero,res);
# endif

      res+=carry;
      carry=(res>Y_MAX);
      if (carry)
        w1[addr(j)]=(y_limb_t)(res)-(y_limb_t)yl;
      else
        w1[addr(j)]=(y_limb_t) res;
      j++;
    }
  while(i<n1);
  /* last limb */
  do
    {
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & ym;
      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxl=0;
        }
# ifdef FFT_CHECK2
      if( j & one)
        add_ssaaaa(smh,sml,smh,sml,zero,res);
      else
        add_ssaaaa(sph,spl,sph,spl,zero,res);
# endif

      res+=carry;
      carry=(res>Y_MAX);
      if (carry)
        w1[addr(j)]=(y_limb_t)(res)-(y_limb_t)yl;
      else
        w1[addr(j)]=(y_limb_t) res;
      j++;
    }
  while (auxl !=0 || carry !=0);
  for(;j<n;j++)
    w1[addr(j)]=0.0;

  /* mod (2^bits_per_limb - 1) */
  add_ssaaaa(sh,sl,bmh,bml,bph,bpl);
  Y_SUM_MOD_NINE(sl,sh);
  sum->_bm = sl;

# ifdef FFT_CHECK2
  /* mod (2^bits_per_double - 1) */
  add_ssaaaa(sh,sl,smh,sml,sph,spl);
  if(sh != 0)
    udiv_qrnnd(sh,sl,sh,sl,ym);
  else
    sl %= ym;
  sum->_sm = sl;
# endif

  /* mod (2^bits_per_limb + 1) */
  Y_SUM_MOD_ELEVEN(bmh,bml);
  Y_SUM_MOD_ELEVEN(bph,bpl);
  sub_ddmmss(bph,bpl,bph,bpl,bmh,bml);
  if( bph > one)
    add_ssaaaa(bph,bpl,bph,bpl,one,one);
  sum->_bph = bph;
  sum->_bpl = bpl;

# ifdef FFT_CHECK2
  /* mod (2^bits_per_double +1) */
  if(sph != 0)
    udiv_qrnnd(sh,spl,sph,spl,yp);
  else
    spl %=  yp;
  if(smh != 0)
    udiv_qrnnd(sh,sml,smh,sml,yp);
  else
    sml %=  yp;
  aux = spl;
  spl -= sml;
  if( spl > aux)
    spl += yp;
  sum->_sp = spl;
# endif

}

/* the randomize version of y_limb_to_workspace()*/

void y_limb_to_workspace_rn_four_tests( y_ptr w1, mp_ptr op, y_size_t n1, y_size_t n, mp_ptr rn, mod_test *sumr)
{
  mp_limb_t auxh,auxl,res,bits=0,yl=Y_LIMB,ym=Y_MASQ,yp;
  mp_limb_t bph,bpl,bmh,bml,zero=0,one=1,sh,sl;
# ifdef FFT_CHECK2

  mp_limb_t sph,spl,smh,sml;
# endif

  int ibit=0,jbit,j=0,i=1;
  long int aux,sum=0,carry=0;

  bph=0;
  bpl=op[0];
  bmh=0;
  bml=0;
# ifdef FFT_CHECK2

  sph=0;
  spl=0;
  smh=0;
  sml=0;
# endif

  yp = Y_LIMB+1;

  auxl=op[0];
  do
    {
      if((j & (BITS_PER_LIMB -1))==0)
        {
          bits = *(rn++);
        }
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & ym;
      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxh=op[i];
          if( i & one )
            add_ssaaaa(bmh,bml,bmh,bml,zero,auxh);
          else
            add_ssaaaa(bph,bpl,bph,bpl,zero,auxh);
          i++;
          res+=(auxh &((1U<<ibit)-1U))<<(BITS_PER_LIMB-jbit);
          auxl=auxh;
        }
# ifdef FFT_CHECK2
      if( j & one)
        add_ssaaaa(smh,sml,smh,sml,zero,res);
      else
        add_ssaaaa(sph,spl,sph,spl,zero,res);
# endif

      aux = res + carry;
      if( bits & 1U )
        {
          if(sum > yl)
            {
              aux -= 2*yl;
              carry = 2;
            }
          else
            {
              aux -= yl;
              carry = 1;
            }
        }
      else
        {
          if(sum < (-yl))
            {
              aux += yl;
              carry = -1;
            }
          else
            carry=0;
        }
      bits >>=1;
      sum+=aux;
      w1[addr(j)]=(y_limb_t)(aux);
      j++;
    }
  while(i<n1);
  /* last limb, we don't randomize last limb */
  do
    {
      jbit=ibit;
      ibit+=Y_BITS;
      res=(auxl>>jbit) & ym;
      if(ibit>=(int)BITS_PER_LIMB)
        {
          ibit-=(int)BITS_PER_LIMB;
          auxl=0;
        }
# ifdef FFT_CHECK2
      if( j & one)
        add_ssaaaa(smh,sml,smh,sml,zero,res);
      else
        add_ssaaaa(sph,spl,sph,spl,zero,res);
# endif

      aux= res + carry;
      carry=( aux >Y_MAX);
      if (carry)
        {
          w1[addr(j)]=(y_limb_t)(aux)-(y_limb_t)yl;
          sum += (aux - yl);
        }
      else
        {
          w1[addr(j)]=(y_limb_t) aux;
          sum +=aux;
        }
      j++;
    }
  while (auxl !=0 || carry !=0);
  for(;j<n;j++)
    w1[addr(j)]=0.0;
  /* mod (2^bits_per_limb - 1) */
  add_ssaaaa(sh,sl,bmh,bml,bph,bpl);
  Y_SUM_MOD_NINE(sl,sh);
  sumr->_bm = sl;


# ifdef FFT_CHECK2
  /* mod (2^bits_per_double - 1) */
  add_ssaaaa(sh,sl,smh,sml,sph,spl);
  if(sh != 0)
    udiv_qrnnd(sh,sl,sh,sl,ym);
  else
    sl %= ym;
  sumr->_sm = sl;
# endif

  /* mod (2^bits_per_limb + 1) */
  Y_SUM_MOD_ELEVEN(bmh,bml);
  Y_SUM_MOD_ELEVEN(bph,bpl);
  sub_ddmmss(bph,bpl,bph,bpl,bmh,bml);
  if( bph > one)
    add_ssaaaa(bph,bpl,bph,bpl,one,one);
  sumr->_bph = bph;
  sumr->_bpl = bpl;

# ifdef FFT_CHECK2
  /* mod (2^bits_per_double +1) */
  if(sph != 0)
    udiv_qrnnd(sh,spl,sph,spl,yp);
  else
    spl %=  yp;
  if(smh != 0)
    udiv_qrnnd(sh,sml,smh,sml,yp);
  else
    sml %=  yp;
  aux = spl;
  spl -= sml;
  if( spl > aux)
    spl += yp;
  sumr->_sp = spl;
# endif

}

/*
   This routine performs the normalization, carry propagation and 
   pass back the float array to a big integer format.
   nt is the size of double array, nc is the size (in limbs) of
   target.
*/
void y_norm_and_carry(mp_ptr opt, y_ptr w, y_size_t nt ,y_size_t nc)
{
  y_limb_t carry=0.0,aux;
# ifdef TRICKY_ROUND

  y_limb_t A=bigA,B=bigB;
# endif

  int i,k,ibit=0,jbit;
  mp_limb_t res=0,yl=Y_LIMB;
  long int ss=0;
# ifdef Y_KILL_BRANCHES

  long int ib;
# endif
# ifdef Y_ROUNDOFF_CHECK

  y_limb_t eau,roundoff;

  Y_ERR=0.0;
# endif
  /* the main loop */
  for(i=0,k=0;(i<nc);k++)
    {
# ifdef Y_ROUNDOFF_CHECK
      eau=w[addr(k)]*Y_INVL;
      roundoff=fabs(eau - RINT(eau));
      if( Y_ERR < roundoff)
        Y_ERR =roundoff;
      carry +=eau;
# else

      carry += (w[addr(k)]*Y_INVL);
# endif

      aux=RINT(carry*Y_TINV);
      Y_FAST_INT(ss,(carry-aux*(y_limb_t)yl));
      carry=aux;
# ifdef Y_KILL_BRANCHES

      ib = ((long int)(ss))>>(BITS_PER_LIMB-1);
      ss += ( ib & yl);
      carry += ib;
# else

      if(ss<0L)
        {
          ss+=(long int)Y_LIMB;
          carry-=1.0;
        }
# endif
      jbit=ibit;
      ibit+=Y_BITS;
      res+=((mp_limb_t)(ss))<<jbit;
      if(ibit>=BITS_PER_LIMB)
        {
          ibit-=BITS_PER_LIMB;
          opt[i++]=res;
          res=((mp_limb_t)(ss))>>(Y_BITS-ibit);
        }
    }
  /*if(carry != 0.0 || res !=0) printf("Oops ! \n");*/
  /*for(;i<nc;i++) opt[i]=0;*/
}

/*
   This routine performs a normalization, carry propagation and convert the
   float array data to a big integer data 
   nt is the size of double array, nc is the size (in limbs) of
   target.
   It also computes the modular residues to be used to check
*/
void y_norm_and_carry_four_tests(mp_ptr opt, y_ptr w, y_size_t nt ,
                                 y_size_t nc, mod_test *sumr)
{
  y_limb_t carry=0.0,aux;
# ifdef TRICKY_ROUND

  y_limb_t A=bigA,B=bigB;
# endif

  int i,k,ibit=0,jbit;
  mp_limb_t res=0,yl=Y_LIMB,yp;
  mp_limb_t bph,bpl,bmh,bml,zero=0,one=1,sh,sl;
# ifdef FFT_CHECK2

  mp_limb_t sph,spl,smh,sml,ym=Y_MASQ;
# endif

  long int ss;
# ifdef Y_ROUNDOFF_CHECK

  y_limb_t eau,roundoff;

  Y_ERR=0.0;
# endif

  bph=0;
  bpl=0;
  bmh=0;
  bml=0;
# ifdef FFT_CHECK2

  sph=0;
  spl=0;
  smh=0;
  sml=0;
# endif

  yp = Y_LIMB+1;



  /* the main loop */
  for(i=0,k=0;(i<nc);k++)
    {
# ifdef Y_ROUNDOFF_CHECK
      eau=w[addr(k)]*Y_INVL;
      roundoff=fabs(eau - RINT(eau));
      if( Y_ERR < roundoff)
        Y_ERR =roundoff;
      carry += eau;
# else

      carry+= (w[addr(k)]*Y_INVL);
# endif

      aux=RINT(carry*Y_TINV);
      /*ss=(long int)(carry-aux*(y_limb_t)Y_LIMB);*/
      Y_FAST_INT(ss,(carry-aux*(y_limb_t)yl));
      carry=aux;
      if(ss<0L)
        {
          ss+=(long int)yl;
          carry-=1.0;
        }
      /*if(ss < 0) printf("Ooooooooops\n");*/
# ifdef FFT_CHECK2
      if( k & one)
        add_ssaaaa(smh,sml,smh,sml,zero,ss);
      else
        add_ssaaaa(sph,spl,sph,spl,zero,ss);
# endif

      jbit=ibit;
      ibit+=Y_BITS;
      res+=((mp_limb_t)(ss))<<jbit;
      if(ibit>=BITS_PER_LIMB)
        {
          ibit-=BITS_PER_LIMB;
          if( i & one )
            add_ssaaaa(bmh,bml,bmh,bml,zero,res);
          else
            add_ssaaaa(bph,bpl,bph,bpl,zero,res);
          opt[i++]=res;
          /*Y_SUM_MOD_NINE(s1,res,sa);*/
          /*
          sa = s1;
          s1 += res;
          if(s1 < sa) s1++;
          */
          res=((mp_limb_t)(ss))>>(Y_BITS-ibit);
        }
    }
  /* mod (2^bits_per_limb - 1) */
  add_ssaaaa(sh,sl,bmh,bml,bph,bpl);
  Y_SUM_MOD_NINE(sl,sh);
  sumr->_bm = sl;

# ifdef FFT_CHECK2
  /* mod (2^bits_per_double - 1) */
  add_ssaaaa(sh,sl,smh,sml,sph,spl);
  if(sh != 0)
    udiv_qrnnd(sh,sl,sh,sl,ym);
  else
    sl %= ym;
  sumr->_sm = sl;
# endif

  /* mod (2^bits_per_limb + 1) */
  Y_SUM_MOD_ELEVEN(bmh,bml);
  Y_SUM_MOD_ELEVEN(bph,bpl);
  sub_ddmmss(bph,bpl,bph,bpl,bmh,bml);
  if( bph > one)
    add_ssaaaa(bph,bpl,bph,bpl,one,one);
  sumr->_bph = bph;
  sumr->_bpl = bpl;

# ifdef FFT_CHECK2
  /* mod (2^bits_per_double +1) */
  if(sph != 0)
    udiv_qrnnd(sh,spl,sph,spl,yp);
  else
    spl %=  yp;
  if(smh != 0)
    udiv_qrnnd(sh,sml,smh,sml,yp);
  else
    sml %=  yp;
  aux = spl;
  spl -= sml;
  if( spl > aux)
    spl += yp;
  sumr->_sp = spl;
# endif

  /*if(carry != 0.0 || res !=0) printf("Oops ! %lX %lX\n",carry,res);*/
  /*for(;i<nc;i++) opt[i]=0;*/
}



/* This routine performs convolution */
void y_convfft(y_ptr w1, y_ptr w2, y_size_t n)
{
  /*
  y_fftf( w1, n);
  y_fftf( w2, n);
  y_conv( w1, w2);
  y_fftb( w1, n );
  */
  y_convolution(w1,w2,n);

}

/* This routine performs square */
void y_squarfft(y_ptr w1, y_size_t n)
{
  /*
  y_fftf( w1, n);
  y_squar(w1);
  y_fftb( w1, n);
  */
  y_auto_convolution(w1,n);
}



/************************** THE INTERFACE ROUTINE *************************
Multiply op1 and op2, every with n1 limbs and n2 limbs, opt has to have 
allocated space enough (n1+n2+1) 
***************************************************************************/
mp_limb_t y_mul(mp_ptr opt, mp_ptr op1, y_size_t n1, mp_ptr op2,
                y_size_t n2)
{
  y_size_t n,needed;
  y_ptr Y_W1,Y_W2=NULL,Y_P1=NULL,Y_P2=NULL;
  int ne;

# if defined(FFT_CHECK)

  mod_test sum1,sum2,sump;
# endif

  ne= (op1 != op2) || (n1 != n2);
  n=(n1>n2) ? n1 :n2;
  /*
      
      Y_BITS is the number of bits we can put in an element of FFT 
      In a 64 bits float IEEE-359 compliant machine uses to be 19 or 20

      A first aproximation, for random number is 
      Y_BITS = (int)(26.0 - 0.4*log(n*2*BITS_PER_LIMB))

  */
  if(n !=last_nsize)
    {
      needed=(2*n*BITS_PER_LIMB);
# ifndef Y_RANDOMIZE_DIGITS

      Y_BITS = (UL)((double)25.00 - 0.4*log((double)needed));
# else

      Y_BITS = (UL)((double)24.00 - 0.5*log((double)needed));
# endif

      Y_MASQ=((1L<<Y_BITS)-1L);
      Y_LIMB=(1L<<Y_BITS);
      Y_MAX=(1L<<(Y_BITS-1));
      Y_TINV=1.0/(y_limb_t)Y_LIMB;

      needed=(needed-1)/(y_size_t)Y_BITS+3;
# if defined(USE_GMP) && defined(Y_RANDOMIZE_DIGITS)

      if (((Y_LENGTH<<1) < needed) )
        random_bits(needed);
# endif

      /*
      inits y-fft and all trig data. The output is the length 
      of fft (in size of complex)
      */
      y_init(needed);

      Y_INVL=1.0/(y_limb_t)(Y_LENGTH);

# ifdef TRICKY_ROUND

      if(bigA==0.0)
        bigA=tricky();
      bigB=bigA;
# endif

    }
  last_nsize=n;
  /* now n is the lengh of FFT (in size of complex) */
  n=Y_LENGTH;

  /*
     Get workspace for fnt. It can be padded. Y_P is the addres of 
     allocated memory, Y_W is the modified addres aligned to 16 
  */

  Y_W1=get_fftworkspace(n<<1,&(Y_P1));
  if(ne)
    Y_W2=get_fftworkspace(n<<1,&(Y_P2));

  /* Init auxiliar vars */
# ifdef Y_RANDOMIZE_DIGITS

  y_limb_to_workspace_rn_four_tests(Y_W1, op1, n1, n<<1,Y_RN,&sum1);
  if (ne)
    y_limb_to_workspace_rn_four_tests(Y_W2, op2, n2, n<<1,Y_RN,&sum2);
# else

  y_limb_to_workspace_four_tests(Y_W1, op1, n1, n<<1,&sum1);
  if (ne)
    y_limb_to_workspace_four_tests(Y_W2, op2, n2, n<<1,&sum2);
# endif

  /* Make the convolution */
  if (ne)
    y_convfft(Y_W1, Y_W2, n);
  else
    y_squarfft(Y_W1, n);

  /* carry and normalize */
  y_norm_and_carry_four_tests(opt, Y_W1, n<<1, n1+n2,&sump);

# if defined(Y_ROUNDOFF_CHECK) && defined(Y_DEBUG)

  printf("Roundoff err: %f\n",Y_ERR);
# endif

  /* Free mem */
  free(Y_P1);
  if(ne)
    free(Y_P2);

# ifdef FFT_CHECK

  if(ne)
    {
      if( y_mod_tests(&sum1,&sum2,&sump)
#  ifdef Y_ROUNDOFF_CHECK
          && (Y_ERR < Y_ROUNDOFF_CHECK)
#  endif
        )
        Y_FFTERR=0;
      else
        Y_FFTERR=1;
    }
  else
    {
      if( y_mod_tests(&sum1,&sum1,&sump)
#  ifdef Y_ROUNDOFF_CHECK
          && (Y_ERR < Y_ROUNDOFF_CHECK)
#  endif
        )
        Y_FFTERR=0;
      else
        Y_FFTERR=1;
    }
# endif
  return opt[n1+n2-1];
}


/************************** THE INTERFACE ROUTINE *************************
Multiply op1 and op2, every with n limbs, opt has to have 
allocated space enough ( 2*n ) 
***************************************************************************/
int y_mpn_mul_fft_full(mp_ptr opt, mp_ptr op1, y_size_t n1,mp_ptr op2, y_size_t n2)
{
  y_size_t n,nsize,needed,i;
  y_ptr Y_W1,Y_W2=NULL,Y_P1=NULL,Y_P2=NULL;
  mp_ptr Y_OP=NULL;
  int save;

# if defined FFT_CHECK

  mod_test sum1,sum2,sump;
# endif

  /*
      
      Y_BITS is the number of bits we can put in an element of FFT 
      In a 64 bits float IEEE-359 compliant machine uses to be 19 or 20

      A first aproximation is 
      Y_BITS = (int)(25.0 - 0.4*log(n*2*BITS_PER_LIMB))

  */
  if( (op1 == op2) && (n1 == n2) )
    return y_mpn_mul_fft_full_sqr(opt, op1, n1);
  nsize= n1 + n2;

  /*See whether all the initialization work con be avoided */
  if(nsize  != last_nsize)
    {
      needed=(nsize*BITS_PER_LIMB);
# ifndef Y_RANDOMIZE_DIGITS

      Y_BITS = (UL)((double)25.00 - 0.4*log((double)needed));
# else

      Y_BITS = (UL)((double)24.00 - 0.5*log((double)needed));
# endif
# ifdef Y_DEBUG

      printf(" Needed : %d BITS: %ld\n",needed,Y_BITS);
# endif

      needed=(needed-1)/(y_size_t)Y_BITS+3;
      Y_MASQ=((1U<<Y_BITS)-1L);
      Y_LIMB=(1U<<Y_BITS);
      Y_MAX=(1U<<(Y_BITS-1));
      Y_TINV=1.0/(y_limb_t)Y_LIMB;

# if defined(USE_GMP) && defined(Y_RANDOMIZE_DIGITS)
      /* Generate random bits */
      if (((Y_LENGTH<<1) < needed) )
        random_bits(needed);
# endif

      /*
      inits y-fft and all trig data. The output is the length 
      of fft (in size of complexs)
      */
      y_init(needed);
      Y_INVL=1.0/(y_limb_t)(Y_LENGTH);
# ifdef TRICKY_ROUND

      if(bigA==0.0)
        bigA=tricky();
      bigB=bigA;
# endif

    }
  last_nsize=nsize;

  /* now n is the lengh of FFT (in size of complex) */
  n=Y_LENGTH;


  /*
     Get workspace for fnt. It can be padded. Y_P is the addres of 
     allocated memory, Y_W is the modified addres aligned to 16 
  */

  Y_W1=get_fftworkspace(n<<1,&(Y_P1));
  Y_W2=get_fftworkspace(n<<1,&(Y_P2));


  /* Init auxiliar vars */
# ifdef Y_RANDOMIZE_DIGITS

  y_limb_to_workspace_rn_four_tests(Y_W1, op1, n1, n<<1, Y_RN, &sum1);
  y_limb_to_workspace_rn_four_tests(Y_W2, op2, n2, n<<1, Y_RN, &sum2);
# else

  y_limb_to_workspace_four_tests(Y_W1, op1, n1, n<<1,&sum1);
  y_limb_to_workspace_four_tests(Y_W2, op2, n2, n<<1,&sum2);
# endif

  /* Make the convolution */
  y_convfft(Y_W1, Y_W2, n);

  /*
     See whether it has to save the results in a provisional memory
     to avoid problems if fails.
     It assumes the operands no overlap partialy
  */
  save = ((opt == op1)||(opt == op2));
  if( save )
    Y_OP=(mp_ptr) malloc (nsize*BYTES_PER_LIMB);

  /* carry and normalize */
  if(save)
    y_norm_and_carry_four_tests(Y_OP, Y_W1, n<<1,nsize,&sump);
  else
    y_norm_and_carry_four_tests(opt, Y_W1, n<<1,nsize,&sump);

# if defined(Y_ROUNDOFF_CHECK) && defined(Y_DEBUG)

  printf("Roundoff err: %f\n",Y_ERR);
# endif

  /* Free mem */
  free(Y_P1);
  free(Y_P2);

# ifdef FFT_CHECK

  if( y_mod_tests(&sum1,&sum2,&sump)
#  ifdef Y_ROUNDOFF_CHECK
      && (Y_ERR < Y_ROUNDOFF_CHECK)
#  endif
    )
    {
      if (save)
        {
          for (i=0;i<nsize;i++)
            opt[i]=Y_OP[i];
          free (Y_OP);
        }
      Y_FFTERR=0;
      return 1;
    }
  else
    {
      if (save)
        free (Y_OP);
      Y_FFTERR=1;
      return 0;
    }
# endif
}

/* the square version of mul. We save the second forward FFT */

int y_mpn_mul_fft_full_sqr(mp_ptr opt, mp_ptr op1, y_size_t nsize)
{
  y_size_t n,needed,i;
  y_ptr Y_W1,Y_P1=NULL;
  mp_ptr Y_OP=NULL;
  int save;

# if defined FFT_CHECK

  mod_test sum1,sump;
# endif

  /*
      
      Y_BITS is the number of bits we can put in an element of FFT 
      In a 64 bits float IEEE-359 compliant machine uses to be 19 or 20

      A first aproximation is 
      Y_BITS = (int)(25.0 - 0.4*log(n*2*BITS_PER_LIMB))

  */
  /*See whether all the initialization work con be avoided */
  if(nsize  != last_nsize)
    {
      needed=(nsize*BITS_PER_LIMB)<<1;
# ifndef Y_RANDOMIZE_DIGITS

      Y_BITS = (UL)((double)25.00 - 0.4*log((double)needed));
# else

      Y_BITS = (UL)((double)24.00 - 0.5*log((double)needed));
# endif
# ifdef Y_DEBUG

      printf(" Needed : %d BITS: %ld\n",needed,Y_BITS);
# endif

      needed=(needed-1)/(y_size_t)Y_BITS+3;
      Y_MASQ=((1U<<Y_BITS)-1L);
      Y_LIMB=(1U<<Y_BITS);
      Y_MAX=(1U<<(Y_BITS-1));
      Y_TINV=1.0/(y_limb_t)Y_LIMB;

# if defined(USE_GMP) && defined(Y_RANDOMIZE_DIGITS)
      /* Generate random bits */
      if (((Y_LENGTH<<1) < needed) )
        random_bits(needed);
# endif

      /*
      inits y-fft and all trig data. The output is the length 
      of fft (in size of complexs)
      */
      y_init(needed);
      Y_INVL=1.0/(y_limb_t)(Y_LENGTH);
# ifdef TRICKY_ROUND

      if(bigA==0.0)
        bigA=tricky();
      bigB=bigA;
# endif

    }
  last_nsize=nsize;

  /* now n is the lengh of FFT (in size of complex) */
  n=Y_LENGTH;


  /*
     Get workspace for fnt. It can be padded. Y_P is the addres of 
     allocated memory, Y_W is the modified addres aligned to 16 
  */

  Y_W1=get_fftworkspace(n<<1,&(Y_P1));

  /* Init auxiliar vars */
# ifdef Y_RANDOMIZE_DIGITS

  y_limb_to_workspace_rn_four_tests(Y_W1, op1, nsize, n<<1, Y_RN, &sum1);
# else

  y_limb_to_workspace_four_tests(Y_W1, op1, nsize, n<<1, &sum1);
# endif

  /* Make the convolution */
  y_squarfft(Y_W1, n);

  /*
     See whether it has to save the results in a provisional memory
     to avoid problems if fails.
     It assumes the operands no overlap partialy
  */
  save = (opt == op1);
  if( save )
    Y_OP=(mp_ptr) malloc (2*nsize*BYTES_PER_LIMB);

  /* carry and normalize */
  if(save)
    y_norm_and_carry_four_tests(Y_OP, Y_W1, n<<1,nsize<<1,&sump);
  else
    y_norm_and_carry_four_tests(opt, Y_W1, n<<1,nsize<<1,&sump);

# if defined(Y_ROUNDOFF_CHECK) && defined(Y_DEBUG)

  printf("Roundoff err: %f\n",Y_ERR);
# endif
  /* Free mem */
  free(Y_P1);

# ifdef FFT_CHECK

  if( y_mod_tests(&sum1,&sum1,&sump)
#  ifdef Y_ROUNDOFF_CHECK
      && (Y_ERR < Y_ROUNDOFF_CHECK)
#  endif
    )
    {
      if (save)
        {
          for (i=0;i<(nsize<<1);i++)
            opt[i]=Y_OP[i];
          free (Y_OP);
        }
      Y_FFTERR=0;
      return 1;
    }
  else
    {
      if (save)
        free (Y_OP);
      Y_FFTERR=1;
      return 0;
    }
# endif
}

#endif

/********************** END CHECK BLOCK ***********************************/


/*$Id$*/
