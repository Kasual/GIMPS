/*$Id$*/
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor
 
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


/*
Make three simultaneous complex muls
it makes _target *= _factor
_ti is the prefix name of target, _fi is the prefix of factors 
*/

#define cplx_three_muls(_t0,_t1,_t2,_f0,_f1,_f2)       \
{                                                      \
   y_limb_t a0,a1,a2,b0,b1,b2;                         \
/* Cycle  1 .mfi .mfi */                               \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 stall. Let's see the compiler solution */  \
/* Cycle  5 */                                         \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  - a0;                     \
   _t0##i = _t0##i * _f0##r  + b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  - a1;                     \
   _t1##i = _t1##i * _f1##r  + b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  - a2;                     \
   _t2##i = _t2##i * _f2##r  + b2;                     \
}


/*
Make four simultaneous complex muls 
it makes _target *= _factor
_ti is the prefix name of target, _fi is the prefix of factors 
*/

#define cplx_four_muls(_t0,_t1,_t2,_t3,_f0,_f1,_f2,_f3)\
{                                                      \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                   \
/* Cycle  1 .mfi .mfi*/                                \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 .mfi .mfi*/                                \
   a3 = _t3##i * _f3##i;                               \
   b3 = _t3##r * _f3##i;                               \
/* Cycle  5  Stall. Schedule work to compiler */       \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  - a0;                     \
   _t0##i = _t0##i * _f0##r  + b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  - a1;                     \
   _t1##i = _t1##i * _f1##r  + b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  - a2;                     \
   _t2##i = _t2##i * _f2##r  + b2;                     \
/* Cycle  9 .mfi .mfi*/                                \
   _t3##r = _t3##r * _f3##r  - a3;                     \
   _t3##i = _t3##i * _f3##r  + b3;                     \
}

/*
Make five simultaneous complex muls.
it makes _target *= _factor
_ti is the prefix name of target, _fi is the prefix of factors 
This is the optimal solution for Itanium, 
just two cycles per complex mul. 
*/

#define cplx_five_muls(_t0,_t1,_t2,_t3,_t4,_f0,_f1,_f2,_f3,_f4) \
{                                                      \
   y_limb_t a0,a1,a2,a3,a4,b0,b1,b2,b3,b4;             \
/* Cycle  1 .mfi .mfi*/                                \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 .mfi .mfi*/                                \
   a3 = _t3##i * _f3##i;                               \
   b3 = _t3##r * _f3##i;                               \
/* Cycle  5 .mfi .mfi*/                                \
   a4 = _t4##i * _f4##i;                               \
   b4 = _t4##r * _f4##i;                               \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  - a0;                     \
   _t0##i = _t0##i * _f0##r  + b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  - a1;                     \
   _t1##i = _t1##i * _f1##r  + b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  - a2;                     \
   _t2##i = _t2##i * _f2##r  + b2;                     \
/* Cycle  9 .mfi .mfi*/                                \
   _t3##r = _t3##r * _f3##r  - a3;                     \
   _t3##i = _t3##i * _f3##r  + b3;                     \
/* Cycle 10 .mfi .mfi*/                                \
   _t4##r = _t4##r * _f4##r  - a4;                     \
   _t4##i = _t4##i * _f4##r  + b4;                     \
}

/*
Make three simultaneous complex divs
it makes _target /= _factor 
_ti is the prefix name of target, _fi is the prefix of factors 
This is valid for trig complex twiddle factors
*/

#define cplx_three_divs(_t0,_t1,_t2,_f0,_f1,_f2)       \
{                                                      \
   y_limb_t a0,a1,a2,b0,b1,b2;                         \
/* Cycle  1 .mfi .mfi */                               \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 stall. Let's see the compiler solution */  \
/* Cycle  5 */                                         \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  + a0;                     \
   _t0##i = _t0##i * _f0##r  - b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  + a1;                     \
   _t1##i = _t1##i * _f1##r  - b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  + a2;                     \
   _t2##i = _t2##i * _f2##r  - b2;                     \
}


/*
Make four simultaneous complex divs 
it makes _target /= _factor
_ti is the prefix name of target, _fi is the prefix of factors 
This is valid for trig complex twiddle factors
*/

#define cplx_four_divs(_t0,_t1,_t2,_t3,_f0,_f1,_f2,_f3)\
{                                                      \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                   \
/* Cycle  1 .mfi .mfi*/                                \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 .mfi .mfi*/                                \
   a3 = _t3##i * _f3##i;                               \
   b3 = _t3##r * _f3##i;                               \
/* Cycle  5  Stall. Schedule work to compiler */       \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  + a0;                     \
   _t0##i = _t0##i * _f0##r  - b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  + a1;                     \
   _t1##i = _t1##i * _f1##r  - b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  + a2;                     \
   _t2##i = _t2##i * _f2##r  - b2;                     \
/* Cycle  9 .mfi .mfi*/                                \
   _t3##r = _t3##r * _f3##r  + a3;                     \
   _t3##i = _t3##i * _f3##r  - b3;                     \
}

/*
Make five simultaneous complex divs.
it makes _target /= _factor
_ti is the prefix name of target, _fi is the prefix of factors 
This is the optimal solution for Itanium, 
just two cycles per complex mul. 
This is valid for trig complex twiddle factors
*/

#define cplx_five_divs(_t0,_t1,_t2,_t3,_t4,_f0,_f1,_f2,_f3,_f4) \
{                                                      \
   y_limb_t a0,a1,a2,a3,a4,b0,b1,b2,b3,b4;             \
/* Cycle  1 .mfi .mfi*/                                \
   a0 = _t0##i * _f0##i;                               \
   b0 = _t0##r * _f0##i;                               \
/* Cycle  2 .mfi .mfi*/                                \
   a1 = _t1##i * _f1##i;                               \
   b1 = _t1##r * _f1##i;                               \
/* Cycle  3 .mfi .mfi*/                                \
   a2 = _t2##i * _f2##i;                               \
   b2 = _t2##r * _f2##i;                               \
/* Cycle  4 .mfi .mfi*/                                \
   a3 = _t3##i * _f3##i;                               \
   b3 = _t3##r * _f3##i;                               \
/* Cycle  5 .mfi .mfi*/                                \
   a4 = _t4##i * _f4##i;                               \
   b4 = _t4##r * _f4##i;                               \
/* Cycle  6 .mfi .mfi*/                                \
   _t0##r = _t0##r * _f0##r  + a0;                     \
   _t0##i = _t0##i * _f0##r  - b0;                     \
/* Cycle  7 .mfi .mfi*/                                \
   _t1##r = _t1##r * _f1##r  + a1;                     \
   _t1##i = _t1##i * _f1##r  - b1;                     \
/* Cycle  8 .mfi .mfi*/                                \
   _t2##r = _t2##r * _f2##r  + a2;                     \
   _t2##i = _t2##i * _f2##r  - b2;                     \
/* Cycle  9 .mfi .mfi*/                                \
   _t3##r = _t3##r * _f3##r  + a3;                     \
   _t3##i = _t3##i * _f3##r  - b3;                     \
/* Cycle 10 .mfi .mfi*/                                \
   _t4##r = _t4##r * _f4##r  + a4;                     \
   _t4##i = _t4##i * _f4##r  - b4;                     \
}



/* This macro makes two divmul of complex trig factor */
/*
  _q0 is the prefix name of first target quotient 
  _p0 is the prefix name of first target product 
  _fb0 is the prefix name of first factor/dividend 
  _fs0 is the prefix name of second factor/divisor 
  _q1, _p1, _fb1, _fs1, _q1 idem for second operation 
  NOTE: outputs cannot overlap inputs
*/
#define cplx_two_divmul(_q0,_p0,_fb0,_fs0,_q1,_p1,_fb1,_fs1) \
{                                                            \
  y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                          \
  /* Cycle 1 .mfi .mfi*/                                     \
  a0 = _fb0##r * _fs0##r;                                    \
  b0 = _fb0##i * _fs0##i;                                    \
  /* Cycle 2 .mfi .mfi*/                                     \
  a1 = _fb1##r * _fs1##r;                                    \
  b1 = _fb1##i * _fs1##i;                                    \
  /* Cycle 3 .mfi .mfi*/                                     \
  a2 = _fb0##i * _fs0##r;                                    \
  b2 = _fb0##r * _fs0##i;                                    \
  /* Cycle 4 .mfi .mfi*/                                     \
  a3 = _fb1##i * _fs1##r;                                    \
  b3 = _fb1##r * _fs1##i;                                    \
  /* Cycle 5 Stall */                                        \
  /* Cycle 6 .mfi .mfi*/                                     \
  _p0##r = a0 - b0;                                          \
  _q0##r = a0 + b0;                                          \
  /* Cycle 7 .mfi .mfi*/                                     \
  _p1##r = a1 - b1;                                          \
  _q1##r = a1 + b1;                                          \
  /* Cycle 8 .mfi .mfi*/                                     \
  _q0##i = a2 - b2;                                          \
  _p0##i = a2 + b2;                                          \
  /* Cycle 9 .mfi .mfi*/                                     \
  _q1##i = a3 - b3;                                          \
  _p1##i = a3 + b3;                                          \
}

/*
This macro completes all the needed twiddle factors
in a radix_6 reduction from powers 1 and 4.
 _t is the prefix name of all target power
 _tp is the prefix name of preloaded powers (1 and 4)
*/
#define cplx_complete_twiddle_6_ia64(_t,_tp) \
{                                            \
   y_limb_t a0,b0;                           \
   /* Cycle 1 .mfi .mfi */                   \
   _t##2i = _tp##1r * _tp##1i;               \
   _t##2r = _tp##1i * _tp##1i;               \
   /* Cycle 2 .mfi .mfi */                   \
   a0 = _tp##4i * _tp##1i;                   \
   b0 = _tp##4r * _tp##1i;                   \
   /* Cycle 3 .mfi .mfi */                   \
   _t##1r = _tp##1r;                         \
   _t##1i = _tp##1i;                         \
   /* Cycle 4 .mfi .mfi */                   \
   _t##4r = _tp##4r;                         \
   _t##4i = _tp##4i;                         \
   /* Cycle 5 Stall*/                        \
   /* Cycle 6 .mfi .mfi */                   \
   _t##2r = _tp##1r * _tp##1r - _t##2r;      \
   _t##2i += _t##2i;                         \
   /* Cycle 7 .mfi .mfi */                   \
   _t##3r = _tp##4r * _tp##1r + a0;          \
   _t##5r = _tp##4r * _tp##1r - a0;          \
   /* Cycle 8 .mfi .mfi */                   \
   _t##3i = _tp##4i * _tp##1r - b0;          \
   _t##5i = _tp##4i * _tp##1r + b0;          \
}

/*
This macro completes all the needed twiddle factors
in a radix_7 reduction from powers 1,4 and 6.
 _t is the prefix name of all target power
 _tp is the prefix name of preloaded powers (1,4 and 6)
*/
#define cplx_complete_twiddle_7_ia64(_t,_tp) \
{                                            \
   y_limb_t a0,b0;                           \
   /* Cycle 1 .mfi .mfi */                   \
   _t##2i = _tp##1r * _tp##1i;               \
   _t##2r = _tp##1i * _tp##1i;               \
   /* Cycle 2 .mfi .mfi */                   \
   a0 = _tp##4i * _tp##1i;                   \
   b0 = _tp##4r * _tp##1i;                   \
   /* Cycle 3 .mfi .mfi */                   \
   _t##1r = _tp##1r;                         \
   _t##1i = _tp##1i;                         \
   /* Cycle 4 .mfi .mfi */                   \
   _t##4r = _tp##4r;                         \
   _t##4i = _tp##4i;                         \
   /* Cycle 5 .mfi .mfi */                   \
   _t##6r = _tp##6r;                         \
   _t##6i = _tp##6i;                         \
   /* Cycle 6 .mfi .mfi */                   \
   _t##2r = _tp##1r * _tp##1r - _t##2r;      \
   _t##2i += _t##2i;                         \
   /* Cycle 7 .mfi .mfi */                   \
   _t##3r = _tp##4r * _tp##1r + a0;          \
   _t##5r = _tp##4r * _tp##1r - a0;          \
   /* Cycle 8 .mfi .mfi */                   \
   _t##3i = _tp##4i * _tp##1r - b0;          \
   _t##5i = _tp##4i * _tp##1r + b0;          \
}

#if !defined(Y_MAXIMUM)
/*
This macro completes all the needed twiddle factors
in a radix_8 reduction from powers 1,2 and 5.
 _t is the prefix name of all target power
 _tp is the prefix name of preloaded powers (1,2 and 5)
*/
# define cplx_complete_twiddle_8_ia64(_t,_tp) \
{                                            \
   y_limb_t a0,b0,a1,b1;                     \
  /* Cycle 1 .mfi .mfi */                    \
  a0 = _tp##5i * _tp##2i;                    \
  b0 = _tp##5i * _tp##1i;                    \
  /* Cycle 2 .mfi .mfi */                    \
  a1 = _tp##5r * _tp##2i;                    \
  b1 = _tp##5r * _tp##1i;                    \
  /* Cycle 3 .mfi .mfi */                    \
  _t##1r = _tp##1r;                          \
  _t##1i = _tp##1i;                          \
  /* Cycle 4 .mfi .mfi */                    \
  _t##2r = _tp##2r;                          \
  _t##2i = _tp##2i;                          \
  /* Cycle 5 .mfi .mfi */                    \
  _t##5r = _tp##5r;                          \
  _t##5i = _tp##5i;                          \
  /* Cycle 6 .mfi .mfi */                    \
  _t##7r = _tp##5r * _tp##2r - a0;           \
  _t##3r = _tp##5r * _tp##2r + a0;           \
  /* Cycle 7 .mfi .mfi */                    \
  _t##6r = _tp##5r * _tp##1r - b0;           \
  _t##4r = _tp##5r * _tp##1r + b0;           \
  /* Cycle 8 .mfi .mfi */                    \
  _t##7i = _tp##5i * _tp##2r + a1;           \
  _t##3i = _tp##5i * _tp##2r - a1;           \
  /* Cycle 9 .mfi .mfi */                    \
  _t##6i = _tp##5i * _tp##1r + b1;           \
  _t##4i = _tp##5i * _tp##1r - b1;           \
}

#else
/*
When defined Y_MAXIMUM there is nothing to
complete, only copy memory.
*/
# define cplx_complete_twiddle_8_ia64(_t,_tp) \
{                                             \
  _t##1r = _tp##1r;                           \
  _t##1i = _tp##1i;                           \
  _t##2r = _tp##2r;                           \
  _t##2i = _tp##2i;                           \
  _t##3r = _tp##3r;                           \
  _t##3i = _tp##3i;                           \
  _t##4r = _tp##4r;                           \
  _t##4i = _tp##4i;                           \
  _t##5r = _tp##5r;                           \
  _t##5i = _tp##5i;                           \
  _t##6r = _tp##6r;                           \
  _t##6i = _tp##6i;                           \
  _t##7r = _tp##7r;                           \
  _t##7i = _tp##7i;                           \
}
#endif

/*
This macro completes all the needed twiddle factors
in a radix_9 reduction from powers 1,2 and 6.
 _t is the prefix name of all target power
 _tp is the prefix name of preloaded powers (1,2 and 6)
*/
#define cplx_complete_twiddle_9_ia64(_t,_tp) \
{                                            \
   y_limb_t a0,b0,a1,b1;                     \
  /* Cycle 1 .mfi .mfi */                    \
  a0 = _tp##6i * _tp##2i;                    \
  b0 = _tp##6i * _tp##1i;                    \
  /* Cycle 2 .mfi .mfi */                    \
  a1 = _tp##6r * _tp##2i;                    \
  b1 = _tp##6r * _tp##1i;                    \
  /* Cycle 3 .mfi .mfi */                    \
  _t##3r = _tp##2i * _tp##1i;                \
  _t##3i = _tp##2i * _tp##1r;                \
  /* Cycle 4 .mfi .mfi */                    \
  _t##1r = _tp##1r;                          \
  _t##1i = _tp##1i;                          \
  /* Cycle 5 .mfi .mfi */                    \
  _t##2r = _tp##2r;                          \
  _t##2i = _tp##2i;                          \
  /* Cycle 6 .mfi .mfi */                    \
  _t##6r = _tp##6r;                          \
  _t##6i = _tp##6i;                          \
  /* Cycle 7 .mfi .mfi */                    \
  _t##8r = _tp##6r * _tp##2r - a0;           \
  _t##4r = _tp##6r * _tp##2r + a0;           \
  /* Cycle 8 .mfi .mfi */                    \
  _t##7r = _tp##6r * _tp##1r - b0;           \
  _t##5r = _tp##6r * _tp##1r + b0;           \
  /* Cycle 9 .mfi .mfi */                    \
  _t##8i = _tp##6i * _tp##2r + a1;           \
  _t##4i = _tp##6i * _tp##2r - a1;           \
  /* Cycle 10 .mfi .mfi */                   \
  _t##7i = _tp##6i * _tp##1r + b1;           \
  _t##5i = _tp##6i * _tp##1r - b1;           \
  /* Cycle 11 .mfi .mfi */                   \
  _t##3r = _tp##2r * _tp##1r - _t##3r;       \
  _t##3i = _tp##2r * _tp##1i + _t##3i;       \
}


/*
This macro performs the following tasks of a radix_n :
  1)Passes three preload pointers to actual pointers.
  2)Computes three new preload pointers
  3)Mul three twiddle factors by the preloaded data.
  4)Preload three complex data.
 
  _t0,_t1,_t2 are the prefix data target name.
  _i0,_i1,_i2 are the index of the preload data and twiddle factors
  _f is the prefix name of twiddle factors
  _ppd is the prefix name of actual pointers
  _ppdl is the prefix name of preload pointers
  _inc is the step between actual and preload pointers.
 
  The scheme is (all data are complex), from j=0 to 2:
 
  _ppd##_i(j) = _ppdl##_i(j); 
  _tj = _s##_i(j) * _f##_i(j);
  _ppdl##_i(j) += _inc;
  _s##_i(j) = *( _ppdl##_i(j)) 
*/
#define cplx_three_muls_fftn_preload(_t0,_t1,_t2,_i0,_i1,_i2,_s,_f,_ppd,_ppdl,_inc)\
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   _t0##r = _s##_i0##i * _f##_i0##i;                               \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   _t0##i = _s##_i0##r * _f##_i0##i;                               \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   _t1##r = _s##_i1##i * _f##_i1##i;                               \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   _t1##i = _s##_i1##r * _f##_i1##i;                               \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   _t2##r = _s##_i2##i * _f##_i2##i;                               \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   _t2##i = _s##_i2##r * _f##_i2##i;                               \
/* Cycle  4 .mfi .mfi */                                           \
   _ppdl##_i0##r += _inc;                                          \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _ppdl##_i1##r += _inc;                                          \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _t0##r = _s##_i0##r * _f##_i0##r  - _t0##r;                     \
   _ppdl##_i2##r += _inc;                                          \
   _t0##i = _s##_i0##i * _f##_i0##r  + _t0##i;                     \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _s##_i1##r * _f##_i1##r  - _t1##r;                     \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _s##_i1##i * _f##_i1##r  + _t1##i;                     \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _s##_i2##r * _f##_i2##r  - _t2##r;                     \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _s##_i2##i * _f##_i2##r  + _t2##i;                     \
/* Cycle  9 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \



/*
This macro performs the following tasks of a radix_n :
1)Passes four preload pointers to actual pointers.
2)Computes four new preload pointers
3)Mul four twiddle factors by the preloaded data.
4)Preload four complex data.

_t0,_t1,_t2,_t3 are the prefix data target name.
_i0,_i1,_i2,_i3 are the index of the preload data and twiddle factors
_f is the prefix name of twiddle factors
_ppd is the prefix name of actual pointers
_ppdl is the prefix name of preload pointers
_inc is the step between actual and preload pointers.

The scheme is (all data are complex), from j=0 to 3:

    _ppd##_i(j) = _ppdl##_i(j);
                  _tj = _s##_i(j) * _f##_i(j);
                        _ppdl##_i(j) += _inc;
                                        _s##_i(j) = *( _ppdl##_i(j))
                                                    */

#define cplx_four_muls_fftn_preload(_t0,_t1,_t2,_t3,_i0,_i1,_i2,_i3,_s,_f,_ppd,_ppdl,_inc)\
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   _t0##r = _s##_i0##i * _f##_i0##i;                               \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   _t0##i = _s##_i0##r * _f##_i0##i;                               \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   _t1##r = _s##_i1##i * _f##_i1##i;                               \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   _t1##i = _s##_i1##r * _f##_i1##i;                               \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   _t2##r = _s##_i2##i * _f##_i2##i;                               \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   _t2##i = _s##_i2##r * _f##_i2##i;                               \
/* Cycle  4 .mfi .mfi */                                           \
   _ppd##_i3##r = _ppdl##_i3##r;                                   \
   _t3##r = _s##_i3##i * _f##_i3##i;                               \
   _ppdl##_i0##r += _inc;                                          \
   _ppd##_i3##i = _ppdl##_i3##i;                                   \
   _t3##i = _s##_i3##r * _f##_i3##i;                               \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mii */                                                \
   _ppdl##_i1##r += _inc;                                          \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _t0##r = _s##_i0##r * _f##_i0##r  - _t0##r;                     \
   _ppdl##_i2##r += _inc;                                          \
   _t0##i = _s##_i0##i * _f##_i0##r  + _t0##i;                     \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _s##_i1##r * _f##_i1##r  - _t1##r;                     \
   _ppdl##_i3##r += _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _s##_i1##i * _f##_i1##r  + _t1##i;                     \
   _ppdl##_i3##i += _inc;                                          \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _s##_i2##r * _f##_i2##r  - _t2##r;                     \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _s##_i2##i * _f##_i2##r  + _t2##i;                     \
/* Cycle  9 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _t3##r = _s##_i3##r * _f##_i3##r  - _t3##r;                     \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
   _t3##i = _s##_i3##i * _f##_i3##r  + _t3##i;                     \
/* Cycle  10 .mfi .mfi */                                          \
   _s##_i3##r = *(_ppdl##_i3##r);                                  \
   _s##_i3##i = *(_ppdl##_i3##i);                                  \

/*
                                                    This macro performs the following tasks of a radix_n :
                                                    1)Passes five preload pointers to actual pointers.
                                                    2)Computes five new preload pointers
                                                    3)Mul five twiddle factors by the preloaded data.
                                                    4)Preload five complex data.

                                                    _t0,_t1,_t2,_t3,_t4 are the prefix data target name.
                                                    _i0,_i1,_i2,_i3,_i4 are the index of the preload data and twiddle factors
                                                    _f is the prefix name of twiddle factors
                                                    _ppd is the prefix name of actual pointers
                                                    _ppdl is the prefix name of preload pointers
                                                    _inc is the step between actual and preload pointers.

                                                    The scheme is (all data are complex),j from 0 to 4:

                                                    _ppd##_i(j) = _ppdl##_i(j);
                                                                  _tj = _s##_i(j) * _f##_i(j);
                                                                        _ppdl##_i(j) += _inc;
                                                                                        _s##_i(j) = *( _ppdl##_i(j))
                                                                                                    */
#define cplx_five_muls_fftn_preload(_t0,_t1,_t2,_t3,_t4,_i0,_i1,_i2,_i3,_i4,_s,_f,_ppd,_ppdl,_inc)\
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   _t0##r = _s##_i0##i * _f##_i0##i;                               \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   _t0##i = _s##_i0##r * _f##_i0##i;                               \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   _t1##r = _s##_i1##i * _f##_i1##i;                               \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   _t1##i = _s##_i1##r * _f##_i1##i;                               \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   _t2##r = _s##_i2##i * _f##_i2##i;                               \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   _t2##i = _s##_i2##r * _f##_i2##i;                               \
/* Cycle  4 .mfi .mfi */                                           \
   _ppd##_i3##r = _ppdl##_i3##r;                                   \
   _t3##r = _s##_i3##i * _f##_i3##i;                               \
   _ppdl##_i0##r += _inc;                                          \
   _ppd##_i3##i = _ppdl##_i3##i;                                   \
   _t3##i = _s##_i3##r * _f##_i3##i;                               \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _ppd##_i4##r = _ppdl##_i4##r;                                   \
   _t4##r = _s##_i4##i * _f##_i4##i;                               \
   _ppdl##_i1##r += _inc;                                          \
   _ppd##_i4##i = _ppdl##_i4##i;                                   \
   _t4##i = _s##_i4##r * _f##_i4##i;                               \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _t0##r = _s##_i0##r * _f##_i0##r  - _t0##r;                     \
   _ppdl##_i2##r += _inc;                                          \
   _t0##i = _s##_i0##i * _f##_i0##r  + _t0##i;                     \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _s##_i1##r * _f##_i1##r  - _t1##r;                     \
   _ppdl##_i3##r += _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _s##_i1##i * _f##_i1##r  + _t1##i;                     \
   _ppdl##_i3##i += _inc;                                          \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _s##_i2##r * _f##_i2##r  - _t2##r;                     \
   _ppdl##_i4##r += _inc;                                          \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _s##_i2##i * _f##_i2##r  + _t2##i;                     \
   _ppdl##_i4##i += _inc;                                          \
/* Cycle  9 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _t3##r = _s##_i3##r * _f##_i3##r  - _t3##r;                     \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
   _t3##i = _s##_i3##i * _f##_i3##r  + _t3##i;                     \
/* Cycle  10 .mfi .mfi */                                          \
   _s##_i3##r = *(_ppdl##_i3##r);                                  \
   _t4##r = _s##_i4##r * _f##_i4##r  - _t4##r;                     \
   _s##_i3##i = *(_ppdl##_i3##i);                                  \
   _t4##i = _s##_i4##i * _f##_i4##r  + _t4##i;                     \
/* Cycle  11 .mfi .mfi */                                          \
   _s##_i4##r = *(_ppdl##_i4##r);                                  \
   _s##_i4##i = *(_ppdl##_i4##i);                                  \


/*
                                                                                                    This macro performs the following tasks of a radix_n :
                                                                                                    1)Passes four preload pointers to actual pointers.
                                                                                                    2)Computes four new preload pointers
                                                                                                    3)Mul four twiddle factors by the data. The twiddle factors are consumed
                                                                                                    after use. So, _t MUST be the same data than _f.
                                                                                                    4)Preload four complex data.

                                                                                                    _t0,_t1,_t2,_t3 are the prefix data target name.
                                                                                                    _i0,_i1,_i2,_i3 are the index of the preload data and twiddle factors
                                                                                                    _f is the prefix name of twiddle factors
                                                                                                    _ppd is the prefix name of actual pointers
                                                                                                    _ppdl is the prefix name of preload pointers
                                                                                                    _inc is the step between actual and preload pointers.

                                                                                                    The scheme is (all data are complex),j from 0 to 4:

                                                                                                    _ppd##_i(j) = _ppdl##_i(j);
                                                                                                                  _tj = _s##_i(j) * _f##_i(j);
                                                                                                                        _ppdl##_i(j) += _inc;
                                                                                                                                        _s##_i(j) = *( _ppdl##_i(j))

                                                                                                                                                    NOTE: BE CAREFUL WITH OVERLAPS. IT IS IMPORTANT THAT
                                                                                                                                                    j = _i(j)
                                                                                                                                                        */
#define cplx_four_muls_fftn_preload_consumed_trigs(_t0,_t1,_t2,_t3,_i0,_i1,_i2,_i3,_s,_f,_ppd,_ppdl,_inc)\
{                                                                  \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                               \
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   a0 = _f##_i0##i * _s##_i0##i;                                   \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   b0 = _f##_i0##r * _s##_i0##i;                                   \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   a1 = _f##_i1##i * _s##_i1##i;                                   \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   b1 = _f##_i1##r * _s##_i1##i;                                   \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   a2 = _f##_i2##i * _s##_i2##i;                                   \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   b2 = _f##_i2##r * _s##_i2##i;                                   \
/* Cycle  4 .mfi .mfi */                                           \
   _ppd##_i3##r = _ppdl##_i3##r;                                   \
   a3 = _f##_i3##i * _s##_i3##i;                                   \
   _ppdl##_i0##r += _inc;                                          \
   _ppd##_i3##i = _ppdl##_i3##i;                                   \
   b3 = _f##_i3##r * _s##_i3##i;                                   \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _ppdl##_i1##r += _inc;                                          \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _t0##r = _f##_i0##r * _s##_i0##r  - a0;                         \
   _ppdl##_i2##r += _inc;                                          \
   _t0##i = _f##_i0##i * _s##_i0##r  + b0;                         \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _f##_i1##r * _s##_i1##r  - a1;                         \
   _ppdl##_i3##r += _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _f##_i1##i * _s##_i1##r  + b1;                         \
   _ppdl##_i3##i += _inc;                                          \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _f##_i2##r * _s##_i2##r  - a2;                         \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _f##_i2##i * _s##_i2##r  + b2;                         \
/* Cycle  9 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _t3##r = _f##_i3##r * _s##_i3##r  - a3;                         \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
   _t3##i = _f##_i3##i * _s##_i3##r  + b3;                         \
/* Cycle  10 .mmf*/                                                \
   _s##_i3##r = *(_ppdl##_i3##r);                                  \
   _s##_i3##i = *(_ppdl##_i3##i);                                  \
}

                                                                                                                                                        /*
                                                                                                                                                        This macro performs the following tasks of a radix_n :
                                                                                                                                                          1)Passes three preload pointers to actual pointers.
                                                                                                                                                          2)Computes three new preload pointers
                                                                                                                                                          3)Mul three twiddle factors by the data. The twiddle factors are consumed
                                                                                                                                                            after use. So, _t MUST be the same data than _f.
                                                                                                                                                          4)Preload three complex data.
                                                                                                                                                         
                                                                                                                                                          _t0,_t1,_t2 are the prefix data target name.
                                                                                                                                                          _i0,_i1,_i2 are the index of the preload data and twiddle factors
                                                                                                                                                          _f is the prefix name of twiddle factors
                                                                                                                                                          _ppd is the prefix name of actual pointers
                                                                                                                                                          _ppdl is the prefix name of preload pointers
                                                                                                                                                          _inc is the step between actual and preload pointers.
                                                                                                                                                         
                                                                                                                                                          The scheme is (all data are complex),j from 0 to 2:
                                                                                                                                                         
                                                                                                                                                          _ppd##_i(j) = _ppdl##_i(j); 
                                                                                                                                                          _tj = _s##_i(j) * _f##_i(j);
                                                                                                                                                          _ppdl##_i(j) += _inc;
                                                                                                                                                          _s##_i(j) = *( _ppdl##_i(j)) 
                                                                                                                                                         
                                                                                                                                                        NOTE: BE CAREFUL WITH OVERLAPS. IT IS IMPORTANT THAT 
                                                                                                                                                          j = _i(j)
                                                                                                                                                        */
#define cplx_three_muls_fftn_preload_consumed_trigs(_t0,_t1,_t2,_i0,_i1,_i2,_s,_f,_ppd,_ppdl,_inc)\
{                                                                  \
   y_limb_t a0,a1,a2,b0,b1,b2;                                     \
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   a0 = _f##_i0##i * _s##_i0##i;                                   \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   b0 = _f##_i0##r * _s##_i0##i;                                   \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   a1 = _f##_i1##i * _s##_i1##i;                                   \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   b1 = _f##_i1##r * _s##_i1##i;                                   \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   a2 = _f##_i2##i * _s##_i2##i;                                   \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   b2 = _f##_i2##r * _s##_i2##i;                                   \
/* Cycle  4 .mfi .mfi */                                           \
   _ppdl##_i0##r += _inc;                                          \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _t0##r = _f##_i0##r * _s##_i0##r  - a0;                         \
   _ppdl##_i1##r += _inc;                                          \
   _t0##i = _f##_i0##i * _s##_i0##r  + b0;                         \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _f##_i1##r * _s##_i1##r  - a1;                         \
   _ppdl##_i2##r += _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _f##_i1##i * _s##_i1##r  + b1;                         \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _f##_i2##r * _s##_i2##r  - a2;                         \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _f##_i2##i * _s##_i2##r  + b2;                         \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
}
                                                                                                                                                        /* This is a variant of former macro to preload backward */
#define cplx_three_muls_fftn_preload_consumed_trigs_down(_t0,_t1,_t2,_i0,_i1,_i2,_s,_f,_ppd,_ppdl,_inc)\
{                                                                  \
   y_limb_t a0,a1,a2,b0,b1,b2;                                     \
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   a0 = _f##_i0##i * _s##_i0##i;                                   \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   b0 = _f##_i0##r * _s##_i0##i;                                   \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   a1 = _f##_i1##i * _s##_i1##i;                                   \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   b1 = _f##_i1##r * _s##_i1##i;                                   \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   a2 = _f##_i2##i * _s##_i2##i;                                   \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   b2 = _f##_i2##r * _s##_i2##i;                                   \
/* Cycle  4 .mfi .mfi */                                           \
   _ppdl##_i0##r -= _inc;                                          \
   _ppdl##_i0##i -= _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _t0##r = _f##_i0##r * _s##_i0##r  - a0;                         \
   _ppdl##_i1##r -= _inc;                                          \
   _t0##i = _f##_i0##i * _s##_i0##r  + b0;                         \
   _ppdl##_i1##i -= _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t1##r = _f##_i1##r * _s##_i1##r  - a1;                         \
   _ppdl##_i2##r -= _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t1##i = _f##_i1##i * _s##_i1##r  + b1;                         \
   _ppdl##_i2##i -= _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   _t2##r = _f##_i2##r * _s##_i2##r  - a2;                         \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   _t2##i = _f##_i2##i * _s##_i2##r  + b2;                         \
/* Cycle  8 .mmf */                                                \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
}

                                                                                                                                                        /*
                                                                                                                                                        This macro performs the following tasks of a last DIT radix_n :
                                                                                                                                                          1)Passes four preload pointers to actual pointers.
                                                                                                                                                          2)Computes four new preload pointers
                                                                                                                                                          3)Mul four twiddle factors by the data. The twiddle factors are consumed
                                                                                                                                                            after use. So, _t HAVE TO be the same data than _f.
                                                                                                                                                            index _i0 and _i1, index _i2 amd _i3 MUST BE PAIRED.
                                                                                                                                                         
                                                                                                                                                          4)Preload four complex data.
                                                                                                                                                         
                                                                                                                                                          _t0,_t1,_t2,_t3 are the prefix data target name.
                                                                                                                                                          _i0,_i1,_i2,_i3 are the index of the preload data and twiddle factors
                                                                                                                                                          _f is the prefix name of twiddle factors
                                                                                                                                                          _ppd is the prefix name of actual pointers
                                                                                                                                                          _ppdl is the prefix name of preload pointers
                                                                                                                                                          _inc is the step between actual and preload pointers.
                                                                                                                                                         
                                                                                                                                                          As example, for paired _i0 and _i1 indexes:
                                                                                                                                                         
                                                                                                                                                          _ppd##_i0 = _ppdl##_i0; 
                                                                                                                                                          _ppd##_i1 = _ppdl##_i1; 
                                                                                                                                                          _t0 = _s##_i1 * _f##_i1;
                                                                                                                                                          _t1 = _s##_i0 * _f##_i0;
                                                                                                                                                          _ppdl##_i0 += _inc;
                                                                                                                                                          _ppdl##_i1 += _inc;
                                                                                                                                                          _s##_i0 = *( _ppdl##_i0); 
                                                                                                                                                          _s##_i1 = *( _ppdl##_i1); 
                                                                                                                                                         
                                                                                                                                                        */
#define cplx_four_muls_fftn_preload_consumed_trigs_paired(_t0,_t1,_t2,_t3,_i0,_i1,_i2,_i3,_s,_f,_ppd,_ppdl,_inc) \
{                                                                  \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                               \
/* Cycle  1 .mfi .mfi */                                           \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   a0 = _f##_i0##i * _s##_i0##i;                                   \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   b0 = _f##_i0##r * _s##_i0##i;                                   \
/* Cycle  2 .mfi .mfi */                                           \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   a1 = _f##_i1##i * _s##_i1##i;                                   \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   b1 = _f##_i1##r * _s##_i1##i;                                   \
/* Cycle  3 .mfi .mfi */                                           \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   a2 = _f##_i2##i * _s##_i2##i;                                   \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   b2 = _f##_i2##r * _s##_i2##i;                                   \
/* Cycle  4 .mfi .mfi */                                           \
   _ppd##_i3##r = _ppdl##_i3##r;                                   \
   a3 = _f##_i3##i * _s##_i3##i;                                   \
   _ppdl##_i0##r += _inc;                                          \
   _ppd##_i3##i = _ppdl##_i3##i;                                   \
   b3 = _f##_i3##r * _s##_i3##i;                                   \
   _ppdl##_i0##i += _inc;                                          \
/* Cycle  5 .mfi .mfi */                                           \
   _ppdl##_i1##r += _inc;                                          \
   _ppdl##_i1##i += _inc;                                          \
/* Cycle  6 .mfi .mfi */                                           \
   a0 = _f##_i0##r * _s##_i0##r  - a0;                             \
   _ppdl##_i2##r += _inc;                                          \
   b0 = _f##_i0##i * _s##_i0##r  + b0;                             \
   _ppdl##_i2##i += _inc;                                          \
/* Cycle  7 .mfi .mfi */                                           \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t0##r = _f##_i1##r * _s##_i1##r  - a1;                         \
   _ppdl##_i3##r += _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t0##i = _f##_i1##i * _s##_i1##r  + b1;                         \
   _ppdl##_i3##i += _inc;                                          \
/* Cycle  8 .mfi .mfi */                                           \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   a2 = _f##_i2##r * _s##_i2##r  - a2;                             \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   b2 = _f##_i2##i * _s##_i2##r  + b2;                             \
/* Cycle  9 .mfi .mfi */                                           \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _t2##r = _f##_i3##r * _s##_i3##r  - a3;                         \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
   _t2##i = _f##_i3##i * _s##_i3##r  + b3;                         \
/* Cycle  10 .mfi .mfi */                                          \
   _s##_i3##r = *(_ppdl##_i3##r);                                  \
   _t1##r = a0;                                                    \
   _s##_i3##i = *(_ppdl##_i3##i);                                  \
   _t1##i = b0;                                                    \
/* Cycle  11 .mfi .mfi */                                          \
   _t3##r = a2;                                                    \
   _t3##i = b2;                                                    \
}
                                                                                                                                                        /* This variant is for backward preload */
#define cplx_four_muls_fftn_preload_consumed_trigs_paired_down(_t0,_t1,_t2,_t3,_i0,_i1,_i2,_i3,_s,_f,_ppd,_ppdl,_inc)\
{                                                                  \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                               \
/* Cycle  1 */                                                     \
   _ppd##_i0##r = _ppdl##_i0##r;                                   \
   a0 = _f##_i0##i * _s##_i0##i;                                   \
   _ppd##_i0##i = _ppdl##_i0##i;                                   \
   b0 = _f##_i0##r * _s##_i0##i;                                   \
/* Cycle  2 */                                                     \
   _ppd##_i1##r = _ppdl##_i1##r;                                   \
   a1 = _f##_i1##i * _s##_i1##i;                                   \
   _ppd##_i1##i = _ppdl##_i1##i;                                   \
   b1 = _f##_i1##r * _s##_i1##i;                                   \
/* Cycle  3 */                                                     \
   _ppd##_i2##r = _ppdl##_i2##r;                                   \
   a2 = _f##_i2##i * _s##_i2##i;                                   \
   _ppd##_i2##i = _ppdl##_i2##i;                                   \
   b2 = _f##_i2##r * _s##_i2##i;                                   \
/* Cycle  4 */                                                     \
   _ppd##_i3##r = _ppdl##_i3##r;                                   \
   a3 = _f##_i3##i * _s##_i3##i;                                   \
   _ppdl##_i0##r -= _inc;                                          \
   _ppd##_i3##i = _ppdl##_i3##i;                                   \
   b3 = _f##_i3##r * _s##_i3##i;                                   \
   _ppdl##_i0##i -= _inc;                                          \
/* Cycle  5 */                                                     \
   _ppdl##_i1##r -= _inc;                                          \
   _ppdl##_i1##i -= _inc;                                          \
/* Cycle  6 */                                                     \
   a0 = _f##_i0##r * _s##_i0##r  - a0;                             \
   _ppdl##_i2##r -= _inc;                                          \
   b0 = _f##_i0##i * _s##_i0##r  + b0;                             \
   _ppdl##_i2##i -= _inc;                                          \
/* Cycle  7 */                                                     \
   _s##_i0##r = *(_ppdl##_i0##r);                                  \
   _t0##r = _f##_i1##r * _s##_i1##r  - a1;                         \
   _ppdl##_i3##r -= _inc;                                          \
   _s##_i0##i = *(_ppdl##_i0##i);                                  \
   _t0##i = _f##_i1##i * _s##_i1##r  + b1;                         \
   _ppdl##_i3##i -= _inc;                                          \
/* Cycle  8 */                                                     \
   _s##_i1##r = *(_ppdl##_i1##r);                                  \
   a2 = _f##_i2##r * _s##_i2##r  - a2;                             \
   _s##_i1##i = *(_ppdl##_i1##i);                                  \
   b2 = _f##_i2##i * _s##_i2##r  + b2;                             \
/* Cycle  9 */                                                     \
   _s##_i2##r = *(_ppdl##_i2##r);                                  \
   _t2##r = _f##_i3##r * _s##_i3##r  - a3;                         \
   _s##_i2##i = *(_ppdl##_i2##i);                                  \
   _t2##i = _f##_i3##i * _s##_i3##r  + b3;                         \
/* Cycle  10 */                                                    \
   _s##_i3##r = *(_ppdl##_i3##r);                                  \
   _t1##r = a0;                                                    \
   _s##_i3##i = *(_ppdl##_i3##i);                                  \
   _t1##i = b0;                                                    \
/* Cycle  11 */                                                    \
   _t3##r = a2;                                                    \
   _t3##i = b2;                                                    \
}


                                                                                                                                                        /*
                                                                                                                                                           This macro does:
                                                                                                                                                          1)Passes eight preload pointers to actual pointers.
                                                                                                                                                          2)Computes eight new preload pointers
                                                                                                                                                          3)Passes eight preload data to actual data
                                                                                                                                                         
                                                                                                                                                          _ti are the prefix name of actual data
                                                                                                                                                          _i0 are the indexes
                                                                                                                                                          _ppd are the prefix name of pointers to actual data
                                                                                                                                                          _ppdl are the prefix name of pointers to preloaded data
                                                                                                                                                           
                                                                                                                                                        */
#define cplx_preload_8_notwd(_t0,_t1,_t2,_t3,_t4,_t5,_t6,_t7,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7,_ppd,_ppdl,_inc)\
   /* Cycle  1 */                     \
   _ppd##_i0##r = ppdl##_i0##r;       \
   _t0##r = _s##_i0##r;               \
   _ppd##_i0##i = ppdl##_i0##i;       \
   _t0##i = _s##_i0##i;               \
   /* Cycle  2 */                     \
   _ppd##_i1##r = ppdl##_i1##r;       \
   _t1##r = _s##_i1##r;               \
   _ppdl##_i0##r +=_inc;              \
   _ppd##_i1##i = ppdl##_i1##i;       \
   _t1##i = _s##_i1##i;               \
   _ppdl##_i0##i +=_inc;              \
   /* Cycle  3 */                     \
   _ppd##_i2##r = ppdl##_i2##r;       \
   _t2##r = _s##_i2##r;               \
   _ppdl##_i1##r +=_inc;              \
   _ppd##_i2##i = ppdl##_i2##i;       \
   _t2##i = _s##_i2##i;               \
   _ppdl##_i1##i +=_inc;              \
   /* Cycle  4 */                     \
   _ppd##_i3##r = ppdl##_i3##r;       \
   _t3##r = _s##_i3##r;               \
   _ppdl##_i2##r +=_inc;              \
   _ppd##_i3##i = ppdl##_i3##i;       \
   _t3##i = _s##_i3##i;               \
   _ppdl##_i2##i +=_inc;              \
   /* Cycle  5 */                     \
   _ppd##_i4##r = ppdl##_i4##r;       \
   _t4##r = _s##_i4##r;               \
   _ppdl##_i3##r +=_inc;              \
   _ppd##_i4##i = ppdl##_i4##i;       \
   _t4##i = _s##_i4##i;               \
   _ppdl##_i3##i +=_inc;              \
   /* Cycle  6 */                     \
   _ppd##_i5##r = ppdl##_i5##r;       \
   _t5##r = _s##_i5##r;               \
   _ppdl##_i4##r +=_inc;              \
   _ppd##_i5##i = ppdl##_i5##i;       \
   _t5##i = _s##_i5##i;               \
   _ppdl##_i4##i +=_inc;              \
   /* Cycle  6 */                     \
   _ppd##_i6##r = ppdl##_i6##r;       \
   _t6##r = _s##_i6##r;               \
   _ppdl##_i5##r +=_inc;              \
   _ppd##_i6##i = ppdl##_i6##i;       \
   _t6##i = _s##_i6##i;               \
   _ppdl##_i5##i +=_inc;              \
   /* Cycle  6 */                     \
   _ppd##_i7##r = ppdl##_i7##r;       \
   _t7##r = _s##_i7##r;               \
   _ppdl##_i6##r +=_inc;              \
   _ppd##_i7##i = ppdl##_i7##i;       \
   _t7##i = _s##_i7##i;               \
   _ppdl##_i6##i +=_inc;              \
   /* Cycle  7 */                     \
   _s##_i0##r = *(_ppdl##_i0##r);     \
   _s##_i0##i = *(_ppdl##_i0##i);     \
   _ppdl##_i7##r +=_inc;              \
   /* Cycle  8 */                     \
   _s##_i1##r = *(_ppdl##_i1##r);     \
   _s##_i1##i = *(_ppdl##_i1##i);     \
   _ppdl##_i7##i +=_inc;              \
   /* Cycle  9 */                     \
   _s##_i2##r = *(_ppdl##_i2##r);     \
   _s##_i2##i = *(_ppdl##_i2##i);     \
   /* Cycle 10 */                     \
   _s##_i3##r = *(_ppdl##_i3##r);     \
   _s##_i3##i = *(_ppdl##_i3##i);     \
   /* Cycle 11 */                     \
   _s##_i4##r = *(_ppdl##_i4##r);     \
   _s##_i4##i = *(_ppdl##_i4##i);     \
   /* Cycle 12 */                     \
   _s##_i5##r = *(_ppdl##_i5##r);     \
   _s##_i5##i = *(_ppdl##_i5##i);     \
   /* Cycle 13 */                     \
   _s##_i6##r = *(_ppdl##_i6##r);     \
   _s##_i6##i = *(_ppdl##_i6##i);     \
   /* Cycle 14 */                     \
   _s##_i7##r = *(_ppdl##_i7##r);     \
   _s##_i7##i = *(_ppdl##_i7##i);


                                                                                                                                                        /*
                                                                                                                                                        This macro performs the following tasks in the last DIT first DIF
                                                                                                                                                        phase in a radix_9 reduction:
                                                                                                                                                          1)init the preload pointers
                                                                                                                                                          2)preload the data
                                                                                                                                                         
                                                                                                                                                          _t is the prefix name of the preloaded data
                                                                                                                                                          _pd is the prefix name of init pointers
                                                                                                                                                          _l is the prefix name of preload pointers
                                                                                                                                                        */
#define cplx_radix_9_dif_dit_init(_t,_pd,_l)    \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _l##1##r = _pd##1;                            \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = _pd##2;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = _pd##3;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = _pd##4;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = _pd##4 + 1;                        \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _l##5##r = _pd##5;                            \
  _t##3##i = *(_l##3##i);                       \
  _l##5##i = _pd##5 + 1;                        \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _l##6##r = _pd##6;                            \
  _t##4##i = *(_l##4##i);                       \
  _l##6##i = _pd##6 + 1;                        \
  /* Cycle 8 */                                 \
  _t##5##r = *(_l##5##r);                       \
  _l##7##r = _pd##7;                            \
  _t##5##i = *(_l##5##i);                       \
  _l##7##i = _pd##7 + 1;                        \
  /* Cycle 9 */                                 \
  _t##6##r = *(_l##6##r);                       \
  _l##8##r = _pd##8;                            \
  _t##6##i = *(_l##6##i);                       \
  _l##8##i = _pd##8 + 1;                        \
  /* Cycle 10 */                                \
  _t##7##r = *(_l##7##r);                       \
  _t##7##i = *(_l##7##i);                       \
  /* Cycle 11 */                                \
  _t##8##r = *(_l##8##r);                       \
  _t##8##i = *(_l##8##i);                       \
}

                                                                                                                                                        /*
                                                                                                                                                        The same than former macro but for radix_8 
                                                                                                                                                        */
#define cplx_radix_8_dif_dit_init(_t,_pd,_l)    \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _l##1##r = _pd##1;                            \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = _pd##2;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = _pd##3;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = _pd##4;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = _pd##4 + 1;                        \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _l##5##r = _pd##5;                            \
  _t##3##i = *(_l##3##i);                       \
  _l##5##i = _pd##5 + 1;                        \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _l##6##r = _pd##6;                            \
  _t##4##i = *(_l##4##i);                       \
  _l##6##i = _pd##6 + 1;                        \
  /* Cycle 8 */                                 \
  _t##5##r = *(_l##5##r);                       \
  _l##7##r = _pd##7;                            \
  _t##5##i = *(_l##5##i);                       \
  _l##7##i = _pd##7 + 1;                        \
  /* Cycle 9 */                                 \
  _t##6##r = *(_l##6##r);                       \
  _t##6##i = *(_l##6##i);                       \
  /* Cycle 10 */                                \
  _t##7##r = *(_l##7##r);                       \
  _t##7##i = *(_l##7##i);                       \
}

                                                                                                                                                        /*
                                                                                                                                                        This is to initiate a radix-4. This inits two fft4.
                                                                                                                                                        See comments on cplx_radix_9_dif_dit_init()
                                                                                                                                                        */
#define cplx_two_radix_4_dif_dit_init(_t,_l,_u,_m,_pd) \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _m##0##r = _pd##0 + 2;                        \
  _m##0##i = _pd##0 + 3;                        \
   /* Cycle 3 */                                \
  _t##0##r = *(_l##0##r);                       \
  _l##1##r = _pd##1;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 4 */                                 \
  _u##0##r = *(_m##0##r);                       \
  _m##1##r = _pd##1 + 2;                        \
  _u##0##i = *(_m##0##i);                       \
  _m##1##i = _pd##1 + 3;                        \
  /* Cycle 5 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##2##r = _pd##2;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 6 */                                 \
  _u##1##r = *(_m##1##r);                       \
  _m##2##r = _pd##2 + 2;                        \
  _u##1##i = *(_m##1##i);                       \
  _m##2##i = _pd##2 + 3;                        \
  /* Cycle 7 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##3##r = _pd##3;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 8 */                                 \
  _u##2##r = *(_m##2##r);                       \
  _m##3##r = _pd##3 + 2;                        \
  _u##2##i = *(_m##2##i);                       \
  _m##3##i = _pd##3 + 3;                        \
  /* Cycle 9 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _t##3##i = *(_l##3##i);                       \
  /* Cycle 10 */                                \
  _u##3##r = *(_m##3##r);                       \
  _u##3##i = *(_m##3##i);                       \
}


                                                                                                                                                        /*
                                                                                                                                                        The same than former macro but for radix_8. Now the pad is 1
                                                                                                                                                        This macro is optimized for use in last DIF first DIT pass 
                                                                                                                                                        */
#define cplx_radix_8_pad1_init(_t,_pd,_jj,_l)   \
{                                               \
  y_ptr ppp;                                    \
  ppp = _pd +  addr((_jj)<<4);                  \
  /* Cycle 1 */                                 \
  _l##0##r = ppp;                               \
  _l##0##i = ppp + 1;                           \
  /* Cycle 2 */                                 \
  _l##1##r = ppp + 2;                           \
  _l##1##i = ppp + 3;                           \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = ppp + 4;                           \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = ppp + 5;                           \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = ppp + 6;                           \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = ppp + 7;                           \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = ppp + 8;                           \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = ppp + 9;                           \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _l##5##r = ppp + 10;                          \
  _t##3##i = *(_l##3##i);                       \
  _l##5##i = ppp + 11;                          \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _l##6##r = ppp + 12;                          \
  _t##4##i = *(_l##4##i);                       \
  _l##6##i = ppp + 13;                          \
  /* Cycle 8 */                                 \
  _t##5##r = *(_l##5##r);                       \
  _l##7##r = ppp + 14;                          \
  _t##5##i = *(_l##5##i);                       \
  _l##7##i = ppp + 15;                          \
  /* Cycle 9 */                                 \
  _t##6##r = *(_l##6##r);                       \
  _t##6##i = *(_l##6##i);                       \
  /* Cycle 10 */                                \
  _t##7##r = *(_l##7##r);                       \
  _t##7##i = *(_l##7##i);                       \
}

                                                                                                                                                        /*
                                                                                                                                                        The same than former macro but for radix_7 
                                                                                                                                                        See comments on cplx_radix_9_dif_dit_init()
                                                                                                                                                        */
#define cplx_radix_7_dif_dit_init(_t,_pd,_l)    \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _l##1##r = _pd##1;                            \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = _pd##2;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = _pd##3;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = _pd##4;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = _pd##4 + 1;                        \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _l##5##r = _pd##5;                            \
  _t##3##i = *(_l##3##i);                       \
  _l##5##i = _pd##5 + 1;                        \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _l##6##r = _pd##6;                            \
  _t##4##i = *(_l##4##i);                       \
  _l##6##i = _pd##6 + 1;                        \
  /* Cycle 8 */                                 \
  _t##5##r = *(_l##5##r);                       \
  _t##5##i = *(_l##5##i);                       \
  /* Cycle 9 */                                 \
  _t##6##r = *(_l##6##r);                       \
  _t##6##i = *(_l##6##i);                       \
}

                                                                                                                                                        /*
                                                                                                                                                        The same than former macro but for radix_6 
                                                                                                                                                        See comments on cplx_radix_9_dif_dit_init()
                                                                                                                                                        */
#define cplx_radix_6_dif_dit_init(_t,_pd,_l)    \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _l##1##r = _pd##1;                            \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = _pd##2;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = _pd##3;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = _pd##4;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = _pd##4 + 1;                        \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _l##5##r = _pd##5;                            \
  _t##3##i = *(_l##3##i);                       \
  _l##5##i = _pd##5 + 1;                        \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _t##4##i = *(_l##4##i);                       \
  /* Cycle 8 */                                 \
  _t##5##r = *(_l##5##r);                       \
  _t##5##i = *(_l##5##i);                       \
}

                                                                                                                                                        /*
                                                                                                                                                        The same than former macro but for radix_5 
                                                                                                                                                        See comments on cplx_radix_9_dif_dit_init()
                                                                                                                                                        */
#define cplx_radix_5_dif_dit_init(_t,_pd,_l)    \
{                                               \
  /* Cycle 1 */                                 \
  _l##0##r = _pd##0;                            \
  _l##0##i = _pd##0 + 1;                        \
  /* Cycle 2 */                                 \
  _l##1##r = _pd##1;                            \
  _l##1##i = _pd##1 + 1;                        \
  /* Cycle 3 */                                 \
  _t##0##r = *(_l##0##r);                       \
  _l##2##r = _pd##2;                            \
  _t##0##i = *(_l##0##i);                       \
  _l##2##i = _pd##2 + 1;                        \
  /* Cycle 4 */                                 \
  _t##1##r = *(_l##1##r);                       \
  _l##3##r = _pd##3;                            \
  _t##1##i = *(_l##1##i);                       \
  _l##3##i = _pd##3 + 1;                        \
  /* Cycle 5 */                                 \
  _t##2##r = *(_l##2##r);                       \
  _l##4##r = _pd##4;                            \
  _t##2##i = *(_l##2##i);                       \
  _l##4##i = _pd##4 + 1;                        \
  /* Cycle 6 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _t##3##i = *(_l##3##i);                       \
  /* Cycle 7 */                                 \
  _t##4##r = *(_l##4##r);                       \
  _t##4##i = *(_l##4##i);                       \
}



                                                                                                                                                        /*
                                                                                                                                                           Macro to make a first radix-8 pass with twiddle factors and prefetch data 
                                                                                                                                                           _t is the locale prefix variable name
                                                                                                                                                           _f is the locale twidle factors prefix variable name
                                                                                                                                                           _q is the pointers prefix variable name
                                                                                                                                                           _i is the index
                                                                                                                                                        */

#define cplx_radix_8_load_all(_t,_i,_q,_l)      \
{                                               \
  unsigned long ii;                             \
  /* Cycle 0 */                                 \
  ii = _i + 1;                                  \
  /* Cycle 1 */                                 \
  _l##0##r = _q##0 + _i;                        \
  _l##0##i = _q##0 + ii;                        \
  /* Cycle 2 */                                 \
  _l##4##r = _q##4 + _i;                        \
  _l##4##i = _q##4 + ii;                        \
  /* Cycle 3 */                                 \
  _l##2##r = _q##2 + _i;                        \
  _t##0##r = *(_l##0##r);                       \
  _l##2##i = _q##2 + ii;                        \
  _t##0##i = *(_l##0##i);                       \
  /* Cycle 4 */                                 \
  _l##6##r = _q##6 + _i;                        \
  _t##4##r = *(_l##4##r);                       \
  _l##6##i = _q##6 + ii;                        \
  _t##4##i = *(_l##4##i);                       \
  /* Cycle 5 */                                 \
  _l##1##r = _q##1 + _i;                        \
  _t##2##r = *(_l##2##r);                       \
  _l##1##i = _q##1 + ii;                        \
  _t##2##i = *(_l##2##i);                       \
  /* Cycle 6 */                                 \
  _l##5##r = _q##5 + _i;                        \
  _t##6##r = *(_l##6##r);                       \
  _l##5##i = _q##5 + ii;                        \
  _t##6##i = *(_l##6##i);                       \
  /* Cycle 7 */                                 \
  _l##3##r = _q##3 + _i;                        \
  _t##1##r = *(_l##1##r);                       \
  _l##3##i = _q##3 + ii;                        \
  _t##1##i = *(_l##1##i);                       \
  /* Cycle 8 */                                 \
  _l##7##r = _q##7 + _i;                        \
  _t##5##r = *(_l##5##r);                       \
  _l##7##i = _q##7 + ii;                        \
  _t##5##i = *(_l##5##i);                       \
  /* Cycle 9 */                                 \
  _t##3##r = *(_l##3##r);                       \
  _t##3##i = *(_l##3##i);                       \
  /* Cycle 10 */                                \
  _t##7##r = *(_l##7##r);                       \
  _t##7##i = *(_l##7##i);                       \
}



                                                                                                                                                        /*
                                                                                                                                                         
                                                                                                                                                        This macro performs the following tasks in a first pass of a radix_8 reduction 
                                                                                                                                                          1) Mul the preloaded data by twiddle factors.
                                                                                                                                                          2) First pass of the last DIT with preloaded data.
                                                                                                                                                          3) Store it in actual data
                                                                                                                                                          4) set preload pointers to actual pointers
                                                                                                                                                          5) Compute new preload pointers
                                                                                                                                                         
                                                                                                                                                        _t is the prefix name of local target values
                                                                                                                                                        _s is the prefix name of local preloaded values
                                                                                                                                                        _f is the prefix name of twiddle factors
                                                                                                                                                        _i is the index 
                                                                                                                                                        _p is the prefix base pointer name 
                                                                                                                                                        _q is the direction prefix name to store in this loop  (loaded in the former)
                                                                                                                                                        _l is the direction prefix name to load in this loop 
                                                                                                                                                        */
#define cplx_radix_8_twd_pp_preload_pass1(_t,_s,_f,_i,_p,_q,_l) \
{                                               \
  y_limb_t a1,a2,a3,a4,a5,a6,a7;                \
  y_limb_t b1,b2,b3,b4,b5,b6,b7;                \
  unsigned long ii;                             \
  /* Cycle 1 */                                 \
  _q##0##r = _l##0##r;                          \
  a1 = _s##4##i * _f##4##i;                     \
  ii = _i + 1;                                  \
  _q##0##i = _l##0##i;                          \
  b1 = _s##4##r * _f##4##i;                     \
  /* Cycle 2 */                                 \
  _q##4##r = _l##4##r;                          \
  a2 = _s##2##i * _f##2##i;                     \
  _q##4##i = _l##4##i;                          \
  b2 = _s##2##r * _f##2##i;                     \
  /* Cycle 3 */                                 \
  _q##2##r = _l##2##r;                          \
  a3 = _s##6##i * _f##6##i;                     \
  _q##2##i = _l##2##i;                          \
  b3 = _s##6##r * _f##6##i;                     \
  /* Cycle 4 */                                 \
  _q##6##r = _l##6##r;                          \
  a4 = _s##1##i * _f##1##i;                     \
  _q##6##i = _l##6##i;                          \
  b4 = _s##1##r * _f##1##i;                     \
  /* Cycle 5 */                                 \
  _q##1##r = _l##1##r;                          \
  a1 = _s##4##r * _f##4##r - a1;                \
  _q##1##i = _l##1##i;                          \
  b1 = _s##4##i * _f##4##r + b1;                \
  /* Cycle 6 */                                 \
  _q##5##r = _l##5##r;                          \
  a2 = _s##2##r * _f##2##r - a2;                \
  _q##5##i = _l##5##i;                          \
  b2 = _s##2##i * _f##2##r + b2;                \
  /* Cycle 7 */                                 \
  _q##3##r = _l##3##r;                          \
  a3 = _s##6##r * _f##6##r - a3;                \
  _q##3##i = _l##3##i;                          \
  b3 = _s##6##i * _f##6##r + b3;                \
  /* Cycle 8 */                                 \
  _q##7##r = _l##7##r;                          \
  a4 = _s##1##r * _f##1##r - a4;                \
  _q##7##i = _l##7##i;                          \
  b4 = _s##1##i * _f##1##r + b4;                \
  /* Cycle 9 */                                 \
  a5 = _s##5##i * _f##5##i;                     \
  _l##0##r = (_p##0 + _i);                      \
  b5 = _s##5##r * _f##5##i;                     \
  _l##0##i = (_p##0 + ii);                      \
  /* Cycle 10 */                                \
  a6 = _s##3##i * _f##3##i;                     \
  _l##4##r = (_p##4 + _i);                      \
  b6 = _s##3##r * _f##3##i;                     \
  _l##4##i = (_p##4 + ii);                      \
  /* Cycle 11 */                                \
  a7 = _s##7##i * _f##7##i;                     \
  _l##2##r = (_p##2 + _i);                      \
  b7 = _s##7##r * _f##7##i;                     \
  _l##2##i = (_p##2 + ii);                      \
  /* Cycle 12 */                                \
  _t##1##r = _s##0##r - a1;                     \
  _l##6##r = (_p##6 + _i);                      \
  _t##1##i = _s##0##i - b1;                     \
  _l##6##i = (_p##6 + ii);                      \
  /* Cycle 13 */                                \
  _t##0##r = _s##0##r + a1;                     \
  _l##1##r = (_p##1 + _i);                      \
  _t##0##i = _s##0##i + b1;                     \
  _l##1##i = (_p##1 + ii);                      \
  /* Cycle 14 */                                \
  a5 = _s##5##r * _f##5##r - a5;                \
  _l##5##r = (_p##5 + _i);                      \
  b5 = _s##5##i * _f##5##r + b5;                \
  _l##5##i = (_p##5 + ii);                      \
  /* Cycle 15 */                                \
  a6 = _s##3##r * _f##3##r - a6;                \
  _l##3##r = (_p##3 + _i);                      \
  b6 = _s##3##i * _f##3##r + b6;                \
  _l##3##i = (_p##3 + ii);                      \
  /* Cycle 16 */                                \
  a7 = _s##7##r * _f##7##r - a7;                \
  _l##7##r = (_p##7 + _i);                      \
  b7 = _s##7##i * _f##7##r + b7;                \
  _l##7##i = (_p##7 + ii);                      \
  /* Cycle 17 */                                \
  _t##2##r = a2 + a3;                           \
  _t##2##i = b2 + b3;                           \
  /* Cycle 18 */                                \
  _t##3##r = a2 - a3;                           \
  _t##3##i = b2 - b3;                           \
  /* Cycle 19 */                                \
  _t##4##r = a4 + a5;                           \
  _t##4##i = b4 + b5;                           \
  /* Cycle 20 */                                \
  _t##5##r = a4 - a5;                           \
  _t##5##i = b4 - b5;                           \
  /* Cycle 19 */                                \
  _t##6##r = a6 + a7;                           \
  _t##6##i = b6 + b7;                           \
  /* Cycle 20 */                                \
  _t##7##r = a6 - a7;                           \
  _t##7##i = b6 - b7;                           \
}


                                                                                                                                                        /*
                                                                                                                                                           This macro does the same than prior macro but without mul by
                                                                                                                                                           twiddle factors.
                                                                                                                                                         
                                                                                                                                                           _t is the locale prefix variable name
                                                                                                                                                           _f is the locale twidle factors prefix variable name
                                                                                                                                                           _q is the pointers prefix variable name
                                                                                                                                                           _i is the index
                                                                                                                                                        */
#define cplx_radix_8_notwd_pp_preload_pass1(_t,_s,_i,_p,_q,_l)\
{                                                         \
  unsigned long ii;                                       \
  /* Cycle  1 */                                          \
  _q##0##r = _l##0##r;                                    \
  _t##0##r = _s##0##r + _s##4##r;                         \
  ii = _i + 1;                                            \
  _q##0##i = _l##0##i;                                    \
  _t##0##i = _s##0##i + _s##4##i;                         \
  /* Cycle  2 */                                          \
  _q##4##r = _l##4##r;                                    \
  _t##1##r = _s##0##r - _s##4##r;                         \
  _l##0##r = (_p##0 + _i);                                \
  _q##4##i = _l##4##i;                                    \
  _t##1##i = _s##0##i - _s##4##i;                         \
  _l##0##i = (_p##0 + ii);                                \
  /* Cycle  3 */                                          \
  _q##2##r = _l##2##r;                                    \
  _t##2##r = _s##2##r + _s##6##r;                         \
  _l##4##r = (_p##4 + _i);                                \
  _q##2##i = _l##2##i;                                    \
  _t##2##i = _s##2##i + _s##6##i;                         \
  _l##4##i = (_p##4 + ii);                                \
  /* Cycle  4 */                                          \
  _q##6##r = _l##6##r;                                    \
  _t##3##r = _s##2##r - _s##6##r;                         \
  _l##2##r = (_p##2 + _i);                                \
  _q##6##i = _l##6##i;                                    \
  _t##3##i = _s##2##i - _s##6##i;                         \
  _l##2##i = (_p##2 + ii);                                \
  /* Cycle  5 */                                          \
  _q##1##r = _l##1##r;                                    \
  _t##4##r = _s##1##r + _s##5##r;                         \
  _l##6##r = (_p##6 + _i);                                \
  _q##1##i = _l##1##i;                                    \
  _t##4##i = _s##1##i + _s##5##i;                         \
  _l##6##i = (_p##6 + ii);                                \
  /* Cycle  6 */                                          \
  _q##5##r = _l##5##r;                                    \
  _t##5##r = _s##1##r - _s##5##r;                         \
  _l##1##r = (_p##1 + _i);                                \
  _q##5##i = _l##5##i;                                    \
  _t##5##i = _s##1##i - _s##5##i;                         \
  _l##1##i = (_p##1 + ii);                                \
  /* Cycle  7 */                                          \
  _q##3##r = _l##3##r;                                    \
  _t##6##r = _s##3##r + _s##7##r;                         \
  _l##5##r = (_p##5 + _i);                                \
  _q##3##i = _l##3##i;                                    \
  _t##6##i = _s##3##i + _s##7##i;                         \
  _l##5##i = (_p##5 + ii);                                \
  /* Cycle  8 */                                          \
  _q##7##r = _l##7##r;                                    \
  _t##7##r = _s##3##r - _s##7##r;                         \
  _l##3##r = (_p##3 + _i);                                \
  _q##7##i = _l##7##i;                                    \
  _t##7##i = _s##3##i - _s##7##i;                         \
  _l##3##i = (_p##3 + ii);                                \
  /* Cycle 9 */                                           \
  _l##7##r = (_p##7 + _i);                                \
  _l##7##i = (_p##7 + ii);                                \
}

                                                                                                                                                        /*
                                                                                                                                                           This macro performs:
                                                                                                                                                           1) the second pass in a DIF radix-8 reduction
                                                                                                                                                           2) preload the data
                                                                                                                                                         
                                                                                                                                                         
                                                                                                                                                        Pass 2 for a radix_8 reduction. Forward 
                                                                                                                                                        _t is the prefix name of local preloaded
                                                                                                                                                        _s is the prefix name of local preloaded   
                                                                                                                                                        _q is the computed direction of preload
                                                                                                                                                        */
#define cplx_radix_8_pp_pass2_preload_F(_t,_s,_q)         \
{                                                         \
  double a0,b0;                                           \
  /* Cycle  1 */                                          \
  a0 = _t##2##r;                                          \
  b0 = _t##2##i;                                          \
  /* Cycle  2 */                                          \
  _s##0##r = *(_q##0##r);                                 \
  _t##2##r = _t##0##r - a0;                               \
  _s##0##i = *(_q##0##i);                                 \
  _t##2##i = _t##0##i - b0;                               \
  /* Cycle  3 */                                          \
  _s##4##r = *(_q##4##r);                                 \
  _t##0##r += a0;                                         \
  _s##4##i = *(_q##4##i);                                 \
  _t##0##i += b0;                                         \
  /* Cycle  4 */                                          \
  a0 = _t##3##r;                                          \
  b0 = _t##3##i;                                          \
  /* Cycle  5 */                                          \
  _s##2##r = *(_q##2##r);                                 \
  _t##3##r = _t##1##r - b0;                               \
  _s##2##i = *(_q##2##i);                                 \
  _t##3##i = _t##1##i + a0;                               \
  /* Cycle  6 */                                          \
  _s##6##r = *(_q##6##r);                                 \
  _t##1##r += b0;                                         \
  _s##6##i = *(_q##6##i);                                 \
  _t##1##i -= a0;                                         \
  /* Cycle  7 */                                          \
  a0 = _t##6##r;                                          \
  b0 = _t##6##i;                                          \
  /* Cycle  8 */                                          \
  _s##1##r = *(_q##1##r);                                 \
  _t##6##r = _t##4##r - a0;                               \
  _s##1##i = *(_q##1##i);                                 \
  _t##6##i = _t##4##i - b0;                               \
  /* Cycle  9 */                                          \
  _s##5##r = *(_q##5##r);                                 \
  _t##4##r += a0;                                         \
  _s##5##i = *(_q##5##i);                                 \
  _t##4##i += b0;                                         \
  /* Cycle 10 */                                          \
  a0 = _t##7##r;                                          \
  b0 = _t##7##i;                                          \
  /* Cycle 11 */                                          \
  _s##3##r = *(_q##3##r);                                 \
  _t##7##r = _t##5##r - b0;                               \
  _s##3##i = *(_q##3##i);                                 \
  _t##7##i = _t##5##i + a0;                               \
  /* Cycle 12 */                                          \
  _s##7##r = *(_q##7##r);                                 \
  _t##5##r += b0;                                         \
  _s##7##i = *(_q##7##i);                                 \
  _t##5##i -= a0;                                         \
}

                                                                                                                                                        /*
                                                                                                                                                           This macro performs:
                                                                                                                                                           1) the second pass in a DIT radix-8 reduction
                                                                                                                                                           2) preload the data
                                                                                                                                                         
                                                                                                                                                         
                                                                                                                                                           _t is the prefix name of local preloaded
                                                                                                                                                           _s is the prefix name of local preloaded   
                                                                                                                                                           _q is the computed direction of preload
                                                                                                                                                        */
#define cplx_radix_8_pp_pass2_preload_B(_t,_s,_q)         \
{                                                         \
  double a0,b0;                                           \
  /* Cycle  1 */                                          \
  a0 = _t##2##r;                                          \
  b0 = _t##2##i;                                          \
  /* Cycle  2 */                                          \
  _s##0##r = *(_q##0##r);                                 \
  _t##2##r = _t##0##r - a0;                               \
  _s##0##i = *(_q##0##i);                                 \
  _t##2##i = _t##0##i - b0;                               \
  /* Cycle  3 */                                          \
  _s##4##r = *(_q##4##r);                                 \
  _t##0##r += a0;                                         \
  _s##4##i = *(_q##4##i);                                 \
  _t##0##i += b0;                                         \
  /* Cycle  4 */                                          \
  a0 = _t##3##r;                                          \
  b0 = _t##3##i;                                          \
  /* Cycle  5 */                                          \
  _s##2##r = *(_q##2##r);                                 \
  _t##3##r = _t##1##r + b0;                               \
  _s##2##i = *(_q##2##i);                                 \
  _t##3##i = _t##1##i - a0;                               \
  /* Cycle  6 */                                          \
  _s##6##r = *(_q##6##r);                                 \
  _t##1##r -= b0;                                         \
  _s##6##i = *(_q##6##i);                                 \
  _t##1##i += a0;                                         \
  /* Cycle  7 */                                          \
  a0 = _t##6##r;                                          \
  b0 = _t##6##i;                                          \
  /* Cycle  8 */                                          \
  _s##1##r = *(_q##1##r);                                 \
  _t##6##r = _t##4##r - a0;                               \
  _s##1##i = *(_q##1##i);                                 \
  _t##6##i = _t##4##i - b0;                               \
  /* Cycle  9 */                                          \
  _s##5##r = *(_q##5##r);                                 \
  _t##4##r += a0;                                         \
  _s##5##i = *(_q##5##i);                                 \
  _t##4##i += b0;                                         \
  /* Cycle 10 */                                          \
  a0 = _t##7##r;                                          \
  b0 = _t##7##i;                                          \
  /* Cycle 11 */                                          \
  _s##3##r = *(_q##3##r);                                 \
  _t##7##r = _t##5##r + b0;                               \
  _s##3##i = *(_q##3##i);                                 \
  _t##7##i = _t##5##i - a0;                               \
  /* Cycle 12 */                                          \
  _s##7##r = *(_q##7##r);                                 \
  _t##5##r -= b0;                                         \
  _s##7##i = *(_q##7##i);                                 \
  _t##5##i += a0;                                         \
}

                                                                                                                                                        /*
                                                                                                                                                           Macro to make the third radix-8 pass and store in a DIF 
                                                                                                                                                           _t is the locale prefix variable name
                                                                                                                                                           _q is the pointers prefix variable name
                                                                                                                                                           _i is the index
                                                                                                                                                        */
#define cplx_radix_8_pp_pass3_store_F(_t,_q)              \
{                                                         \
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;                       \
  y_limb_t b0,b1,b2,b3,b4,b5,b6,b7;                       \
  /* Cycle  1 */                                          \
  a0 = _t##5##r + _t##5##i;                               \
  b0 = _t##5##i - _t##5##r;                               \
  /* Cycle  2 */                                          \
  a1 = _t##7##r - _t##7##i;                               \
  b1 = _t##7##i + _t##7##r;                               \
  /* Cycle  3 */                                          \
  a2 = _t##0##r + _t##4##r;                               \
  b2 = _t##0##i + _t##4##i;                               \
  /* Cycle  4 */                                          \
  a3 = _t##0##r - _t##4##r;                               \
  b3 = _t##0##i - _t##4##i;                               \
  /* Cycle  5 */                                          \
  a4 = _t##2##r + _t##6##i;                               \
  b4 = _t##2##i - _t##6##r;                               \
  /* Cycle  6 */                                          \
  a5 = _t##2##r - _t##6##i;                               \
  b5 = _t##2##i + _t##6##r;                               \
  /* Cycle  7 */                                          \
  a6 = a0 * fr + _t##1##r;                                \
  b6 = b0 * fr + _t##1##i;                                \
  /* Cycle  8 */                                          \
  a0 = a0 * fi + _t##1##r;                                \
  b0 = b0 * fi + _t##1##i;                                \
  /* Cycle  9 */                                          \
  *(_q##0##r) = a2;                                       \
  *(_q##0##i) = b2;                                       \
  a7 = a1 * fi + _t##3##r;                                \
  /* Cycle 10 */                                          \
  *(_q##4##r) = a3;                                       \
  *(_q##4##i) = b3;                                       \
  b7 = b1 * fi + _t##3##i;                                \
  /* Cycle 11 */                                          \
  *(_q##2##r) = a4;                                       \
  *(_q##2##i) = b4;                                       \
  a1 = a1 * fr + _t##3##r;                                \
  /* Cycle  12 */                                         \
  *(_q##6##r) = a5;                                       \
  *(_q##6##i) = b5;                                       \
  b1 = b1 * fr + _t##3##i;                                \
  /* Cycle  13 */                                         \
  *(_q##1##r) = a6;                                       \
  *(_q##1##i) = b6;                                       \
  /* Cycle  14 */                                         \
  *(_q##5##r) = a0;                                       \
  *(_q##5##i) = b0;                                       \
  /* Cycle  15 */                                         \
  *(_q##3##r) = a7;                                       \
  *(_q##3##i) = b7;                                       \
  /* Cycle  16 */                                         \
  *(_q##7##r) = a1;                                       \
  *(_q##7##i) = b1;                                       \
}

                                                                                                                                                        /*
                                                                                                                                                           Macro to make the third radix-8 pass and store in a DIF 
                                                                                                                                                           _t is the locale prefix variable name
                                                                                                                                                           _q is the pointers prefix variable name
                                                                                                                                                           _i is the index
                                                                                                                                                        */
#define cplx_radix_8_pp_pass3_store_B(_t,_q)              \
{                                                         \
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;                       \
  y_limb_t b0,b1,b2,b3,b4,b5,b6,b7;                       \
  /* Cycle  1 */                                          \
  a0 = _t##5##r - _t##5##i;                               \
  b0 = _t##5##i + _t##5##r;                               \
  /* Cycle  2 */                                          \
  a1 = _t##7##r + _t##7##i;                               \
  b1 = _t##7##i - _t##7##r;                               \
  /* Cycle  3 */                                          \
  a2 = _t##0##r + _t##4##r;                               \
  b2 = _t##0##i + _t##4##i;                               \
  /* Cycle  4 */                                          \
  a3 = _t##0##r - _t##4##r;                               \
  b3 = _t##0##i - _t##4##i;                               \
  /* Cycle  5 */                                          \
  a4 = _t##2##r - _t##6##i;                               \
  b4 = _t##2##i + _t##6##r;                               \
  /* Cycle  6 */                                          \
  a5 = _t##2##r + _t##6##i;                               \
  b5 = _t##2##i - _t##6##r;                               \
  /* Cycle  7 */                                          \
  a6 = a0 * fr + _t##1##r;                                \
  b6 = b0 * fr + _t##1##i;                                \
  /* Cycle  8 */                                          \
  a0 = a0 * fi + _t##1##r;                                \
  b0 = b0 * fi + _t##1##i;                                \
  /* Cycle  9 */                                          \
  *(_q##0##r) = a2;                                       \
  *(_q##0##i) = b2;                                       \
  a7 = a1 * fi + _t##3##r;                                \
  /* Cycle 10 */                                          \
  *(_q##4##r) = a3;                                       \
  *(_q##4##i) = b3;                                       \
  b7 = b1 * fi + _t##3##i;                                \
  /* Cycle 11 */                                          \
  *(_q##2##r) = a4;                                       \
  *(_q##2##i) = b4;                                       \
  a1 = a1 * fr + _t##3##r;                                \
  /* Cycle  12 */                                         \
  *(_q##6##r) = a5;                                       \
  *(_q##6##i) = b5;                                       \
  b1 = b1 * fr + _t##3##i;                                \
  /* Cycle  13 */                                         \
  *(_q##1##r) = a6;                                       \
  *(_q##1##i) = b6;                                       \
  /* Cycle  14 */                                         \
  *(_q##5##r) = a0;                                       \
  *(_q##5##i) = b0;                                       \
  /* Cycle  15 */                                         \
  *(_q##3##r) = a7;                                       \
  *(_q##3##i) = b7;                                       \
  /* Cycle  16 */                                         \
  *(_q##7##r) = a1;                                       \
  *(_q##7##i) = b1;                                       \
}


#define cplx_radix_8_pp_last_dit_pass3(_t)                \
{                                                         \
  y_limb_t a0,a1,a2,a3,a4,a5;                             \
  y_limb_t b0,b1,b2,b3,b4,b5;                             \
  /* Cycle  1 */                                          \
  a0 = _t##5##r - _t##5##i;                               \
  b0 = _t##5##i + _t##5##r;                               \
  /* Cycle  2 */                                          \
  a1 = _t##7##r + _t##7##i;                               \
  b1 = _t##7##i - _t##7##r;                               \
  /* Cycle  3 */                                          \
  a2 = _t##0##r + _t##4##r;                               \
  b2 = _t##0##i + _t##4##i;                               \
  /* Cycle  4 */                                          \
  _t##4##r = _t##0##r - _t##4##r;                         \
  _t##4##i = _t##0##i - _t##4##i;                         \
  /* Cycle  5 */                                          \
  a3 = _t##2##r - _t##6##i;                               \
  b3 = _t##2##i + _t##6##r;                               \
  /* Cycle  6 */                                          \
  a4 = _t##2##r + _t##6##i;                               \
  b4 = _t##2##i - _t##6##r;                               \
  /* Cycle  7 */                                          \
  a5 = a0 * fr + _t##1##r;                                \
  b5 = b0 * fr + _t##1##i;                                \
  /* Cycle  8 */                                          \
  _t##5##r = a0 * fi + _t##1##r;                          \
  _t##5##i = b0 * fi + _t##1##i;                          \
  /* Cycle  9 */                                          \
  _t##0##r = a2;                                          \
  _t##0##i = b2;                                          \
  /* Cycle 10 */                                          \
  a0 = a1 * fi + _t##3##r;                                \
  b0 = b1 * fi + _t##3##i;                                \
  /* Cycle 11 */                                          \
  _t##2##r = a3;                                          \
  _t##2##i = b3;                                          \
  /* Cycle  12 */                                         \
  _t##6##r = a4;                                          \
  _t##6##i = b4;                                          \
  /* Cycle  13 */                                         \
  _t##7##r = a1 * fr + _t##3##r;                          \
  _t##7##i = b1 * fr + _t##3##i;                          \
  /* Cycle  14 */                                         \
  _t##1##r = a5;                                          \
  _t##1##i = b5;                                          \
  /* Cycle  15 */                                         \
  _t##3##r = a0;                                          \
  _t##3##i = b0;                                          \
}

#define cplx_radix_8_pp_first_dif_pass1(_t)               \
{                                                         \
  y_limb_t a0,a1,a2,a3;                                   \
  y_limb_t b0,b1,b2,b3;                                   \
  /* Cycle  1 */                                          \
  a0 = _t##0##r;                                          \
  b0 = _t##0##i;                                          \
  /* Cycle  2 */                                          \
  a2 = _t##2##r;                                          \
  b2 = _t##2##i;                                          \
  /* Cycle  3 */                                          \
  a1 = _t##1##r;                                          \
  b1 = _t##1##i;                                          \
  /* Cycle  4 */                                          \
  a3 = _t##3##r;                                          \
  b3 = _t##3##i;                                          \
  /* Cycle  5 */                                          \
  _t##0##r += t##4##r;                                    \
  _t##0##i += t##4##i;                                    \
  /* Cycle  6 */                                          \
  _t##4##r = a0 - t##4##r;                                \
  _t##4##i = b0 - t##4##i;                                \
  /* Cycle  7 */                                          \
  _t##2##r += t##6##r;                                    \
  _t##2##i += t##6##i;                                    \
  /* Cycle  8 */                                          \
  _t##6##r = a2 - t##6##r;                                \
  _t##6##i = b2 - t##6##i;                                \
  /* Cycle  9 */                                          \
  _t##1##r += t##5##r;                                    \
  _t##1##i += t##5##i;                                    \
  /* Cycle 10 */                                          \
  _t##5##r = a1 - t##5##r;                                \
  _t##5##i = b1 - t##5##i;                                \
  /* Cycle 11 */                                          \
  _t##3##r += t##7##r;                                    \
  _t##3##i += t##7##i;                                    \
  /* Cycle 10 */                                          \
  _t##7##r = a3 - t##7##r;                                \
  _t##7##i = b3 - t##7##i;                                \
}

#define cplx_radix_8_pp_first_dif_pass2(_t)               \
{                                                         \
  y_limb_t a0,a1,a2,a3;                                   \
  y_limb_t b0,b1,b2,b3;                                   \
  /* Cycle  1 */                                          \
  a0 = _t##2##r;                                          \
  b0 = _t##2##i;                                          \
  /* Cycle  2 */                                          \
  a1 = _t##6##r;                                          \
  b1 = _t##6##i;                                          \
  /* Cycle  3 */                                          \
  a2 = _t##3##r;                                          \
  b2 = _t##3##i;                                          \
  /* Cycle  4 */                                          \
  a3 = _t##7##r;                                          \
  b3 = _t##7##i;                                          \
  /* Cycle  5 */                                          \
  _t##2##r = _t##0##r - a0;                               \
  _t##2##i = _t##0##i - b0;                               \
  /* Cycle  6 */                                          \
  _t##0##r += a0;                                         \
  _t##0##i += b0;                                         \
  /* Cycle  7 */                                          \
  _t##6##r = _t##4##r - b1;                               \
  _t##6##i = _t##4##i + a1;                               \
  /* Cycle  8 */                                          \
  _t##4##r += b1;                                         \
  _t##4##i -= a1;                                         \
  /* Cycle  9 */                                          \
  _t##3##r = _t##1##r - a2;                               \
  _t##3##i = _t##1##i - b2;                               \
  /* Cycle 10 */                                          \
  _t##1##r += a2;                                         \
  _t##1##i += b2;                                         \
  /* Cycle 11 */                                          \
  _t##7##r = _t##5##r - b3;                               \
  _t##7##i = _t##5##i + a3;                               \
  /* Cycle 12 */                                          \
  _t##5##r += b3;                                         \
  _t##5##i -= a3;                                         \
}


                                                                                                                                                        /*
                                                                                                                                                           Macro to make the third radix-8 pass and store in a DIF 
                                                                                                                                                           _t is the locale prefix variable name
                                                                                                                                                           _q is the pointers prefix variable name
                                                                                                                                                           _i is the index
                                                                                                                                                        */
#define cplx_radix_8_pp_first_dif_pass3(_t,_q)            \
{                                                         \
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;                       \
  y_limb_t b0,b1,b2,b3,b4,b5,b6,b7;                       \
  /* Cycle  1 */                                          \
  a0 = _t##5##r + _t##5##i;                               \
  b0 = _t##5##i - _t##5##r;                               \
  /* Cycle  2 */                                          \
  a1 = _t##7##r - _t##7##i;                               \
  b1 = _t##7##i + _t##7##r;                               \
  /* Cycle  3 */                                          \
  a2 = _t##0##r + _t##1##r;                               \
  b2 = _t##0##i + _t##1##i;                               \
  /* Cycle  4 */                                          \
  a3 = _t##0##r - _t##1##r;                               \
  b3 = _t##0##i - _t##1##i;                               \
  /* Cycle  5 */                                          \
  a4 = _t##2##r + _t##3##i;                               \
  b4 = _t##2##i - _t##3##r;                               \
  /* Cycle  6 */                                          \
  a5 = _t##2##r - _t##3##i;                               \
  b5 = _t##2##i + _t##3##r;                               \
  /* Cycle  7 */                                          \
  a6 = a0 * fr + _t##4##r;                                \
  b6 = b0 * fr + _t##4##i;                                \
  /* Cycle  8 */                                          \
  a0 = a0 * fi + _t##4##r;                                \
  b0 = b0 * fi + _t##4##i;                                \
  /* Cycle  9 */                                          \
  *(_q##0##r) = a2;                                       \
  *(_q##0##i) = b2;                                       \
  a7 = a1 * fi + _t##6##r;                                \
  /* Cycle 10 */                                          \
  *(_q##4##r) = a3;                                       \
  *(_q##4##i) = b3;                                       \
  b7 = b1 * fi + _t##6##i;                                \
  /* Cycle 11 */                                          \
  *(_q##2##r) = a4;                                       \
  *(_q##2##i) = b4;                                       \
  a1 = a1 * fr + _t##6##r;                                \
  /* Cycle  12 */                                         \
  *(_q##6##r) = a5;                                       \
  *(_q##6##i) = b5;                                       \
  b1 = b1 * fr + _t##6##i;                                \
  /* Cycle  13 */                                         \
  *(_q##1##r) = a6;                                       \
  *(_q##1##i) = b6;                                       \
  /* Cycle  14 */                                         \
  *(_q##5##r) = a0;                                       \
  *(_q##5##i) = b0;                                       \
  /* Cycle  15 */                                         \
  *(_q##3##r) = a7;                                       \
  *(_q##3##i) = b7;                                       \
  /* Cycle  16 */                                         \
  *(_q##7##r) = a1;                                       \
  *(_q##7##i) = b1;                                       \
}



                                                                                                                                                        /*********************************************************************
                                                                                                                                                        This is the dyadic square for nested-compls representation
                                                                                                                                                         xk =(xk + xmk* ) * (xk + xmk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
                                                                                                                                                             G^(-k) * (xk - xmk* ) * (xk - xmk* );
                                                                                                                                                         xkm =(xkm + xk* ) * (xmk + xk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
                                                                                                                                                             G^(k) * (xmk - xk* ) * (xmk - xk* );
                                                                                                                                                         where all are pseudo-complex vars. 
                                                                                                                                                        **********************************************************************/
                                                                                                                                                        /* 36 fadds  24  Fmuls*/
#define square_nested_0_0_1_4(xk,xmk,yk,ymk,gkr,gki) \
{                                                    \
   y_limb_t _r0,_r1,_r2,_r3,_r4,_r5;                 \
   y_limb_t _s0,_s1,_s2,_s3,_s4,_s5;                 \
   /* Cycle  1 */                                    \
   _r4=(xk##r - xmk##r);                             \
   _r5=(xk##i + xmk##i);                             \
   /* Cycle  2 */                                    \
   _s4=(yk##r - ymk##r);                             \
   _s5=(yk##i + ymk##i);                             \
   /* Cycle  3 */                                    \
   _r2= xk##i * xk##i;                               \
   _s2= yk##i * yk##i;                               \
   /* Cycle  4 */                                    \
   xk##i *= xk##r;                                   \
   yk##i *= yk##r;                                   \
   /* Cycle  5 */                                    \
   _r0 = _r4 * 0.5;                                  \
   _r1 = _r5 * 0.5;                                  \
   /* Cycle  6 */                                    \
   _s0 = _s4 * 0.5;                                  \
   _s1 = _s5 * 0.5;                                  \
   /* Cycle  7 */                                    \
   _r4 = _r1 * _r1;                                  \
   _s4 = _s1 * _s1;                                  \
   /* Cycle  8 */                                    \
   xk##r = (xk##r * xk##r) - _r2;                    \
   yk##r = (yk##r * yk##r) - _s2;                    \
   /* Cycle  9 */                                    \
   _r1 *= _r0;                                       \
   _s1 *= _s0;                                       \
   /* Cycle 10 */                                    \
   _r5 = xmk##i * xmk##i ;                           \
   _s5 = ymk##i * ymk##i ;                           \
   /* Cycle 11 */                                    \
   xmk##i *= xmk##r;                                 \
   ymk##i *= ymk##r;                                 \
   /* Cycle 12 */                                    \
   _r1 +=  _r1;                                      \
   _s1 +=  _s1;                                      \
   /* Cycle 13 */                                    \
   _r2 = (_r0  * _r0) - _r4;                         \
   _s2 = (_s0  * _s0) - _s4;                         \
   /* Cycle 14 */                                    \
   _r4 = gkr + 1.0;                                  \
   _s4 = gki + 1.0;                                  \
   /* Cycle 15 */                                    \
   xmk##r =(xmk##r * xmk##r) - _r5;                  \
   ymk##r =(ymk##r * ymk##r) - _s5;                  \
   /* Cycle 16 */                                    \
   _r3 = gki * _r1;                                  \
   _s3 = gkr * _s1;                                  \
   /* Cycle 17 */                                    \
   _r5 = gki * _r2;                                  \
   _s5 = gkr * _s2;                                  \
   /* Cycle 18 */                                    \
   xk##i += xk##i;                                   \
   yk##i += yk##i;                                   \
   /* Cycle 19 */                                    \
   /* Cycle 20 */                                    \
   xmk##i += xmk##i;                                 \
   ymk##i += ymk##i;                                 \
   /* Cycle 21 */                                    \
   _r0= (_r4 * _r2) - _r3;                           \
   _s0= (_s4 * _s2) + _s3;                           \
   /* Cycle 22 */                                    \
   _r1= (_r4 * _r1) + _r5;                           \
   _s1= (_s4 * _s1) - _s5;                           \
   /* Cycle 23 */                                    \
   /* Cycle 24 */                                    \
   /* Cycle 25 */                                    \
   xk##r -= _r0;                                     \
   xmk##r -= _r0;                                    \
   /* Cycle 26 */                                    \
   yk##r -= _s0;                                     \
   ymk##r -= _s0;                                    \
   /* Cycle 27 */                                    \
   xk##i -= _r1;                                     \
   xmk##i += _r1;                                    \
   /* Cycle 28 */                                    \
   yk##i -= _s1;                                     \
   ymk##i += _s1;                                    \
}

#define square_nested_1_2_3_4(xk,xmk,yk,ymk,gkr,gki) \
{                                                    \
   y_limb_t _r0,_r1,_r2,_r3,_r4,_r5;                 \
   y_limb_t _s0,_s1,_s2,_s3,_s4,_s5;                 \
   /* Cycle  1 */                                    \
   _r4=(xk##r - xmk##r);                             \
   _r5=(xk##i + xmk##i);                             \
   /* Cycle  2 */                                    \
   _s4=(yk##r - ymk##r);                             \
   _s5=(yk##i + ymk##i);                             \
   /* Cycle  3 */                                    \
   _r2= xk##i * xk##i;                               \
   _s2= yk##i * yk##i;                               \
   /* Cycle  4 */                                    \
   xk##i *= xk##r;                                   \
   yk##i *= yk##r;                                   \
   /* Cycle  5 */                                    \
   _r0 = _r4 * 0.5;                                  \
   _r1 = _r5 * 0.5;                                  \
   /* Cycle  6 */                                    \
   _s0 = _s4 * 0.5;                                  \
   _s1 = _s5 * 0.5;                                  \
   /* Cycle  7 */                                    \
   _r4 = _r1 * _r1;                                  \
   _s4 = _s1 * _s1;                                  \
   /* Cycle  8 */                                    \
   xk##r = (xk##r * xk##r) - _r2;                    \
   yk##r = (yk##r * yk##r) - _s2;                    \
   /* Cycle  9 */                                    \
   _r1 *= _r0;                                       \
   _s1 *= _s0;                                       \
   /* Cycle 10 */                                    \
   _r5 = xmk##i * xmk##i ;                           \
   _s5 = ymk##i * ymk##i ;                           \
   /* Cycle 11 */                                    \
   xmk##i *= xmk##r;                                 \
   ymk##i *= ymk##r;                                 \
   /* Cycle 12 */                                    \
   _r1 +=  _r1;                                      \
   _s1 +=  _s1;                                      \
   /* Cycle 13 */                                    \
   _r2 = (_r0  * _r0) - _r4;                         \
   _s2 = (_s0  * _s0) - _s4;                         \
   /* Cycle 14 */                                    \
   _r4 = gkr - 1.0;                                  \
   _s4 = gki - 1.0;                                  \
   /* Cycle 15 */                                    \
   xmk##r =(xmk##r * xmk##r) - _r5;                  \
   ymk##r =(ymk##r * ymk##r) - _s5;                  \
   /* Cycle 16 */                                    \
   _r3 = gki * _r1;                                  \
   _s3 = gkr * _s1;                                  \
   /* Cycle 17 */                                    \
   _r5 = gki * _r2;                                  \
   _s5 = gkr * _s2;                                  \
   /* Cycle 18 */                                    \
   xk##i += xk##i;                                   \
   yk##i += yk##i;                                   \
   /* Cycle 19 */                                    \
   /* Cycle 20 */                                    \
   xmk##i += xmk##i;                                 \
   ymk##i += ymk##i;                                 \
   /* Cycle 21 */                                    \
   _r0= (_r4 * _r2) - _r3;                           \
   _s0= (_s4 * _s2) + _s3;                           \
   /* Cycle 22 */                                    \
   _r1= (_r4 * _r1) + _r5;                           \
   _s1= (_s4 * _s1) - _s5;                           \
   /* Cycle 23 */                                    \
   /* Cycle 24 */                                    \
   /* Cycle 25 */                                    \
   xk##r += _r0;                                     \
   xmk##r += _r0;                                    \
   /* Cycle 26 */                                    \
   yk##r += _s0;                                     \
   ymk##r += _s0;                                    \
   /* Cycle 27 */                                    \
   xk##i += _r1;                                     \
   xmk##i -= _r1;                                    \
   /* Cycle 28 */                                    \
   yk##i += _s1;                                     \
   ymk##i -= _s1;                                    \
}

#if !defined(Y_MAXIMUM)
                                                                                                                                                        /* This Macro is to init the loop in a last DIF first DIT pass
                                                                                                                                                         1) Preload trigs factors
                                                                                                                                                         2) Computes preload pointers
                                                                                                                                                         3) Preload data
                                                                                                                                                         _jj is the index of the upward block
                                                                                                                                                         _kk is the index of the backward block
                                                                                                                                                        */
# define last_dif_4_before_loop_preload(_jj,_kk)           \
{                                                         \
   y_limb_t a2,a3;                                        \
   y_limb_t b2,b3;                                        \
   y_size_t _j,_k;                                        \
   tpw1r = *(pxtr);                                       \
   tpw1i = *(pxti);                                       \
   dw1r = *(pxdr);                                        \
   dw1i = *(pxdi);                                        \
   _j= addr(_jj<<3);                                      \
   _k= addr(_kk<<3);                                      \
   a3 = tpw1i * tpw1r;                                    \
   b3 = dw1i * dw1r;                                      \
   a2 = tpw1i * tpw1i;                                    \
   b2 = dw1i * dw1i;                                      \
   tw2i = a3 + a3;                                        \
   dw2i = b3 + b3;                                        \
   tw2r = tpw1r * tpw1r - a2;                             \
   dw2r = dw1r * dw1r - b2;                               \
   a2 = tpw1i * tw2i;                                     \
   b2 = dw1i * dw2i;                                      \
   a3 = tpw1i * tw2r;                                     \
   b3 = dw1i * dw2r;                                      \
   tw3r = tpw1r * tw2r - a2;                              \
   dw3r = dw1r * dw2r - b2;                               \
   tw3i = tpw1r * tw2i + a3;                              \
   dw3i = dw1r * dw2i + b3;                               \
   ptl0r = p0r + _j;                                      \
   pdl0r = p0r + _k;                                      \
   ptl0i = p0i + _j;                                      \
   pdl0i = p0i + _k;                                      \
   ptl1r = p1r + _j;                                      \
   pdl1r = p1r + _k;                                      \
   ptl1i = p1i + _j;                                      \
   pdl1i = p1i + _k;                                      \
   ptl2r = p2r + _j;                                      \
   pdl2r = p2r + _k;                                      \
   ptl2i = p2i + _j;                                      \
   pdl2i = p2i + _k;                                      \
   ptl3r = p3r + _j;                                      \
   pdl3r = p3r + _k;                                      \
   ptl3i = p3i + _j;                                      \
   pdl3i = p3i + _k;                                      \
   tp0r = *(ptl0r);                                       \
   sp0r = *(pdl0r);                                       \
   tp0i = *(ptl0i);                                       \
   sp0i = *(pdl0i);                                       \
   tp1r = *(ptl1r);                                       \
   sp1r = *(pdl1r);                                       \
   tp1i = *(ptl1i);                                       \
   sp1i = *(pdl1i);                                       \
   tp2r = *(ptl2r);                                       \
   sp2r = *(pdl2r);                                       \
   tp2i = *(ptl2i);                                       \
   sp2i = *(pdl2i);                                       \
   tp3r = *(ptl3r);                                       \
   sp3r = *(pdl3r);                                       \
   tp3i = *(ptl3i);                                       \
   sp3i = *(pdl3i);                                       \
}
#else

                                                                                                                                                        /*Version for Y_MAXIMUM */
# define last_dif_4_before_loop_preload(_jj,_kk)           \
{                                                         \
   y_size_t _j,_k;                                        \
   tpw1r = *(pxtr);                                       \
   tpw1i = *(pxti);                                       \
   dw1r = *(pxdr);                                        \
   dw1i = *(pxdi);                                        \
   tw2r = *(pxtr+2);                                      \
   tw2i = *(pxti+2);                                      \
   dw2r = *(pxdr+2);                                      \
   dw2i = *(pxdi+2);                                      \
   tw3r = *(pxtr+4);                                      \
   tw3i = *(pxti+4);                                      \
   dw3r = *(pxdr+4);                                      \
   dw3i = *(pxdi+4);                                      \
   _j= addr(_jj<<3);                                      \
   _k= addr(_kk<<3);                                      \
   ptl0r = p0r + _j;                                      \
   pdl0r = p0r + _k;                                      \
   ptl0i = p0i + _j;                                      \
   pdl0i = p0i + _k;                                      \
   ptl1r = p1r + _j;                                      \
   pdl1r = p1r + _k;                                      \
   ptl1i = p1i + _j;                                      \
   pdl1i = p1i + _k;                                      \
   ptl2r = p2r + _j;                                      \
   pdl2r = p2r + _k;                                      \
   ptl2i = p2i + _j;                                      \
   pdl2i = p2i + _k;                                      \
   ptl3r = p3r + _j;                                      \
   pdl3r = p3r + _k;                                      \
   ptl3i = p3i + _j;                                      \
   pdl3i = p3i + _k;                                      \
   tp0r = *(ptl0r);                                       \
   sp0r = *(pdl0r);                                       \
   tp0i = *(ptl0i);                                       \
   sp0i = *(pdl0i);                                       \
   tp1r = *(ptl1r);                                       \
   sp1r = *(pdl1r);                                       \
   tp1i = *(ptl1i);                                       \
   sp1i = *(pdl1i);                                       \
   tp2r = *(ptl2r);                                       \
   sp2r = *(pdl2r);                                       \
   tp2i = *(ptl2i);                                       \
   sp2i = *(pdl2i);                                       \
   tp3r = *(ptl3r);                                       \
   sp3r = *(pdl3r);                                       \
   tp3i = *(ptl3i);                                       \
   sp3i = *(pdl3i);                                       \
}
#endif

                                                                                                                                                        /* This macro make all last dif-4 work with preloaded data
                                                                                                                                                           It makes both upward and downward block transform */
#define big_macro_last_dif_4_preload(_jj,_kk)             \
{                                                         \
   y_limb_t a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5;          \
   y_size_t _j,_k;                                        \
   _j= addr(((_jj + 1)<<3));                              \
   _k= addr(((_kk - 1)<<3));                              \
   tw1r = tpw1r;                                          \
   tw1i = tpw1i;                                          \
   /* cycle 1 */                                          \
   pt0r = ptl0r;                                          \
   a0 = tpw1i * tp1i;                                     \
   pd0r = pdl0r;                                          \
   b0 = dw1i * sp1i;                                      \
   /* cycle 2 */                                          \
   pt0i = ptl0i;                                          \
   a1 = tpw1i * tp1r;                                     \
   pd0i = pdl0i;                                          \
   b1 = dw1i * sp1r;                                      \
   /* cycle 3 */                                          \
   pt1r = ptl1r;                                          \
   a2 = tw2i * tp2i;                                      \
   pd1r = pdl1r;                                          \
   b2 = dw2i * sp2i;                                      \
   /* cycle 4 */                                          \
   pt1i = ptl1i;                                          \
   a3 = tw2i * tp2r;                                      \
   pd1i = pdl1i;                                          \
   b3 = dw2i * sp2r;                                      \
   /* cycle 5 */                                          \
   pt2r = ptl2r;                                          \
   a4 = tw3i * tp3i;                                      \
   pd2r = pdl2r;                                          \
   b4 = dw3i * sp3i;                                      \
   /* cycle 6 */                                          \
   pt2i = ptl2i;                                          \
   a5 = tw3i * tp3r;                                      \
   pd2i = pdl2i;                                          \
   b5 = dw3i * sp3r;                                      \
   /* cycle 7 */                                          \
   pxtr += Y_STEP;                                        \
   tp1r = tpw1r * tp1r - a0;                              \
   pxti += Y_STEP;                                        \
   sp1r = dw1r * sp1r - b0;                               \
   /* cycle 8 */                                          \
   pxdr -= Y_STEP;                                        \
   tp1i = tpw1r * tp1i + a1;                              \
   pxdi -= Y_STEP;                                        \
   sp1i = dw1r * sp1i + b1;                               \
   /* cycle 7 */                                          \
   pt3r = ptl3r;                                          \
   tp2r = tw2r * tp2r - a2;                               \
   pd3r = pdl3r;                                          \
   sp2r = dw2r * sp2r - b2;                               \
   /* cycle 8 */                                          \
   pt3i = ptl3i;                                          \
   tp2i = tw2r * tp2i + a3;                               \
   pd3i = pdl3i;                                          \
   sp2i = dw2r * sp2i + b3;                               \
   /* cycle 9 */                                          \
   ptl0r = p0r + _j;                                      \
   tp3r = tw3r * tp3r - a4;                               \
   pdl0r = p0r + _k;                                      \
   sp3r = dw3r * sp3r - b4;                               \
   /* cycle 10 */                                         \
   ptl0i = p0i + _j;                                      \
   tp3i = tw3r * tp3i + a5;                               \
   pdl0i = p0i + _k;                                      \
   sp3i = dw3r * sp3i + b5;                               \
   /* Cycle 11 */                                         \
   /* Cycle 12 */                                         \
   t0r = tp0r + tp2r;                                     \
   ptl1r = p1r + _j;                                      \
   u0r = sp0r + sp2r;                                     \
   pdl1r = p1r + _k;                                      \
   /* cycle 11 */                                         \
   t1r = tp0r - tp2r;                                     \
   ptl1i = p1i + _j;                                      \
   u1r = sp0r - sp2r;                                     \
   pdl1i = p1i + _k;                                      \
   /* cycle 12 */                                         \
   t0i = tp0i + tp2i;                                     \
   ptl2r = p2r + _j;                                      \
   u0i = sp0i + sp2i;                                     \
   pdl2r = p2r + _k;                                      \
    /* cycle 13 */                                        \
   tpw1i = *(pxti);                                       \
   t1i = tp0i - tp2i;                                     \
   ptl2i = p2i + _j;                                      \
   dw1i = *(pxdi);                                        \
   u1i = sp0i - sp2i;                                     \
   pdl2i = p2i + _k;                                      \
   /* cycle 14 */                                         \
   tp0r = *(ptl0r);                                       \
   t2r = tp1r + tp3r;                                     \
   ptl3r = p3r + _j;                                      \
   sp0r = *(pdl0r);                                       \
   u2r = sp1r + sp3r;                                     \
   pdl3r = p3r + _k;                                      \
   /* Cycle 15 */                                         \
   tp0i = *(ptl0i);                                       \
   t3r = tp1r - tp3r;                                     \
   ptl3i = p3i + _j;                                      \
   sp0i = *(pdl0i);                                       \
   u3r = sp1r - sp3r;                                     \
   pdl3i = p3i + _k;                                      \
   /* Cycle 16 */                                         \
   tp1r = *(ptl1r);                                       \
   t2i = tp1i + tp3i;                                     \
   sp1r = *(pdl1r);                                       \
   u2i = sp1i + sp3i;                                     \
   /* cycle 17 */                                         \
   tpw1r = *(pxtr);                                       \
   t3i = tp1i - tp3i;                                     \
   dw1r = *(pxdr);                                        \
   u3i = sp1i - sp3i;                                     \
   /* cycle 18 */                                         \
   tp1i = *(ptl1i);                                       \
   a0 = t0r;                                              \
   sp1i = *(pdl1i);                                       \
   b0 = u0r;                                              \
   /* cycle 19 */                                         \
   tp2r = *(ptl2r);                                       \
   a1 = t0i;                                              \
   sp2r = *(pdl2r);                                       \
   b1 = u0i;                                              \
   /* cycle 20 */                                         \
   tp2i = *(ptl2i);                                       \
   t0r += t2r;                                            \
   sp2i = *(pdl2i);                                       \
   u0r += u2r;                                            \
   /* cycle 21 */                                         \
   tp3r = *(ptl3r);                                       \
   t2r = a0 - t2r;                                        \
   sp3r = *(pdl3r);                                       \
   u2r = b0 - u2r;                                        \
   /* cycle 22 */                                         \
   tp3i = *(ptl3i);                                       \
   t0i += t2i;                                            \
   sp3i = *(pdl3i);                                       \
   u0i += u2i;                                            \
   /* cycle 23 */                                         \
   t2i = a1 - t2i;                                        \
   u2i = b1 - u2i;                                        \
   /* cycle 24 */                                         \
   a0 = t3i;                                              \
   b0 = u3i;                                              \
   /* cycle 25 */                                         \
   a1 = t3r;                                              \
   b1 = u3r;                                              \
   /* cycle 26 */                                         \
   t3r = t1r - a0;                                        \
   u3r = u1r - b0;                                        \
   /* cycle 27 */                                         \
   t1r += a0;                                             \
   u1r += b0;                                             \
   /* cycle 28 */                                         \
   t3i = t1i + a1;                                        \
   u3i = u1i + b1;                                        \
   /* cycle 29 */                                         \
   t1i -= a1;                                             \
   u1i -= b1;                                             \
}

                                                                                                                                                        /*
                                                                                                                                                        This macro make both downward and upward first DIT-4
                                                                                                                                                        block, and store the results in memory. 
                                                                                                                                                        */
#define big_macro_first_dit_4_preload()                   \
{                                                         \
   y_limb_t a0,a1,a2,a3;                                  \
   y_limb_t b0,b1,b2,b3;                                  \
   /* cycle 1  */                                         \
   a3 = tpw1i * tpw1r;                                    \
   b3 = dw1i * dw1r;                                      \
   /* cycle 2  */                                         \
   a2 = tpw1i * tpw1i;                                    \
   b2 = dw1i * dw1i;                                      \
   /* cycle 3  */                                         \
   a0 = t0r;                                              \
   b0 = u0r;                                              \
   /* cycle 4  */                                         \
   a1 = t0i;                                              \
   b1 = u0i;                                              \
   /* cycle 5  */                                         \
   t0r += t2r;                                            \
   u0r += u2r;                                            \
   /* cycle 6  */                                         \
   t2r = a0 - t2r;                                        \
   u2r = b0 - u2r;                                        \
   /* cycle 7  */                                         \
   tw2i = a3 + a3;                                        \
   dw2i = b3 + b3;                                        \
   /* cycle 8  */                                         \
   tw2r = tpw1r * tpw1r - a2;                             \
   dw2r = dw1r * dw1r - b2;                               \
   /* cycle 9  */                                         \
   t0i += t2i;                                            \
   u0i += u2i;                                            \
   /* cycle 10 */                                         \
   t2i = a1 - t2i;                                        \
   u2i = b1 - u2i;                                        \
   /* cycle 11 */                                         \
   a0 = t1r;                                              \
   b0 = u1r;                                              \
   /* cycle 12 */                                         \
   a1 = t1i;                                              \
   b1 = u1i;                                              \
   /* cycle 13 */                                         \
   a2 = tpw1i * tw2i;                                     \
   b2 = dw1i * dw2i;                                      \
   /* cycle 14 */                                         \
   a3 = tpw1i * tw2r;                                     \
   b3 = dw1i * dw2r;                                      \
   /* cycle 15 */                                         \
   t1r += t3r;                                            \
   u1r += u3r;                                            \
   /* cycle 16 */                                         \
   t3r = a0 - t3r;                                        \
   u3r = b0 - u3r;                                        \
   /* cycle 17 */                                         \
   t1i += t3i;                                            \
   u1i += u3i;                                            \
   /* cycle 18 */                                         \
   t3i = a1 - t3i;                                        \
   u3i = b1 - u3i;                                        \
   /* cycle 19 */                                         \
   tw3r = tpw1r * tw2r - a2;                              \
   dw3r = dw1r * dw2r - b2;                               \
   /* cycle 20 */                                         \
   tw3i = tpw1r * tw2i + a3;                              \
   dw3i = dw1r * dw2i + b3;                               \
   /* cycle 21 */                                         \
   a0 = t0r + t1r;                                        \
   b0 = u0r + u1r;                                        \
   /* cycle 22 */                                         \
   a1 = t0r - t1r;                                        \
   b1 = u0r - u1r;                                        \
   /* cycle 23 */                                         \
   a2 = t0i + t1i;                                        \
   b2 = u0i + u1i;                                        \
   /* cycle 24 */                                         \
   *(pt0r) = a0;                                          \
   a3 = t0i - t1i;                                        \
   *(pd0r) = b0;                                          \
   b3 = u0i - u1i;                                        \
   /* cycle 25 */                                         \
   *(pt2r) = a1;                                          \
   a0 = t2r - t3i;                                        \
   *(pd2r) = b1;                                          \
   b0 = u2r - u3i;                                        \
   /* cycle 26 */                                         \
   *(pt0i) = a2;                                          \
   a1 = t2r + t3i;                                        \
   *(pd0i) = b2;                                          \
   b1 = u2r + u3i;                                        \
   /* cycle 27 */                                         \
   *(pt2i) = a3;                                          \
   a2 = t2i + t3r;                                        \
   *(pd2i) = b3;                                          \
   b2 = u2i + u3r;                                        \
   /* cycle 28 */                                         \
   *(pt1r) = a0;                                          \
   a3 = t2i - t3r;                                        \
   *(pd1r) = b0;                                          \
   b3 = u2i - u3r;                                        \
   /* cycle 29 */                                         \
   *(pt3r) = a1;                                          \
   *(pd3r) = b1;                                          \
   /* cycle 30 */                                         \
   *(pt1i) = a2;                                          \
   *(pd1i) = b2;                                          \
   /* cycle 31 */                                         \
   *(pt3i) = a3;                                          \
   *(pd3i) = b3;                                          \
}

                                                                                                                                                        /*
                                                                                                                                                         Radix_8 FFT reduction. Forward (DIF)                              
                                                                                                                                                          22 fadd, 22 fsub, 4 fused muladd                      
                                                                                                                                                          4 register changes, 8 aux registers                   
                                                                                                                                                         _i are scrambled indexes                              
                                                                                                                                                         _t is the prefix name of data                          */
#define cplx_fft8_ia64_F(_t,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   a0 = _t##_i0##r + _t##_i1##r;                          \
   b0 = _t##_i0##i + _t##_i1##i;                          \
  /* Cycle  2 */                                          \
   _t##_i1##r = _t##_i0##r - _t##_i1##r;                  \
   _t##_i1##i = _t##_i0##i - _t##_i1##i;                  \
  /* Cycle  3 */                                          \
   a1 = _t##_i2##r + _t##_i3##r;                          \
   b1 = _t##_i2##i + _t##_i3##i;                          \
  /* Cycle  4 */                                          \
   _t##_i3##r = _t##_i2##r - _t##_i3##r;                  \
   _t##_i3##i = _t##_i2##i - _t##_i3##i;                  \
  /* Cycle  5 */                                          \
   a2 = _t##_i4##r + _t##_i5##r;                          \
   b2 = _t##_i4##i + _t##_i5##i;                          \
  /* Cycle  6 */                                          \
   _t##_i5##r = _t##_i4##r - _t##_i5##r;                  \
   _t##_i5##i = _t##_i4##i - _t##_i5##i;                  \
  /* Cycle  7 */                                          \
   a3 = _t##_i6##r + _t##_i7##r;                          \
   b3 = _t##_i6##i + _t##_i7##i;                          \
  /* Cycle  8 */                                          \
  _t##_i7##r = _t##_i6##r - _t##_i7##r;                   \
  _t##_i7##i = _t##_i6##i - _t##_i7##i;                   \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _t##_i0##r = a0 + a1;                                   \
  _t##_i0##i = b0 + b1;                                   \
  /* Cycle 10 */                                          \
  _t##_i2##r = a0 - a1;                                   \
  _t##_i2##i = b0 - b1;                                   \
  /* Cycle 11 */                                          \
  a0 = _t##_i1##r - _t##_i3##i;/* t3 */                   \
  b0 = _t##_i1##i + _t##_i3##r;/* t3 */                   \
  /* Cycle 12 */                                          \
  _t##_i4##r = a2 + a3;                                   \
  _t##_i4##i = b2 + b3;                                   \
  /* Cycle 13 */                                          \
  _t##_i6##r = a2 - a3;                                   \
  _t##_i6##i = b2 - b3;                                   \
  /* Cycle 14 */                                          \
  a2 = _t##_i5##r - _t##_i7##i;/* t7 */                   \
  b2 = _t##_i5##i + _t##_i7##r;/* t7 */                   \
  /* Cycle 15 */                                          \
  _t##_i5##r +=  _t##_i7##i;                              \
  _t##_i5##i -=  _t##_i7##r;                              \
  /* Cycle 16 */                                          \
  _t##_i1##r +=  _t##_i3##i;                              \
  _t##_i1##i -=  _t##_i3##r;                              \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  a3 = _t##_i0##r + _t##_i4##r;/* t0 */                   \
  b3 = _t##_i0##i + _t##_i4##i;/* t0 */                   \
 /* Cycle  18 */                                          \
  _t##_i4##r = _t##_i0##r - _t##_i4##r;                   \
  _t##_i4##i = _t##_i0##i - _t##_i4##i;                   \
  /* Cycle 19 */                                          \
  _t##_i7##r = a2 - b2;                                   \
  _t##_i7##i = b2 + a2;                                   \
  /* Cycle 20 */                                          \
  a2 = _t##_i5##r + _t##_i5##i;                           \
  b2 = _t##_i5##i - _t##_i5##r;                           \
  /* Cycle 21 */                                          \
  a1 = _t##_i2##r - _t##_i6##i;/* t6 */                   \
  b1 = _t##_i2##i + _t##_i6##r;/* t6 */                   \
  /* Cycle 22 */                                          \
  _t##_i2##r += _t##_i6##i;                               \
  _t##_i2##i -= _t##_i6##r;                               \
  /* Cycle 23 */                                          \
  _t##_i0##r = a3;                                        \
  _t##_i0##i = b3;                                        \
  /* Cycle 24 */                                          \
  _t##_i3##r = _t##_i7##r * fi + a0;                      \
  _t##_i3##i = _t##_i7##i * fi + b0;                      \
  /* Cycle 25 */                                          \
  _t##_i7##r = _t##_i7##r * fr + a0;                      \
  _t##_i7##i = _t##_i7##i * fr + b0;                      \
  /* Cycle 26 */                                          \
  _t##_i5##r = a2 * fi + _t##_i1##r;                      \
  _t##_i5##i = b2 * fi + _t##_i1##i;                      \
  /* Cycle 27 */                                          \
  _t##_i1##r = a2 * fr + _t##_i1##r;                      \
  _t##_i1##i = b2 * fr + _t##_i1##i;                      \
  /* Cycle 28 */                                          \
  _t##_i6##r = a1;                                        \
  _t##_i6##i = b1;                                        \
}



                                                                                                                                                        /*
                                                                                                                                                         Radix_8 FFT reduction. Backward (DIT)                             
                                                                                                                                                         22 fadd, 22 fsub, 4 fused muladd                      
                                                                                                                                                         4 register changes, 8 aux registers                   
                                                                                                                                                         _i are scrambled indexes                              
                                                                                                                                                         _t is the prefix name of data
                                                                                                                                                        */
#define cplx_fft8_ia64_B(_t,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   a0 = _t##_i0##r + _t##_i1##r;                          \
   b0 = _t##_i0##i + _t##_i1##i;                          \
  /* Cycle  2 */                                          \
   _t##_i1##r = _t##_i0##r - _t##_i1##r;                  \
   _t##_i1##i = _t##_i0##i - _t##_i1##i;                  \
  /* Cycle  3 */                                          \
   a1 = _t##_i2##r + _t##_i3##r;                          \
   b1 = _t##_i2##i + _t##_i3##i;                          \
  /* Cycle  4 */                                          \
   _t##_i3##r = _t##_i2##r - _t##_i3##r;                  \
   _t##_i3##i = _t##_i2##i - _t##_i3##i;                  \
  /* Cycle  5 */                                          \
   a2 = _t##_i4##r + _t##_i5##r;                          \
   b2 = _t##_i4##i + _t##_i5##i;                          \
  /* Cycle  6 */                                          \
   _t##_i5##r = _t##_i4##r - _t##_i5##r;                  \
   _t##_i5##i = _t##_i4##i - _t##_i5##i;                  \
  /* Cycle  7 */                                          \
   a3 = _t##_i6##r + _t##_i7##r;                          \
   b3 = _t##_i6##i + _t##_i7##i;                          \
  /* Cycle  8 */                                          \
  _t##_i7##r = _t##_i6##r - _t##_i7##r;                   \
  _t##_i7##i = _t##_i6##i - _t##_i7##i;                   \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _t##_i0##r = a0 + a1;                                   \
  _t##_i0##i = b0 + b1;                                   \
  /* Cycle 10 */                                          \
  _t##_i2##r = a0 - a1;                                   \
  _t##_i2##i = b0 - b1;                                   \
  /* Cycle 11 */                                          \
  a0 = _t##_i1##r + _t##_i3##i;/* t3 */                   \
  b0 = _t##_i1##i - _t##_i3##r;/* t3 */                   \
  /* Cycle 12 */                                          \
  _t##_i4##r = a2 + a3;                                   \
  _t##_i4##i = b2 + b3;                                   \
  /* Cycle 13 */                                          \
  _t##_i6##r = a2 - a3;                                   \
  _t##_i6##i = b2 - b3;                                   \
  /* Cycle 14 */                                          \
  a2 = _t##_i5##r + _t##_i7##i;/* t7 */                   \
  b2 = _t##_i5##i - _t##_i7##r;/* t7 */                   \
  /* Cycle 15 */                                          \
  _t##_i5##r -=  _t##_i7##i;                              \
  _t##_i5##i +=  _t##_i7##r;                              \
  /* Cycle 16 */                                          \
  _t##_i1##r -=  _t##_i3##i;                              \
  _t##_i1##i +=  _t##_i3##r;                              \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  a3 = _t##_i0##r + _t##_i4##r;/* t0 */                   \
  b3 = _t##_i0##i + _t##_i4##i;/* t0 */                   \
 /* Cycle  18 */                                          \
  _t##_i4##r = _t##_i0##r - _t##_i4##r;                   \
  _t##_i4##i = _t##_i0##i - _t##_i4##i;                   \
  /* Cycle 19 */                                          \
  _t##_i7##r = a2 + b2;                                   \
  _t##_i7##i = b2 - a2;                                   \
  /* Cycle 20 */                                          \
  a2 = _t##_i5##r - _t##_i5##i;                           \
  b2 = _t##_i5##i + _t##_i5##r;                           \
  /* Cycle 21 */                                          \
  a1 = _t##_i2##r + _t##_i6##i;/* t6 */                   \
  b1 = _t##_i2##i - _t##_i6##r;/* t6 */                   \
  /* Cycle 22 */                                          \
  _t##_i2##r -= _t##_i6##i;                               \
  _t##_i2##i += _t##_i6##r;                               \
  /* Cycle 23 */                                          \
  _t##_i0##r = a3;                                        \
  _t##_i0##i = b3;                                        \
  /* Cycle 24 */                                          \
  _t##_i3##r = _t##_i7##r * fi + a0;                      \
  _t##_i3##i = _t##_i7##i * fi + b0;                      \
  /* Cycle 25 */                                          \
  _t##_i7##r = _t##_i7##r * fr + a0;                      \
  _t##_i7##i = _t##_i7##i * fr + b0;                      \
  /* Cycle 26 */                                          \
  _t##_i5##r = a2 * fi + _t##_i1##r;                      \
  _t##_i5##i = b2 * fi + _t##_i1##i;                      \
  /* Cycle 27 */                                          \
  _t##_i1##r = a2 * fr + _t##_i1##r;                      \
  _t##_i1##i = b2 * fr + _t##_i1##i;                      \
  /* Cycle 28 */                                          \
  _t##_i6##r = a1;                                        \
  _t##_i6##i = b1;                                        \
}



                                                                                                                                                        /*
                                                                                                                                                         Radix_8 FFT reduction. Forward (DIF)                              
                                                                                                                                                         22 fadd, 22 fsub, 4 fused muladd                      
                                                                                                                                                         4 register changes, 8 aux registers                   
                                                                                                                                                         _i are scrambled indexes                              
                                                                                                                                                         _t is the prefix name of data 
                                                                                                                                                         
                                                                                                                                                        The results are stored in memory:    *(_pdj) = _t##_i(j);
                                                                                                                                                         
                                                                                                                                                        */
#define cplx_fft8_ia64_store_F(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_pd7,_t,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   a0 = _t##_i0##r + _t##_i1##r;                          \
   b0 = _t##_i0##i + _t##_i1##i;                          \
  /* Cycle  2 */                                          \
   _t##_i1##r = _t##_i0##r - _t##_i1##r;                  \
   _t##_i1##i = _t##_i0##i - _t##_i1##i;                  \
  /* Cycle  3 */                                          \
   a1 = _t##_i2##r + _t##_i3##r;                          \
   b1 = _t##_i2##i + _t##_i3##i;                          \
  /* Cycle  4 */                                          \
   _t##_i3##r = _t##_i2##r - _t##_i3##r;                  \
   _t##_i3##i = _t##_i2##i - _t##_i3##i;                  \
  /* Cycle  5 */                                          \
   a2 = _t##_i4##r + _t##_i5##r;                          \
   b2 = _t##_i4##i + _t##_i5##i;                          \
  /* Cycle  6 */                                          \
   _t##_i5##r = _t##_i4##r - _t##_i5##r;                  \
   _t##_i5##i = _t##_i4##i - _t##_i5##i;                  \
  /* Cycle  7 */                                          \
   a3 = _t##_i6##r + _t##_i7##r;                          \
   b3 = _t##_i6##i + _t##_i7##i;                          \
  /* Cycle  8 */                                          \
  _t##_i7##r = _t##_i6##r - _t##_i7##r;                   \
  _t##_i7##i = _t##_i6##i - _t##_i7##i;                   \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _t##_i0##r = a0 + a1;                                   \
  _t##_i0##i = b0 + b1;                                   \
  /* Cycle 10 */                                          \
  _t##_i2##r = a0 - a1;                                   \
  _t##_i2##i = b0 - b1;                                   \
  /* Cycle 11 */                                          \
  a0 = _t##_i1##r - _t##_i3##i;/* t3 */                   \
  b0 = _t##_i1##i + _t##_i3##r;/* t3 */                   \
  /* Cycle 12 */                                          \
  _t##_i4##r = a2 + a3;                                   \
  _t##_i4##i = b2 + b3;                                   \
  /* Cycle 13 */                                          \
  _t##_i6##r = a2 - a3;                                   \
  _t##_i6##i = b2 - b3;                                   \
  /* Cycle 14 */                                          \
  a2 = _t##_i5##r - _t##_i7##i;/* t7 */                   \
  b2 = _t##_i5##i + _t##_i7##r;/* t7 */                   \
  /* Cycle 15 */                                          \
  _t##_i5##r +=  _t##_i7##i;                              \
  _t##_i5##i -=  _t##_i7##r;                              \
  /* Cycle 16 */                                          \
  _t##_i1##r +=  _t##_i3##i;                              \
  _t##_i1##i -=  _t##_i3##r;                              \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  a3 = _t##_i0##r + _t##_i4##r;/* t0 */                   \
  b3 = _t##_i0##i + _t##_i4##i;/* t0 */                   \
 /* Cycle  18 */                                          \
  _t##_i4##r = _t##_i0##r - _t##_i4##r;                   \
  _t##_i4##i = _t##_i0##i - _t##_i4##i;                   \
  /* Cycle 19 */                                          \
  _t##_i7##r = a2 - b2;                                   \
  _t##_i7##i = b2 + a2;                                   \
  /* Cycle 20 */                                          \
  a2 = _t##_i5##r + _t##_i5##i;                           \
  b2 = _t##_i5##i - _t##_i5##r;                           \
  /* Cycle 21 */                                          \
  a1 = _t##_i2##r - _t##_i6##i;/* t6 */                   \
  b1 = _t##_i2##i + _t##_i6##r;/* t6 */                   \
  /* Cycle 22 */                                          \
  _t##_i2##r += _t##_i6##i;                               \
  _t##_i2##i -= _t##_i6##r;                               \
  /* Cycle 23 */                                          \
  *(_pd0##r) = a3;                                        \
  *(_pd0##i) = b3;                                        \
  /* Cycle 24 */                                          \
  *(_pd4##r) = _t##_i4##r;                                \
  _t##_i3##r = _t##_i7##r * fi + a0;                      \
  *(_pd4##i) = _t##_i4##i;                                \
  _t##_i3##i = _t##_i7##i * fi + b0;                      \
  /* Cycle 25 */                                          \
  _t##_i7##r = _t##_i7##r * fr + a0;                      \
  _t##_i7##i = _t##_i7##i * fr + b0;                      \
  /* Cycle 26 */                                          \
  *(_pd2##r) = _t##_i2##r;                                \
  _t##_i5##r = a2 * fi + _t##_i1##r;                      \
  *(_pd2##i) = _t##_i2##i;                                \
  _t##_i5##i = b2 * fi + _t##_i1##i;                      \
  /* Cycle 27 */                                          \
  *( _pd6##r ) = a1;                                      \
  _t##_i1##r = a2 * fr + _t##_i1##r;                      \
  *( _pd6##i ) = b1;                                      \
  _t##_i1##i = b2 * fr + _t##_i1##i;                      \
  /* Cycle 28 */                                          \
  /* Cycle 29 */                                          \
  *(_pd3##r) = _t##_i3##r;                                \
  *(_pd3##i) = _t##_i3##i;                                \
  /* Cycle 30 */                                          \
  *(_pd7##r) = _t##_i7##r;                                \
  *(_pd7##i) = _t##_i7##i;                                \
  /* Cycle 31 */                                          \
  *(_pd5##r) = _t##_i5##r;                                \
  *(_pd5##i) = _t##_i5##i;                                \
  /* Cycle 32 */                                          \
  *(_pd1##r) = _t##_i1##r;                                \
  *(_pd1##i) = _t##_i1##i;                                \
}

                                                                                                                                                        /* This is for backward DIT */
#define cplx_fft8_ia64_store_B(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_pd7,_t,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   a0 = _t##_i0##r + _t##_i1##r;                          \
   b0 = _t##_i0##i + _t##_i1##i;                          \
  /* Cycle  2 */                                          \
   _t##_i1##r = _t##_i0##r - _t##_i1##r;                  \
   _t##_i1##i = _t##_i0##i - _t##_i1##i;                  \
  /* Cycle  3 */                                          \
   a1 = _t##_i2##r + _t##_i3##r;                          \
   b1 = _t##_i2##i + _t##_i3##i;                          \
  /* Cycle  4 */                                          \
   _t##_i3##r = _t##_i2##r - _t##_i3##r;                  \
   _t##_i3##i = _t##_i2##i - _t##_i3##i;                  \
  /* Cycle  5 */                                          \
   a2 = _t##_i4##r + _t##_i5##r;                          \
   b2 = _t##_i4##i + _t##_i5##i;                          \
  /* Cycle  6 */                                          \
   _t##_i5##r = _t##_i4##r - _t##_i5##r;                  \
   _t##_i5##i = _t##_i4##i - _t##_i5##i;                  \
  /* Cycle  7 */                                          \
   a3 = _t##_i6##r + _t##_i7##r;                          \
   b3 = _t##_i6##i + _t##_i7##i;                          \
  /* Cycle  8 */                                          \
  _t##_i7##r = _t##_i6##r - _t##_i7##r;                   \
  _t##_i7##i = _t##_i6##i - _t##_i7##i;                   \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _t##_i0##r = a0 + a1;                                   \
  _t##_i0##i = b0 + b1;                                   \
  /* Cycle 10 */                                          \
  _t##_i2##r = a0 - a1;                                   \
  _t##_i2##i = b0 - b1;                                   \
  /* Cycle 11 */                                          \
  a0 = _t##_i1##r + _t##_i3##i;/* t3 */                   \
  b0 = _t##_i1##i - _t##_i3##r;/* t3 */                   \
  /* Cycle 12 */                                          \
  _t##_i4##r = a2 + a3;                                   \
  _t##_i4##i = b2 + b3;                                   \
  /* Cycle 13 */                                          \
  _t##_i6##r = a2 - a3;                                   \
  _t##_i6##i = b2 - b3;                                   \
  /* Cycle 14 */                                          \
  a2 = _t##_i5##r + _t##_i7##i;/* t7 */                   \
  b2 = _t##_i5##i - _t##_i7##r;/* t7 */                   \
  /* Cycle 15 */                                          \
  _t##_i5##r -=  _t##_i7##i;                              \
  _t##_i5##i +=  _t##_i7##r;                              \
  /* Cycle 16 */                                          \
  _t##_i1##r -=  _t##_i3##i;                              \
  _t##_i1##i +=  _t##_i3##r;                              \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  a3 = _t##_i0##r + _t##_i4##r;/* t0 */                   \
  b3 = _t##_i0##i + _t##_i4##i;/* t0 */                   \
 /* Cycle  18 */                                          \
  _t##_i4##r = _t##_i0##r - _t##_i4##r;                   \
  _t##_i4##i = _t##_i0##i - _t##_i4##i;                   \
  /* Cycle 19 */                                          \
  _t##_i7##r = a2 + b2;                                   \
  _t##_i7##i = b2 - a2;                                   \
  /* Cycle 20 */                                          \
  a2 = _t##_i5##r - _t##_i5##i;                           \
  b2 = _t##_i5##i + _t##_i5##r;                           \
  /* Cycle 21 */                                          \
  a1 = _t##_i2##r + _t##_i6##i;/* t6 */                   \
  b1 = _t##_i2##i - _t##_i6##r;/* t6 */                   \
  /* Cycle 22 */                                          \
  _t##_i2##r -= _t##_i6##i;                               \
  _t##_i2##i += _t##_i6##r;                               \
  /* Cycle 23 */                                          \
  *(_pd0##r) = a3;                                        \
  *(_pd0##i) = b3;                                        \
  /* Cycle 24 */                                          \
  *(_pd4##r) = _t##_i4##r;                                \
  _t##_i3##r = _t##_i7##r * fi + a0;                      \
  *(_pd4##i) = _t##_i4##i;                                \
  _t##_i3##i = _t##_i7##i * fi + b0;                      \
  /* Cycle 25 */                                          \
  _t##_i7##r = _t##_i7##r * fr + a0;                      \
  _t##_i7##i = _t##_i7##i * fr + b0;                      \
  /* Cycle 26 */                                          \
  *(_pd2##r) = _t##_i2##r;                                \
  _t##_i5##r = a2 * fi + _t##_i1##r;                      \
  *(_pd2##i) = _t##_i2##i;                                \
  _t##_i5##i = b2 * fi + _t##_i1##i;                      \
  /* Cycle 27 */                                          \
  *( _pd6##r ) = a1;                                      \
  _t##_i1##r = a2 * fr + _t##_i1##r;                      \
  *( _pd6##i ) = b1;                                      \
  _t##_i1##i = b2 * fr + _t##_i1##i;                      \
  /* Cycle 28 */                                          \
  /* Cycle 29 */                                          \
  *(_pd3##r) = _t##_i3##r;                                \
  *(_pd3##i) = _t##_i3##i;                                \
  /* Cycle 30 */                                          \
  *(_pd7##r) = _t##_i7##r;                                \
  *(_pd7##i) = _t##_i7##i;                                \
  /* Cycle 31 */                                          \
  *(_pd5##r) = _t##_i5##r;                                \
  *(_pd5##i) = _t##_i5##i;                                \
  /* Cycle 32 */                                          \
  *(_pd1##r) = _t##_i1##r;                                \
  *(_pd1##i) = _t##_i1##i;                                \
}



                                                                                                                                                        /*
                                                                                                                                                        This macro performs a radix_8 DIF reduction of data and store to memory.
                                                                                                                                                        1) sets the preload pointers to actual pointers
                                                                                                                                                        2) computes the new preload pointers
                                                                                                                                                        3) Makes a fft-8
                                                                                                                                                        4) preload data for next loop
                                                                                                                                                        5) stores data to memory
                                                                                                                                                         
                                                                                                                                                        _pdi are the prefixes of pointers to store the data in memory
                                                                                                                                                        _t are the local actual data
                                                                                                                                                        _s are the preloaded data
                                                                                                                                                        _ij are the indexes of preloaded data
                                                                                                                                                        _ppd are the pointers of actual data
                                                                                                                                                        _ppdl are the pointers of preloaded data
                                                                                                                                                        _inc is the increment of preload pointers in this call
                                                                                                                                                         
                                                                                                                                                        The scheme:
                                                                                                                                                         
                                                                                                                                                          _t##j = _s##_i(j);
                                                                                                                                                          make_ a_fft8_with_t;
                                                                                                                                                          *(pd##j) = _t##j;
                                                                                                                                                         
                                                                                                                                                        */
#define cplx_fft8_preload_store_F(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_pd7,_t,_s,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7,_ppd,_ppdl,_inc) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   _ppd##_i0##r = _ppdl##_i0##r;                          \
   _ppd##_i0##i = _ppdl##_i0##i;                          \
   _t##0r = _s##_i0##r + _s##_i1##r;                      \
   _ppd##_i1##r = _ppdl##_i1##r;                          \
   _ppd##_i1##i = _ppdl##_i1##i;                          \
   _t##0i = _s##_i0##i + _s##_i1##i;                      \
  /* Cycle  2 */                                          \
   _ppd##_i2##r = _ppdl##_i2##r;                          \
   _ppd##_i2##i = _ppdl##_i2##i;                          \
   _t##1r = _s##_i0##r - _s##_i1##r;                      \
   _ppd##_i3##r = _ppdl##_i3##r;                          \
   _ppd##_i3##i = _ppdl##_i3##i;                          \
   _t##1i = _s##_i0##i - _s##_i1##i;                      \
  /* Cycle  3 */                                          \
   _ppd##_i4##r = _ppdl##_i4##r;                          \
   _ppd##_i4##i = _ppdl##_i4##i;                          \
   _t##2r = _s##_i2##r + _s##_i3##r;                      \
   _ppd##_i5##r = _ppdl##_i5##r;                          \
   _ppd##_i5##i = _ppdl##_i5##i;                          \
   _t##2i = _s##_i2##i + _s##_i3##i;                      \
  /* Cycle  4 */                                          \
   _ppd##_i6##r = _ppdl##_i6##r;                          \
   _ppd##_i6##i = _ppdl##_i6##i;                          \
   _t##3r = _s##_i2##r - _s##_i3##r;                      \
   _ppd##_i7##r = _ppdl##_i7##r;                          \
   _ppd##_i7##i = _ppdl##_i7##i;                          \
   _t##3i = _s##_i2##i - _s##_i3##i;                      \
  /* Cycle  5 */                                          \
   _ppdl##_i0##r += _inc;                                 \
   _ppdl##_i0##i += _inc;                                 \
   _t##4r = _s##_i4##r + _s##_i5##r;                      \
   _ppdl##_i1##r += _inc;                                 \
   _ppdl##_i1##i += _inc;                                 \
   _t##4i = _s##_i4##i + _s##_i5##i;                      \
  /* Cycle  6 */                                          \
   _ppdl##_i2##r += _inc;                                 \
   _ppdl##_i2##i += _inc;                                 \
   _t##5r = _s##_i4##r - _s##_i5##r;                      \
   _ppdl##_i3##r += _inc;                                 \
   _ppdl##_i3##i += _inc;                                 \
   _t##5i = _s##_i4##i - _s##_i5##i;                      \
  /* Cycle  7 */                                          \
   _t##6r = _s##_i6##r + _s##_i7##r;                      \
   _ppdl##_i4##r += _inc;                                 \
   _ppdl##_i4##i += _inc;                                 \
   _t##6i = _s##_i6##i + _s##_i7##i;                      \
   _ppdl##_i5##r += _inc;                                 \
   _ppdl##_i5##i += _inc;                                 \
  /* Cycle  8 */                                          \
   _t##7r = _s##_i6##r - _s##_i7##r;                      \
   _ppdl##_i6##r += _inc;                                 \
   _ppdl##_i6##i += _inc;                                 \
   _t##7i = _s##_i6##i - _s##_i7##i;                      \
   _ppdl##_i7##r += _inc;                                 \
   _ppdl##_i7##i += _inc;                                 \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _s##_i0##r = *( _ppdl##_i0##r );                        \
  a0 = _t##0r + _t##2r; /* t0 */                          \
  _s##_i0##i = *( _ppdl##_i0##i );                        \
  b0 = _t##0i + _t##2i; /* t0 */                          \
  /* Cycle 10 */                                          \
  _s##_i1##r = *( _ppdl##_i1##r );                        \
  _t##2r = _t##0r - _t##2r;                               \
  _s##_i1##i = *( _ppdl##_i1##i );                        \
  _t##2i = _t##0i - _t##2i;                               \
  /* Cycle 11 */                                          \
  _s##_i2##r = *( _ppdl##_i2##r );                        \
  a1 = _t##1r - _t##3i;        /* t3 */                   \
  _s##_i2##i = *( _ppdl##_i2##i );                        \
  b1 = _t##1i + _t##3r;        /* t3 */                   \
  /* Cycle 12 */                                          \
  _s##_i3##r = *( _ppdl##_i3##r );                        \
  a2 = _t##4r - _t##6r;/* t6 */                           \
  _s##_i3##i = *( _ppdl##_i3##i );                        \
  b2 = _t##4i - _t##6i;/* t6 */                           \
  /* Cycle 13 */                                          \
  _s##_i4##r = *( _ppdl##_i4##r );                        \
  _t##4r += _t##6r;                                       \
  _s##_i4##i = *( _ppdl##_i4##i );                        \
  _t##4i += _t##6i;                                       \
  /* Cycle 14 */                                          \
  _s##_i5##r = *( _ppdl##_i5##r );                        \
  a3 = _t##5r - _t##7i;        /* t7 */                   \
  _s##_i5##i = *( _ppdl##_i5##i );                        \
  b3 = _t##5i + _t##7r;        /* t7 */                   \
  /* Cycle 15 */                                          \
  _s##_i6##r = *( _ppdl##_i6##r );                        \
  _t##5r +=  _t##7i;                                      \
  _s##_i6##i = *( _ppdl##_i6##i );                        \
  _t##5i -=  _t##7r;                                      \
  /* Cycle 16 */                                          \
  _s##_i7##r = *( _ppdl##_i7##r );                        \
  _t##1r +=  _t##3i;                                      \
  _s##_i7##i = *( _ppdl##_i7##i );                        \
  _t##1i -=  _t##3r;                                      \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  _t##0r = a0 + _t##4r;                                   \
  _t##0i = b0 + _t##4i;                                   \
 /* Cycle  18 */                                          \
  _t##4r = a0 - _t##4r;                                   \
  _t##4i = b0 - _t##4i;                                   \
  /* Cycle 19 */                                          \
  _t##7r = a3 - b3;                                       \
  _t##7i = b3 + a3;                                       \
  /* Cycle 20 */                                          \
  a0 = _t##5r + _t##5i;                                   \
  b0 = _t##5i - _t##5r;                                   \
  /* Cycle 21 */                                          \
  _t##6r = _t##2r - b2;                                   \
  _t##6i = _t##2i + a2;                                   \
  /* Cycle 22 */                                          \
  _t##2r += b2;                                           \
  _t##2i -= a2;                                           \
  /* Cycle 23 */                                          \
  *(_pd0##r) = _t##0r;                                    \
  *(_pd0##i) = _t##0i;                                    \
  /* Cycle 24 */                                          \
  *(_pd4##r) = _t##4r;                                    \
  _t##3r = _t##7r * fi + a1;                              \
  *(_pd4##i) = _t##4i;                                    \
  _t##3i = _t##7i * fi + b1;                              \
  /* Cycle 25 */                                          \
  _t##7r = _t##7r * fr + a1;                              \
  _t##7i = _t##7i * fr + b1;                              \
  /* Cycle 26 */                                          \
  *( _pd2##r ) = _t##2r;                                  \
  _t##5r = a0 * fi + _t##1r;                              \
  *( _pd2##i ) = _t##2i;                                  \
  _t##5i = b0 * fi + _t##1i;                              \
  /* Cycle 27 */                                          \
  *( _pd6##r ) = _t##6r;                                  \
  _t##1r = a0 * fr + _t##1r;                              \
  *( _pd6##i ) = _t##6i;                                  \
  _t##1i = b0 * fr + _t##1i;                              \
  /* Cycle 28 */                                          \
  /* Cycle 29 */                                          \
  *(_pd3##r) = _t##3r;                                    \
  *(_pd3##i) = _t##3i;                                    \
  /* Cycle 30 */                                          \
  *(_pd7##r) = _t##7r;                                    \
  *(_pd7##i) = _t##7i;                                    \
  /* Cycle 31 */                                          \
  *(_pd5##r) = _t##5r;                                    \
  *(_pd5##i) = _t##5i;                                    \
  /* Cycle 32 */                                          \
  *(_pd1##r) = _t##1r;                                    \
  *(_pd1##i) = _t##1i;                                    \
}

                                                                                                                                                        /* This is for bakward DIT */
#define cplx_fft8_preload_store_B(_pd0,_pd1,_pd2,_pd3,_pd4,_pd5,_pd6,_pd7,_t,_s,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7,_ppd,_ppdl,_inc) \
{                                                         \
   y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                      \
  /******************** Pass1 ******************/         \
  /* Cycle  1 */                                          \
   _ppd##_i0##r = _ppdl##_i0##r;                          \
   _ppd##_i0##i = _ppdl##_i0##i;                          \
   _t##0r = _s##_i0##r + _s##_i1##r;                      \
   _ppd##_i1##r = _ppdl##_i1##r;                          \
   _ppd##_i1##i = _ppdl##_i1##i;                          \
   _t##0i = _s##_i0##i + _s##_i1##i;                      \
  /* Cycle  2 */                                          \
   _ppd##_i2##r = _ppdl##_i2##r;                          \
   _ppd##_i2##i = _ppdl##_i2##i;                          \
   _t##1r = _s##_i0##r - _s##_i1##r;                      \
   _ppd##_i3##r = _ppdl##_i3##r;                          \
   _ppd##_i3##i = _ppdl##_i3##i;                          \
   _t##1i = _s##_i0##i - _s##_i1##i;                      \
  /* Cycle  3 */                                          \
   _ppd##_i4##r = _ppdl##_i4##r;                          \
   _ppd##_i4##i = _ppdl##_i4##i;                          \
   _t##2r = _s##_i2##r + _s##_i3##r;                      \
   _ppd##_i5##r = _ppdl##_i5##r;                          \
   _ppd##_i5##i = _ppdl##_i5##i;                          \
   _t##2i = _s##_i2##i + _s##_i3##i;                      \
  /* Cycle  4 */                                          \
   _ppd##_i6##r = _ppdl##_i6##r;                          \
   _ppd##_i6##i = _ppdl##_i6##i;                          \
   _t##3r = _s##_i2##r - _s##_i3##r;                      \
   _ppd##_i7##r = _ppdl##_i7##r;                          \
   _ppd##_i7##i = _ppdl##_i7##i;                          \
   _t##3i = _s##_i2##i - _s##_i3##i;                      \
  /* Cycle  5 */                                          \
   _ppdl##_i0##r += _inc;                                 \
   _ppdl##_i0##i += _inc;                                 \
   _t##4r = _s##_i4##r + _s##_i5##r;                      \
   _ppdl##_i1##r += _inc;                                 \
   _ppdl##_i1##i += _inc;                                 \
   _t##4i = _s##_i4##i + _s##_i5##i;                      \
  /* Cycle  6 */                                          \
   _ppdl##_i2##r += _inc;                                 \
   _ppdl##_i2##i += _inc;                                 \
   _t##5r = _s##_i4##r - _s##_i5##r;                      \
   _ppdl##_i3##r += _inc;                                 \
   _ppdl##_i3##i += _inc;                                 \
   _t##5i = _s##_i4##i - _s##_i5##i;                      \
  /* Cycle  7 */                                          \
   _t##6r = _s##_i6##r + _s##_i7##r;                      \
   _ppdl##_i4##r += _inc;                                 \
   _ppdl##_i4##i += _inc;                                 \
   _t##6i = _s##_i6##i + _s##_i7##i;                      \
   _ppdl##_i5##r += _inc;                                 \
   _ppdl##_i5##i += _inc;                                 \
  /* Cycle  8 */                                          \
   _t##7r = _s##_i6##r - _s##_i7##r;                      \
   _ppdl##_i6##r += _inc;                                 \
   _ppdl##_i6##i += _inc;                                 \
   _t##7i = _s##_i6##i - _s##_i7##i;                      \
   _ppdl##_i7##r += _inc;                                 \
   _ppdl##_i7##i += _inc;                                 \
  /**************************** Pass2    **********/      \
  /* Cycle  9 */                                          \
  _s##_i0##r = *( _ppdl##_i0##r );                        \
  a0 = _t##0r + _t##2r; /* t0 */                          \
  _s##_i0##i = *( _ppdl##_i0##i );                        \
  b0 = _t##0i + _t##2i; /* t0 */                          \
  /* Cycle 10 */                                          \
  _s##_i1##r = *( _ppdl##_i1##r );                        \
  _t##2r = _t##0r - _t##2r;                               \
  _s##_i1##i = *( _ppdl##_i1##i );                        \
  _t##2i = _t##0i - _t##2i;                               \
  /* Cycle 11 */                                          \
  _s##_i2##r = *( _ppdl##_i2##r );                        \
  a1 = _t##1r + _t##3i;        /* t3 */                   \
  _s##_i2##i = *( _ppdl##_i2##i );                        \
  b1 = _t##1i - _t##3r;        /* t3 */                   \
  /* Cycle 12 */                                          \
  _s##_i3##r = *( _ppdl##_i3##r );                        \
  a2 = _t##4r - _t##6r;/* t6 */                           \
  _s##_i3##i = *( _ppdl##_i3##i );                        \
  b2 = _t##4i - _t##6i;/* t6 */                           \
  /* Cycle 13 */                                          \
  _s##_i4##r = *( _ppdl##_i4##r );                        \
  _t##4r += _t##6r;                                       \
  _s##_i4##i = *( _ppdl##_i4##i );                        \
  _t##4i += _t##6i;                                       \
  /* Cycle 14 */                                          \
  _s##_i5##r = *( _ppdl##_i5##r );                        \
  a3 = _t##5r + _t##7i;        /* t7 */                   \
  _s##_i5##i = *( _ppdl##_i5##i );                        \
  b3 = _t##5i - _t##7r;        /* t7 */                   \
  /* Cycle 15 */                                          \
  _s##_i6##r = *( _ppdl##_i6##r );                        \
  _t##5r -=  _t##7i;                                      \
  _s##_i6##i = *( _ppdl##_i6##i );                        \
  _t##5i +=  _t##7r;                                      \
  /* Cycle 16 */                                          \
  _s##_i7##r = *( _ppdl##_i7##r );                        \
  _t##1r -=  _t##3i;                                      \
  _s##_i7##i = *( _ppdl##_i7##i );                        \
  _t##1i +=  _t##3r;                                      \
  /**************************** Pass3    **********/      \
  /* Cycle 17 */                                          \
  _t##0r = a0 + _t##4r;                                   \
  _t##0i = b0 + _t##4i;                                   \
 /* Cycle  18 */                                          \
  _t##4r = a0 - _t##4r;                                   \
  _t##4i = b0 - _t##4i;                                   \
  /* Cycle 19 */                                          \
  _t##7r = a3 + b3;                                       \
  _t##7i = b3 - a3;                                       \
  /* Cycle 20 */                                          \
  a0 = _t##5r - _t##5i;                                   \
  b0 = _t##5i + _t##5r;                                   \
  /* Cycle 21 */                                          \
  _t##6r = _t##2r + b2;                                   \
  _t##6i = _t##2i - a2;                                   \
  /* Cycle 22 */                                          \
  _t##2r -= b2;                                           \
  _t##2i += a2;                                           \
  /* Cycle 23 */                                          \
  *(_pd0##r) = _t##0r;                                    \
  *(_pd0##i) = _t##0i;                                    \
  /* Cycle 24 */                                          \
  *(_pd4##r) = _t##4r;                                    \
  _t##3r = _t##7r * fi + a1;                              \
  *(_pd4##i) = _t##4i;                                    \
  _t##3i = _t##7i * fi + b1;                              \
  /* Cycle 25 */                                          \
  _t##7r = _t##7r * fr + a1;                              \
  _t##7i = _t##7i * fr + b1;                              \
  /* Cycle 26 */                                          \
  *( _pd2##r ) = _t##2r;                                  \
  _t##5r = a0 * fi + _t##1r;                              \
  *( _pd2##i ) = _t##2i;                                  \
  _t##5i = b0 * fi + _t##1i;                              \
  /* Cycle 27 */                                          \
  *( _pd6##r ) = _t##6r;                                  \
  _t##1r = a0 * fr + _t##1r;                              \
  *( _pd6##i ) = _t##6i;                                  \
  _t##1i = b0 * fr + _t##1i;                              \
  /* Cycle 28 */                                          \
  /* Cycle 29 */                                          \
  *(_pd3##r) = _t##3r;                                    \
  *(_pd3##i) = _t##3i;                                    \
  /* Cycle 30 */                                          \
  *(_pd7##r) = _t##7r;                                    \
  *(_pd7##i) = _t##7i;                                    \
  /* Cycle 31 */                                          \
  *(_pd5##r) = _t##5r;                                    \
  *(_pd5##i) = _t##5i;                                    \
  /* Cycle 32 */                                          \
  *(_pd1##r) = _t##1r;                                    \
  *(_pd1##i) = _t##1i;                                    \
}


                                                                                                                                                        /*
                                                                                                                                                        This macro performs two radix_4 DIF reduction of data and store to memory.
                                                                                                                                                        1) Makes two fft-4
                                                                                                                                                        2) stores data to memory
                                                                                                                                                         
                                                                                                                                                        _pj,_qj are the prefixes of pointers to store the data in memory
                                                                                                                                                        _t,_u are the local actual data
                                                                                                                                                        _ij are the indexes of preloaded data
                                                                                                                                                         
                                                                                                                                                        The scheme:
                                                                                                                                                         
                                                                                                                                                          _t##j = _s##_i(j);
                                                                                                                                                          make_ a_fft8_with_t;
                                                                                                                                                          *(pd##j) = _t##j;
                                                                                                                                                         
                                                                                                                                                        */
#define cplx_two_fft4_preload_store_F(_p0,_p1,_p2,_p3,_q0,_q1,_q2,_q3,_t,_u,_i0,_i1,_i2,_i3) \
{                                                         \
  y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                       \
  /* Pass 1 *******************/                          \
  /* Cycle  1 */                                          \
  a0 = _t##_i0##r - _t##_i1##r;/*t1*/                     \
  b0 = _t##_i0##i - _t##_i1##i;                           \
  /* Cycle  2 */                                          \
  a1 = _u##_i0##r - _u##_i1##r;/*u1*/                     \
  b1 = _u##_i0##i - _u##_i1##i;                           \
  /* Cycle  3 */                                          \
  a2 = _t##_i2##r - _t##_i3##r;/*t3*/                     \
  b2 = _t##_i2##i - _t##_i3##i;/*t3*/                     \
  /* Cycle  4 */                                          \
  a3 = _u##_i2##r - _u##_i3##r;/*u3*/                     \
  b3 = _u##_i2##i - _u##_i3##i;/*u3*/                     \
  /* Cycle  5 */                                          \
  _t##_i0##r += _t##_i1##r;                               \
  _t##_i0##i += _t##_i1##i;                               \
  /* Cycle  6 */                                          \
  _u##_i0##r += _u##_i1##r;                               \
  _u##_i0##i += _u##_i1##i;                               \
  /* Cycle  7 */                                          \
  _t##_i2##r += _t##_i3##r;                               \
  _t##_i2##i += _t##_i3##i;                               \
  /* Cycle  8 */                                          \
  _u##_i2##r += _u##_i3##r;                               \
  _u##_i2##i += _u##_i3##i;                               \
  /* Pass 2 *******************/                          \
  /* Cycle  9 */                                          \
  _t##_i3##r = a0 - b2;                                   \
  _t##_i3##i = b0 + a2;                                   \
  /* Cycle 10 */                                          \
  _t##_i1##r = a0 + b2;                                   \
  _t##_i1##i = b0 - a2;                                   \
  /* Cycle 11 */                                          \
  _u##_i3##r = a1 - b3;                                   \
  _u##_i3##i = b1 + a3;                                   \
  /* Cycle 12 */                                          \
  _u##_i1##r = a1 + b3;                                   \
  _u##_i1##i = b1 - a3;                                   \
  /* Cycle 13 */                                          \
  a0 = _t##_i0##r - _t##_i2##r;/* t2 */                   \
  b0 = _t##_i0##i - _t##_i2##i;/* t2 */                   \
  /* Cycle 14 */                                          \
  *(_p3##r) = _t##_i3##r;                                 \
  _t##_i0##r += _t##_i2##r;                               \
  *(_p3##i) = _t##_i3##i;                                 \
  _t##_i0##i += _t##_i2##i;                               \
  /* Cycle 15 */                                          \
  *(_p1##r) = _t##_i1##r;                                 \
  a1 = _u##_i0##r - _u##_i2##r;/* u2 */                   \
  *(_p1##i) = _t##_i1##i;                                 \
  b1 = _u##_i0##i - _u##_i2##i;/* u2 */                   \
  /* Cycle 16 */                                          \
  *(_q3##r) = _u##_i3##r;                                 \
  _u##_i0##r += _u##_i2##r;                               \
  *(_q3##i) = _u##_i3##i;                                 \
  _u##_i0##i += _u##_i2##i;                               \
  /* Cycle 17 */                                          \
  *(_q1##r) = _u##_i1##r;                                 \
  *(_q1##i) = _u##_i1##i;                                 \
  /* Cycle 18 */                                          \
  *(_p2##r) = a0;                                         \
  *(_p2##i) = b0;                                         \
  /* Cycle 19 */                                          \
  *(_p0##r) = _t##_i0##r;                                 \
  *(_p0##i) = _t##_i0##i;                                 \
  /* Cycle 20 */                                          \
  *(_q2##r) = a1;                                         \
  *(_q2##i) = b1;                                         \
  /* Cycle 21 */                                          \
  *(_q0##r) = _u##_i0##r;                                 \
  *(_q0##i) = _u##_i0##i;                                 \
}

                                                                                                                                                        /* Idem for backward DIT */
#define cplx_two_fft4_preload_store_B(_p0,_p1,_p2,_p3,_q0,_q1,_q2,_q3,_t,_u,_i0,_i1,_i2,_i3) \
{                                                         \
  y_limb_t a0,a1,a2,a3,b0,b1,b2,b3;                       \
  /* Pass 1 *******************/                          \
  /* Cycle  1 */                                          \
  a0 = _t##_i0##r - _t##_i1##r;/*t1*/                     \
  b0 = _t##_i0##i - _t##_i1##i;                           \
  /* Cycle  2 */                                          \
  a1 = _u##_i0##r - _u##_i1##r;/*u1*/                     \
  b1 = _u##_i0##i - _u##_i1##i;                           \
  /* Cycle  3 */                                          \
  a2 = _t##_i2##r - _t##_i3##r;/*t3*/                     \
  b2 = _t##_i2##i - _t##_i3##i;/*t3*/                     \
  /* Cycle  4 */                                          \
  a3 = _u##_i2##r - _u##_i3##r;/*u3*/                     \
  b3 = _u##_i2##i - _u##_i3##i;/*u3*/                     \
  /* Cycle  5 */                                          \
  _t##_i0##r += _t##_i1##r;                               \
  _t##_i0##i += _t##_i1##i;                               \
  /* Cycle  6 */                                          \
  _u##_i0##r += _u##_i1##r;                               \
  _u##_i0##i += _u##_i1##i;                               \
  /* Cycle  7 */                                          \
  _t##_i2##r += _t##_i3##r;                               \
  _t##_i2##i += _t##_i3##i;                               \
  /* Cycle  8 */                                          \
  _u##_i2##r += _u##_i3##r;                               \
  _u##_i2##i += _u##_i3##i;                               \
  /* Pass 2 *******************/                          \
  /* Cycle  9 */                                          \
  _t##_i3##r = a0 + b2;                                   \
  _t##_i3##i = b0 - a2;                                   \
  /* Cycle 10 */                                          \
  _t##_i1##r = a0 - b2;                                   \
  _t##_i1##i = b0 + a2;                                   \
  /* Cycle 11 */                                          \
  _u##_i3##r = a1 + b3;                                   \
  _u##_i3##i = b1 - a3;                                   \
  /* Cycle 12 */                                          \
  _u##_i1##r = a1 - b3;                                   \
  _u##_i1##i = b1 + a3;                                   \
  /* Cycle 13 */                                          \
  a0 = _t##_i0##r - _t##_i2##r;/* t2 */                   \
  b0 = _t##_i0##i - _t##_i2##i;/* t2 */                   \
  /* Cycle 14 */                                          \
  *(_p3##r) = _t##_i3##r;                                 \
  _t##_i0##r += _t##_i2##r;                               \
  *(_p3##i) = _t##_i3##i;                                 \
  _t##_i0##i += _t##_i2##i;                               \
  /* Cycle 15 */                                          \
  *(_p1##r) = _t##_i1##r;                                 \
  a1 = _u##_i0##r - _u##_i2##r;/* u2 */                   \
  *(_p1##i) = _t##_i1##i;                                 \
  b1 = _u##_i0##i - _u##_i2##i;/* u2 */                   \
  /* Cycle 16 */                                          \
  *(_q3##r) = _u##_i3##r;                                 \
  _u##_i0##r += _u##_i2##r;                               \
  *(_q3##i) = _u##_i3##i;                                 \
  _u##_i0##i += _u##_i2##i;                               \
  /* Cycle 17 */                                          \
  *(_q1##r) = _u##_i1##r;                                 \
  *(_q1##i) = _u##_i1##i;                                 \
  /* Cycle 18 */                                          \
  *(_p2##r) = a0;                                         \
  *(_p2##i) = b0;                                         \
  /* Cycle 19 */                                          \
  *(_p0##r) = _t##_i0##r;                                 \
  *(_p0##i) = _t##_i0##i;                                 \
  /* Cycle 20 */                                          \
  *(_q2##r) = a1;                                         \
  *(_q2##i) = b1;                                         \
  /* Cycle 21 */                                          \
  *(_q0##r) = _u##_i0##r;                                 \
  *(_q0##i) = _u##_i0##i;                                 \
}

                                                                                                                                                        /*
                                                                                                                                                        This macro performs two radix_4 DIF reduction of data and store to memory.
                                                                                                                                                        1) sets the preload pointers to actual pointers
                                                                                                                                                        2) computes the new preload pointers
                                                                                                                                                        3) Makes two fft-4
                                                                                                                                                        4) preload data for next loop
                                                                                                                                                        5) stores data to memory
                                                                                                                                                         
                                                                                                                                                        _pi,_qi are the prefixes of pointers to store the data in memory
                                                                                                                                                        _t,_u are the local actual data
                                                                                                                                                        _s,_v are the preloaded data
                                                                                                                                                        _ij are the indexes of preloaded data
                                                                                                                                                        _pd,_qd are the pointers of actual data
                                                                                                                                                        _pdl,_qdl are the pointers of preloaded data
                                                                                                                                                        _inc is the increment of preload pointers in this call
                                                                                                                                                         
                                                                                                                                                        The scheme:
                                                                                                                                                         
                                                                                                                                                          _t##j = _s##_i(j);
                                                                                                                                                          _u##j = _v##_i(j);
                                                                                                                                                          make_two_fft4_with_t_s;
                                                                                                                                                          *(_pd##j) = _t##j;
                                                                                                                                                          *(_qd##j) = _u##j;
                                                                                                                                                         
                                                                                                                                                        */

#define cplx_two_fft4_notwd_preload_store_F(_p0,_p1,_p2,_p3,_q0,_q1,_q2,_q3,_t,_s,_u,_v,_i0,_i1,_i2,_i3,_pd,_pdl,_qd,_qdl,_inc) \
{                                                         \
  y_limb_t a0,a1,b0,b1;                                   \
  /* Pass 1 *******************/                          \
  /* Cycle  1 */                                          \
  _pd##_i0##r = _pdl##_i0##r;                             \
  _pd##_i0##i = _pdl##_i0##i;                             \
  _t##1r = _s##_i0##r - _s##_i1##r;                       \
  _pd##_i1##r = _pdl##_i1##r;                             \
  _pd##_i1##i = _pdl##_i1##i;                             \
  _t##1i = _s##_i0##i - _s##_i1##i;                       \
  /* Cycle  2 */                                          \
  _qd##_i0##r = _qdl##_i0##r;                             \
  _qd##_i0##i = _qdl##_i0##i;                             \
  _u##1r = _v##_i0##r - _v##_i1##r;                       \
  _qd##_i1##r = _qdl##_i1##r;                             \
  _qd##_i1##i = _qdl##_i1##i;                             \
  _u##1i = _v##_i0##i - _v##_i1##i;                       \
  /* Cycle  3 */                                          \
  _pd##_i2##r = _pdl##_i2##r;                             \
  _pd##_i2##i = _pdl##_i2##i;                             \
  a0 = _s##_i2##r - _s##_i3##r;/*t3*/                     \
  _pd##_i3##r = _pdl##_i3##r;                             \
  _pd##_i3##i = _pdl##_i3##i;                             \
  b0 = _s##_i2##i - _s##_i3##i;/*t3*/                     \
  /* Cycle  4 */                                          \
  _qd##_i2##r = _qdl##_i2##r;                             \
  _qd##_i2##i = _qdl##_i2##i;                             \
  a1 = _v##_i2##r - _v##_i3##r;/*u3*/                     \
  _qd##_i3##r = _qdl##_i3##r;                             \
  _qd##_i3##i = _qdl##_i3##i;                             \
  b1 = _v##_i2##i - _v##_i3##i;/*u3*/                     \
  /* Cycle  5 */                                          \
  _pdl##_i0##r += _inc;                                   \
  _pdl##_i0##i += _inc;                                   \
  _t##0r = _s##_i0##r + _s##_i1##r;                       \
  _qdl##_i0##r += _inc;                                   \
  _qdl##_i0##i += _inc;                                   \
  _t##0i = _s##_i0##i + _s##_i1##i;                       \
  /* Cycle  6 */                                          \
  _pdl##_i1##r += _inc;                                   \
  _pdl##_i1##i += _inc;                                   \
  _u##0r = _v##_i0##r + _v##_i1##r;                       \
  _qdl##_i1##r += _inc;                                   \
  _qdl##_i1##i += _inc;                                   \
  _u##0i = _v##_i0##i + _v##_i1##i;                       \
  /* Cycle  7 */                                          \
  _pdl##_i2##r += _inc;                                   \
  _pdl##_i2##i += _inc;                                   \
  _t##2r = _s##_i2##r + _s##_i3##r;                       \
  _qdl##_i2##r += _inc;                                   \
  _qdl##_i2##i += _inc;                                   \
  _t##2i = _s##_i2##i + _s##_i3##i;                       \
  /* Cycle  8 */                                          \
  _pdl##_i3##r += _inc;                                   \
  _pdl##_i3##i += _inc;                                   \
  _u##2r = _v##_i2##r + _v##_i3##r;                       \
  _qdl##_i3##r += _inc;                                   \
  _qdl##_i3##i += _inc;                                   \
  _u##2i = _v##_i2##i + _v##_i3##i;                       \
  /* Pass 2 *******************/                          \
  /* Cycle  9 */                                          \
  _s##_i0##r = *(_pdl##_i0##r);                           \
  _t##3r = _t##1r - b0;                                   \
  _s##_i0##i = *(_pdl##_i0##i);                           \
  _t##3i = _t##1i + a0;                                   \
  /* Cycle 10 */                                          \
  _v##_i0##r = *(_qdl##_i0##r);                           \
  _t##1r += b0;                                           \
  _v##_i0##i = *(_qdl##_i0##i);                           \
  _t##1i -= a0;                                           \
  /* Cycle 11 */                                          \
  _s##_i1##r = *(_pdl##_i1##r);                           \
  _u##3r = _u##1r - b1;                                   \
  _s##_i1##i = *(_pdl##_i1##i);                           \
  _u##3i = _u##1i + a1;                                   \
  /* Cycle 12 */                                          \
  _v##_i1##r = *(_qdl##_i1##r);                           \
  _u##1r += b1;                                           \
  _v##_i1##i = *(_qdl##_i1##i);                           \
  _u##1i -= a1;                                           \
  /* Cycle 13 */                                          \
  _s##_i2##r = *(_pdl##_i2##r);                           \
  _s##_i2##i = *(_pdl##_i2##i);                           \
  a0 = _t##0r - _t##2r;/* t2 */                           \
  /* Cycle 14 */                                          \
  _v##_i2##r = *(_qdl##_i2##r);                           \
  _v##_i2##i = *(_qdl##_i2##i);                           \
  b0 = _t##0i - _t##2i;/* t2 */                           \
  /* Cycle 15 */                                          \
  _s##_i3##r = *(_pdl##_i3##r);                           \
  _s##_i3##i = *(_pdl##_i3##i);                           \
  _t##0r += _t##2r;                                       \
  /* Cycle 16 */                                          \
  _v##_i3##r = *(_qdl##_i3##r);                           \
  _v##_i3##i = *(_qdl##_i3##i);                           \
  _t##0i += _t##2i;                                       \
  /* Cycle 17 */                                          \
  *(_p3##r) = _t##3r;                                     \
  *(_p3##i) = _t##3i;                                     \
  /* Cycle 18 */                                          \
  *(_p1##r) = _t##1r;                                     \
  a1 = _u##0r - _u##2r;/* u2 */                           \
  *(_p1##i) = _t##1i;                                     \
  b1 = _u##0i - _u##2i;/* u2 */                           \
  /* Cycle 19 */                                          \
  *(_q3##r) = _u##3r;                                     \
  _u##0r += _u##2r;                                       \
  *(_q3##i) = _u##3i;                                     \
  _u##0i += _u##2i;                                       \
  /* Cycle 20 */                                          \
  *(_q1##r) = _u##1r;                                     \
  *(_q1##i) = _u##1i;                                     \
  /* Cycle 21 */                                          \
  *(_p2##r) = a0;                                         \
  *(_p2##i) = b0;                                         \
  /* Cycle 22 */                                          \
  *(_p0##r) = _t##0r;                                     \
  *(_p0##i) = _t##0i;                                     \
  /* Cycle 23 */                                          \
  *(_q2##r) = a1;                                         \
  *(_q2##i) = b1;                                         \
  /* Cycle 24 */                                          \
  *(_q0##r) = _u##0r;                                     \
  *(_q0##i) = _u##0i;                                     \
}








