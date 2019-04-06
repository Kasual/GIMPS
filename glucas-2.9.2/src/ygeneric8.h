/* $Id$ */
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
   Macro to make a first radix-8 pass with twiddle factors and prefetch data 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_twd_pp_pass1(_t,_f,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;\
  a0 = *(_pd##4 + _i);\
  a1 = *(_pd##4 + _i + 1);\
  a2 = a0 * (_f##4##r);\
  prefetch_p(_pd##4 + _i);\
  a3 = a1 * (_f##4##i);\
  a4 = *(_pd##0 + _i);\
  a5 = *(_pd##0 + _i + 1);\
  a6 = a0 * (_f##4##i);\
  prefetch_p(_pd##0 + _i);\
  a7 = a1 * (_f##4##r);\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 -= a3;\
  prefetch_p(_pd##2 + _i);\
  a3 = a0 * (_f##2##r);\
  a6 += a7;\
  _t##0##r = a4 + a2;\
  a7 = a1 * (_f##2##i);\
  _t##1##r = a4 - a2;\
  a4 = a0 * (_f##2##i);\
  _t##0##i = a5 + a6;\
  a2 = a1 * (_f##2##r);\
  _t##1##i = a5 - a6;\
  a0 = *(_pd##6 + _i);\
  a1 = *(_pd##6 + _i + 1);\
  a5 = a0 * (_f##6##r);\
  a3 -= a7;\
  a6 = a1 * (_f##6##i);\
  prefetch_p(_pd##6 + _i);\
  a2 += a4;\
  a4 = a0 * (_f##6##i);\
  a5 -= a6;\
  a6 = a1 * (_f##6##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  prefetch_p(_pd##1 + _i);\
  a4 += a6;\
  a6 = a0 * (_f##1##r);\
  _t##2##r = a3 + a5;\
  a7 = a1 * (_f##1##i);\
  _t##3##r = a3 - a5;\
  a3 = a0 * (_f##1##i);\
  _t##2##i = a2 + a4;\
  a5 = a1 * (_f##1##r);\
  _t##3##i = a2 - a4;\
  a0 = *(_pd##5 + _i);\
  a1 = *(_pd##5 + _i + 1);\
  a2 = a0 * (_f##5##r);\
  a6 -= a7;\
  a4 = a1 * (_f##5##i);\
  a3 += a5;\
  prefetch_p(_pd##5 + _i);\
  a7 = a0 * (_f##5##i);\
  a2 -= a4;\
  a5 = a1 * (_f##5##r);\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  prefetch_p(_pd##3 + _i);\
  a5 += a7;\
  a4 = a0 * (_f##3##r);\
  _t##4##r = a6 + a2;\
  a7 = a1 * (_f##3##i);\
  _t##5##r = a6 - a2;\
  a6 = a0 * (_f##3##i);\
  _t##4##i = a3 + a5;\
  a2 = a1 * (_f##3##r);\
  _t##5##i = a3 - a5;\
  a0 = *(_pd##7 + _i);\
  a1 = *(_pd##7 + _i + 1);\
  a3 = a0 * (_f##7##r);\
  a4 -= a7;\
  a5 = a1 * (_f##7##i);\
  a2 += a6;\
  a6 = a0 * (_f##7##i);\
  a3 -= a5;\
  a5 = a1 * (_f##7##r);\
  prefetch_p(_pd##7 + _i);\
  _t##6##r = a4 + a3;\
  a6 += a5;\
  _t##7##r = a4 - a3;\
  _t##6##i = a2 + a6;\
  _t##7##i = a2 - a6;\
}

/*
   Macro to make a first radix-8 pass with twiddle factors 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_twd_pp_pass1_no_fetch(_t,_f,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;\
  a0 = *(_pd##4 + _i);\
  a1 = *(_pd##4 + _i + 1);\
  a2 = a0 * (_f##4##r);\
  a3 = a1 * (_f##4##i);\
  a4 = *(_pd##0 + _i);\
  a5 = *(_pd##0 + _i + 1);\
  a6 = a0 * (_f##4##i);\
  a7 = a1 * (_f##4##r);\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 -= a3;\
  a3 = a0 * (_f##2##r);\
  a6 += a7;\
  _t##0##r = a4 + a2;\
  a7 = a1 * (_f##2##i);\
  _t##1##r = a4 - a2;\
  a4 = a0 * (_f##2##i);\
  _t##0##i = a5 + a6;\
  a2 = a1 * (_f##2##r);\
  _t##1##i = a5 - a6;\
  a0 = *(_pd##6 + _i);\
  a1 = *(_pd##6 + _i + 1);\
  a5 = a0 * (_f##6##r);\
  a3 -= a7;\
  a6 = a1 * (_f##6##i);\
  a2 += a4;\
  a4 = a0 * (_f##6##i);\
  a5 -= a6;\
  a6 = a1 * (_f##6##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  a4 += a6;\
  a6 = a0 * (_f##1##r);\
  _t##2##r = a3 + a5;\
  a7 = a1 * (_f##1##i);\
  _t##3##r = a3 - a5;\
  a3 = a0 * (_f##1##i);\
  _t##2##i = a2 + a4;\
  a5 = a1 * (_f##1##r);\
  _t##3##i = a2 - a4;\
  a0 = *(_pd##5 + _i);\
  a1 = *(_pd##5 + _i + 1);\
  a2 = a0 * (_f##5##r);\
  a6 -= a7;\
  a4 = a1 * (_f##5##i);\
  a3 += a5;\
  a7 = a0 * (_f##5##i);\
  a2 -= a4;\
  a5 = a1 * (_f##5##r);\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  a5 += a7;\
  a4 = a0 * (_f##3##r);\
  _t##4##r = a6 + a2;\
  a7 = a1 * (_f##3##i);\
  _t##5##r = a6 - a2;\
  a6 = a0 * (_f##3##i);\
  _t##4##i = a3 + a5;\
  a2 = a1 * (_f##3##r);\
  _t##5##i = a3 - a5;\
  a0 = *(_pd##7 + _i);\
  a1 = *(_pd##7 + _i + 1);\
  a3 = a0 * (_f##7##r);\
  a4 -= a7;\
  a5 = a1 * (_f##7##i);\
  a2 += a6;\
  a6 = a0 * (_f##7##i);\
  a3 -= a5;\
  a5 = a1 * (_f##7##r);\
  _t##6##r = a4 + a3;\
  a6 += a5;\
  _t##7##r = a4 - a3;\
  _t##6##i = a2 + a6;\
  _t##7##i = a2 - a6;\
}


/*
   Macro to make a first radix-8 pass without twiddle factors and 
   with prefetch data 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_notwd_pp_pass1(_t,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##0 + _i);\
  a1 = *(_pd##0 + _i + 1);\
  prefetch_p(_pd##0 + _i);\
  a2 = *(_pd##4 + _i);\
  a3 = *(_pd##4 + _i + 1);\
  prefetch_p(_pd##4 + _i);\
  _t##0##r = a0 + a2;\
  _t##1##r = a0 - a2;\
  _t##0##i = a1 + a3;\
  _t##1##i = a1 - a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  prefetch_p(_pd##2 + _i);\
  a2 = *(_pd##6 + _i);\
  a3 = *(_pd##6 + _i + 1);\
  prefetch_p(_pd##6 + _i);\
  _t##2##r = a0 + a2;\
  _t##3##r = a0 - a2;\
  _t##2##i = a1 + a3;\
  _t##3##i = a1 - a3;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  prefetch_p(_pd##1 + _i);\
  a2 = *(_pd##5 + _i);\
  a3 = *(_pd##5 + _i + 1);\
  prefetch_p(_pd##5 + _i);\
  _t##4##r = a0 + a2;\
  _t##5##r = a0 - a2;\
  _t##4##i = a1 + a3;\
  _t##5##i = a1 - a3;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  prefetch_p(_pd##3 + _i);\
  a2 = *(_pd##7 + _i);\
  a3 = *(_pd##7 + _i + 1);\
  prefetch_p(_pd##7 + _i);\
  _t##6##r = a0 + a2;\
  _t##7##r = a0 - a2;\
  _t##6##i = a1 + a3;\
  _t##7##i = a1 - a3;\
}

/*
   Macro to make a first radix-8 pass without twiddle factors and 
   without prefetch data 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_notwd_pp_pass1_no_fetch(_t,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##0 + _i);\
  a1 = *(_pd##0 + _i + 1);\
  a2 = *(_pd##4 + _i);\
  a3 = *(_pd##4 + _i + 1);\
  _t##0##r = a0 + a2;\
  _t##1##r = a0 - a2;\
  _t##0##i = a1 + a3;\
  _t##1##i = a1 - a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = *(_pd##6 + _i);\
  a3 = *(_pd##6 + _i + 1);\
  _t##2##r = a0 + a2;\
  _t##3##r = a0 - a2;\
  _t##2##i = a1 + a3;\
  _t##3##i = a1 - a3;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  a2 = *(_pd##5 + _i);\
  a3 = *(_pd##5 + _i + 1);\
  _t##4##r = a0 + a2;\
  _t##5##r = a0 - a2;\
  _t##4##i = a1 + a3;\
  _t##5##i = a1 - a3;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  a2 = *(_pd##7 + _i);\
  a3 = *(_pd##7 + _i + 1);\
  _t##6##r = a0 + a2;\
  _t##7##r = a0 - a2;\
  _t##6##i = a1 + a3;\
  _t##7##i = a1 - a3;\
}

/*
   Macro to make a second radix-8 pass in a DIF  
   _t is the locale prefix variable name
*/
#define cplx_radix_8_pp_pass2_F(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##2##r;\
  a1 = _t##2##i;\
  _t##2##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##2##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##3##r;\
  a1 = _t##3##i;\
  _t##3##r = _t##1##r - a1;\
  _t##1##r += a1;\
  _t##3##i = _t##1##i + a0;\
  _t##1##i -= a0;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##4##r - a0;\
  _t##4##r += a0;\
  _t##6##i = _t##4##i - a1;\
  _t##4##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##5##r - a1;\
  _t##5##r += a1;\
  _t##7##i = _t##5##i + a0;\
  _t##5##i -= a0;\
}

/*
   Macro to make a second radix-8 pass in a DIT  
   _t is the locale prefix variable name
*/
#define cplx_radix_8_pp_pass2_B(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##2##r;\
  a1 = _t##2##i;\
  _t##2##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##2##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##3##r;\
  a1 = _t##3##i;\
  _t##3##r = _t##1##r + a1;\
  _t##1##r -= a1;\
  _t##3##i = _t##1##i - a0;\
  _t##1##i += a0;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##4##r - a0;\
  _t##4##r += a0;\
  _t##6##i = _t##4##i - a1;\
  _t##4##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##5##r +a1;\
  _t##5##r -= a1;\
  _t##7##i = _t##5##i - a0;\
  _t##5##i += a0;\
}


/*
   Macro to make the third radix-8 pass and store in a DIF 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_pp_pass3_store_F(_t,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r + _t##5##i;\
  a1 = _t##5##i - _t##5##r;\
  a2 = _t##0##r + _t##4##r;\
  a3 = _t##0##i + _t##4##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##4##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##4##i;\
  *( _pd##0 + _i) = a2;\
  *( _pd##0 + _i + 1) =a3;\
  a2 = _t##1##r + a0;\
  a3 = _t##1##i + a1;\
  *( _pd##4 + _i) = a4;\
  *( _pd##4 + _i + 1) =a5;\
  a4 = _t##1##r - a0;\
  a5 = _t##1##i - a1;\
  *( _pd##1 + _i) = a2;\
  *( _pd##1 + _i + 1) =a3;\
  a2 = _t##2##r + _t##6##i;\
  a3 = _t##2##i - _t##6##r;\
  *( _pd##5 + _i) = a4;\
  *( _pd##5 + _i + 1) =a5;\
  a0 = _t##7##r - _t##7##i;\
  a1 = _t##7##i + _t##7##r;\
  *( _pd##2+ _i) = a2;\
  *( _pd##2+ _i + 1) =a3;\
  a4= _t##2##r - _t##6##i;\
  a0 *= F_1_8i;\
  a5= _t##2##i + _t##6##r;\
  a1 *= F_1_8i;\
  *( _pd##6+ _i) = a4;\
  *( _pd##6+ _i + 1) =a5;\
  a4 = _t##3##r + a0;\
  a5 = _t##3##i + a1;\
  a2 = _t##3##r - a0;\
  a3 = _t##3##i - a1;\
  *( _pd##3+ _i) = a4;\
  *( _pd##3+ _i + 1) =a5;\
  *( _pd##7+ _i) = a2;\
  *( _pd##7+ _i + 1) =a3;\
}

/*
   Macro to make the third radix-8 pass and store in a DIT 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_8_pp_pass3_store_B(_t,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r - _t##5##i;\
  a1 = _t##5##i + _t##5##r;\
  a2 = _t##0##r + _t##4##r;\
  a3 = _t##0##i + _t##4##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##4##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##4##i;\
  *( _pd##0 + _i) = a2;\
  *( _pd##0 + _i + 1) =a3;\
  a2 = _t##1##r + a0;\
  a3 = _t##1##i + a1;\
  *( _pd##4 + _i) = a4;\
  *( _pd##4 + _i + 1) =a5;\
  a4 = _t##1##r - a0;\
  a5 = _t##1##i - a1;\
  *( _pd##1 + _i) = a2;\
  *( _pd##1 + _i + 1) =a3;\
  a2 = _t##2##r - _t##6##i;\
  a3 = _t##2##i + _t##6##r;\
  *( _pd##5 + _i) = a4;\
  *( _pd##5 + _i + 1) =a5;\
  a0 = _t##7##r + _t##7##i;\
  a1 = _t##7##i - _t##7##r;\
  *( _pd##2+ _i) = a2;\
  *( _pd##2+ _i + 1) =a3;\
  a4= _t##2##r + _t##6##i;\
  a0 *= F_1_8i;\
  a5= _t##2##i - _t##6##r;\
  a1 *= F_1_8i;\
  *( _pd##6+ _i) = a4;\
  *( _pd##6+ _i + 1) =a5;\
  a4 = _t##3##r + a0;\
  a5 = _t##3##i + a1;\
  a2 = _t##3##r - a0;\
  a3 = _t##3##i - a1;\
  *( _pd##3+ _i) = a4;\
  *( _pd##3+ _i + 1) =a5;\
  *( _pd##7+ _i) = a2;\
  *( _pd##7+ _i + 1) =a3;\
}

/*
   Macro to make the first radix-8 pass in the last DIF step with prefetch 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _f is the tiwddle factors locale prefix name
*/
#define cplx_radix_8_twd_pp_last_dif_pass1(_t,_f,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;\
  a0 = *(_pd + 8);\
  a1 = *(_pd + 9);\
  a2 = a0 * (_f##4##r);\
  prefetch_data(_pd, 24);\
  a3 = a1 * (_f##4##i);\
  a4 = *(_pd );\
  a5 = *(_pd + 1);\
  a6 = a0 * (_f##4##i);\
  prefetch_data(_pd, 16);\
  a7 = a1 * (_f##4##r);\
  a0 = *(_pd + 4);\
  a1 = *(_pd + 5);\
  a2 -= a3;\
  a6 += a7;\
  prefetch_data(_pd, 20);\
  a3 = a0 * (_f##2##r);\
  _t##0##r = a4 + a2;\
  a7 = a1 * (_f##2##i);\
  _t##1##r = a4 - a2;\
  a4 = a0 * (_f##2##i);\
  _t##0##i = a5 + a6;\
  a2 = a1 * (_f##2##r);\
  _t##1##i = a5 - a6;\
  a0 = *(_pd + 12);\
  a1 = *(_pd + 13);\
  a5 = a0 * (_f##6##r);\
  a3 -= a7;\
  a6 = a1 * (_f##6##i);\
  prefetch_data(_pd, 28);\
  a2 += a4;\
  a4 = a0 * (_f##6##i);\
  a5 -= a6;\
  a6 = a1 * (_f##6##r);\
  a0 = *(_pd + 2);\
  a1 = *(_pd + 3);\
  a4 += a6;\
  a6 = a0 * (_f##1##r);\
  _t##2##r = a3 + a5;\
  a7 = a1 * (_f##1##i);\
  _t##3##r = a3 - a5;\
  a3 = a0 * (_f##1##i);\
  _t##2##i = a2 + a4;\
  a5 = a1 * (_f##1##r);\
  _t##3##i = a2 - a4;\
  a0 = *(_pd + 10);\
  a1 = *(_pd + 11);\
  a2 = a0 * (_f##5##r);\
  a6 -= a7;\
  a4 = a1 * (_f##5##i);\
  a3 += a5;\
  a7 = a0 * (_f##5##i);\
  a2 -= a4;\
  a5 = a1 * (_f##5##r);\
  a0 = *(_pd + 6);\
  a1 = *(_pd + 7);\
  a5 += a7;\
  a4 = a0 * (_f##3##r);\
  _t##4##r = a6 + a2;\
  a7 = a1 * (_f##3##i);\
  _t##5##r = a6 - a2;\
  a6 = a0 * (_f##3##i);\
  _t##4##i = a3 + a5;\
  a2 = a1 * (_f##3##r);\
  _t##5##i = a3 - a5;\
  a0 = *(_pd + 14);\
  a1 = *(_pd + 15);\
  a3 = a0 * (_f##7##r);\
  a4 -= a7;\
  a5 = a1 * (_f##7##i);\
  a2 += a6;\
  a6 = a0 * (_f##7##i);\
  a3 -= a5;\
  a5 = a1 * (_f##7##r);\
  _t##6##r = a4 + a3;\
  a6 += a5;\
  _t##7##r = a4 - a3;\
  _t##6##i = a2 + a6;\
  _t##7##i = a2 - a6;\
}

/*
   Macro to make the first radix-8 pass in the last DIF step 
   The prefetch is backward.
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _f is the tiwddle factors locale prefix name
*/
#define cplx_radix_8_twd_pp_last_dif_pass1_down(_t,_f,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;\
  a0 = *(_pd + 8);\
  a1 = *(_pd + 9);\
  a2 = a0 * (_f##4##r);\
  prefetch_data(_pd, (-8));\
  a3 = a1 * (_f##4##i);\
  a4 = *(_pd );\
  a5 = *(_pd + 1);\
  a6 = a0 * (_f##4##i);\
  prefetch_data(_pd, (-16));\
  a2 -= a3;\
  a7 = a1 * (_f##4##r);\
  a0 = *(_pd + 4);\
  a1 = *(_pd + 5);\
  a6 += a7;\
  prefetch_data(_pd, (-12));\
  a3 = a0 * (_f##2##r);\
  _t##0##r = a4 + a2;\
  a7 = a1 * (_f##2##i);\
  _t##1##r = a4 - a2;\
  a4 = a0 * (_f##2##i);\
  _t##0##i = a5 + a6;\
  a2 = a1 * (_f##2##r);\
  _t##1##i = a5 - a6;\
  a0 = *(_pd + 12);\
  a1 = *(_pd + 13);\
  a5 = a0 * (_f##6##r);\
  a3 -= a7;\
  a6 = a1 * (_f##6##i);\
  prefetch_data(_pd, (-4));\
  a2 += a4;\
  a4 = a0 * (_f##6##i);\
  a5 -= a6;\
  a6 = a1 * (_f##6##r);\
  a0 = *(_pd + 2);\
  a1 = *(_pd + 3);\
  a4 += a6;\
  a6 = a0 * (_f##1##r);\
  _t##2##r = a3 + a5;\
  a7 = a1 * (_f##1##i);\
  _t##3##r = a3 - a5;\
  a3 = a0 * (_f##1##i);\
  _t##2##i = a2 + a4;\
  a5 = a1 * (_f##1##r);\
  _t##3##i = a2 - a4;\
  a0 = *(_pd + 10);\
  a1 = *(_pd + 11);\
  a2 = a0 * (_f##5##r);\
  a6 -= a7;\
  a4 = a1 * (_f##5##i);\
  a3 += a5;\
  a7 = a0 * (_f##5##i);\
  a2 -= a4;\
  a5 = a1 * (_f##5##r);\
  a0 = *(_pd + 6);\
  a1 = *(_pd + 7);\
  a5 += a7;\
  a4 = a0 * (_f##3##r);\
  _t##4##r = a6 + a2;\
  a7 = a1 * (_f##3##i);\
  _t##5##r = a6 - a2;\
  a6 = a0 * (_f##3##i);\
  _t##4##i = a3 + a5;\
  a2 = a1 * (_f##3##r);\
  _t##5##i = a3 - a5;\
  a0 = *(_pd + 14);\
  a1 = *(_pd + 15);\
  a3 = a0 * (_f##7##r);\
  a4 -= a7;\
  a5 = a1 * (_f##7##i);\
  a2 += a6;\
  a6 = a0 * (_f##7##i);\
  a3 -= a5;\
  a5 = a1 * (_f##7##r);\
  _t##6##r = a4 + a3;\
  a6 += a5;\
  _t##7##r = a4 - a3;\
  _t##6##i = a2 + a6;\
  _t##7##i = a2 - a6;\
}


/*
   Macro to make the last radix-8 pass in the last DIF step. No store 
   _t is the locale prefix variable name
*/
#define cplx_radix_8_pp_last_dif_pass3(_t)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r + _t##5##i;\
  a1 = _t##5##i - _t##5##r;\
  a2 = _t##0##r + _t##4##r;\
  a3 = _t##0##i + _t##4##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##4##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##4##i;\
  _t##0##r = a2;\
  _t##0##i = a3;\
  a2 = _t##1##r + a0;\
  a3 = _t##1##i + a1;\
  _t##4##r = a4;\
  _t##4##i = a5;\
  a4 = _t##1##r - a0;\
  a5 = _t##1##i - a1;\
  _t##1##r = a2;\
  _t##1##i = a3;\
  a2 = _t##2##r + _t##6##i;\
  a3 = _t##2##i - _t##6##r;\
  _t##5##r = a4;\
  _t##5##i = a5;\
  a0 = _t##7##r - _t##7##i;\
  a1 = _t##7##i + _t##7##r;\
  a4 = _t##2##r - _t##6##i;\
  a5 = _t##2##i + _t##6##r;\
  _t##2##r = a2;\
  _t##2##i = a3;\
  a0 *= F_1_8i;\
  a1 *= F_1_8i;\
  _t##6##r = a4;\
  _t##6##i = a5;\
  _t##7##r = _t##3##r - a0;\
  _t##7##i = _t##3##i - a1;\
  _t##3##r += a0;\
  _t##3##i += a1;\
}

/*
   Macro to make the first radix-8 pass in the first DIT step 
   _t is the locale prefix variable name
*/
#define cplx_radix_8_pp_first_dit_pass1(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##4##r;\
  a1 = _t##4##i;\
  _t##4##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##4##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##5##r;\
  a1 = _t##5##i;\
  _t##5##r = _t##1##r - a0;\
  _t##1##r += a0;\
  _t##5##i = _t##1##i - a1;\
  _t##1##i += a1;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##2##r - a0;\
  _t##2##r += a0;\
  _t##6##i = _t##2##i - a1;\
  _t##2##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##3##r - a0;\
  _t##3##r += a0;\
  _t##7##i = _t##3##i - a1;\
  _t##3##i += a1;\
}

/*
   Macro to make the second radix-8 pass in the first DIT step 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _f is the tiwddle factors locale prefix name
*/
#define cplx_radix_8_pp_first_dit_pass2(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##2##r;\
  a1 = _t##2##i;\
  _t##2##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##2##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##4##r + a1;\
  _t##4##r -= a1;\
  _t##6##i = _t##4##i - a0;\
  _t##4##i += a0;\
  a0 = _t##3##r;\
  a1 = _t##3##i;\
  _t##3##r = _t##1##r - a0;\
  _t##1##r += a0;\
  _t##3##i = _t##1##i - a1;\
  _t##1##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##5##r + a1;\
  _t##5##r -= a1;\
  _t##7##i = _t##5##i - a0;\
  _t##5##i += a0;\
}

/*
   Macro to make the third radix-8 pass in the first DIT pass 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
*/
#define cplx_radix_8_pp_first_dit_pass3(_t,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r - _t##5##i;\
  a1 = _t##5##i + _t##5##r;\
  a2 = _t##0##r + _t##1##r;\
  a3 = _t##0##i + _t##1##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##1##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##1##i;\
  *( _pd ) = a2;\
  *( _pd + 1) =a3;\
  a2 = _t##4##r + a0;\
  a3 = _t##4##i + a1;\
  *( _pd + 8) = a4;\
  *( _pd + 9) = a5;\
  a4 = _t##4##r - a0;\
  a5 = _t##4##i - a1;\
  *( _pd + 2) = a2;\
  *( _pd + 3) = a3;\
  a2 = _t##2##r - _t##3##i;\
  a3 = _t##2##i + _t##3##r;\
  *( _pd + 10) = a4;\
  *( _pd + 11) = a5;\
  a0 = _t##7##r + _t##7##i;\
  a1 = _t##7##i - _t##7##r;\
  *( _pd + 4) = a2;\
  *( _pd + 5) = a3;\
  a4= _t##2##r + _t##3##i;\
  a0 *= F_1_8i;\
  a5= _t##2##i - _t##3##r;\
  a1 *= F_1_8i;\
  *( _pd + 12) = a4;\
  *( _pd + 13) = a5;\
  a4 = _t##6##r + a0;\
  a5 = _t##6##i + a1;\
  a2 = _t##6##r - a0;\
  a3 = _t##6##i - a1;\
  *( _pd + 6) = a4;\
  *( _pd + 7) = a5;\
  *( _pd + 14) = a2;\
  *( _pd + 15) = a3;\
}

#define cplx_radix_8_twd_pp_last_dit_pass1(_t,_f,_i,_pd)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5,a6,a7;\
  a0 = *(_pd##4 + _i);\
  a1 = *(_pd##4 + _i + 1);\
  a2 = a0 * (_f##4##r);\
  prefetch_p(_pd##4 + _i);\
  a3 = a1 * (_f##4##i);\
  a4 = *(_pd##0 + _i);\
  a5 = *(_pd##0 + _i + 1);\
  a6 = a0 * (_f##4##i);\
  prefetch_p(_pd##0 + _i);\
  a2 -= a3;\
  a7 = a1 * (_f##4##r);\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a6 += a7;\
  prefetch_p(_pd##2 + _i);\
  a3 = a0 * (_f##2##r);\
  _t##0##r = a4 + a2;\
  a7 = a1 * (_f##2##i);\
  _t##1##r = a4 - a2;\
  a4 = a0 * (_f##2##i);\
  _t##0##i = a5 + a6;\
  a2 = a1 * (_f##2##r);\
  _t##1##i = a5 - a6;\
  a0 = *(_pd##6 + _i);\
  a1 = *(_pd##6 + _i + 1);\
  a5 = a0 * (_f##6##r);\
  a3 -= a7;\
  a6 = a1 * (_f##6##i);\
  prefetch_p(_pd##6 + _i);\
  a2 += a4;\
  a4 = a0 * (_f##6##i);\
  a5 -= a6;\
  a6 = a1 * (_f##6##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  prefetch_p(_pd##1 + _i);\
  a4 += a6;\
  a6 = a0 * (_f##1##r);\
  _t##2##r = a3 + a5;\
  a7 = a1 * (_f##1##i);\
  _t##3##r = a3 - a5;\
  a3 = a0 * (_f##1##i);\
  _t##2##i = a2 + a4;\
  a5 = a1 * (_f##1##r);\
  _t##3##i = a2 - a4;\
  a0 = *(_pd##5 + _i);\
  a1 = *(_pd##5 + _i + 1);\
  prefetch_p(_pd##5 + _i);\
  a2 = a0 * (_f##5##r);\
  a6 -= a7;\
  a4 = a1 * (_f##5##i);\
  a3 += a5;\
  a7 = a0 * (_f##5##i);\
  a2 -= a4;\
  a5 = a1 * (_f##5##r);\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  prefetch_p(_pd##3 + _i);\
  a5 += a7;\
  a4 = a0 * (_f##3##r);\
  _t##4##r = a6 + a2;\
  a7 = a1 * (_f##3##i);\
  _t##5##r = a6 - a2;\
  a6 = a0 * (_f##3##i);\
  _t##4##i = a3 + a5;\
  a2 = a1 * (_f##3##r);\
  _t##5##i = a3 - a5;\
  a0 = *(_pd##7 + _i);\
  a1 = *(_pd##7 + _i + 1);\
  prefetch_p(_pd##7 + _i);\
  a3 = a0 * (_f##7##r);\
  a4 -= a7;\
  a5 = a1 * (_f##7##i);\
  a2 += a6;\
  a6 = a0 * (_f##7##i);\
  a3 -= a5;\
  a5 = a1 * (_f##7##r);\
  _t##6##r = a4 + a3;\
  a6 += a5;\
  _t##7##r = a4 - a3;\
  _t##6##i = a2 + a6;\
  _t##7##i = a2 - a6;\
}

#define cplx_radix_8_pp_last_dit_pass3(_t)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r - _t##5##i;\
  a1 = _t##5##i + _t##5##r;\
  a2 = _t##0##r + _t##4##r;\
  a3 = _t##0##i + _t##4##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##4##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##4##i;\
  _t##0##r = a2;\
  _t##0##i = a3;\
  a2 = _t##1##r + a0;\
  a3 = _t##1##i + a1;\
  _t##4##r = a4;\
  _t##4##i = a5;\
  a4 = _t##1##r - a0;\
  a5 = _t##1##i - a1;\
  _t##1##r = a2;\
  _t##1##i = a3;\
  a2 = _t##2##r - _t##6##i;\
  a3 = _t##2##i + _t##6##r;\
  _t##5##r = a4;\
  _t##5##i = a5;\
  a0 = _t##7##r + _t##7##i;\
  a1 = _t##7##i - _t##7##r;\
  a4 = _t##2##r + _t##6##i;\
  a5 = _t##2##i - _t##6##r;\
  _t##2##r = a2;\
  _t##2##i = a3;\
  a0 *= F_1_8i;\
  a1 *= F_1_8i;\
  _t##6##r = a4;\
  _t##6##i = a5;\
  _t##7##r = _t##3##r - a0;\
  _t##7##i = _t##3##i - a1;\
  _t##3##r += a0;\
  _t##3##i += a1;\
}

#define cplx_radix_8_pp_first_dif_pass1(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##4##r;\
  a1 = _t##4##i;\
  _t##4##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##4##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##5##r;\
  a1 = _t##5##i;\
  _t##5##r = _t##1##r - a0;\
  _t##1##r += a0;\
  _t##5##i = _t##1##i - a1;\
  _t##1##i += a1;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##2##r - a0;\
  _t##2##r += a0;\
  _t##6##i = _t##2##i - a1;\
  _t##2##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##3##r - a0;\
  _t##3##r += a0;\
  _t##7##i = _t##3##i - a1;\
  _t##3##i += a1;\
}

#define cplx_radix_8_pp_first_dif_pass2(_t)\
{\
  y_limb_t a0,a1;\
  a0 = _t##2##r;\
  a1 = _t##2##i;\
  _t##2##r = _t##0##r - a0;\
  _t##0##r += a0;\
  _t##2##i = _t##0##i - a1;\
  _t##0##i += a1;\
  a0 = _t##6##r;\
  a1 = _t##6##i;\
  _t##6##r = _t##4##r - a1;\
  _t##4##r += a1;\
  _t##6##i = _t##4##i + a0;\
  _t##4##i -= a0;\
  a0 = _t##3##r;\
  a1 = _t##3##i;\
  _t##3##r = _t##1##r - a0;\
  _t##1##r += a0;\
  _t##3##i = _t##1##i - a1;\
  _t##1##i += a1;\
  a0 = _t##7##r;\
  a1 = _t##7##i;\
  _t##7##r = _t##5##r - a1;\
  _t##5##r += a1;\
  _t##7##i = _t##5##i + a0;\
  _t##5##i -= a0;\
}

#define cplx_radix_8_pp_first_dif_pass3(_t,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3,a4,a5;\
  a0 = _t##5##r + _t##5##i;\
  a1 = _t##5##i - _t##5##r;\
  a2 = _t##0##r + _t##1##r;\
  a3 = _t##0##i + _t##1##i;\
  a0 *= F_1_8r;\
  a4 = _t##0##r - _t##1##r;\
  a1 *= F_1_8r;\
  a5 = _t##0##i - _t##1##i;\
  *( _pd##0 + _i ) = a2;\
  *( _pd##0 + _i + 1) =a3;\
  a2 = _t##4##r + a0;\
  a3 = _t##4##i + a1;\
  *( _pd##4 + _i) = a4;\
  *( _pd##4 + _i + 1) = a5;\
  a4 = _t##4##r - a0;\
  a5 = _t##4##i - a1;\
  *( _pd##1 + _i) = a2;\
  *( _pd##1 + _i + 1) = a3;\
  a2 = _t##2##r + _t##3##i;\
  a3 = _t##2##i - _t##3##r;\
  *( _pd##5 + _i) = a4;\
  *( _pd##5 + _i + 1) = a5;\
  a0 = _t##7##r - _t##7##i;\
  a1 = _t##7##i + _t##7##r;\
  *( _pd##2 + _i) = a2;\
  *( _pd##2 + _i  + 1) = a3;\
  a4= _t##2##r - _t##3##i;\
  a0 *= F_1_8i;\
  a5= _t##2##i + _t##3##r;\
  a1 *= F_1_8i;\
  *( _pd##6 + _i) = a4;\
  *( _pd##6 + _i + 1) = a5;\
  a4 = _t##6##r + a0;\
  a5 = _t##6##i + a1;\
  a2 = _t##6##r - a0;\
  a3 = _t##6##i - a1;\
  *( _pd##3 + _i) = a4;\
  *( _pd##3 + _i + 1) = a5;\
  *( _pd##7 + _i) = a2;\
  *( _pd##7 + _i + 1) = a3;\
}
