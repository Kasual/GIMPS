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
   Macro to make a radix-4 pass with twiddle factors and prefetch data 
   in a DIF 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_twd_pp_passes_F(_t,_f,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = a0 * (_f##2##r);\
  prefetch_p(_pd##2 + _i);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##3##r = a0 * (_f##2##i);\
  prefetch_p(_pd##0 + _i);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  _t##3##r += _t##3##i;\
  prefetch_p(_pd##1 + _i);\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  prefetch_p(_pd##3 + _i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i + _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  a2 = _t##1##i - _t##3##r;\
  a1 = _t##1##r - _t##3##i;\
  a0 = _t##1##r + _t##3##i;\
  *(_pd##3 + _i) = a1;\
  *(_pd##3 + _i + 1) =a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  *(_pd##1 + _i) = a0;\
  *(_pd##1 + _i + 1) =a2;\
  a0 = _t##0##r + _t##2##r;\
  a2 = _t##0##i + _t##2##i;\
  *(_pd##2 + _i) = a1;\
  *(_pd##2 + _i + 1) =a3;\
  *(_pd##0 + _i) = a0;\
  *(_pd##0 + _i + 1) =a2;\
}


/*
   Macro to make a radix-4 pass with twiddle factors without prefetch data in 
   a DIF 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointer prefix variable name
   _i is the index
*/
#define cplx_radix_4_twd_pp_passes_no_fetch_F(_t,_f,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = a0 * (_f##2##r);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##3##r = a0 * (_f##2##i);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  _t##3##r += _t##3##i;\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i + _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  a2 = _t##1##i - _t##3##r;\
  a1 = _t##1##r - _t##3##i;\
  a0 = _t##1##r + _t##3##i;\
  *(_pd##3 + _i) = a1;\
  *(_pd##3 + _i + 1) =a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  *(_pd##1 + _i) = a0;\
  *(_pd##1 + _i + 1) =a2;\
  a0 = _t##0##r + _t##2##r;\
  a2 = _t##0##i + _t##2##i;\
  *(_pd##2 + _i) = a1;\
  *(_pd##2 + _i + 1) =a3;\
  *(_pd##0 + _i) = a0;\
  *(_pd##0 + _i + 1) =a2;\
}


/*
   Macro to make a radix-4 pass with twiddle factors and prefetch data 
   in a DIT 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointer prefix variable name
   _i is the index
*/
#define cplx_radix_4_twd_pp_passes_B(_t,_f,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = a0 * (_f##2##r);\
  prefetch_p(_pd##2 + _i);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##3##r = a0 * (_f##2##i);\
  prefetch_p(_pd##0 + _i);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  _t##3##r += _t##3##i;\
  prefetch_p(_pd##1 + _i);\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  prefetch_p(_pd##3 + _i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i - _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  a2 = _t##1##i + _t##3##r;\
  a1 = _t##1##r + _t##3##i;\
  a0 = _t##1##r - _t##3##i;\
  *(_pd##3 + _i) = a1;\
  *(_pd##3 + _i + 1) =a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  *(_pd##1 + _i) = a0;\
  *(_pd##1 + _i + 1) =a2;\
  a0 = _t##0##r + _t##2##r;\
  a2 = _t##0##i + _t##2##i;\
  *(_pd##2 + _i) = a1;\
  *(_pd##2 + _i + 1) =a3;\
  *(_pd##0 + _i) = a0;\
  *(_pd##0 + _i + 1) =a2;\
}


/*
   Macro to make a radix-4 pass with twiddle factors without prefetch data in 
   a DIT 
   _t is the locale prefix variable name
   _f is the locale twidle factors prefix variable name
   _pd is the pointer prefix variable name
   _i is the index
*/
#define cplx_radix_4_twd_pp_passes_no_fetch_B(_t,_f,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = a0 * (_f##2##r);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##3##r = a0 * (_f##2##i);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  _t##3##r += _t##3##i;\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i - _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  a2 = _t##1##i + _t##3##r;\
  a1 = _t##1##r + _t##3##i;\
  a0 = _t##1##r - _t##3##i;\
  *(_pd##3 + _i) = a1;\
  *(_pd##3 + _i + 1) =a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  *(_pd##1 + _i) = a0;\
  *(_pd##1 + _i + 1) =a2;\
  a0 = _t##0##r + _t##2##r;\
  a2 = _t##0##i + _t##2##i;\
  *(_pd##2 + _i) = a1;\
  *(_pd##2 + _i + 1) =a3;\
  *(_pd##0 + _i) = a0;\
  *(_pd##0 + _i + 1) =a2;\
}


/*
   Macro to make a radix-4 pass with prefetch data in a DIF 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_notwd_pp_passes_F(_t,_pd,_i)\
{\
  y_limb_t a0,a1;\
  _t##3##r = *(_pd##2 + _i);\
  _t##3##i = *(_pd##2 + _i + 1);\
  prefetch_p (_pd##2 + _i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  prefetch_p (_pd##0 + _i);\
  _t##0##r =  _t##2##r + _t##3##r;\
  _t##1##r =  _t##2##r - _t##3##r;\
  _t##0##i =  _t##2##i + _t##3##i;\
  _t##1##i =  _t##2##i - _t##3##i;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##3 + _i);\
  prefetch_p (_pd##3 + _i);\
  _t##2##r = a0 + a1;\
  _t##3##r = a0 - a1;\
  a1 = *(_pd##3 + _i + 1);\
  a0 = *(_pd##1 + _i + 1);\
  prefetch_p (_pd##1 + _i);\
  _t##2##i = a0 + a1;\
  _t##3##i = a0 - a1;\
  \
  a0 = _t##0##r - _t##2##r;\
  a1 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  *( _pd##2 + _i) = a0;\
  *( _pd##2 + _i + 1) = a1;\
  a0 = _t##1##r - _t##3##i;\
  a1 = _t##1##i + _t##3##r;\
  t##1##r += _t##3##i;\
  t##1##i -= _t##3##r;\
  *( _pd##0 + _i) = _t##0##r;\
  *( _pd##0 + _i + 1) = _t##0##i;\
  *( _pd##3 + _i) = a0;\
  *( _pd##3 + _i + 1) = a1;\
  *( _pd##1 + _i) = _t##1##r;\
  *( _pd##1 + _i + 1) = _t##1##i;\
}


/*
   Macro to make a radix-4 pass with prefetch data in a DIT 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_notwd_pp_passes_B(_t,_pd,_i)\
{\
  y_limb_t a0,a1;\
  _t##3##r = *(_pd##2 + _i);\
  _t##3##i = *(_pd##2 + _i + 1);\
  prefetch_p (_pd##2 + _i);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  prefetch_p (_pd##0 + _i);\
  _t##0##r =  _t##2##r + _t##3##r;\
  _t##1##r =  _t##2##r - _t##3##r;\
  _t##0##i =  _t##2##i + _t##3##i;\
  _t##1##i =  _t##2##i - _t##3##i;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##3 + _i);\
  prefetch_p (_pd##3 + _i);\
  _t##2##r = a0 + a1;\
  _t##3##r = a0 - a1;\
  a1 = *(_pd##3 + _i + 1);\
  a0 = *(_pd##1 + _i + 1);\
  prefetch_p (_pd##1 + _i);\
  _t##2##i = a0 + a1;\
  _t##3##i = a0 - a1;\
  \
  a0 = _t##0##r - _t##2##r;\
  a1 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  *( _pd##2 + _i) = a0;\
  *( _pd##2 + _i + 1) = a1;\
  a0 = _t##1##r + _t##3##i;\
  a1 = _t##1##i - _t##3##r;\
  t##1##r -= _t##3##i;\
  t##1##i += _t##3##r;\
  *( _pd##0 + _i) = _t##0##r;\
  *( _pd##0 + _i + 1) = _t##0##i;\
  *( _pd##3 + _i) = a0;\
  *( _pd##3 + _i + 1) = a1;\
  *( _pd##1 + _i) = _t##1##r;\
  *( _pd##1 + _i + 1) = _t##1##i;\
}

/*
   Macro to make a radix-4 pass without prefetch data in a DIF 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_notwd_pp_passes_no_fetch_F(_t,_pd,_i)\
{\
  y_limb_t a0,a1;\
  _t##3##r = *(_pd##2 + _i);\
  _t##3##i = *(_pd##2 + _i + 1);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##0##r =  _t##2##r + _t##3##r;\
  _t##1##r =  _t##2##r - _t##3##r;\
  _t##0##i =  _t##2##i + _t##3##i;\
  _t##1##i =  _t##2##i - _t##3##i;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##3 + _i);\
  _t##2##r = a0 + a1;\
  _t##3##r = a0 - a1;\
  a1 = *(_pd##3 + _i + 1);\
  a0 = *(_pd##1 + _i + 1);\
  _t##2##i = a0 + a1;\
  _t##3##i = a0 - a1;\
  \
  a0 = _t##0##r - _t##2##r;\
  a1 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  *( _pd##2 + _i) = a0;\
  *( _pd##2 + _i + 1) = a1;\
  a0 = _t##1##r - _t##3##i;\
  a1 = _t##1##i + _t##3##r;\
  t##1##r += _t##3##i;\
  t##1##i -= _t##3##r;\
  *( _pd##0 + _i) = _t##0##r;\
  *( _pd##0 + _i + 1) = _t##0##i;\
  *( _pd##3 + _i) = a0;\
  *( _pd##3 + _i + 1) = a1;\
  *( _pd##1 + _i) = _t##1##r;\
  *( _pd##1 + _i + 1) = _t##1##i;\
}


/*
   Macro to make a radix-4 pass without prefetch data in a DIF 
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_notwd_pp_passes_no_fetch_B(_t,_pd,_i)\
{\
  y_limb_t a0,a1;\
  _t##3##r = *(_pd##2 + _i);\
  _t##3##i = *(_pd##2 + _i + 1);\
  _t##2##r = *(_pd##0 + _i);\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##0##r =  _t##2##r + _t##3##r;\
  _t##1##r =  _t##2##r - _t##3##r;\
  _t##0##i =  _t##2##i + _t##3##i;\
  _t##1##i =  _t##2##i - _t##3##i;\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##3 + _i);\
  _t##2##r = a0 + a1;\
  _t##3##r = a0 - a1;\
  a1 = *(_pd##3 + _i + 1);\
  a0 = *(_pd##1 + _i + 1);\
  _t##2##i = a0 + a1;\
  _t##3##i = a0 - a1;\
  \
  a0 = _t##0##r - _t##2##r;\
  a1 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  *( _pd##2 + _i) = a0;\
  *( _pd##2 + _i + 1) = a1;\
  a0 = _t##1##r + _t##3##i;\
  a1 = _t##1##i - _t##3##r;\
  t##1##r -= _t##3##i;\
  t##1##i += _t##3##r;\
  *( _pd##0 + _i) = _t##0##r;\
  *( _pd##0 + _i + 1) = _t##0##i;\
  *( _pd##3 + _i) = a0;\
  *( _pd##3 + _i + 1) = a1;\
  *( _pd##1 + _i) = _t##1##r;\
  *( _pd##1 + _i + 1) = _t##1##i;\
}

/*
   Macro to make the last radix-4 DIF pass with prefetch data 
   _t is the locale prefix variable name
   _f is the locale prefix twiddle factor name
   _pd is the pointers prefix variable name
*/
#define cplx_radix_4_pp_last_dif_passes(_t,_f,_pd)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd + 4);\
  a1 = *(_pd + 5);\
  a2 = a0 * (_f##2##r);\
  prefetch_data(_pd , 12);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd );\
  _t##2##i = *(_pd + 1);\
  _t##3##r = a0 * (_f##2##i);\
  prefetch_data(_pd , 8);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd + 2);\
  a1 = *(_pd + 3);\
  _t##3##r += _t##3##i;\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd + 6);\
  a1 = *(_pd + 7);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i + _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  _t##1##i -= _t##3##r;\
  _t##3##r = _t##1##r - _t##3##i;\
  _t##1##r += _t##3##i;\
  _t##3##i = a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  _t##2##r = a1;\
  _t##2##i = a3;\
}

/*
   Macro to make the last radix-4 DIF pass with backward prefetch data 
   _t is the locale prefix variable name
   _f is the locale prefix twiddle factor name
   _pd is the pointers prefix variable name
*/
#define cplx_radix_4_pp_last_dif_passes_down(_t,_f,_pd)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd + 4);\
  a1 = *(_pd + 5);\
  a2 = a0 * (_f##2##r);\
  prefetch_data(_pd , (-4));\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd );\
  _t##2##i = *(_pd + 1);\
  _t##3##r = a0 * (_f##2##i);\
  prefetch_data(_pd , (-8));\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd + 2);\
  a1 = *(_pd + 3);\
  _t##3##r += _t##3##i;\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd + 6);\
  a1 = *(_pd + 7);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i + _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  _t##1##i -= _t##3##r;\
  _t##3##r = _t##1##r - _t##3##i;\
  _t##1##r += _t##3##i;\
  _t##3##i = a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  _t##2##r = a1;\
  _t##2##i = a3;\
}

/*
   Macro to make the first radix-4 DIT pass  
   _t is the locale prefix variable name
   _pd is the pointers prefix variable name
*/

#define cplx_radix_4_pp_first_dit_passes(_t,_pd)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = _t##0##r - _t##2##r;\
  a1 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  a2 = _t##1##r - _t##3##r;\
  a3 = _t##1##i - _t##3##i;\
  _t##1##r += _t##3##r;\
  _t##1##i += _t##3##i;\
  _t##2##r = a0 - a3;\
  _t##2##i = a1 + a2;\
  _t##3##r = a0 + a3;\
  _t##3##i = a1 - a2;\
  a0 = _t##0##r - _t##1##r;\
  a1 = _t##0##i - _t##1##i;\
  _t##0##r += _t##1##r;\
  _t##0##i += _t##1##i;\
  *( _pd + 2) = _t##2##r;\
  *( _pd + 3) = _t##2##i;\
  *( _pd + 6) = _t##3##r;\
  *( _pd + 7) = _t##3##i;\
  *( _pd ) = _t##0##r;\
  *( _pd + 1) = _t##0##i;\
  *( _pd + 4) = a0;\
  *( _pd + 5) = a1;\
}

/*
   Macro to make the last radix-4 DIT pass with prefetch data 
   _t is the locale prefix variable name
   _f is the locale prefix twiddle factor name
   _pd is the pointers prefix variable name
   _i is the index
*/
#define cplx_radix_4_pp_last_dit_passes(_t,_f,_pd,_i)\
{\
  y_limb_t a0,a1,a2,a3;\
  a0 = *(_pd##2 + _i);\
  a1 = *(_pd##2 + _i + 1);\
  a2 = a0 * (_f##2##r);\
  prefetch_p(_pd##2 + _i);\
  a3 = a1 * (_f##2##i);\
  _t##2##r = *(_pd##0 + _i );\
  _t##2##i = *(_pd##0 + _i + 1);\
  _t##3##r = a0 * (_f##2##i);\
  prefetch_p(_pd##0 + _i);\
  a2 -= a3;\
  _t##3##i = a1 * (_f##2##r);\
  a0 = *(_pd##1 + _i);\
  a1 = *(_pd##1 + _i + 1);\
  prefetch_p(_pd##1 + _i);\
  _t##3##r += _t##3##i;\
  a3 = a0 * (_f##1##r);\
  _t##0##r = _t##2##r + a2;\
  _t##3##i = a1 * (_f##1##i);\
  _t##1##r = _t##2##r - a2;\
  _t##2##r = a0 * (_f##1##i);\
  _t##0##i = _t##2##i + _t##3##r;\
  a2 = a1 * (_f##1##r);\
  _t##1##i = _t##2##i - _t##3##r;\
  a0 = *(_pd##3 + _i);\
  a1 = *(_pd##3 + _i + 1);\
  prefetch_p(_pd##3 + _i);\
  _t##2##i = a0 * (_f##3##r);\
  a3 -= _t##3##i;\
  _t##3##r = a1 * (_f##3##i);\
  a2 += _t##2##r;\
  a0 *= (_f##3##i);\
  _t##2##i -= _t##3##r;\
  a1 *= (_f##3##r);\
  _t##2##r = a3 + _t##2##i;\
  _t##3##r = a3 - _t##2##i;\
  a0 += a1;\
  \
  a3 = _t##1##i - _t##3##r;\
  _t##3##i = a2 - a0;\
  _t##2##i = a2 + a0;\
  _t##1##i += _t##3##r;\
  _t##3##r = _t##1##r + _t##3##i;\
  _t##1##r -= _t##3##i;\
  _t##3##i = a3;\
  a1 = _t##0##r - _t##2##r;\
  a3 = _t##0##i - _t##2##i;\
  _t##0##r += _t##2##r;\
  _t##0##i += _t##2##i;\
  _t##2##r = a1;\
  _t##2##i = a3;\
}

