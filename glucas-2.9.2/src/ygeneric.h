/*$Id$*/
/*  This file is part of
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
/* This are macros to be used with pseudo complex FFT                       */

/***************  SINGLE LOAD AND STORE MACROS ***************/


/* to read some trigonomeric factors */
#define cplx_data_trig_load_2(_t0,_t1,_p) \
   _t0##r = *( _p );\
   _t0##i = *( _p + 1);\
   _t1##r = *( _p + 2);\
   _t1##i = *( _p + 3);

#define cplx_data_trig_load_3(_t0,_t1,_t2,_p) \
   _t0##r = *( _p );\
   _t0##i = *( _p + 1);\
   _t1##r = *( _p + 2);\
   _t1##i = *( _p + 3);\
   _t2##r = *( _p + 4);\
   _t2##i = *( _p + 5);


/* To load from real padded array of data as a complex array */
#define cplx_data_to_local(_t,_array,_k) \
	_t##r = _array[addr(( _k )<<1)]; \
	_t##i = _array[addr(( _k )<<1)+1];


#define cplx_data_to_local_p(_t,_pd,_array,_k) \
	_pd = _array + addr(( _k)<<1);\
	_t##r = *(_pd ); \
	_t##i = *( _pd + 1);\
        prefetch_p(_pd);


#define cplx_data_to_local_pp(_t,_index,_pd) \
	_t##r = *( _pd + _index); \
	_t##i = *( _pd + _index + 1);\
        prefetch_p( _pd + _index);


/* To store on read padded array from local data as complex */

#define cplx_local_to_data(_array,_k,_t) \
	_array[addr(( _k )<<1)]=_t##r ; \
	_array[addr(( _k )<<1)+1]=_t##i;


#define cplx_local_to_data_p(_pd,_t) \
	*( _pd ) = _t##r ; \
	*( _pd + 1) = _t##i;\
        prefetch_p(_pd);


/* To load from real non-padded array of data as a complex array */
#define cplx_mem_to_local(_t,_array,_k) \
	_t##r =_array[( _k )<<1]; \
	_t##i =_array[(( _k )<<1)+1];

#define cplx_mem_to_local_p(_t,_pd) \
	_t##r =  *(_pd);\
	_t##i =  *(_pd + 1); \
         prefetch_p (_pd);


/* To strore to a non-padded array from local */
#define cplx_local_to_mem(_array,_k,_t) \
	_array[( _k )<<1]=_t##r; \
	_array[(( _k )<<1)+1]=_t##i;


/**************** PSEUDO_COMPLEX GENERAL AUXILIAR MACROS **************/
/**********************************************************************
   cplx_muladdsub macro :
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex and _t is temporal aux var.
***********************************************************************/

# define cplx_muladdsub(_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai,_a=_t0##r; \
		_ar = (_t1##r ) * (_f##r ) - (_t1##i ) * (_f##i ); \
		_ai = (_t1##r ) * (_f##i ) + (_t1##i ) * (_f##r ); \
		_t0##r += _ar;\
		_t1##r = _a - _ar;\
		_a = _t0##i;\
		_t0##i += _ai;\
		_t1##i = _a - _ai;\
	}


/**********************************************************************
   cplx_load_muladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
 where all the variables are pseudo-complex and _t is temporal aux var.
**********************************************************************/
#define cplx_load_muladdsub(_t0,_t1,_f,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_ar = (_br) * (_f##r) - (_bi) * (_f##i); \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r); \
		_t0##r = _array[addr(_k0<<1)] + _ar; \
		_t1##r = _array[addr(_k0<<1)] - _ar; \
		_t0##i = _array[addr(_k0<<1)+1] + _ai; \
		_t1##i = _array[addr(_k0<<1)+1] - _ai; \
                prefetch_p(&_array[addr(_k0<<1)]);\
	}


#define cplx_load_muladdsub_p(_t0,_t1,_pd0,_pd1,_f,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_pd1= _array + addr((_k1)<<1);\
		_br= *( _pd1 );\
		_bi= *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_ar = (_br) * (_f##r) - (_bi) * (_f##i); \
		_pd0 = _array + addr((_k0)<<1);\
		_ai = (_br) * (_f##i) + (_bi) * (_f##r); \
		_t0##r = *( _pd0 ) + _ar; \
		_t1##r = *( _pd0 ) - _ar; \
		_t0##i = *( _pd0 + 1) + _ai; \
		_t1##i = *( _pd0 + 1) - _ai; \
                prefetch_p(_pd0);\
	}

#define cplx_load_muladdsub_pp(_t0,_t1,_f,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_ar = (_br) * (_f##r) - (_bi) * (_f##i); \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r); \
	        _t0##r = *( _pd0 + _index) + _ar; \
		_t1##r = *( _pd0 + _index) - _ar; \
		_t0##i = *( _pd0 + _index + 1) + _ai; \
		_t1##i = *( _pd0 + _index + 1) - _ai; \
                prefetch_p( _pd0 + _index);\
	}

#define cplx_load_muladdsub_pp_ind(_t0,_t1,_pf,_i,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_ar = (_br) * (_pf[ 2 * _i - 2]) - (_bi) * (_pf[2 * _i - 1]); \
		_ai = (_br) * (_pf[ 2 * _i - 1]) + (_bi) * (_pf[2 * _i - 2]); \
	        _t0##r = *( _pd0 + _index) + _ar; \
		_t1##r = *( _pd0 + _index) - _ar; \
		_t0##i = *( _pd0 + _index + 1) + _ai; \
		_t1##i = *( _pd0 + _index + 1) - _ai; \
                prefetch_p( _pd0 + _index);\
	}

#define cplx_load_muladdsub_pp_no_fetch(_t0,_t1,_f,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
		_ar = (_br) * (_f##r) - (_bi) * (_f##i); \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r); \
		_t0##r = *( _pd0 + _index) + _ar; \
		_t1##r = *( _pd0 + _index) - _ar; \
		_t0##i = *( _pd0 + _index + 1) + _ai; \
		_t1##i = *( _pd0 + _index + 1) - _ai; \
	}

#define cplx_preload_muladdsub(_t0,_t1,_f,_index,_pd0,_pd1,_tp0,_tp1) \
	{                                                             \
		y_limb_t _ar,_ai;                                     \
		_ar = (_tp1##r) * (_f##r) - (_tp1##i) * (_f##i);      \
		_ai = (_tp1##r) * (_f##i) + (_tp1##i) * (_f##r);      \
		_t0##r = _tp0##r + _ar;                               \
		_t1##r = _tp0##r - _ar;                               \
		_t0##i = _tp0##i + _ai;                               \
		_t1##i = _tp0##i - _ai;                               \
                _tp0##r = *( _pd0 + _index);                          \
                _tp0##i = *( _pd0 + _index + 1);                      \
                _tp1##r = *( _pd1 + _index );                         \
                _tp1##i = *( _pd1 + _index + 1);                      \
	}


/**********************************************************************
   cplx_muladdsub_store macro:
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
		_array( _k0) = _t0;
		_array( _k1) = _t1;
   where all the variables are pseudo-complex and _t is temporal aux var.
**********************************************************************/
#define cplx_muladdsub_store(_array,_k0,_k1,_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i); \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r); \
		_array[addr(_k0<<1)]= _t0##r + _ar; \
                _array[addr(_k1<<1)]= _t0##r - _ar; \
		_array[addr(_k0<<1)+1]= _t0##i + _ai; \
                _array[addr(_k1<<1)+1]= _t0##i - _ai; \
        }


#define cplx_muladdsub_store_p(_index,_pd0,_pd1,_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i); \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r); \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
	}




/**********************************************************************
   cplx_mulmuladdsub macro:
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/
# define cplx_mulmuladdsub(_t0,_t1,_f0,_f1) \
	{ \
	        y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##r) * (_f1##r) - (_t1##i) * (_f1##i); \
		_br = (_t0##r) * (_f0##r) - (_t0##i) * (_f0##i); \
		_ai = (_t1##r) * (_f1##i) + (_t1##i) * (_f1##r); \
		_bi = (_t0##r) * (_f0##i) + (_t0##i) * (_f0##r); \
                _t0##r = _br + _ar;\
                _t1##r = _br - _ar;\
                _t0##i = _bi + _ai;\
		_t1##i = _bi - _ai;\
	}


#define cplx_mulmul(_t0,_t1,_f0,_f1)                                \
	{                                                           \
		y_limb_t _ar,_br;                                   \
		_ar = (_t1##r) * (_f1##r) - (_t1##i) * (_f1##i);    \
		_br = (_t0##r) * (_f0##r) - (_t0##i) * (_f0##i);    \
		_t1##i = (_t1##r) * (_f1##i) + (_t1##i) * (_f1##r); \
		_t0##i = (_t0##r) * (_f0##i) + (_t0##i) * (_f0##r); \
		_t1##r = _ar;\
		_t0##r = _br;\
	}

/* The factors here are supossed trig factors */
#define cplx_divdiv(_t0,_t1,_f0,_f1)                                \
	{                                                           \
		y_limb_t _ar,_br;                                   \
		_ar = (_t1##r) * (_f1##r) + (_t1##i) * (_f1##i);    \
		_br = (_t0##r) * (_f0##r) + (_t0##i) * (_f0##i);    \
		_t1##i = (_t1##i) * (_f1##r) - (_t1##r) * (_f1##i); \
		_t0##i = (_t0##i) * (_f0##r) - (_t0##r) * (_f0##i); \
		_t1##r = _ar;\
		_t0##r = _br;\
	}


/**********************************************************************
   cplx_load_mulmuladdsub macro:
		_t0 = _array( _k0);
		_t1 = _array( _k1);
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/
#define cplx_load_mulmuladdsub(_t0,_t1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= _array[addr(_k0<<1)];\
		_bi= _array[addr(_k0<<1)+1];\
                prefetch_p(&_array[addr(_k0<<1)]);\
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i); \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r); \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i); \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r); \
		_t0##r = _ar + _cr; \
		_t1##r = _ar - _cr; \
		_t0##i = _ai + _ci; \
		_t1##i = _ai - _ci; \
	}


# define cplx_load_mulmuladdsub_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_pd0 = _array + addr((_k0)<<1);\
		_br= *( _pd0 );\
		_bi= *( _pd0 + 1);\
                prefetch_p(_pd0);\
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i); \
		_pd1 = _array + addr((_k1)<<1);\
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r); \
		_br= *( _pd1 );\
		_bi= *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i); \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r); \
		_t0##r = _ar + _cr; \
		_t1##r = _ar - _cr; \
		_t0##i = _ai + _ci; \
		_t1##i = _ai - _ci; \
	}

# define cplx_load_mulmuladdsub_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= *( _pd0 + _index);\
		_bi= *( _pd0 + _index + 1);\
                prefetch_p( _pd0 + _index);\
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i); \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r); \
		_br= *( _pd1 + _index);\
		_bi= *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i); \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r); \
		_t0##r = _ar + _cr; \
		_t1##r = _ar - _cr; \
		_t0##i = _ai + _ci; \
		_t1##i = _ai - _ci; \
	}

# define cplx_load_mulmuladdsub_pp_ind(_t0,_t1,_pf,_i0,_i1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= *( _pd0 + _index);\
		_bi= *( _pd0 + _index + 1);\
                prefetch_p( _pd0 + _index);\
		_ar = (_br) * (_pf[2 * _i0 - 2]) - (_bi) * (_pf[2 * _i0 - 1]); \
		_ai = (_br) * (_pf[2 * _i0 - 1]) + (_bi) * (_pf[2 * _i0 - 2]); \
		_br= *( _pd1 + _index);\
		_bi= *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_cr = (_br) * (_pf[2 * _i1 - 2]) - (_bi) * (_pf[2 * _i1 - 1]); \
		_ci = (_br) * (_pf[2 * _i1 - 1]) + (_bi) * (_pf[2 * _i1 - 2]); \
		_t0##r = _ar + _cr; \
		_t1##r = _ar - _cr; \
		_t0##i = _ai + _ci; \
		_t1##i = _ai - _ci; \
	}

# define cplx_load_mulmuladdsub_pp_no_fetch(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= *( _pd0 + _index);\
		_bi= *( _pd0 + _index + 1);\
		_ar = (_br) * (_f0##r) - (_bi) * (_f0##i); \
		_ai = (_br) * (_f0##i) + (_bi) * (_f0##r); \
		_br= *( _pd1 + _index);\
		_bi= *( _pd1 + _index + 1);\
		_cr = (_br) * (_f1##r) - (_bi) * (_f1##i); \
		_ci = (_br) * (_f1##i) + (_bi) * (_f1##r); \
		_t0##r = _ar + _cr; \
		_t1##r = _ar - _cr; \
		_t0##i = _ai + _ci; \
		_t1##i = _ai - _ci; \
	}

# define cplx_preload_mulmuladdsub(_t0,_t1,_f0,_f1,_index,_pd0,_pd1,_tp0,_tp1)\
	{                                                               \
		y_limb_t _ar,_ai,_br,_bi;                               \
		_ar = (_tp0##r) * (_f0##r) - (_tp0##i) * (_f0##i);      \
		_br = (_tp1##r) * (_f1##r) - (_tp1##i) * (_f1##i);      \
		_ai = (_tp0##r) * (_f0##i) + (_tp0##i) * (_f0##r);      \
		_bi = (_tp1##r) * (_f1##i) + (_tp1##i) * (_f1##r);      \
		_t0##r = _ar + _br;                                     \
		_t1##r = _ar - _br;                                     \
		_t0##i = _ai + _bi;                                     \
		_t1##i = _ai - _bi;                                     \
                _tp0##r = *( _pd0 + _index) ;                           \
                _tp0##i = *( _pd0 + _index + 1);                        \
                _tp1##r = *( _pd1 + _index );                           \
                _tp1##i = *( _pd1 + _index + 1);                        \
	}

#define cplx_load_mul_p(_t,_pd,_f,_array,_k)\
{\
	y_limb_t _ar,_ai;\
	_pd = _array + addr((_k)<<1);\
	_ar = *( _pd );\
	_ai = *( _pd + 1);\
         prefetch_p(_pd);\
	_t##r = _ar * _f##r - _ai * _f##i;\
	_t##i = _ai * _f##r + _ar * _f##i;\
}

#define cplx_load_mul_pp(_t,_f,_index,_pd)\
{\
	y_limb_t _ar,_ai;\
	_ar = *( _pd + _index );\
	_ai = *( _pd + _index + 1);\
        prefetch_p( _pd + _index);\
	_t##r = _ar * _f##r - _ai * _f##i;\
	_t##i = _ai * _f##r + _ar * _f##i;\
}

#define cplx_preloaded_mul(_t,_f,_index,_pd,_tp)  \
	_t##r = _tp##r * _f##r - _tp##i * _f##i;  \
	_t##i = _tp##r * _f##i + _tp##i * _f##r;  \
        _tp##r = *( _pd + _index) ;               \
        _tp##i = *( _pd + _index + 1);            \



#define cplx_load_mulmul_p(_t0,_t1,_pd0,_pd1,_f0,_f1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_b,_cr,_ci; \
		_pd0 = _array + addr((_k0)<<1);\
		_b = *( _pd0 );\
		_ai= *( _pd0 + 1);\
                prefetch_p(_pd0);\
		_ar = (_b) * (_f0##r) - (_ai) * (_f0##i); \
		_pd1 = _array + addr((_k1)<<1);\
		_ai = (_b) * (_f0##i) + (_ai) * (_f0##r); \
		_b = *( _pd1 );\
		_ci = *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_cr = (_b) * (_f1##r) - (_ci) * (_f1##i); \
		_ci = (_b) * (_f1##i) + (_ci) * (_f1##r); \
		_t0##r = _ar; \
		_t0##i = _ai; \
		_t1##r = _cr; \
		_t1##i = _ci; \
	}

#define cplx_load_mulmul_pp(_t0,_t1,_f0,_f1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_b,_cr,_ci; \
		_b = *( _pd0 + _index );\
		_ai= *( _pd0 + _index + 1);\
                prefetch_p( _pd0 + _index);\
		_ar = (_b) * (_f0##r) - (_ai) * (_f0##i); \
		_ai = (_b) * (_f0##i) + (_ai) * (_f0##r); \
		_b = *( _pd1 + _index );\
		_ci = *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_cr = (_b) * (_f1##r) - (_ci) * (_f1##i); \
		_ci = (_b) * (_f1##i) + (_ci) * (_f1##r); \
		_t0##r = _ar; \
		_t0##i = _ai; \
		_t1##r = _cr; \
		_t1##i = _ci; \
	}

#define cplx_preloaded_mulmul(_t0,_t1,_f0,_f1,_index,_pd0,_pd1,_tp0,_tp1)  \
	_t0##r = _tp0##r * _f0##r - _tp0##i * _f0##i;     \
	_t0##i = _tp0##r * _f0##i + _tp0##i * _f0##r;     \
	_t1##r = _tp1##r * _f1##r - _tp1##i * _f1##i;     \
	_t1##i = _tp1##r * _f1##i + _tp1##i * _f1##r;     \
        _tp0##r = *( _pd0 + _index) ;                     \
        _tp0##i = *( _pd0 + _index + 1);                  \
        _tp1##r = *( _pd1 + _index );                     \
        _tp1##i = *( _pd1 + _index + 1);


/**********************************************************************
   cplx_divmul macro:
		_t0= _f0 / _f1;
		_t1= _f0 * _f1;
   where all the variables  are pseudo-complex. _f0 and _f1 are module
   1 pseudo-complex (like trig. factros)
**********************************************************************/
#define cplx_divmul(_t0,_t1,_f0,_f1)\
	{\
		y_limb_t _w0,_w1,_w2,_w3;\
		_w0 = _f0##r * _f1##r;\
		_w1 = _f0##i * _f1##i;\
		_w2 = _f0##i * _f1##r;\
		_w3 = _f0##r * _f1##i;\
		_t0##r = _w0 + _w1;\
		_t1##r = _w0 - _w1;\
		_t0##i = _w2 - _w3;\
		_t1##i = _w2 + _w3;\
	}


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/
# define cplx_addsub(_t0,_t1) \
	{ \
		y_limb_t _a=_t0##r; \
		_t0##r += (_t1##r); \
		_t1##r = _a - (_t1##r); \
		_a = (_t0##i); \
		_t0##i += (_t1##i); \
		_t1##i = _a - (_t1##i); \
		}


/**********************************************************************
   cplx_load_addsub macro:
		_t0 = _array( _k0) + _array (_k1);
		_t1 = _array( _k0) - _array (_k1);
 where all the variables are pseudo-complex .
**********************************************************************/

#define cplx_load_addsub(_t0,_t1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = _array[addr(_k1<<1)];\
		_ai = _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_t0##r = _array[addr(_k0<<1)] + _ar; \
		_t1##r = _array[addr(_k0<<1)] - _ar; \
		_t0##i = _array[addr(_k0<<1)+1] + _ai; \
		_t1##i = _array[addr(_k0<<1)+1] - _ai; \
                prefetch_p(&_array[addr(_k0<<1)]);\
	}

#define cplx_load_addsub_p(_t0,_t1,_pd0,_pd1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_pd1 = _array + addr((_k1)<<1);\
		_ar = *( _pd1 );\
		_ai = *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_pd0 = _array + addr((_k0)<<1);\
                prefetch_p(_pd0);\
		_t0##r = *(_pd0 ) + _ar; \
		_t1##r = *(_pd0 ) - _ar; \
		_t0##i = *(_pd0 + 1) + _ai; \
		_t1##i = *(_pd0 + 1) - _ai; \
	}

#define cplx_load_addsub_pp(_t0,_t1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = *( _pd1 + _index);\
		_ai = *( _pd1 + _index + 1);\
                prefetch_p( _pd1 + _index);\
		_t0##r = *( _pd0 + _index) + _ar; \
		_t1##r = *( _pd0 + _index) - _ar; \
		_t0##i = *( _pd0 + _index + 1) + _ai; \
		_t1##i = *( _pd0 + _index + 1) - _ai; \
                prefetch_p( _pd0 + _index);\
	}

#define cplx_load_addsub_pp_no_fetch(_t0,_t1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = *( _pd1 + _index);\
		_ai = *( _pd1 + _index + 1);\
		_t0##r = *( _pd0 + _index) + _ar; \
		_t1##r = *( _pd0 + _index) - _ar; \
		_t0##i = *( _pd0 + _index + 1) + _ai; \
		_t1##i = *( _pd0 + _index + 1) - _ai; \
	}

#define cplx_preload_addsub(_t0,_t1,_index,_pd0,_pd1,_tp0,_tp1)    \
        _t0##r = _tp0##r + _tp1##r;                                \
        _t1##r = _tp0##r - _tp1##r;                                \
        _t0##i = _tp0##i + _tp1##i;                                \
        _t1##i = _tp0##i - _tp1##i;                                \
        _tp0##r = *( _pd0 + _index) ;                              \
        _tp0##i = *( _pd0 + _index + 1);                           \
        _tp1##r = *( _pd1 + _index );                              \
        _tp1##i = *( _pd1 + _index + 1);


#if defined(SUM_CHECK) && defined(YNORM)
/* if defined sumcheck and is the norm-carry step then t0 is used
   to sum the group */

# define cplx_addsub_store(_k0, _k1, _t0, _t1)   \
	_d[addr(_k1 << 1)] = _t0##r - _t1##r;    \
	_d[addr(_k1 << 1) + 1] = _t0##i - _t1##i;\
        _t0##r += _t1##r;                        \
        _t0##i += _t1##i;                        \
	_d[addr(_k0 << 1)] = _t0##r;             \
	_d[addr(_k0 << 1) + 1] = _t0##i;


# define cplx_addsub_store_p(_index, _pd0, _pd1, _t0, _t1)\
	*( _pd1 + _index) = _t0##r - _t1##r;     \
	*( _pd1 + _index + 1) = _t0##i - _t1##i; \
        _t0##r += _t1##r;                        \
        _t0##i += _t1##i;                        \
	*( _pd0 + _index) = _t0##r;              \
	*( _pd0 + _index + 1) = _t0##i;

#else

# define cplx_addsub_store(_k0, _k1, _t0, _t1)   \
	_d[addr(_k0 << 1)] = _t0##r + _t1##r;    \
	_d[addr(_k0 << 1) + 1] = _t0##i + _t1##i;\
	_d[addr(_k1 << 1)] = _t0##r - _t1##r;    \
	_d[addr(_k1 << 1) + 1] = _t0##i - _t1##i;


# define cplx_addsub_store_p(_index,_pd0,_pd1,_t0,_t1)\
	*( _pd0 + _index) = _t0##r + _t1##r;     \
	*( _pd0 + _index + 1) = _t0##i + _t1##i; \
	*( _pd1 + _index) = _t0##r - _t1##r;     \
	*( _pd1 + _index + 1) = _t0##i - _t1##i;

#endif /* SUM_CHECK */
/**********************************************************************
   cplx_mul macro:
		_t= _f0*f1;
   where all the variables are pseudo complex. 
**********************************************************************/
#define cplx_mul(_t,_f0,_f1) \
	{ \
		y_limb_t _r; \
		_r= ((_f0##r) * (_f1##r)) - ((_f0##i) * (_f1##i)); \
		_t##i= ((_f0##r) * (_f1##i)) + ((_f0##i) * (_f1##r)); \
		_t##r = _r; \
	}

/**********************************************************************
   cplx_squar macro:
		_t= _f0 * _f0;
   where all the variables are pseudo-complex. 
***********************************************************************/
#define cplx_squar(_t,_f0) \
	{ \
		y_limb_t _r; \
		_r= (_f0##r +_f0##i ) * (_f0##r - _f0##i);\
		_t##i= 2.0 * (_f0##r) * (_f0##i) ;\
		_t##r = _r; \
	}

/***********************************************************************
   cplx_add macro:
		_t= _s0 + _s1;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_add(_t,_s0,_s1) \
	_t##r= (_s0##r) + (_s1##r); \
	_t##i= (_s0##i) + (_s1##i);

/***********************************************************************
   cplx_sub macro:
		_t= _m - _s;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_sub(_t,_m,_s) \
	_t##r= (_m##r) - (_s##r); \
	_t##i= (_m##i) - (_s##i);


/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
/*
# define cplx_mul_1_4_F_addsub(_t0,_t1) \
{\
   y_limb_t _a0=_t1##r,_a1=_t1##i;\
   _t1##r = _t0##r - _a1;\
   _t0##r += _a1;\
   _t1##i = _t0##i + _a0;\
   _t0##i -= _a0;\
}
*/

# define cplx_mul_1_4_F_addsub(_t0,_t1) \
{ \
	y_limb_t _ar= _t0##r, _ai = _t0##i ,_a = _t1##r; \
	_t0##r += (_t1##i); \
	_t0##i -= _a; \
	_t1##r = _ar - (_t1##i); \
	_t1##i = _ai + _a; \
}


/***********************************************************************
   cplx_mul_1_4_F_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
	_d[addr( _k0<<1 )]= _t0##r + _t1##i; \
	_d[addr( _k1<<1 )]= _t0##r - _t1##i; \
	_d[addr( _k0<<1 )+ 1]= _t0##i - _t1##r; \
	_d[addr( _k1<<1 )+ 1]= _t0##i + _t1##r;

#define cplx_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	*(_pd0 + _index) = _t0##r + _t1##i; \
	*(_pd1 + _index) = _t0##r - _t1##i; \
	*( _pd0 + _index + 1) = _t0##i - _t1##r; \
	*( _pd1 + _index + 1) = _t0##i + _t1##r;



/***********************************************************************
   cplx_mul_1_4_B macro:
		_t = _t1 * G^(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
/*
# define cplx_mul_1_4_B_addsub(_t0,_t1) \
{\
   y_limb_t _a0=_t1##r,_a1=_t1##i;\
   _t1##r = _t0##r + _a1;\
   _t0##r -= _a1;\
   _t1##i = _t0##i - _a0;\
   _t0##i += _a0;\
}
*/

#define cplx_mul_1_4_B_addsub(_t0,_t1) \
{ \
	y_limb_t _ar=_t0##r, _ai=_t0##i ,_a= _t1##r; \
	_t0##r -= (_t1##i); \
	_t0##i += _a; \
	_t1##r = _ar + (_t1##i); \
	_t1##i = _ai - _a; \
}


/***********************************************************************
   cplx_mul_1_4_B_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
	_d[addr( _k0<<1 )]= _t0##r - _t1##i; \
	_d[addr( _k1<<1 )]= _t0##r + _t1##i; \
	_d[addr( _k0<<1 )+ 1]= _t0##i + _t1##r; \
	_d[addr( _k1<<1 )+ 1]= _t0##i - _t1##r;

#define cplx_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	*( _pd0 + _index) = _t0##r - _t1##i; \
	*( _pd1 + _index) = _t0##r + _t1##i; \
	*( _pd0 + _index + 1) = _t0##i + _t1##r; \
	*( _pd1 + _index + 1) = _t0##i - _t1##r;


/***********************************************************************
   cplx_mul_1_8_F macro:
		_t = _t1 * G^-(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_8_F_addsub(_tt0,_tt1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##r + _tt1##i )*F_1_8r; \
		_ai = (_tt1##i - _tt1##r )*F_1_8r; \
		_tt1##r = _tt0##r - _ar;\
		_tt0##r += _ar;\
		_tt1##i = _tt0##i - _ai;\
		_tt0##i += _ai;\
	}

/***********************************************************************
   cplx_mul_1_8_F_addsub_store macro:
		_t = _t1 * G^-(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r + _t1##i )*F_1_8r; \
		_ai = (_t1##i - _t1##r )*F_1_8r; \
		_d[addr( _k0<<1 )]= _t0##r + _ar; \
		_d[addr( _k1<<1 )]= _t0##r - _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai; \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai; \
	}

#define cplx_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r + _t1##i )*F_1_8r; \
		_ai = (_t1##i - _t1##r )*F_1_8r; \
		*( _pd0 + _index)= _t0##r + _ar; \
		*( _pd1 + _index)= _t0##r - _ar; \
		*( _pd0 + _index + 1)= _t0##i + _ai; \
		*( _pd1 + _index + 1)= _t0##i - _ai; \
	}

/***********************************************************************
   cplx_mul_1_8_B macro:
		_t = _t1 * G^(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_8_B_addsub(_tt0,_tt1) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##r - _tt1##i )*F_1_8r; \
		_ai = (_tt1##i + _tt1##r )*F_1_8r; \
		_tt1##r = _tt0##r - _ar;\
		_tt0##r += _ar;\
		_tt1##i = _tt0##i - _ai;\
		_tt0##i += _ai;\
	}

/***********************************************************************
   cplx_mul_1_8_B_addsub_store macro:
		_t = _t1 * G^(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r - _t1##i )*F_1_8r; \
		_ai = (_t1##i + _t1##r )*F_1_8r; \
		_d[addr( _k0<<1 )]= _t0##r + _ar; \
		_d[addr( _k1<<1 )]= _t0##r - _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai; \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai; \
	}

#define cplx_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r - _t1##i )*F_1_8r; \
		_ai = (_t1##i + _t1##r )*F_1_8r; \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
	}


/***********************************************************************
   cplx_mul_3_8_F macro:
		_t = _t1 * G^-(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_3_8_F_addsub(_tt0,_tt1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##r - _tt1##i )*F_1_8i; \
		_ai = (_tt1##i + _tt1##r )*F_1_8i; \
		_tt1##r = _tt0##r - _ar;\
		_tt0##r += _ar;\
		_tt1##i = _tt0##i - _ai;\
		_tt0##i += _ai;\
	}

/***********************************************************************
   cplx_mul_3_8_F_addsub_store macro:
		_t = _t1 * G^-(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_F_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r - _t1##i )*F_1_8i; \
		_ai = (_t1##i + _t1##r )*F_1_8i; \
		_d[addr( _k0<<1 )]= _t0##r + _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai; \
		_d[addr( _k1<<1 )]= _t0##r - _ar; \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai; \
	}

#define cplx_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##i - _t1##r )*F_1_8r; \
		_ai = (_t1##i + _t1##r )*F_1_8r; \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd0 + _index + 1) = _t0##i - _ai; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd1 + _index + 1) = _t0##i + _ai; \
	}

/***********************************************************************
   cplx_mul_3_8_B macro:
		_t = _t1 * G^(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_B_addsub(_tt0,_tt1) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1##r + _tt1##i )*F_1_8i; \
		_ai = (_tt1##i - _tt1##r )*F_1_8i; \
		_tt1##r = _tt0##r - _ar;\
		_tt0##r += _ar;\
		_tt1##i = _tt0##i - _ai;\
		_tt0##i += _ai;\
	}

/***********************************************************************
   cplx_mul_3_8_B_addsub_store macro:
		_t = _t1 * G^(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_B_addsub_store(_d,_k0,_k1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r + _t1##i )*F_1_8i; \
		_ai = (_t1##i - _t1##r )*F_1_8i; \
		_d[addr( _k0<<1 )]= _t0##r + _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0##i + _ai; \
		_d[addr( _k1<<1 )]= _t0##r - _ar; \
		_d[addr( _k1<<1 )+ 1]= _t0##i - _ai; \
	}

#define cplx_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0,_t1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r + _t1##i )*F_1_8r; \
		_ai = (_t1##r - _t1##i )*F_1_8r; \
		*( _pd0 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index) = _t0##r + _ar; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
	}



/***********************************************************************
 This is the dyadic mul for nested complex representation 
 xk =(xk + xmk* ) * (yk + ymk* ) + 2*((xk * yk) - (xmk* * ymk* ) -
     G^(-k) * (xk - xmk* ) * (yk - ymk* );
 xkm =(xkm + xk* ) * (ymk + yk* ) + 2*((xmk * xk*) - (ymk * yk* ) -
     G^(k) * (xmk - xk* ) * (ymk - yk* );
 where all are pseudo-complex vars. 
***********************************************************************/
/* 17 Fadds 18 FMULS */
#define conv_nested(xk,xmk,yk,ymk,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;\
   _y0i=(xk##i + xmk##i)*0.25;\
   _y1r=(yk##r - ymk##r);\
   _y1i=(yk##i + ymk##i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);\
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);\
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);\
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);\
   _y3r= ((1.0 +(gkr)) * _y2r - (gki) * _y2i);\
   _y3i= ((1.0 +(gkr)) * _y2i + (gki) * _y2r);\
   xk##r = _y0r - _y3r;\
   xk##i = _y0i - _y3i;\
   xmk##r = _y1r - _y3r;\
   xmk##i = _y1i + _y3i;\
}

/*********************************************************************
This is the dyadic square for nested-compls representation
 xk =(xk + xmk* ) * (xk + xmk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(-k) * (xk - xmk* ) * (xk - xmk* );
 xkm =(xkm + xk* ) * (xmk + xk* ) + 2*((xk * xk) - (xmk* * xmk* ) -
     G^(k) * (xmk - xk* ) * (xmk - xk* );
 where all are pseudo-complex vars. 
**********************************************************************/
/* 18 fadds  12  Fmuls*/
#define square_nested(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= (((gkr)+ 1.0) * _r2 - (gki) * _r1);\
   _r1= (((gkr)+ 1.0) * _r1 + (gki) * _r2);\
   xk##r -= _r0;\
   xmk##r -= _r0;\
   xk##i -= _r1;\
   xmk##i += _r1;\
}

#define square_nested_eq(xk,gkr) \
{\
   y_limb_t _y0r=xk##r,_y0i=xk##i;\
   xk##r =(_y0r * _y0r + (gkr) * _y0i * _y0i);\
   xk##i = 2.0 * _y0r * _y0i;\
}

/* These are version to exploid better the symmetries */
#define conv_nested_1_4(xk,xmk,yk,ymk,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;\
   _y0i=(xk##i + xmk##i)*0.25;\
   _y1r=(yk##r - ymk##r);\
   _y1i=(yk##i + ymk##i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);\
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);\
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);\
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);\
   _y3r= ((gki +1.0) * _y2r + gkr * _y2i);\
   _y3i= ((gki +1.0) * _y2i - gkr * _y2r);\
   xk##r = _y0r - _y3r;\
   xk##i = _y0i - _y3i;\
   xmk##r = _y1r - _y3r;\
   xmk##i = _y1i + _y3i;\
}

#define conv_nested_1_2(xk,xmk,yk,ymk,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;\
   _y0i=(xk##i + xmk##i)*0.25;\
   _y1r=(yk##r - ymk##r);\
   _y1i=(yk##i + ymk##i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);\
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);\
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);\
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);\
   _y3r= ((gkr - 1.0) * _y2r - gki * _y2i);\
   _y3i= ((gkr - 1.0) * _y2i + gki * _y2r);\
   xk##r = _y0r + _y3r;\
   xk##i = _y0i + _y3i;\
   xmk##r = _y1r + _y3r;\
   xmk##i = _y1i - _y3i;\
}

#define conv_nested_3_4(xk,xmk,yk,ymk,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xk##r - xmk##r)*0.25;\
   _y0i=(xk##i + xmk##i)*0.25;\
   _y1r=(yk##r - ymk##r);\
   _y1i=(yk##i + ymk##i);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xk##r * yk##r) - (xk##i * yk##i);\
   _y0i=(xk##r * yk##i) + (xk##i * yk##r);\
   _y1r=(xmk##r * ymk##r) - (xmk##i * ymk##i);\
   _y1i=(xmk##r * ymk##i) + (xmk##i * ymk##r);\
   _y3r= ((gki - 1.0) * _y2r + gkr * _y2i);\
   _y3i= ((gki - 1.0) * _y2i - gkr * _y2r);\
   xk##r = _y0r + _y3r;\
   xk##i = _y0i + _y3i;\
   xmk##r = _y1r + _y3r;\
   xmk##i = _y1i - _y3i;\
}

#define square_nested_1_4(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gki + 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki + 1.0) * _r1 - gkr * _r2);\
   xk##r -= _r0;\
   xmk##r -= _r0;\
   xk##i -= _r1;\
   xmk##i += _r1;\
}

#define square_nested_1_2(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gkr - 1.0) * _r2 - gki * _r1);\
   _r1= ((gkr - 1.0) * _r1 + gki * _r2);\
   xk##r += _r0;\
   xmk##r += _r0;\
   xk##i += _r1;\
   xmk##i -= _r1;\
}

#define square_nested_3_4(xk,xmk,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xk##r - xmk##r)*0.5;\
   _r1=(xk##i + xmk##i)*0.5;\
   _r2=(xk##r + xk##i) * (xk##r - xk##i);\
   xk##i *= xk##r;\
   xk##r = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmk##r + xmk##i) * (xmk##r - xmk##i);\
   xmk##i *= xmk##r;\
   xmk##r = _r3;\
   _r1 +=  _r1;\
   xk##i += xk##i;\
   xmk##i += xmk##i;\
   _r0= ((gki - 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki - 1.0) * _r1 - gkr * _r2);\
   xk##r += _r0;\
   xmk##r += _r0;\
   xk##i += _r1;\
   xmk##i -= _r1;\
}


/*$Id$*/










