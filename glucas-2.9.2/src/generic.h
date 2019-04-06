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
#ifndef Y_CACHE_LINE
# define Y_CACHE_LINE 4
#endif
/* Prefetch is not defined for generic code */
/*  # define prefetch_data(_array,_k)*/ /* */

/* # define prefetch_p(_pd) */ /* */


/* to read some trigonomeric factors */
#define cplx_data_trig_load_2(_t0r,_t0i,_t1r,_ti1,_p) \
   _t0r = *( _p );\
   _t0i = *( _p + 1);\
   _t1r = *( _p + 2);\
   _t1i = *( _p + 3);

#define cplx_data_trig_load_3(_t0r,_t0i,_t1r,_t1i,_t2r,_t2i,_p) \
   _t0r = *( _p );\
   _t0i = *( _p + 1);\
   _t1r = *( _p + 2);\
   _t1i = *( _p + 3);\
   _t2r = *( _p + 4);\
   _t2i = *( _p + 5);


/* To load from real padded array of data as a complex array */
#define cplx_data_to_local(_tr,_ti,_array,_k) \
	_tr = _array[addr(( _k )<<1)]; \
	_ti = _array[addr(( _k )<<1)+1];


#define cplx_data_to_local_p(_tr,_ti,_pd,_array,_k) \
	_pd = _array + addr(( _k)<<1);\
	_tr = *(_pd ); \
	_ti = *( _pd + 1);\
        prefetch_p(_pd);


#define cplx_data_to_local_pp(_tr,_ti,_index,_pd) \
	_tr = *( _pd + _index); \
	_ti = *( _pd + _index + 1);


/* To store on read padded array from local data as complex */

#define cplx_local_to_data(_array,_k,_tr,_ti) \
	_array[addr(( _k )<<1)]=_tr ; \
	_array[addr(( _k )<<1)+1]=_ti;


#define cplx_local_to_data_p(_pd,_tr,_ti) \
	*( _pd ) = _tr ; \
	*( _pd + 1) = _ti;


/* To load from real non-padded array of data as a complex array */
#define cplx_mem_to_local(_tr,_ti,_array,_k) \
	_tr =_array[( _k )<<1]; \
	_ti =_array[(( _k )<<1)+1];

#define cplx_mem_to_local_p(_tr,_ti,_pd) \
	_tr =  *(_pd);\
	_ti =  *(_pd + 1); \
         prefetch_data (_pd + 2);


/* To strore to a non-padded array from local */
#define cplx_local_to_mem(_array,_k,_tr,_ti) \
	_array[( _k )<<1]=_tr; \
	_array[(( _k )<<1)+1]=_ti;


/**************** PSEUDO_COMPLEX GENERAL AUXILIAR MACROS **************/
/**********************************************************************
   cplx_muladdsub macro :
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex and _t is temporal aux var.
***********************************************************************/

# define cplx_muladdsub(_t0r,_t0i,_t1r,_t1i,_fr,_fi) \
	{ \
		y_limb_t _ar,_ai,_a=_t0r; \
		_ar = (_t1r ) * (_fr ) - (_t1i ) * (_fi ); \
		_ai = (_t1r ) * (_fi ) + (_t1i ) * (_fr ); \
		_t0r += _ar;\
		_t1r = _a - _ar;\
		_a = _t0i;\
		_t0i += _ai;\
		_t1i = _a - _ai;\
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
#define cplx_load_muladdsub(_t0r,_t0i,_t1r,_t1i,_fr,_fi,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_ar = (_br) * (_fr) - (_bi) * (_fi); \
		_ai = (_br) * (_fi) + (_bi) * (_fr); \
		_t0r = _array[addr(_k0<<1)] + _ar; \
		_t1r = _array[addr(_k0<<1)] - _ar; \
		_t0i = _array[addr(_k0<<1)+1] + _ai; \
		_t1i = _array[addr(_k0<<1)+1] - _ai; \
                prefetch_p(&_array[addr(_k0<<1)]);\
	}


#define cplx_load_muladdsub_p(_t0r,_t0i,_t1r,_t1i,_pd0,_pd1,_fr,_fi,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_pd1= _array + addr((_k1)<<1);\
		_br= *( _pd1 );\
		_bi= *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_ar = (_br) * (_fr) - (_bi) * (_fi); \
		_pd0 = _array + addr((_k0)<<1);\
		_ai = (_br) * (_fi) + (_bi) * (_fr); \
		_t0r = *( _pd0 ) + _ar; \
		_t1r = *( _pd0 ) - _ar; \
		_t0i = *( _pd0 + 1) + _ai; \
		_t1i = *( _pd0 + 1) - _ai; \
                prefetch_p(_pd0);\
	}

#define cplx_load_muladdsub_pp(_t0r,_t0i,_t1r,_t1i,_fr,_fi,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
		_ar = (_br) * (_fr) - (_bi) * (_fi); \
		_ai = (_br) * (_fi) + (_bi) * (_fr); \
		_t0r = *( _pd0 + _index) + _ar; \
		_t1r = *( _pd0 + _index) - _ar; \
		_t0i = *( _pd0 + _index + 1) + _ai; \
		_t1i = *( _pd0 + _index + 1) - _ai; \
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
#define cplx_muladdsub_store(_array,_k0,_k1,_t0r,_t0i,_t1r,_t1i,_fr,_fi) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r) * (_fr) - (_t1i) * (_fi); \
		_ai = (_t1r) * (_fi) + (_t1i) * (_fr); \
		_array[addr(_k0<<1)]= _t0r + _ar; \
                _array[addr(_k1<<1)]= _t0r - _ar; \
		_array[addr(_k0<<1)+1]= _t0i + _ai; \
                _array[addr(_k1<<1)+1]= _t0i - _ai; \
        }


#define cplx_muladdsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i,_fr,_fi) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r) * (_fr) - (_t1i) * (_fi); \
		_ai = (_t1r) * (_fi) + (_t1i) * (_fr); \
		*( _pd0 + _index) = _t0r + _ar; \
		*( _pd1 + _index) = _t0r - _ar; \
		*( _pd0 + _index + 1) = _t0i + _ai; \
		*( _pd1 + _index + 1) = _t0i - _ai; \
	}




/**********************************************************************
   cplx_mulmuladdsub macro:
		_a= _t0 * _f0;
		_b= _t1 * _f1;
		_t0= _a + _b;
		_t1= _a - _b;
   where all the variables  are pseudo-complex. _a and _b are temp. aux.
**********************************************************************/
# define cplx_mulmuladdsub(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i) \
	{ \
	        y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1r) * (_f1r) - (_t1i) * (_f1i); \
		_br = (_t0r) * (_f0r) - (_t0i) * (_f0i); \
		_ai = (_t1r) * (_f1i) + (_t1i) * (_f1r); \
		_bi = (_t0r) * (_f0i) + (_t0i) * (_f0r); \
                _t0r = _br + _ar;\
                _t1r = _br - _ar;\
                _t0i = _bi + _ai;\
		_t1i = _bi - _ai;\
	}


#define cplx_mulmul(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i) \
	{ \
		y_limb_t _ar,_br; \
		_ar = (_t1r) * (_f1r) - (_t1i) * (_f1i); \
		_br = (_t0r) * (_f0r) - (_t0i) * (_f0i); \
		_t1i = (_t1r) * (_f1i) + (_t1i) * (_f1r); \
		_t0i = (_t0r) * (_f0i) + (_t0i) * (_f0r); \
		_t1r = _ar;\
		_t0r = _br;\
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
#define cplx_load_mulmuladdsub(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_br= _array[addr(_k0<<1)];\
		_bi= _array[addr(_k0<<1)+1];\
                prefetch_p(&_array[addr(_k0<<1)]);\
		_ar = (_br) * (_f0r) - (_bi) * (_f0i); \
		_ai = (_br) * (_f0i) + (_bi) * (_f0r); \
		_br= _array[addr(_k1<<1)];\
		_bi= _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_cr = (_br) * (_f1r) - (_bi) * (_f1i); \
		_ci = (_br) * (_f1i) + (_bi) * (_f1r); \
		_t0r = _ar + _cr; \
		_t1r = _ar - _cr; \
		_t0i = _ai + _ci; \
		_t1i = _ai - _ci; \
	}


# define cplx_load_mulmuladdsub_p(_t0r,_t0i,_t1r,_t1i,_pd0,_pd1,_f0r,_f0i,_f1r,_f1i,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci; \
		_pd0 = _array + addr((_k0)<<1);\
		_br= *( _pd0 );\
		_bi= *( _pd0 + 1);\
                prefetch_p(_pd0);\
		_ar = (_br) * (_f0r) - (_bi) * (_f0i); \
		_pd1 = _array + addr((_k1)<<1);\
		_ai = (_br) * (_f0i) + (_bi) * (_f0r); \
		_br= *( _pd1 );\
		_bi= *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_cr = (_br) * (_f1r) - (_bi) * (_f1i); \
		_ci = (_br) * (_f1i) + (_bi) * (_f1r); \
		_t0r = _ar + _cr; \
		_t1r = _ar - _cr; \
		_t0i = _ai + _ci; \
		_t1i = _ai - _ci; \
	}

# define cplx_load_mulmuladdsub_pp(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci;      \
		_br= *( _pd0 + _index);                \
		_bi= *( _pd0 + _index + 1);            \
		_ar = (_br) * (_f0r) - (_bi) * (_f0i); \
		_ai = (_br) * (_f0i) + (_bi) * (_f0r); \
                prefetch_p( _pd0 + _index);            \
		_br= *( _pd1 + _index);                \
		_bi= *( _pd1 + _index + 1);            \
		_cr = (_br) * (_f1r) - (_bi) * (_f1i); \
		_ci = (_br) * (_f1i) + (_bi) * (_f1r); \
                prefetch_p( _pd1 + _index);            \
		_t0r = _ar + _cr;                      \
		_t1r = _ar - _cr;                      \
		_t0i = _ai + _ci;                      \
		_t1i = _ai - _ci;                      \
	}

#define cplx_load_mul_p(_tr,_ti,_pd,_fr,_fi,_array,_k)\
{\
	y_limb_t _ar,_ai;\
	_pd = _array + addr((_k)<<1);\
	_ar = *( _pd );\
	_ai = *( _pd + 1);\
         prefetch_p(_pd);\
	_tr = _ar * _fr - _ai * _fi;\
	_ti = _ai * _fr + _ar * _fi;\
}

#define cplx_load_mul_pp(_tr,_ti,_fr,_fi,_index,_pd)\
{\
	y_limb_t _ar,_ai;\
	_ar = *( _pd + _index );\
	_ai = *( _pd + _index + 1);\
	_tr = _ar * _fr - _ai * _fi;\
	_ti = _ai * _fr + _ar * _fi;\
}

#define cplx_load_mulmul_p(_t0r,_t0i,_t1r,_t1i,_pd0,_pd1,_f0r,_f0i,_f1r,_f1i,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai,_b,_cr,_ci; \
		_pd0 = _array + addr((_k0)<<1);\
		_b = *( _pd0 );\
		_ai= *( _pd0 + 1);\
                prefetch_p(_pd0);\
		_ar = (_b) * (_f0r) - (_ai) * (_f0i); \
		_pd1 = _array + addr((_k1)<<1);\
		_ai = (_b) * (_f0i) + (_ai) * (_f0r); \
		_b = *( _pd1 );\
		_ci = *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_cr = (_b) * (_f1r) - (_ci) * (_f1i); \
		_ci = (_b) * (_f1i) + (_ci) * (_f1r); \
		_t0r = _ar; \
		_t0i = _ai; \
		_t1r = _cr; \
		_t1i = _ci; \
	}

#define cplx_load_mulmul_pp(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_b,_cr,_ci; \
		_b = *( _pd0 + _index );\
		_ai= *( _pd0 + _index + 1);\
		_ar = (_b) * (_f0r) - (_ai) * (_f0i); \
		_ai = (_b) * (_f0i) + (_ai) * (_f0r); \
		_b = *( _pd1 + _index );\
		_ci = *( _pd1 + _index + 1);\
		_cr = (_b) * (_f1r) - (_ci) * (_f1i); \
		_ci = (_b) * (_f1i) + (_ci) * (_f1r); \
		_t0r = _ar; \
		_t0i = _ai; \
		_t1r = _cr; \
		_t1i = _ci; \
	}


/**********************************************************************
   cplx_divmul macro:
		_t0= _f0 / _f1;
		_t1= _f0 * _f1;
   where all the variables  are pseudo-complex. _f0 and _f1 are module
   1 pseudo-complex (like trig. factros)
**********************************************************************/
#define cplx_divmul(_t0r,_t0i,_t1r,_t1i,_f0r,_f0i,_f1r,_f1i)\
	{\
		y_limb_t _w0,_w1,_w2,_w3;\
		_w0 = _f0r * _f1r;\
		_w1 = _f0i * _f1i;\
		_w2 = _f0i * _f1r;\
		_w3 = _f0r * _f1i;\
		_t0r = _w0 + _w1;\
		_t1r = _w0 - _w1;\
		_t0i = _w2 - _w3;\
		_t1i = _w2 + _w3;\
	}


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/
# define cplx_addsub(_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _a=_t0r; \
		_t0r += (_t1r); \
		_t1r = _a - (_t1r); \
		_a = (_t0i); \
		_t0i += (_t1i); \
		_t1i = _a - (_t1i); \
		}


/**********************************************************************
   cplx_load_addsub macro:
		_t0 = _array( _k0) + _array (_k1);
		_t1 = _array( _k0) - _array (_k1);
 where all the variables are pseudo-complex .
**********************************************************************/

#define cplx_load_addsub(_t0r,_t0i,_t1r,_t1i,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = _array[addr(_k1<<1)];\
		_ai = _array[addr(_k1<<1)+1];\
                prefetch_p(&_array[addr(_k1<<1)]);\
		_t0r = _array[addr(_k0<<1)] + _ar; \
		_t1r = _array[addr(_k0<<1)] - _ar; \
		_t0i = _array[addr(_k0<<1)+1] + _ai; \
		_t1i = _array[addr(_k0<<1)+1] - _ai; \
                prefetch_p(&_array[addr(_k0<<1)]);\
	}

#define cplx_load_addsub_p(_t0r,_t0i,_t1r,_t1i,_pd0,_pd1,_array,_k0,_k1) \
	{ \
		y_limb_t _ar,_ai; \
		_pd1 = _array + addr((_k1)<<1);\
		_ar = *( _pd1 );\
		_ai = *( _pd1 + 1);\
                prefetch_p(_pd1);\
		_pd0 = _array + addr((_k0)<<1);\
                prefetch_p(_pd0);\
		_t0r = *(_pd0 ) + _ar; \
		_t1r = *(_pd0 ) - _ar; \
		_t0i = *(_pd0 + 1) + _ai; \
		_t1i = *(_pd0 + 1) - _ai; \
	}

#define cplx_load_addsub_pp(_t0r,_t0i,_t1r,_t1i,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = *( _pd1 + _index);\
		_ai = *( _pd1 + _index + 1);\
		_t0r = *( _pd0 + _index) + _ar; \
		_t1r = *( _pd0 + _index) - _ar; \
		_t0i = *( _pd0 + _index + 1) + _ai; \
		_t1i = *( _pd0 + _index + 1) - _ai; \
	}


#define cplx_addsub_store(_k0,_k1,_t0r,_t0i,_t1r,_t1i)\
	_d[addr(_k0<<1)] = _t0r + _t1r;\
	_d[addr(_k0<<1)+1] = _t0i + _t1i;\
	_d[addr(_k1<<1)] = _t0r - _t1r;\
	_d[addr(_k1<<1)+1] = _t0i - _t1i;


#define cplx_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i)\
	*( _pd0 + _index) = _t0r + _t1r;\
	*( _pd0 + _index + 1) = _t0i + _t1i;\
	*( _pd1 + _index) = _t0r - _t1r;\
	*( _pd1 + _index + 1) = _t0i - _t1i;


/**********************************************************************
   cplx_mul macro:
		_t= _f0*f1;
   where all the variables are pseudo complex. 
**********************************************************************/
#define cplx_mul(_tr,_ti,_f0r,_f0i,_f1r,_f1i) \
	{ \
		y_limb_t _r; \
		_r= ((_f0r) * (_f1r)) - ((_f0i) * (_f1i)); \
		_ti= ((_f0r) * (_f1i)) + ((_f0i) * (_f1r)); \
		_tr = _r; \
	}

/**********************************************************************
   cplx_squar macro:
		_t= _f0 * _f0;
   where all the variables are pseudo-complex. 
***********************************************************************/
#define cplx_squar(_tr,_ti,_fr,_fi) \
	{ \
		y_limb_t _r; \
		_r= (_fr +_fi ) * (_fr - _fi);\
		_ti= 2.0 * (_fr) * (_fi) ;\
		_tr = _r; \
	}

/***********************************************************************
   cplx_add macro:
		_t= _s0 + _s1;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_add(_tr,_ti,_s0r,_s0i,_s1r,_s1i) \
	_tr= (_s0r) + (_s1r); \
	_ti= (_s0i) + (_s1i);

/***********************************************************************
   cplx_sub macro:
		_t= _m - _s;
   where all the variables are pseudo-complex. 
************************************************************************/
#define cplx_sub(_tr,_ti,_mr,_mi,_sr,_si) \
	_tr= (_mr) - (_sr); \
	_ti= (_mi) - (_si);


/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_4_F_addsub(_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar= _t0r, _ai = _t0i ,_a = _t1r; \
		_t0r += (_t1i); \
		_t0i -= _a; \
		_t1r = _ar - (_t1i); \
		_t1i = _ai + _a; \
	}

/***********************************************************************
   cplx_mul_1_4_F_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_F_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	_d[addr( _k0<<1 )]= _t0r + _t1i; \
	_d[addr( _k1<<1 )]= _t0r - _t1i; \
	_d[addr( _k0<<1 )+ 1]= _t0i - _t1r; \
	_d[addr( _k1<<1 )+ 1]= _t0i + _t1r;

#define cplx_mul_1_4_F_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	*(_pd0 + _index) = _t0r + _t1i; \
	*(_pd1 + _index) = _t0r - _t1i; \
	*( _pd0 + _index + 1) = _t0i - _t1r; \
	*( _pd1 + _index + 1) = _t0i + _t1r;



/***********************************************************************
   cplx_mul_1_4_B macro:
		_t = _t1 * G^(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub(_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar=_t0r, _ai=_t0i ,_a= _t1r; \
		_t0r -= (_t1i); \
		_t0i += _a; \
		_t1r = _ar + (_t1i); \
		_t1i = _ai - _a; \
	}

/***********************************************************************
   cplx_mul_1_4_B_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	_d[addr( _k0<<1 )]= _t0r - _t1i; \
	_d[addr( _k1<<1 )]= _t0r + _t1i; \
	_d[addr( _k0<<1 )+ 1]= _t0i + _t1r; \
	_d[addr( _k1<<1 )+ 1]= _t0i - _t1r;

#define cplx_mul_1_4_B_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	*( _pd0 + _index) = _t0r - _t1i; \
	*( _pd1 + _index) = _t0r + _t1i; \
	*( _pd0 + _index + 1) = _t0i + _t1r; \
	*( _pd1 + _index + 1) = _t0i - _t1r;


/***********************************************************************
   cplx_mul_1_8_F macro:
		_t = _t1 * G^-(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_8_F_addsub(_tt0r,_tt0i,_tt1r,_tt1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1r + _tt1i )*F_1_8r; \
		_ai = (_tt1i - _tt1r )*F_1_8r; \
		_tt1r = _tt0r - _ar;\
		_tt0r += _ar;\
		_tt1i = _tt0i - _ai;\
		_tt0i += _ai;\
	}

/***********************************************************************
   cplx_mul_1_8_F_addsub_store macro:
		_t = _t1 * G^-(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_F_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r + _t1i )*F_1_8r; \
		_ai = (_t1i - _t1r )*F_1_8r; \
		_d[addr( _k0<<1 )]= _t0r + _ar; \
		_d[addr( _k1<<1 )]= _t0r - _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0i + _ai; \
		_d[addr( _k1<<1 )+ 1]= _t0i - _ai; \
	}

#define cplx_mul_1_8_F_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r + _t1i )*F_1_8r; \
		_ai = (_t1i - _t1r )*F_1_8r; \
		*( _pd0 + _index)= _t0r + _ar; \
		*( _pd1 + _index)= _t0r - _ar; \
		*( _pd0 + _index + 1)= _t0i + _ai; \
		*( _pd1 + _index + 1)= _t0i - _ai; \
	}

/***********************************************************************
   cplx_mul_1_8_B macro:
		_t = _t1 * G^(1/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_8_B_addsub(_tt0r,_tt0i,_tt1r,_tt1i) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1r - _tt1i )*F_1_8r; \
		_ai = (_tt1i + _tt1r )*F_1_8r; \
		_tt1r = _tt0r - _ar;\
		_tt0r += _ar;\
		_tt1i = _tt0i - _ai;\
		_tt0i += _ai;\
	}

/***********************************************************************
   cplx_mul_1_8_B_addsub_store macro:
		_t = _t1 * G^(1/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_8_B_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r - _t1i )*F_1_8r; \
		_ai = (_t1i + _t1r )*F_1_8r; \
		_d[addr( _k0<<1 )]= _t0r + _ar; \
		_d[addr( _k1<<1 )]= _t0r - _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0i + _ai; \
		_d[addr( _k1<<1 )+ 1]= _t0i - _ai; \
	}

#define cplx_mul_1_8_B_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r - _t1i )*F_1_8r; \
		_ai = (_t1i + _t1r )*F_1_8r; \
		*( _pd0 + _index) = _t0r + _ar; \
		*( _pd1 + _index) = _t0r - _ar; \
		*( _pd0 + _index + 1) = _t0i + _ai; \
		*( _pd1 + _index + 1) = _t0i - _ai; \
	}


/***********************************************************************
   cplx_mul_3_8_F macro:
		_t = _t1 * G^-(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_3_8_F_addsub(_tt0r,_tt0i,_tt1r,_tt1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_tt1r - _tt1i )*F_1_8i; \
		_ai = (_tt1i + _tt1r )*F_1_8i; \
		_tt1r = _tt0r - _ar;\
		_tt0r += _ar;\
		_tt1i = _tt0i - _ai;\
		_tt0i += _ai;\
	}

/***********************************************************************
   cplx_mul_3_8_F_addsub_store macro:
		_t = _t1 * G^-(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_F_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r - _t1i )*F_1_8i; \
		_ai = (_t1i + _t1r )*F_1_8i; \
		_d[addr( _k0<<1 )]= _t0r + _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0i + _ai; \
		_d[addr( _k1<<1 )]= _t0r - _ar; \
		_d[addr( _k1<<1 )+ 1]= _t0i - _ai; \
	}

#define cplx_mul_3_8_F_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r - _t1i )*F_1_8i; \
		_ai = (_t1i + _t1r )*F_1_8i; \
		*( _pd0 + _index) = _t0r + _ar; \
		*( _pd0 + _index + 1) = _t0i + _ai; \
		*( _pd1 + _index) = _t0r - _ar; \
		*( _pd1 + _index + 1) = _t0i - _ai; \
	}

/***********************************************************************
   cplx_mul_3_8_B macro:
		_t = _t1 * G^(3/8);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_B_addsub(_tt0r,_tt0i,_tt1r,_tt1i) \
        { \
		y_limb_t _ar,_ai; \
		_ar = (_tt1r + _tt1i )*F_1_8i; \
		_ai = (_tt1i - _tt1r )*F_1_8i; \
		_tt1r = _tt0r - _ar;\
		_tt0r += _ar;\
		_tt1i = _tt0i - _ai;\
		_tt0i += _ai;\
	}

/***********************************************************************
   cplx_mul_3_8_B_addsub_store macro:
		_t = _t1 * G^(3/8);
		_d[k0]= _t0 + _t;
		_t[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_3_8_B_addsub_store(_d,_k0,_k1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r + _t1i )*F_1_8i; \
		_ai = (_t1i - _t1r )*F_1_8i; \
		_d[addr( _k0<<1 )]= _t0r + _ar; \
		_d[addr( _k0<<1 )+ 1]= _t0i + _ai; \
		_d[addr( _k1<<1 )]= _t0r - _ar; \
		_d[addr( _k1<<1 )+ 1]= _t0i - _ai; \
	}

#define cplx_mul_3_8_B_addsub_store_p(_index,_pd0,_pd1,_t0r,_t0i,_t1r,_t1i) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1r + _t1i )*F_1_8i; \
		_ai = (_t1i - _t1r )*F_1_8i; \
		*( _pd0 + _index) = _t0r + _ar; \
		*( _pd0 + _index + 1) = _t0i + _ai; \
		*( _pd1 + _index) = _t0r - _ar; \
		*( _pd1 + _index + 1) = _t0i - _ai; \
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
#define conv_nested(xkr,xki,xmkr,xmki,ykr,yki,ymkr,ymki,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xkr - xmkr)*0.25;\
   _y0i=(xki + xmki)*0.25;\
   _y1r=(ykr - ymkr);\
   _y1i=(yki + ymki);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xkr * ykr) - (xki * yki);\
   _y0i=(xkr * yki) + (xki * ykr);\
   _y1r=(xmkr * ymkr) - (xmki * ymki);\
   _y1i=(xmkr * ymki) + (xmki * ymkr);\
   _y3r= ((1.0 +(gkr)) * _y2r - (gki) * _y2i);\
   _y3i= ((1.0 +(gkr)) * _y2i + (gki) * _y2r);\
   xkr = _y0r - _y3r;\
   xki = _y0i - _y3i;\
   xmkr = _y1r - _y3r;\
   xmki = _y1i + _y3i;\
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
#define square_nested(xkr,xki,xmkr,xmki,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xkr - xmkr)*0.5;\
   _r1=(xki + xmki)*0.5;\
   _r2=(xkr + xki) * (xkr - xki);\
   xki *= xkr;\
   xkr = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmkr + xmki) * (xmkr - xmki);\
   xmki *= xmkr;\
   xmkr = _r3;\
   _r1 +=  _r1;\
   xki += xki;\
   xmki += xmki;\
   _r0= (((gkr)+ 1.0) * _r2 - (gki) * _r1);\
   _r1= (((gkr)+ 1.0) * _r1 + (gki) * _r2);\
   xkr -= _r0;\
   xmkr -= _r0;\
   xki -= _r1;\
   xmki += _r1;\
}

#define square_nested_eq(xkr,xki,gkr) \
{\
   y_limb_t _y0r=xkr,_y0i=xki;\
   xkr =(_y0r * _y0r + (gkr) * _y0i * _y0i);\
   xki = 2.0 * _y0r * _y0i;\
}

/* These are version to exploid better the symmetries */
#define conv_nested_1_4(xkr,xki,xmkr,xmki,ykr,yki,ymkr,ymki,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xkr - xmkr)*0.25;\
   _y0i=(xki + xmki)*0.25;\
   _y1r=(ykr - ymkr);\
   _y1i=(yki + ymki);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xkr * ykr) - (xki * yki);\
   _y0i=(xkr * yki) + (xki * ykr);\
   _y1r=(xmkr * ymkr) - (xmki * ymki);\
   _y1i=(xmkr * ymki) + (xmki * ymkr);\
   _y3r= ((gki +1.0) * _y2r + gkr * _y2i);\
   _y3i= ((gki +1.0) * _y2i - gkr * _y2r);\
   xkr = _y0r - _y3r;\
   xki = _y0i - _y3i;\
   xmkr = _y1r - _y3r;\
   xmki = _y1i + _y3i;\
}

#define conv_nested_1_2(xkr,xki,xmkr,xmki,ykr,yki,ymkr,ymki,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xkr - xmkr)*0.25;\
   _y0i=(xki + xmki)*0.25;\
   _y1r=(ykr - ymkr);\
   _y1i=(yki + ymki);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xkr * ykr) - (xki * yki);\
   _y0i=(xkr * yki) + (xki * ykr);\
   _y1r=(xmkr * ymkr) - (xmki * ymki);\
   _y1i=(xmkr * ymki) + (xmki * ymkr);\
   _y3r= ((gkr - 1.0) * _y2r - gki * _y2i);\
   _y3i= ((gkr - 1.0) * _y2i + gki * _y2r);\
   xkr = _y0r + _y3r;\
   xki = _y0i + _y3i;\
   xmkr = _y1r + _y3r;\
   xmki = _y1i - _y3i;\
}

#define conv_nested_3_4(xkr,xki,xmkr,xmki,ykr,yki,ymkr,ymki,gkr,gki) \
{ \
   y_limb_t _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i;\
   _y0r=(xkr - xmkr)*0.25;\
   _y0i=(xki + xmki)*0.25;\
   _y1r=(ykr - ymkr);\
   _y1i=(yki + ymki);\
   _y2r=_y0r*_y1r - _y0i*_y1i;\
   _y2i=_y0r*_y1i + _y0i*_y1r;\
   _y0r=(xkr * ykr) - (xki * yki);\
   _y0i=(xkr * yki) + (xki * ykr);\
   _y1r=(xmkr * ymkr) - (xmki * ymki);\
   _y1i=(xmkr * ymki) + (xmki * ymkr);\
   _y3r= ((gki - 1.0) * _y2r + gkr * _y2i);\
   _y3i= ((gki - 1.0) * _y2i - gkr * _y2r);\
   xkr = _y0r + _y3r;\
   xki = _y0i + _y3i;\
   xmkr = _y1r + _y3r;\
   xmki = _y1i - _y3i;\
}

#define square_nested_1_4(xkr,xki,xmkr,xmki,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xkr - xmkr)*0.5;\
   _r1=(xki + xmki)*0.5;\
   _r2=(xkr + xki) * (xkr - xki);\
   xki *= xkr;\
   xkr = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmkr + xmki) * (xmkr - xmki);\
   xmki *= xmkr;\
   xmkr = _r3;\
   _r1 +=  _r1;\
   xki += xki;\
   xmki += xmki;\
   _r0= ((gki + 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki + 1.0) * _r1 - gkr * _r2);\
   xkr -= _r0;\
   xmkr -= _r0;\
   xki -= _r1;\
   xmki += _r1;\
}

#define square_nested_1_2(xkr,xki,xmkr,xmki,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xkr - xmkr)*0.5;\
   _r1=(xki + xmki)*0.5;\
   _r2=(xkr + xki) * (xkr - xki);\
   xki *= xkr;\
   xkr = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmkr + xmki) * (xmkr - xmki);\
   xmki *= xmkr;\
   xmkr = _r3;\
   _r1 +=  _r1;\
   xki += xki;\
   xmki += xmki;\
   _r0= ((gkr - 1.0) * _r2 - gki * _r1);\
   _r1= ((gkr - 1.0) * _r1 + gki * _r2);\
   xkr += _r0;\
   xmkr += _r0;\
   xki += _r1;\
   xmki -= _r1;\
}

#define square_nested_3_4(xkr,xki,xmkr,xmki,gkr,gki) \
{\
   y_limb_t _r0,_r1,_r2,_r3;\
   _r0=(xkr - xmkr)*0.5;\
   _r1=(xki + xmki)*0.5;\
   _r2=(xkr + xki) * (xkr - xki);\
   xki *= xkr;\
   xkr = _r2;\
   _r2 =(_r0 + _r1) * (_r0 - _r1);\
   _r1 *= _r0;\
   _r3 =(xmkr + xmki) * (xmkr - xmki);\
   xmki *= xmkr;\
   xmkr = _r3;\
   _r1 +=  _r1;\
   xki += xki;\
   xmki += xmki;\
   _r0= ((gki - 1.0) * _r2 + gkr * _r1);\
   _r1= ((gki - 1.0) * _r1 - gkr * _r2);\
   xkr += _r0;\
   xmkr += _r0;\
   xki += _r1;\
   xmki -= _r1;\
}


/*$Id$*/










