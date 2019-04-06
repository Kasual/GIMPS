/* $Id$ */
/*  This file is part of
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
/* This are macros to be used with pseudo complex FFT                       */

/* THIS ARE VECTORIZED VERSIONS */

/***************  SINGLE LOAD AND STORE MACROS ***************/




/* To load from real padded array of data as a complex array */
#define cplx_data_to_local_vector(_t,_u,_array,_k) \
	_t##r = _array[addr(( _k )<<1)]; \
	_t##i = _array[addr(( _k )<<1)+1];\
	_u##r = _array[addr(( _k )<<1)+2];\
	_u##i = _array[addr(( _k )<<1)+3];


#define cplx_data_to_local_p_vector(_t,_u,_pd,_array,_k) \
	_pd = _array + addr(( _k)<<1);\
	_t##r = *( _pd ); \
	_t##i = *( _pd + 1);\
	_u##r = *( _pd + 2); \
	_u##i = *( _pd + 3);\
        prefetch_p(_pd);


#define cplx_data_to_local_pp_vector(_t,_u,_index,_pd) \
	_t##r = *( _pd + _index); \
	_t##i = *( _pd + _index + 1);\
	_u##r = *( _pd + _index + 2); \
	_u##i = *( _pd + _index + 3);\
        prefetch_p( _pd + _index);


/* To store on read padded array from local data as complex */

#define cplx_local_to_data_vector(_array,_k,_t,_u) \
	_array[addr(( _k )<<1)]=_t##r ; \
	_array[addr(( _k )<<1)+1]=_t##i; \
	_array[addr(( _k )<<1)+2]=_u##r ; \
	_array[addr(( _k )<<1)+3]=_u##i; \


#define cplx_local_to_data_p_vector(_pd,_t,_u) \
	*( _pd ) = _t##r ; \
	*( _pd + 1) = _t##i;\
	*( _pd + 2) = _u##r ; \
	*( _pd + 3) = _u##i;\
        prefetch_p(_pd);


/* To load from real non-padded array of data as a complex array */
#define cplx_mem_to_local_vector(_t,_u,_array,_k) \
	_t##r =_array[( _k )<<1]; \
	_t##i =_array[(( _k )<<1)+1];\
	_u##r =_array[( _k )<<1]; \
	_u##i =_array[(( _k )<<1)+1];

#define cplx_mem_to_local_p_vector(_t,_u,_pd) \
	_t##r =  *(_pd);\
	_t##i =  *(_pd + 1); \
	_u##r =  *(_pd + 2); \
	_u##i =  *(_pd + 3); \
         prefetch_p (_pd);


/* To strore to a non-padded array from local */
#define cplx_local_to_mem_vector(_array,_k,_t,_u) \
	_array[( _k )<<1]=_t##r; \
	_array[(( _k )<<1)+1]=_t##i;\
	_array[(( _k )<<1)+2]=_u##r; \
	_array[(( _k )<<1)+3]=_u##i;\


/**************** PSEUDO_COMPLEX GENERAL AUXILIAR MACROS **************/
/**********************************************************************
   cplx_muladdsub macro :
		_t = _t1 * _f;
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex and _t is temporal aux var.
***********************************************************************/

# define cplx_muladdsub_vector(_t0,_t1,_u0,_u1,_f1,_f2) \
	{ \
		y_limb_t _ar,_br,_ai,_bi,_a=_t0##r,_b=_u0##r; \
		_ar = (_t1##r ) * (_f1##r ) - (_t1##i ) * (_f1##i ); \
		_br = (_u1##r ) * (_f2##r ) - (_u1##i ) * (_f2##i ); \
		_ai = (_t1##r ) * (_f1##i ) + (_t1##i ) * (_f1##r ); \
		_bi = (_u1##r ) * (_f2##i ) + (_u1##i ) * (_f2##r ); \
		_t0##r += _ar;\
		_u0##r += _br;\
		_t1##r = _a - _ar;\
		_u1##r = _b - _br;\
		_a = _t0##i;\
		_b = _u0##i;\
		_t0##i += _ai;\
		_u0##i += _bi;\
		_t1##i = _a - _ai;\
		_u1##i = _b - _bi;\
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
#define cplx_load_muladdsub_pp_vector_F(_t0,_t1,_u0,_u1,_f,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci,_dr,_di; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
		_dr= *( _pd1 + _index + 2);\
		_di= *( _pd1 + _index + 3);\
                prefetch_p( _pd1 + _index);\
		_ar = (_br) * (_f##r) - (_bi) * (_f##i); \
		_cr = (_dr) * (_f##r) - (_di) * (_f##i); \
		_ai = (_br) * (_f##i) + (_bi) * (_f##r); \
		_ci = (_dr) * (_f##i) + (_di) * (_f##r); \
		_br= *( _pd0 + _index );\
		_bi= *( _pd0 + _index + 1);\
		_dr= *( _pd0 + _index + 2);\
		_di= *( _pd0 + _index + 3);\
                prefetch_p( _pd0 + _index);\
		_t0##r = _br + _ar; \
		_t1##r = _br - _ar; \
		_u0##r = _dr + _cr; \
		_u1##r = _dr - _cr; \
		_t0##i = _bi + _ai; \
		_t1##i = _bi - _ai; \
		_u0##i = _di + _ci; \
		_u1##i = _di - _ci; \
	}


#define cplx_load_muladdsub_pp_vector_B(_t0,_t1,_u0,_u1,_f1,_f2,_index,_pd0,_pd1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi,_cr,_ci,_dr,_di; \
		_br= *( _pd1 + _index );\
		_bi= *( _pd1 + _index + 1);\
		_dr= *( _pd1 + _index + 2);\
		_di= *( _pd1 + _index + 3);\
                prefetch_p( _pd1 + _index);\
		_ar = (_br) * (_f1##r) - (_bi) * (_f1##i); \
		_cr = (_dr) * (_f2##r) - (_di) * (_f2##i); \
		_ai = (_br) * (_f1##i) + (_bi) * (_f1##r); \
		_ci = (_dr) * (_f2##i) + (_di) * (_f2##r); \
		_br= *( _pd0 + _index );\
		_bi= *( _pd0 + _index + 1);\
		_dr= *( _pd0 + _index + 2);\
		_di= *( _pd0 + _index + 3);\
                prefetch_p( _pd0 + _index);\
		_t0##r = _br + _ar; \
		_t1##r = _br - _ar; \
		_u0##r = _dr + _cr; \
		_u1##r = _dr - _cr; \
		_t0##i = _bi + _ai; \
		_t1##i = _bi - _ai; \
		_u0##i = _di + _ci; \
		_u1##i = _di - _ci; \
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
#define cplx_muladdsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1,_f1,_f2) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##r) * (_f1##r) - (_t1##i) * (_f1##i); \
		_br = (_u1##r) * (_f2##r) - (_u1##i) * (_f2##i); \
		_ai = (_t1##r) * (_f1##i) + (_t1##i) * (_f1##r); \
		_bi = (_u1##r) * (_f2##i) + (_u1##i) * (_f2##r); \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 2) = _u0##r + _br; \
		*( _pd1 + _index + 2) = _u0##r - _br; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
		*( _pd0 + _index + 3) = _u0##i + _bi; \
		*( _pd1 + _index + 3) = _u0##i - _bi; \
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
# define cplx_load_mulmuladdsub_pp_vector_F(_t0,_t1,_u0,_u1,_f0,_f1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _a0,_b0,_a1,_b1,_a2,_b2,_a3,_b3,_a4,_b4,_a5,_b5; \
		_a2= *( _pd0 + _index);\
		_a3= *( _pd0 + _index + 1);\
		_b2= *( _pd0 + _index + 2);\
		_b3= *( _pd0 + _index + 3);\
                prefetch_p( _pd0 + _index);\
		_a0 = (_a2) * (_f0##r) - (_a3) * (_f0##i); \
		_b0 = (_b2) * (_f0##r) - (_b3) * (_f0##i); \
		_a1 = (_a2) * (_f0##i) + (_a3) * (_f0##r); \
		_b1 = (_b2) * (_f0##i) + (_b3) * (_f0##r); \
		_a2= *( _pd1 + _index);\
		_a3= *( _pd1 + _index + 1);\
		_b2= *( _pd1 + _index + 2);\
		_b3= *( _pd1 + _index + 3);\
                prefetch_p( _pd1 + _index);\
		_a4 = (_a2) * (_f1##r) - (_a3) * (_f1##i); \
		_b4 = (_b2) * (_f1##r) - (_b3) * (_f1##i); \
		_a5 = (_a2) * (_f1##i) + (_a3) * (_f1##r); \
		_b5 = (_b2) * (_f1##i) + (_b3) * (_f1##r); \
		_t0##r = _a0 + _a4; \
		_t1##r = _a0 - _a4; \
		_u0##r = _b0 + _b4; \
		_u1##r = _b0 - _b4; \
		_t0##i = _a1 + _a5; \
		_t1##i = _a1 - _a5; \
		_u0##i = _b1 + _b5; \
		_u1##i = _b1 - _b5; \
	}


# define cplx_load_mulmuladdsub_pp_vector_B(_t0,_t1,_u0,_u1,_f0,_f1,_g0,_g1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _a0,_b0,_a1,_b1,_a2,_b2,_a3,_b3,_a4,_b4,_a5,_b5; \
		_a2= *( _pd0 + _index);\
		_a3= *( _pd0 + _index + 1);\
		_b2= *( _pd0 + _index + 2);\
		_b3= *( _pd0 + _index + 3);\
                prefetch_p( _pd0 + _index);\
		_a0 = (_a2) * (_f0##r) - (_a3) * (_f0##i); \
		_b0 = (_b2) * (_g0##r) - (_b3) * (_g0##i); \
		_a1 = (_a2) * (_f0##i) + (_a3) * (_f0##r); \
		_b1 = (_b2) * (_g0##i) + (_b3) * (_g0##r); \
		_a2= *( _pd1 + _index);\
		_a3= *( _pd1 + _index + 1);\
		_b2= *( _pd1 + _index + 2);\
		_b3= *( _pd1 + _index + 3);\
                prefetch_p( _pd1 + _index);\
		_a4 = (_a2) * (_f1##r) - (_a3) * (_f1##i); \
		_b4 = (_b2) * (_g1##r) - (_b3) * (_g1##i); \
		_a5 = (_a2) * (_f1##i) + (_a3) * (_f1##r); \
		_b5 = (_b2) * (_g1##i) + (_b3) * (_g1##r); \
		_t0##r = _a0 + _a4; \
		_t1##r = _a0 - _a4; \
		_u0##r = _b0 + _b4; \
		_u1##r = _b0 - _b4; \
		_t0##i = _a1 + _a5; \
		_t1##i = _a1 - _a5; \
		_u0##i = _b1 + _b5; \
		_u1##i = _b1 - _b5; \
	}


/**********************************************************************
   cplx_divmul macro:
		_t0= _f0 / _f1;
		_t1= _f0 * _f1;
   where all the variables  are pseudo-complex. _f0 and _f1 are module
   1 pseudo-complex (like trig. factros)
**********************************************************************/
#define cplx_divmul_vector(_t0,_t1,_u0,_u1,_f0,_f1,_g0,_g1)\
	{\
		y_limb_t _w0,_w1,_w2,_w3,_x0,_x1,_x2,_x3;\
		_w0 = _f0##r * _f1##r;\
		_w1 = _f0##i * _f1##i;\
		_w2 = _f0##i * _f1##r;\
		_w3 = _f0##r * _f1##i;\
		_x0 = _g0##r * _g1##r;\
		_t0##r = _w0 + _w1;\
		_x1 = _g0##i * _g1##i;\
		_t1##r = _w0 - _w1;\
		_x2 = _g0##i * _g1##r;\
		_t0##i = _w2 - _w3;\
		_x3 = _g0##r * _g1##i;\
		_t1##i = _w2 + _w3;\
		_u0##r = _x0 + _x1;\
		_u1##r = _x0 - _x1;\
		_u0##i = _x2 - _x3;\
		_u1##i = _x2 + _x3;\
	}


/**********************************************************************
   cplx_addsub:
		_a= _t0;
		_t0 = _t0 + t1;
		_t1= _a - t1;
   where all the variables are pseudo-complex. _a is temporal aux.
**********************************************************************/
# define cplx_addsub_vector(_t0,_t1,_u0,_u1) \
	{ \
		y_limb_t _a=_t0##r,_b=u0##r; \
		_t0##r += (_t1##r); \
		_u0##r += (_u1##r); \
		_t1##r = _a - (_t1##r); \
		_u1##r = _b - (_u1##r); \
		_a = (_t0##i); \
		_b = (_u0##i); \
		_t0##i += (_t1##i); \
		_u0##i += (_u1##i); \
		_t1##i = _a - (_t1##i); \
		_u1##i = _b - (_u1##i); \
	}


/**********************************************************************
   cplx_load_addsub macro:
		_t0 = _array( _k0) + _array (_k1);
		_t1 = _array( _k0) - _array (_k1);
   where all the variables are pseudo-complex .
**********************************************************************/

#define cplx_load_addsub_pp_vector(_t0,_t1,_u0,_u1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _a0,_a1,_a2,_a3,_b0,_b1,_b2,_b3; \
		_a2 = *( _pd1 + _index);\
		_a3 = *( _pd1 + _index + 1);\
		_b2 = *( _pd1 + _index + 2);\
		_b3 = *( _pd1 + _index + 3);\
                prefetch_p( _pd1 + _index);\
		_a0 = *( _pd0 + _index);\
		_a1 = *( _pd0 + _index + 1);\
		_b0 = *( _pd0 + _index + 2);\
		_b1 = *( _pd0 + _index + 3);\
                prefetch_p( _pd0 + _index);\
		_t0##r = _a0  + _a2; \
		_t1##r = _a0  - _a2; \
		_u0##r = _b0  + _b2; \
		_u1##r = _b0  - _b2; \
		_t0##i = _a1  + _a3; \
		_t1##i = _a1  - _a3; \
		_u0##i = _b1  + _b3; \
		_u1##i = _b1  - _b3; \
	}

#define cplx_load_addsub_pp_vector_no_fetch(_t0,_t1,_u0,_u1,_index,_pd0,_pd1) \
	{ \
		y_limb_t _a0,_a1,_a2,_a3,_b0,_b1,_b2,_b3; \
		_a2 = *( _pd1 + _index);\
		_a3 = *( _pd1 + _index + 1);\
		_b2 = *( _pd1 + _index + 2);\
		_b3 = *( _pd1 + _index + 3);\
		_a0 = *( _pd0 + _index);\
		_a1 = *( _pd0 + _index + 1);\
		_b0 = *( _pd0 + _index + 2);\
		_b1 = *( _pd0 + _index + 3);\
		_t0##r = _a0  + _a2; \
		_t1##r = _a0  - _a2; \
		_u0##r = _b0  + _b2; \
		_u1##r = _b0  - _b2; \
		_t0##i = _a1  + _a3; \
		_t1##i = _a1  - _a3; \
		_u0##i = _b1  + _b3; \
		_u1##i = _b1  - _b3; \
	}



#define cplx_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1)\
	*( _pd0 + _index) = _t0##r + _t1##r;\
	*( _pd0 + _index + 1) = _t0##i + _t1##i;\
	*( _pd0 + _index + 2) = _u0##r + _u1##r;\
	*( _pd0 + _index + 3) = _u0##i + _u1##i;\
	*( _pd1 + _index) = _t0##r - _t1##r;\
	*( _pd1 + _index + 1) = _t0##i - _t1##i;\
	*( _pd1 + _index + 2) = _u0##r - _u1##r;\
	*( _pd1 + _index + 3) = _u0##i - _u1##i;


/***********************************************************************
   cplx_mul_1_4_F macro:
		_t = _t1 * G^-(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
# define cplx_mul_1_4_F_addsub_vector(_t0,_t1,_u0,_u1) \
{ \
	y_limb_t _ar=_t0##r,_ai=_t0##i,_a=_t1##r,_br=_u0##r,_bi=_u0##i,_b=_u1##r; \
	_t0##r += (_t1##i); \
	_u0##r += (_u1##i); \
	_t0##i -= _a; \
	_u0##i -= _b; \
	_t1##r = _ar - (_t1##i); \
	_u1##r = _br - (_u1##i); \
	_t1##i = _ai + _a; \
	_u1##i = _bi + _b; \
}


/***********************************************************************
   cplx_mul_1_4_F_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_F_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	*(_pd0 + _index) = _t0##r + _t1##i; \
	*(_pd0 + _index + 1) = _t0##i - _t1##r; \
	*(_pd0 + _index + 2) = _u0##r + _u1##i; \
	*(_pd0 + _index + 3) = _u0##i - _u1##r; \
	*(_pd1 + _index) = _t0##r - _t1##i; \
	*(_pd1 + _index + 1) = _t0##i + _t1##r;\
	*(_pd1 + _index + 2) = _u0##r - _u1##i; \
	*(_pd1 + _index + 3) = _u0##i + _u1##r;



/***********************************************************************
   cplx_mul_1_4_B macro:
		_t = _t1 * G^(1/4);
		_t0= _t0 + _t;
		_t1= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub_vector(_t0,_t1,_u0,_u1) \
{ \
	y_limb_t _ar=_t0##r, _ai=_t0##i ,_a= _t1##r; \
	y_limb_t _br=_u0##r, _bi=_u0##i ,_b= _u1##r; \
	_t0##r -= (_t1##i); \
	_u0##r -= (_u1##i); \
	_t0##i += _a; \
	_u0##i += _b; \
	_t1##r = _ar + (_t1##i); \
	_u1##r = _br + (_u1##i); \
	_t1##i = _ai - _a; \
	_u1##i = _bi - _b; \
}


/***********************************************************************
   cplx_mul_1_4_B_store macro:
		_t = _t1 * G^-(1/4);
		_d[k0]= _t0 + _t;
		_d[k1]= _t0 - _t;
   where all the variables are pseudo-complex. _t is a temporal aux var and
   G is unit G=exp(2*pi*i), i^2=-1; 
***********************************************************************/
#define cplx_mul_1_4_B_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	*( _pd0 + _index) = _t0##r - _t1##i; \
	*( _pd0 + _index + 1) = _t0##i + _t1##r; \
	*( _pd0 + _index + 2) = _u0##r - _u1##i; \
	*( _pd0 + _index + 3) = _u0##i + _u1##r; \
	*( _pd1 + _index) = _t0##r + _t1##i; \
	*( _pd1 + _index + 1) = _t0##i - _t1##r;\
	*( _pd1 + _index + 2) = _u0##r + _u1##i; \
	*( _pd1 + _index + 3) = _u0##i - _u1##r;


#define cplx_mul_1_8_F_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##r + _t1##i )*F_1_8r; \
		_br = (_u1##r + _u1##i )*F_1_8r; \
		_ai = (_t1##i - _t1##r )*F_1_8r; \
		_bi = (_u1##i - _u1##r )*F_1_8r; \
		*( _pd0 + _index)= _t0##r + _ar; \
		*( _pd1 + _index)= _t0##r - _ar; \
		*( _pd0 + _index + 2)= _u0##r + _br; \
		*( _pd1 + _index + 2)= _u0##r - _br; \
		*( _pd0 + _index + 1)= _t0##i + _ai; \
		*( _pd1 + _index + 1)= _t0##i - _ai; \
		*( _pd0 + _index + 3)= _u0##i + _bi; \
		*( _pd1 + _index + 3)= _u0##i - _bi; \
	}

#define cplx_mul_1_8_B_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##r - _t1##i )*F_1_8r; \
		_br = (_u1##r - _u1##i )*F_1_8r; \
		_ai = (_t1##i + _t1##r )*F_1_8r; \
		_bi = (_u1##i + _u1##r )*F_1_8r; \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 2) = _u0##r + _br; \
		*( _pd1 + _index + 2) = _u0##r - _br; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
		*( _pd0 + _index + 3) = _u0##i + _bi; \
		*( _pd1 + _index + 3) = _u0##i - _bi; \
	}

#define cplx_mul_3_8_F_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##i - _t1##r )*F_1_8r; \
		_br = (_u1##i - _u1##r )*F_1_8r; \
		_ai = (_t1##i + _t1##r )*F_1_8r; \
		_bi = (_u1##i + _u1##r )*F_1_8r; \
		*( _pd0 + _index) = _t0##r + _ar; \
		*( _pd0 + _index + 1) = _t0##i - _ai; \
		*( _pd0 + _index + 2) = _u0##r + _br; \
		*( _pd0 + _index + 3) = _u0##i - _bi; \
		*( _pd1 + _index) = _t0##r - _ar; \
		*( _pd1 + _index + 1) = _t0##i + _ai; \
		*( _pd1 + _index + 2) = _u0##r - _br; \
		*( _pd1 + _index + 3) = _u0##i + _bi; \
	}

#define cplx_mul_3_8_B_addsub_store_p_vector(_index,_pd0,_pd1,_t0,_t1,_u0,_u1) \
	{ \
		y_limb_t _ar,_ai,_br,_bi; \
		_ar = (_t1##r + _t1##i )*F_1_8r; \
		_br = (_u1##r + _u1##i )*F_1_8r; \
		_ai = (_t1##r - _t1##i )*F_1_8r; \
		_bi = (_u1##r - _u1##i )*F_1_8r; \
		*( _pd0 + _index) = _t0##r - _ar; \
		*( _pd0 + _index + 1) = _t0##i + _ai; \
		*( _pd0 + _index + 2) = _u0##r - _br; \
		*( _pd0 + _index + 3) = _u0##i + _bi; \
		*( _pd1 + _index) = _t0##r + _ar; \
		*( _pd1 + _index + 1) = _t0##i - _ai; \
		*( _pd1 + _index + 2) = _u0##r + _br; \
		*( _pd1 + _index + 3) = _u0##i - _bi; \
	}


/* $Id$ */










