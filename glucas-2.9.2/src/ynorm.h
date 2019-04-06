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
/* To store with add on read padded array from local data as complex */
#if (Y_TARGET == 1 || Y_TARGET == 11 || Y_TARGET == 12 || Y_TARGET == 13)
# ifndef Y_KILL_BRANCHES
#  define Y_KILL_BRANCHES
# endif
# if Y_TARGET == 1
#  include "ynormx86.h"
# else
#  include "ynormx86p3.h"
# endif
#else

/*#define Y_ALT_MAX*/
/* generic target */
#ifdef Y_ALT_MAX
# define check_err(_maxerr,_err) \
  _err = (_err - _maxerr) * 0.5 ;\
  _maxerr += _err + fabs(_err);\

#else

# define check_err(_maxerr,_err) _maxerr = (_maxerr >= _err) ? _maxerr : _err;

#endif /* Y_ALT_MAX */

# define cplx_local_to_data_add(_array, _k, _t)   \
	_array[addr(( _k ) << 1)] += _t##r ;      \
	_array[addr(( _k ) << 1) + 1] += _t##i;

# define cplx_local_to_data_add_p(_array, _pd, _t) \
	_array[_pd ] += _t##r ;                    \
	_array[_pd + 1] += _t##i;


# define cplx_muladdsub_store_add(_array, _k0, _k1, _t0, _t1, _f) \
	{                                                         \
		y_limb_t _ar,_ai;                                 \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i);    \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r);    \
		_array[addr(_k0<<1)]+= _t0##r + _ar;              \
                _array[addr(_k1<<1)]+= _t0##r - _ar;              \
		_array[addr(_k0<<1)+1]+= _t0##i + _ai;            \
                _array[addr(_k1<<1)+1]+= _t0##i - _ai;            \
        }

#if defined(_OPENMP) || defined(_PTHREADS) || defined(_SUNMP)
# define cplx_muladdsub_store_add_mp(_d0, _d1, _t0, _t1, _f)   \
	{                                                      \
		y_limb_t _ar,_ai;                              \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i); \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r); \
		*(_d0) += _t0##r + _ar;                        \
                *(_d1) += _t0##r - _ar;                        \
		*(_d0 + 1) += _t0##i + _ai;                    \
                *(_d1 + 1) += _t0##i - _ai;                    \
        }
#endif

# ifndef Y_KILL_BRANCHES

/* This macro makes the carry and norm task for a set _j */
/*
#  define cplx_carry_norm_check(_d, _j)      \
{                                            \
  BIG_DOUBLE temp0, tempErr, err;            \
  temp0 = _d##_j##r;                         \
  tempErr = RINT( temp0 * ttmp##_j );        \
  err = fabs( temp0 * ttmp##_j - tempErr);   \
  check_err (maxerr, err);                   \
  temp0 = tempErr + carry##_j;               \
  if ( bj##_j >= c)                          \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT (temp0);              \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##r = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
  temp0 = _d##_j##i;                         \
  tempErr = RINT( temp0 * ttmp##_j );        \
  err = fabs( temp0 * ttmp##_j - tempErr);   \
  check_err (maxerr, err);                   \
   temp0 = tempErr + carry##_j ;             \
  if ( bj##_j >= c)                          \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##i = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
}
*/

#  define cplx_carry_norm_check(_d, _j)          \
{                                                \
  BIG_DOUBLE tempErr, err;                       \
  tempErr = carry##_j + _d##_j##r * ttmp##_j;    \
  _d##_j##r = RINT( tempErr);                    \
  err = fabs( _d##_j##r - tempErr);              \
  if ( bj##_j >= c)                              \
    {                                            \
      _d##_j##r *= hiinv;                        \
      ttp##_j += ttp##_j;                        \
      carry##_j = RINT (_d##_j##r);              \
      bj##_j -= c;                               \
      ttmp##_j *= ttmpBig;                       \
    }                                            \
  else                                           \
    {                                            \
      _d##_j##r *= loinv;                        \
      carry##_j = RINT(_d##_j##r);               \
      bj##_j += b;                               \
      ttmp##_j *= ttmpSmall;                     \
    }                                            \
  check_err (maxerr, err);                       \
  _d##_j##r = (_d##_j##r - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                           \
  tempErr = carry##_j + _d##_j##i * ttmp##_j;    \
  _d##_j##i = RINT( tempErr);                    \
  err = fabs( _d##_j##i - tempErr);              \
  if ( bj##_j >= c)                              \
    {                                            \
      _d##_j##i *= hiinv;                        \
      ttp##_j += ttp##_j;                        \
      carry##_j = RINT(_d##_j##i);               \
      bj##_j -= c;                               \
      ttmp##_j *= ttmpBig;                       \
    }                                            \
  else                                           \
    {                                            \
      _d##_j##i *= loinv;                        \
      carry##_j = RINT(_d##_j##i);               \
      bj##_j += b;                               \
      ttmp##_j *= ttmpSmall;                     \
    }                                            \
  check_err (maxerr, err);                       \
  _d##_j##i = (_d##_j##i - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                           \
}

#  define cplx_carry_norm_check_vector(_d, _j, _k)              \
{                                                               \
  BIG_DOUBLE tempErr0, err0, tempErr1, err1;                    \
  tempErr0 = RINT( _d##_j##r * ttmp##_j );                      \
  tempErr1 = RINT( _d##_k##r * ttmp##_k );                      \
  err0 = fabs( _d##_j##r * ttmp##_j - tempErr0);                \
  err1 = fabs( _d##_k##r * ttmp##_k - tempErr1);                \
  if (err0 > maxerr) maxerr = err0;                             \
  if (err1 > maxerr) maxerr = err1;                             \
  _d##_j##r = tempErr0 + carry##_j;                             \
  _d##_k##r = tempErr1 + carry##_k;                             \
  if ( bj##_j >= c) _d##_j##r *= hiinv; else _d##_j##r *= loinv;\
  if ( bj##_k >= c) _d##_k##r *= hiinv; else _d##_k##r *= loinv;\
  if ( bj##_j >= c) ttp##_j += ttp##_j;                         \
  if ( bj##_k >= c) ttp##_k += ttp##_k;                         \
  carry##_j = RINT (_d##_j##r);                                 \
  carry##_k = RINT (_d##_k##r);                                 \
  if ( bj##_j >= c) ttmp##_j *= ttmpBig;                        \
  else ttmp##_j *= ttmpSmall;                                   \
  if ( bj##_k >= c) ttmp##_k *= ttmpBig;                        \
  else ttmp##_k *= ttmpSmall;                                   \
  if ( bj##_j >= c) bj##_j -= c; else bj##_j += b;              \
  if ( bj##_k >= c) bj##_k -= c; else bj##_k += b;              \
  _d##_j##r = (_d##_j##r - carry##_j ) * ttp##_j;               \
  _d##_k##r = (_d##_k##r - carry##_k ) * ttp##_k;               \
  ttp##_j *= ttpSmall;                                          \
  ttp##_k *= ttpSmall;                                          \
  tempErr0 = RINT( _d##_j##i * ttmp##_j );                      \
  tempErr1 = RINT( _d##_k##i * ttmp##_k );                      \
  err0 = fabs( _d##_j##i * ttmp##_j - tempErr0);                \
  err1 = fabs( _d##_k##i * ttmp##_k - tempErr1);                \
  if (err0>maxerr) maxerr=err0;                                 \
  if (err1>maxerr) maxerr=err1;                                 \
  _d##_j##i = tempErr0 + carry##_j;                             \
  _d##_k##i = tempErr1 + carry##_k;                             \
  if ( bj##_j >= c) _d##_j##i *= hiinv; else _d##_j##i *= loinv;\
  if ( bj##_k >= c) _d##_k##i *= hiinv; else _d##_k##i *= loinv;\
  if ( bj##_j >= c) ttp##_j += ttp##_j;                         \
  if ( bj##_k >= c) ttp##_k += ttp##_k;                         \
  carry##_j = RINT (_d##_j##i);                                 \
  carry##_k = RINT (_d##_k##i);                                 \
  if ( bj##_j >= c) ttmp##_j *= ttmpBig;                        \
  else ttmp##_j *= ttmpSmall;                                   \
  if ( bj##_k >= c) ttmp##_k *= ttmpBig;                        \
  else ttmp##_k *= ttmpSmall;                                   \
  if ( bj##_j >= c) bj##_j -= c; else bj##_j += b;              \
  if ( bj##_k >= c) bj##_k -= c; else bj##_k += b;              \
  _d##_j##i = (_d##_j##i - carry##_j ) * ttp##_j;               \
  _d##_k##i = (_d##_k##i - carry##_k ) * ttp##_k;               \
  ttp##_j *= ttpSmall;                                          \
  ttp##_k *= ttpSmall;                                          \
}


#  define cplx_carry_norm(_d, _j)                \
  _d##_j##r *= ttmp##_j;                         \
  _d##_j##r = carry##_j + RINT( _d##_j##r );     \
  if ( bj##_j >= c)                              \
    {                                            \
      _d##_j##r *= hiinv;                        \
      ttp##_j += ttp##_j;                        \
      carry##_j = RINT(_d##_j##r);               \
      bj##_j -= c;                               \
      ttmp##_j *= ttmpBig;                       \
    }                                            \
  else                                           \
    {                                            \
      _d##_j##r *= loinv;                        \
      carry##_j = RINT(_d##_j##r);               \
      bj##_j += b;                               \
      ttmp##_j *= ttmpSmall;                     \
    }                                            \
  _d##_j##r = (_d##_j##r - carry##_j ) * ttp##_j;\
  _d##_j##i *= ttmp##_j;                         \
  ttp##_j *= ttpSmall;                           \
  _d##_j##i = carry##_j + RINT( _d##_j##i );     \
  if ( bj##_j >= c)                              \
    {                                            \
      _d##_j##i *= hiinv;                        \
      carry##_j = RINT(_d##_j##i);               \
      ttp##_j += ttp##_j;                        \
      bj##_j -= c;                               \
      ttmp##_j *= ttmpBig;                       \
    }                                            \
  else                                           \
    {                                            \
      _d##_j##i *= loinv;                        \
      carry##_j = RINT(_d##_j##i);               \
      bj##_j += b;                               \
      ttmp##_j *= ttmpSmall;                     \
    }                                            \
  _d##_j##i = (_d##_j##i - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;


#  define cplx_carry_norm_vector(_d, _j, _k)                    \
  _d##_j##r =  carry##_j + RINT( _d##_j##r * ttmp##_j );        \
  _d##_k##r =  carry##_k + RINT( _d##_k##r * ttmp##_k );        \
  if ( bj##_j >= c) _d##_j##r *= hiinv; else _d##_j##r *= loinv;\
  if ( bj##_k >= c) _d##_k##r *= hiinv; else _d##_k##r *= loinv;\
  if ( bj##_j >= c) ttp##_j += ttp##_j;                         \
  if ( bj##_k >= c) ttp##_k += ttp##_k;                         \
  carry##_j = RINT (_d##_j##r);                                 \
  carry##_k = RINT (_d##_k##r);                                 \
  if ( bj##_j >= c) ttmp##_j *= ttmpBig;                        \
  else ttmp##_j *= ttmpSmall;                                   \
  if ( bj##_k >= c) ttmp##_k *= ttmpBig;                        \
  else ttmp##_k *= ttmpSmall;                                   \
  if ( bj##_j >= c) bj##_j -= c; else bj##_j += b;              \
  if ( bj##_k >= c) bj##_k -= c; else bj##_k += b;              \
  _d##_j##r = (_d##_j##r - carry##_j ) * ttp##_j;               \
  _d##_k##r = (_d##_k##r - carry##_k ) * ttp##_k;               \
  ttp##_j *= ttpSmall;                                          \
  ttp##_k *= ttpSmall;                                          \
  _d##_j##i =  carry##_j + RINT( _d##_j##i * ttmp##_j );        \
  _d##_k##i =  carry##_k + RINT( _d##_k##i * ttmp##_k );        \
  if ( bj##_j >= c) _d##_j##i *= hiinv; else _d##_j##i *= loinv;\
  if ( bj##_k >= c) _d##_k##i *= hiinv; else _d##_k##i *= loinv;\
  if ( bj##_j >= c) ttp##_j += ttp##_j;                         \
  if ( bj##_k >= c) ttp##_k += ttp##_k;                         \
  carry##_j = RINT (_d##_j##i);                                 \
  carry##_k = RINT (_d##_k##i);                                 \
  if ( bj##_j >= c) ttmp##_j *= ttmpBig;                        \
  else ttmp##_j *= ttmpSmall;                                   \
  if ( bj##_k >= c) ttmp##_k *= ttmpBig;                        \
  else ttmp##_k *= ttmpSmall;                                   \
  if ( bj##_j >= c) bj##_j -= c; else bj##_j += b;              \
  if ( bj##_k >= c) bj##_k -= c; else bj##_k += b;              \
  _d##_j##i = (_d##_j##i - carry##_j ) * ttp##_j;               \
  _d##_k##i = (_d##_k##i - carry##_k ) * ttp##_k;               \
  ttp##_j *= ttpSmall;                                          \
  ttp##_k *= ttpSmall;



/* This macro is for last carrie phase */
#  define cplx_carry( _d, _j)                \
{                                            \
  BIG_DOUBLE temp0;                          \
  temp0 = _d##_j##r;                         \
  if ( bj##_j >= c )                         \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
    }                                        \
  _d##_j##r = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
  temp0 = carry##_j ;                        \
  if ( bj##_j >= c )                         \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
    }                                        \
  _d##_j##i = (temp0 - carry##_j ) * ttp##_j;\
  if(carry##_j != 0.0)                       \
    {                                        \
	printf(" Wrong carries: %f \n", carry##_j); \
	exit(1);                             \
    }                                        \
}

/* This macro only has diferent with respect cplx_carry_norm because it
take into account if l==lastloop */
#  define cplx_carry_norm_last_check(_d,_j)  \
{                                            \
  BIG_DOUBLE temp0,tempErr,err;              \
  temp0= _d##_j##r ;                         \
  tempErr = RINT( temp0 * ttmp##_j );        \
  err = fabs( temp0 * ttmp##_j - tempErr);   \
  if (err>maxerr) maxerr=err;                \
  temp0 = tempErr + carry##_j;               \
  if ( bj##_j >=c )                          \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##r = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
  if(l == (pad-1)) bj##_j = 0;               \
  temp0= _d##_j##i;                          \
  tempErr = RINT( temp0 * ttmp##_j );        \
  err = fabs( temp0 * ttmp##_j - tempErr);   \
  if (err>maxerr) maxerr=err;                \
  temp0 = tempErr + carry##_j ;              \
  if ( bj##_j >= c )                         \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##i = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
}

#  define cplx_carry_norm_last(_d, _j)       \
{                                            \
  BIG_DOUBLE temp0, tempErr;                 \
  tempErr = RINT( _d##_j##r * ttmp##_j );    \
  temp0 = tempErr + carry##_j;               \
  if ( bj##_j >=c )                          \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##r = (temp0 - carry##_j ) * ttp##_j;\
  if(l == (pad-1)) bj##_j = 0;               \
  ttp##_j *= ttpSmall;                       \
  tempErr = RINT( _d##_j##i * ttmp##_j );    \
  temp0 = tempErr + carry##_j ;              \
  if ( bj##_j >= c )                         \
    {                                        \
      temp0 *= hiinv;                        \
      ttp##_j += ttp##_j;                    \
      carry##_j = RINT(temp0);               \
      bj##_j -= c;                           \
      ttmp##_j *= ttmpBig;                   \
    }                                        \
  else                                       \
    {                                        \
      temp0 *= loinv;                        \
      carry##_j = RINT(temp0);               \
      bj##_j += b;                           \
      ttmp##_j *= ttmpSmall;                 \
    }                                        \
  _d##_j##i = (temp0 - carry##_j ) * ttp##_j;\
  ttp##_j *= ttpSmall;                       \
}

# else


/* This macro makes the carry and norm task for a set _j */
#  define cplx_carry_norm_check(_d, _j)            \
{                                                  \
  BIG_DOUBLE temp0, err;                           \
  temp0 = _d##_j##r * tt##_j[1] + carry##_j;       \
  _d##_j##r = RINT( temp0 );                       \
  err = fabs( temp0  - _d##_j##r);                 \
  check_err(maxerr,err);                           \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);       \
  _d##_j##r *= aux[issmall];                       \
  bj##_j += ib[issmall];                           \
  tt##_j[issmall] += tt##_j[issmall];              \
  carry##_j = RINT( _d##_j##r);                    \
  tt##_j[1] *= aux[2];                             \
  _d##_j##r = (_d##_j##r - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                             \
                                                   \
  temp0 = _d##_j##i * tt##_j[1] + carry##_j;       \
  _d##_j##i = RINT( temp0 );                       \
  err = fabs( temp0  - _d##_j##i);                 \
  check_err(maxerr,err);                           \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);       \
  _d##_j##i *= aux[issmall];                       \
  bj##_j += ib[issmall];                           \
  tt##_j[issmall] += tt##_j[issmall];              \
  carry##_j = RINT( _d##_j##i);                    \
  tt##_j[1] *= aux[2];                             \
  _d##_j##i = (_d##_j##i - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                             \
}


/* This macro makes the carry and norm task for a set of two _j's */
#  define cplx_carry_norm_check_vector(_d, _j0, _j1)   \
{                                                      \
  BIG_DOUBLE temp0, temp1, err;                        \
  temp0 = _d##_j0##r * tt##_j0[1] + carry##_j0;        \
  temp1 = _d##_j1##r * tt##_j1[1] + carry##_j1;        \
  _d##_j0##r = RINT( temp0 );                          \
  _d##_j1##r = RINT( temp1 );                          \
  err = fabs( temp0 - _d##_j0##r);                     \
  check_err(maxerr,err);                               \
  err = fabs( temp1 - _d##_j1##r);                     \
  check_err(maxerr,err);                               \
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);          \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);         \
  _d##_j0##r *= aux[issmall];                          \
  _d##_j1##r *= aux[issmall1];                         \
  bj##_j0 += ib[issmall];                              \
  bj##_j1 += ib[issmall1];                             \
  tt##_j0[issmall] += tt##_j0[issmall];                \
  tt##_j1[issmall1] += tt##_j1[issmall1];              \
  carry##_j0 = RINT(_d##_j0##r);                       \
  carry##_j1 = RINT(_d##_j1##r);                       \
  tt##_j0[1] *= aux[2];                                \
  tt##_j1[1] *= aux[2];                                \
  _d##_j0##r = (_d##_j0##r - carry##_j0 ) * tt##_j0[0];\
  _d##_j1##r = (_d##_j1##r - carry##_j1 ) * tt##_j1[0];\
  tt##_j0[0] *= aux[3];                                \
  tt##_j1[0] *= aux[3];                                \
                                                       \
  temp0 = _d##_j0##i * tt##_j0[1] + carry##_j0;        \
  temp1 = _d##_j1##i * tt##_j1[1] + carry##_j1;        \
  _d##_j0##i = RINT( temp0 );                          \
  _d##_j1##i = RINT( temp1 );                          \
  err = fabs( temp0 - _d##_j0##i);                     \
  check_err(maxerr,err);                               \
  err = fabs( temp1 - _d##_j1##i);                     \
  check_err(maxerr,err);                               \
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);          \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);         \
  _d##_j0##i *= aux[issmall];                          \
  _d##_j1##i *= aux[issmall1];                         \
  bj##_j0 += ib[issmall];                              \
  bj##_j1 += ib[issmall1];                             \
  tt##_j0[issmall] += tt##_j0[issmall];                \
  tt##_j1[issmall1] += tt##_j1[issmall1];              \
  carry##_j0 = RINT(_d##_j0##i);                       \
  carry##_j1 = RINT(_d##_j1##i);                       \
  tt##_j0[1] *= aux[2];                                \
  tt##_j1[1] *= aux[2];                                \
  _d##_j0##i = (_d##_j0##i - carry##_j0 ) * tt##_j0[0];\
  _d##_j1##i = (_d##_j1##i - carry##_j1 ) * tt##_j1[0];\
  tt##_j0[0] *= aux[3];                                \
  tt##_j1[0] *= aux[3];                                \
}



#  define cplx_carry_norm(_d, _j)              \
{                                              \
  BIG_DOUBLE temp0, tempErr;                   \
  tempErr = RINT( _d##_j##r * tt##_j[1] );     \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
                                               \
  tempErr = RINT( _d##_j##i * tt##_j[1]);      \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
}

#  define cplx_carry_norm_vector(_d,_j0,_j1)      \
{                                                 \
  BIG_DOUBLE temp0, temp1, tempErr0, tempErr1;    \
  tempErr0 = RINT( _d##_j0##r * tt##_j0[1] );     \
  tempErr1 = RINT( _d##_j1##r * tt##_j1[1] );     \
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);     \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);    \
  temp0 = (tempErr0 + carry##_j0) * aux[issmall]; \
  temp1 = (tempErr1 + carry##_j1) * aux[issmall1];\
  bj##_j0 += ib[issmall];                         \
  bj##_j1 += ib[issmall1];                        \
  tt##_j0[issmall] += tt##_j0[issmall];           \
  tt##_j1[issmall1] += tt##_j1[issmall1];         \
  carry##_j0 = RINT(temp0);                       \
  carry##_j1 = RINT(temp1);                       \
  tt##_j0[1] *= aux[2];                           \
  tt##_j1[1] *= aux[2];                           \
  _d##_j0##r = (temp0 - carry##_j0 ) * tt##_j0[0];\
  _d##_j1##r = (temp1 - carry##_j1 ) * tt##_j1[0];\
  tt##_j0[0] *= aux[3];                           \
  tt##_j1[0] *= aux[3];                           \
                                                  \
  tempErr0 = RINT( _d##_j0##i * tt##_j0[1]);      \
  tempErr1 = RINT( _d##_j1##i * tt##_j1[1]);      \
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);     \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);    \
  temp0 = (tempErr0 + carry##_j0) * aux[issmall]; \
  temp1 = (tempErr1 + carry##_j1) * aux[issmall1];\
  bj##_j0 += ib[issmall];                         \
  bj##_j1 += ib[issmall1];                        \
  tt##_j0[issmall] += tt##_j0[issmall];           \
  tt##_j1[issmall1] += tt##_j1[issmall1];         \
  carry##_j0 = RINT(temp0);                       \
  carry##_j1 = RINT(temp1);                       \
  tt##_j0[1] *= aux[2];                           \
  tt##_j1[1] *= aux[2];                           \
  _d##_j0##i = (temp0 - carry##_j0 ) * tt##_j0[0];\
  _d##_j1##i = (temp1 - carry##_j1 ) * tt##_j1[0];\
  tt##_j0[0] *= aux[3];                           \
  tt##_j1[0] *= aux[3];                           \
}

/*
#  define cplx_carry_norm(_d, _j)                       \
  _d##_j##r = RINT( carry##_j + _d##_j##r * tt##_j[1] );\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);            \
  _d##_j##r *= aux[issmall];                            \
  bj##_j += ib[issmall];                                \
  tt##_j[issmall] += tt##_j[issmall];                   \
  carry##_j = RINT( _d##_j##r);                         \
  tt##_j[1] *= aux[2];                                  \
  _d##_j##r = (_d##_j##r - carry##_j ) * tt##_j[0];     \
  tt##_j[0] *= aux[3];                                  \
                                                        \
  _d##_j##i = RINT( carry##_j + _d##_j##i * tt##_j[1] );\
  issmall = (bj##_j - c) >> (BITS_PER_UL - 1);          \
  _d##_j##i *= aux[issmall];                            \
  bj##_j += ib[issmall];                                \
  tt##_j[issmall] += tt##_j[issmall];                   \
  carry##_j = RINT( _d##_j##i);                         \
  tt##_j[1] *= aux[2];                                  \
  _d##_j##i = (_d##_j##i - carry##_j ) * tt##_j[0];     \
  tt##_j[0] *= aux[3];                                  \


#  define cplx_carry_norm_vector(_d,_j0,_j1)                \
  _d##_j0##r = RINT( carry##_j0 + _d##_j0##r * tt##_j0[1] );\
  _d##_j1##r = RINT( carry##_j1 + _d##_j1##r * tt##_j1[1] );\
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);               \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);              \
  _d##_j0##r *= aux[issmall];                               \
  _d##_j1##r *= aux[issmall1];                              \
  bj##_j0 += ib[issmall];                                   \
  bj##_j1 += ib[issmall1];                                  \
  tt##_j0[issmall] += tt##_j0[issmall];                     \
  tt##_j1[issmall1] += tt##_j1[issmall1];                   \
  carry##_j0 = RINT( _d##_j0##r);                           \
  carry##_j1 = RINT( _d##_j1##r);                           \
  tt##_j0[1] *= aux[2];                                     \
  tt##_j1[1] *= aux[2];                                     \
  _d##_j0##r = (_d##_j0##r - carry##_j0 ) * tt##_j0[0];     \
  _d##_j1##r = (_d##_j1##r - carry##_j1 ) * tt##_j1[0];     \
  tt##_j0[0] *= aux[3];                                     \
  tt##_j1[0] *= aux[3];                                     \
                                                            \
  _d##_j0##i = RINT( carry##_j0 + _d##_j0##i * tt##_j0[1] );\
  _d##_j1##i = RINT( carry##_j1 + _d##_j1##i * tt##_j1[1] );\
  issmall=(bj##_j0 - c) >> (BITS_PER_UL - 1);               \
  issmall1=(bj##_j1 - c) >> (BITS_PER_UL - 1);              \
  _d##_j0##i *= aux[issmall];                               \
  _d##_j1##i *= aux[issmall1];                              \
  bj##_j0 += ib[issmall];                                   \
  bj##_j1 += ib[issmall1];                                  \
  tt##_j0[issmall] += tt##_j0[issmall];                     \
  tt##_j1[issmall1] += tt##_j1[issmall1];                   \
  carry##_j0 = RINT( _d##_j0##i);                           \
  carry##_j1 = RINT( _d##_j1##i);                           \
  tt##_j0[1] *= aux[2];                                     \
  tt##_j1[1] *= aux[2];                                     \
  _d##_j0##i = (_d##_j0##i - carry##_j0 ) * tt##_j0[0];     \
  _d##_j1##i = (_d##_j1##i - carry##_j1 ) * tt##_j1[0];     \
  tt##_j0[0] *= aux[3];                                     \
  tt##_j1[0] *= aux[3];                                     \
                                                            \
*/

/* This macro is for last carrie phase */
#  define cplx_carry(_d, _j)                       \
{                                                  \
  BIG_DOUBLE temp0;                                \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);       \
  temp0 = _d##_j##r * aux[issmall];                \
  bj##_j += ib[issmall];                           \
  tt##_j[issmall] += tt##_j[issmall];              \
  carry##_j = RINT(temp0);                         \
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];    \
  tt##_j[0] *= aux[3];                             \
                                                   \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);       \
  temp0 = carry##_j * aux[issmall];                \
  carry##_j = RINT(temp0);                         \
  tt##_j[issmall] += tt##_j[issmall];              \
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];    \
  if(carry##_j != 0.0)                             \
    {                                              \
	printf(" Carries too high final pass \n"); \
	exit(1);                                   \
    }                                              \
}

/* This macro only has diferent with respect cplx_carry_norm because it
take into account if l==lastloop */
#  define cplx_carry_norm_last_check(_d, _j)   \
{                                              \
  BIG_DOUBLE temp0, tempErr, err;              \
  temp0 = _d##_j##r;                           \
  tempErr = RINT( temp0 * tt##_j[1]);          \
  err = fabs( temp0 * tt##_j[1] - tempErr);    \
  check_err(maxerr,err);                       \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
                                               \
  temp0 = _d##_j##i;                           \
  tempErr = RINT( temp0 * tt##_j[1]);          \
  err = fabs( temp0 * tt##_j[1] - tempErr);    \
  check_err(maxerr,err);                       \
  bj##_j &= nolast;                            \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
}

#  define cplx_carry_norm_last(_d, _j)         \
{                                              \
  BIG_DOUBLE temp0,tempErr;                    \
  tempErr = RINT( _d##_j##r * tt##_j[1] );     \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
                                               \
  tempErr = RINT( _d##_j##i * tt##_j[1] );     \
  bj##_j &= nolast;                            \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);   \
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];                       \
  tt##_j[issmall] += tt##_j[issmall];          \
  carry##_j = RINT(temp0);                     \
  tt##_j[1] *= aux[2];                         \
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];                         \
}

# endif

#endif

/*$Id$*/






