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
/*
This includes the norm and carry phase macros for ia64 
 
We supossed compiler uses predication features of this processor 
 
This code is written with Intel ia64 in mind. The cycles in the comments
are, of course, only aproximated. These schedule pattern is what we would like
the compiler to do. But, sure, the compiler can find better solutions 
*/

#ifndef fmax
# define fmax(_x1,_x2) (_x1 > _x2) ? _x1 : _x2
#endif

/* this macro makes carry and norm pass for 2 simul elements */
#  define cplx_carry_norm_check_vector_ia64(_d,_r,_i0,_i1)      \
{                                                               \
  BIG_DOUBLE err0,err1;                                         \
  err0 = _d##_i0##_r * ttmp##_i0 + carry##_i0;                  \
  err1 = _d##_i1##_r * ttmp##_i1 + carry##_i1;                  \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  _d##_i0##_r = RINT(err0);                                     \
  _d##_i1##_r = RINT(err1);                                     \
  err0 -= _d##_i0##_r;                                          \
  err1 -= _d##_i1##_r;                                          \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  err0 = fabs(err0);                                            \
  err1 = fabs(err1);                                            \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  err0 = fmax(err0,err1);                                       \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  maxerr=fmax(maxerr,err0);                                     \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
}

#  define cplx_carry_norm_check_vector_3_ia64(_d,_r,_i0,_i1,_i2)\
{                                                               \
  BIG_DOUBLE err0,err1,err2;                                    \
  err0 = _d##_i0##_r * ttmp##_i0 + carry##_i0;                  \
  err1 = _d##_i1##_r * ttmp##_i1 + carry##_i1;                  \
  err2 = _d##_i2##_r * ttmp##_i2 + carry##_i2;                  \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  _d##_i0##_r = RINT(err0);                                     \
  _d##_i1##_r = RINT(err1);                                     \
  _d##_i2##_r = RINT(err2);                                     \
  err0 -= _d##_i0##_r;                                          \
  err1 -= _d##_i1##_r;                                          \
  err2 -= _d##_i2##_r;                                          \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  err0 = fabs(err0);                                            \
  err1 = fabs(err1);                                            \
  err2 = fabs(err2);                                            \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                                \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  err0 = fmax(err0,err1);                                       \
  maxerr=fmax(maxerr,err2);                                     \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  maxerr = fmax(maxerr,err0);                                   \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
}


# define cplx_carry_norm_check_vector_4_ia64(_d,_r,_i0,_i1,_i2,_i3)  \
{                                                               \
  BIG_DOUBLE err0,err1,err2,err3;                               \
  err0 = _d##_i0##_r * ttmp##_i0 + carry##_i0;                  \
  err1 = _d##_i1##_r * ttmp##_i1 + carry##_i1;                  \
  err2 = _d##_i2##_r * ttmp##_i2 + carry##_i2;                  \
  err3 = _d##_i3##_r * ttmp##_i3 + carry##_i3;                  \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  _d##_i0##_r = RINT(err0);                                     \
  _d##_i1##_r = RINT(err1);                                     \
  _d##_i2##_r = RINT(err2);                                     \
  _d##_i3##_r = RINT(err3);                                     \
  err0 -= _d##_i0##_r;                                          \
  err1 -= _d##_i1##_r;                                          \
  err2 -= _d##_i2##_r;                                          \
  err3 -= _d##_i3##_r;                                          \
  if ( bj##_i0 >= c) _d##_i0##_r *= hiinv;                      \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *= hiinv;                      \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *= hiinv;                      \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *= hiinv;                      \
  else _d##_i3##_r *= loinv;                                    \
  err0 = fabs(err0);                                            \
  err1 = fabs(err1);                                            \
  err2 = fabs(err2);                                            \
  err3 = fabs(err3);                                            \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  err0=fmax(err0,err1);                                         \
  err2=fmax(err2,err3);                                         \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  err0=fmax(err0,err2);                                         \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  maxerr = fmax(maxerr,err0);                                   \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;                                         \
  if (err3>maxerr) maxerr=err3;                                 \
}



#  define cplx_carry_norm_vector_ia64(_d,_r,_i0,_i1)            \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \

#  define cplx_carry_norm_vector_3_ia64(_d,_r,_i0,_i1,_i2)      \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;


#  define cplx_carry_norm_vector_4_ia64(_d,_r,_i0,_i1,_i2,_i3)  \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  _d##_i3##_r = (_d##_i3##_r * ttmp##_i3 + carry##_i3);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  _d##_i3##_r = RINT(_d##_i3##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *=  hiinv;                     \
  else _d##_i3##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;


#  define cplx_carry_norm_vector_5_ia64(_d,_r,_i0,_i1,_i2,_i3,_i4)  \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  _d##_i3##_r = (_d##_i3##_r * ttmp##_i3 + carry##_i3);         \
  _d##_i4##_r = (_d##_i4##_r * ttmp##_i4 + carry##_i4);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  if ( bj##_i4 >= c) ttp##_i4 += ttp##_i4;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  _d##_i3##_r = RINT(_d##_i3##_r);                              \
  _d##_i4##_r = RINT(_d##_i4##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *=  hiinv;                     \
  else _d##_i3##_r *= loinv;                                    \
  if ( bj##_i4 >= c) _d##_i4##_r *=  hiinv;                     \
  else _d##_i4##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  carry##_i4 = RINT(_d##_i4##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  if ( bj##_i4 >= c) ttmp##_i4 *= ttmpBig;                      \
  else ttmp##_i4 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  if ( bj##_i4 >= c) bj##_i4 -= c; else bj##_i4 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  _d##_i4##_r -= carry##_i4;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  _d##_i4##_r *= ttp##_i4;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;                                         \
  ttp##_i4 *= ttpSmall;


#  define cplx_carry_norm_vector_6_ia64(_d,_r,_i0,_i1,_i2,_i3,_i4,_i5) \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  _d##_i3##_r = (_d##_i3##_r * ttmp##_i3 + carry##_i3);         \
  _d##_i4##_r = (_d##_i4##_r * ttmp##_i4 + carry##_i4);         \
  _d##_i5##_r = (_d##_i5##_r * ttmp##_i5 + carry##_i5);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  if ( bj##_i4 >= c) ttp##_i4 += ttp##_i4;                      \
  if ( bj##_i5 >= c) ttp##_i5 += ttp##_i5;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  _d##_i3##_r = RINT(_d##_i3##_r);                              \
  _d##_i4##_r = RINT(_d##_i4##_r);                              \
  _d##_i5##_r = RINT(_d##_i5##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *=  hiinv;                     \
  else _d##_i3##_r *= loinv;                                    \
  if ( bj##_i4 >= c) _d##_i4##_r *=  hiinv;                     \
  else _d##_i4##_r *= loinv;                                    \
  if ( bj##_i5 >= c) _d##_i5##_r *=  hiinv;                     \
  else _d##_i5##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  carry##_i4 = RINT(_d##_i4##_r);                               \
  carry##_i5 = RINT(_d##_i5##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  if ( bj##_i4 >= c) ttmp##_i4 *= ttmpBig;                      \
  else ttmp##_i4 *= ttmpSmall;                                  \
  if ( bj##_i5 >= c) ttmp##_i5 *= ttmpBig;                      \
  else ttmp##_i5 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  if ( bj##_i4 >= c) bj##_i4 -= c; else bj##_i4 += b;           \
  if ( bj##_i5 >= c) bj##_i5 -= c; else bj##_i5 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  _d##_i4##_r -= carry##_i4;                                    \
  _d##_i5##_r -= carry##_i5;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  _d##_i4##_r *= ttp##_i4;                                      \
  _d##_i5##_r *= ttp##_i5;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;                                         \
  ttp##_i4 *= ttpSmall;                                         \
  ttp##_i5 *= ttpSmall;

#  define cplx_carry_norm_vector_7_ia64(_d,_r,_i0,_i1,_i2,_i3,_i4,_i5,_i6) \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  _d##_i3##_r = (_d##_i3##_r * ttmp##_i3 + carry##_i3);         \
  _d##_i4##_r = (_d##_i4##_r * ttmp##_i4 + carry##_i4);         \
  _d##_i5##_r = (_d##_i5##_r * ttmp##_i5 + carry##_i5);         \
  _d##_i6##_r = (_d##_i6##_r * ttmp##_i6 + carry##_i6);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  if ( bj##_i4 >= c) ttp##_i4 += ttp##_i4;                      \
  if ( bj##_i5 >= c) ttp##_i5 += ttp##_i5;                      \
  if ( bj##_i6 >= c) ttp##_i6 += ttp##_i6;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  _d##_i3##_r = RINT(_d##_i3##_r);                              \
  _d##_i4##_r = RINT(_d##_i4##_r);                              \
  _d##_i5##_r = RINT(_d##_i5##_r);                              \
  _d##_i6##_r = RINT(_d##_i6##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *=  hiinv;                     \
  else _d##_i3##_r *= loinv;                                    \
  if ( bj##_i4 >= c) _d##_i4##_r *=  hiinv;                     \
  else _d##_i4##_r *= loinv;                                    \
  if ( bj##_i5 >= c) _d##_i5##_r *=  hiinv;                     \
  else _d##_i5##_r *= loinv;                                    \
  if ( bj##_i6 >= c) _d##_i6##_r *=  hiinv;                     \
  else _d##_i6##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  carry##_i4 = RINT(_d##_i4##_r);                               \
  carry##_i5 = RINT(_d##_i5##_r);                               \
  carry##_i6 = RINT(_d##_i6##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  if ( bj##_i4 >= c) ttmp##_i4 *= ttmpBig;                      \
  else ttmp##_i4 *= ttmpSmall;                                  \
  if ( bj##_i5 >= c) ttmp##_i5 *= ttmpBig;                      \
  else ttmp##_i5 *= ttmpSmall;                                  \
  if ( bj##_i6 >= c) ttmp##_i6 *= ttmpBig;                      \
  else ttmp##_i6 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  if ( bj##_i4 >= c) bj##_i4 -= c; else bj##_i4 += b;           \
  if ( bj##_i5 >= c) bj##_i5 -= c; else bj##_i5 += b;           \
  if ( bj##_i6 >= c) bj##_i6 -= c; else bj##_i6 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  _d##_i4##_r -= carry##_i4;                                    \
  _d##_i5##_r -= carry##_i5;                                    \
  _d##_i6##_r -= carry##_i6;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  _d##_i4##_r *= ttp##_i4;                                      \
  _d##_i5##_r *= ttp##_i5;                                      \
  _d##_i6##_r *= ttp##_i6;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;                                         \
  ttp##_i4 *= ttpSmall;                                         \
  ttp##_i5 *= ttpSmall;                                         \
  ttp##_i6 *= ttpSmall;

#  define cplx_carry_norm_vector_8_ia64(_d,_r,_i0,_i1,_i2,_i3,_i4,_i5,_i6,_i7) \
  _d##_i0##_r = (_d##_i0##_r * ttmp##_i0 + carry##_i0);         \
  _d##_i1##_r = (_d##_i1##_r * ttmp##_i1 + carry##_i1);         \
  _d##_i2##_r = (_d##_i2##_r * ttmp##_i2 + carry##_i2);         \
  _d##_i3##_r = (_d##_i3##_r * ttmp##_i3 + carry##_i3);         \
  _d##_i4##_r = (_d##_i4##_r * ttmp##_i4 + carry##_i4);         \
  _d##_i5##_r = (_d##_i5##_r * ttmp##_i5 + carry##_i5);         \
  _d##_i6##_r = (_d##_i6##_r * ttmp##_i6 + carry##_i6);         \
  _d##_i7##_r = (_d##_i7##_r * ttmp##_i7 + carry##_i7);         \
  if ( bj##_i0 >= c) ttp##_i0 += ttp##_i0;                      \
  if ( bj##_i1 >= c) ttp##_i1 += ttp##_i1;                      \
  if ( bj##_i2 >= c) ttp##_i2 += ttp##_i2;                      \
  if ( bj##_i3 >= c) ttp##_i3 += ttp##_i3;                      \
  if ( bj##_i4 >= c) ttp##_i4 += ttp##_i4;                      \
  if ( bj##_i5 >= c) ttp##_i5 += ttp##_i5;                      \
  if ( bj##_i6 >= c) ttp##_i6 += ttp##_i6;                      \
  if ( bj##_i7 >= c) ttp##_i7 += ttp##_i7;                      \
  _d##_i0##_r = RINT(_d##_i0##_r);                              \
  _d##_i1##_r = RINT(_d##_i1##_r);                              \
  _d##_i2##_r = RINT(_d##_i2##_r);                              \
  _d##_i3##_r = RINT(_d##_i3##_r);                              \
  _d##_i4##_r = RINT(_d##_i4##_r);                              \
  _d##_i5##_r = RINT(_d##_i5##_r);                              \
  _d##_i6##_r = RINT(_d##_i6##_r);                              \
  _d##_i7##_r = RINT(_d##_i7##_r);                              \
  if ( bj##_i0 >= c) _d##_i0##_r *=  hiinv;                     \
  else _d##_i0##_r *= loinv;                                    \
  if ( bj##_i1 >= c) _d##_i1##_r *=  hiinv;                     \
  else _d##_i1##_r *= loinv;                                    \
  if ( bj##_i2 >= c) _d##_i2##_r *=  hiinv;                     \
  else _d##_i2##_r *= loinv;                                    \
  if ( bj##_i3 >= c) _d##_i3##_r *=  hiinv;                     \
  else _d##_i3##_r *= loinv;                                    \
  if ( bj##_i4 >= c) _d##_i4##_r *=  hiinv;                     \
  else _d##_i4##_r *= loinv;                                    \
  if ( bj##_i5 >= c) _d##_i5##_r *=  hiinv;                     \
  else _d##_i5##_r *= loinv;                                    \
  if ( bj##_i6 >= c) _d##_i6##_r *=  hiinv;                     \
  else _d##_i6##_r *= loinv;                                    \
  if ( bj##_i7 >= c) _d##_i7##_r *=  hiinv;                     \
  else _d##_i7##_r *= loinv;                                    \
  carry##_i0 = RINT(_d##_i0##_r);                               \
  carry##_i1 = RINT(_d##_i1##_r);                               \
  carry##_i2 = RINT(_d##_i2##_r);                               \
  carry##_i3 = RINT(_d##_i3##_r);                               \
  carry##_i4 = RINT(_d##_i4##_r);                               \
  carry##_i5 = RINT(_d##_i5##_r);                               \
  carry##_i6 = RINT(_d##_i6##_r);                               \
  carry##_i7 = RINT(_d##_i7##_r);                               \
  if ( bj##_i0 >= c) ttmp##_i0 *= ttmpBig;                      \
  else ttmp##_i0 *= ttmpSmall;                                  \
  if ( bj##_i1 >= c) ttmp##_i1 *= ttmpBig;                      \
  else ttmp##_i1 *= ttmpSmall;                                  \
  if ( bj##_i2 >= c) ttmp##_i2 *= ttmpBig;                      \
  else ttmp##_i2 *= ttmpSmall;                                  \
  if ( bj##_i3 >= c) ttmp##_i3 *= ttmpBig;                      \
  else ttmp##_i3 *= ttmpSmall;                                  \
  if ( bj##_i4 >= c) ttmp##_i4 *= ttmpBig;                      \
  else ttmp##_i4 *= ttmpSmall;                                  \
  if ( bj##_i5 >= c) ttmp##_i5 *= ttmpBig;                      \
  else ttmp##_i5 *= ttmpSmall;                                  \
  if ( bj##_i6 >= c) ttmp##_i6 *= ttmpBig;                      \
  else ttmp##_i6 *= ttmpSmall;                                  \
  if ( bj##_i7 >= c) ttmp##_i7 *= ttmpBig;                      \
  else ttmp##_i7 *= ttmpSmall;                                  \
  if ( bj##_i0 >= c) bj##_i0 -= c; else bj##_i0 += b;           \
  if ( bj##_i1 >= c) bj##_i1 -= c; else bj##_i1 += b;           \
  if ( bj##_i2 >= c) bj##_i2 -= c; else bj##_i2 += b;           \
  if ( bj##_i3 >= c) bj##_i3 -= c; else bj##_i3 += b;           \
  if ( bj##_i4 >= c) bj##_i4 -= c; else bj##_i4 += b;           \
  if ( bj##_i5 >= c) bj##_i5 -= c; else bj##_i5 += b;           \
  if ( bj##_i6 >= c) bj##_i6 -= c; else bj##_i6 += b;           \
  if ( bj##_i7 >= c) bj##_i7 -= c; else bj##_i7 += b;           \
  _d##_i0##_r -= carry##_i0;                                    \
  _d##_i1##_r -= carry##_i1;                                    \
  _d##_i2##_r -= carry##_i2;                                    \
  _d##_i3##_r -= carry##_i3;                                    \
  _d##_i4##_r -= carry##_i4;                                    \
  _d##_i5##_r -= carry##_i5;                                    \
  _d##_i6##_r -= carry##_i6;                                    \
  _d##_i7##_r -= carry##_i7;                                    \
  _d##_i0##_r *= ttp##_i0;                                      \
  _d##_i1##_r *= ttp##_i1;                                      \
  _d##_i2##_r *= ttp##_i2;                                      \
  _d##_i3##_r *= ttp##_i3;                                      \
  _d##_i4##_r *= ttp##_i4;                                      \
  _d##_i5##_r *= ttp##_i5;                                      \
  _d##_i6##_r *= ttp##_i6;                                      \
  _d##_i7##_r *= ttp##_i7;                                      \
  ttp##_i0 *= ttpSmall;                                         \
  ttp##_i1 *= ttpSmall;                                         \
  ttp##_i2 *= ttpSmall;                                         \
  ttp##_i3 *= ttpSmall;                                         \
  ttp##_i4 *= ttpSmall;                                         \
  ttp##_i5 *= ttpSmall;                                         \
  ttp##_i6 *= ttpSmall;                                         \
  ttp##_i7 *= ttpSmall;


