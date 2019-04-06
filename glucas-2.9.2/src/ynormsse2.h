/*$Id$*/
/*
   (c) 2003-2006 Guillermo Ballester Valor 
   
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
 
   This include SSE2 macros
*/
/* To store with add on read padded array from local data as complex */
/* generic target */

#define sse2_local_to_data_add( _array, _k, _t)                      \
  {                                                                  \
     Y__M128D _auxr, _auxi;                                          \
        Y_MM_LOAD_PAIR_NEXT ( _auxr,_auxi, _array + addr((_k)<<1) ); \
        Y_MM_ADD_PD( _auxr, _auxr, _t##r);                           \
        Y_MM_ADD_PD( _auxi, _auxi, _t##i);                           \
        Y_MM_STORE_PAIR_NEXT ( _array + addr((_k)<<1), _auxr, _auxi);\
  }

#define sse2_local_to_data_inter_add( _array, _k, _t)                \
  {                                                                  \
     Y__M128D auxr, auxi;                                             \
        Y_MM_LOAD_INTER_PAIR_NEXT ( auxr, auxi, _array + addr((_k)<<1) );\
        Y_MM_ADD_PD( auxr, auxr, _t##r);                           \
        Y_MM_ADD_PD( auxi, auxi, _t##i);                           \
        Y_MM_STORE_INTER_PAIR_NEXT ( _array + addr((_k)<<1), auxr, auxi);\
  }

#define sse2_local_to_data_add_p( _array, _pd, _t)      \
  {                                                      \
     Y__M128D auxr, auxi;                                \
        Y_MM_LOAD_PAIR_NEXT ( auxr, auxi, _array + _pd );\
        Y_MM_ADD_PD( auxr, auxr, _t##r);                 \
        Y_MM_ADD_PD( auxi, auxi, _t##i);                 \
        Y_MM_STORE_PAIR_NEXT ( _array + _pd, auxr, auxi);\
  }

#if defined(SUM_CHECK)

# define cplx_carry_norm_check_sse2(_d,_j)                         \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}


#  define cplx_carry_norm_check_sse2_last(_d,_j)                   \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0HI);      \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_check_sse2_last_odd(_d,_j)               \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0LO);      \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

# define cplx_carry_norm_check_sse2_vector(_dj,_dk,_j,_k)          \
{                                                                  \
  Y__M128D maskj, maskk;                                           \
  int imaskj, imaskk;                                              \
  Y_MM_ADD_PD ( summ, summ, _dj);                                  \
  Y_MM_MUL_PD ( maskj, _dj, ttmp##_j );                            \
  Y_MM_MUL_PD ( maskk, _dk, ttmp##_k );                            \
  Y_MM_ADD_PD ( summ, summ, _dk);                                  \
  Y_MM_TWO_RINT_PD ( _dj, _dk, maskj, maskk);                      \
  Y_MM_SUB_PD ( maskj, maskj, _dj);                                \
  Y_MM_SUB_PD ( maskk, maskk, _dk);                                \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_ABS_PD ( maskk, maskk);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _dj, _dj, carry##_j);                              \
  Y_MM_MAX_PD ( maxerr, maxerr, maskk);                            \
  Y_MM_ADD_PD ( _dk, _dk, carry##_k);                              \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_CMPGT_PD ( maskk , bj##_k, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MOVEMASK_PD ( imaskk, maskk);                               \
  Y_MM_MUL_PD ( _dj, _dj, MM_inv[imaskj]);                         \
  Y_MM_MUL_PD ( _dk, _dk, MM_inv[imaskk]);                         \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_AND_PD ( maskk, maskk, ttp##_k);                            \
  Y_MM_TWO_RINT_PD ( carry##_j, carry##_k, _dj, _dk);              \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_ADD_PD ( ttp##_k, ttp##_k, maskk);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_MUL_PD ( ttmp##_k, ttmp##_k, MM_auxt[imaskk]);              \
  Y_MM_SUB_PD ( _dj, _dj, carry##_j );                             \
  Y_MM_SUB_PD ( _dk, _dk, carry##_k );                             \
  Y_MM_MUL_PD ( _dj, _dj, ttp##_j );                               \
  Y_MM_MUL_PD ( _dk, _dk, ttp##_k );                               \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_ADD_PD ( bj##_k, bj##_k, MM_bc[imaskk]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
  Y_MM_MUL_PD ( ttp##_k, ttp##_k, MM_Hsmall );                     \
}

#else

# define cplx_carry_norm_check_sse2(_d,_j)                         \
{                                                                  \
  Y__M128D maskj, auxj;                                            \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_CMPGT_PD ( auxj , bj##_j, MM_c );                           \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_MOVEMASK_PD ( imaskj, auxj);                                \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_AND_PD ( auxj, auxj, ttp##_j);                              \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, auxj);                           \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}



#  define cplx_carry_norm_check_sse2_last(_d,_j)                   \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0HI);      \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}


#  define cplx_carry_norm_check_sse2_last_odd(_d,_j)               \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_PD ( _d, maskj);                                       \
  Y_MM_SUB_PD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0LO);      \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

# define cplx_carry_norm_check_sse2_vector(_dj,_dk,_j,_k)          \
{                                                                  \
  Y__M128D maskj, maskk;                                           \
  int imaskj, imaskk;                                              \
  Y_MM_MUL_PD ( maskj, _dj, ttmp##_j );                            \
  Y_MM_MUL_PD ( maskk, _dk, ttmp##_k );                            \
  Y_MM_TWO_RINT_PD ( _dj, _dk, maskj, maskk);                      \
  Y_MM_SUB_PD ( maskj, maskj, _dj);                                \
  Y_MM_SUB_PD ( maskk, maskk, _dk);                                \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_ABS_PD ( maskk, maskk);                                     \
  Y_MM_MAX_PD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_PD ( _dj, _dj, carry##_j);                              \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MAX_PD ( maxerr, maxerr, maskk);                            \
  Y_MM_ADD_PD ( _dk, _dk, carry##_k);                              \
  Y_MM_CMPGT_PD ( maskk , bj##_k, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MOVEMASK_PD ( imaskk, maskk);                               \
  Y_MM_MUL_PD ( _dj, _dj, MM_inv[imaskj]);                         \
  Y_MM_MUL_PD ( _dk, _dk, MM_inv[imaskk]);                         \
  Y_MM_TWO_RINT_PD ( carry##_j, carry##_k, _dj, _dk);              \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_AND_PD ( maskk, maskk, ttp##_k);                            \
  Y_MM_SUB_PD ( _dj, _dj, carry##_j );                             \
  Y_MM_SUB_PD ( _dk, _dk, carry##_k );                             \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_ADD_PD ( ttp##_k, ttp##_k, maskk);                          \
  Y_MM_MUL_PD ( ttmp##_k, ttmp##_k, MM_auxt[imaskk]);              \
  Y_MM_MUL_PD ( _dj, _dj, ttp##_j );                               \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( _dk, _dk, ttp##_k );                               \
  Y_MM_ADD_PD ( bj##_k, bj##_k, MM_bc[imaskk]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
  Y_MM_MUL_PD ( ttp##_k, ttp##_k, MM_Hsmall );                     \
}

#endif /*sumcheck*/

#if defined(SUM_CHECK)

# define cplx_carry_norm_sse2(_d,_j)                               \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_PD ( _d, _d);                                          \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last(_d,_j)                         \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_PD ( _d, _d);                                          \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0HI);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last_odd(_d,_j)                     \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_PD ( _d, _d);                                          \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0LO);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_vector(_dj,_dk,_j,_k)               \
{                                                                  \
  Y__M128D maskj, maskk;                                           \
  int imaskj, imaskk;                                              \
  Y_MM_ADD_PD ( summ, summ, _dj);                                  \
  Y_MM_MUL_PD ( maskj, _dj, ttmp##_j );                            \
  Y_MM_ADD_PD ( summ, summ, _dk);                                  \
  Y_MM_MUL_PD ( maskk, _dk, ttmp##_k );                            \
  Y_MM_TWO_RINT_PD ( _dj, _dk, maskj, maskk);                      \
  Y_MM_ADD_PD ( _dj, _dj, carry##_j);                              \
  Y_MM_ADD_PD ( _dk, _dk, carry##_k);                              \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_CMPGT_PD ( maskk , bj##_k, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MOVEMASK_PD ( imaskk, maskk);                               \
  Y_MM_MUL_PD ( _dj, _dj, MM_inv[imaskj]);                         \
  Y_MM_MUL_PD ( _dk, _dk, MM_inv[imaskk]);                         \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_AND_PD ( maskk, maskk, ttp##_k);                            \
  Y_MM_TWO_RINT_PD ( carry##_j, carry##_k, _dj, _dk);              \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_ADD_PD ( ttp##_k, ttp##_k, maskk);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_MUL_PD ( ttmp##_k, ttmp##_k, MM_auxt[imaskk]);              \
  Y_MM_SUB_PD ( _dj, _dj, carry##_j );                             \
  Y_MM_SUB_PD ( _dk, _dk, carry##_k );                             \
  Y_MM_MUL_PD ( _dj, _dj, ttp##_j );                               \
  Y_MM_MUL_PD ( _dk, _dk, ttp##_k );                               \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_ADD_PD ( bj##_k, bj##_k, MM_bc[imaskk]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
  Y_MM_MUL_PD ( ttp##_k, ttp##_k, MM_Hsmall );                     \
}

#else

# define cplx_carry_norm_sse2(_d,_j)                               \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_RINT_EQ_PD ( _d);                                           \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last(_d,_j)                         \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_EQ_PD ( _d);                                           \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0HI);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last_odd(_d,_j)                     \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_PD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_EQ_PD ( _d);                                           \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0LO);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  if defined(Y_AMD64)
#   define cplx_carry_norm_sse2_vector(_dj,_dk,_j,_k)               \
{                                                                  \
  Y__M128D maskj, maskk;                                           \
  int imaskj, imaskk;                                              \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_CMPGT_PD ( maskk , bj##_k, MM_c );                          \
  Y_MM_MUL_PD ( _dj, _dj, ttmp##_j );                              \
  Y_MM_MUL_PD ( _dk, _dk, ttmp##_k );                              \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MOVEMASK_PD ( imaskk, maskk);                               \
  Y_MM_ADD_PD ( _dj, _dj, carry##_j);                              \
  Y_MM_ADD_PD ( _dk, _dk, carry##_k);                              \
  Y_MM_TWO_RINT_EQ_PD ( _dj, _dk);                                 \
  Y_MM_MUL_PD ( _dj, _dj, MM_inv[imaskj]);                         \
  Y_MM_MUL_PD ( _dk, _dk, MM_inv[imaskk]);                         \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_AND_PD ( maskk, maskk, ttp##_k);                            \
  Y_MM_TWO_RINT_PD ( carry##_j, carry##_k, _dj, _dk);              \
  Y_MM_SUB_PD ( _dj, _dj, carry##_j );                             \
  Y_MM_SUB_PD ( _dk, _dk, carry##_k );                             \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_ADD_PD ( ttp##_k, ttp##_k, maskk);                          \
  Y_MM_MUL_PD ( _dj, _dj, ttp##_j );                               \
  Y_MM_MUL_PD ( _dk, _dk, ttp##_k );                               \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_MUL_PD ( ttmp##_k, ttmp##_k, MM_auxt[imaskk]);              \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
  Y_MM_ADD_PD ( bj##_k, bj##_k, MM_bc[imaskk]);                    \
  Y_MM_MUL_PD ( ttp##_k, ttp##_k, MM_Hsmall );                     \
}

#  else

#   define cplx_carry_norm_sse2_vector(_dj,_dk,_j,_k)              \
{                                                                  \
  Y__M128D maskj, maskk;                                           \
  int imaskj, imaskk;                                              \
  Y_MM_MUL_PD ( maskj, _dj, ttmp##_j );                            \
  Y_MM_MUL_PD ( maskk, _dk, ttmp##_k );                            \
  Y_MM_TWO_RINT_PD ( _dj, _dk, maskj, maskk);                      \
  Y_MM_ADD_PD ( _dj, _dj, carry##_j);                              \
  Y_MM_ADD_PD ( _dk, _dk, carry##_k);                              \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_CMPGT_PD ( maskk , bj##_k, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MOVEMASK_PD ( imaskk, maskk);                               \
  Y_MM_MUL_PD ( _dj, _dj, MM_inv[imaskj]);                         \
  Y_MM_MUL_PD ( _dk, _dk, MM_inv[imaskk]);                         \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_AND_PD ( maskk, maskk, ttp##_k);                            \
  Y_MM_TWO_RINT_PD ( carry##_j, carry##_k, _dj, _dk);              \
  Y_MM_SUB_PD ( _dj, _dj, carry##_j );                             \
  Y_MM_SUB_PD ( _dk, _dk, carry##_k );                             \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_ADD_PD ( ttp##_k, ttp##_k, maskk);                          \
  Y_MM_MUL_PD ( _dj, _dj, ttp##_j );                               \
  Y_MM_MUL_PD ( _dk, _dk, ttp##_k );                               \
  Y_MM_MUL_PD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_MUL_PD ( ttmp##_k, ttmp##_k, MM_auxt[imaskk]);              \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_ADD_PD ( bj##_k, bj##_k, MM_bc[imaskk]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
  Y_MM_MUL_PD ( ttp##_k, ttp##_k, MM_Hsmall );                     \
}
# endif /* Y_AMD64 */
#endif /*sumcheck */

/* SD version for ODD radices */
#if defined(SUM_CHECK)

# define cplx_carry_norm_check_sse2_sd(_d,_j)                      \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_SD ( summ, summ, _d);                                   \
  Y_MM_MUL_SD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_SD ( _d, maskj);                                       \
  Y_MM_SUB_SD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_SD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_sd(_d,_j)                           \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_SD ( summ, summ, _d);                                   \
  Y_MM_MUL_SD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_SD ( _d, _d);                                          \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_check_sse2_last_odd_sd(_d,_j)            \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_SD ( summ, summ, _d);                                   \
  Y_MM_MUL_SD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_SD ( _d, maskj);                                       \
  Y_MM_SUB_SD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_SD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0LO);      \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last_odd_sd(_d,_j)                  \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_ADD_PD ( summ, summ, _d);                                   \
  Y_MM_MUL_SD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_SD ( _d, _d);                                          \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0LO);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#else

# define cplx_carry_norm_check_sse2_sd(_d,_j)                      \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_SD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_SD ( _d, maskj);                                       \
  Y_MM_SUB_SD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_SD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_sd(_d,_j)                           \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_SD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_SD ( _d, _d);                                          \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_check_sse2_last_odd_sd(_d,_j)            \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_SD ( maskj, _d, ttmp##_j );                             \
  Y_MM_RINT_SD ( _d, maskj);                                       \
  Y_MM_SUB_SD ( maskj, maskj, _d);                                 \
  Y_MM_ABS_PD ( maskj, maskj);                                     \
  Y_MM_MAX_SD ( maxerr, maxerr, maskj);                            \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  if( l == (pad - 2) ) Y_MM_AND_PD( bj##_j, bj##_j, MM_Y0LO);      \
  Y_MM_CMPGT_SD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#  define cplx_carry_norm_sse2_last_odd_sd(_d,_j)                  \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_MUL_SD ( _d, _d, ttmp##_j );                                \
  Y_MM_RINT_SD ( _d, _d);                                          \
  Y_MM_ADD_SD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  if( l == (pad - 2) ) Y_MM_AND_PD( maskj, maskj, MM_Y0LO);        \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_SD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_SD ( carry##_j, _d);                                   \
  Y_MM_ADD_SD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_MUL_SD ( ttmp##_j, ttmp##_j, MM_auxt[imaskj]);              \
  Y_MM_SUB_SD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_SD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_SD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_SD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

#endif /*sumcheck */

#  define cplx_carry_sse2(_d,_j)                                   \
{                                                                  \
  Y__M128D maskj;                                                  \
  int imaskj;                                                      \
  Y_MM_RINT_PD ( _d, _d);                                          \
  Y_MM_ADD_PD ( _d, _d, carry##_j);                                \
  Y_MM_CMPGT_PD ( maskj , bj##_j, MM_c );                          \
  Y_MM_MOVEMASK_PD ( imaskj, maskj);                               \
  Y_MM_MUL_PD ( _d, _d, MM_inv[imaskj]);                           \
  Y_MM_AND_PD ( maskj, maskj, ttp##_j);                            \
  Y_MM_RINT_PD ( carry##_j, _d);                                   \
  Y_MM_ADD_PD ( ttp##_j, ttp##_j, maskj);                          \
  Y_MM_SUB_PD ( _d, _d, carry##_j );                               \
  Y_MM_MUL_PD ( _d, _d, ttp##_j );                                 \
  Y_MM_ADD_PD ( bj##_j, bj##_j, MM_bc[imaskj]);                    \
  Y_MM_MUL_PD ( ttp##_j, ttp##_j, MM_Hsmall );                     \
}

/*$Id$*/
