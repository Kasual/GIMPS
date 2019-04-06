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

#if defined(USE_ASM) && defined(__CYGWIN__)
/* Cygwin's GCC uses the variable exactly as given in asm, but exports with _ prefix in C */
#define MEMORY_VARIABLE_NAME(name) "_" #name
#else
#define MEMORY_VARIABLE_NAME(name) #name
#endif

/* To store with add on read padded array from local data as complex */

#ifdef USE_ASM

# define asm_real_carry_norm_check(_bj,_data,_carry,_tt)                \
__asm__ volatile("                                          \n"         \
"  fldl %2			#1	temp0=datareal      \n"         \
"  fmull 8(%%edi)		#2-4    TEMP0*TT[1]=TEMP0   \n"         \
"  movl %%eax,%%edx		#3                          \n"         \
"  subl " MEMORY_VARIABLE_NAME(c) ",%%edx                       #4                          \n"         \
"  fld  %%st(0)			#5	temp0,temp0         \n"         \
"  faddl 32(%%esi)              #6-8    TEMP0+BIGA,temp0    \n"         \
"  shrl $31,%%edx		#7                          \n"         \
"  addl (%%ebx,%%edx,4),%%eax   #8                          \n"         \
"  fsubl 32(%%esi)		#9-11	TEMPERR,temp0       \n"         \
"  fldl (%%edi,%%edx,8)         #10	tt[ismall],TEMPERR,temp0 \n"    \
"  fadd %%st(0)			#11-13  2*TT[issmall],TEMPERR,temp0 \n" \
"  fxch %%st(1)			#11	TEMPERR,2*TT[issmall],temp0 \n" \
"  fsub %%st(0),%%st(2)		#12-14  temperr,2*TT[issmall],temperr-temp0 \n"   \
"  faddl %3			#13-15  TEMPERR+CARRY,2*TT[issmall],TEMPERR-TEMP0 \n"\
"  fxch %%st(1)			#14	2*tt[issmall],TEMPERR+CARRY,TEMPERR-TEMP0 \n"\
"  fstpl (%%edi,%%edx,8)        #15     TEMP+CARRY,temperr-temp0 \n"    \
"  fmull (%%esi,%%edx,8)	#16-18	INV*(temp+carry),temperr-temp0 \n" \
"  fxch %%st(1)			#16     temperr-temp0,NTMP0 \n"         \
"  fabs				#17     err,NTEMP0          \n"         \
"  fsubl %4			#18-20  ERR-MAXERR=ERR,NTEMP0 \n"       \
"  fld  %%st(1)			#19     temp0,ERR,temp0     \n"         \
"  faddl 32(%%esi)		#20-22	TEMP0+BIGA,ERR,temp0 \n"        \
"  fxch %%st(1)			#20	ERR,TEMP0+BIGA,temp0 \n"        \
"  fmull 40(%%esi)		#21-23  0.5*ERR=ERR,TEMP0+BIGA,temp0 \n"\
"  fldl 8(%%edi)		#22	TT[1],ERR,TEMP0+BIGA,temp0 \n"  \
"  fmull 16(%%esi)		#23-25  ttmpBig*TT[1]=TT[1],ERR,temp0+biga,temp0 \n"\
"  fxch %%st(2)			#23	temp0+biga,ERR,TT[1],temp0 \n"  \
"  fsubl 32(%%esi)		#24-26	NEWCARRY,ERR,TT[1],temp0 \n"    \
"  fld  %%st(1)			#25     err,NEWCARRY,err,TT[1],temp0 \n"\
"  fabs 			#26     abs(err),NEWCARRY,err,tt[1],temp0 \n" \
"  faddp %%st(0),%%st(2)	#27-29  newcarry,ERR,tt[1],temp0 \n"    \
"  fstl %3			#28	newcarry,ERR,tt[1],temp0 \n"    \
"  fsubrp %%st(0),%%st(3)	#29-31  ERR,tt[1],TEMP0-NEWCARRY \n"    \
"  faddl %4                     #30-32  MAXERR,tt[1],TEMP0-NEWCARRY \n" \
"  fldl (%%edi)			#31	tt[0],MAXERR,tt[1],TEMP0-NEWCARRY \n"\
"  fmul %%st(0),%%st(3)         #32-34  TT[0],MAXERR,tt[1],NEWDATA  \n" \
"                               #33     stall                    \n"    \
"  fmull 24(%%esi)		#34-36  TT[0]*ttpSmall,MAXERR,tt[1],NEWDATA \n" \
"  fxch %%st(2)			#35	tt[1],maxerr,TT[0],NEWDATA \n"  \
"  fstpl 8(%%edi)		#36	maxerr,tt[0],newdata \n"        \
"  fstpl %4			#37     tt[0],newdata        \n"        \
"  fstpl (%%edi)		#38	newdata              \n"        \
"  fstpl %2			#39 	void                 \n"        \
  : "=&a" (_bj) : "0" (_bj),"m"(_data),"m"(_carry),                     \
 "m"(maxerr),"D"(&_tt[0]),"S"(&aux[0]),"b"(&ib[0]):"memory","dx")

# define asm_real_carry_norm(_bj,_data,_carry,_tt)                      \
__asm__ volatile("                                           \n"        \
"  fldl	%2			#1	temp0=data           \n"        \
"  fmull 8(%%edi)		#2-4    TEMP0*TT[1]=TEMP0    \n"        \
"  movl %%eax,%%edx		#3                           \n"        \
"  subl " MEMORY_VARIABLE_NAME(c) ",%%edx                       #4                           \n"        \
"  faddl 32(%%esi)		#5-7	TEMP0+BIGA           \n"        \
"  shrl $31,%%edx          	#6                           \n"        \
"  addl (%%ebx,%%edx,4),%%eax   #7                           \n"        \
"  fsubl 32(%%esi)		#8-10   TEMPERR              \n"        \
"  fldl (%%edi,%%edx,8)		#9	tt[issmall],TEMPERR  \n"        \
"  fadd  %%st(0)		#10-12  2*TT[issmall]=TT[issmall],TEMPERR\n"\
"  fxch  %%st(1)		#10	TEMPERR,TT[issmall]  \n"        \
"  faddl %3			#11-13	TEMPERR+CARRY,TT[issmall] \n"   \
"  fxch %%st(1)			#12	TT[issmall],TEMPERR+CARRY \n"   \
"  fstpl (%%edi,%%edx,8)	#13     TEMPERR+CARRY        \n"        \
"  fmull (%%esi,%%edx,8)	#14-16	inv*(temperr+CARRY)=TEMP0 \n"   \
"  fldl	8(%%edi)		#15	tt[1],TEMP0          \n"        \
"  fmull 16(%%esi)		#16-18  ttmpBig*tt[1],TEMP0  \n"        \
"  fld  %%st(1)			#17     temp0,TT[1],temp0    \n"        \
"  faddl 32(%%esi)		#18-20	TEMP0+BIGA,TT[1],temp0 \n"      \
"  fldl (%%edi)			#19	tt[0],TEMP0+BIGa,TT[1],temp0 \n"\
"  fxch %%st(2)			#20	tt[1],TEMP0+BIGA,tt[0],temp0 \n"\
"  fstpl 8(%%edi)		#21	temp0+biga,tt[0],temp0 \n"      \
"  fsubl 32(%%esi)		#22-24	NEWCARRY,tt[0],temp0 \n"        \
"  fxch %%st(1)			#22	tt[0],NEWCARRY,temp0 \n"        \
"  fmull 24(%%esi)		#23-25	NEWTT[0],NEWCARRY,temp0 \n"     \
"  fxch %%st(2)			#23	temp0,NEWCARRY,NEWTT[0] \n"     \
"  fsub %%st(1)			#24	TEMP0-NEWCARRY,NEWCARRY,NEWTT[0]\n"\
"  fxch %%st(1)                 #25     NEWCARRY,TEMP0-NEWCARRY,NEWTT[0]\n"\
"  fstpl %3                     #26     TEMP0-NEWCARRY,NEWTT[0] \n"     \
"  fmull (%%edi)		#25-27	NEWDATA,NEWtt[0]     \n"        \
"  fxch %%st(1)			#26     newtt[0],NEWDATA     \n"        \
"  fstpl (%%edi)		#27     NEWDATA              \n"        \
"  fstpl %2			#28     void                 \n"        \
     :"=&a" (_bj) :"0" (_bj),"m"(_data),"m"(_carry),                    \
     "D"(&_tt[0]),"S"(&aux[0]),"b"(&ib[0]):"memory","dx")

#define asm_real_carry(_bj,_data,_carry,_tt)                          \
__asm__ volatile("                                              \n"   \
"movl %%eax,%%edx                 #1                             \n"  \
"subl " MEMORY_VARIABLE_NAME(c) ", %%edx                    #2                             \n"  \
"shrl $31,%%edx                   #3                             \n"  \
"fldl %2                          #4       datare                \n"  \
"fmull (%%esi,%%edx,8)            #5-8     datare*xinv[is]       \n"  \
"addl (%%ebx,%%edx,4),%%eax       #6                             \n"  \
"fldl (%%edi,%%edx,8)             #7       tt[is],datare*xinv[is] \n" \
"fadd %%st(0)                     #8-11    2*tt[is],datare*xinv[is]\n"\
"fld  %%st(1)                     #9       temp0,2*tt[is],temp0  \n"  \
"faddl 32(%%esi)                  #10-12   temp0+bigA,ntt[is],temp0 \n"\
"fxch %%st(1)                     #11      ntt[is],temp0+bigA,temp0 \n"\
"fstpl (%%edi,%%edx,8)            #12      temp0+bigA,temp0      \n"  \
"fsubl 32(%%esi)                  #13-15   carry,temp0           \n"  \
"fldl (%%edi)                     #14      tt[0],carry,temp0     \n"  \
"fxch %%st(2)                     #14      temp0,carry,tt[0]     \n"  \
"                                 #15      STALL                 \n"  \
"fsub %%st(1)                     #16-18   temp0-carry,carry,tt[0]\n" \
"fxch %%st(1)                     #17      carry,temp0-carry,tt[0]\n" \
"fstpl %3                         #18      temp0-carry,tt[0]     \n"  \
"fmul %%st(1)                     #19-21   newdata,tt[0]         \n"  \
"fxch %%st(1)                     #20      tt[0],newdata         \n"  \
"fmull 24(%%esi)                  #21-22   ntt[0],newdata        \n"  \
"fxch %%st(1)                     #21      newdata,ntt[0]        \n"  \
"fstpl %2                         #22      ntt[0]                \n"  \
"fstpl (%%edi)                    #23                            \n"  \
     :"=&a"(_bj) :"0"(_bj),"m"(_data),"m"(_carry),                    \
     "D"(&_tt[0]),"S"(&aux[0]),"b"(&ib[0]):"memory","dx");

#endif


#define cplx_local_to_data_add(_array,_k,_t) \
	_array[addr(( _k )<<1)] += _t##r ; \
	_array[addr(( _k )<<1)+1] += _t##i;

#define cplx_local_to_data_add_p(_array,_pd,_t) \
	_array[_pd ] += _t##r ; \
	_array[_pd + 1] += _t##i;


#define cplx_muladdsub_store_add(_array,_k0,_k1,_t0,_t1,_f) \
	{ \
		y_limb_t _ar,_ai; \
		_ar = (_t1##r) * (_f##r) - (_t1##i) * (_f##i); \
		_ai = (_t1##r) * (_f##i) + (_t1##i) * (_f##r); \
		_array[addr(_k0<<1)]+= _t0##r + _ar; \
                _array[addr(_k1<<1)]+= _t0##r - _ar; \
		_array[addr(_k0<<1)+1]+= _t0##i + _ai; \
                _array[addr(_k1<<1)+1]+= _t0##i - _ai; \
        }

#ifdef Y_ALTERNATIVE_MAX
# define check_err(_maxerr,_err)\
  _err = (_err - _maxerr) * aux[5];\
  _maxerr += (_err + fabs(_err));
#else
# define check_err(_maxerr,_err) _maxerr = (_maxerr >= _err) ? _maxerr : _err;
#endif


/* This macro makes the carry and norm task for a set _j */
#ifndef USE_ASM

# define cplx_carry_norm_check(_d,_j)\
{\
  BIG_DOUBLE temp0,tempErr,err;\
  temp0 = _d##_j##r;\
  tempErr = RINT( temp0 * tt##_j[1]); \
  err = fabs( temp0 * tt##_j[1] - tempErr);\
  check_err(maxerr,err);\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
  \
  temp0 = _d##_j##i;\
  tempErr = RINT( temp0 * tt##_j[1]); \
  err = fabs( temp0 * tt##_j[1] - tempErr);\
  check_err(maxerr,err);\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
}

# define cplx_carry_norm_check_vector(_d,_j0,_j1)\
  cplx_carry_norm_check(_d,_j0);\
  cplx_carry_norm_check(_d,_j1);

#else

# define cplx_carry_norm_check(_d,_j)\
{\
	asm_real_carry_norm_check(bj##_j ,_d##_j##r ,carry##_j ,tt##_j );\
	asm_real_carry_norm_check(bj##_j ,_d##_j##i ,carry##_j ,tt##_j );\
}

#endif

#ifndef USE_ASM
# define cplx_carry_norm(_d,_j)\
{\
  BIG_DOUBLE temp0,tempErr;\
  tempErr = RINT( _d##_j##r * tt##_j[1] ); \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
  \
  tempErr = RINT( _d##_j##i * tt##_j[1]); \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
}
#else
# define cplx_carry_norm(_d,_j)\
{\
	asm_real_carry_norm(bj##_j ,_d##_j##r ,carry##_j ,tt##_j );\
	asm_real_carry_norm(bj##_j ,_d##_j##i ,carry##_j ,tt##_j );\
}
#endif

#ifndef USE_ASM

/* This macro is for last carrie phase */

#define cplx_carry(_d,_j)\
{\
  BIG_DOUBLE temp0;\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = _d##_j##r * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = carry##_j * aux[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[issmall] += tt##_j[issmall];\
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  if(carry##_j != 0.0) \
    {\
	printf(" Carries too high final pass \n"); \
	exit(1);\
    }\
}

#else

# define cplx_carry(_d,_j)\
{\
  asm_real_carry(bj##_j ,_d##_j##r ,carry##_j ,tt##_j);\
  _d##_j##i = carry##_j;\
  asm_real_carry(bj##_j ,_d##_j##i ,carry##_j ,tt##_j);\
  if(carry##_j != 0.0) \
    {\
	printf(" Carries too high final pass \n"); \
	exit(1);\
    }\
}

#endif



/* This macro only has diferent with respect cplx_carry_norm because it
   take into account if l==lastloop */
#ifndef USE_ASM
# define cplx_carry_norm_last_check(_d,_j)\
{\
  BIG_DOUBLE temp0,tempErr,err;\
  temp0 = _d##_j##r;\
  tempErr = RINT( temp0 * tt##_j[1]); \
  err = fabs( temp0 * tt##_j[1] - tempErr);\
  check_err(maxerr,err);\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
  \
  temp0 = _d##_j##i;\
  tempErr = RINT( temp0 * tt##_j[1]); \
  err = fabs( temp0 * tt##_j[1] - tempErr);\
  check_err(maxerr,err);\
  if( l == (pad-1))  bj##_j = 0;\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
}

#else

# define cplx_carry_norm_last_check(_d,_j)\
	asm_real_carry_norm_check(bj##_j,_d##_j##r,carry##_j,tt##_j);\
        bj##_j &= nolast;\
	asm_real_carry_norm_check(bj##_j,_d##_j##i,carry##_j,tt##_j);
#endif

#ifndef USE_ASM

# define cplx_carry_norm_last(_d,_j)\
{\
  BIG_DOUBLE temp0,tempErr;\
  tempErr = RINT( _d##_j##r * tt##_j[1] ); \
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##r = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
  \
  tempErr = RINT( _d##_j##i * tt##_j[1] ); \
  if( l == (pad-1))  bj##_j = 0;\
  issmall=(bj##_j - c) >> (BITS_PER_UL - 1);\
  temp0 = (tempErr + carry##_j) * aux[issmall];\
  bj##_j += ib[issmall];\
  tt##_j[issmall] += tt##_j[issmall];\
  carry##_j = RINT(temp0);\
  tt##_j[1] *= aux[2];\
  _d##_j##i = (temp0 - carry##_j ) * tt##_j[0];\
  tt##_j[0] *= aux[3];\
}

#else

# define cplx_carry_norm_last(_d,_j)\
	asm_real_carry_norm(bj##_j,_d##_j##r,carry##_j,tt##_j);\
        bj##_j &= nolast;\
	asm_real_carry_norm(bj##_j,_d##_j##i,carry##_j,tt##_j);
#endif

# define cplx_carry_norm_check_vector(_d,_j0,_j1)\
  cplx_carry_norm_check(_d,_j0);\
  cplx_carry_norm_check(_d,_j1);

# define cplx_carry_norm_vector(_d,_j0,_j1)\
  cplx_carry_norm(_d,_j0);\
  cplx_carry_norm(_d,_j1);


/*$Id$*/
















