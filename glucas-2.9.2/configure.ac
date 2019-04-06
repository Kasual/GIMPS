dnl $Id$
dnl
dnl Process this file with autoconf to produce a configure script.
dnl (c) 2000-2013
dnl Guillermo Ballester valor gbv@oxixares.com
dnl Klaus Kastens kiste@bawue.com
dnl

AC_INIT(Glucas, 2.9.2)

AC_REVISION($Rev$)

AC_CONFIG_SRCDIR([src/glucas.c])
AM_CONFIG_HEADER([config.h])

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE

AC_LANG_C
AC_PROG_CC
dnl AM_PROG_LIBTOOL


dnl ecc is Intel compiler for ia64
dnl ccc is Compaq-C compiler  for linux-alpha
dnl gcc is gnu compiler
dnl cc is default c-compiler on most system
AC_PROG_CC([ecc ccc gcc cc])
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_RANLIB

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_CHECK_LIB(m, floor)

dnl Checks for header files.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([limits.h unistd.h sys/resource.h sys/time.h sys/utsname.h arpa/inet.h sys/socket.h netinet/in.h netdb.h])

AC_CHECK_FUNCS([strchr getopt atoi strlen gethrtime getrusage gettimeofday chdir floor pow gethostbyname gethrtime gettimeofday memset modf socket strchr strstr uname])

AC_CHECK_DECLS(__builtin_prefetch)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
dnl AC_FUNC_MALLOC 
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE

dnl Checks for library functions.
AC_TYPE_SIGNAL
AC_CHECK_SIZEOF(unsigned long int)
AC_CHECK_SIZEOF(double)
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_STRFTIME


dnl
dnl This block finds 
dnl Default values for CFLAGS CTHREAD and CSTATIC
dnl 
{ 
CFLAGS1=""
case "${CC}" in
	gcc*)
		CFLAGS0=""
		CFLAGS="-O3"
		CTHREAD=""
		CSTATIC="-static"
		;;
	icc)
		CTHREAD=""
		CSTATIC="-static"

dnl
dnl Is actually an IA64 system 
dnl
echo -n "checking whether is actualy an IA64... " >&6
AC_TRY_COMPILE(,{
#if defined(__ia64) 
}
#endif
,IS_IA64="yes"; echo yes, IS_IA64="no"; echo no)

dnl
dnl Intel CC compiler >= 8.0 deprecates ecc, so sometimes 'icc' is old 'ecc' 
dnl
echo -n "checking whether Intel CC version >= 8 and IA64... " >&6
AC_TRY_COMPILE(,{
#if __INTEL_COMPILER >= 800 && defined(__ia64) 
}
#endif
,CFLAGS="-O2"; CFLAGS0="-O0"; echo yes, 
CFLAGS="-O3 -rcd -align"; CFLAGS0="-O0 -rcd -align"; echo no)
		;;

	ecc)
		CFLAGS0="-O0"
		CFLAGS="-O2"
		CTHREAD=""
		CSTATIC="-static"
		;;
	ccc)
		CFLAGS0=""
		CFLAGS="-O5"
		CTHREAD=""
		CSTATIC="-non_shared"
		;;
	*)
		CFLAGS0="-O0"
		CFLAGS="-O3"
		CTHREAD=""
		CSTATIC="-static"
		;;
esac
}


dnl
dnl test to see whether GCC is v.3 or higher
dnl If GCC >= 3.0 then the old problem about 
dnl stack alignement is gone and -mpreferred-stack-bounday
dnl supressed. It is defined -DNO_HACK_ALIGN  
dnl
{
echo -n "checking whether GCC version >= 3 ... " >&6
AC_TRY_COMPILE(,{
#if __GNUC__ >= 3
}
#endif
,echo yes, [echo no; if test ! "$CC" = "icc"; then PREFERRED_STACK="2";CFLAGS="$CFLAGS -malign-double"];fi)
}

dnl
dnl test to see whether GCC is v.3.1 or higher
dnl If GCC >= 3.1 then the some default options for pentium3 
dnl are different
{
echo -n "checking whether GCC version >= 3.1 ... " >&6
AC_TRY_COMPILE(,{
#if (( __GNUC__ == 3) && ( __GNUC_MINOR__ > 0 )) || ( __GNUC__ > 3)
}
#endif
,echo yes; GCC_3_1="yes" , echo no )
}

dnl
dnl test to see whether GCC is v.3.3 or higher
dnl If GCC >= 3.3 then the some options default options for pentium3 
dnl are different
{
echo -n "checking whether GCC version >= 3.3 ... " >&6
AC_TRY_COMPILE(,{
#if (( __GNUC__ == 3) && ( __GNUC_MINOR__ >= 3 )) || ( __GNUC__ > 3)
}
#endif
,echo yes; GCC_3_3="yes" , echo no )
}

dnl
dnl test to see whether GCC is v.3.4 or higher
dnl If GCC >= 3.4 Then some target options can be different
dnl are different
{
echo -n "checking whether GCC version >= 3.4 ... " >&6
AC_TRY_COMPILE(,{
#if (( __GNUC__ == 3) && ( __GNUC_MINOR__ >= 4 )) || ( __GNUC__ > 3)
}
#endif
,echo yes; GCC_3_4="yes" , echo no )
}

dnl
dnl test to see whether GCC is v.4.0 or higher
dnl If GCC >= 4.0 Then some target options can be different
dnl are different
{
echo -n "checking whether GCC version >= 4.0 ... " >&6
AC_TRY_COMPILE(,{
#if (( __GNUC__ == 4) && ( __GNUC_MINOR__ >= 0 )) || ( __GNUC__ > 4)
}
#endif
,echo yes; GCC_4_0="yes" , echo no )
}


dnl
dnl Block to see the compiler options 
dnl depending on the target
dnl 
{
echo "configure:$LINENO: Selecting $CC compiler flags for $target" >&5
echo "Selecting $CC compiler flags for $target ... " >&6

dnl
dnl Block to supress -mpreferred-stack-boundary if needed
dnl and set compiler options  on x86 family
dnl 
if test -z "$GLUCAS_FLAGS"; then GLUCAS_FLAGS=""; fi

case "${target}" in
	i386* | i486* | i586* | i686* | k6* | pentium* | athlon* | x86_64* | core2* )
	if test ! "$CC" = "icc" ; then
	   if test -n "$PREFERRED_STACK" ; then
		OLDFLAGS="$CFLAGS"
		echo -n "checking whether the -mpreferred-stack-boundary works ... "
		CFLAGS="$CFLAGS -fno-defer-pop -mpreferred-stack-boundary=$PREFERRED_STACK"
		AC_TRY_COMPILE(,/* Test program */,echo yes,echo no; CFLAGS="$OLDFLAGS")
		if test "$PREFERRED_STACK" != "2"; then
		   CFLAGS="$CFLAGS -DNO_HACK_ALIGN"
		fi;
	   else
		GLUCAS_FLAGS="$GLUCAS_FLAGS -DNO_HACK_ALIGN"
	   fi;

dnl default options
dnl
dnl enable for profiler
dnl
	AC_ARG_ENABLE(profiler,
		[  --enable-profiler             enable profiler use, default no],
		[if test -n "$enable_profiler"; then 
		if test "$enable_profiler" != "no" ; then 
		   CFLAGS="$CFLAGS -Wall -g -pg -malign-double "
		else
		   CFLAGS="$CFLAGS -Wall -fomit-frame-pointer -malign-double"
		fi
	fi] )


dnl target specific options	
	case "${target}" in
        i386*) 
           OLDFLAGS="$CFLAGS"
	   CPU="i386"
	   DEF_CPU="i386"
	   ;;

	i486*)
           OLDFLAGS="$CFLAGS"
	   CPU="i486"
	   DEF_CPU="i486"
	   ;;

	i586*)
           OLDFLAGS="$CFLAGS"
	   CPU="i586"
	   DEF_CPU="i586"
	   ;;

	pentiumm* )
           OLDFLAGS="$CFLAGS"
	   CPU="pentium-mmx"
	   DEF_CPU="i586"
   	   ;;

	pentiumpro* )
           OLDFLAGS="$CFLAGS"
	   CPU="i686"
	   DEF_CPU="i686"
  	   ;;

	pentium2* )
           OLDFLAGS="$CFLAGS"
	   CPU="pentium2"
	   DEF_CPU="i686"
	   ;;
 
	pentium3* )
           OLDFLAGS="$CFLAGS"
	   CPU="pentium3"
	   DEF_CPU="i686"
	   ;;

	pentium4* )
           OLDFLAGS="$CFLAGS"
	   CPU="pentium4"
	   DEF_CPU="i686"
	   ;;

	athlon* )
           OLDFLAGS="$CFLAGS"
	   CPU="athlon"
	   DEF_CPU="i686"
	   ;;

	x86_64* )
           OLDFLAGS="$CFLAGS"
	   CPU="k8"
	   DEF_CPU="k8"
	   ;;

	core2* )
           OLDFLAGS="$CFLAGS"
	   CPU="pentium4"
	   DEF_CPU="i686"
	   ;;


	k62* )
           OLDFLAGS="$CFLAGS"
	   CPU="k6-2"
	   DEF_CPU="k6"
	   ;;

	k63* )
           OLDFLAGS="$CFLAGS"
	   CPU="k6-3"
	   DEF_CPU="k6"
	   ;;

	k6* )
	   CFLAGS="$CFLAGS -mcpu=k6"
	   ;;

	esac
dnl
dnl	now the test about cpu flags
dnl
	if test "$GCC_3_4" = "yes"; then
	   CPUFLAG="-mcpu=$CPU -mtune=$CPU"
	else
	   CPUFLAG="-mcpu=$CPU -march=$CPU"
	fi	     
        echo -n "checking whether the $CPUFLAG works ... "
	CFLAGS="$CFLAGS $CPUFLAG"
	AC_TRY_COMPILE(,/* Test program */,echo yes,echo no; CFLAGS="$OLDFLAGS -mcpu=i686");


      else
dnl
dnl	this is for intel icc/ecc compiler
dnl

dnl
dnl enable for profiler
dnl
	AC_ARG_ENABLE(profiler,
		[  --enable-profiler             enable profiler use, default no],
		[if test -n "$enable_profiler"; then 
		if test "$enable_profiler" != "no" ; then 
		   CFLAGS="$CFLAGS -g -fp -pg"
		fi
	fi] )

	
     fi   	
     ;; 

dnl
dnl	case for AMD64 
dnl
	x86_64*)
	
	if test ! "$CC" = "icc" ; then
	
		AC_ARG_ENABLE(profiler,
			[  --enable-profiler             enable profiler use, default no],
			[if test -n "$enable_profiler"; then 
			if test "$enable_profiler" != "no" ; then 
			  CFLAGS="$CFLAGS -Wall -g -pg "
			else
			  CFLAGS="$CFLAGS -Wall -fomit-frame-pointer "
			fi
			fi] )

	else

		AC_ARG_ENABLE(profiler,
			[  --enable-profiler             enable profiler use, default no],
			[if test -n "$enable_profiler"; then 
			if test "$enable_profiler" != "no" ; then 
			  CFLAGS="$CFLAGS -g -fp -pg"
			fi
			fi] )
	fi
	;;

esac


dnl
dnl checking whether Openmp flag -omp or -openmp works
dnl
    if test "$CC" != "gcc" ; then
	OLDFLAGS="$CFLAGS"
        echo -n "checking whether the OpenMP -omp, -openmp, +Oopenmp or -xopenmp works ... "
	CFLAGS="$CFLAGS -omp"
	AC_TRY_COMPILE(,/* Test program */,echo yes; OMP="YES" ; omp="-omp" ,echo no) 
	CFLAGS="$OLDFLAGS"
	if test "$OMP" != "YES"; then 
 	  OLDFLAGS="$CFLAGS"
  	  CFLAGS="$CFLAGS -openmp"
 	  AC_TRY_COMPILE(,/* Test program */,echo yes; OMP="YES" ; omp="-openmp" , echo no) 
 	  CFLAGS="$OLDFLAGS"
	fi
	if test "$OMP" != "YES"; then 
 	  OLDFLAGS="$CFLAGS"
  	  CFLAGS="$CFLAGS +Oopenmp"
 	  AC_TRY_COMPILE(,/* Test program */,echo yes; OMP="YES" ; omp="+Oopenmp" , echo no) 
 	  CFLAGS="$OLDFLAGS"
	fi
	if test "$OMP" != "YES"; then 
 	  OLDFLAGS="$CFLAGS"
  	  CFLAGS="$CFLAGS -xopenmp"
 	  AC_TRY_COMPILE(,/* Test program */,echo yes; OMP="YES" ; omp="-xopenmp" , echo no) 
 	  CFLAGS="$OLDFLAGS"
	fi
    fi
        


dnl
dnl enable for ytarget
dnl
AC_ARG_ENABLE(ytarget,
	[  --enable-ytarget=n            configure to Y_TARGET=n, default 0 
                                see info files],
	[if test -n "$enable_ytarget" -a "$enable_ytarget" != "no" ; then 
		Y_TARGET="$enable_ytarget"; fi
	if test "$enable_ytarget" = "yes" ; then Y_TARGET=0;  fi] )

dnl
dnl enable for yaval
dnl
AC_ARG_ENABLE(yaval,
	[  --enable-yaval=n              configure to Y_AVAL=n, default 3], 
	[if test -n "$enable_yaval" -a "$enable_yaval" != "no" ; then 
		Y_AVAL="$enable_yaval"; fi
	if test "$enable_yaval" = "yes" ; then Y_AVAL=3;  fi] )

dnl
dnl enable for ymem_threshold
dnl
AC_ARG_ENABLE(ymem_threshold,
	[  --enable-ymem_threshold=n     Set Y_MEM_THRESHOLD=n, default 2048],
	[if test -n "$enable_ymem_threshold" -a "$enable_ymem_threshold" != "no" ; then 
		Y_MEM_THRESHOLD="$enable_ymem_threshold"; fi
	if test "$enable_ymem_threshold" = "yes" ; then Y_MEM_THRESHOLD=2048 ;  fi])

dnl
dnl enable for yblocksize
dnl
AC_ARG_ENABLE(yblocksize,
	[  --enable-yblocksize=n         Set Y_BLOCKSIZE=n, default 512],
	[if test -n "$enable_yblocksize" -a "$enable_yblocksize" != "no" ; then 
		Y_BLOCKSIZE="$enable_yblocksize"; fi
	if test "$enable_yblocksize" = "yes" ; then Y_BLOCKSIZE=512 ;  fi])

dnl
dnl enable for yshift
dnl
AC_ARG_ENABLE(yshift,
	[  --enable-yshift=n             Set Y_SHIFT=n, default 8],
	[if test -n "$enable_yshift" -a "$enable_yshift" != "no" ; then 
		Y_SHIFT="$enable_yshift"; fi
	if test "$enable_yshift" = "yes" ; then Y_SHIFT=8 ;  fi] )


dnl
dnl enable for ypadding_level
dnl
AC_ARG_ENABLE(ypadding_level,
	[  --enable-ypadding_level=n     Set Y_PADDING_LEVEL=n, default 1],
	[if test -n "$enable_ypadding_level" -a "$enable_ypadding_level" != "no" ; then 
		Y_PADDING_LEVEL="$enable_ypadding_level"; fi
	if test "$enable_ypadding_level" = "yes" ; then Y_PADDING_LEVEL=1 ;  fi] )

dnl
dnl enable for ycache_line
dnl
AC_ARG_ENABLE(ycache_line,
	[  --enable-ycache_line=n        Set Y_CACHE_LINE=n, default 4],
	[if test -n "$enable_ycache_line" -a "$enable_ycache_line" != "no" ; then 
		Y_CACHE_LINE="$enable_ycache_line"; fi
	if test "$enable_yshift" = "yes" ; then Y_CACHE_LINE=4 ;  fi] )

dnl
dnl enable for ysecure
dnl
AC_ARG_ENABLE(ysecure,
	[  --enable-ysecure              define Y_SECURE, default no],
	[if test -n "$enable_ysecure"; then 
	   if test "$enable_ysecure" != "no" ; then 
		Y_SECURE=1   
	   else 
                Y_SECURE=0
           fi
	fi] )

dnl
dnl enable for ysse2
dnl
AC_ARG_ENABLE(ysse2,
	[  --enable-ysse2                define Y_USE_SSE2, default no],
	[if test -n "$enable_ysse2"; then 
	   if test "$enable_ysse2" != "no" ; then 
		Y_USE_SSE2=1   
	   else 
                Y_USE_SSE2=0
           fi
	fi] )

dnl
dnl enable for yminimum
dnl
AC_ARG_ENABLE(yminimum,
	[  --enable-yminimum             define Y_MINIMUM, default no],
	[if test -n "$enable_yminimum"; then 
	   if test "$enable_yminimum" != "no" ; then 
		Y_MINIMUM=1   
	   else 
                Y_MINIMUM=0
           fi
	fi] )

dnl
dnl enable for ymaximum
dnl
AC_ARG_ENABLE(ymaximum,
	[  --enable-ymaximum             define Y_MAXIMUM, default no],
	[if test -n "$enable_ymaximum"; then 
	   if test "$enable_ymaximum" != "no" ; then 
		Y_MAXIMUM=1   
	   else 
                Y_MAXIMUM=0
           fi
	fi] )

dnl
dnl enable for ykill_branches
dnl
AC_ARG_ENABLE(ykill_branches,
	[  --enable-ykill_branches       define Y_KILL_BRANCHES, default no],
	[if test -n "$enable_ykill_branches"; then 
	   if test "$enable_ykill_branches" != "no" ; then 
		Y_KILL_BRANCHES=1   
	   else 
                Y_KILL_BRANCHES=0
           fi
	fi] )

dnl
dnl enable for ymany_registers
dnl
AC_ARG_ENABLE(ymany_registers,
	[  --enable-ymany_registers      define Y_MANY_REGISTERS, default no],
	[if test -n "$enable_ymany_registers"; then 
	   if test "$enable_ymany_registers" != "no" ; then 
		Y_MANY_REGISTERS=1   
	   else 
                Y_MANY_REGISTERS=0
           fi
	fi] )

dnl
dnl enable for yvectorize
dnl
AC_ARG_ENABLE(yvectorize,
	[  --enable-yvectorize           define Y_VECTORIZE, default no],
	[if test -n "$enable_yvectorize"; then 
	   if test "$enable_yvectorize" != "no" ; then 
		Y_VECTORIZE=1   
	   else 
                Y_VECTORIZE=0
           fi
	fi] )

dnl
dnl enable for yvectorize2
dnl
AC_ARG_ENABLE(yvectorize2,
	[  --enable-yvectorize2          define Y_VECTORIZE2, default no],
	[if test -n "$enable_yvectorize2"; then 
	   if test "$enable_yvectorize2" != "no" ; then 
		Y_VECTORIZE2=1   
	   else 
                Y_VECTORIZE2=0
           fi
	fi] )

dnl
dnl enable for ylong_macros
dnl
AC_ARG_ENABLE(ylong_macros,
	[  --enable-ylong_macros         define Y_LONG_MACROS, default no],
	[if test -n "$enable_ylong_macros"; then 
	   if test "$enable_ylong_macros" != "no" ; then 
		Y_LONG_MACROS=1   
	   else 
                Y_LONG_MACROS=0
           fi
	fi] )

dnl
dnl enable for ydefines_file
dnl
AC_ARG_ENABLE(ydefines_file,
	[  --enable-ydefines_file        define Y_DEFINES_FILE, default no],
	[if test -n "$enable_ydefines_file"; then 
	   if test "$enable_ydefines_file" != "no" ; then 
		Y_DEFINES_FILE=1   
	   else 
                Y_DEFINES_FILE=0
           fi
	fi] )

dnl
dnl enable for yprefetch_expensive
dnl
AC_ARG_ENABLE(yprefetch_expensive,
	[  --enable-yprefetch_expensive  define Y_PREFETCH_EXPENSIVE, default no],
	[if test -n "$enable_yprefetch_expensive"; then 
	   if test "$enable_yprefetch_expensive" != "no" ; then 
		Y_PREFETCH_EXPENSIVE=1   
	   else 
                Y_PREFETCH_EXPENSIVE=0
           fi
	fi] )

dnl
dnl enable for build_yeafft
dnl
AC_ARG_ENABLE(ybuild_yeafft,
	[  --enable-ybuild_yeafft        builds libyeafft.a in adition to Glucas, default no],
	[if test -n "$enable_ybuild_yeafft"; then 
	   if test "$enable_ybuild_yeafft" != "no" ; then 
		Y_BUILD_YEAFFT=true   
	   else 
                Y_BUILD_YEAFFT=false
           fi
	fi] )

AM_CONDITIONAL(BUILD_YEAFFT, test x$Y_BUILD_YEAFFT = xtrue )

case "${target}" in
	i386*-*-*gnu* | i386*-*-*bsd* )
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	;;

	i486*-*-*gnu* | i486*-*-*bsd* )
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	;;
	
	k6*-*-*gnu* | k6*-*-*bsd*)
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then Y_MINIMUM=1 ; fi 
	;;

	i586*-*-*gnu* | pentiumm*-linux-*gnu* )
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp5"
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	else 
		if test -z "$Y_TARGET"; then Y_TARGET=1 ; fi 
		if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then Y_MINIMUM=1 ; fi 
	;;

	pentiumpro*-linux-gnu* | pentium2*-linux-gnu* )
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp6"
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	else 
		if test -z "$Y_TARGET"; then Y_TARGET=1 ; fi 
	fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then Y_MINIMUM=1 ; fi 
	;;


	i586*-*-*bsd* | pentiumm*-*-*bsd* | pentiumpro*-*-*bsd* | pentium2*-*-*bsd*)
	if test -z "$Y_TARGET"; then Y_TARGET=1 ; fi 
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then Y_MINIMUM=1 ; fi 
	;;

	i686*-*-*gnu* | pentium3*-linux-*gnu* )
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp6 -xK"
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	else
		if test -z "$Y_TARGET"; then
			if test "$GCC_3_1" = "yes" ; then 
			   Y_TARGET=0
			else 
			   Y_TARGET=11 
			fi 
		fi
	fi
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then 
		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
	fi
	;;

	pentium4*-linux-*gnu | pentium4*-pc-cygwin | core2*-*-* )
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp7 -xW "
		CFLAGS0="$CFLAGS0 -tpp7 -xW "
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
		if test -z "$Y_USE_SSE2"; then Y_USE_SSE2=1 ; fi 
		Y_PENTIUM4=1;		
	else
		if test -z "$Y_TARGET"; then
			if test "$GCC_3_1" = "yes" ; then 
			   Y_TARGET=0
			else 
			   Y_TARGET=13 
			fi 
			if test "$GCC_3_3" = "yes" ; then 
			   CFLAGS="$CFLAGS -msse2"
			   CFLAGS0="$CFLAGS0 -msse2"
			   Y_PENTIUM4=1
			   if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
			   if test -z "$Y_USE_SSE2"; then Y_USE_SSE2=1 ; fi 
			else
			   CFLAGS="$CFLAGS -Wall -malign-double -mcpu=i686 -march=i686"
			fi			   
		else
		   CFLAGS="$CFLAGS -Wall -malign-double -mcpu=i686 -march=i686"
		fi
        fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=16 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=1024 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=7 ; fi 
dnl	if test -z "$Y_MINIMUM"; then 
dnl		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
dnl	fi
	;;

	x86_64*-linux-*gnu | amd64*-linux-*gnu | x86_64*-*-*bsd* | amd64*-*-*bsd*)
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp7 -xW"
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
		if test -z "$Y_USE_SSE2"; then Y_USE_SSE2=1 ; fi 
		Y_PENTIUM4=1;		
	else
		if test -z "$Y_TARGET"; then
			if test "$GCC_3_1" = "yes" ; then 
			   Y_TARGET=0
			else 
			   Y_TARGET=13 
			fi 
			if test "$GCC_3_3" = "yes" ; then 
			   CFLAGS="$CFLAGS -msse2"
			   CFLAGS0="$CFLAGS0 -msse2"
			   Y_AMD64=1
			   if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
			   if test -z "$Y_USE_SSE2"; then Y_USE_SSE2=1 ; fi 
			   if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
			else
			   CFLAGS="$CFLAGS -Wall -msse2"
			fi			   
		else
		   CFLAGS="$CFLAGS -Wall -msse2"
		fi
        fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=16 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=128 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=6 ; fi
	if test -z "$Y_PADDING_LEVEL"; then Y_PADDING_LEVEL=4 ; fi
 
dnl	if test -z "$Y_MINIMUM"; then 
dnl		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
dnl	fi
	;;

	i686*-*-*bsd* | pentium3*-*-*bsd* )
	if test -z "$Y_TARGET"; then
		if test "$GCC_3_1" = "yes" ; then 
		   Y_TARGET=0
		else 
		   Y_TARGET=11 
		fi 
	fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then 
		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
	fi
	;;

	pentium4*-*-*bsd* )
	if test -z "$Y_TARGET"; then
		if test "$GCC_3_1" = "yes" ; then 
		   Y_TARGET=0
		else 
		   Y_TARGET=13 
		fi 
	fi
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=16 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then 
		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
	fi
	;;

	athlon*-linux-*gnu*  | athlon*-*-*bsd* )
	if test "$CC" = "icc" ; then
		CFLAGS="$CFLAGS -tpp6"
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	else
		if test -z "$Y_TARGET"; then
			if test "$GCC_3_1" = "yes" ; then 
			   Y_TARGET=0
			else 
			   Y_TARGET=12 
			fi
		fi 
	fi
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=128 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=6 ; fi 
	if test -z "$Y_PADDING_LEVEL"; then Y_PADDING_LEVEL=3 ; fi 
	if test -z "$Y_AVAL"; then Y_AVAL=4 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_SECURE"; then Y_SECURE=1 ; fi 
	if test -z "$Y_MINIMUM"; then 
		if test "$GCC_3_1" != "yes" ; then Y_MINIMUM=1 ; fi 
	fi
	;;

	ia64*-*-*gnu* | ia64*-*-*bsd* )
	if test "$CC" = "ecc" -o "$CC" = "icc" ; then 
	   if test "IS_IA64" = "yes" ; then
	      CFLAGS="-O2"
	   else
	      CFLAGS="-O3 -rcd"
	   fi
	else
	   CFLAGS="-O3"
	fi
	GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_ITANIUM"
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	;;

	ia64*-*-*hpux* )
	if test "$CC" = "cc" ; then 
	   CFLAGS="-O4 -fast"
	else
	   CFLAGS="-O3"
	fi
	GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_ITANIUM"
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	;;

	alphaev56*-dec-*)
	CFLAGS="-O5 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
	CTHREAD="-pthread"
	CSTATIC="-non_shared"
	if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	alphaev6-dec-*)
	CFLAGS="-O5 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
	CTHREAD="-pthread"
	CSTATIC="-non_shared"
	if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;


	alphaev67-dec-* | alphaev68-dec-* | alpha*-dec-* )
	CFLAGS="-O5 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
	CTHREAD="-pthread"
	CSTATIC="-non_shared"
	if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;


	alphaev56*-*-*gnu* )
	if test "$CC" = "ccc" ; then
		CFLAGS="-O4 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
		CSTATIC="-non_shared"
		if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
		if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	else
		CFLAGS="$CFLAGS -mcpu=ev56"
		if test -z "$Y_TARGET"; then Y_TARGET=32 ; fi 
	fi
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	alphaev6-*-*gnu* )
	if test "$CC" = "ccc" ; then
		CFLAGS="-O4 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
		CSTATIC="-non_shared"
		if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	else
		CFLAGS="$CFLAGS -mcpu=ev6"
		if test -z "$Y_TARGET"; then Y_TARGET=32 ; fi 
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	fi
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	alphaev67-*-*gnu* | alphaev68-*-*gnu* )
	if test "$CC" = "ccc" ; then
		CFLAGS="-O4 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
		if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
        	CSTATIC="-non_shared"
	else
		CFLAGS="$CFLAGS -mcpu=ev6"
		if test -z "$Y_TARGET"; then Y_TARGET=32 ; fi 
		if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	fi
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	alpha*-*-*gnu* )
	if test "$CC" = "ccc" ; then
		CFLAGS="-O4 -fast -tune host -arch host -assume accuracy_sensitive -unroll 1"
		CSTATIC="-non_shared"
		if test -z "$Y_TARGET"; then Y_TARGET=31 ; fi 
	else
		CFLAGS="$CFLAGS -mcpu=ev5"
		if test -z "$Y_TARGET"; then Y_TARGET=32 ; fi 
		if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	fi
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=8 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	sparcv9-*-* | ultrasparc*-*-* )
	if test "$CC" = "gcc" ; then 
	  CFLAGS="-O3 -Wall -fomit-frame-pointer -mcpu=ultrasparc"
	else
	  CFLAGS="-fast -Xa"
	  CSTATIC="-Bstatic"
 	  OLDFLAGS="$CFLAGS"
          echo -n "checking whether the SunMP -xexplicitpar works ... "
	  CFLAGS="$CFLAGS -xexplicitpar"
	  AC_TRY_COMPILE(,/* Test program */,echo yes; SUNMP="YES" , echo no) 
	  CFLAGS="$OLDFLAGS"
	fi
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	;;

	powerpc970-*-darwin* )
	if test "$GCC_4_0" = "yes" ; then 
		CFLAGS="-O3 -fomit-frame-pointer -mcpu=powerpc -falign-loops=16 -mtune=G5 -mcpu=970 -ffast-math"
		enable_static=no
	else
		CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -mcpu=powerpc"
	fi
	if test -z "$Y_TARGET"; then Y_TARGET=24 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=6 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	if test -z "$Y_CACHE_LINE"; then Y_CACHE_LINE=16 ; fi 
	;;

	powerpc*-*-darwin* )
	CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -mcpu=powerpc"
	if test -z "$Y_TARGET"; then Y_TARGET=23 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=6 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	;;

	powerpc*-*-*gnu* )
	CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -mcpu=powerpc"
	if test -z "$Y_TARGET"; then Y_TARGET=23 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=6 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	;;

	powerpc*-ibm-aix* )
	if test "$CC" = "gcc" ; then 
	  CFLAGS="-O3 -Wall"
	  if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	else
	  CFLAGS="-O3 -qstrict"
	fi	  
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	;;

	hppa64*-*-*gnu* | hppa2*-*-*gnu* )
	CFLAGS="-O3"
	if test -z "$Y_TARGET"; then Y_TARGET=61 ; fi 
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=32768 ; fi 
	if test -z "$Y_BLOCKSIZE"; then Y_BLOCKSIZE=4096 ; fi 
	if test -z "$Y_SHIFT"; then Y_SHIFT=9 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	if test -z "$Y_LONG_MACROS"; then Y_LONG_MACROS=1 ; fi 
	if test -z "$Y_VECTORIZE2"; then Y_VECTORIZE2=1 ; fi 
	;;

	s390-*-*-*gnu* )
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=4096 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	;;


	alpha*-*-*  | sparc*-*-* | powerpc*-*-* )
	if test -z "$Y_MEM_THRESHOLD"; then Y_MEM_THRESHOLD=8192 ; fi 
	if test -z "$Y_KILL_BRANCHES"; then Y_KILL_BRANCHES=1 ; fi 
	if test -z "$Y_VECTORIZE"; then Y_VECTORIZE=1 ; fi 
	;;


	*-*-*) 
	echo " ****WARNING************************************ "
	echo " Your system has not been checked for Glucas yet."
        echo " It identifies itself as being of type ${target}."
        echo " Use --enable-ytarget=<number> to permit Glucas to"
        echo " attempt target-specific optimisations if a similar"
        echo " known target exists. Otherwise, or in addition,"
        echo " use the simple 'makefile' file and try your compiler"
        echo " and Glucas options (see doc files). With the aid of"
	echo " provided tests, your goal is to get a Glucas program as"
	echo " fast as possible. If you have done so, please send your" 
	echo " complier options to <gbv@oxixares.com> to be included" 
	echo " in future configure releases."
	echo " THANKS FOR YOUR HELP."
	;;
esac 
}

dnl concatenate partial results
if test -n "$Y_DEFINES_FILE" -a "$Y_DEFINES_FILE" = "1" ; then
   GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_DEFINES_FILE"
else 
 if test -n "$Y_TARGET" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_TARGET=$Y_TARGET"
 fi

 if test -n "$Y_AVAL" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_AVAL=$Y_AVAL"
 fi

 if test -n "$Y_MEM_THRESHOLD" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_MEM_THRESHOLD=$Y_MEM_THRESHOLD"
 fi

 if test -n "$Y_BLOCKSIZE" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_BLOCKSIZE=$Y_BLOCKSIZE"
 fi

 if test -n "$Y_SHIFT" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_SHIFT=$Y_SHIFT" 
 fi

 if test -n "$Y_PADDING_LEVEL" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_PADDING_LEVEL=$Y_PADDING_LEVEL"
 fi

 if test -n "$Y_CACHE_LINE" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_CACHE_LINE=$Y_CACHE_LINE"
 fi

 if test -n "$Y_SECURE" -a "$Y_SECURE" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_SECURE"
 fi 

 if test -n "$Y_MINIMUM" -a "$Y_MINIMUM" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_MINIMUM"
 fi

 if test -n "$Y_MAXIMUM" -a "$Y_MAXIMUM" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_MAXIMUM"
 fi
 
 if test -n "$Y_KILL_BRANCHES" -a "$Y_KILL_BRANCHES" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_KILL_BRANCHES"
 fi

 if test -n "$Y_MANY_REGISTERS" -a "$Y_MANY_REGISTERS" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_MANY_REGISTERS"
 fi

 if test -n "$Y_VECTORIZE" -a "$Y_VECTORIZE" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_VECTORIZE"
 fi

 if test -n "$Y_VECTORIZE2" -a "$Y_VECTORIZE2" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_VECTORIZE2"
 fi

 if test -n "$Y_LONG_MACROS" -a "$Y_LONG_MACROS" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_LONG_MACROS"
 fi

 if test -n "$Y_PREFETCH_EXPENSIVE" -a "$Y_PREFETCH_EXPENSIVE" = "1" ; then 
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_PREFETCH_EXPENSIVE"
 fi

 if test -n "$Y_USE_SSE2" -a "$Y_USE_SSE2" = "1" ; then
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_USE_SSE2"
 fi 

 if test -n "$Y_PENTIUM4" -a "$Y_PENTIUM4" = "1" ; then
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_PENTIUM4"
 fi 

 if test -n "$Y_AMD64" -a "$Y_AMD64" = "1" ; then
    GLUCAS_FLAGS="$GLUCAS_FLAGS -DY_AMD64"
 fi 
fi

CFLAGS="$CFLAGS $GLUCAS_FLAGS"
CFLAGS0="$CFLAGS0 $GLUCAS_FLAGS"

dnl
dnl enable for posix thread
dnl
AC_ARG_ENABLE(pthread,
	[  --enable-pthread=n            configure to run Glucas with n Posix-threads ],
	[if test -n "$enable_pthread" -a "$enable_pthread" != "no"; then
		LIBS="$LIBS -lpthread"
		if test "$enable_pthread" = "yes"; then enable_pthread=2; fi
		CFLAGS="$CFLAGS $CTHREAD -D_PTHREADS=$enable_pthread"
		CFLAGS0="$CFLAGS0 -D_PTHREADS=$enable_pthread"
	fi]
) 


dnl
dnl enable for OpenMP directives
dnl
AC_ARG_ENABLE(openmp,
	[  --enable-openmp=n             configure to run Glucas with OpenMP directives ],
	[if test -n "$enable_openmp" -a "$enable_openmp" != "no"; then
	   if test "$OMP" != "YES"; then
		echo "Warning: $CC compiler does not accept OpenMP directives"
	   else 
		if test "$enable_openmp" != "yes" ; then
		    CFLAGS="$CFLAGS -DY_NUM_THREADS=$enable_openmp" 
		    CFLAGS0="$CFLAGS0 -DY_NUM_THREADS=$enable_openmp" 
		fi
		CFLAGS="$omp $CFLAGS"
	   fi
	fi]
) 

dnl
dnl enable for SunMP directives
dnl
AC_ARG_ENABLE(sunmp,
	[  --enable-sunmp=n              configure to run Glucas with SunMP directives ],
	[if test -n "$enable_sunmp" -a "$enable_sunmp" != "no"; then
	   if test "$SUNMP" != "YES"; then
		echo "Warning: $CC compiler does not accept SunMP directives"
	   else 
		CFLAGS="-xexplicitpar $CFLAGS -D_SUNMP"  
		if test "$enable_sunmp" != "yes" ; then
		    CFLAGS="$CFLAGS -DY_NUM_THREADS=$enable_sunmp" 
		else
		    CFLAGS="$CFLAGS -DY_NUM_THREADS=2" 
		    CFLAGS0="$CFLAGS0 -DY_NUM_THREADS=2" 
		fi
	   fi
	fi]
) 

dnl enable for build static linked binaries
AC_ARG_ENABLE(static,
	[  --enable-static               configure to build statically linked binaries ],
	[if test -n "$enable_static" -a "$enable_static" != "no"; then
		CFLAGS="$CSTATIC $CFLAGS"
	fi]
) 



echo "configure:$LINENO: Minimum optimization flags for tricky.c:" >&5  
echo "Minimum optimization flags for tricky.c:" >&6  
echo "configure:$LINENO: CFLAGS0=$CFLAGS0" >&5
echo "CFLAGS0=$CFLAGS0"  >&6
echo "configure:$LINENO: Full optimization flags:" >&5
echo "Full optimization flags:" >&6
echo "configure:$LINENO: CFLAGS=$CFLAGS" >&5
echo "CFLAGS=$CFLAGS" >&6

if test -n "$Y_DEFINES_FILE" -a "$Y_DEFINES_FILE" = "1" ; then

echo "/* This file was generated by configure script */" > src/yeafft_defines.h

if test -n "$Y_TARGET" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_TARGET $Y_TARGET
_ACEOF
fi

if test -n "$Y_AVAL" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_AVAL $Y_AVAL
_ACEOF
fi

if test -n "$Y_MEM_THRESHOLD" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_MEM_THRESHOLD $Y_MEM_THRESHOLD
_ACEOF
fi

if test -n "$Y_BLOCKSIZE" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_BLOCKSIZE $Y_BLOCKSIZE
_ACEOF
fi

if test -n "$Y_SHIFT" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_SHIFT $Y_SHIFT
_ACEOF
fi

if test -n "$Y_PADDING_LEVEL" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_PADDING_LEVEL $Y_PADDING_LEVEL
_ACEOF
fi

if test -n "$Y_CACHE_LINE" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_CACHE_LINE $Y_CACHE_LINE
_ACEOF
fi

if test -n "$Y_SECURE" -a "$Y_SECURE" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_SECURE
_ACEOF
fi

if test -n "$Y_MINIMUM" -a "$Y_MINIMUM" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_MINIMUM
_ACEOF
fi

if test -n "$Y_MAXIMUM" -a "$Y_MAXIMUM" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_MAXIMUM
_ACEOF
fi

if test -n "$Y_KILL_BRANCHES" -a "$Y_KILL_BRANCHES" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_KILL_BRANCHES
_ACEOF
fi

if test -n "$Y_MANY_REGISTERS" -a "$Y_MANY_REGISTERS" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_MANY_REGISTERS
_ACEOF
fi

if test -n "$Y_VECTORIZE" -a "$Y_VECTORIZE" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_VECTORIZE
_ACEOF
fi

if test -n "$Y_VECTORIZE2" -a "$Y_VECTORIZE2" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_VECTORIZE2
_ACEOF
fi

if test -n "$Y_LONG_MACROS" -a "$Y_LONG_MACROS" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_LONG_MACROS
_ACEOF
fi

if test -n "$Y_PREFETCH_EXPENSIVE" -a "$Y_PREFETCH_EXPENSIVE" = "1" ;  
then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_PREFETCH_EXPENSIVE
_ACEOF
fi

if test -n "$Y_USE_SSE2" -a "$Y_USE_SSE2" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_USE_SSE2
_ACEOF
fi

if test -n "$Y_PENTIUM4" -a "$Y_PENTIUM4" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_PENTIUM4
_ACEOF
fi

if test -n "$Y_AMD64" -a "$Y_AMD64" = "1" ; then
cat >> src/yeafft_defines.h <<_ACEOF
#define Y_AMD64
_ACEOF
fi
echo "configure:$LINENO:" >&5
echo " "
cat src/yeafft_defines.h  >&5
cat src/yeafft_defines.h >&6

fi

AC_SUBST(CFLAGS0)

AC_OUTPUT(Makefile src/Makefile doc/Makefile)

