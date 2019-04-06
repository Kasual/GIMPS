# $Id$

bin_PROGRAMS = Glucas

# set the include path found by configure
INCLUDES = -I$(top_srcdir)/. -I$(top_srcdir)/src $(all_includes)

# the library search path.

.c.s:
	$(CC) -S -DHAVE_CONFIG_H -I. -I.. $(CFLAGS) $(CPPFLAGS) -o $@ $<

tricky.o: $(top_srcdir)/src/tricky.c
	$(CC) -c $(CFLAGS0) $(top_srcdir)/src/tricky.c 


LDADD = tricky.o  

bin_SCRIPTS = sumGlucas cmpGlucas Glucas_test

EXTRA_DIST = tricky.c Makefile.Glucas ${bin_SCRIPTS}
 
CLEANFILES = g1* g2* g3* g4* g5* g6* g7* g8* g9* s1* s2* s3* s4* s5* s6* s7* s8* s9* t1* t2* t3* t4* t5* t6* t7* t8* t9*
Glucas_SOURCES = compat.c difdit_16.c difdit_2.c difdit_4.c difdit_8.c \
		difditmm_4.c difditmm_8.c dyadic.c dyadic32.c gio.c glucas.c gsetup.c gstd.c \
		gthreads.c iteration.c lucasini.c radix_10.c radix_12.c radix_14.c radix_16.c \
		radix_2.c radix_32.c radix_4.c radix_5.c radix_6.c radix_7.c radix_8.c radix_9.c \
		radixmm_16.c radixmm_4.c radixmm_5.c radixmm_6.c radixmm_7.c radixmm_8.c radixmm_9.c \
		radixmp_16.c radixmp_32.c radixmp_4.c radixmp_8.c selftest.c version.c yeafft.c \
		yeafft1.c yeainit.c yealucas.c yealucmp.c ynorm.c ynorm_10.c ynorm_12.c \
		ynorm_14.c ynorm_16.c ynorm_4.c ynorm_5.c ynorm_6.c ynorm_7.c ynorm_8.c ynorm_9.c \
		ynormm_5.c ynormm_6.c ynormm_7.c ynormm_8.c ynormm_9.c

if BUILD_YEAFFT
lib_LIBRARIES = libyeafft.a
libyeafft_a_SOURCES = difdit_16.c difdit_2.c difdit_4.c difdit_8.c difditmm_4.c \
	difditmm_8.c dyadic32.c dyadic.c radix_10.c radix_12.c radix_14.c radix_16.c \
	radix_2.c radix_32.c radix_4.c radix_5.c radix_6.c radix_7.c radix_8.c radix_9.c \
	radixmm_16.c radixmm_4.c radixmm_5.c radixmm_6.c radixmm_7.c radixmm_8.c radixmm_9.c \
	radixmp_16.c radixmp_32.c radixmp_4.c radixmp_8.c tricky.c yeafft1.c yeafft.c \
	yeainit.c yealucas.c yealucmp.c ynorm_10.c ynorm_12.c ynorm_14.c ynorm_16.c \
	ynorm_4.c ynorm_5.c ynorm_6.c ynorm_7.c ynorm_8.c ynorm_9.c ynorm.c ynormm_5.c \
	ynormm_6.c ynormm_7.c ynormm_8.c ynormm_9.c
endif

noinst_HEADERS = fft3.h fft3ia64.h fft3sse2.h fft3x86.h fft5.h fft5ia64.h \
	fft5sse2.h fft5x86.h fft7.h fft7ia64.h fft7sse2.h fft7x86.h generic.h glucas.h \
	gmp_yea.h gsetup.h gthreads.h macos-codewarrior-prefix.h mccomp.h mccomplx.h \
	prefetch.h round.h version.h ycx86.h ydebug.h yeafft.h ygeneric4.h ygeneric8.h \
	ygeneric.h ygensse2.h ygenvect.h yia64.h yminimum.h ynorm.h ynormia64.h \
	ynormsse2.h ynormx86.h ynormx86p3.h ysse2.h yx86ath.h yx86comm.h yx86.h yx86p3.h \
	yx86p4.h
