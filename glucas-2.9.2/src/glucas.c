/* glucas.c
   
   $Id$
   (c) Guillermo Ballester Valor   gbv@oxixares.com
   (c) Klaus Kastens               kiste@bawue.de
 
   2000-2006
 
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

Revision history:
   v0.0 15/apr/2000 
	The starting point. It only uses the convolution provided by
   my convolver package YEAFFT (Yet Another Fast Fourier Transform). 
   It was the first YEAFFT version. It makes an in-place forward Decimation 
   in Frequency transform (DIF), a dyadic multiplication (using the 
   nested-complex representation of a real float array), and a backward 
   decimation in time (DIT). Only tested in my Pentium-166 MMX. Timings
   10% higher than FFTW. Not bad at all for a first version. Not released.
 
   v.0.1 20/may/2000
	Uses the new dyadic mul. features of YEAFFT package. to avoid 
   memory accesses when possible, YEAFFT makes the last forward DIF pass,
   the dyadic mul and first backward DIT pass in a single pass. 10% average
   gain of performance. Timings similar to FFTW on small exponents (in a 
   Pentium). For big FFT lengths YEAFFT is better. Now, tested in my Pentium 
   and in Ernst Mayer MIPS and Alpha machines. It seems really good. About 
   only 10% slower than E.Mayer Mlucas fortran-90 code. Not released.
 
   v.1.0 05/jun/2000
	A lot of small improvements. Over all, Glucas now makes the last
   backward DIT pass, carry and normalization, and first DIF pass using only 
   a single memory access. 10%-15% performance improvement. Results similar
   to Mlucas. Now the program adjusts its accuracy dynamically.    
 
   v.1.1 10/jun/2000
        A hidden bug in y_norm*.c fixed. This bug only worked when testing 
   very big exponents and randomly. It was in the init bjs on carry_norm 
   process. There was a lack of precision for 32 bits unsigned integers. 
 
   V.1.2 12/jun/2000
        I introduced a bug in radix_12.c a bad day trying to write something
   better. Fortunately my script test has detected it and I fixed. The main
   header file yeafft.h has changed to let define some important defs in
   compile time using the -Dmacro=def facility of most C compilers. I wrote
   a little README.file
 
   v.1.3 18/jun/2000
        New radix-9 introduced. If Y_AVAL > 3, YEAFFT can now use radix-9 
   reduction and then gain 10% performance in some exponents.
 
   v.1.4 12/jul/2000
        Best non_power_of_two radix-reduction. Now Glucas uses Nussbaumer 
   method with some improvements from E.W. MAYER. 5% faster.
        Improved check-results-security. Now Glucas makes a roundoff error
   check every iteration. During first init_iterations Glucas assures the 
   roundoff error is less than init_kErrLimit, and dynamically adjust the
   precision to achieve that. In the following iterations the roundoff 
   error must be less than kErrLimit. See the code to the actual defined
   values.
        A semi_bug fixed. Now Glucas saves the update rate of weight factors 
   during the adjust_of_accuracy phase. When restarting, Glucas reads the 
   optimal values found.
        New features: 
	-A verbose mode can be activated via glucas.ini file. See
   README.Glucas file to the details. Glucas can output, every established
   iterations, the last roundoff error, the user and system time, and the
   real time between two outputs. 
        -For test/check proposes, an only_check mode can be activated from
   glucas.ini file. See again README.Glucas file.
 
   v.1.10 18/aug/2000:
	This is the first stand alone version for glucas. No more use of MERS
   routines, but still there is a lot of code lines from it, in special the
   setup routines. Not affected the core and/or speed. Deep changes on usage
   and input/output features:
	-Now Glucas does the work queued in a queue file we can select. 
   Using some command arguments  options, we can select how Glucas manage
   the queue. There are several forms and input file formats we can add work
   to  queue (see the documentation files to information).
	-We can select, using command input arguments, the 'ini' file, 
   the 'queue' file, and the 'results' file.
	-We can enter a no_nice option. (see the README.Glucas.html, it can 
   slow down other jobs, be careful).
	Other important improvements:
	-The response time to a break signal has been reduced. In earlier 
   versions we had to wait about two iterations, now we wait less than
   a iteration.
	-There is now two save files formats we can choose. The default
   stores the residue of a iteration in integer mode but we can select the
   old style (compatible with MERS save format) which saves the residue 
   using the internal floats. The integer format is only a bit slower but 
   the save files are more than three times shorter, there is almost a 
   universal compatibility, and Glucas can now resize the FFT length at run
   time taking the results obtained with other FFT run length.
 
 
   v.2.0 01/dec/2000
        -New Interchangeable Mersenne Residue Format is introduced.
	-A lot of small other improvements
 
   v.2.1 06/dec/2000
        -Now Glucas can be compiled by Codewarrior on iMac
 
   v.2.2 16/dec/2000
        -A hidden bug affecting little exponents fixed
	-The info also displays the sec/iter.
 
   v.2.3 20/jan/2001 (not released)
        -Fixed a bug in displaying the sec/iter.
 
   v.2.4 28/jan/2001 (not released)
        -A bug increasing the memory requirements every time a save file 
	is written fixed.
   
   v.2.5 06/feb/2001 
        -Added the command line option -W to work with other than actual directory.
 
   v.2.6 28/feb/2001
        This version has included:
        -The possibility to set the roundoff error check in the 'inifile'
	Roundoff_check=1. Glucas will make this check all the iterations
	Roundoff_check=0. The roundoff_check is disabled.
 
	-The threshold is less conservative. By default the roundoff check works as
	follows:
	1)Is activated the first 2^17=131092 iterations.
	2)Is activated the first 2^10=1024 iterations after a restart.
	3)If no 1) and 2) then there is a roundoff check every 64 iterations.
	The threshold in 1) is 0.40. The threshold in 2 and 3 is 0.45.
	
	Fixed bugs:
	-A bug found by Brian J. Beesley when is defined Y_MANY_REGISTERS. Now Glucas
	does not search the unwritten routine 'ynorm_18'
    
	-A semi-bug in selftest 25 and 26. The check exponents are changed and so Glucas 
	does not jump to the next FFT-length.
 
    v.2.7 12/mar/2001
 
	Some QA facilities has been added, and some new features:
 
	-The default format for i/o save files is the 'interchangeable format'. It will be so
    for the next versions.
 
	-The default file 'glucas.ini' will be created whether there is no _ini_ file in 
    the working directory.
 
	-When Glucas start, if there is no 'Machine_id=n' line in the ini file, it will 
    append a Machine Identifier line.
 
	-The g<exponent> file now will have two more lines, the Identifier machine and the
    first iteration with roundcheck done in that machine. If one use a save file from other
    machine, then there will be not a g<> file and then, as prevention, it will make roundoff
    check during 131092 next iterations. In addition, if the Machine_id and the identifier 
    read in g<exponent> file  are not the same then is supposed a 'Glucas_environment_change':
        1) g<exponent> file will be updated with the actual iteration and the Machine_id read
    from _ini_ file
	2) A roundoff check error will be made next 131092 iterations.
 
	-When Glucas writes a save file, if Verbose_flag is activated it also display the Res64
    of the iteration being saved. 
 
	-For QA proposes, it is possible to write Interim files every n iterations including 
    the line
 
	QA_interim_file=n
	
    in the _ini_ file. The Res64 will be also displayed.
 
	
	Fixed bugs:
       
	-It has been reported coredumps when there is no _ini_ file. Now there is no problem
    with that.
 
	-When the working directory was other than current and the path was too large, a core
    dump was probable. It is fixed now.
 
        -The max number of exponent in a queue file was 50. Now is 128, I think it is enough.
 
    v.2.7b 25/mar/2001
    
	-The display iterations info now prints the percent of work completed
    instead of the final iteration.	
 
	Fixed bugs:
 
	-For some C-compilers, open the file with "r" is not the same than "rb". In GNU-c
    compiler the "b" does not has any effect. Related with this, Klaus Kastens found a problem
    with Mac Clients. In compatibility mode, save files from Mac and others platforms are not
    compatible. Now there is fixed.
 
        -In some circumstances, the selftest were not done properly. Now is fixed.
 
    V.2.8a  14/jun/2001
	
	-The arrays with factors for DWT now has the same padding than the
	main array. It is to avoid some thrashing memory.
 
	-Some tricks has been included to avoid to compute the true address
	from logical address so often. Hope there will be a gain of about 
	3-5%. 
	
	-For some systems, there is save of time computing almost all Trig.
	factors when needed instead of read from memory. This is activated
	defining Y_MINIMUM. Some new code make this option faster, but is 
	less accurate.
		
	-For x86 machines some new assembler macros has been written.
 
	-For Pentium3 and Athlon, there is some hints to prefetch data. There
	is a gain about 15% of performance. B.J. Beesley has found a big 
	improvement for Athlons tunning properly the prefetch hints.
 
	-Selftest now has only activated Roundoff check error during first
	50 iterations. So one can see the effects of disabling this check.
 
	v.2.8b  07/Aug/2001
 
	-Great Prefetch working progress has made for Alphas. B.J.Beesley 
	found the way to insert assembler prefetch hints in Compaq-c code. The 
	improvement is about 30% or more for ev6 and ev67.!
 
	-Other prefetch hints has been coded for other platforms. At the moment
	there is no significant gain for other than x86 and powerpc. 
 
	-Good news for Mac OS users, both with classic MacOS and Mac OS X:
	Big performance improvement for powerpc family (10%, about 3% adding
	prefetch hints and 7% tuning the parameters). Klaus kastens and 
	Tom Cage did the job. 
	
	-Binaries for Itanium IA-64 has improved a lot, but now the credits
	are for GNU/gcc team. With gcc 3.0 now Glucas is almost twice faster.
 
	-Long macros has been coded for radices 4 and 8. It could take more
	advantages of prefetch and help to less clever compilers. It can be 
	activated with -DY_LONG_MACROS compiler flag. The gain is from 5% to
	-1% .
 
	-Some routines has been recoded to hide some dependency stalls and to
	make easy to vectorize with instructions like altivec G4+ or SSE2.
	We can activate it -DY_VECTORIZE. We have to set -DY_KILL_BRANCHES to
	do any effect. About 1% gain in most cases.
 
	-Radix 4 can be recoded to use other vectorized macros, doing two
	radix-4 transform in a single loop pass. It can be useful sometimes.
	it is equivalent to unroll those loops. To activate -DY_VECTORIZE2
 
	-When the prefetch hints have a lot of cost, we can try the flag
	-DY_VECTORIZE_EXPENSIVE. Radices 4 and 8 will be unrolled to save a
	lot of prefetch calls. It is still an experimental feature.
	
	-Selftest now outputs all the Glucas flags active. It is useful
	in tuning and developing tasks
 
	-If option Alternative_output_flag == 2, now the output is driven
	both to stdout (console) and file set with option Output_file.
 
	-Compiler time has been reduced a lot in most cases. The routines 
	are trivial when we no need them. It will make the developer work 
	easy.
 
	v.2.8b1  08/Sep/2001
 
	Glucas 2.8b doesn't compile under Mac OS X. (Missing <sys/time.h>
	in gio.c)
	Problem reported by Jonas Hansbo.
 
	v.2.8c 10/Oct/2001
 
	- Most of the FFT code has been written specially for IA64 (anyway 
	standard C code) and the use of Intel compiler has doubled the speed 
	for that platform. Now is the fastest Glucas runner (per clock). 
	The improvement is not only because of a rewritten code, a new preload
	scheme has been introduced. Basically it loads in a loop cycle what we
	need in the next one. It is only efficient in processors with lot of 
	registers as Itanium's (128 FP registers). To active this code in
	IA64 architecture compile with -DY_ITANIUM
 
	-Klaus Kastens has improved the timing routines. The elapsed real 
	time now is measured with a millisecond precision. 
 
	-Klaus also has written a new output format for verbose information. 
	Now there is information about the percentage of processor time spent
	by Glucas.
 
	-Some small changes for other platforms. No significant performance 
	improvements.
	
	-Fixed a rare memory bug affecting to huge selftest under FreeBSD
	(reported by Gregory Matus). Better memory management.
 
	-Fixed signal code for FreeBSD/x86 and Mac OS X.
 
 
	v 2.8d 02/Jan/2001
 
	-Glucas has been partially recoded to get multihtreading capabilities
	using Posix threads, OpenMP (a collection of libraries and pragma 
	directives is emerging as standar), and Sun WorkShop MP C (the Sun 
	multiprocessing version). No changes in single processor code.
 
 
	v 2.9.0 31/may/2002
 
	-Better task management. The queue files now includes more information
	about task to do, not only the exponent.
 
	-Introduced initial random shifts capabilities. Now Glucas can 
	double-check exponents already checked with Glucas (with different
	intial random shifts). For first Lucas-Lehmer test, the random shift
	is 0. For a double-check task, there is non zero random shift. 
	The command line options '-d' and '-D=n' has been added.
 
	-The range of exponents has been doubled. Now Glucas can work with 
	exponents up to 156000000. It will take a lot to complete this range :)
	Buildin seftests from 4608 to 8192 K FFT-runlength has been added 
	(selftest 27 to 31)
	
	-Improved configure script. Now it brings fully support for Intel 
	C++ compiler.
 
	-Improved documentation files. Now it includes texinfo files. It can be
	transformed to 'dvi', 'pdf' or 'html' format easily.
 
	-Fixed some small bugs.
	*********************************************************************
 
 
   ACKNOWLDEGMENTS: (in chronological order)
 
   This main file is an adaptation of Richard Crandall lucdwt.c, Sweeney 
   MacLucasUNIX.c and Will Edgington code. There are few things mine own in 
   it. The FFT routines and some special Lucas-Lehmer test routines it calls 
   are all from my hand. All the interface code is written by Will Edgington,
   I think.  
 
   Thanks to Erns W. Mayer to suggest me to write my own FFT package, to do 
   the first timings, to give me a telnet user account in his machines, and 
   by some ideas I borrow from his Mlucas code.
 
   Thanks to George Woltman, the core of all that Mersenne Prime search. 
   He has written the fastest code for Pentium I ever seen. After all, my 
   high language code only reach half performance than George's (tuned 
   carefully in assembler). 
 
   The alignment of doubles on stack, using GNU/GCC compiler on intel x86
   has got me a headache. It is almost the most important thing which 
   affect to performance. Thanks to FFTW developer team to point me on the
   solution. Actually, Glucas uses its tricks to align the stack.
     
   B.J. Beesley has made a great QA work. He also told me the way we can
   include assembler prefetch hints in Compaq-c compiler for alphas, and 
   reported some small bugs. 
 
   Klaus Kastens has made a deep work for ppc processors. His suggestions
   help me to improve some parts of code. Now he is a developer of glucas
   project. He is really an invaluable help.
 
   Tom Cage has made an excellent work of beta tester and porting Glucas
   to iMac.
 
   Thomas Perrier also helped tuning and building binaries for ppc processors.
 
   Gregory Matus, who reported memory problems under OpenBSD.  
 
   Tony Reix, by some improvement in Posix Threads code.
 
   And, over all, thanks to my wife Ana. Her patience has made possible 
   to write this in my spare time.
    
   THE PROGRAM:
   This program is designed to make the Lucas-Lehmer primality test to 
   Mersenne numbers. Is all written in C, using intensively the macro
   facilities, and being able to run in any system with a C-compiler.
   Actually, it runs in a Msdos machine compiled with the ancient  
   BORLAND Turbo  C++ (nevertheless, for small exponents).
    
   The memory requirements is about q bytes, where q is the exponent of
   mersenne number being tested. (m(q)).
 
   To see more about it, read the comments in yeafft.h file. There a lot
   of features it can be tuned for a particular system. You also must read 
   the file README.Glucas.htm
 
   Please, let me know any problem.
 
   Guillermo Ballester Valor
   gbv@oxixares.com
 
   Klaus Kastens
   kiste@bawue.de
*/

/* Include Files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <errno.h>

/* To use CodeWarrior C++ in Macintosh, as Tom Cage suggested */
#ifdef macintosh
# include <console.h>
#endif

/* some includes needed */
#include "version.h"
#include "gsetup.h"

#if (defined(HAVE_UNISTD_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
# include <unistd.h>
#endif


/* include needed to use yeafft package */
#include "yeafft.h"
#include "ydebug.h"

/* glucas include */
#include "glucas.h"

/* some definitions needed to Edgington Mers. package */
#define MASK Y_MASK
#define SHIFT Y_SHIFT
#define kErrLimit 0.40
#define init_kErrLimit 0.35
#define init_iterations 131072
#define starting_iterations (Y_SELFTEST_ITERS >> 1)
/* kErrChkFreq need to be a power of two */
#ifdef Y_SECURE
# define kErrChkFreq 1
#else
# define kErrChkFreq 64
#endif
#define kErrChk 1
#define kLast 2


/************************ definitions ************************************/
const char RCSglucas_c[] = "$Id$";
const char program_name[] = "Glucas"; /* for perror() and similar */
char version[80];

/* global variables needed */
BIG_DOUBLE high, low, highinv, lowinv;
BIG_DOUBLE SumLast = 0.0, SumIn = 0.0, SumOut = 0.0, SumNorm, ErrLimit;
BIG_DOUBLE Gsmall, Gbig, Hsmall, Hbig, Err, Y_XBITS;
BIG_DOUBLE *two_to_phi, *two_to_minusphi, *pttp=NULL, *pttmp=NULL, *PX=NULL;
UL b, c, UPDATE, SHIFT_UPDATE, save_iterations, QA_save, smode = 2,
    last_inf, first_giter;

char guser[SIZE_PATH]="Anonymous", chkpnt_s[SIZE_PATH], chkpnt_t[SIZE_PATH],
                      bits[64], chkpnt_u[SIZE_PATH], resfile[SIZE_PATH], queuefile[SIZE_PATH],
                      wdirectory[SIZE_PATH], UserID[SIZE_PATH], UserPWD[SIZE_PATH],
                      ComputerID[SIZE_PATH], ProxyHost[SIZE_PATH], DaysofWork[16],speed[16],
                      hours[16], WorktodoFile[SIZE_PATH], Machine_idr[SIZE_PATH];

FILE *glucasout;
/* Flags */

int Verbose_flag, Alternative_output_flag, Iteration_output,
Last_error_flag, Time_flag, Only_check_flag, Check_iteration,
Selftest_flag = 0, User_info,terminated = 0, ecf=kErrChk, ecfn = 0,
                /*Use_primenet,*/ CPUHours = 0, CPUSpeed = 0, pos;

#ifdef Y_DEBUG
/* debug global vars */
unsigned int TP_TIMEH = 0, TP_TIMEL = 0, TP_COUNT = 0, TP_TIMEMAX = 0,
                                      TP_COUNTMIN = 0, TP_TIMEMIN=0xFFFFFFFF;
#endif
#if defined(__GNUC__)
int BAD_STACK = 0;
#endif

/* Current iteration, Number of exponents queued and array of exponents */
UL Y_ITER = 0, Y_EXPONENT, Y_NQ = 0;
struct gtask Y_QD[Y_MAX_QD];

/* Current shift bit */
UL Y_SBIT0 = 0, Y_SBIT, Y_SFORCE = 0;

/*
   Y_KILL is to force a kill every Y_KILL iterations when Y_KILL != 0
   This is only to check i/o and save files
*/
int Y_KILL = 0;

/*
   Y_RESIDUE_EQU_TWO is 1 when a residue is 2 
*/
int Y_RESIDUE_EQU_TWO = 0;

/*
  Y_PASS2 marks the begining of pass2. It is set by 
  y_check_threshold_pass2()
*/
int Y_PASS2;

/* The last saved residue mod (2^32 - 1) */
UL Y_SUMCHECK;

/*
  Force next FFT runlength. To be save
*/
int FORCE_NEXT_FFT = 0;

/* The path of some control files */
char inifile[SIZE_PATH], inputfile[SIZE_PATH], name[SIZE_PATH],
/*pnetfile[SIZE_PATH],sdate[32],*/ gbuf[8 * SIZE_PATH];

/********  The TRICK to round to nearest *******************************/
/*
This plays with the internal hardward round to nearest when we add a
small number to a very big number, called bigA number. For intel fpu
this is 3*2^62. Then you have to sustract the same number (named different 
to confuse the clever compilers). This makes rint() very fast.
 
If you want to use this, define TRICKY_ROUND.
*/

/*#define TRICKY_ROUND*/

#if defined(TRICKY_ROUND)
BIG_DOUBLE bigA, bigB;
# if defined(Y_USE_SSE2)
Y__M128D MM_bigA, MM_bigB;
# endif /* Y_USE_SSE2 */
#endif /* TRICKY_ROUND */

# if defined(Y_USE_SSE2)
Y__M128D MM_bc[4], MM_auxt[4], MM_inv[4], MM_c, MM_Hsmall;
#endif /* Y_USE_SSE2 */


/* Macro to assign the proper name to ini files */

#define assign_name_inifiles(_filename, _defaultname)\
  if (_filename[0] == '\0')                          \
    {                                                \
      if(wdirectory[0] != '\0')                      \
	{                                            \
	  strcpy (_filename, wdirectory);            \
	  strcat (_filename, DIR_SLASH);             \
	}                                            \
      strcat (_filename, _defaultname);              \
    }                                                \
  else if(wdirectory[0] != '\0')                     \
    {                                                \
      strcpy (name, _filename);                      \
      strcpy (_filename, wdirectory);                \
      strcat (_filename, DIR_SLASH);                 \
      strcat (_filename, name);                      \
    }


/**************************************************************
 *
 *      Main Function
 *
 **************************************************************/

int main(int argc, char *argv[])
{
  UL q = 0L, n, n2=0, j = 1L, first=1L, last = 2L, last_iter_big = 0/*,
                                        mask = MASK, shift = SHIFT*/;
  struct gtask *qd;
  size_t size, k;
  BIG_DOUBLE  *x =NULL, /* *px=NULL, */ bits_per_double,
                  last_big_err = 0.0;
#if defined(SUM_CHECK)

  BIG_DOUBLE sumerr, last_sumerr = 0.0;
  UL last_itererr = 0;
#endif

  double timetmp;
  int restarting = 0, newfft =0;
  int mode = 0, nonice=0, iopt, nq=0, i, lmode, save=0, flag;


  /* Patch for stack alignement in some x86 GCC versions */
  HACK_CHECK_ESP();

  /* set stdout line buffered */
  setvbuf( stdout, NULL, _IOLBF, BUFSIZ);

  /*inits random seed */
  srand(time(0));

  UPDATE = Y_UPDATE;
  SHIFT_UPDATE = Y_SHIFT_UPDATE;
  two_to_phi = NULL;
  two_to_minusphi = NULL;
  qd = &Y_QD[0];

#if defined(TRICKY_ROUND)
  /*
     if defined SSE2 code the constants are set when calling tricky()
     note than double bigA could not be the same than the vector-SSE2 
     value.
  */
  bigA = tricky();
  bigB = bigA;
#endif /* TRICKY_ROUND */

  inifile[0] = '\0';
  resfile[0] = '\0';
  queuefile[0] = '\0';
  wdirectory[0] = '\0';
  ProxyHost[0] = '\0';
  WorktodoFile[0] = '\0';
  ComputerID[0] = '\0';
#ifdef macintosh

  argc = ccommand(&argv);
#endif

  /* read the arguments using getopt or a simplified version of mine*/
#ifdef _PTHREADS

  while ((iopt = getopt(argc, argv, "fanocdhvzi:q:r:s:S:w:W:D:T:N:K:"))!=-1)
#else

  while ((iopt = getopt(argc, argv, "fanocdhvzi:q:r:s:S:w:W:D:N:K:"))!=-1)
#endif

    switch (iopt)
      {
      case 'a':	/*The mode how the work is inserted in queue is append*/
        mode = 0;
        break;

      case 'f': /*The mode is insert: Made it just now! */
        mode = 1;
        break;

      case 'n': /* No nice run_mode */
        nonice = 1;
        break;

      case 'o': /* Old style save files format */
        smode = 0;
        break;

      case 'c': /* Universal compatible style save files format */
        smode = 2;
        break;

      case 'd': /* Is a double check, use random shift mode */
        Y_SBIT0 = 1;
        break;

      case 's': /* Brief self test (100 iters / exponent)*/
        Y_SELFTEST_ITERS = 100;
        if((strlen (optarg) > 0) && (selftest (optarg) >= 0) && self_inifile())
          {
            Selftest_flag = 1;
            nonice = 1;
            smode = 2;
            strcpy (queuefile, "selftest.que");
            strcpy (inifile, "selftest.ini");
            strcpy (resfile, "selftest.res");
          }
        break;

      case 'S': /* Larger self test (1000 iters / exponent)*/
        Y_SELFTEST_ITERS = 1000;
        if((strlen (optarg) > 0) && (selftest (optarg) >= 0) && self_inifile())
          {
            Selftest_flag = 1;
            nonice = 1;
            smode = 2;
            strcpy (queuefile, "selftest.que");
            strcpy (inifile, "selftest.ini");
            strcpy (resfile, "selftest.res");
          }
        break;

      case 'i': /* The ini file is ... */
        if(strlen (optarg) >0)
          strcpy (inifile, optarg);
        break;

      case 'q': /* The queue file is ... */
        if(strlen (optarg) >0)
          strcpy (queuefile, optarg);
        break;

      case 'r': /* The result file is ... */
        if(strlen (optarg) >0)
          strcpy (resfile, optarg);
        break;

      case 'w': /* The worktodo file is ... */
        if(strlen (optarg) >0)
          strcpy (WorktodoFile, optarg);
        break;

      case 'W': /* The working directory is ... */
        if(strlen (optarg) >0)
          strcpy (wdirectory, optarg);
        break;

      case 'D': /* force initial shift to .... */
        if(strlen (optarg) >0)
          Y_SFORCE = atoi (optarg);
        Y_SBIT0 = 1;
        break;

      case 'K': /* force to kill every Y_KILL iterations .... */
        if(strlen (optarg) >0)
          Y_KILL = atoi (optarg);
        break;

      case 'N': /* force Nice priority to .... */
        if(strlen (optarg) >0)
          Nice = atoi (optarg);
        break;

#ifdef _PTHREADS

      case 'T': /* Use T threads */
        if(strlen (optarg) >0)
          Y_NTHREADS = atoi (optarg);
        break;
#endif

      case 'v': /* Print version and exits */
        printf ("%s %s %s\n", program_name, version_string, build_string);
        return 0;

      case 'h': /* Help flag */
        print_usage ();
        return 0;

      case 'z': /* use next FFT */
        FORCE_NEXT_FFT = 1;
        break;

      case '?':
        print_usage ();
        return 0;
      }

#if defined(HAVE_CONFIG_H) && defined(HAVE_CHDIR)
  if (wdirectory[0] != '\0' && chdir( wdirectory))
    {
      fprintf (stderr, "Glucas: Cannot move to directory %s\n", wdirectory);
      exit (EXIT_FAILURE);
    }
  wdirectory[0] = '\0'; /* So it is like no -W flag */
#endif

  /* Set the default values if not in args */
  assign_name_inifiles (inifile, "glucas.ini");
  assign_name_inifiles (queuefile, "glucas.que");
  assign_name_inifiles (resfile, "results");
  assign_name_inifiles (WorktodoFile, "worktodo.ini");

  /*
    See whether is more work file to add to Glucas queue 
    Fill qd[] with the queued exponents and tasks to test         
    nq is the number of exponents in queue               
    qd[0].exponent will be the first exponent we need to test     
  */
  if (argc == optind)
    {
      inputfile[0] = '\0';
      nq = manage_queue (mode, inputfile, qd);
    }
  else
    {
      iopt=argc - 1;
      while(iopt >= optind)
        {
          nq = manage_queue (mode, argv[iopt], qd);
          iopt--;
        }
    }

  /* Set the program nice (if nonice==0) and the handler signals */
  setup (nonice);

  /* Set other options with the optional ini file */
  read_inifile (inifile);

  /* To set the info  version */
  sprintf (version, "%s %s %s", program_name, version_string, guser);

  /* Read the Primenet parameters */
  /*if(Use_primenet)
    { 
      if(wdirectory[0] != '\0')
  {
   strcpy (name, "primenet.ini");
   strcpy (pnetfile, wdirectory);
   strcat (pnetfile, DIR_SLASH);
   strcat (pnetfile, name);
  }
      else strcpy(pnetfile, "primenet.ini");
      read_primenet_inifile (pnetfile);
  */
  /* Get exponents whether no enough work */
  /*  if((remaining_work()/86400.0) < (double) atoi(DaysofWork))
  {
   if(Verbose_flag) 
     {
       fprintf(glucasout,"Getting exponents from server\n");
       if(Alternative_output_flag==2)
  {
    printf("Getting exponents from server\n");
  }
     }
   if(primenet_checkout(WorktodoFile))
     fprintf(stderr,"Glucas : Warning : Unable to get exponents from Primenet\n");
   else nq = manage_queue( mode, WorktodoFile, qd);
  }
    }
  */
  /* Initialize timer */
  if(Time_flag)
    timetmp = sec_elapsed ();

  /*
    This is the big loop. 
    Terminate is set to one with a signal. Then the Lucas-Lehmer loop 
    ends setting 'terminated' to one, saving the intermediate files, and
    exiting. 
    There are some return points in it when some errors are detected  
  */
  while (!terminate)
    {
#if defined(_PTHREADS)
      PTHREAD_trace(fprintf(stderr, "Main Thread: %p: in big loop (while !terminate)\n", (void *) pthread_self()););
#endif

      do
        /*
           while (restarting) 
           restarting is set to one when begins again because of an error
           round off error exceed, in other case, the do {} ends when a job
           finishes normaly (and restarting remains 0). 
           Then it delete the save_files and look for new jobs with a 
           new entry to this do {}.
           If no more work in queue, ginput()  returns EXIT_SUCCESS and the 
           program exits normally.
        */
        {
#if defined(_PTHREADS)
          PTHREAD_trace(fprintf(stderr, "Main Thread: %p: in loop (do ... while restarting)\n", (void *) pthread_self()););
#endif

          /*
             read the next job to do or the last save file if exists 
             If a save file exists and n2 != 0 when calling ginput, then 
             the data is allocate using a x[n2] float array. It is  useful to
             change the FFT run length at run time if errors are too high
          */
          switch ( ginput (nq, qd, &q, &n2, &j, &Err, &x))
            {
              /* no more input */
            case 0:
              sprintf (gbuf, "Terminated all the queued job in file: %s.\n",
                       queuefile);
              write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
              return (EXIT_SUCCESS);

              /*something wrong; error message, if any, already printed */
            case 1:
            default:
              print_date ();
              return (EXIT_FAILURE);

              /* continuing work from a partial result */
            case 2:
              if(Y_SBIT0 != qd[0].parm1)
                {
                  qd[0].parm1 = Y_SBIT0;
                  pos = sprintf (gbuf, "Warning: Used save file are not generated from %s:\n", queuefile);
                  sprintf(gbuf + pos, "   Different initial shift in queue file and savefile\n");
                  write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                }
              i = read_update (q);
              Y_EXPONENT = q;
              Y_SBIT = compute_sbit (Y_SBIT0, j, q);
              if(i > 3)
                {
                  if(strcmp (ComputerID, Machine_idr))
                    {
                      /*
                      The starting is from other machine results, 
                      we must to activate the roundoff check errors 
                      via First_giter
                      */
                      pos = sprintf (gbuf, "Detected a change in Glucas environment.\n");
                      sprintf (gbuf + pos,"Roundoff check will be activated during next 2^17=131072 iterations.\n");
                      write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);

                      first_giter = j;
                      strcpy (Machine_idr, ComputerID);
                      write_update (q);
                    }
                }
              else
                {
                  /*
                     The starting is from other machine results,
                     we didnt't have a g* file, 
                     we must to activate the roundoff check errors via 
                     First_giter
                  */
                  pos = sprintf (gbuf, "Detected a change in Glucas environment.\n");
                  sprintf (gbuf + pos,"Roundoff check will be activated during next 2^17=131072 iterations.\n");
                  write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                  first_giter = j;
                  strcpy (Machine_idr, ComputerID);
                  write_update (q);
                }

              /*restarting = 1;*/
              /* not the usual sense of restarting (FFT errors too high) */
              n = y_init (n2);
              size = addr (n2);
              newfft = 0;
              pos = sprintf (gbuf, "Restarting from iteration %ld .Exponent %ld.", j, q);

              if(Y_SBIT0)
                pos += sprintf (gbuf + pos, " Shifted %ld\n", Y_SBIT0);
              else
                pos += sprintf (gbuf + pos, "\n");

              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
              break;

              /* Start from iteration 1. No save file */
            case 3:
              read_update (q);
              Y_EXPONENT = q;
              strcpy(Machine_idr, ComputerID);
              /* Set shift mode if task is DoubleCheck */
              if (qd[0].task == DoubleCheck)
                Y_SBIT0 = qd[0].parm1;
              else
                Y_SBIT0 = 0;
              if( !newfft ) /* Changed FFT-Length ? */
                {
                  /* Regresion taken from EW Mayer's timing page */
                  bits_per_double = (BIG_DOUBLE)25.33 - 0.37 * log ((BIG_DOUBLE)q);
                  n2 =(UL) (((BIG_DOUBLE)q-1.0)/bits_per_double + 1.0);
                  if (FORCE_NEXT_FFT)
                    n2 += (n2 >> 2);
                }
              else
                newfft = 0;

              n = y_init (n2);

              /* now n is the fft length (in size of complex) */
              n2 = n + n;
              size = addr (n2);

              if (PX != NULL)
                free((char *) PX);

              PX = ALLOC_DOUBLES (size);

              if( PX == NULL)
                {
                  fprintf (stderr, "glucas: No memory to allocate data.\n");
                  exit (EXIT_FAILURE);
                }

              x = ALIGN_DOUBLES (PX);
              terminated = 0;

              /* is the exponent too small ? */
              if(Y_NRADICES < 2)
                {
                  fprintf (stderr, "Exponent too small to this version.\n");
                  exit (EXIT_FAILURE);
                }

              /* see whether UPDATE is consistent with FFT plan */
              if( (Y_LRIGHT[1] << 1) < UPDATE)
                {
                  UPDATE = (Y_LRIGHT[1] << 1);
                  SHIFT_UPDATE = (Y_K + 1);
                }
              write_update (q);

              /* Iteration 0 is 4, 1 is 14 and so on */
              for (k = 0; k < size; k++)
                x[k] = 0.0;

              set_L0 (x, q, n2);
              j = (UL)0;
              Err = 0.0;

              if(Selftest_flag)
                {
                  self_info();
                  pos = sprintf(gbuf, "Selftest %d (%d K FFT-runlength).", nself(q), nalter(q));
                  sprintf (gbuf + pos, " %d iterations for M%ld...\n", Y_SELFTEST_ITERS, q);
                  fputs (gbuf, glucasout);
                  fputs (gbuf, stdout);
                }
              pos = sprintf (gbuf, "Going to work with exponent %ld\n", q);
              pos += sprintf (gbuf + pos, "Starting from iteration 1. Exponent %ld.", q);

              if(Y_SBIT0)
                pos += sprintf (gbuf + pos," Initial shift %ld.\n", Y_SBIT0);
              else
                pos += sprintf (gbuf + pos, "\n");

              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
              break;
            }
          restarting = 0;

          /* Init the Lucas-Lehmer variables */
          Y_WORK = x;
          init_lucas (q, n2);
#if defined(SUM_CHECK)

          SumLast = 0.0;
          SumOut = 0.0;
#endif

          /* The first iteration done, readed from a save file if j!=0 */
          first = j;
          last_inf = j;
          Y_ITER = j;

          /* Init the timers */
          sec_elapsed ();

#if defined(_PTHREADS) && defined(linux)

          if(y_user != NULL)
            {
              for(k = 0; k < Y_NTHREADS; k++)
                y_user[k] = 0.0;
            }
#endif

          sec_elapsed_user();

          /* the last iteration to do in the primary loop */
          if(Only_check_flag)
            {
              last = ((UL)Check_iteration < q - (UL)2) ?
                     (UL)Check_iteration : q - (UL)2 ;
            }
          else
            {
              last = q - (UL)2;
              if( !primeq(q) )
                {
                  pos = sprintf (gbuf, "M%ld should not be queued for a complete Lucas-Lehmer test.\n", q);
                  sprintf(gbuf + pos, "%ld is not prime. Skipping job!\n", q);
                  write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                  sprintf (gbuf, "M%ld is not prime. Exponent not prime!\n", q);
                  write_resultsfile (gbuf, resfile, 1);
                  last = 0;
                }
            }

          /*
             Lucas Lehmer test loop. No end up to last or a signal interrupt.
             'terminated' is set to one when a save file has written with 
             success because of a break signal (which sets terminate to 1).
             'restarting' is set to one when an error exceed occurs.
             'j' is the last iteration done in the loop.  
          */
          for ( ; !terminated && !restarting && j < last; )
            {
#if (kErrChkFreq > 1)
              if ((j & ((UL) kErrChkFreq - (UL)1)) == (UL)0 || (j - first_giter) <  init_iterations)
#endif

                flag = (UL) ecf;
#if (kErrChkFreq > 1)

              else
                flag = (UL) ecfn;
#endif
              /* with this line there will be at least starting_iterations
              roundoff checks. It will avoid potential problems when the 
              hardware contest changes after a restart */
              if((j - first) < (UL)starting_iterations)
                flag = (UL) 1;

              /* patch to kill every Y_KILL iterations (and check the save files process) */
              if (Y_KILL && j != first && (j % Y_KILL) == 0)
                terminate = 1;

              /*
              See the kind of iteration to do:
              Normal. lmode = 0. A normal iteration.
              First. lmode = -1. The first iteration in the loop.
              Last. lmode = 1. The last iteration in the loop.
              save. lmode = 2. Save a residue an do a normal iteration.
              break. lmode = 9. Save a  residue inmediately without doing
              a new iteration  
              */
              if( j == first)
                lmode = -1;
              else if( j == (last - 1))
                lmode = 1;
              else if( terminate )
                lmode = 9;
              else if( save )
                lmode = 2;
              else
                lmode = 0;

              /*
              Do an iteration.
              if i==1 the iteration has actually done
              if i==0 it only prepared the x[] array to be saved
              if i==-1 an error has occurred
              */
              i = (UL)lucas_square (q, x, n, j, lmode, flag);
              if (i<0)
                return (errno);
              j += (UL) i;
              Y_ITER = j;

              /* See whether it has to write a save file */
              save = (((save_iterations != 0) &&
                       ((j % save_iterations) == (UL) 0 ) &&
                       (j < last) &&  (j>0))) ||
                     ((QA_save !=0) && (j % QA_save == (UL) 0) && (j>0));

              /* See whether it has to output information */
              if ((Verbose_flag) && ( j % Iteration_output == 0))
                iteration_inf (j, last);

#if defined(SUM_CHECK)

              if (SumLast != 0.0 && !terminated)
                {
                  sumerr = fabs ( SumLast * SumLast - SumOut * SumNorm);
                  /*printf("iter=%ld Sumerr = %lf\n", j, sumerr);*/
                  if (sumerr > ErrLimit)
                    {
                      if (last_sumerr != sumerr || j != last_itererr )
                        {
                          /* a sumcheck error, new */
                          pos = sprintf (gbuf, "Glucas: Sumerr = %lf exceeding limits",sumerr);
                          pos += sprintf (gbuf + pos, " at iteration %lu.\n", j);
                          pos += sprintf (gbuf + pos, "Restarting from last save file\n");
                          write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                          last_sumerr = sumerr;
                          last_itererr = j;
                        }
                      else if (last_sumerr == sumerr && j == last_itererr)
                        {
                          pos = sprintf (gbuf, "Glucas Error!: Sumerr = %lf exceeding limits", sumerr);
                          pos += sprintf (gbuf + pos, " at iteration %lu.\n", j);
                          pos += sprintf (gbuf + pos, "This is a reproducible error\n");
                          sprintf(gbuf + pos, "Glucas Sumcheck error!. Suiciding at iteration %ld\n",j);
                          write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                          fprintf (stderr, "Glucas Sumcheck error!. Suiciding at iteration %ld...\n",j);
                          exit(EXIT_FAILURE);
                        }
                      restarting = 1;
                    }
                }
              SumLast = SumIn;
              SumIn = 0.0;
              SumOut = 0.0;
#endif /*SUM_CHECK*/

              /*
              If error greater than limit, change the parameters.
              Prior be sure it is a reproducible sotfware error.
              First, try to increase the UPDATE rate. If no good 
              results, then increase the FFT length. 
              */
              if (((Err > kErrLimit) && (j > init_iterations)) ||
                  ((Err > init_kErrLimit) && (j <= init_iterations)))
                {
                  if(j < (UL) 10) /* If  error is too high in first iterations abort */
                    {
                      fprintf (stderr, "Glucas panic!. Something is wrong. Suiciding at iteration %ld...\n",j);
                      exit(EXIT_FAILURE);
                    }

                  if (last_iter_big == 0)
                    {
                      pos = sprintf (gbuf, "Glucas: Max. roundoff error exceeded");

                      if (j> init_iterations)
                        pos += sprintf (gbuf + pos, " at iteration %lu. (%f > %f)\n", j, Err, kErrLimit);
                      else
                        pos += sprintf (gbuf + pos, " at iteration %lu. (%f > %f)\n", j, Err, init_kErrLimit);
                      sprintf (gbuf + pos, "Restarting from last save file if it exists.\n");
                      write_glucasout (gbuf, Alternative_output_flag, 0, Verbose_flag);
                    }

                  if (last_big_err == Err || j == last_iter_big )
                    {
                      /* Is a lack of software precision */
                      {
                        write_glucasout ("Last error was a reproducible error.\n",
                                         Alternative_output_flag, 0, Verbose_flag);
                        if(UPDATE > (UL) 4)
                          {
                            UPDATE >>= 1;
                            SHIFT_UPDATE--;
                            write_update (q);
                            write_glucasout ("Restarting. Accuracy increased.\n",
                                             Alternative_output_flag, 0, Verbose_flag);
                            n2 = 0;
                          }
                        else
                          {
                            UPDATE = Y_UPDATE;
                            SHIFT_UPDATE = Y_SHIFT_UPDATE;
                            write_update (q);
                            write_glucasout ("Restarting. Changed FFT length.\n",
                                             Alternative_output_flag, 0, Verbose_flag);
                            newfft = 1;
                            n2 += 2;
                            n = y_init(n2);
                            n2 = n + n;
                          }
                        last_big_err = 0.0;
                        last_iter_big = 0;
                      }
                    }
                  else if (last_big_err == 0.0)
                    {
                      last_big_err = Err;
                      last_iter_big = j;
                    }
                  else
                    {
                      write_glucasout ( "WARNING: it was a non-reproducible (hardware) error\n",
                                        Alternative_output_flag, 0, Verbose_flag);
                      last_big_err = 0.0;
                      last_iter_big = 0;
                    }
#if defined(_PTHREADS)
                  y_cancel_all_threads();
#endif

                  restarting = 1;
                }
            }/*end of a normal L-L cycle. It will perform another iter only if
          	       !terminated && !restarting && j < last .Otherwise continue in the following
          	       code */

          if ( last_iter_big && j > last_iter_big)
            {
              write_glucasout ( "WARNING: it was a non-reproducible (hardware) error\n",
                                Alternative_output_flag, 0, Verbose_flag);
              last_big_err = 0.0;
              last_iter_big = 0;
            }

          if (restarting)
            continue;
          else if (j < last)/* not done, but need to quit, for example,
            				by a SIGTERM signal*/
            {
              sprintf (gbuf, "Terminating at iteration %ld by break signal.\n",j);
              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);

              return((write_check_point (q, n2, j, Err, x, smode)
                      <= 0) ? errno : 0);
            }
        }
      while (restarting);

      if(iszero (x, n2))
        {

          if(isnewprime (q))
            {
              /* NEW PRIME  DISCOVERED !!!! */
              /* print_date2(sdate); */
              /*if(Use_primenet)
              {
              if(ComputerID[0]!='\0') 
              pos = sprintf(gbuf,"%sUID: %s/%s, M%ld is prime!\n",
              sdate,UserID,ComputerID,q);
              else
              pos = sprintf(gbuf,"%sUID: %s, M%ld is prime!\n",
              sdate,UserID,q);
              }
              else*/
              pos = sprintf (gbuf, "UID: %s/%s, M%ld is prime! G29: 00000000,%ld\n",
                             UserID, ComputerID, q, Y_SBIT);
              pos += sprintf (gbuf + pos, "NEW MERSENNE PRIME DISCOVERED !!!!\n");
              pos += sprintf (gbuf + pos, "PLEASE, send email to woltman@alum.mit.edu AND\n");
              pos += sprintf (gbuf + pos, "gbv@oxixares.com with the file \"%s_last\" attached.\n",chkpnt_s);
              pos += sprintf (gbuf + pos, "(You can find it in your working directory)\n");
              write_glucasout(gbuf, Alternative_output_flag, 1, Verbose_flag);
              write_resultsfile(gbuf, resfile, 1);
              /* copy last save_file */
              copy_last_save_file ();
            }
          else
            {
              pos = sprintf (gbuf, "M%ld is a known Mersenne prime!\n", q);
              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
              write_resultsfile (gbuf, resfile, 1);
            }
        }
      else if (res64 (x, q, n2, bits) ==  16  )
        {
          bits[16] = '\0';

          if(last ==  (q - (UL) 2))
            {
              UL aux;
              aux = residue64_sumcheck_shift (q, n2, x, &Y_SUMCHECK);
              /*print_date2(sdate);*/
              /*if(Use_primenet)
              {
              sprintf(gbuf,"M%ld is not prime. Res64: %s. %s\r\n",
              q,bits,version);
              }
              else*/
#if ULONG_MAX > 0xFFFFFFFF

              sprintf (gbuf, "M%ld is not prime. Res64: %s. G29: %.16lX,%ld\n",
                       q, bits, Y_SUMCHECK, Y_SBIT);
#else

              sprintf (gbuf, "M%ld is not prime. Res64: %s. G29: %.8lX%.8lX,%ld\n",
                       q, bits, aux, Y_SUMCHECK, Y_SBIT);
#endif
              /*fputs(sdate,glucasout);*/
              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
              write_resultsfile (gbuf, resfile, 1);
              /* Send the results to primenet */
              /*if(Use_primenet)
              {
              if(Verbose_flag) 
              {
               fprintf(glucasout,"Sending results to server\n");
               if(Alternative_output_flag==2)
              {
              printf("Sending results to server\n");
              }
              }
              primenet_checkin(gbuf);
              }*/
            }
          else
            {
              pos = sprintf (gbuf, "M%ld. Iteration %ld. Res64: %s. %s\n", q, last,
                             bits, version);
              write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
              write_resultsfile (gbuf, resfile, 1);

              if(Selftest_flag)
                {
                  if(self_success (q, bits))
                    pos = sprintf (gbuf, "Selftest %d success!\n", nself(q));
                  else
                    pos = sprintf (gbuf, "***** Selftest %d failed!\n", nself(q));

                  write_glucasout (gbuf, Alternative_output_flag, 1, Verbose_flag);
                  write_resultsfile (gbuf, resfile, 1);
                  fputs (gbuf,stdout);
#ifdef Y_DEBUG

                  printf (" DEBUG RESULTS :\n Max=%ld Min=%ld at %ld of %ld\n", TP_TIMEMAX, TP_TIMEMIN, TP_COUNTMIN, TP_COUNT);

                  if(TP_COUNT)
                    printf (" AVER =%6.0f\n",((double) TP_TIMEH * 65536.0 * 65536.0 + (double) TP_TIMEL)/(double) TP_COUNT);
#endif

                }
            }
        }
      /*print_date();*/
      delete_update (q);
      nq = delete_first_job (qd);
      /*if(Use_primenet)
      {
      if(delete_worktodo_exponent( WorktodoFile , q))
      exit(EXIT_FAILURE);
      Y_ITER=0;
      if((remaining_work()/86400.0) < (double) atoi(DaysofWork))
      {
       if(Verbose_flag) 
      {
      fprintf(glucasout,"Getting exponents from server\n");
      if(Alternative_output_flag==2)
      {
        printf("Getting exponents from server\n");
      }
      }
       if(primenet_checkout(WorktodoFile))
      fprintf(stderr,"Glucas : Warning : Unable to get exponents from Primenet\n");
       else nq = manage_queue( mode, WorktodoFile, qd);
      }
      }
      */
      remove
        (chkpnt_s);
      remove
        (chkpnt_t);
#if defined(_PTHREADS)

      y_cancel_all_threads();
#endif

      free ((char *) PX);
      PX = NULL;
      n2 = 0;
      terminated = 0;
      Err = 0.0;
    }
  return (EXIT_SUCCESS);
}
/* $Id$ */

