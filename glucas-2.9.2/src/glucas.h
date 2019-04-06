/*$Id$*/
/*
 
   (c) 2000-2006 Guillermo Ballester Valor, Klaus Kastens

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

/**************************************************************************
 This part of file has the prototypes for special pass used in Lucas-Lehmer
test where last dit pass, carry and normalize, and first dif pass are made
in a single pass 
***************************************************************************/
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Mathematical constants in <math.h> only defined in Unix98 compatible */
#ifndef M_LN2
# define M_LN2 0.69314718055994530942
#endif

#if (!defined(HAVE_GETOPT) && defined(HAVE_CONFIG_H)) || defined(macintosh) || (defined(__MWERKS__) && defined(__INTEL__))
int getopt(int , char **, char *);
extern int optind, goptindex;
extern char  *optarg;
#endif

#if (!defined(HAVE_STRCHR) && defined(HAVE_CONFIG_H))
const char *strchr(const char *string, int ch);
#endif

#if defined(HAVE_LIMITS_H) || !defined(HAVE_CONFIG_H)
# include <limits.h>
/* Patch suggested by B.J.Beesley */
# if defined(SUN_V9_GCC) || defined(__sparc)
#  undef ULONG_MAX
#  define ULONG_MAX 0xFFFFFFFF
# endif
#endif

/* struct to define the task to do in a exponent */
struct gtask
  {
    UL  exponent; /* the exponent */
    int  task; /* the task, see the former macro defs */
    UL  parm1; /*
    		parameter 1 optional needed for task
    		Initial shift bit in a L-L test.
    	     */
    UL  parm2; /* parameter 2 optional needed for task */
  };

/* Max. number of exponents in queue */
#define Y_MAX_QD 128

/* some define macros with tasks */

#define Test  0
#define DoubleCheck 1

/* The size of strings to store the complete path of files */
#define SIZE_PATH 256

void init_lucas (UL , UL);

void set_L0 (y_ptr, UL q, UL n);

void generate_shift0 (struct gtask *);

int lucas_square (UL , y_ptr, UL , UL , int, int);

int primeq (UL);

void print_usage (void);

int y_fftf_pass1_lucas (y_ptr , y_size_t, y_size_t);

int y_fftf_pass1_lucas_mp (y_ptr , y_size_t, y_size_t, int);

void y_fftf_pass1_lucas_first (y_ptr, y_size_t);

void y_fftb_pass1_lucas_last (y_ptr, y_size_t);

void y_fftb_pass1_lucas (y_ptr , y_size_t , y_size_t, int);

void y_fftb_pass1_lucas_mp (y_ptr , y_size_t , y_size_t, int);

void lucas_auto_convolution (y_ptr , y_size_t);

void lucas_last_auto_convolution (y_ptr , y_size_t);

void lucas_first_auto_convolution (y_ptr , y_size_t);

void normalize (y_ptr ,UL ,UL);

void first_normalize (y_ptr, UL);

void last_normalize (y_ptr ,UL ,UL);

void y_normalize (y_ptr ,UL ,UL);

void y2_normalize (y_ptr ,UL , UL);

void y3_normalize (y_ptr ,UL , UL);

void y_last_normalize (y_ptr ,UL , UL);

int isnewprime (UL);

int iszero (y_ptr, UL);

void check_balanced (y_ptr, UL);

void balancetosdrep (y_ptr, UL);

int res64 (y_ptr, UL, UL, char *);

UL residue_char (y_ptr, UL, UL,unsigned  char *, UL, UL);

int write_update (UL);

int read_update (UL);

int delete_update (UL);

int print_date (void);

time_t print_date2 (char *);

time_t print_date3 (char *);

int write_glucasout (char *, int, int, int);

int write_resultsfile (char *, char *, int);

int validate_inifile_string (char *);

int read_inifile (char *);

int read_primenet_inifile (char *);

int primenet_checkin (char *);

int primenet_checkout (char *);

double sec_elapsed (void);

double sec_elapsed_user (void);

size_t inf_time (UL, char *, size_t);

int iteration_inf ( UL, UL );

int read_queue (struct gtask *, int);

int write_queue (struct gtask *, int);

int delete_first_job (struct gtask *);

int manage_queue (int, char *, struct gtask *);

int write_standar (UL, UL, y_ptr, UL *, UL *);

int read_standar (UL, UL, y_ptr, UL *, UL *);

int residue_sumcheck (UL, UL, y_ptr, UL *);

UL residue64_sumcheck_shift (UL, UL, y_ptr, UL *);

int read_check_point (UL, UL *, UL *, UL *, y_ptr, y_ptr *);

int write_check_point (UL, UL, UL, BIG_DOUBLE, y_ptr, int);

int write_interim_file (UL);

int copy_last_save_file (void);

int ginput (int, struct gtask *, UL *, UL *, UL *, y_ptr, y_ptr *);

int selftest (char *);

double selftest_parser (char *);

void virtual_CPU_parameters (double);

double remaining_work (void);

char * parse_checkout_response (char *);

char * parse_checkin_response (char *);

int delete_worktodo_exponent (char *, UL);

int self_inifile (void);

int self_success (UL, char *);

int self_info (void);

int nself (UL);

int nalter (UL);

size_t fwrite_UL_pad (UL, FILE *);

size_t fread_UL_pad (UL *, FILE *);

size_t fwrite_UL (UL, FILE *);

size_t fread_UL (UL *, FILE *);

size_t fwrite_twin (UL, UL, FILE *);

size_t fread_twin (UL *, UL *, FILE *);

size_t fwrite_UL_array (UL *u, size_t, FILE *);

size_t fread_UL_array (UL *u, size_t, FILE *);

void add_sumcheck (UL * , UL );

void add_sumcheck2 (UL * , UL *, UL, UL);

UL sumcheck32 (UL);

void lucas_dit_carry_norm_dif (y_ptr, UL, UL);

void dit_carry_norm_dif_2 (y_ptr, UL, UL);

void dit_carry_norm_dif_4 (y_ptr, UL, UL);

void dit_carry_norm_dif_5 (y_ptr, UL, UL);

void dit_carry_norm_dif_5_sse2 (y_ptr, UL, UL);

void dit_carry_norm_dif_6 (y_ptr, UL, UL);

void dit_carry_norm_dif_6_sse2 (y_ptr, UL, UL);

void dit_carry_norm_dif_7 (y_ptr, UL, UL);

void dit_carry_norm_dif_7_sse2 (y_ptr, UL, UL);

void dit_carry_norm_dif_8 (y_ptr, UL, UL);

void dit_carry_norm_dif_8_sse2 (y_ptr, UL, UL);

void dit_carry_norm_dif_9 (y_ptr, UL, UL);

void dit_carry_norm_dif_9_sse2 (y_ptr, UL, UL);

void dit_carry_norm_dif_10 (y_ptr, UL, UL);

void dit_carry_norm_dif_12 (y_ptr, UL, UL);

void dit_carry_norm_dif_14 (y_ptr, UL, UL);

void dit_carry_norm_dif_16 (y_ptr, UL, UL);

void substract_two (y_ptr, UL);

void substract_two_5 (y_ptr, UL);

void substract_two_6 (y_ptr, UL);

void substract_two_7 (y_ptr, UL);

void substract_two_8 (y_ptr, UL);

void substract_two_9 (y_ptr, UL);

void substract_two_10 (y_ptr, UL);

void substract_two_12 (y_ptr, UL);

void substract_two_14 (y_ptr, UL);

void substract_two_16 (y_ptr, UL);

#if defined(Y_USE_SSE2)
void substract_two_5_sse2 (y_ptr, UL);

void substract_two_6_sse2 (y_ptr, UL);

void substract_two_7_sse2 (y_ptr, UL);

void substract_two_8_sse2 (y_ptr, UL);

void substract_two_9_sse2 (y_ptr, UL);
#endif

void void_carry_norm_5_sse2 (void);

void void_carry_norm_6_sse2 (void);

void void_carry_norm_7_sse2 (void);

void void_carry_norm_8_sse2 (void);

void void_carry_norm_9_sse2 (void);

UL compute_sbit (UL, UL, UL);


/**************************************************************************
    THIS FILE INCLUDES MACROS USED IN SPECIAL LUCAS-LEHMER TEST          
    WHERE LAST BACKWARD DIF, CARRY AND NORMALIZATION, AND FIRST FORWARD
    DIF IS MADE WITH ONLY ONE MEMORY-PASS 
***************************************************************************/
/*************************************************
  Some global definitions .
 
  See original Dr. Crandall lucdwt.c code to learn more.
   high =2^bits_per_big_word
   low  =2^bits_per_small_word , high=low*2
   highinv = 1/high
   lowinv  = 1/low
   b = q % N
   c = N - b       
*/
extern BIG_DOUBLE high, highinv, low, lowinv, Gbig, Gsmall, Hbig, Hsmall, Err;
extern BIG_DOUBLE SumLast, SumIn, SumOut, SumNorm, ErrLimit;
extern BIG_DOUBLE *two_to_phi, *two_to_minusphi, *pttp, *pttmp, *PX;
extern UL b, c, UPDATE, SHIFT_UPDATE, save_iterations, QA_save, last_inf,
  first_giter, smode;
extern int Nice;
extern int Verbose_flag, Alternative_output_flag, Iteration_output,
  Last_error_flag, Time_flag, Only_check_flag, Check_iteration, Selftest_flag,
  User_info, terminated, ecf, ecfn, Use_primenet;
extern char guser[], chkpnt_s[], chkpnt_t[], chkpnt_u[], bits[], Machine_idr[],
  version[], resfile[], queuefile[], wdirectory[], WorktodoFile[];
extern const char program_name[];
extern FILE *glucasout;
extern UL qtest[];
extern char *qres[];
extern int qfl[];
extern const char *const inistring[];
extern char UserID[], UserPWD[], ComputerID[], ProxyHost[], speed[], hours[],
  DaysofWork[];
extern int CPUSpeed, CPUHours;
extern double Y_ALIGNED(16) tref[], Y_XBITS;
extern int Y_SELFTEST_ITERS;

#if defined(Y_USE_SSE2)
extern Y__M128D MM_bc[4], MM_auxt[4], MM_inv[4], MM_c, MM_Hsmall;
#endif /*Y_USE_SSE2*/


/* The rounding trick */
#include "round.h"
#define BE_TRICKY TRICKY_ROUND

extern UL Y_ITER, Y_EXPONENT, Y_NQ;
extern struct gtask Y_QD[];
extern UL Y_SBIT0, Y_SBIT, Y_SFORCE;
extern UL Y_SUMCHECK;
extern int Y_KILL, Y_RESIDUE_EQU_TWO;

/*************************************************************************
  the path of some control files
************************************************************************/
extern char inifile[SIZE_PATH], inputfile[SIZE_PATH], name[SIZE_PATH],
  /*pnetfile[SIZE_PATH],sdate[32],*/ gbuf[8 * SIZE_PATH];

/**************************************************************************
 To reduce memory comsumption and increase the speed, the arrays
two_to_phi[] and two_to_minusphi[] are stored only partialy.
**************************************************************************/
#define Y_UPDATE 32
#define Y_SHIFT_UPDATE 5

/*************************************************************************
 Some util macros    
*************************************************************************/
/* This is the character in directory tree: '\' or '/' */
#define DIR_SLASH "/"

/* to open files with some control */
#define fopen_file(_fd,_namefile,_option)                            \
    if(NULL==( _fd = fopen( _namefile, _option)))                    \
    {                                                                \
      perror(program_name);                                          \
      fprintf(stderr,"\t%s:%d\n", __FILE__, __LINE__);               \
      fprintf(stderr,"Unable to open %s file with option \"%s\"\n",  \
             _namefile,_option);                                     \
      exit(EXIT_FAILURE);                                            \
    }

/* to concatenate string without overflow the target string */
#define concatenate_string(_to,_from,_limit)                          \
  if((strlen( _to ) + strlen( _from ))<= _limit ) strcat( _to, _from);\
  else                                                                \
    {                                                                 \
      fprintf(stderr,"Glucas : overflow concatenating strings");      \
      exit(EXIT_FAILURE);                                             \
    }

/*
   this routine defines modular mul for 32 bits machine, to avoid 
   catastrophics overflows.
 
   Is not fast, but sure for exponents up to 2^30 
*/

#define modmul_32(_tar, _f1, _f2, _q) \
  {                                   \
    UL _ax = _f1, _ay = _f2;          \
    _tar = 0;                         \
    do                                \
      {                               \
	if(_ax & 1U) _tar += _ay;     \
	if(_tar >= _q) _tar -= _q;    \
	_ay += _ay;                   \
	if(_ay >= _q) _ay -= _q;      \
	_ax >>=1;                     \
      } while (_ax);                  \
  }
/*
#define modmul_32(_tar, _i, _b, _n)      \
{                                        \
  UL _j, aux = (_b);                     \
  _tar = 0;                              \
  for (_j = 1U; _j <= (_i); _j <<= 1)    \
  {                                      \
      if (_j & (_i))                     \
      {                                  \
	  _tar += aux;                   \
	  if (_tar >= (_n)) _tar -= (_n);\
      }                                  \
      aux <<= 1;                         \
      if (aux >= _n) aux -= (_n);        \
  }                                      \
}
*/



/*$Id$*/
