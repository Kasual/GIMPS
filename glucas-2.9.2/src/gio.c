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
/*
 
   This file contains the input/output enhanced features for Glucas
   not included in MERS package.
   
   The most important information wich is not saved in intermediate files
   created by input() is UPDATE. To speed it up, Glucas, like MacLucasUNIX and
   other programs, stores only two_to_phi[] and two_to_minusphi[] factors 
   every UPDATE elements, so it reduces the needed memory and the access to it.
   
   UPDATE is asigned at runtime depending on the errors detected, and can vary 
   from 32 to 4. The file gnnn, where nnn is the exponent, has this 
   information.
 
*/
#include "gsetup.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#if defined(HAVE_LIMITS_H) || !defined(HAVE_CONFIG_H)
# include <limits.h>
/* Patch suggested by B.J.Beesley */
# if defined(SUN_V9_GCC) || defined(__sparc)
#  undef ULONG_MAX
#  define ULONG_MAX 0xFFFFFFFF
# endif
#endif

#if (defined(HAVE_SYS_TIME_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
# include <sys/time.h>
#elif defined(__APPLE__)
# include <sys/time.h>
#endif

#if (defined(HAVE_SYS_RESOURCE_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
# include <sys/resource.h>
#endif

#if defined(__MWERKS__) && defined(macintosh)
# include <Timer.h>
#endif

#if (defined(HAVE_UNISTD_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
# include <unistd.h>
#endif


/* include needed to use yeafft package */
#include "yeafft.h"
#include "ydebug.h"

#if ULONG_MAX > 0xFFFFFFFF
# define MAGIC_NUMBER_1 (UL)(0x7776757473727170)
# define MAGIC_NUMBER_2 (UL)(0x0706050403020100)
#define FMTHEX "%16.16lX"
#else
# define MAGIC_NUMBER_1 (UL)(0x73727170)
# define MAGIC_NUMBER_2 (UL)(0x03020100)
#define FMTHEX "%8.8lX"
#endif

/* glucas include */
#include "glucas.h"

void print_usage(void)
{
  printf("Usage: \n      Glucas [-options] [file1] [file2]...\n");
  printf("  Valid options :\n");
  printf("-c . Universal compatible style save file format.\n");
  printf("-d . Force random shift for initial value.\n");
  printf("-D <initial_shift>. Force initial shift to 'initial_shift'\n");
  printf("-f . The work in files are inserted in queue file\n");
  printf("     will be done first. The default is append.\n");
  printf("-i <inifile> . The ini file is 'inifile', otherwise is\n");
  printf("    'glucas.ini'.\n");
  printf("-h . Print this help output and exit\n");
  printf("-n . No nice mode. The program will run at normal priority.\n");
  printf("     The default is nice.\n");
  printf("-N <nice_adjustement>. Sets the priority. Default is 40.\n");
  printf("-o . Old style save file format. The default is new style.\n");
  printf("-q <queuefile> . The queue file is 'queuefile', otherwise is\n");
  printf("    'glucas.que'.\n");
  printf("-r <resultfile> .The result file is 'resultfile', otherwise\n");
  printf("    is 'result'.\n");
  printf("-s <n | small | big | huge | enormous | primenet | all>. Selftest.\n");
#ifdef _PTHREADS
  printf("-T <n>. Glucas will run using n threads.\n");
#endif
  printf("-S <n | small | big | huge | enormous | primenet | all>. Longer selftest.\n");
  printf("-v . Print version and exit\n");
  printf("-W <work_directory>. All files are read from /write to work_directory.\n");
  printf("-z . Force to use bigger than necessary FFT runlength.\n");
  printf("  ABOUT FILES: A file in arguments will be first taken as file.\n");
  printf("  If does not exist, then as an exponent to test\n");
  printf("  (In this case the work will be appended).\n");
  printf("  Files taken from PRIMENET, like worktodo.ini files, are\n");
  printf("  read and Lucas-Lehmer test work queued. Glucas does not\n");
  printf("  perform any factorization work.\n");
  printf("  YOU CAN SET MORE OPTIONS EDITING AN 'INI' FILE. See documentation.\n");
}

/*
 * Only used to verify exponent is prime, fast enough for p < 80M
 */
int primeq(unsigned long p)
{
  unsigned long j;

  if (p % 2 == 0 || p < 3)
    return (p == 2 ? 1 : 0);

  for (j = 3; (j * j) <= p; j += 2)
    if (p % j == 0)
      return 0;

  return 1;
}

int iszero(BIG_DOUBLE *x, UL N)
{
#ifdef TRICKY_ROUND
  BIG_DOUBLE A = bigA, B = bigB;
#endif

  UL j;

  for(j = 0; j < N; ++j)
    {
      if ( RINT(x[addr(j)]))
        return 0;
    }
  return 1;
}

/* See whether a Mersenne Prime is actually new */
int isnewprime(UL q)
{
  if( (q !=(UL)3 ) && (q != (UL)5 ) && (q != (UL)7) &&
      (q != (UL)13 ) && (q != (UL)17) && (q != (UL)19) &&
      (q != (UL)31 ) && (q != (UL)61) && (q != (UL)89) &&
      (q != (UL)107) && (q != (UL)127) && (q != (UL)521) &&
      (q != (UL)607) && (q != (UL)1279) &&
      (q != (UL)2203 ) && (q != (UL)2281 ) &&
      (q != (UL)3217 ) && (q != (UL)4253 ) &&
      (q != (UL)4423 ) && (q != (UL)9689 ) &&
      (q != (UL)9941 ) && (q != (UL)11213 ) &&
      (q != (UL)19937 ) && (q != (UL)21701 ) &&
      (q != (UL)23209 ) && (q != (UL)44497 ) &&
      (q != (UL)86243 ) && (q != (UL)110503 ) &&
      (q != (UL)132049 ) && (q != (UL)216091 ) &&
      (q != (UL)756839 ) && (q != (UL)859433 ) &&
      (q != (UL)1257787 ) && (q != (UL)1398269 ) &&
      (q != (UL)2976221 ) && (q != (UL)3021377 ) &&
      (q != (UL)6972593 ) && (q != (UL)13466917 ) &&
      (q != (UL)20996011 ) && (q != (UL)24036583 ) &&
      (q != (UL)25964951 ) && (q != (UL)30402457 ) &&
      (q != (UL)32582657 ) && (q != (UL)37156667 ) &&
      (q != (UL)42643801 ) && (q != (UL)43112609 ) &&
      (q != (UL)57885161 ) && (q != (UL)74207281 ))
    return 1;
  return 0;
}


void check_balanced(BIG_DOUBLE *x,UL N)
{
  UL   j,bj = 0,NminusOne = N-1;
  BIG_DOUBLE limit, hilim,lolim;

  hilim = high*0.5;
  lolim = low*0.5;
  for (j=0; j<NminusOne ; ++j)
    {
      if (bj >= c)
        {
          limit = hilim;
          bj-=c;
        }
      else
        {
          limit = lolim;
          bj+=b;
        }
      assert ((x[addr(j)]<=limit) && (x[addr(j)]>=-limit));
    }
  assert ((x[addr(NminusOne)]<=lolim) && (x[addr(NminusOne)]>=-lolim));

}

void balancetosdrep(  BIG_DOUBLE  *x, UL N)
{
  BIG_DOUBLE carry=0.0;
  UL   bj = N, j=0, NminusOne=N-1;
  while (j<NminusOne)
    {
      x[addr(j)]+=carry;
      if( bj >= c)
        {
          bj -= c;
          if(x[addr(j)] < 0.0)
            {
              carry = -1.0;
              x[addr(j)] += high;
            }
          else
            carry = 0.0;
        }
      else
        {
          bj += b;
          if(x[addr(j)] < 0.0)
            {
              carry = -1.0;
              x[addr(j)] += low;
            }
          else
            carry = 0.0;
        }
      j++;
    }

  x[addr(j)]+= carry;
  if(x[addr(j)] < 0.0)
    {
      carry = -1.0;
      x[addr(j)] += low;
    }
  else
    carry = 0.0;
  bj=N;

  j=0;
  while (carry != 0.0)
    {
      x[addr(j)]+=carry;
      if( bj >= c)
        {
          bj -= c;
          if(x[addr(j)] < 0.0)
            {
              carry = -1.0;
              x[addr(j)] += high;
            }
          else
            carry = 0.0;
        }
      else
        {
          bj += b;
          if(x[addr(j)] < 0.0)
            {
              carry = -1.0;
              x[addr(j)] += low;
            }
          else
            carry = 0.0;
        }
      j++;
    }
  return;
}

/*
   This version of res64 uses Y_SBIT to read the bits in a shifted 
   environment 
*/
int res64(BIG_DOUBLE *x, UL q, UL N , char *bits)
{
  UL i = 15L, j, k, k1 = 0L, word = 0L, totalbits = 64L, s, ib;

  /* the float array element where is the initial bit*/
  j = (UL)floor((double)Y_SBIT / Y_XBITS);
  /* the bit in the element */
  ib = (UL)floor((double)Y_SBIT - (double)j * Y_XBITS);
  do
    {
      k = (UL)( ceil((BIG_DOUBLE)(q) * (j + 1L) / N) -
                ceil((BIG_DOUBLE)(q) * j / N));
      k -= ib;
      if (k > totalbits)
        k = totalbits;
      totalbits -= k;
      word += (((UL)floor(x[addr(j)] + 0.5)) >> ib) << k1;
      ib = 0;
      j++;
      if(j == N)
        j = 0; /* See whether we've reached the end of float array */
      k += k1;
      while(k > 3L)
        {
          s= word & 15L;
          if(s > 9L)
            s+=7L;
          bits[i--] = (char)('0' + s);
          word >>= 4L;
          k -= 4L;
        }
      k1 = k;
    }
  while (totalbits);
  if(k1 > 0)
    bits[i--]= (char)('0' + (word & 15L));
  return (int)16;
}

/*
   This version of uses shift to read the bits in a shifted 
   environment, it returns unsigned chars (8 bits per unsigned).
   The caller routine must to assure there are enough chars 
   allocated
*/
UL residue_char(BIG_DOUBLE *x, UL q, UL N ,unsigned  char *bits, UL nbits, UL shift)
{
  UL i, j, k, k1=0L, word=0L, totalbits = nbits, s,ib;

  /* Number of unsigned chars to allocate minus one */
  i = (nbits - 1) >> 3;
  /* the float array element where is the initial bit */
  if( shift >= q)
    shift -= q;
  j = (UL)floor((double)shift / Y_XBITS);
  /* the bit in the element */
  ib = (UL)floor((double)shift - (double)j * Y_XBITS);
  do
    {
      k = (UL)( ceil((BIG_DOUBLE)(q) * (j + 1L) / N) - ceil((BIG_DOUBLE)(q) * j / N));
      k -= ib;
      if (k > totalbits)
        k = totalbits;
      totalbits -= k; /* now totalbits are the remaining bits to read */
      word += (((UL)floor(x[addr(j)] + 0.5)) >> ib) << k1;
      ib = 0; /* set ib to 0 after first read */
      j++; /* next float array index */
      if(j == N)
        j=0; /* See whether we've reached the end of float array */
      k += k1; /* now k are the bits in word to move to uchars */
      while(k > 7L)
        {
          s= word & 255U;
          bits[i--] = (unsigned char)(s);
          word >>= 8;
          k-=8;
        }
      k1 = k; /* k1 are the bits still no placed in uchar */
    }
  while (totalbits);
  /* move the final bits */
  if(k1 > 0)
    bits[i--]= (unsigned char)(word & ((1U << k) - 1U));
  return ((nbits-1) >> 3) + 1;
}


/*
   This routines write the gNNNNN file, which includes some information
   about ratio computing of DWT factors 
*/
int write_update(UL q)
{
  int i;
  FILE *iofich;
  char name[SIZE_PATH], num[16];
  if(wdirectory[0]=='\0')
    sprintf(name,"g%ld",q);
  else
    {
      sprintf (num, "%ld", q);
      strcpy (name, wdirectory);
      strcat (name, DIR_SLASH);
      strcat (name, num);
    }
  fopen_file (iofich, name, "w");
  setvbuf (iofich, NULL, _IOLBF , BUFSIZ);
  i=fprintf (iofich, "%ld\n%ld\n%ld\n%s\n%ld\n", q, UPDATE, SHIFT_UPDATE,
             Machine_idr, first_giter);
  fclose (iofich);
  return i;
}

/*
  Read the gNNNN file
  Returns 0 if no file, otherwise returns the number
  of readed integers
*/
int read_update(UL q)
{
  UL qr,u,su;
  int i;
  FILE *iofich;
  char name[SIZE_PATH];

  if(wdirectory[0] == '\0')
    sprintf (name, "g%ld", q);
  else
    sprintf(name, "%s/g%ld", wdirectory, q);

  iofich = fopen (name, "r");
  if(iofich == NULL)
    {
      UPDATE = Y_UPDATE;
      SHIFT_UPDATE = Y_SHIFT_UPDATE;
      return 0;
    }
  i=fscanf (iofich, "%ld %ld %ld %s %ld", &qr, &u, &su, Machine_idr,
            &first_giter);
  if(i >= 3 && q == qr )
    {
      UPDATE = u;
      SHIFT_UPDATE = su;
    }
  else
    {
      UPDATE = Y_UPDATE;
      SHIFT_UPDATE = Y_SHIFT_UPDATE;
    }
  fclose(iofich);
  if((UL)1 << SHIFT_UPDATE != UPDATE)
    {
      fprintf (stderr, "Inconsistent UPDATE/SHIFT_UPDATE\n");
      exit (EXIT_FAILURE);
    }
  return i;
}

/* to remove gNNNN file after finish l-l test */
int delete_update( UL q)
{
  char name[SIZE_PATH];

  if (wdirectory[0] == '\0')
    sprintf (name, "g%ld", q);
  else
    sprintf (name, "%s/g%ld", wdirectory, q);
  return remove
           (name);
}

/* This routine prints the date */
int print_date(void)
{
  time_t curtime;
  struct tm *loctime;
  char s[32],*p;
  curtime = time(NULL);
  loctime = localtime( &curtime);
  strcpy(s,"[");
  strcat(s, asctime (loctime));
  p = strchr (s, '\n');
  *p = '\0';
  strcat (s, "]\n");
  if (Alternative_output_flag == 2)
    fputs (s, stdout);
  return fputs (s, glucasout);
}

/*
   this routine prints in the time string as example
   "[Sun Jan 27 12:32:14 2002]\n" 
*/
time_t print_date2(char *s)
{
  time_t curtime;
  struct tm *loctime;
  char *p;

  curtime = time(NULL);
  loctime = localtime (&curtime);
  strcpy (s, "[");
  strcat (s, asctime(loctime));
  p = strchr(s, '\n');
  *p = '\0';
  strcat (s, "]\n");
  return curtime ;
}

/*
   this routine prints the abbreviate date,
   "[Jan 27 12:32:54]" 
*/
time_t print_date3(char *s)
{
  time_t curtime;
  struct tm *loctime;

  curtime = time(NULL);
  loctime = localtime( &curtime);
  strftime( s, 20, "[%b %d %H:%M:%S] ", loctime);
  return curtime ;
}

int write_glucasout(char *s, int alternative, int timestamp, int verbose)
{
  char buf[8 * SIZE_PATH];

  if(verbose)
    {
      if (timestamp == 1) /* complete date */
        {
          print_date2 (buf);
          strncat(buf, s, SIZE_PATH);
          s = buf;
        }
      else if (timestamp == 2) /* short date */
        {
          print_date3 (buf);
          strncat(buf, s, SIZE_PATH);
          s = buf;
        }
      fprintf(glucasout, "%s", s);
      fflush(glucasout);
      if(alternative == 2)
        printf("%s", s);
    }
  return 0;
}


int write_resultsfile(char *s, char *path, int timestamp)
{
  char buf[8 * SIZE_PATH];
  FILE *res;
  static time_t lasttime;
  time_t thistime;

  if (timestamp)
    {
      thistime = print_date2(buf);
      if (difftime(thistime, lasttime) > 60.0)
        {
          strncat(buf, s, SIZE_PATH);
          s = buf;
        }
      lasttime = thistime;
    }
  fopen_file(res, path, "a");
  setvbuf( res, NULL, _IOLBF, BUFSIZ );
  fprintf(res, "%s", s);
  fclose(res);
  return 0;
}


#if (!defined(HAVE_STRCHR) && defined(HAVE_CONFIG_H))
/* same as the BSD index() and the SYSV strchr() */
const char *strchr(const char *string, int ch)
{
  if (string == NULL)
    return(NULL);
  while (*string != '\0')
    if (*string == ch)
      return(string);
    else
      ++string;
  return(NULL);
}
#endif


#if (!defined(HAVE_UNISTD_H) && defined(HAVE_CONFIG_H)) || defined(pccompiler) || defined(macintosh)
/* This is an alternative to systems which have not getopt facilities */
/* This is a very simplified routine, please, avoid to use this for other
   proposes */

int optind = 1, goptindex=1;
char  *optarg;

int getopt( int gargc, char **gargv, char *optstring)
{
  int ch,l;
  char *aux;

  while(optind < gargc)
    {
      if (gargv[optind][0] =='-') /* is this an opt argument ?*/
        {
          l=strlen(gargv[optind]);
          ch=gargv[optind][goptindex];
          if((aux = strchr(optstring,ch)) != NULL)
            {
              if( aux[1]==':') /* An argument? */
                {
                  optind++;
                  if(optind < gargc )
                    optarg=gargv[optind++];
                  else
                    return -1;
                  goptindex=1;
                  return ch;
                }
              else
                optarg=NULL;
              goptindex++;
              if(goptindex == l)
                {
                  goptindex =1;
                  optind++;
                }
              return ch;
            }
          else
            {
              goptindex++;
              if(goptindex == l)
                {
                  goptindex =1;
                  optind++;
                }
              return '?';
            }
        }
      else
        return -1;
    }
  return -1;
}

#endif

/* This is the array of valid inputs in 'ini' files */
const char *const inistring[]=
  {
    "Verbose_flag", "Alternative_output_flag", "Iteration_output",
    "Output_file", "File_output", "Last_error_flag", "Time_flag",
    "Only_check_flag", "Check_iteration", "Check_iterations",
    "Save_iterations", "QA_interim_file", "ComputerID",
    "User_information", "Roundoff_check", "Use_primenet", "UserID",
    "UserPWD", "CPUSpeed", "CPUHours", "ProxyHost", "DaysofWork",
    "Machine_id", "User_info"
  };


/* returns the ordinal of valid string. Returns -1 if not match*/
/* Validate strings in 'ini' files */
int validate_inifile_string(char *p)
{
  unsigned int i;

  for (i = 0; i < sizeof(inistring) / sizeof(inistring[0]); i++)
    if (strcmp(inistring[i], p) == 0)
      return (int)i;
  return -1;
}

/* Read optional glucas.ini file */
int read_inifile( char *inifile )
{
  char cad[SIZE_PATH], *p, File_output[SIZE_PATH], name[SIZE_PATH],
  outbuf[SIZE_PATH];
  FILE *fich;
  int i, n = 0, ic;

  /* set the defaults */

  Verbose_flag = 0;
  Alternative_output_flag = 0;
  Iteration_output = 1000;
  Last_error_flag = 1;
  Time_flag = 1;
  Only_check_flag = 0;
  Check_iteration = 100;
  File_output[0] = '\0';
  save_iterations = 5000;
  QA_save = 0;
  /*Use_primenet=0;*/

  fich = fopen (inifile,"r");
  /* See whether the file exists, if not, Glucas automatically writes the default one */
  if (fich == NULL)
    {
      fopen_file (fich, inifile, "w");
      setvbuf (fich, NULL, _IOLBF, BUFSIZ );
#if defined(__MWERKS__) && defined(macintosh)
      /* Verbose_flag=1 + Alternative_output_flag=0 to open the console window */
      fprintf (fich, "Verbose_flag=1\n");
      fprintf (fich, "Last_error_flag=1\n");
      fprintf (fich, "Time_flag=1\n");
      fprintf (fich, "Alternative_output_flag=0\n");
      fprintf (fich, "Iteration_output=100\n");
#else

      fprintf (fich, "Verbose_flag=1\n");
      fprintf (fich, "Last_error_flag=1\n");
      fprintf (fich, "Time_flag=1\n");
      fprintf (fich, "Alternative_output_flag=1\n");
      fprintf (fich, "Iteration_output=20000\n");
#endif

      fprintf (fich, "Only_check_flag=0\n");
      fprintf (fich, "Check_iterations=100\n");
      fprintf (fich, "File_output=glucas.out\n");
      fprintf (fich, "Save_iterations=5000\n");
      fprintf (fich, "QA_interim_file=1000000\n");
      srand (time (NULL));
      fprintf (fich, "ComputerID=%d\n", rand ());
      fclose (fich);
    }
  else
    fclose (fich);

  fopen_file (fich, inifile, "r");
  i = fscanf (fich, "%s ", cad);
  do
    {
      p = strchr(cad, '=');
      *p = '\0';
      if( (ic = validate_inifile_string (cad)) < 0)
        fprintf (stderr, "Warning: Unknown entry %s in glucas.ini\n", cad);
      switch (ic)
        {
        case 0: /* Verbose flag */
          Verbose_flag = atoi (p + 1);
          break;
        case 1: /* Alternative output flag */
          Alternative_output_flag = atoi (p + 1);
          break;
        case 2: /* Iteration output */
          Iteration_output = atoi (p + 1);
          break;
        case 3: /* "Output_file" error in older documentation */
        case 4: /* File output */
          sprintf(File_output,"%s", (p + 1));
          break;
        case 5: /* Last error flag */
          Last_error_flag = atoi (p + 1);
          break;
        case 6: /* Time Flag */
          Time_flag = atoi (p + 1);
          break;
        case 7: /* Only check flag */
          Only_check_flag = atoi (p + 1);
          break;
        case 8: /* "Check_iteration" backwards compatibility */
        case 9: /* Check iterations */
          Check_iteration = atoi (p + 1);
          break;
        case 10: /* Save iterations */
          save_iterations = (UL) atoi (p + 1);
          break;
        case 11: /* Interim files */
          QA_save = atoi (p + 1);
          break;
        case 12: /* Computer id */
          strcpy (ComputerID, p + 1);
          break;
        case 13: /* User information */
          sprintf (guser, "%s", p + 1);
          break;
        case 14: /* Roundoff check flag */
          ecf = atoi (p + 1);
          ecfn = ecf;
          break;
        case 16: /* Primenet User id */
          strcpy (UserID, p + 1);
          break;
          /*case 15:  Use primenet
          Use_primenet=atoi(p+1);
          break;*/
        case 21: /* Days of work */
          strcpy (DaysofWork, p + 1);
          break;
        default:
          break;
        }
      i=fscanf (fich, "%s ",cad);
    }
  while (i > 0);
  fclose(fich);

  if(DaysofWork[0] == '\0')
    strcpy (DaysofWork, "30");

  if(glucasout != NULL && Alternative_output_flag != 0)
    fclose(glucasout);

  if(Alternative_output_flag)
    {
      if(!strlen (File_output))
        if(wdirectory[0] == '\0')
          {
            fopen_file (glucasout, "glucas.out", "a");
          }
        else
          {
            strcpy (name, wdirectory);
            strcat (name, DIR_SLASH);
            strcat (name, "glucas.out");
            fopen_file (glucasout, name, "a");
          }
      else
        if(wdirectory[0] == '\0')
          {
            fopen_file (glucasout, File_output, "a");
          }
        else
          {
            strcpy (name, wdirectory);
            strcat (name, DIR_SLASH);
            strcat (name, File_output);
            fopen_file (glucasout, name, "a");
          }
      setvbuf (glucasout, NULL, _IOLBF, BUFSIZ );
    }
  else
    glucasout = stdout;

  /*Write a machine Identifier if not readed */
  if (ComputerID[0] == '\0')
    {
      fich = fopen (inifile, "a");
      if( fich == NULL)
        {
          sprintf (outbuf, "Warning: Glucas can't write the Machine Id. in %s file.\n",
                   inifile);
          write_glucasout (outbuf, Alternative_output_flag, 0, 1);
          return n;
        }
      srand (time (NULL));
      sprintf (ComputerID, "%d", rand ());
      fprintf (fich, "\nComputerID=%s\n", ComputerID);
      sprintf (outbuf, "Appended ComputerID=%s in file %s.\n", ComputerID, inifile);
      write_glucasout (outbuf, Alternative_output_flag, 0, Verbose_flag);
      fclose (fich);
    }
  return n;
}

/* Return seconds elapsed since last call */
double sec_elapsed(void)
{
  double delta;
#if (defined(HAVE_GETHRTIME) || !defined(HAVE_CONFIG_H)) && defined(__sun)
  /*
    On some versions of SunOS 'gettimeofday()' returns incorrect values
    in the tv_usec field. Therefore we use 'gethrtime()' instead.
  */
  static hrtime_t start;
  hrtime_t end;

  end = gethrtime();
  delta = 1E-9 * (end - start);
  start = end;
#elif (defined(HAVE_GETTIMEOFDAY) || !defined(HAVE_CONFIG_H)) && !defined(vms) && !defined(pccompiler) && !defined(macintosh)

  static struct timeval start;
  struct timeval end;

  if (gettimeofday(&end,NULL) == 0)
    {
      /*  delta = timeval_subtract(&end, &start);  */
      delta = end.tv_sec + 1E-6 * end.tv_usec - start.tv_sec - 1E-6 * start.tv_usec;
      start = end;
    }
  else
    delta = 0;
#elif defined(__MWERKS__) && defined(macintosh)

  static double start;
  double end;
  struct UnsignedWide tmp;

  Microseconds(&tmp);
  end = 1E-6 * (4294967296.0 * tmp.hi + tmp.lo);
  delta = end - start;
  start = end;
#else

  static time_t start;
  time_t end;

  end = time(NULL);
  if (end != (time_t)-1)
    {
      delta = end - start;
      start = end;
    }
  else
    delta = 0;
#endif

  return delta;
}

/* Return seconds of user time elapsed from last call */
/* If getrusage is not defined it returns 0.0         */
double sec_elapsed_user(void)
{
  double user_time=0.0;
#if (defined(HAVE_GETRUSAGE) || !defined(HAVE_CONFIG_H)) && !defined(vms) && !defined(pccompiler) && !defined(macintosh)

  static double begin_user_time;
  double end_user_time;
  struct rusage ru;

  getrusage(RUSAGE_SELF, &ru);
  end_user_time = (BIG_DOUBLE)(ru.ru_utime.tv_sec) +
                  1e-6*(BIG_DOUBLE)(ru.ru_utime.tv_usec);
#if defined(_PTHREADS) && defined(linux)
  /* this computes the time used by the child threads */
  {
    int i;
    for(i = 0; i < Y_NTHREADS; i++)
      end_user_time += y_user[i];
  }
#endif
  user_time = end_user_time - begin_user_time;
  begin_user_time = end_user_time;
#endif

  return user_time;
}

/* this is a modification of W. Edgington code */
size_t inf_time (UL j, char *buf, size_t pos)
{
  double elapsed;
#if (defined(HAVE_GETRUSAGE) || !defined(HAVE_CONFIG_H)) && !defined(vms) && !defined(pccompiler) && !defined(macintosh)

  double user_time, percent_cpu;
#else

  double hh, mm, ss;
#endif

  elapsed = sec_elapsed();
#if (defined(HAVE_GETRUSAGE) || !defined(HAVE_CONFIG_H)) && !defined(vms) && !defined(pccompiler) && !defined(macintosh)

  user_time = sec_elapsed_user();
  /*
    On most systems (Linux/{ppc,sunos,x86}, SunOS) timings from getrusage only
    have 10ms precision (1 jiffie = 10ms). This can result into strange (>100%)
    values for %CPU for small user_time values.
  */
#if !defined(_OPENMP) && !defined(_SUNMP) && !defined(_PTHREADS)

  if (user_time < elapsed)
    percent_cpu = 100 * user_time/elapsed;
  else
    percent_cpu = 100;
#else

  percent_cpu = 100 * user_time/elapsed;
#endif

  pos += sprintf(buf + pos, "%.2f user ", user_time);
  pos += sprintf(buf + pos, "%3.0f%% CPU ", percent_cpu);
#else

  ss = 60 * modf(60 * modf(elapsed/3600, &hh), &mm);
  if (hh > 0)
    pos += sprintf(buf + pos, "%.0f:%02.0f:%02.0f real ", hh, mm, ss);
  else
    pos += sprintf(buf + pos, "%.0f:%02.0f real ", mm, ss);
#endif

  if((elapsed > 0) && ( j > last_inf))
    pos += sprintf(buf + pos, "(%.4f sec/iter).", elapsed/(j-last_inf));
  last_inf=j;
  return pos;
}


/* This routine prints some information */
int iteration_inf( UL j, UL last)
{
  double completed;
  size_t pos = 0;
  char buf[128], tmp[20];

  completed = ((double)(j) * 100.0) / last;
  pos += sprintf(buf + pos, "Iter. %*lu (%6.2f%%), ", sprintf(tmp, "%lu", last), j, completed);
  if (Last_error_flag)
    pos += sprintf(buf + pos, "Err= %.3f, ",Err);
  if (Time_flag)
    pos = inf_time(j, buf, pos);
  pos += sprintf(buf + pos, "\n");
  if ( Selftest_flag)
    write_glucasout(buf, Alternative_output_flag, 0, 1);
  else
    write_glucasout(buf, Alternative_output_flag, 2, 1);
  return 0;
}

/* This routine returns the number of exponents readed from queuefile */
/* It stores the values since position 'n' in array qd[] */
int read_queue(struct gtask *qd, int n)
{
  FILE *queue;
  UL qaux;
  int i = 0, j;
  UL aux[5];
  char caux[128], *ca;

  /*
     first check the version of queue file format, if less than four items
     or second item is bigger than 10 then is version 1, otherwise
     is version 2
  */
  queue = fopen(queuefile,"r");
  if(queue == NULL)
    return 0;

  ca = fgets (caux, 128, queue);
  fclose (queue);
  if ((j = sscanf( caux, "%ld %ld %ld %ld", &aux[0],
                   &aux[1], &aux[2], &aux[3])) < 4)
    { /* Old format */
      queue = fopen (queuefile, "r");
      if(queue == NULL)
        return 0;
      i = 0;
      while (fscanf (queue, "%ld ", &qaux) != EOF)
        {
          if((n + i) < Y_MAX_QD)
            {
              qd[n + i].exponent = qaux;
              qd[n + i].task = 0;
              qd[n + i].parm1 = 0;
              qd[n + i].parm2 = 0;
              i++;
              if (qaux < 64)
                {
                  fprintf (stderr, "Glucas: Error. Bad exponent in queue file %s.\n", queuefile);
                  exit (EXIT_FAILURE);
                }
            }
          else
            fprintf(stderr,"Warning: Too much exponents in queue file, exponent %ld will not be checked\n",qaux);
        }
      fclose(queue);
    }
  else if (j == 4)
    { /* new format */
      queue = fopen (queuefile, "r");
      if(queue == NULL)
        return 0;
      i = 0;
      while ((n + i) < Y_MAX_QD && fgets(caux, 128, queue) != NULL)
        {
          if (sscanf(caux, "%ld %d %ld %ld ", &qd[n + i].exponent,
                     &qd[n + i].task, &qd[n + i].parm1, &qd[n + i].parm2) != 4)
            { /* malformed line */
              fprintf (stderr, "Glucas: Error. Wrong line in queue file %s.\n", queuefile);
              exit (EXIT_FAILURE);
            }

          if (qd[n + i].exponent < 64)
            {
              fprintf (stderr, "Glucas: Error. Bad exponent in queue file %s.\n", queuefile);
              exit (EXIT_FAILURE);
            }

          /* See whether the shift bit is greater than exponent */
          if (qd[n + i].exponent <= qd[n +i].parm1)
            qd[n + i].parm1 %= qd[n + i].exponent;

          /* Check the task */
          if ((qd[n + i].task != Test) && (qd[n + i].task != DoubleCheck))
            {
              fprintf (stderr, "Glucas: Error. Unknown task in queue file %s.\n", queuefile);
              exit (EXIT_FAILURE);
            }
          i++;
        }
      while (fgets(caux, 128, queue) != NULL)
        {
          sscanf (caux, "%ld %ld %ld %ld ", &aux[1], &aux[2], &aux[3], &aux[4]);
          fprintf (stderr, "Warning: Too much exponents in queue file, exponent %ld will be removed from queue\n", aux[1]);
          i++;
        }
      fclose (queue);
    }
  else
    { /* wrong format */
      fprintf (stderr, "Glucas: Error. Unknown format for queue file %s.\n", queuefile);
      exit (EXIT_FAILURE);
    }
  if(i <= Y_MAX_QD)
    return i;
  else
    return Y_MAX_QD;
}

/*
   This routine writes the exponents stored in struct array qd[] to queuefile 
*/
int write_queue(struct gtask *qd, int n)
{
  FILE *queue;
  int i;

  fopen_file(queue,queuefile,"w");
  setvbuf( queue, NULL, _IOLBF, BUFSIZ );

  for(i = 0; i < n; i++)
    fprintf (queue, "%ld %d %ld %ld\n", qd[i].exponent, qd[i].task,
             qd[i].parm1, qd[i].parm2);
  fclose (queue);
  return i;
}

/* This routine deletes first line in queuefile */
int delete_first_job(struct gtask *qd)
{
  int j, i;
  i = read_queue (qd , 0);
  if(i > 1)
    {
      for(j = 0; j < (i - 1); j++)
        {
          qd[j].exponent = qd[j + 1].exponent;
          qd[j].task = qd[j + 1].task;
          qd[j].parm1 = qd[j + 1].parm1;
          qd[j].parm2 = qd[j + 1].parm2;
        }
      write_queue (qd, i - 1);
      qd[i - 1].exponent = 0;
      qd[i - 1].task = 0;
      qd[i - 1].parm1 = 0;
      qd[i - 1].parm2 = 0;
    }
  if(i == 1)
    remove
      (queuefile);
  Y_NQ = i - 1;
  return (i - 1);
}

/*
   This routine generates an initial shift bit if necessary and stores 
   it in qd.parm1
*/
void generate_shift0(struct gtask *gt)
{
  if(gt->task == Test)
    {
      gt->parm1 = 0;
      return;
    }
  if(Y_SFORCE)
    gt->parm1 = (Y_SFORCE % gt->exponent);
  else
    {
      do
        {
          gt->parm1 = rand() % gt->exponent;
        }
      while (gt->parm1 == 0);
    }
}

/*
   This manage the queue file. It returns the number of exponents in queue 
   (and in qd[] array) after adding the job (if any) from inputfile 
*/
int manage_queue( int mode, char *inputfile, struct gtask *qd)
{
  FILE *fich;
  char *p1, *p2, cad[60], line[60], name[SIZE_PATH];
  UL q = 0;
  int i, j, k = 0,qv = 0,it = 0;

  i = read_queue (qd, 0);
  Y_NQ = i;
  if(inputfile[0] == '\0')
    return i;
  if(wdirectory[0] == '\0')
    fich = fopen (inputfile,"r");
  else
    {
      strcpy (name, wdirectory);
      strcat (name, DIR_SLASH);
      strcat (name, inputfile);
      fich = fopen(name, "r");
    }
  if (fich == NULL)
    {
      /* very inmediate file?. Better than the name of a file is the
      exponent itself */
      sscanf (inputfile, "%ld", &q);
      {
        /* be sure the exponent is not already in queue */
        if (q == 0 || q > 156000000)
          qv = 0;
        else
          qv = 1;
        for (j = 0; j < i; j++)
          if(q == qd[j].exponent)
            qv = 0;
        if ( qv )
          {
            if(i < Y_MAX_QD)
              {
                qd[i].exponent = q; /* append */
                if(Y_SBIT0)
                  qd[i].task = DoubleCheck;/* make a double check*/
                else
                  qd[i].task = Test; /* L-L test */
                generate_shift0 (&qd[i]);
                qd[i].parm2 = 0;
                i++;
                write_queue (qd, i);
              }
            else
              fprintf (stderr, "Warning: Too much exponents in queue, exponent %ld not queued\n", q);
            Y_NQ = i;
            return i;
          }
      }
      Y_NQ = i;
      return i;
    }
  j = fscanf (fich, "%s ", line);
  do
    {
      strcpy (cad, line);
      p1 = strchr(cad, '=');
      if(p1 == NULL)
        q = (UL) atoi (cad);
      else
        {
          *p1 = '\0';
          p2 = strchr(p1 + 1,',');
          if(p2 != NULL)
            *p2 = '\0';
          switch (cad[0])
            {
            case 'T':
              q = (UL) atoi (p1 + 1);
              it = Test;
              break;
            case 'D':
              q = (UL) atoi (p1 + 1);
              it = DoubleCheck;
              break;
            default:
              fprintf (stderr, "Warning: Unknown entry %s in %s\n", line, inputfile);
            }
        }
      if (q == 0 || q >156000000)
        qv = 0;
      else
        qv = 1;
      /* It is a different q? */
      for (j = 0; j < i; j++)
        if(q == qd[j].exponent)
          qv =0;
      if(qv)
        {
          if(i < Y_MAX_QD)
            {
              if(mode == 0)
                {
                  qd[i].exponent = q; /* append */
                  qd[i].task = it;
                  generate_shift0 (&qd[i]);
                  qd[i].parm2 = 0;
                  i++;
                }
              else
                {
                  /* insert */
                  for(j = i; j > k; j--)
                    {
                      qd[j].exponent = qd[j - 1].exponent;
                      qd[j].task = qd[j - 1].task;
                      qd[j].parm1 = qd[j - 1].parm1;
                      qd[j].parm2 = qd[j - 1].parm2;
                    }
                  i++;
                  qd[k].exponent = q;
                  qd[k].task = it;
                  generate_shift0 (&qd[k]);
                  qd[k].parm2 = 0;
                  k++;
                }
            }
          else
            fprintf (stderr, "Warning: Too much exponents in queue, exponent %ld not queued\n", q);
        }
    }
  while (fscanf (fich, "%s ", line) > 0);
  write_queue (qd, i);
  fclose (fich);
  Y_NQ = i;
  return i;
}

/*
   This is the routine to control read/write from/to save files 
   It calls read_check_point() and write_check_point(). 
*/
int ginput(int nq, struct gtask *qd, UL *q, UL *n, UL *j, BIG_DOUBLE *err,
           BIG_DOUBLE **x)
{
  UL qr, nr, nrr;

  if(nq == 0 )
    return 0;
  qr = qd[0].exponent;
  nr = read_check_point( qr, q, n, j, err, x);

  switch(nr)
    {
    case 0:
      /* Try to read t-file, if exists */
      if(wdirectory[0] =='\0')
        {
          sprintf (chkpnt_s, "s%ld", qr);
          sprintf (chkpnt_t, "t%ld", qr);
        }
      else
        {
          sprintf (chkpnt_s, "%s/s%ld", wdirectory, qr);
          sprintf (chkpnt_t, "%s/t%ld", wdirectory, qr);
        }
      if(rename (chkpnt_t, chkpnt_s) != -1)
        {
          nrr = read_check_point( qr, q, n, j, err, x);
          switch(nrr)
            {
            case 0:/* no t-file */
              if(wdirectory[0] == '\0')
                {
                  sprintf (chkpnt_s, "s%ld", *q);
                  sprintf (chkpnt_t, "t%ld", *q);
                  sprintf (chkpnt_u, "u%ld", *q);
                }
              else
                {
                  sprintf (chkpnt_s, "%s/s%ld", wdirectory, *q);
                  sprintf (chkpnt_t, "%s/t%ld", wdirectory, *q);
                  sprintf (chkpnt_u, "%s/u%ld", wdirectory, *q);
                }
              return 3;
            case 1:/* t-file readed with success */
              write_glucasout ("No s-save file. Recovered from t-file\n", Alternative_output_flag, 0, 1);
              return 2;
            case -1:/* Also problems with t-file */
              rename (chkpnt_s, chkpnt_t);
              return 1;
            default:
              fprintf(stderr, "ginput: Unspected read_check_point() return value\n");
              exit (EXIT_FAILURE);
            }
        }
      else
        {
          /* No s-file nor t-file */
          if(wdirectory[0] == '\0')
            {
              sprintf (chkpnt_s, "s%ld", *q);
              sprintf (chkpnt_t, "t%ld", *q);
              sprintf (chkpnt_u, "u%ld", *q);
            }
          else
            {
              sprintf (chkpnt_s, "%s/s%ld", wdirectory, *q);
              sprintf (chkpnt_t, "%s/t%ld", wdirectory, *q);
              sprintf (chkpnt_u, "%s/u%ld", wdirectory, *q);
            }
          return 3;
        }
    case 1:
      return 2;
    case -1:
      /* Try to read t-file, if exists */
      if(wdirectory[0] == '\0')
        {
          sprintf (chkpnt_s, "s%ld", qr);
          sprintf (chkpnt_t, "t%ld", qr);
          sprintf (chkpnt_u, "u%ld", qr);
        }
      else
        {
          sprintf (chkpnt_s, "%s/s%ld", wdirectory, qr);
          sprintf (chkpnt_t, "%s/t%ld", wdirectory, qr);
          sprintf (chkpnt_t, "%s/u%ld", wdirectory, qr);
        }
      rename (chkpnt_s, chkpnt_u);
      if(rename (chkpnt_t, chkpnt_s) != -1)
        {
          nrr = read_check_point( qr, q, n, j, err, x);
          switch(nrr)
            {
            case 0:/* no t-file */
              rename (chkpnt_s, chkpnt_t);
              rename (chkpnt_u, chkpnt_s);
              return 1;
            case 1:/* t-file readed with success */
              write_glucasout ("Error in save file. Recovered from t-file\n", Alternative_output_flag, 0, 1);
              remove
                (chkpnt_u);
              return 2;
            case -1:/* Also problems with t-file */
              rename (chkpnt_s, chkpnt_t);
              rename (chkpnt_u, chkpnt_s);
              return 1;
            }
        }
      else
        {
          rename (chkpnt_u, chkpnt_s);
          return 1;
        }
    default:
      fprintf(stderr, "ginput: Unspected read_check_point() return value\n");
      exit (EXIT_FAILURE);
    }
  return 1;
}

/* It returns 0 if no checkpoint is present, so we need to start the exponent for iteration 1 */
/* It returns 1 if success, it returns -1 if any problem */
/* qr ir the request q.  n == 0 in entry point then the n will be
   the stored in the file, other way *n will remain unaltered */

int read_check_point(UL qr, UL *q, UL *n, UL *j, BIG_DOUBLE *err,
                     BIG_DOUBLE **x)
{
  int retval, mode=0; /* mode=0 is for binary-read  storage */

  /*fver is the mformat version. Currently 1 or 2 */
  UL tmp, tmp1, *ultmp=NULL, sumcheck=(UL)0, fver;
  FILE *chkpnt_fp;


  if(wdirectory[0] == '\0')
    sprintf (chkpnt_s, "s%ld", qr);
  else
    sprintf(chkpnt_s, "%s/s%ld", wdirectory,qr);

  if ((chkpnt_fp = fopen (chkpnt_s, "rb")) == NULL)
    {
      /* no checkpoint to read, we need to start */
      *q = qr;
      *err = 0.0;
      *j = 0;
      return 0;
    }

  if ((retval = fread (&tmp, sizeof(tmp), 1, chkpnt_fp)) == sizeof(tmp)
      || retval == 1 )
    {
      switch (tmp)
        {
        case MAGIC_NUMBER_1:
          mode = 0; /* binary/real mode */
          break;
        case MAGIC_NUMBER_2:
          mode = 1; /* binary/integer mode */
          break;
        default:
          fclose (chkpnt_fp);
          /* try now the Interchangeable format */
          fopen_file (chkpnt_fp, chkpnt_s, "rb");
          fread_UL (&tmp, chkpnt_fp);
          if ((tmp & 0xFFFFFFFF) == 0x006A64B1)
            {
              mode = 2;
            }
          else
            {
              fprintf (stderr, "%s: Unknown format file.\n", program_name);
              fclose (chkpnt_fp);
              return (-1);
            }
        }
    }
  else
    {
      fprintf (stderr, "%s: Unable to read first block.\n", program_name);
      fclose (chkpnt_fp);
      return (-1);
    }
  add_sumcheck (&sumcheck, tmp);

  if (mode == 2)
    {
#if ULONG_MAX == 0XFFFFFFFF
      /* read the format version block  */
      fread_UL (&fver, chkpnt_fp);
      if((fver & 0x000000FF) > 2)
        {
          fprintf (stderr,"%s: Unspected format file version.\n", program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
      add_sumcheck (&sumcheck, fver);
#else

      if((tmp & 0x000000FF00000000) > 0x200000000)
        {
          fprintf (stderr, "%s: Unspected format file version.\n", program_name);
          fclose(chkpnt_fp);
          return (-1);
        }
      else
        fver= (tmp >> 32) & 255U;
#endif


      /*Read the program/version/task block */
      if ((retval = fread_UL (&tmp, chkpnt_fp)) == sizeof(tmp)
          || retval == 1 )
        add_sumcheck(&sumcheck, tmp);
      else
        {
          fprintf (stderr, "%s: Unable to read second block.\n", program_name);
          fclose(chkpnt_fp);
          return (-1);
        }

      /*
      Patch to older than v2.7a version for Macs. This bug was found by Klaus Kastens. 
      For gcc-compiler, open a file with "r" is the same than "rb". For other compilers
      is not. We used "r" for prior versions we now will use "rb" and "wb" 
      */

#if ULONG_MAX ==0xFFFFFFFF
      if((tmp & 0x00FFFF00) <= 0x00020700)
        {
          fclose (chkpnt_fp);
          sumcheck = 0;
          fopen_file (chkpnt_fp, chkpnt_s, "r");
          fread_UL (&tmp, chkpnt_fp);
          add_sumcheck (&sumcheck, tmp);
          fread_UL (&tmp1, chkpnt_fp);
          add_sumcheck (&sumcheck, tmp1);
          fread_UL (&tmp, chkpnt_fp);
          add_sumcheck (&sumcheck, tmp);
        }
#else
      if((tmp & 0x0000000000FFFF00) <= 0x0000000000020700)
        {
          fclose(chkpnt_fp);
          sumcheck = 0;
          fopen_file (chkpnt_fp,chkpnt_s, "r");
          fread_UL (&tmp, chkpnt_fp);
          add_sumcheck (&sumcheck, tmp);
          fread_UL (&tmp, chkpnt_fp);
          add_sumcheck (&sumcheck, tmp);
        }
#endif


#if ULONG_MAX == 0XFFFFFFFF
      /* read the second part of header. The task */
      fread_UL (&tmp1, chkpnt_fp);

      if (tmp1 & 0x000000FF)
        {
          fprintf (stderr, "%s: Unknown task. At the moment Glucas only make Lucas-Lehmer tests.\n",
                   program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
      add_sumcheck (&sumcheck, tmp1);
#else

      if(tmp & 0x000000FF00000000)
        {
          fprintf (stderr, "%s: Unknown task. At the moment Glucas only make Lucas-Lehmer tests.\n",
                   program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
#endif

    }

  /*
     version 2 is compatible with version 1, i.e., it can read version 1 formats
     without any problem. In this case Y_SBIT0=0 
  */
  if( mode == 2)
    retval = fread_twin (q, &Y_SBIT0, chkpnt_fp);
  else
    retval = fread (q, sizeof(*q), 1, chkpnt_fp);
  if( retval != 1 && retval != sizeof(*q))
    {
      fprintf (stderr, "%s : Cannot read q\n", program_name);
      fclose (chkpnt_fp);
      return (-1);
    }
  if (*q != qr)
    {
      fprintf (stderr, "%s : q readed = %ld is not the same than request %ld\n",
               program_name, *q, qr);
      fclose (chkpnt_fp);
      return (-1);
    }

  add_sumcheck (&sumcheck, *q);
  add_sumcheck (&sumcheck, Y_SBIT0);

  if( mode == 2)
    retval = fread_UL_pad (&tmp1, chkpnt_fp);
  else
    retval = fread (&tmp1, sizeof(tmp1), 1, chkpnt_fp);
  if( retval != 1 && retval != sizeof(tmp1))
    {
      fprintf (stderr,"%s : Cannot read n\n", program_name);
      fclose (chkpnt_fp);
      return (-1);
    }
  add_sumcheck (&sumcheck, tmp1);

  if(mode == 0)
    *n = tmp1;
  else if( (*n) == 0)
    *n = tmp1;

#if defined(__MWERKS__)
#pragma warn_possunwant off
#endif

  for (tmp1 = *n; (tmp1 & (UL)1) == 0; tmp1 >>= 1)
    ;
#if defined(__MWERKS__)
#pragma warn_possunwant reset
#endif

  if  ( *n > (((UL)1 << 31) + 1L) ||
        (tmp1 != (UL)1 && tmp1 != (UL)3 && tmp1 != (UL)5 && tmp1 != (UL)7 && tmp1 != (UL)9))
    {
      fprintf(stderr,"%s : Bad FFT run length\n",program_name);
      fclose(chkpnt_fp);
      return(-1);
    }

  if( mode == 2)
    retval = fread_UL_pad (j, chkpnt_fp);
  else
    retval = fread (j, sizeof(*j), 1, chkpnt_fp);
  if( retval != 1 && retval != sizeof(*j))
    {
      fprintf (stderr, "%s : Cannot read j\n", program_name);
      fclose (chkpnt_fp);
      return (-1);
    }
  if( *j > (*q - (UL)2))
    {
      fprintf (stderr, "%s : Invalid iteration\n", program_name);
      return (-1);
    }

  add_sumcheck (&sumcheck, *j);

  if ( mode == 0)
    {
      if ((retval = fread(err, sizeof(*err), 1, chkpnt_fp)) != sizeof(*err)
          && retval != 1 )
        {
          fprintf (stderr, "%s : Cannot read err\n", program_name);
          fclose(chkpnt_fp);
          return (-1);
        }
    }
  else /* read the error stored in integer form */
    {
      if(mode == 2)
        retval = fread_UL_pad (&tmp1, chkpnt_fp);
      else
        retval = fread(&tmp1, sizeof(tmp1), 1, chkpnt_fp);
      if( retval != 1 && retval != sizeof(tmp1))
        {
          fprintf (stderr, "%s : Cannot read err\n", program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
      *err = (BIG_DOUBLE)0.000001 * tmp1;
    }

  add_sumcheck(&sumcheck, tmp1);

  if( *err < 0.0 || *err > 0.49)
    {
      fprintf (stderr, "%s : Invalid readed err\n", program_name);
      fclose (chkpnt_fp);
      return (-1);
    }

  if ((PX) != NULL)
    free((char *)(PX));
  if ((PX = CALLOC_DOUBLES((addr((*n))))) == NULL)
    {
      fprintf (stderr, "%s : Cannot allocate space for data\n", program_name);
      fclose (chkpnt_fp);
      return (-2);
    }
  else
    (*x) = ALIGN_DOUBLES(PX);

  if (mode  ==  0)
    {
      if ((retval = (int)fread (*x, sizeof(x[0][0]), (size_t)*n, chkpnt_fp)) !=
          (int)((*n)*sizeof(x[0][0])) && retval != (int)*n )
        {
          fprintf (stderr, "%s : Cannot read data in real mode\n", program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
      for (tmp = *n - 1; (addr(tmp) - (tmp)) != 0; tmp--)
        x[0][addr(tmp)] = x[0][tmp];
    }
  else
    {
      if ((ultmp = (UL *)calloc((size_t)((*q - (UL)1) / BITS_PER_UL + 1),
                                sizeof(UL))) == NULL)
        {
          fprintf (stderr, "%s : Cannot allocate space for temporal integer data\n",
                   program_name);
          fclose (chkpnt_fp);
          return (-2);
        }
      tmp1 = (*q - (UL)1)/(UL)BITS_PER_UL + (UL)1;
      if( mode == 2)
        {
          retval = fread_UL_array (ultmp, (size_t)tmp1, chkpnt_fp);
#if ULONG_MAX > 0xFFFFFFFF

          retval -= tmp1;
#else

          retval = (retval << 1) - tmp1 - (tmp1 & (UL)1);
#endif

          if( retval != 0)
            {
              fprintf (stderr, "%s : Cannot read data in integer mode. Read %d items. Expected %ld\n",
                       program_name, retval, tmp1);
            }
        }
      else if (( retval = (int)fread (ultmp, sizeof(UL), (size_t)tmp1, chkpnt_fp)) !=
               (int)(tmp1 * sizeof(UL)) && retval != (int)tmp1)
        {
          fprintf (stderr, "%s : Cannot read data in integer mode. Read %d items. Expected %ld\n",
                   program_name, retval, tmp1);
        }
      if( (tmp1 = (UL) read_standar (*q, *n, x[0], ultmp, &sumcheck)) != *n )
        {
          fprintf (stderr,"%s : Error reading integer data. Readed %ld items\n",
                   program_name, tmp1);
          fclose (chkpnt_fp);
          return (-1);
        }
      free((char *) ultmp);

      /* read last carry if any*/
      if( mode == 2)
        retval = fread_UL_pad (&tmp1, chkpnt_fp);
      else
        retval = fread (&tmp1, sizeof(tmp1), 1, chkpnt_fp);
      if( retval != 1 && retval != sizeof(tmp1))
        {
          fprintf (stderr,"%s : Error reading last carry.\n", program_name);
          fclose (chkpnt_fp);
          return (-1);
        }
      add_sumcheck (&sumcheck , tmp1);
      /* sumcheck += tmp1; */
      if( mode == 2)
        retval = fread_UL_pad (&tmp1, chkpnt_fp);
      else
        retval = fread (&tmp1, sizeof(tmp1), 1, chkpnt_fp);
      if( retval != 1 && retval != sizeof(tmp1))
        {
          fprintf (stderr, "%s : Error reading sumcheck.\n", program_name);
          fclose (chkpnt_fp);
          return (-1);
        }

#if ULONG_MAX > 0xFFFFFFFF
      sumcheck = sumcheck32 ( sumcheck );
#endif

      if(sumcheck != tmp1)
        {
          fprintf (stderr, "%s : Sumckeck error failed!. Read " FMTHEX " and computed " FMTHEX "\n",
                   program_name, tmp1, sumcheck);
          fclose(chkpnt_fp);
          return (-1);
        }
    }
  fclose(chkpnt_fp);
  if(wdirectory[0] == '\0')
    {
      sprintf (chkpnt_s, "s%ld", *q);
      sprintf (chkpnt_t, "t%ld", *q);
      sprintf (chkpnt_u, "u%ld", *q);
    }
  else
    {
      sprintf (chkpnt_s, "%s/s%ld", wdirectory, *q);
      sprintf (chkpnt_t, "%s/t%ld", wdirectory, *q);
      sprintf (chkpnt_u, "%s/u%ld", wdirectory, *q);
    }
  return(1);
}


/* THIS ROUTINE HAS A LOT OF LINES FROM WILL EDGINGTON MERS's */
/* print what we've done so far to a (local) checkpoint file */
/* If mode==0  then it writes the real/binary form: it is faster but
       files are biger and machine dependent 
   If mode ==1 then it writes the integer/binary form: it is slow but
       files are smaller and the format 'almost' universal
   If mode ==2 then it writes the universal integer/binary form.
  if All is OK. then it returns 1, else it returns 0 */

int write_check_point(UL q, UL n, UL j, BIG_DOUBLE err, BIG_DOUBLE *x, int mode)
{
  FILE *chkpnt_fp;
  /* magic number ala the file(1) command; note that it's byte order dependent, */
  /* which we want so that we don't try to use a save file from a different */
  /* architecture or byte order */
  UL k, tmp1=0, *ultmp=NULL, sumcheck=(UL)0, sumcheck0=(UL)0;
  int retval; /* return value from fwrite(); Linux man page and libc.a do not agree on correct value */

  /* refuse to write checkpoint files for small exponents */
  /* Anyway, YEAFFT, does not run in such little exponents */
  if (mode == 0)
    k=MAGIC_NUMBER_1;
  else if (mode == 1)
    k=MAGIC_NUMBER_2;
  /*
    0x00000033 will be Glucas firm,
    reserved 0x00000011 for mprime
    and 0x00000022 for Mlucas
    The three higher bytes are reserved for version information
  */
  else
    {
      k = (UL) 0x00020933; /* Version 2.9 . Program : Glucas */
      tmp1 = (UL) 0x006A64B1; /* The signature 6972593 */
    }

  if (q < (UL)1000)
    return(1);
  if ((chkpnt_fp = fopen(chkpnt_s, "r")) != NULL)
    {
      (void)fclose(chkpnt_fp);
      /* this line is a diff of W. Edgington. Prior to delete t-file, rename it to u-file */
      /* when all is OK, then it will delete u-file */
      (void)rename(chkpnt_t, chkpnt_u);
      (void)rename(chkpnt_s, chkpnt_t);
    }
  if ((chkpnt_fp = fopen(chkpnt_s, "wb")) == NULL)
    {
      perror(program_name);
      fprintf(stderr, "%s: could not open '%s' for writing\n", program_name, chkpnt_s);
      return(0);
    }

  /* 'tag' the file as coming from this program, computer architecture, etc. */
  /* For the version of format, we will set version 1 when Y_SBIT0==0, so it will be backward
     compatible for all version prior to Glucas 2.9.0 releases */
  if (mode == 2)
    {
#if ULONG_MAX > 0xFFFFFFFF
      if(Y_SBIT0==0)
        tmp1 += 0x100000000U; /* version 1 of format */
      else
        tmp1 += 0x200000000U; /* version 2 of format */
      retval = fwrite_UL(tmp1,chkpnt_fp);
#else

      retval =fwrite_UL(tmp1,chkpnt_fp);
      add_sumcheck (&sumcheck, tmp1);
      if(Y_SBIT0==0)
        tmp1= 1U;
      else
        tmp1 = 2U;
      retval =fwrite_UL(tmp1,chkpnt_fp);
#endif

      if (retval != sizeof(tmp1) && retval != 1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot format checkpoint info (k = %ld)\n", program_name, k);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          return(0);
        }
      add_sumcheck (&sumcheck, tmp1);
    }

  if (mode == 2)
    retval = fwrite_UL_pad(k, chkpnt_fp);
  else
    retval = fwrite(&k, sizeof(k), 1, chkpnt_fp);
  if (retval != sizeof(k) && retval != 1)
    {
      perror(program_name);
      fprintf(stderr, "%s: cannot write checkpoint info (k = %ld)\n", program_name, k);
      fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
      return(0);
    }
  add_sumcheck (&sumcheck, k);

  if (mode == 2)
    retval = fwrite_twin(q, Y_SBIT0, chkpnt_fp);
  else
    retval = fwrite(&q, sizeof(q), 1, chkpnt_fp);
  if (retval != sizeof(q) && retval != 1)
    {
      perror(program_name);
      fprintf(stderr, "%s: cannot write checkpoint info (q = %ld)\n", program_name, q);
      fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
      return(0);
    }
  add_sumcheck (&sumcheck, q);
  add_sumcheck (&sumcheck, Y_SBIT0);

  if (mode == 2)
    retval = fwrite_UL_pad(n, chkpnt_fp);
  else
    retval = fwrite(&n, sizeof(n), 1, chkpnt_fp);
  if (retval != sizeof(n) && retval != 1)
    {
      perror(program_name);
      fprintf(stderr, "%s: cannot write checkpoint info (n = %ld , q = %ld )\n",
              program_name, n, q);
      fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
      return(0);
    }
  add_sumcheck (&sumcheck, n);

  if (mode == 2)
    retval = fwrite_UL_pad(j, chkpnt_fp);
  else
    retval = fwrite(&j, sizeof(j), 1, chkpnt_fp);
  if (retval != sizeof(j) && retval != 1)
    {
      perror(program_name);
      fprintf(stderr, "%s: cannot write checkpoint info (j = %ld , q = %ld )\n",
              program_name, j, q);
      fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
      return(0);
    }
  add_sumcheck (&sumcheck, j);

  /* if mode is 1 skip to write err in float form*/
  if(mode == 0)
    {
      if ((retval = fwrite(&err, sizeof(err), 1, chkpnt_fp)) != sizeof(err) && retval != 1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write checkpoint info (err, q = %ld )\n", program_name, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          return(0);
        }
    }
  else
    {
      tmp1=err*(BIG_DOUBLE)1000000.0;
      if (mode == 2)
        retval = fwrite_UL_pad(tmp1, chkpnt_fp);
      else
        retval = fwrite(&tmp1, sizeof(tmp1), 1, chkpnt_fp);
      if (retval != sizeof(tmp1) && retval != 1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write checkpoint info (err, q = %ld )\n", program_name, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          return(0);
        }
      add_sumcheck (&sumcheck, tmp1);
    }

  if(mode == 0)
    {
      for (k = 0; k < n; k++)
        if ((retval = fwrite(&(x[addr(k)]), sizeof(x[0]), 1, chkpnt_fp)) != sizeof(x[0]) && retval != 1)
          {
            perror(program_name);
            fprintf(stderr, "%s: cannot write checkpoint info (x, q = %ld )\n", program_name, q);
            fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
            return(0);
          }
    }
  else
    {
      if (((ultmp) = (UL *)calloc((size_t)((q-(UL)1)/BITS_PER_UL + (UL)1),
                                  sizeof(UL))) == NULL)
        {
          fprintf(stderr,"%s : Cannot allocate space for temporal integer data\n",program_name);
          return (0);
        }
      if( (tmp1 = (UL)write_standar( q, n, x, ultmp, &sumcheck0)) != ((q-(UL)1)/(UL)BITS_PER_UL + (UL)1))
        {
          fprintf(stderr,"%s : Error writing integer data. Written %ld items of %ld\n",
                  program_name,tmp1,(q-(UL)1)/(UL)BITS_PER_UL+ (UL)1);
          free ((char *) ultmp);
          return (-1);
        }
      Y_SUMCHECK=sumcheck0;
      add_sumcheck(&sumcheck, sumcheck0);
      if( mode == 2)
        {
          retval = fwrite_UL_array(ultmp, (size_t)tmp1, chkpnt_fp);
#if ULONG_MAX > 0xFFFFFFFF

          retval -= tmp1;
#else

          retval = (retval<<1) - tmp1 - (tmp1 & (UL)1);
#endif

          if (retval !=0)
            {
              perror(program_name);
              fprintf(stderr, "%s: cannot write integer checkpoint info (x, q = %ld )\n", program_name, q);
              fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
              free ((char *) ultmp);
              return(0);
            }
        }
      else if ((UL)(retval = fwrite(&(ultmp[0]), sizeof(UL),(size_t) tmp1, chkpnt_fp)) != (UL)(sizeof(UL)*tmp1)
               && (UL)retval != tmp1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write integer checkpoint info (x, q = %ld )\n", program_name, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          free (ultmp);
          return(0);
        }

      tmp1 = (UL)0; /* Last carry, 0 in this version */
      if (mode == 2)
        retval = fwrite_UL_pad(tmp1, chkpnt_fp);
      else
        retval = fwrite(&tmp1, sizeof(tmp1), 1, chkpnt_fp);
      if (retval != sizeof(tmp1) && retval != 1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write last carry (j = %ld , q = %ld )\n",
                  program_name, j, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          return(0);
        }
      free ((char *) ultmp);
#if ULONG_MAX > 0xFFFFFFFF
      /* pass to 32 bits nine proof */
      sumcheck = sumcheck32( sumcheck );
#endif

      if (mode == 2)
        retval = fwrite_UL_pad(sumcheck, chkpnt_fp);
      else
        retval = fwrite(&sumcheck, sizeof(sumcheck), 1, chkpnt_fp);
      if (retval != sizeof(tmp1) && retval != 1)
        {
          perror(program_name);
          fprintf(stderr, "%s: cannot write sumcheck (j = %ld , q = %ld )\n",
                  program_name, j, q);
          fprintf(stderr, "%s: fwrite returned %d\n", program_name, retval);
          return(0);
        }
    }

  fflush(chkpnt_fp);
#ifdef HAVE_UNISTD_H

  fsync(fileno(chkpnt_fp));
#endif
  /* now delete the u-file */
  remove
    (chkpnt_u);
  return(!fclose(chkpnt_fp));
}

/* This routine copy the s-save file to interim file */
int write_interim_file(UL iter)
{
  char interim_file[32], buf[8];
  FILE *inp, *outp;
  size_t i;

  interim_file[0]='\0';
  sprintf(interim_file,"%s_%8.8ld",chkpnt_s,iter);
  fopen_file(inp,chkpnt_s,"rb");
  fopen_file(outp,interim_file,"wb");

  do
    {
      i = fread(buf,sizeof(char),8,inp);
      if(i !=0)
        fwrite(buf,sizeof(char),i,outp);
    }
  while (i != 0);
  fclose(inp);
  fflush(outp);
#ifdef HAVE_UNISTD_H

  fsync(fileno(outp));
#endif

  fclose(outp);
  return 1;
}

/*
   This routine copy the s-save file to interim file
   It is used as a Test proposes in case of the discovery
   of a New Prime!
*/
int copy_last_save_file( void )
{
  char interim_file[32], buf[8];
  FILE *inp, *outp;
  size_t i;

  interim_file[0]='\0';
  sprintf(interim_file,"%s_last",chkpnt_s);
  fopen_file(inp,chkpnt_s,"rb");
  fopen_file(outp,interim_file,"wb");

  do
    {
      i = fread(buf,sizeof(char),8,inp);
      if(i !=0)
        fwrite(buf,sizeof(char),i,outp);
    }
  while (i != 0);
  fclose(inp);
  fflush(outp);
#ifdef HAVE_UNISTD_H

  fsync(fileno(outp));
#endif

  fclose(outp);
  return 1;
}

/*$Id$*/




