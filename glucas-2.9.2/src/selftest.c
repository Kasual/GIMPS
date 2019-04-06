/*$Id$*/
/*
    selftest.c. An collection of routines to make fast selftest
 
    Copyright (C) 2000-2006 Guillermo Ballester Valor, Klaus Kastens
 
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
/* This file includes some self-test routines */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "yeafft.h"
#include "gsetup.h"
#include "version.h"
#include "glucas.h"
#include "prefetch.h"

#if (defined(HAVE_SYS_UTSNAME_H) || !defined(HAVE_CONFIG_H)) && !defined(pccompiler) && !defined(macintosh)
# include <sys/utsname.h>
#elif defined(__MWERKS__) && (defined(macintosh) || defined(__INTEL__))
# include <utsname.h>
#endif

#define N_TESTS 31

UL qtest[N_TESTS]={ 2457601, 2752513, 2981887, 3276799, 4128767,
                    4816897, 5301301, 6029311, 6684673, 7798783,
                    9043967, 10223617, 11927551, 13631489, 16515073,
                    18874367, 20971521, 23592959, 26738687, 33030143,
                    36700159, 40370177, 46137343, 52428799, 60293119,
                    69730303, 81593293, 90123133, 108000043, 123000011,
                    145000033};

int qfl[N_TESTS] = {128, 144, 160, 192, 224, 256, 288, 320, 384, 448,
                    512, 576, 640, 768, 896, 1024, 1152, 1280, 1536, 1792,
                    2048, 2304, 2560, 3072, 3584, 4096, 4608, 5120, 6144, 7168,
                    8192};

char *qres[N_TESTS] = {"9D4603AF2919D905", "3C6D78B43A3C383A",
                       "F030A87BAE7E3979", "5B3B4C979ECF2E65",
                       "13A3AF2C0817B0E5", "FFE7F14D21B5E23C",
                       "ED892A998FBFBF65", "6B71D9119E4ED7CC",
                       "DEAECA0423F160C1", "378C4887FD13514D",
                       "60D5A9FE57A28B27", "0181C857137FF671",
                       "BDDAD6CB117EB5F6", "FCDD6D7816D0CDEF",
                       "54945CC36411836B", "339DB12CC64DFA7C",
                       "B93CF91BF334F0D6", "A0E97160575C10F3",
                       "21084004CDD83DDF", "1A479CDE7A34184C",
                       "3A368A8683127C76", "28A78C781B1F4E72",
                       "BD04190F60B26096", "E4477F3E938ADF11",
                       "FF945C59464BD00B", "59773B81752D65DC",
                       "3FEC67C5FA5789F8", "6197221020E383EA",
                       "5844DC56AF9D671D", "D6E014BCC0E96D34",
                       "2E0662F113E1871A"};

char *qres1[N_TESTS]= {"C0DD4583ECECDFAD", "EA1C78FD7D67E05A",
                       "2FCE731B3CB87BF2", "AABB54AE0CF01F01",
                       "9F95B93871BC3BEF", "CA6ACD6EE7E60EF1",
                       "BAA4BB15828564CD", "79364F07F7D16DF4",
                       "F6A616BB428EB94C", "940D6B7F2482B8EA",
                       "B296342A047E7D14", "11EA7223602C3445",
                       "5CA99F3FC2C66C89", "15FBD2AA1C8F932D",
                       "43F9EFA0880B2CA5", "3FC59A02F92C16F0",
                       "996AB2B97DEDA7B2", "CA1AAE5F2C4CDC30",
                       "4A14529456E2DF72", "787260E734015486",
                       "EB017F393EE04A34", "6C64D913016B0FDD",
                       "5DBF7116E20938A5", "CB22FF0389E65A9E",
                       "265326339B94C4CF", "944C9AEAC2D90D3D",
                       "4D99F2B0E4A3D4A8", "59EB7B9FC20E8B0D",
                       "12303E9E826AC73F", "4A635FCAEFF251B6",
                       "A5B21C2F2FCE7427"};

int Y_SELFTEST_ITERS;

int selftest(char *a)
{
  UL shift = 0;
  int ntest = 0, i, it;
  char name[SIZE_PATH];
  FILE *self;

  /* Set shift mode */
  if (Y_KILL)
    it = Test;
  else
    it = DoubleCheck;


  if(wdirectory[0] == '\0')
    strcpy(name, "selftest.que");
  else
    {
      strcpy (name, wdirectory);
      strcat (name, DIR_SLASH);
      strcat (name, "selftest.que");
    }
  fopen_file (self, name, "w");
  setvbuf (self, NULL, _IOLBF, BUFSIZ);
  switch (a[0])
    {
    case('a'): /* all the test */
            for(i = 0; i < N_TESTS; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    case('p'): /* "PrimeNet" from 256 to 2048 k*/
            for(i = 5; i < 21; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    case('s'): /* from 128 to 512 k*/
            for(i = 0; i < 11; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    case('b'): /* from 512 to 2048 k*/
            for(i = 10; i < 21; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    case('h'): /* from 2048 to 4096 k*/
            for(i = 20; i < 26; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    case('e'): /* from 4096 to 8192 k*/
            for(i = 25; i < N_TESTS; i++)
          {
            if(it)
              do
                {
                  shift = rand() % qtest[i];
                }
              while (shift == 0);
            fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
          }
      fclose (self);
      return 0;
    default:
      ntest = atoi (a);
      if(ntest >0 && ntest <= N_TESTS)
        {
          if(it)
            do
              {
                shift = rand() % qtest[ntest-1];
              }
            while (shift == 0);
          fprintf(self, "%ld %d %ld 0\n", qtest[ntest-1], it, shift);
          fclose (self);
          return ntest;
        }
      for(i = 0; i < N_TESTS; i++)
        {
          if(ntest == qfl[i])
            {
              if(it)
                do
                  {
                    shift = rand() % qtest[i];
                  }
                while (shift == 0);
              fprintf (self, "%ld %d %ld 0\n", qtest[i], it, shift);
              fclose (self);
              return (i+1);
            }
        }
      return -1;
    }
}

int self_inifile( void )
{
  FILE *f;
  char  name[SIZE_PATH];
  int iter_out = Y_SELFTEST_ITERS / 10;

  if(wdirectory[0]=='\0')
    strcpy(name,"selftest.ini");
  else
    {
      strcpy (name, wdirectory);
      strcat (name, DIR_SLASH);
      strcat (name, "selftest.ini");
    }


  if((f = fopen(name,"w")) != NULL)
    {
      setvbuf( f, NULL, _IOLBF, BUFSIZ );
      fprintf (f, "Verbose_flag=1\n");
      fprintf (f, "Alternative_output_flag=1\n");
      fprintf (f, "Iteration_output=%d\n", iter_out);
      fprintf (f, "File_output=selftest.out\n");
      fprintf (f, "Last_error_flag=1\n");
      fprintf (f, "Time_flag=1\n");
      fprintf (f, "Only_check_flag=1\n");
      fprintf (f, "Check_iterations=%d\n", Y_SELFTEST_ITERS);
      fprintf (f, "Roundoff_check=0\n");
      fprintf (f, "User_information=selftest\n");
      fprintf (f, "ComputerID=13466917\n");
      fclose (f);
      return  1;
    }
  else
    return 0;
}

int self_success( UL q, char *res)
{
  int i;

  for(i = 0; i < N_TESTS; i++)
    {
      if(q == qtest[i])
        {
          if(Y_SELFTEST_ITERS == 100)
            {
              if(strcmp (qres[i], res) == 0)
                return 1;
              else
                return 0;
            }
          else
            {
              if(strcmp (qres1[i], res) == 0)
                return 1;
              else
                return 0;
            }
        }
    }
  return 0;
}


int self_info(void)
{
  int i=1;
#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)

  int nt;
#endif
#if (defined(HAVE_SYS_UTSNAME_H) || !defined(HAVE_CONFIG_H)) && (!defined(pccompiler) || defined(__MWERKS__))

  struct utsname INFO;
#endif

  fprintf(glucasout,"\nSELFTEST INFORMATION\n");
  if(Alternative_output_flag == 2)
    printf("\nSELFTEST INFORMATION\n");
#if (defined(HAVE_SYS_UTSNAME_H) || !defined(HAVE_CONFIG_H)) && (!defined(pccompiler) || defined(__MWERKS__))

  if(uname( &INFO )!= -1)
    {
      fprintf(glucasout,"Host: %s.\nOS: %s. Release: %s. Version: %s\n",
              INFO.nodename,INFO.sysname,INFO.release,INFO.version);
      fprintf(glucasout,"Machine: %s \n",INFO.machine);
      if(Alternative_output_flag == 2)
        {
          printf("Host: %s.\nOS: %s. Release: %s. Version: %s\n",
                 INFO.nodename,INFO.sysname,INFO.release,INFO.version);
          printf("Machine: %s \n",INFO.machine);
        }
    }
#endif

  /* Print version and build */
  fprintf (glucasout, "%s %s-%s\n", program_name, version_string,
           build_string);
  if(Alternative_output_flag == 2)
    printf ("%s %s-%s\n", program_name, version_string, build_string);

  fprintf (glucasout,"-DY_AVAL=%d ",(int)Y_AVAL);
  fprintf (glucasout,"-DY_MEM_THRESHOLD=%d ",(int)Y_MEM_THRESHOLD);
  if(Alternative_output_flag == 2)
    {
      printf ("-DY_AVAL=%d ", (int)Y_AVAL);
      printf ("-DY_MEM_THRESHOLD=%d ", (int)Y_MEM_THRESHOLD);
    }
#ifdef Y_PADDING
  fprintf (glucasout, "-DY_BLOCKSIZE=%d -DY_SHIFT=%d", (int)Y_BLOCKSIZE,
           (int)Y_SHIFT);
  if (Alternative_output_flag == 2)
    printf("-DY_BLOCKSIZE=%d -DY_SHIFT=%d",(int)Y_BLOCKSIZE,(int)Y_SHIFT);
#endif

  fprintf(glucasout,"\n");
  fprintf(glucasout,"-DY_TARGET=%d ",(int)Y_TARGET);
  if (Alternative_output_flag == 2)
    {
      printf("\n");
      printf("-DY_TARGET=%d ",(int)Y_TARGET);
    }
#if defined(_OPENMP) || defined(_SUNMP) || defined(_PTHREADS)
  if(!Y_NTHREADS)
    {
# if defined(_OPENMP)
#  ifdef Y_NUM_THREADS
#   pragma omp parallel default(shared)

      {
        nt=Y_NUM_THREADS;
      }
#  else
#   pragma omp parallel default(shared)

      {
        nt = omp_get_num_procs();
      }
#  endif
# elif defined(_SUNMP)
      nt = Y_NUM_THREADS;
# else

      if(Y_NTHREADS == 0)
        nt=_PTHREADS;
# endif

    }
  else
    nt=Y_NTHREADS;
#endif
#ifdef _OPENMP

  fprintf(glucasout,"-D_OPENMP=%i ",nt);
  if(Alternative_output_flag==2)
    printf("-D_OPENMP=%i ",nt);
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef _SUNMP
  fprintf(glucasout,"-D_SUNMP=%i ",nt);
  if(Alternative_output_flag==2)
    printf("-D_SUNMP=%i ",nt);
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef _PTHREADS
  fprintf(glucasout,"-D_PTHREADS=%i ",nt);
  if(Alternative_output_flag==2)
    printf("-D_PTHREADS=%i ",nt);
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_MANY_REGISTERS
  fprintf(glucasout,"-DY_MANY_REGISTERS ");
  if(Alternative_output_flag==2)
    printf("-DY_MANY_REGISTERS ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_KILL_BRANCHES
  fprintf(glucasout,"-DY_KILL_BRANCHES ");
  if(Alternative_output_flag==2)
    printf("-DY_KILL_BRANCHES ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_USE_SSE2
  fprintf(glucasout,"-DY_USE_SSE2 ");
  if(Alternative_output_flag==2)
    printf("-DY_USE_SSE2 ");
  i++;
  if (i == 4 || i == 8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag == 2)
        printf("\n");
    }
#endif
#ifdef Y_MINIMUM
  fprintf(glucasout,"-DY_MINIMUM ");
  if(Alternative_output_flag==2)
    printf("-DY_MINIMUM ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_MAXIMUM
  fprintf(glucasout,"-DY_MAXIMUM ");
  if(Alternative_output_flag==2)
    printf("-DY_MAXIMUM ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_LONG_MACROS
  fprintf(glucasout,"-DY_LONG_MACROS ");
  if(Alternative_output_flag==2)
    printf("-DY_LONG_MACROS ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_VECTORIZE
  fprintf(glucasout,"-DY_VECTORIZE ");
  if(Alternative_output_flag==2)
    printf("-DY_VECTORIZE ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_VECTORIZE2
  fprintf(glucasout,"-DY_VECTORIZE2 ");
  if(Alternative_output_flag==2)
    printf("-DY_VECTORIZE2 ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_PREFETCH_EXPENSIVE
  fprintf(glucasout,"-DY_PREFETCH_EXPENSIVE ");
  if(Alternative_output_flag==2)
    printf("-DY_PREFETCH_EXPENSIVE ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#ifdef Y_ITANIUM
  fprintf(glucasout,"-DY_ITANIUM ");
  if(Alternative_output_flag==2)
    printf("-DY_ITANIUM ");
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
#if defined(Y_CACHE_LINE) && Y_TARGET!=0
  fprintf(glucasout,"-DY_CACHE_LINE=%d ", (int)Y_CACHE_LINE);
  if(Alternative_output_flag==2)
    printf("-DY_CACHE_LINE=%d ", (int)Y_CACHE_LINE);
  i++;
  if (i==4 || i==8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
#endif
  if (i!=4 && i!=8)
    {
      fprintf(glucasout,"\n");
      if(Alternative_output_flag==2)
        printf("\n");
    }
  return 0;
}


int nself( UL q)
{
  int i;
  for(i=0;i<N_TESTS;i++)
    if(q == qtest[i])
      return (i+1);
  return 0;
}

int nalter( UL q)
{
  int i;
  for(i=0;i<N_TESTS;i++)
    if(q == qtest[i])
      return qfl[i];
  return 0;
}

