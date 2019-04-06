/* $Id$ */
/*
  (c) 2000-2005  Guillermo Ballester Valor

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

  This file includes the code to control the Lucas-Lehmer iteration
 
*/
#include <stdlib.h>
#include <stdio.h>
#include "yeafft.h"
#include "glucas.h"
#include "gsetup.h"


/****************************************************************************
 *           Lucas Test - specific routines                                 *
 ***************************************************************************/


/*
  The Lucas square routine. it returns the number of iterations actually
  done, i.e. 0 or 1. When returns -1 then an error has occurred in writing a
  save file. See main routine comments to the meaning of lmode arg. 
*/
int lucas_square(UL q, BIG_DOUBLE *x, UL N, UL iter, int lmode,
                 int error_log)
{
  char outbuf[SIZE_PATH];

  HACK_ALIGN_STACK_EVEN();

  switch (lmode)
    {
    case -1:
      /*
      First. The first iteration in the loop.
      It is not multithreaded
      */
      first_normalize (x, N << 1);
      lucas_first_auto_convolution (x, N);
      return 1;

    case 1:
      /*
      Last. The last iteration in the loop.
      */
#if defined(_PTHREADS)

      Y_ERR_FLAG = error_log;

      /* this tell threads they have to finish after backward pass 1 */
      Y_THREADS_TERMINATE = 1;

      /* Start the loop ! */
      y_barrier ((Y_NTHREADS + 1), &gbarrier_start);

      /* Wait after finish the loop ! */
      y_barrier ((Y_NTHREADS + 1), &gbarrier_finished);

      /* cancel the threads */
      y_join_and_freemem_all_threads();

#else /* _PTHREADS */

      lucas_dit_carry_norm_dif (x, N << 1, error_log);
      lucas_auto_convolution (x, N);

#endif /* _PTHREADS */

      y_fftb_pass1_lucas_last (x, N);
      last_normalize (x, N << 1, error_log);
      terminated = 1;
      return 1;

    case 9:
      /* Save a  residue inmediately without doing an iteration*/
      y_fftb_pass1_lucas_last (x, N);
      last_normalize (x, N << 1, error_log);
      terminate = 0;
      terminated = 1;
      return 0;

    case 2:
      /*Save. Save a residue an do a normal iteration.*/
      y_fftb_pass1_lucas_last (x, N);
      last_normalize (x, N << 1, error_log);

      if(write_check_point (q, N << 1, iter, Err, x, smode) <= 0)
        return -1;
      res64 (x, q, N << 1, bits);

      /* Check whether we have 0 residues. Something is wrong */
      if ( (strcmp(bits, "0000000000000000") == 0) && iszero ( x, N << 1))
        {
          sprintf( outbuf, "M%ld. Unespected zero residues at iteration %ld. Exiting", q,
                   iter);
          write_glucasout (outbuf, Alternative_output_flag, 2, Verbose_flag);
          exit(EXIT_FAILURE);
        }

      /* is a suspicious 2 residue? */
      if ( strcmp(bits, "0000000000000002") == 0)
        {
          sprintf( outbuf, "M%ld. Suspicious residue at iteration %ld", q,
                   iter);
          write_glucasout (outbuf, Alternative_output_flag, 2, Verbose_flag);
          if (Y_RESIDUE_EQU_TWO)
            {
              sprintf( outbuf, "M%ld. Is a reincident case. Exiting ", q);
              write_glucasout (outbuf, Alternative_output_flag, 2, Verbose_flag);
            }
          else
            Y_RESIDUE_EQU_TWO = 1;
        }
      else
        Y_RESIDUE_EQU_TWO = 0;

      /* Read inifile, is a hot read */
      read_inifile (inifile);

      if((QA_save != 0) && (iter % QA_save == (UL)0))
        {
          write_interim_file (iter);
          sprintf (outbuf, "M%ld. Saved Interim file at iteration %ld. Res64: %s.\n",
                   q, iter, bits);
          write_glucasout (outbuf, Alternative_output_flag, 2, Verbose_flag);
          write_resultsfile (outbuf, resfile, 1);
        }
      else
        {
          sprintf (outbuf, "M%ld. Saved file at iteration %ld. Res64: %s.\n", q, iter, bits);
          write_glucasout (outbuf, Alternative_output_flag, 2, Verbose_flag);
        }
      first_normalize (x, N << 1);
      lucas_first_auto_convolution (x, N);
      return 1;

    default:
      /*Normal.  A normal iteration.*/
#if defined(_PTHREADS)

      Y_ERR_FLAG = error_log;

      /* Start the loop ! */
      y_barrier ((Y_NTHREADS + 1), &gbarrier_start);

      /* Wait after finish the loop ! */
      y_barrier ((Y_NTHREADS + 1), &gbarrier_finished);
#else

      lucas_dit_carry_norm_dif (x, N << 1, error_log);

      lucas_auto_convolution (x, N);
#endif

      return 1;
    }
}
