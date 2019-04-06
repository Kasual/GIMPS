/*$Id$*/
/*  This file is part of
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
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
/* This include the files with macros depending on Y_TARGET */
#include "prefetch.h"

#ifdef Y_TARGET
# if (Y_TARGET == 0) || (Y_TARGET == 16) || (Y_TARGET == 17)
#   include "ygeneric.h"
# elif Y_TARGET == 1
#   include "yx86.h"
# elif Y_TARGET == 11
#   include "yx86p3.h"
# elif Y_TARGET == 12
#   include "yx86ath.h"
# elif Y_TARGET == 13
#   include "yx86p4.h"
# elif (Y_TARGET == 21) || (Y_TARGET == 23) || (Y_TARGET == 24)
#   include "ygeneric.h"
# elif (Y_TARGET == 32) || (Y_TARGET == 31)
#   include "ygeneric.h"
# elif Y_TARGET == 41
#   include "ygeneric.h"
# elif Y_TARGET == 51
#   include "ygeneric.h"
# elif Y_TARGET == 61
#   include "ygeneric.h"
# elif Y_TARGET == 99
#   include "ygeneric.h"
# endif
#else
# include "ygeneric.h"
#endif

/*$Id$*/















