/*$Id$*/
/*
    YEAFFT. A library to make real convolutions using Fast Fourier
	    Transforms. 
    Copyright (C) 2000-2006  Guillermo Ballester Valor
 
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
/*  */
#define cplx_trig_min_load(_t0,_p) \
    _t0##r = *(_p);\
    _t0##i = *(_p + 1);

#define cplx_trig_min_square(_t,_s)\
    _t##r = (_s##r - _s##i)*(_s##r + _s##i);\
    _t##i = 2.0 * _s##r * _s##i ;

#define cplx_trig_min_mul(_t,_s0,_s1)\
    _t##r = (_s0##r * _s1##r) - (_s0##i * _s1##i);\
    _t##i = (_s0##r * _s1##i) + (_s0##i * _s1##r);

#define cplx_trig_min_div(_t,_s0,_s1)\
    _t##r = (_s0##r * _s1##r) + (_s0##i * _s1##i);\
    _t##i = (_s0##i * _s1##r) - (_s0##r * _s1##i);

#define cplx_trig_min_muldiv(_tm,_td,_s0,_s1)\
    _tm##r = (_s0##r * _s1##r) - (_s0##i * _s1##i);\
    _td##r = (_s0##r * _s1##r) + (_s0##i * _s1##i);\
    _tm##i = (_s0##r * _s1##i) + (_s0##i * _s1##r);\
    _td##i = (_s0##i * _s1##r) - (_s0##r * _s1##i);



