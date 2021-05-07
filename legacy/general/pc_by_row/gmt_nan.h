/*--------------------------------------------------------------------
 *   The GMT-system:	gmt_nan.h	[Automatically Generated]
 *
 *   Copyright (c) 1991-2004 by P. Wessel and W. H. F. Smith
 *   See COPYING file for copying and redistribution conditions.
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 of the License.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * Machine-dependent macros for generation of NaNs on this system
 * This file was created by the program gmt_nan_init.c
 */

#ifndef _GMT_NAN_H
#define _GMT_NAN_H
#ifdef NO_IEEE
#define GMT_make_fnan(x) (x = MAX_FLT)
#else
#define GMT_make_fnan(x) (((unsigned int *) &x)[0] = 0x7fffffff)
#endif
#ifdef NO_IEEE
#define GMT_make_dnan(x) (x = MAX_DBL)
#else
#define GMT_make_dnan(x) (((unsigned int *) &x)[0] = 0xffffffff, ((unsigned int *) &x)[1] = 0x7fffffff)

#endif
/* Conditionally define the last resort (poor man's isnan)
   in case they are not supplied in gmt_math.h */

#ifndef GMT_is_fnan
#define GMT_is_fnan(x) ((x) != (x))
#endif
#ifndef GMT_is_dnan
#define GMT_is_dnan(x) ((x) != (x))
#endif

#endif /* _GMT_NAN_H */
