/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                                 constants.h
 *
 * Some constants often used in numeric computing.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef	CONSTANTS_H
#define CONSTANTS_H


namespace splab
{

    const double	EPS		    = 2.220446049250313e-016;

    const double	PI		    = 3.141592653589793;
    const double	TWOPI	    = 6.283185307179586;
    const double	HALFPI	    = 1.570796326794897;
    const double	D2R 	    = 0.017453292519943;

    const double	EXP		    = 2.718281828459046;

    const double	RT2		    = 1.414213562373095;
    const double	IRT2	    = 0.707106781186547;

    const int		FORWARD	    = 1;
    const int		INVERSE	    = 0;

    const int		MAXTERM	    = 20;
    const int       INITSIZE    = 20;
    const int       EXTFACTOR   = 2;

}
// namespace splab


#endif
// CONSTANTS_H
