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
 *                                    bwt.h
 *
 * Dyadic Wavelet Transform.
 *
 * These routines are designed for computing the dyadic wavelet transform
 * and it's inverse transform using quadratic spline wavelet.
 *
 * To distinguish with the "dwt" (discrete wavelet transform), we call this
 * file as "bwt", but in fact, it should be dyadic wavelet transform.                                                                                 *
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef BWT_H
#define BWT_H


#include <vector.h>
#include <utilities.h>


namespace splab
{

    template<typename Type> Vector< Vector<Type> > bwt( const Vector<Type>&,
                                                        int );
    template<typename Type> Vector<Type> ibwt( const Vector< Vector<Type> >&,
                                               int );


    #include <bwt-impl.h>

}
// namespace splab


#endif
// BWT_H
