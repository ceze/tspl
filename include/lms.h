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
 *                                   lms.h
 *
 * Least Mean Square Adaptive Filter.
 *
 * Least mean squares (LMS) algorithms are a class of adaptive filter used
 * to mimic a desired filter by finding the filter coefficients that relate
 * to producing the least mean squares of the error signal (difference
 * between the desired and the actual signal). It is a stochastic gradient
 * descent method in that the filter is only adapted based on the error at
 * the current time.
 *
 * This file implement three types of the LMS algorithm: conventional LMS,
 * algorithm, LMS-Newton algorhm and normalized LMS algorithm.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef LMS_H
#define LMS_H


#include <vector.h>
#include <matrix.h>


namespace splab
{

    template<typename Type>
    Type lms( const Type&, const Type&, Vector<Type>&, const Type& );

    template<typename Type>
    Type lmsNewton( const Type&, const Type&, Vector<Type>&,
                    const Type&, const Type&, const Type& );

    template<typename Type>
    Type lmsNormalize( const Type&, const Type&, Vector<Type>&,
                       const Type&, const Type& );


    #include <lms-impl.h>

}
// namespace splab


#endif
// LMS_H
