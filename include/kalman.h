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
 *                                   kalman.h
 *
 * Kalman Filter.
 *
 * The Kalman filter is an efficient recursive filter that estimates the
 * internal state of a linear dynamic system from a series of noisy
 * measurements. In most applications, the internal state is much larger
 * (more degrees of freedom) than the few "observable" parameters which are
 * measured. However, by combining a series of measurements, the Kalman
 * filter can estimate the entire internal state.
 *
 * A wide variety of Kalman filters have now been developed, from Kalman's
 * original formulation, now called the simple Kalman filter, the Kalman-Bucy
 * filter, Schmidt's extended filter, the information filter, and a variety
 * of square-root filters that were developed by Bierman, Thornton and so on.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef KALMAN_H
#define KALMAN_H


#include <vector.h>
#include <matrix.h>
#include <inverse.h>


namespace splab
{

    template<typename Type>
    Vector<Type> kalman( const Matrix<Type>&, const Matrix<Type>&,
                         const Matrix<Type>&, const Matrix<Type>&,
                         const Vector<Type>&, const Vector<Type>&,
                         const Vector<Type>& );


    #include <kalman-impl.h>

}
// namespace splab


#endif
// KALMAN_H
