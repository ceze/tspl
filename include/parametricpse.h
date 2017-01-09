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
 *                               parametricpse.h
 *
 * Parametric Power Spectrum Estimation Mothods.
 *
 * Techniques for spectrum estimation can generally be divided into parametric
 * (such as classical spectrum estimation) and non-parametric methods. The
 * parametric approaches assume that the underlying stationary stochastic
 * process has a certain structure which can be described using a small number
 * of parameters (for example, auto-regressive model). In these approaches,
 * the task is to estimate the parameters of the model that describes the
 * stochastic process.
 *
 * The widely used model is AR model, so this file provides three subroutines
 * to estimate the parameter of AR model, they are Yule-Walker method, Burg's
 * recursive mothod and forward-and-backward linear prediction least square
 * method.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef PARAMETRICPSE_H
#define PARAMETRICPSE_H


#include <vector.h>
#include <fft.h>
#include <correlation.h>
#include <levinson.h>
#include <linequs1.h>


namespace splab
{

    template<typename Type> Vector<Type> yulewalkerPSE( const Vector<Type>&,
                                                        int, Type& );
    template<typename Type> Vector<Type> burgPSE( const Vector<Type>&,
                                                  int, Type& );
    template<typename Type> Vector<Type> fblplsPSE( const Vector<Type>&,
                                                    int, Type& );

    template<typename Type> Vector<Type> armaPSD( const Vector<Type>&,
                                                  const Vector<Type>&,
                                                  const Type&, int );


    #include <parametricpse-impl.h>

}
// namespace splab


#endif
// PARAMETRICPSE_H
