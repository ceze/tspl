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
 *                                 nleroot.h
 *
 * Solution for nonlinear equation.
 *
 * This file includes three routines for finding root of nonlinear equation.
 * The bisection method don't need compute the function's gradient, but it
 * has a slow convergence rate, linear convergence. Newton method has a
 * quadratic convergence, however, it has to compute the gradient of function.
 * The secant method don't need compute the gradient, in adition, it has a
 * superlinear convergence rate( the order is 1.618 ), so it's a practical
 * method in many cases.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef NLEROOT_H
#define NLEROOT_H


#include <cmath>
#include <constants.h>
#include <nlfunc.h>


using namespace std;


namespace splab
{

    template<typename Type> Type bisection( NLFunc<Type> &f, Type a, Type b,
                                            Type tol=Type(1.0e-6) );

    template<typename Type> Type newton( NLFunc<Type> &f, Type x0,
                                         const Type tol=Type(1.0e-6),
                                         int maxItr=MAXTERM );

    template<typename Type> Type secant( NLFunc<Type> &f, Type x1, Type x2,
                                         const Type tol=Type(1.0e-6),
                                         int maxItr=MAXTERM );


    #include <nleroot-impl.h>

}
// namespace splab


#endif
// NLEROOT_H
