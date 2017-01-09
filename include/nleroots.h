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
 *                                 nleroots.h
 *
 * Solution for nonlinear equations.
 *
 * This file includes two routines for finding roots of nonlinear equations.
 * The Seidel method is and improvement of Fixed-Point iteration mithod,
 * which don't need compute the Jacobi matrix, but with a slow convergence
 * rate. The other is Newton method, which can provide a faster rate of
 * convergence, but at the const of computing Jacobi matrix and its inverse
 * to get the iterative increment.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef NLEROOTS_H
#define NLEROOTS_H


#include <nlfuncs.h>
#include <linequs1.h>


namespace splab
{

    template<typename Type>
    Vector<Type> seidel( NLEqus<Type> &G, const Vector<Type> &X0,
                         const Type tol=Type(1.0e-6),
                         int maxItr=MAXTERM );

    template<typename Type>
    Vector<Type> newton( NLFuncs<Type> &F, const Vector<Type> &X0,
                         const Type tol=Type(1.0e-6),
                         const Type eps=Type(1.0e-9),
                         int maxItr=MAXTERM );


    #include <nleroots-impl.h>

}
// namespace splab


#endif
// NLEROOTS_H
