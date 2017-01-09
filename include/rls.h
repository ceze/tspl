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
 *                                   rls.h
 *
 * Recursive Least Square Filter.
 *
 * The RLS adaptive filter recursively finds the filter coefficients that
 * minimize a weighted linear least squares cost function relating to the
 * input signals. This in contrast to other algorithms such as the LMS that
 * aim to reduce the mean square error.
 *
 * The input signal of RLS adaptive filter is considered deterministic,
 * while for the LMS and similar algorithm it is considered stochastic.
 * Compared to most of its competitors, the RLS exhibits extremely fast
 * convergence. However, this benefit comes at the cost of high computational
 * complexity, and potentially poor tracking performance when the filter to
 * be estimated changes.
 *
 * This file includes five types usually used RLS algorithms, they are:
 * conventional RLS (rls),       stabilised fast transversal RLS (sftrls),
 * lattice RLS (lrls),           error feedblck lattice RLS (eflrls),
 * QR based RLS (qrrls).
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef RLS_H
#define RLS_H


#include <vector.h>
#include <matrix.h>


namespace splab
{

    template<typename Type>
    Type rls( const Type&, const Type&, Vector<Type>&,
              const Type&, const Type& );

    template<typename Type>
    Type sftrls( const Type&, const Type&, Vector<Type>&,
                 const Type&, const Type&, const string& );

    template<typename Type>
    Type lrls( const Type&, const Type&, Vector<Type>&,
               const Type&, const Type&, const string& );

    template<typename Type>
    Type eflrls( const Type&, const Type&, Vector<Type>&,
                 const Type&, const Type&, const string& );

    template<typename Type>
    Type qrrls( const Type&, const Type&, Vector<Type>&,
                const Type&, const string& );


    #include <rls-impl.h>

}
// namespace splab


#endif
// RLS_H
