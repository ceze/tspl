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
 *                               pseudoinverse.h
 *
 * Matrix Pseudoinverse
 *
 * If 'A' is an m-by-n (where m != n) matrix, or n-by-n but nor full rank,
 * then we can't compute the ordinary inverse. In these cases, we should
 * compute the pseudoinverse. And here we use SVD to compute the pseudoinverse,
 * which is consistent with Matlab's "pinv".
 *
 * The algorithms provided in this file can be applied both for REAL matrix
 * or COMPLEX matrix.
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef PSEUDOINVERSE_H
#define PSEUDOINVERSE_H


#include <matrix.h>
#include <svd.h>
#include <csvd.h>


namespace splab
{

    template<typename Real>
    Matrix<Real> pinv( const Matrix<Real>&, Real tol=Real(-1.0) );

    template<typename Type>
    Matrix<complex<Type> > pinv( const Matrix<complex<Type> >&,
                                 Type tol=Type(-1.0) );


	#include <pseudoinverse-impl.h>
}
// namespace splab


#endif
// PSEUDOINVERSE_H
