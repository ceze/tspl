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
 *                                linequs3.h
 *
 * Function template for solving Rank Defect linear equations.
 *
 * For a m-by-n coefficient matrix A and m-by-1 constant vector b, if A is a
 * Rank Deficient Matrix, then the conventional methods will fail to solve
 * Ax = b. Three methods for solving such problem are provided in this file,
 * they are: Truncated SVD, Dampted SVD and Tikhonov Regularization.
 *
 * These methods adapt to both REAL or COMPLEX linear equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef RANKDEFLINEQUS_H
#define RANKDEFLINEQUS_H


#include <svd.h>
#include <csvd.h>


namespace splab
{

	template<typename Real>
	Vector<Real> tsvd( const Matrix<Real>&, const Vector<Real>&,
                       Real tol=Real(-1.0) );
	template<typename Type>
	Vector<complex<Type> > tsvd( const Matrix<complex<Type> >&,
                                 const Vector<complex<Type> >&, Type tol=Type(-1.0) );

	template<typename Real>
	Vector<Real> dsvd( const Matrix<Real>&, const Vector<Real>&, Real& );
	template<typename Type>
	Vector<complex<Type> > dsvd( const Matrix<complex<Type> >&,
                                 const Vector<complex<Type> >&, Type& );

	template<typename Real>
	Vector<Real> tikhonov( const Matrix<Real>&, const Vector<Real>&, Real& );
	template<typename Type>
	Vector<complex<Type> > tikhonov( const Matrix<complex<Type> >&,
                                     const Vector<complex<Type> >&, Type& );


	#include <linequs3-impl.h>

}
// namespace splab


#endif
// RANKDEFLINEQUS_H
