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
 *                                linequs2.h
 *
 * Function template for solving linear equations.
 *
 * For a m-by-n (m!=n) coefficient matrix A and m-by-1 constant vector b, if
 * m>n, the exact solution of Ax=b is not existent, but we can find the least
 * square solution that minimize norm of Ax-b; if m<n, there are infinite
 * many solutions, but we can find the minimum norm sulution that minimize
 * norm of x.
 *
 * These functions can be used for both REAL or COMPLEX linear equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef UNDETLINEQUS_H
#define UNDETLINEQUS_H


#include <qrd.h>
#include <svd.h>
#include <cqrd.h>
#include <csvd.h>
#include <linequs1.h>


namespace splab
{

	template<typename Type>
	Vector<Type> lsSolver( const Matrix<Type>&, const Vector<Type>& );
	template<typename Real>
	Vector<Real> qrLsSolver( const Matrix<Real>&, const Vector<Real>& );
	template<typename Real>
	Vector<Real> svdLsSolver( const Matrix<Real>&, const Vector<Real>& );
	template<typename Type>
	Vector<complex<Type> > qrLsSolver( const Matrix<complex<Type> >&,
                                       const Vector<complex<Type> >& );
	template<typename Type>
	Vector<complex<Type> > svdLsSolver( const Matrix<complex<Type> >&,
                                        const Vector<complex<Type> >& );

	template<typename Type>
	Vector<Type> lnSolver( const Matrix<Type>&, const Vector<Type>& );
	template<typename Real>
	Vector<Real> qrLnSolver( const Matrix<Real>&, const Vector<Real>& );
	template<typename Real>
	Vector<Real> svdLnSolver( const Matrix<Real>&, const Vector<Real>& );
	template<typename Type>
	Vector<complex<Type> > qrLnSolver( const Matrix<complex<Type> >&,
                                       const Vector<complex<Type> >& );
	template<typename Type>
	Vector<complex<Type> > svdLnSolver( const Matrix<complex<Type> >&,
                                        const Vector<complex<Type> >& );


	#include <linequs2-impl.h>

}
// namespace splab


#endif
// UNDETLINEQUS_H
