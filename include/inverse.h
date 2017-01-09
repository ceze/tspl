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
 *                                  inverse.h
 *
 * Matrix Inverse
 *
 * Matrix's inverse is a very important algorithm. If 'A' is an n-by-n full
 * rank square matrix, then it's inverse 'invA' is also an n-by-n full rank
 * square matrix. The general method to compute such matrix's inverse is
 * use LU decomposition to solve A*invA=I, if 'A' is SPD, then Cholesky
 * decomposition will be more efficiency.
 *
 * The algorithms provided in this file can be applied both for REAL matrix
 * or COMPLEX matrix.
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef INVERSE_H
#define INVERSE_H


#include <string>
#include <matrix.h>
#include <cholesky.h>
#include <lud.h>


namespace splab
{

	template<typename Type> Matrix<Type> inv( const Matrix<Type>&,
                                              const string &type="nspd" );
    template<typename Type>
    Matrix<complex<Type> > cinv( const Matrix<complex<Type> >&,
                                 const string &type="nspd" );

	template<typename Type> Matrix<Type> colPivInv( const Matrix<Type>& );
    template<typename Type> Matrix<Type> cmpPivInv( const Matrix<Type>& );


	#include <inverse-impl.h>
}
// namespace splab


#endif
// INVERSE_H
