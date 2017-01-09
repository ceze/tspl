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
 *                               linequs1.h
 *
 * Function template for solving deterministic linear equations.
 *
 * For a n-by-n coefficient matrix A and n-by-1 constant vector b, Gauss
 * elimination method can solve the equations Ax=b. But the decomposition
 * methods are more commonly used. if A is not SPD, the LU decomposition
 * can be use to solve the equations; if A is SPD, the Cholesky decomposition
 * is the better choice; if A is a tridiagonal matrix, then the Forward
 * Elimination and Backward Substitution maybe the best choice.
 *
 * These functions can be used for both REAL or COMPLEX linear equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef DETLINEQUS_H
#define DETLINEQUS_H


#include <vector.h>
#include <matrix.h>
#include <cholesky.h>
#include <lud.h>


namespace splab
{
    template<typename Type>
    Vector<Type> gaussSolver( const Matrix<Type>&, const Vector<Type>& );
    template<typename Type>
    void gaussSolver( Matrix<Type>&, Matrix<Type>& );

    template<typename Type>
    Vector<Type> luSolver( const Matrix<Type>&, const Vector<Type>& );
    template<typename Type>
    Matrix<Type> luSolver( const Matrix<Type>&, const Matrix<Type>& );

    template<typename Type>
    Vector<Type> choleskySolver( const Matrix<Type>&, const Vector<Type>& );
	template<typename Type>
	Matrix<Type> choleskySolver( const Matrix<Type>&, const Matrix<Type>& );

    template<typename Type>
    Vector<Type> utSolver( const Matrix<Type>&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> ltSolver( const Matrix<Type>&, const Vector<Type>& );

	template<typename Type>
	Vector<Type> febsSolver( const Vector<Type>&, const Vector<Type>&,
                             const Vector<Type>&, const Vector<Type>& );


	#include <linequs1-impl.h>

}
// namespace splab


#endif
// DETLINEQUS_H
