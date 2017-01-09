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
 *                                    lud.h
 *
 * Class template of LU decomposition.
 *
 * For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
 * unit lower triangular matrix L, an n-by-n upper triangular matrix U, and
 * a permutation vector piv of length m so that A(piv,:) = L*U. If m < n,
 * then L is m-by-m and U is m-by-n.
 *
 * The matrix can be both REAL or COMPLEX.
 *
 * The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail. The primary use of the LU
 * decomposition is in the solution of square systems of simultaneouslinear
 * equations. This will fail if isNonsingular() returns false.
 *
 * Adapted from Template Numerical Toolkit.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef LUD_H
#define LUD_H


#include <matrix.h>


namespace splab
{

	template <typename Type>
	class LUD
	{

	public :

        LUD();
		~LUD();

		void dec( const Matrix<Type> &A );
		Matrix<Type> getL();
		Matrix<Type> getU();
		Vector<int>  getPivot() const;

		Type det();
		bool isNonsingular();

		Vector<Type> solve( const Vector<Type> &b );
		Matrix<Type> solve( const Matrix<Type> &B );

    private:

        // Array for internal storage of decomposition.
		Matrix<Type> LU;
		int m;
		int n;

		int pivsign;
		Vector<int> piv;

		Vector<Type> permuteCopy( const Vector<Type> &A,
                                  const Vector<int> &piv );
		Matrix<Type> permuteCopy( const Matrix<Type> &A,
		                          const Vector<int> &piv, int j0, int j1 );

	};
	// class LUD


	#include <lud-impl.h>

}
// namespace splab


#endif
// LUD_H
