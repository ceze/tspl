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
 *                                    svd.h
 *
 * Class template of Singular Value Decomposition.
 *
 * For an m-by-n matrix A, the singular value decomposition is an m-by-m
 * orthogonal matrix U, an m-by-n diagonal matrix S, and an n-by-n orthogonal
 * matrix V so that A = U*S*V^T.
 *
 * The singular values, sigma[k] = S[k][k], are ordered so that sigma[0] >=
 * sigma[1] >= ... >= sigma[n-1].
 *
 * For economy size, denotes p = min(m,n), then U is m-by-p, S is p-by-p,
 * and V is n-by-p, this file provides the economic decomposition format.
 *
 * The singular value decompostion always exists, so the constructor will
 * never fail. The matrix condition number and the effective numerical rank
 * can be computed from this decomposition.
 *
 * Adapted from Template Numerical Toolkit.
 *
 * Zhang Ming, 2010-01 (revised 2010-08), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef SVD_H
#define SVD_H


#include <matrix.h>


namespace splab
{

	template <typename Real>
	class SVD
	{

	public:

        SVD();
		~SVD();

        void dec( const Matrix<Real> &A );
		Matrix<Real> getU() const;
		Matrix<Real> getV() const;
		Matrix<Real> getSM();
		Vector<Real> getSV() const;

		Real norm2() const;
		Real cond() const;
		int  rank();

    private:

		// the orthogonal matrix and singular value vector
		Matrix<Real> U;
		Matrix<Real> V;
		Vector<Real> S;

        // docomposition for matrix with rows >= columns
		void decomposition( Matrix<Real>&, Matrix<Real>&,
                            Vector<Real>&, Matrix<Real>& );

	};
	// class SVD


    #include <svd-impl.h>

}
// namespace splab


#endif
// SVD_H
