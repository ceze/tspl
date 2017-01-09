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
 *                                    evd.h
 *
 * Class template of eigenvalues and eigenvectors decomposition.
 *
 * For a real matrix A, we have A*V = V*D, where the eigenvalue matrix D is
 * diagonal and the eigenvector matrix V is linear independence. That is the
 * kth diagonal value of D is the eigenvalue and the kth column of V
 * represents the corresponding eigenvector of D[k][k]. If A is symmetric,
 * then V is a orthogonal matrix, which means A = V*D*V', and eigenvalues
 * are all real numbers.
 *
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal
 * with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 * a+i*b, in 2-by-2 blocks [a,b; -b,a]. That is, if the complex eigenvalues
 * look like
 *
 *    u + iv     .        .          .      .    .
 *      .      u - iv     .          .      .    .
 *      .        .      a + ib       .      .    .
 *      .        .        .        a - ib   .    .
 *      .        .        .          .      x    .
 *      .        .        .          .      .    y
 *
 * then D looks like
 *
 *      u        v        .          .      .    .
 *     -v        u        .          .      .    .
 *      .        .        a          b      .    .
 *      .        .       -b          a      .    .
 *      .        .        .          .      x    .
 *      .        .        .          .      .    y
 *
 * This keeps V a real matrix in both symmetric and non-symmetric cases, and
 * A*V = V*D.
 *
 * The matrix V may be badly conditioned, or even singular, so the validity
 * of the equation A=V*D*inverse(V) depends upon the condition number of V.
 *
 * Adapted from Template Numerical Toolkit.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef EVD_H
#define EVD_H


#include <matrix.h>


namespace splab
{

	template <typename Real>
	class EVD
	{

    public:

        EVD();
		~EVD();

        // decomposition
		void dec( const Matrix<Real> &A );

		// the eigenvalues are real or complex
		bool isSymmetric() const;
		bool isComplex( Real tol=Real(EPS) );

        // get eigenvectors
		Matrix<Real> getV() const;
		Matrix<complex<Real> > getCV();

        // get eigenvalues
        Vector<Real> getD() const;
        Vector<complex<Real> > getCD();
//        Matrix<Real> getDM();
//        Matrix<complex<Real> > getCDM();

    private:

		int     n;
		bool    symmetric;
        Real    cdivr,
                cdivi;

		// eigenvalues and its real and image part
		Matrix<Real> V;
		Vector<Real> d;
		Vector<Real> e;

		// temporary storage for internal variables
		Vector<Real> ort;
		Matrix<Real> H;

		void tred2();
		void tql2();
		void cdiv( Real xr, Real xi, Real yr, Real yi );
		void hqr2();
		void others();
		void normalized();

	};
	// class EVD


    #include <evd-impl.h>

}
// namespace splab


#endif
// EVD_H
