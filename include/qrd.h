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
 *                                    qrd.h
 *
 * Class template of QR decomposition for real matrix.
 *
 * For an m-by-n matrix A, the QR decomposition is an m-by-m orthogonal
 * matrix Q and an m-by-n upper triangular matrix R so that A = Q*R.
 *
 * For economy size, denotes p = min(m,n), then Q is m-by-p, and R is n-by-p,
 * this file provides the economic decomposition format.
 *
 * The QR decompostion always exists, even if the matrix does not have full
 * rank, so the constructor will never fail. The Q and R factors can be
 * retrived via the getQ() and getR() methods. Furthermore, a solve() method
 * is provided to find the least squares solution of Ax=b or AX=B using the
 * QR factors.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef QRD_H
#define QRD_H


#include <matrix.h>


namespace splab
{

	template <typename Real>
	class QRD
	{

	public:

		QRD();
		~QRD();

		void dec( const Matrix<Real> &A );
		bool isFullRank() const;

		Matrix<Real> getQ();
		Matrix<Real> getR();
		Matrix<Real> getH();

		Vector<Real> solve( const Vector<Real> &b );
		Matrix<Real> solve( const Matrix<Real> &B );

    private:

		// internal storage of QR
		Matrix<Real> QR;

		// diagonal of R.
		Vector<Real> RDiag;

	};
	// class QRD


    #include <qrd-impl.h>

}
// namespace splab


#endif
// QRD_H
