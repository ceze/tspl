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
 *                                    cevd.h
 *
 * Class template of eigenvalues decomposition for complex matrix.
 *
 * For a complex matrix A, we have A*V = V*D, where the eigenvalue matrix D
 * is diagonal and the eigenvector matrix V is linear independence. That is,
 * the diagonal values of D are the eigenvalues and the columns of V represent
 * the corresponding eigenvectors of D. If A is Hermitian, then V is a unitary
 * matrix, which means A = V*D*V', and eigenvalues are all real numbers.
 *
 * The matrix V may be badly conditioned, or even singular, so the validity
 * of the equation A=V*D*inverse(V) depends upon the condition number of V.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CEVD_H
#define CEVD_H


#include <evd.h>


namespace splab
{

	template <typename Type>
	class CEVD
	{

    public:

        CEVD();
		~CEVD();

        // decomposition
		void dec( const Matrix<complex<Type> > &A );

        // the eigenvalues are real or complex
		bool isHertimian() const;

        // get eigenvectors and
		Matrix<complex<Type> > getV() const;
		Vector<complex<Type> > getD() const;
		Vector<Type> getRD() const;

    private:

        bool hermitian;

		// eigenvectors and eigenvalues
		Matrix<complex<Type> > V;
		Vector<complex<Type> > d;
		Vector<Type> rd;

	};
	// class CEVD


    #include <cevd-impl.h>

}
// namespace splab


#endif
// CEVD_H
