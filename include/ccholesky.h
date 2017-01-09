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
 *                               ccholesky.h
 *
 * Class template of complex matrix Cholesky factorization.
 *
 * For a conjugate symmetric, positive definite matrix A, this class computes
 * the Cholesky factorization, i.e. it computes a lower triangular matrix L
 * such that A = L*L^H. If the matrix is not conjugate symmetric or positive
 * definite, the class computes only a partial decomposition. This can be
 * tested with the isSpd() flag.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CCHOLESKY_H
#define CCHOLESKY_H


#include <matrix.h>


namespace splab
{

    template <typename Type>
    class CCholesky
    {

    public:

        CCholesky();
        ~CCholesky();

        bool isSpd() const;
        void dec( const Matrix<Type> &A );
        Matrix<Type> getL() const;

        Vector<Type> solve( const Vector<Type> &b );
        Matrix<Type> solve( const Matrix<Type> &B );

    private:

        bool spd;
        Matrix<Type> L;

    };
    //	class CCholesky


    #include <ccholesky-impl.h>

}
// namespace splab


#endif
// CCHOLESKY_H
