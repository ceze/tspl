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
 *                               cholesky.h
 *
 * Class template of Cholesky decomposition.
 *
 * For a symmetric, positive definite matrix A, this function computes the
 * Cholesky factorization, i.e. it computes a lower triangular matrix L such
 * that A = L*L'. If the matrix is not symmetric or positive definite, the
 * function computes only a partial decomposition. This can be tested with
 * the isSpd() flag.
 *
 * This class also supports factorization of complex matrix by specializing
 * some member functions.
 *
 * Adapted from Template Numerical Toolkit.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CHOLESKY_H
#define CHOLESKY_H


#include <matrix.h>


namespace splab
{

    template <typename Type>
    class Cholesky
    {

    public:

        Cholesky();
        ~Cholesky();

        bool isSpd() const;
        void dec( const Matrix<Type> &A );
        Matrix<Type> getL() const;

        Vector<Type> solve( const Vector<Type> &b );
        Matrix<Type> solve( const Matrix<Type> &B );

    private:

        bool spd;

        Matrix<Type> L;

    };
    //	class Cholesky


    #include <cholesky-impl.h>

}
// namespace splab


#endif
// CHOLESKY_H
