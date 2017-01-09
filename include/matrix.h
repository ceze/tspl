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
 *                                 matrix.h
 *
 * Class template of matrix which is designed for basic linear algebra
 * operations such as:
 *              A + x    x + A    A += x    A1 + A2   A1 += A2
 *              A - x    x - A    A -= x    A1 - A2   A1 -= A2
 *              A * x    x * A    A *= x    A1 * A2
 *              A / x    x / A    A /= x
 *              mum,     min,     max       mean      swap
 *              eye,     diag,    trT       trH       norm
 *              trMult   multTr   elemMult  elemMultEq
 *              elemDivd elemDivdEq
 * These operators and functions can be applied to both real matrix and
 * complex matrix.
 *
 * The class also provides the basic math functions such as:
 *              cos    sin    tan    acos   asin   atan
 *              abs    exp    log    log10  sqrt   pow
 * This should include "matrixmath.h" file.
 *
 * When debugging, use #define BOUNDS_CHECK above your "#include matrix.h"
 * line. When done debugging, comment out #define BOUNDS_CHECK for better
 * performance.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef MATRIX_H
#define MATRIX_H


#include <vector.h>


namespace splab
{

    template <typename Type>
    class Matrix
    {

    public:

        // constructors and destructor
        Matrix();
        Matrix( const Matrix<Type> &A );
        Matrix( int rows, int columns, const Type &x = Type(0) );
        Matrix( int rows, int columns, const Type *v );
        ~Matrix();

        // assignments
        Matrix<Type>& operator=( const Matrix<Type> &A );
        Matrix<Type>& operator=( const Type &x );

        // accessors
        Type* operator[]( int i );
        const Type* operator[]( int i ) const;
        Type& operator()( int row, int column );
        const Type& operator()( int row, int column ) const;

        // type conversion
        operator Type*();
        operator const Type*() const;
//        operator Type**();
//        operator const Type**() const;

        // others
        long size() const;
        int dim( int dimension ) const;
        int rows() const;
        int cols() const;
        Matrix<Type>& resize( int rows, int columns );
        Vector<Type> getRow( int row ) const;
        Vector<Type> getColumn( int column ) const;
        void setRow( const Vector<Type> &v, int row );
        void setColumn( const Vector<Type> &v, int column );

        // computed assignment
        Matrix<Type>& operator+=( const Type& );
        Matrix<Type>& operator+=( const Matrix<Type>& );
        Matrix<Type>& operator-=( const Type& );
        Matrix<Type>& operator-=( const Matrix<Type>& );
        Matrix<Type>& operator*=( const Type& );
        // WARNING: element-by-element
        Matrix<Type>& operator*=( const Matrix<Type>& );
        Matrix<Type>& operator/=( const Type& );
        // WARNING: element-by-element
        Matrix<Type>& operator/=( const Matrix<Type>& );

    private:

        // 0-based and 1-based data pointer
        Type *pv0, *pv1;

        // 0-based and 1-based row pointer's pointer
        Type **prow0, **prow1;

        // row number, column number and total number
        int	 nRow;
        int	 nColumn;
        long nTotal;

        void init( int rows, int columns );
        void copyFromArray( const Type *v );
        void setByScalar( const Type &x );
        void destroy();

    };
    // class Matrix

    // input and output
    template<typename Type>
    ostream& operator<<( ostream&, const Matrix<Type>& );
    template<typename Type>
    istream& operator>>( istream&, Matrix<Type>& );

    // arithmetic operators
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator+( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator+( const Type&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator+( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator-( const Type&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator-( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator*( const Type&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> operator*( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Vector<Type> operator*( const Matrix<Type>&, const Vector<Type>& );
    template<typename Type>
    Matrix<Type> operator/( const Matrix<Type>&, const Type& );
    template<typename Type>
    Matrix<Type> operator/( const Type&, const Matrix<Type>& );

    // optimizied version of multiplication
    template<typename Type>
    Matrix<Type>& optMult( const Matrix<Type>&, const Matrix<Type>&,
                           Matrix<Type>& );
    template<typename Type>
    Vector<Type>& optMult( const Matrix<Type>&, const Vector<Type>&,
                           Vector<Type>& );

    // element-by-element multiplication and division
    template<typename Type>
    Matrix<Type> elemMult( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> elemDivd( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type>& elemMultEq( Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type>& elemDivdEq( Matrix<Type>&, const Matrix<Type>& );

    // transpose and conjugate transpose
    template<typename Type> Matrix<Type> trT( const Matrix<Type>& );
    template<typename Type> Matrix<Type> trH( const Matrix<Type>& );

    // transpose multiplication
    template<typename Type>
    Matrix<Type> trMult( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Vector<Type> trMult( const Matrix<Type>&, const Vector<Type>& );
    template<typename Type>
    Matrix<Type> multTr( const Matrix<Type>&, const Matrix<Type>& );
    template<typename Type>
    Matrix<Type> multTr( const Vector<Type>&, const Vector<Type>& );
    template<typename Type>
    Matrix<complex<Type> > trMult( const Matrix<complex<Type> >&,
                                   const Matrix<complex<Type> >& );
    template<typename Type>
    Vector<complex<Type> > trMult( const Matrix<complex<Type> >&,
                                   const Vector<complex<Type> >& );
    template<typename Type>
    Matrix<complex<Type> > multTr( const Matrix<complex<Type> >&,
                                   const Matrix<complex<Type> >& );
    template<typename Type>
    Matrix<complex<Type> > multTr( const Vector<complex<Type> >&,
                                   const Vector<complex<Type> >& );

    // unit and diagonal matrix
    template<typename Type> Matrix<Type> eye( int, const Type& );
    template<typename Type> Vector<Type> diag( const Matrix<Type>& );
    template<typename Type> Matrix<Type> diag( const Vector<Type>& );

    // utilities
    template<typename Type> Type norm( const Matrix<Type>& );
    template<typename Type> Type norm( const Matrix<complex<Type> >& );
    template<typename Type> void swap( Matrix<Type>&, Matrix<Type>& );
    template<typename Type> Vector<Type> sum( const Matrix<Type>& );
    template<typename Type> Vector<Type> min( const Matrix<Type>& );
    template<typename Type> Vector<Type> max( const Matrix<Type>& );
    template<typename Type> Vector<Type> mean( const Matrix<Type>& );
    template<typename Type> Matrix<Type> abs( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> arg( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> real( const Matrix<complex<Type> >& );
    template<typename Type> Matrix<Type> imag( const Matrix<complex<Type> >& );
    template<typename Type>
    Matrix<complex<Type> > complexMatrix( const Matrix<Type>& );
    template<typename Type>
    Matrix<complex<Type> > complexMatrix( const Matrix<Type>&,
                                          const Matrix<Type>& );


    #include <matrix-impl.h>

}
// namespace splab


#endif
// MATRIX_H
