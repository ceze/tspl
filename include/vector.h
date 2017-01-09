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
 *                                 vector.h
 *
 * Class template of vector which is designed for basic linear algebra
 * operations such as:
 *              v + k    k + v    v += k    v1 + v2    v1 += v2
 *              v - k    k - v    v -= k    v1 - v2    v1 -= v2
 *              v * k    k * v    v *= k    v1 * v2    v1 *= v2
 *              v / k    k / v    v /= k    v1 / v2    v1 /= v2
 *              mum,     min,     max       swap       reverse
 *              norm     dotProd
 * These operators and functions can be applied to both real vector and
 * complex vector.
 *
 * The class also provides the basic math functions such as:
 *              cos    sin    tan    acos   asin   atan
 *              abs    exp    log    log10  sqrt   pow
 * This should include "matrixmath.h" file.
 *
 * When debugging, use #define BOUNDS_CHECK above your "#include vector.h"
 * line. When done debugging, comment out #define BOUNDS_CHECK for better
 * performance.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef VECTOR_H
#define VECTOR_H


#include <iostream>
#include <cassert>
#include <cmath>
#include <complex>
#include <usingdeclare.h>
#include <constants.h>


namespace splab
{

    template <typename Type>
    class Vector
    {

    public:

        typedef         Type*   iterator;
        typedef const   Type*   const_iterator;

        // constructors and destructor
        Vector();
        Vector( const Vector<Type> &v );
        Vector( int length, const Type &x = Type(0) );
        Vector( int length, const Type *array );
        ~Vector();

        // assignments
        Vector<Type>& operator=( const Vector<Type> &v );
        Vector<Type>& operator=( const Type &x );

        // accessors
        Type& operator[]( int i );
        const Type& operator[]( int i ) const;
        Type& operator()( int i );
        const Type& operator()( int i ) const;

        // iterators
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;

        // type conversion
        operator Type*();
        operator const Type*() const;

        // others
        int size() const;
        int dim() const;
        Vector<Type>& resize( int length );

        // computed assignment
        Vector<Type>& operator+=( const Type& );
        Vector<Type>& operator-=( const Type& );
        Vector<Type>& operator*=( const Type& );
        Vector<Type>& operator/=( const Type& );
        Vector<Type>& operator+=( const Vector<Type>& );
        Vector<Type>& operator-=( const Vector<Type>& );
        Vector<Type>& operator*=( const Vector<Type>& );
        Vector<Type>& operator/=( const Vector<Type>& );

    private:

        // data pointer for 0-offset indexing
        Type *pv0;

        // data pointer for 1-offset indexing
        Type *pv1;

        // the row number of vector
        int	 nRow;

        void init( int length );
        void copyFromArray( const Type *v );
        void setByScalar( const Type &x );
        void destroy();

    };
    // class Vector


    // input and output
    template<typename Type>
    ostream& operator<<( ostream&, const Vector<Type>& );
    template<typename Type>
    istream& operator>>( istream&, Vector<Type>& );

    // arithmetic operators
    template<typename Type>
    Vector<Type> operator-( const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator+( const Vector<Type>&, const Type& );
    template<typename Type>
    Vector<Type> operator+( const Type&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator+( const Vector<Type>&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator-( const Vector<Type>&, const Type& );
    template<typename Type>
    Vector<Type> operator-( const Type&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator-( const Vector<Type>&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator*( const Vector<Type>&, const Type& );
    template<typename Type>
    Vector<Type> operator*( const Type&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator*( const Vector<Type>&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator/( const Vector<Type>&, const Type& );
    template<typename Type>
    Vector<Type> operator/( const Type&, const Vector<Type>& );
    template<typename Type>
    Vector<Type> operator/( const Vector<Type>&, const Vector<Type>& );

    // dot product
    template<typename Type>
    Type dotProd( const Vector<Type>&, const Vector<Type>& );
    template<typename Type> complex<Type>
    dotProd( const Vector<complex<Type> >&, const Vector<complex<Type> >& );

    // utilities
    template<typename Type> Type sum( const Vector<Type>& );
    template<typename Type> Type min( const Vector<Type>& );
    template<typename Type> Type max( const Vector<Type>& );
    template<typename Type> Type norm( const Vector<Type>& );
    template<typename Type> Type norm( const Vector<complex<Type> >& );
    template<typename Type> void swap( Vector<Type>&, Vector<Type>& );
    template<typename Type> Vector<Type> linspace( Type, Type, int );
    template<typename Type> Vector<Type> abs( const Vector<complex<Type> >& );
    template<typename Type> Vector<Type> arg( const Vector<complex<Type> >& );
    template<typename Type> Vector<Type> real( const Vector<complex<Type> >& );
    template<typename Type> Vector<Type> imag( const Vector<complex<Type> >& );
    template<typename Type>
    Vector<complex<Type> > complexVector( const Vector<Type>& );
    template<typename Type>
    Vector<complex<Type> > complexVector( const Vector<Type>&,
                                          const Vector<Type>&  );


    #include <vector-impl.h>

}
// namespace splab


#endif
// VECTOR_H
