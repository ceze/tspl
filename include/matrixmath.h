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
 *                               matrixmath.h
 *
 * This file provides the basic math functions such as:
 *              cos    sin    tan    acos   asin   atan
 *              abs    exp    log    log10  sqrt   pow
 *
 * When debugging, use #define BOUNDS_CHECK above your "#include matrix.h"
 * line. When done debugging, comment out #define BOUNDS_CHECK for better
 * performance.
 *
 * Zhang Ming, 2010-01 (revised 2010-08), Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef MATRIXMATH_H
#define MATRIXMATH_H


#include <matrix.h>


namespace splab
{

    template<typename Type> Matrix<Type> abs( const Matrix<Type>& );
    template<typename Type> Matrix<Type> cos( const Matrix<Type>& );
    template<typename Type> Matrix<Type> sin( const Matrix<Type>& );
    template<typename Type> Matrix<Type> tan( const Matrix<Type>& );
    template<typename Type> Matrix<Type> acos( const Matrix<Type>& );
    template<typename Type> Matrix<Type> asin( const Matrix<Type>& );
    template<typename Type> Matrix<Type> atan( const Matrix<Type>& );

    template<typename Type> Matrix<Type> exp( const Matrix<Type>& );
    template<typename Type> Matrix<Type> log( const Matrix<Type>& );
    template<typename Type> Matrix<Type> log10( const Matrix<Type>& );

    template<typename Type> Matrix<Type> sqrt( const Matrix<Type>& );
    template<typename Type> Matrix<Type> pow( const Matrix<Type>&,
                                              const Matrix<Type>& );
    template<typename Type> Matrix<Type> pow( const Matrix<Type>&,
                                              const Type& );
    template<typename Type> Matrix<Type> pow( const Type&,
                                              const Matrix<Type>& );


    #include <matrixmath-impl.h>

}
// namespace splab


#endif
// MATRIXMATH_H
