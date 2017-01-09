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
 *                                vectormath.h
 *
 * This file provides the basic math functions such as:
 *              cos    sin    tan    acos   asin   atan
 *              abs    exp    log    log10  sqrt   pow
 *
 * Zhang Ming, 2010-08, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef VECTORMATH_H
#define VECTORMATH_H


#include <vector.h>


namespace splab
{

    template<typename Type> Vector<Type> abs( const Vector<Type>& );
    template<typename Type> Vector<Type> cos( const Vector<Type>& );
    template<typename Type> Vector<Type> sin( const Vector<Type>& );
    template<typename Type> Vector<Type> tan( const Vector<Type>& );
    template<typename Type> Vector<Type> acos( const Vector<Type>& );
    template<typename Type> Vector<Type> asin( const Vector<Type>& );
    template<typename Type> Vector<Type> atan( const Vector<Type>& );

    template<typename Type> Vector<Type> exp( const Vector<Type>& );
    template<typename Type> Vector<Type> log( const Vector<Type>& );
    template<typename Type> Vector<Type> log10( const Vector<Type>& );

    template<typename Type> Vector<Type> sqrt( const Vector<Type>& );
    template<typename Type> Vector<Type> pow( const Vector<Type>&,
                                              const Vector<Type>& );
    template<typename Type> Vector<Type> pow( const Vector<Type>&,
                                              const Type& );
    template<typename Type> Vector<Type> pow( const Type&,
                                              const Vector<Type>& );
    template<typename Type>
    Vector<Type> gauss( const Vector<Type>&, const Type&, const Type& );


    #include <vectormath-impl.h>

}
// namespace splab


#endif
// VECTORMATH_H
