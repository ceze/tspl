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
 *                                 statistics.h
 *
 * Some generally used routines about probability and statistics, such as
 * median, mean, variance, standard variance, skew, kurtosis and the
 * PDF(probability density function).
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef STATISTICS_H
#define STATISTICS_H


#include <vector.h>


namespace splab
{

    template<typename Type> Type mid( const Vector<Type>& );
    template<typename Type> Type mean( const Vector<Type>& );

    template<typename Type> Type var( const Vector<Type>& );
    template<typename Type> Type stdVar( const Vector<Type>& );
    template<typename Type> Vector<Type> standard( const Vector<Type>& );

    template<typename Type> Type skew( const Vector<Type>& );
    template<typename Type> Type kurt( const Vector<Type>& );

    template<typename Type> Vector<Type> pdf( Vector<Type>&,
                                              const Type lambda=Type(1.0) );


    #include <statistics-impl.h>

}
// namespace splab


#endif
// STATISTICS_H
