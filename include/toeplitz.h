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
 *                                  toeplitz.h
 *
 * A Toeplitz matrix is defined by one row and one column. A symmetric
 * Toeplitz matrix is defined by just one row. Toeplitz generates Toeplitz
 * matrices given just the row or row and column description.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef TOEPLITZ_H
#define TOEPLITZ_H


#include <vector.h>
#include <matrix.h>


namespace splab
{
    template<typename Type> Matrix<Type> toeplitz( const Vector<Type>&,
                                                   const Vector<Type>& );
    template<typename Type> Matrix<Type> toeplitz( const Vector<Type>& );


    #include <toeplitz-impl.h>

}
// namespace splab


#endif
// TOEPLITZ_H
