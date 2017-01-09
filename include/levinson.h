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
 *                                  levinson.h
 *
 * If the coefficient matrix of a linear equations is Toeplitz, then it can
 * be solved in a high computational efficiency way through Levinson-Durbin
 * algorithm. The subroutiones in this file will be used for solving Wiener
 * -Hopf equeations in Wiener filtring and Yule- Walker equations in
 * parametric spectrum estimation, respectively.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef LEVINSON_H
#define LEVINSON_H


#include <vector.h>


namespace splab
{

    template<typename Type>
    Vector<Type> levinson( const Vector<Type>&, const Vector<Type>& );

    template<typename Type>
    Vector<Type> levinson( const Vector<Type>&, Type& );


    #include <levinson-impl.h>

}
// namespace splab


#endif
// LEVINSON_H
