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
 *                                  nlfuncs.h
 *
 * Nonlinear equations and functions.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef NLFUNCS_H
#define NLFUNCS_H


#include <vector.h>
#include <matrix.h>


namespace splab
{

    template <typename Type>
    class NLEqus
    {

    public:

        /**
         * Compute the values of equations at point X.
         */
        Vector<Type> operator()( const Vector<Type> &X )
        {
            Vector<Type> XNew( X.dim() );

            XNew(1) = ( X(1)*X(1) - X(2) + Type(0.5) ) / 2;
            XNew(2) = ( -X(1)*X(1) - 4*X(2)*X(2) + 8*X(2) + 4 ) / 8;

            return XNew;
        }

    };
    // class NLEqus


    template <typename Type>
    class NLFuncs
    {

    public:

        /**
         * Compute the values of functions at point X.
         */
        Vector<Type> operator()( Vector<Type> &X )
        {
            Vector<Type> FX( X.dim() );

            FX(1) = X(1)*X(1) - 2*X(1) - X(2) + Type(0.5);
            FX(2) = X(1)*X(1) + 4*X(2)*X(2) - 4;

            return FX;
        }

        /**
         * Compute the Jacobian-Matrix at point X.
         */
        Matrix<Type> jacobi( Vector<Type> &X )
        {
            Matrix<Type> JX( X.dim(), X.dim() );

            JX(1,1) = 2*X(1) - 2;    JX(1,2) = Type(-1);
            JX(2,1) = 2*X(1);        JX(2,2) = 8*X(2);

            return JX;
        }

    };
    // class NLFuncs

}
// namespace splab


#endif
// NLFUNCS_H
