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
 *                                spline3interp.h
 *
 * Cubic splines interpolation method.
 *
 * For a given points set "Pn" {xn,yn}, this class can find a cubic polynomial
 * in each interval [x_i, x_i+1], such that both of the polynomials and their
 * first order derivative are continuous at the bound of each interval.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef SPLINE3INTERP_H
#define SPLINE3INTERP_H


#include <matrix.h>
#include <interpolation.h>


namespace splab
{

    template <typename Type>
    class Spline3Interp : public Interpolation<Type>
    {

     public:

        using Interpolation<Type>::xi;
        using Interpolation<Type>::yi;

        Spline3Interp( const Vector<Type> &xn, const Vector<Type> &yn,
                       Type d2l=Type(0), Type d2r=Type(0) );
        ~Spline3Interp();

        void calcCoefs();
        Type evaluate( Type x );
        Matrix<Type> getCoefs() const;

    private:

        void derivative2( Vector<Type> &dx, Vector<Type> &d1,
                          Vector<Type> &d2 );

        Type M0,
             Mn;
        Matrix<Type> coefs;

    };
    // class Spline3Interp


    #include <spline3interp-impl.h>

}
// namespace splab


#endif
// SPLINE3INTERP_H
