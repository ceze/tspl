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
 *                                newtoninterp.h
 *
 * Newton interpolation method.
 *
 * For a given points set "Pn" {xn,yn}, this class can find a polynomial which
 * cross all points of "Pn".
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef NEWTONINTERP_H
#define NEWTONINTERP_H


#include <interpolation.h>


namespace splab
{

    template <typename Type>
    class NewtonInterp : public Interpolation<Type>
    {

     public:

        using Interpolation<Type>::xi;
        using Interpolation<Type>::yi;

        NewtonInterp( const Vector<Type> &xn, const Vector<Type> &yn );
        ~NewtonInterp();

        void calcCoefs();
        Type evaluate( Type x );
        Vector<Type> getCoefs() const;

    private:

        Vector<Type> coefs;

    };
    // class NewtonInterp


    #include <newtoninterp-impl.h>

}
// namespace splab


#endif
// NEWTONINTERP_H
