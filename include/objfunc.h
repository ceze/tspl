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
 *                                objfunc.h
 *
 * Function Object.
 *
 * We can use this functor computing the value of objective function and it's
 * gradient vector. The objective function is supposed to be multidimensional,
 * one dimention is the special case of "vector.dim()=1".
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef OBJFUNC_H
#define OBJFUNC_H


#include <vector.h>


namespace splab
{

    template <typename Type>
    class ObjFunc
    {

    public:

        /**
         * Initialize the parameters
         */
        ObjFunc( Type aa, Type bb, Type cc ) : a(aa), b(bb), c(cc)
        { }

        /**
         * Compute the value of objective function at point x.
         */
        Type operator()( Vector<Type> &x )
        {
            return ( a*x(1) * exp(b*x(1)*x(1)+c*x(2)*x(2)) );
        }

        /**
         * Compute the gradient of objective function at point x.
         */
        Vector<Type> grad( Vector<Type> &x )
        {
            Type expTerm = exp(a*x(1)*x(1)+b*x(2)*x(2));
            Vector<Type> df(x.dim());
            df(1) = (a+2*a*b*x(1)*x(1)) * expTerm;
            df(2) = 2*a*c*x(1)*x(2) * expTerm;

            return df;
        }

    private:

        // parameters
        Type a,
             b,
             c;

    };
    // class ObjFunc

}
// namespace splab


#endif
// OBJFUNC_H
