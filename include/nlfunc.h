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
 *                                  nlfunc.h
 *
 * Nonlinear function.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef NLFUNC_H
#define NLFUNC_H


namespace splab
{

    template <typename Type>
    class NLFunc
    {

    public:

        /**
         * Initialize the parameters
         */
        NLFunc( Type aa, Type bb, Type cc ) : a(aa), b(bb), c(cc)
        { }

        /**
         * Compute the value of objective function at point x.
         */
        Type operator()( Type &x )
        {
            return ( a*x*x + b*x + c );
        }

        /**
         * Compute the gradient of objective function at point x.
         */
        Type grad( Type &x )
        {
            return ( 2*a*x + b );
        }

    private:

        // parameters
        Type a,
             b,
             c;

    };
    // class NLFunc

}
// namespace splab


#endif
// NLFUNC_H
