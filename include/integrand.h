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
 *                               integrand.h
 *
 * Integrand used in numerical integral.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef INTEGRAND_H
#define INTEGRAND_H


#include <cmath>
#include <constants.h>


namespace splab
{

    /**
     * Integrand, a function objective.
     */
    template <typename Type>
    class Func
    {

    public:

        /**
         * Initialize the parameters
         */
        Func( Type aa, Type bb ) : a(aa), b(bb)
        { }

        /**
         * Compute the value of objective function at point x.
         */
        Type operator()( Type x )
        {
            return a + b*sin(x);
        }

    private:

        Type a,
             b;

    };
    // class Func

}
// namespace splab


#endif
// INTEGRAND_H
