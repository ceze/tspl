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
 *                               fitcurves.h
 *
 * Fitted functions for least square fitting.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef FITCURVES_H
#define FITCURVES_H


#include <iostream>
#include <constants.h>


namespace splab
{

    template <typename Type>
    class Funcs
    {

    public:

        const static int dim = 4;

        // Compute the value of functions at point x.
        Type operator()( int i, const Type &x )
        {
            switch( i )
            {
                case 1:
                {
                    return 1;
                    break;
                }
                case 2:
                {
                    return log(max(x,EPS));
                    break;
                }
                case 3:
                {
                    return x;
                    break;
                }
                case 4:
                {
                    return x*x;
                    break;
                }
                default:
                {
                    std::cerr << "The dimension 'i' exceed the bound!"
                              << std::endl;
                    return 0;
                    break;
                }
            }
        }

    };
    // class ObjFunc

}
// namespace splab


#endif
// FITCURVES_H
