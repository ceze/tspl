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
 *                                 integral.h
 *
 * Implementation for Romberg method.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


template <typename Type>
Type romberg( Func<Type> &f, const Type &a, const Type &b, const Type tol )
{
    int n = 5,
        m = 2*n;
    Type x = 0,
         tmp = 0,
         result = 0;

    for( int k=0; k<MAXTERM; ++k )
    {
        Type h = (b-a)/m;
        for( int i=0; i<n; ++i )
        {
            x = a + (2*i+1)*h;
            tmp += 4*f(x);
        }
        for( int i=1; i<n; ++i )
        {
            x = a + 2*i*h;
            tmp += 2*f(x);
        }
        tmp += f(a)+f(b);
        tmp *= h/3;

        if( abs(tmp-result) < tol )
        {
            result = tmp;
            return result;
        }
        else
        {
            result = tmp;
            tmp = 0;
            n = m;
            m = 2*n;
        }
    }

    std::cerr << "Dodn't get the specified precision value!" << std::endl;
    return result;
}
