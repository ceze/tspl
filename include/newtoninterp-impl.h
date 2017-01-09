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
 *                               newtoninterp-impl.h
 *
 * Implementation for NewtonInterp class.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Type>
NewtonInterp<Type>::NewtonInterp( const Vector<Type> &xn,
                                  const Vector<Type> &yn )
                    : Interpolation<Type>( xn, yn ), coefs(yn)
{
}

template <typename Type>
NewtonInterp<Type>::~NewtonInterp()
{
}


/**
 * Compute polynomial' coefsficients.
 */
template <typename Type>
void NewtonInterp<Type>::calcCoefs()
{
    int N = xi.size();
	for( int j=1; j<N; ++j )
		for( int i=N-1; i>=j; --i )
			coefs[i] = (coefs[i]-coefs[i-1]) / (xi[i]-xi[i-j]);
}


/**
 * Compute the value of polynomial at given "x".
 */
template <typename Type>
Type NewtonInterp<Type>::evaluate( Type x )
{
    int N = xi.size();
	Type y = 0,
         tmp = 0;

	for( int j=0; j<N; ++j )
	{
		tmp = coefs[j];
		for( int i=0; i<j; ++i )
			tmp *= x-xi[i];

		y += tmp;
	}

	return y;
}


/**
 * Get polynomial' coefsficients.
 */
template <typename Type>
inline Vector<Type> NewtonInterp<Type>::getCoefs() const
{
	return coefs;
}
