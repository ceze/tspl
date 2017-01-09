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
 *                               lsfitting-impl.h
 *
 * Implementation for LSFitting class.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Type>
LSFitting<Type>::LSFitting( const Vector<Type> &xn, const Vector<Type> &yn,
                            Funcs<Type> &f )
                 : Interpolation<Type>( xn, yn ), phi(f)
{
    coefs.resize(phi.dim);
}

template <typename Type>
LSFitting<Type>::~LSFitting()
{
}


/**
 * Compute fitted parameters.
 */
template <typename Type>
void LSFitting<Type>::calcCoefs()
{
    int N = xi.size(),
        M = coefs.dim();

    Type tmp;
    Vector<Type> b(M);
    Matrix<Type> C( M, M );

	for( int i=1; i<=M; ++i )
	{
	    tmp = 0;
	    for( int k=0; k<N; ++k )
            tmp += phi(i,xi[k]) * phi(i,xi[k]);
        C(i,i) = tmp;

	    for( int j=i+1; j<=M; ++j )
	    {
	        tmp = 0;
	        for( int k=0; k<N; ++k )
                tmp += phi(i,xi[k]) * phi(j,xi[k]);
            C(i,j) = C(j,i) = tmp;
	    }

	    tmp = 0;
	    for( int k=0; k<N; ++k )
            tmp += phi(i,xi[k]) * yi[k];
        b(i) = tmp;
	}

    coefs = choleskySolver( C, b );
}


/**
 * Compute the value of fitted function at given "x".
 */
template <typename Type>
Type LSFitting<Type>::evaluate( Type x )
{
	Type y = 0;
	for( int j=0; j<coefs.size(); ++j )
		y += coefs[j] * phi( j, x );

	return y;
}


/**
 * Get the fitted parameters.
 */
template <typename Type>
inline Vector<Type> LSFitting<Type>::getCoefs() const
{
	return coefs;
}
