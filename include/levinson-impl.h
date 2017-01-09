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
 *                               levinson-impl.h
 *
 * Implementationfor Levinson-Durbin alogrithm.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Levinson algorithm for solving Toeplitz equations.
 * t    : t(0), t(1), ..., t(n-1) of Toeplitz coefficient matrix
 * b    : constant vector
 */
template <typename Type>
Vector<Type> levinson( const Vector<Type> &t, const Vector<Type> &b )
{
    assert( t.size() == b.size() );

    int n = t.size();
	Type alpha, beta, q, c, omega;
	Vector<Type> y(n), yy(n), x(n);

	alpha = t[0];
	if( abs(alpha) < EPS )
	{
		cerr << "The matrix is ill-conditioned!" << endl;
		return x;
	}
	y[0] = 1;
	x[0] = b[0] / alpha;

	for( int k=1; k<n; ++k )
	{
		q = 0;
		beta = 0;
		for( int j=0; j<k; ++j )
		{
			q += x[j] * t[k-j];
			beta += y[j] * t[j+1];
		}
		c = -beta / alpha;

		yy[0] = c * y[k-1];
		y[k] = y[k-1];
		for( int i=1; i<k; ++i )
			yy[i] = y[i-1] + c*y[k-i-1];
		yy[k] = y[k-1];

		alpha += c*beta;
		if( abs(alpha) < EPS )
		{
			cerr << "The matrix is ill-conditioned!" << endl;
			return x;
		}

		omega = (b[k]-q) / alpha;
		for( int i=0; i<k; ++i )
		{
			x[i] += omega*yy[i];
			y[i] = yy[i];
		}
		x[k] = omega*y[k];
	}

	return x;
}


/**
 * Levinson-Durbin algorithm for solving Youle-Walker equations.
 * rn       : r(0), r(1), ..., r(p)
 * sigma2   : the variance of exciting white noise
 */
template <typename Type>
Vector<Type> levinson( const Vector<Type> &rn, Type &sigma2 )
{
    int p = rn.size()-1;
    Type tmp;
    Vector<Type> ak(p+1), akPrev(p+1);

    ak[0] = Type(1.0);
    sigma2 = rn[0];
    ak[1] = -rn[1]/sigma2;
    sigma2 *= 1 - ak[1]*ak[1];

    for( int k=2; k<=p; ++k )
    {
        tmp = 0;
        for( int i=0; i<k; ++i )
            tmp += ak[i]*rn[k-i];
        ak[k] = -tmp/sigma2;

        for( int i=1; i<k; ++i )
            akPrev[i] = ak[i] + ak[k]*ak[k-i];
        for( int i=1; i<k; ++i )
            ak[i] = akPrev[i];

        sigma2 *= 1 - ak[k]*ak[k];
    }

	return ak;
}
