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
 *                               spline3interp-impl.h
 *
 * Implementation for Spline3Interp class.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Type>
Spline3Interp<Type>::Spline3Interp( const Vector<Type> &xn,
                                    const Vector<Type> &yn,
                                    Type d2l, Type d2r )
                     : Interpolation<Type>( xn, yn ), M0(d2l), Mn(d2r)
{
}

template <typename Type>
Spline3Interp<Type>::~Spline3Interp()
{
}


/**
 * Compute the second derivative at interpolated points.
 */
template <typename Type>
void Spline3Interp<Type>::derivative2( Vector<Type> &dx, Vector<Type> &d1,
                                  Vector<Type> &d2 )
{
    int N = xi.size(),
        M = N-1;
    Vector<Type> b(M),
                 v(M),
                 y(M),
                 alpha(M),
                 beta(M-1);

	for( int i=1; i<M; ++i )
		b[i] = 2 * (dx[i]+dx[i-1]);

	v[1] = 6*(d1[1]-d1[0]) - dx[0]*d2[0];
	for( int i=1; i<M-1; ++i )
		v[i] = 6 * (d1[i]-d1[i-1]);
	v[M-1] = 6*(d1[M-1]-d1[M-2]) - dx[M-1]*d2[M];

	alpha[1] = b[1];
	for( int i=2; i<M; ++i )
		alpha[i] = b[i] - dx[i]*dx[i-1]/alpha[i-1];

	for( int i=1; i<M-1; ++i )
		beta[i] = dx[i]/alpha[i];

	y[1] = v[1]/alpha[1];
	for( int i=2; i<M; ++i )
		y[i] = (v[i]-dx[i]*y[i-1]) / alpha[i];

	d2[M-1] = y[M-1];
	for( int i=M-2; i>0; --i )
		d2[i] = y[i] - beta[i]*d2[i+1];
}


/**
 * Compute the polynomial' coefsficient in each interval.
 */
template <typename Type>
void Spline3Interp<Type>::calcCoefs()
{
    int N = xi.size(),
        M = N-1;

    Vector<Type> m(N),
                 h(M),
                 d(M);

    m[0] = M0;
    m[M] = Mn;
    for( int i=0; i<M; ++i )
    {
        h[i] = xi[i+1]-xi[i];
        d[i] = (yi[i+1]-yi[i]) / h[i];
    }

    derivative2( h, d, m );

	coefs.resize( M, 4 );
	for( int i=0; i<M; ++i )
	{
		coefs[i][0] = yi[i];
		coefs[i][1] = d[i] - h[i]*(2*m[i]+m[i+1])/6;
		coefs[i][2] = m[i] / 2;
		coefs[i][3] = (m[i+1]-m[i]) / (6*h[i]);
	}
}


/**
 * Compute the value of polynomial at given "x".
 */
template <typename Type>
Type Spline3Interp<Type>::evaluate( Type x )
{
    int k = -1,
        N = xi.size(),
        M = N-1;

	Type dx,
         y;

	for( int i=0; i<M; ++i )
	{
		if( (xi[i]<=x) && (xi[i+1]>=x) )
		{
			k = i;
			dx = x-xi[i];
			break;
		}
	}
	if(k!=-1)
	{
		y = ( ( coefs[k][3]*dx + coefs[k][2] ) * dx + coefs[k][1] ) * dx
		  + coefs[k][0];
		return y;
	}
	else
	{
		cerr << "The value is out of range!" << endl;
		return Type(0);
	}
}


/**
 * Get polynomial' coefsficients.
 */
template <typename Type>
inline Matrix<Type> Spline3Interp<Type>::getCoefs() const
{
	return coefs;
}
