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
 *                               nleroot-impl.h
 *
 * Implementation for nonlinear equation rooting.
 *
 * Zhang Ming, 2010-04, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Bisection method for finding function root.
 */
template <typename Type>
Type bisection( NLFunc<Type> &f, Type a, Type b, Type tol )
{
	Type    c = 0,
            fc = 0,
            fa = f(a),
            fb = f(b);

    if( a > b )
    {
        c = a;
        a = b;
        b = c;
    }
    if( fa*fb > 0 )
    {
        cerr << " Invalid interval!" << endl;
        return Type(0);
    }

    int maxItr = (int) ceil( (log(b-a)-log(tol)) / log(2.0) );
	for( int i=0; i<maxItr; ++i )
	{
	    c = ( a+ b) / 2;
	    fc = f(c);

	    if( abs(fc) < EPS )
            return c;

	    if( fb*fc > 0 )
	    {
	        b = c;
	        fb = fc;
	    }
	    else
	    {
	        a = c;
	        fa = fc;
	    }

	    if( (b-a) < tol )
            return (b+a)/2;
	}

	cout << "No solution for the specified tolerance!" << endl;
	return c;
}


/**
 * Newton method for finding function root.
 */
template <typename Type>
Type newton( NLFunc<Type> &f, Type x0, Type tol, int maxItr )
{
	Type    xkOld = x0,
            xkNew,
            fkGrad;

	for( int i=0; i<maxItr; ++i )
	{
	    fkGrad = f.grad(xkOld);
	    while( abs(fkGrad) < EPS )
	    {
	        xkOld *= 1.2;
	        fkGrad = f.grad(xkOld);
	    }
		xkNew = xkOld - f(xkOld)/fkGrad;

		if( abs(f(xkNew)) < tol )
			return xkNew;
		if( abs(xkOld-xkNew) <tol )
			return xkNew;

		xkOld = xkNew;
	}

	cout << "No solution for the specified tolerance!" << endl;
	return xkNew;
}


/**
 * Secant method for finding function root.
 */
template <typename Type>
Type secant( NLFunc<Type> &f, Type x1, Type x2, Type tol, int maxItr )
{
	Type    x_k,
            x_k1 = x1,
            x_k2 = x2,
            f_k1 = f(x_k1),
            f_k2;

	for( int i=0; i<maxItr; ++i )
	{
		f_k2 = f(x_k2);
		if( abs(f_k2-f_k1) < EPS )
        {
            x_k2 = (x_k1+x_k2) / 2;
            f_k2 = f(x_k2);
        }

		x_k = x_k2 - (x_k2-x_k1)*f_k2/(f_k2-f_k1);

		if( abs(f(x_k)) < tol )
			return x_k;
		if( abs(x_k2-x_k) < tol )
			return x_k;

		f_k1 = f_k2;
		x_k1 = x_k2;
		x_k2 = x_k;
	}

	cout << "No solution for the specified tolerance!" << endl;
	return x_k;
}
