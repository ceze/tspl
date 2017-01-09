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
 *                               nleroots-impl.h
 *
 * Implementation for nonlinear equations rooting.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Bisection method for finding function root.
 */
template <typename Type>
Vector<Type> seidel( NLEqus<Type> &G, const Vector<Type> &X0,
                     const Type tol, int maxItr )
{
    int N = X0.dim();
    Vector<Type> X(X0), XNew(X0), tmp(N);

    for( int k=0; k<maxItr; ++k )
    {
        // update X by the Seidel iteration method
        for( int i=1; i<=N; ++i )
        {
            tmp = G(XNew);
            XNew(i) = tmp(i);
        }

        // conditional judgement for stopping iteration
        Type err = norm( XNew-X ),
             relErr = err / ( norm(XNew) + EPS );
        if( (err < tol) || (relErr < tol) )
            return XNew;

        X = XNew;
    }

    cout << "No solution for the specified tolerance!" << endl;
	return XNew;
}


/**
 * Newton method for finding function root.
 */
template <typename Type>
Vector<Type> newton( NLFuncs<Type> &F, const Vector<Type> &X0,
                     const Type tol, const Type eps, int maxItr )
{
    Vector<Type> X(X0), XNew( X0.dim() );

    for( int k=0; k<maxItr; ++k )
    {
        // update X by the Newton iteration method
        XNew = X - luSolver( F.jacobi(X), F(X) );

        // conditional judgement for stopping iteration
        Type err = norm( XNew-X ),
             relErr = err / ( norm(XNew) + EPS );
        if( (err < tol) || (relErr < tol) || (norm(F(XNew)) < eps) )
            return XNew;

        X = XNew;
    }

    cout << "No solution for the specified tolerance!" << endl;
	return XNew;
}
