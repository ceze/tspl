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
 *                               advmath-impl.h
 *
 * Implementationfor advance math functions.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Calcutes the inverse hyperbolic cosine of x.
 */
inline double acosh( double x )
{
    return log( x + sqrt(x*x-1) );
}


/**
 * Calculates the inverse hyperbolic sine of x.
 */
inline double asinh( double x )
{
    return log( x + sqrt(x*x + 1) );
}


/**
 * Calculates Jacobian elliptic arctan function using
 * arithmetic-geometric mean method.
 */
double arcsc( double u, double k )
{
    int     i, L;
    double  A, B, BT, Y;

    A = 1;
    B = k;
    Y = 1.0 / u;
    L = 0;

    for( i=0; i<MAXTERM; ++i )
    {
        BT = A * B;
        A = A + B;
        B = 2 * sqrt(BT);
        Y = Y - BT / Y;

        if( Y == 0 )
            Y = sqrt(BT) * EPS;
        if( abs(A-B) < (A*EPS) )
            break;

        L = 2 * L;
        if( Y < 0 )
            L++;
    }

    if( Y < 0 )
        L++;

    return ( (atan(A/Y) + PI * L) / A );
}


/**
 * Calculates Jacobian elliptic functions using
 * arithmetic-geometric mean method.
 */
void ellipticFun( double u, double k,
                  double *sn, double *cn, double *dn )
{
    int     i, imax;
    double  A[MAXTERM], B[MAXTERM],
            C[MAXTERM], P[MAXTERM];

    k = k * k;
    A[0] = 1;
    B[0] = sqrt(1-k);
    C[0] = sqrt(k);

    for( i=1; i<MAXTERM; ++i )
    {
        A[i] = (A[i-1] + B[i-1]) / 2;
        B[i] = sqrt( A[i-1] * B[i-1] );
        C[i] = (A[i-1] - B[i-1]) / 2;
        if( C[i] < EPS )
            break;
    }

    if( i == MAXTERM )
        imax = i - 1;
    else
        imax = i;

    P[imax] = pow(2.0,imax) * A[imax] * u;
    for( i=imax; i>0; i-- )
        P[i-1] = (asin(C[i]*sin(P[i])/A[i]) + P[i]) / 2;

    *sn = sin(P[0]);
    *cn = cos(P[0]);
    *dn = sqrt( 1 - k*(*sn)*(*sn) );
}


/**
 * Calculates complete elliptic integral using
 * arithmetic-geometric mean method.
 */
double ellipticIntegral( double k )
{
    int     i;
    double  A[MAXTERM], B[MAXTERM], C[MAXTERM];

    k = k*k;
    A[0] = 1;
    B[0] = sqrt(1-k);
    C[0] = sqrt(k);

    for( i=1; i<MAXTERM; ++i )
    {
        A[i] = (A[i-1] + B[i-1]) / 2;
        B[i] = sqrt(A[i-1]*B[i-1]);
        C[i] = (A[i-1] - B[i-1]) / 2;
        if( C[i] < EPS )
            break;
    }

    return PI / (2*A[i]);
}


/**
 * Factoring quadratic equation with complex coefficients.
 * Equation form is a*x^2 + b*x + c, solutions will be placed
 * in complex numbers pointers d and e.
 */
void quadradicRoot( complex<double> a, complex<double> b, complex<double> c,
                    complex<double> &x1, complex<double> &x2 )
{
    complex<double> a2,
                    ac4,
                    sq;

    a2  = 2.0 * a;
    sq  = sqrt( b*b - 4.0*a*c );
    x1  = ( -b + sq ) / a2;
    x2  = ( -b - sq ) / a2;
}
