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
 *                               inverse-impl.h
 *
 * Implementation for matrix inverse
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Cpmpute the inverse of square matrix.
 */
template <typename Type>
Matrix<Type> inv( const Matrix<Type> &A, const string &type )
{
    int n = A.rows();
    assert( n == A.cols() );

    if( type == "spd" )
    {
        Cholesky<Type> cho;
        cho.dec(A);
        if( cho.isSpd() )
            return cho.solve( eye(n,Type(1)) );
        else
        {
            cerr << "The matrix is not symmetric!" << endl;
            return A;
        }
    }
    else
    {
        LUD<Type> lu;
        lu.dec(A);
        return lu.solve( eye(n,Type(1)) );
    }
}


/**
 * Gauss-jordan column pivot elimination for computing matrix's inverse.
 * The matrix can be both REAL or COMPLEX.
 */
template <typename Type>
Matrix<Type> colPivInv( const Matrix<Type> &A )
{
    int rows = A.rows();
    int clumns = A.cols();

    assert( rows == clumns );

    Matrix<Type> invA(A);
    Vector<int> index( rows );
    int i, j, k;
    Type tmp = 0;

    for( k=0; k<rows; ++k )
    {
        //Findint pivot and exchange if necessary.
        index[k] = k;
        Type mvl = invA[k][k];
        for( i=k+1; i<rows; ++i )
        {
            tmp = abs(invA[i][k]);
            if( abs(tmp) > abs(mvl) )
            {
                mvl = tmp;
                index[k] = i;
            }
        }
        if( abs(mvl) < EPS )
        {
            cerr << "\n" << "A is a singular matrix." << "\n";
            return Matrix<Type>(0,0);
        }

        if( index[k] != k )
        {
            tmp = 0;
            for( j=0; j<rows; ++j )
            {
                tmp = invA[k][j];
                invA[k][j] = invA[index[k]][j];
                invA[index[k]][j] = tmp;
            }
        }

        // Calculating the kth column.
        invA[k][k] = Type(1) / invA[k][k];
        for( i=0; i<rows; ++i )
            if( i != k )
                invA[i][k] = - invA[k][k]*invA[i][k];

        // Calculating all elements excptint the kth row and column.
        for( i=0; i<rows; ++i )
            if( i != k )
                for( j=0; j<rows; ++j )
                    if( j != k )
                        invA[i][j] += invA[i][k] * invA[k][j];

        // Calculating the kth row.
        for( j=0; j<rows; ++j )
            if( j != k )
                invA[k][j] *= invA[k][k];
    }

    //Exchanging back.
    for( k=rows-1; k>=0; --k )
        if( index[k] != k )
        {
            tmp = 0;
            for( i=0; i<rows; ++i )
            {
                tmp = invA[i][k];
                invA[i][k] = invA[i][index[k]];
                invA[i][index[k]] = tmp;
            }
        }

    return invA;
}


/**
 * Gauss-jordan complete pivot elimination for computing matrix's inverse.
 * The matrix can be both REAL or COMPLEX.
 */
template <typename Type>
Matrix<Type> cmpPivInv( const Matrix<Type> &A )
{
    int n = A.rows();
    assert( n == A.cols() );
    Matrix<Type> invA(A);

    int k;
    Type tmp = 0,
         mvl = 0;
    Vector<int> rowIndex(n),
                colIndex(n);

    for( k=0; k<n; ++k )
    {
        //finding pivot
        mvl = 0;
        for( int i=k; i<n; ++i )
            for( int j=k; j<n; ++j )
            {
                tmp = abs(invA[i][j]);
                if( abs(invA[i][j]) > abs(mvl) )
                {
                    mvl = tmp;
                    rowIndex[k] = i;
                    colIndex[k] = j;
                }
            }

        if( abs(mvl) < EPS )
        {
           cerr << endl << "A is a singular matrix." << endl;
            return invA;
        }

        // row exchange
        if( rowIndex[k] != k )
            for( int i=0; i<n; ++i )
                swap( invA[k][i], invA[rowIndex[k]][i] );

        // column exchange
        if( colIndex[k] != k )
            for( int j=0; j<n; ++j )
                swap( invA[j][k], invA[j][rowIndex[k]] );

        // Calculating the kth column.
        invA[k][k] = Type(1) / invA[k][k];
        for( int j=0; j<n; ++j )
            if( j != k )
                invA[k][j] = invA[k][j]*invA[k][k];

        // Calculating all elements excptint the kth row and column.
        for( int i=0; i<n; ++i )
            if( i != k )
                for( int j=0; j<n; ++j )
                    if( j != k )
                        invA[i][j] -= invA[i][k]*invA[k][j];

        for( int i=0; i<n; ++i )
            if( i != k )
                invA[i][k] = -invA[i][k]*invA[k][k];
    }

    //Exchanging back.
    for( k=n-1; k>=0; --k )
    {
        if( colIndex[k] != k )
            for( int j=0; j<n; ++j )
                swap( invA[k][j], invA[colIndex[k]][j] );
        if( rowIndex[k] != k )
            for( int i=0; i<n; ++i )
                swap( invA[i][k], invA[i][rowIndex[k]] );
    }

    return invA;
}


/**
 * Cpmpute the inverse of square complex matrix. This algorithm translate
 * complex matrix inverse to real matrix inverse, but the amount of
 * calculation is very great.
 */
template <typename Type>
Matrix<complex<Type> > cinv( const Matrix<complex<Type> > &M,
                             const string &type )
{
    int n = M.rows();
    assert( n == M.cols() );

    Matrix<Type> A = real(M),
                 B = imag(M),
                 C = inv( A + B*inv(A,type)*B ),
                 D = -inv( B + A*inv(B,type)*A );

    Matrix<complex<Type> > IM(n,n);
    for( int i=0; i<n; ++i )
        for( int j=0; j<n; ++j )
            IM[i][j] = complex<Type>( C[i][j], D[i][j] );

    return IM;
}
