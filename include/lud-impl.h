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
 *                               lud-impl.h
 *
 * Implementation for LUD class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
LUD<Type>::LUD()
{
}

template<typename Type>
LUD<Type>::~LUD()
{
}


/**
 * permute copy
 */
template <typename Type>
Vector<Type> LUD<Type>::permuteCopy( const Vector<Type> &A,
                                     const Vector<int> &piv )
{
    int pivLength = piv.dim();
    if( pivLength != A.dim() )
        return Vector<Type>();

    Vector<Type> x(pivLength);
    for( int i=0; i<pivLength; ++i )
        x[i] = A[piv[i]];

    return x;
}

template <typename Type>
Matrix<Type> LUD<Type>::permuteCopy( const Matrix<Type> &A,
                                     const Vector<int> &piv, int j0, int j1 )
{
    int pivLength = piv.dim();
    Matrix<Type> X( pivLength, j1-j0+1 );

    for( int i=0; i<pivLength; ++i )
        for( int j=j0; j<=j1; ++j )
            X[i][j-j0] = A[piv[i]][j];

    return X;
}


/**
 * Return pivot permutation vector
 */
template<typename Type>
inline Vector<int> LUD<Type>::getPivot() const
{
    return piv;
}


/**
 * LU Decomposition
 */
template <typename Type>
void LUD<Type>::dec(const Matrix<Type> &A)
{
    m = A.rows();
    n = A.cols();
    piv.resize(m);
    LU = A;

    // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
    for( int i=0; i<m; ++i )
        piv[i] = i;

    pivsign = 1;
    Type *LUrowi = 0;
    Vector<Type> LUcolj(m);

    // outer loop
    for( int j=0; j<n; ++j )
    {
        // Make a copy of the j-th column to localize references.
        for( int i=0; i<m; ++i )
            LUcolj[i] = LU[i][j];

        // Apply previous transformations.
        for( int i=0; i<m; ++i )
        {
            LUrowi = LU[i];

            // Most of the time is spent in the following dot product.
            int kmax = (i < j)? i : j;
            Type s = 0;

            for( int k=0; k<kmax; ++k )
                s += LUrowi[k]*LUcolj[k];

            LUrowi[j] = LUcolj[i] -= s;
        }

        // Find pivot and exchange if necessary.
        int p = j;
        for( int i=j+1; i<m; ++i )
            if( abs(LUcolj[i]) > abs(LUcolj[p]) )
                p = i;

        if( p != j )
        {
            int k=0;
            for( k=0; k<n; ++k )
                swap( LU[p][k], LU[j][k] );

            swap( piv[p], piv[j] );
            pivsign = -pivsign;
        }

        // compute multipliers
        if( (j < m) && ( abs(LU[j][j]) != 0 ) )
            for( int i=j+1; i<m; ++i )
                LU[i][j] /= LU[j][j];
    }
}


/**
 * Return lower triangular matrix L.
 */
template <typename Type>
Matrix<Type> LUD<Type>::getL()
{
    int p = min( m, n );
    Matrix<Type> tmp( m, p );

    for( int i=0; i<m; ++i )
        for( int j=0; j<i && j<n; ++j )
            tmp[i][j] = LU[i][j];
   for( int i=0; i<p; ++i )
        tmp[i][i] = 1;

    return tmp;
}


/**
 * Return upper triangular matrix U.
 */
template <typename Type>
Matrix<Type>LUD<Type>::getU()
{
    int p = min( m, n );
    Matrix<Type> tmp( p, n );

    for( int i=0; i<m; ++i )
        for( int j=i; j<n; ++j )
            tmp[i][j] = LU[i][j];

    return tmp;
}


/**
 * Compute determinant using LU factors.
 */
template <typename Type>
Type LUD<Type>::det()
{
    if( m != n )
        return 0;

    Type d = Type(pivsign);
    for( int j=0; j<n; ++j )
        d *= LU[j][j];

    return d;
}


/**
 * true if upper triangular factor U is nonsingular, 0 otherwise.
 */
template <typename Type>
inline bool LUD<Type>::isNonsingular()
{
    for( int j=0; j<n; ++j )
        if( abs(LU[j][j]) == 0 )
            return false;

    return true;
}


/**
 * Solve A*x = b, where x and b are vectors of length equal
 * to the number of rows in A.
 */
template <typename Type>
Vector<Type> LUD<Type>::solve( const Vector<Type> &b )
{
    // dimensions: A is mxn, X is nxk, B is mxk
    if( b.dim() != m )
        return Vector<Type>();

    if( !isNonsingular() )
        return Vector<Type>();

    Vector<Type> x = permuteCopy( b, piv );

    // solve L*Y = B(piv)
    for( int k=0; k<n; ++k )
        for( int i=k+1; i<n; ++i )
            x[i] -= x[k]*LU[i][k];

    // solve U*x = y;
    for( int k=n-1; k>=0; --k )
    {
        x[k] /= LU[k][k];
        for( int i=0; i<k; ++i )
            x[i] -= x[k]*LU[i][k];
    }

    return x;
}


/**
 * Solve A*X = B
 */
template <typename Type>
Matrix<Type> LUD<Type>::solve( const Matrix<Type> &B )
{
    // dimensions: A is mxn, X is nxk, B is mxk
    if( B.rows() != m )
        return Matrix<Type>(0,0);

    if( !isNonsingular() )
        return Matrix<Type>(0,0);

    // copy right hand side with pivoting
    int nx = B.cols();
    Matrix<Type> X = permuteCopy( B, piv, 0, nx-1 );

    // solve L*Y = B(piv,:)
    for( int k=0; k<n; ++k )
        for( int i=k+1; i<n; ++i )
            for( int j=0; j<nx; ++j )
                X[i][j] -= X[k][j]*LU[i][k];

    // solve U*X = Y;
    for( int k=n-1; k>=0; --k )
    {
        for( int j=0; j<nx; ++j )
            X[k][j] /= LU[k][k];

        for( int i=0; i<k; ++i )
            for( int j=0; j<nx; ++j )
                X[i][j] -= X[k][j]*LU[i][k];
    }

    return X;
}
