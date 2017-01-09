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
 *                               csvd-impl.h
 *
 * Implementation for CSVD class.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
CSVD<Type>::CSVD()
{
}

template<typename Type>
CSVD<Type>::~CSVD()
{
}


/**
 * Making singular decomposition.
 */
template <typename Type>
void CSVD<Type>::dec( const Matrix< complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols(),
        p = min( m, n );

    U.resize( m, p );
    V.resize( n, p );
    S.resize( p );
    if( m >= n )
    {
        Matrix< complex<Type> > B(A);
        decomposition( B, U, S, V );
    }
    else
    {
        Matrix< complex<Type> > B( trH( A ) );
        decomposition( B, V, S, U );
    }
}


/**
 * The singular value decomposition of a complex M by N (N<=M)
 * matrix A has the form A = U S V*, where:
 *     U is an M by N unitary matrix (for economy size),
 *     S is an M by N diagonal matrix,
 *     V is an N by N unitary matrix.
 *
 * Input:
 *     a ---> the complex matrix to be decomposition.
 *     m ---> the number of rows in a.
 *     n ---> the number of columns in a.
 *
 * Output:
 *     u ---> the left singular vectors.
 *     s ---> the singular values.
 *     v ---> the right singular vectors.
 */
template <typename Type>
void CSVD<Type>::decomposition( Matrix< complex<Type> > &B,
                                Matrix< complex<Type> > &U,
                                Vector<Type> &S,
                                Matrix< complex<Type> > &V )
{
    int k, k1,
        L, L1,
        nM1,
        m = B.rows(),
        n = B.cols();

    Type    tol, eta;
	Type    w, x, y, z;
	Type    cs, sn, f, g, h;
    complex<Type>   q;

    Vector<Type>    b(m), c(m), t(m);

    L = 0;
    nM1 = n - 1;
	eta = 2.8E-16;      // the relative machine precision

	// Householder Reduction
	c[0] = 0;
	k = 0;
	while( 1 )
	{
		k1 = k + 1;

		// elimination of a[i][k], i = k, ..., m-1
		z = 0;
		for( int i=k; i<m; ++i )
			z += norm(B[i][k]);
		b[k] = 0;

		if( z > EPS )
		{
			z = sqrt(z);
			b[k] = z;
			w = abs(B[k][k]);
			q = 1;
			if( w != Type(0) )
                q = B[k][k] / w;

            B[k][k] = q * ( z + w );

			if( k != n - 1 )
                for( int j=k1; j<n; ++j )
				{
                    q = 0;
					for( int i=k; i<m; ++i )
                        q += conj(B[i][k]) * B[i][j];
                    q /= z * ( z + w );

					for( int i=k; i<m; ++i )
                        B[i][j] -= q * B[i][k];
				}

			// phase transformation
			q = -conj(B[k][k]) / abs(B[k][k]);
			for( int j=k1; j<n; ++j )
				B[k][j] *= q;
		}

		// elimination of a[k][j], j = k+2, ..., n-1
		if( k == nM1 )
            break;
		z = 0;
		for( int j=k1; j<n; ++j )
            z += norm( B[k][j] );
		c[k1] = 0;

		if( z > EPS )
		{
			z = sqrt(z);
			c[k1] = z;
			w = abs( B[k][k1] );
			q = 1;
			if( w != Type(0) )
				q = B[k][k1] / w;

			B[k][k1] = q * ( z + w );
			for( int i=k1; i<m; ++i )
			{
				q = 0;
				for( int j=k1; j<n; ++j )
					q += conj(B[k][j]) * B[i][j];

				q /= z * ( z + w );
				for( int j=k1; j<n; ++j )
					B[i][j] -= q * B[k][j];
			}

			// phase transformation
			q = -conj(B[k][k1]) / abs(B[k][k1]);
			for( int i=k1; i<m; ++i )
                B[i][k1] *= q;
		}
		k = k1;
	}

	// tolerance for negligible elements
	tol = 0;
	for( k=0; k<n; ++k )
	{
		S[k] = b[k];
		t[k] = c[k];
		if( S[k]+t[k] > tol )
            tol = S[k] + t[k];
	}
	tol *= eta;

	// initialization fo U and V
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<m; ++i )
            U[i][j] = 0;
        U[j][j] = 1;

        for( int i=0; i<n; ++i )
            V[i][j] = 0;
        V[j][j] = 1;
    }

	// QR diagonallization
	for( k = nM1; k >= 0; --k )
	{
		// test for split
		while( 1 )
		{
			for( L=k; L>=0; --L )
            {
				if( abs(t[L]) <= tol )
                    goto Test;
				if( abs(S[L - 1]) <= tol )
                    break;
			}

			// cancellation of E(L)
			cs = 0;
			sn = 1;
			L1 = L - 1;
			for( int i=L; i<=k; ++i )
			{
				f = sn * t[i];
				t[i] *= cs;
				if( abs(f) <= tol )
                    goto Test;

				h = S[i];
				w = sqrt( f*f + h*h );
				S[i] = w;
				cs = h / w;
				sn = -f / w;

                for( int j=0; j<n; ++j )
                {
                    x = real( U[j][L1] );
                    y = real( U[j][i] );
                    U[j][L1] = x*cs + y*sn;
                    U[j][i]  = y*cs - x*sn;
                }
			}

			// test for convergence
	Test:	w = S[k];
			if( L == k )
                break;

			// origin shift
			x = S[L];
			y = S[k-1];
			g = t[k-1];
			h = t[k];
			f = ( (y-w)*(y+w) + (g-h)*(g+h) ) / ( 2*h*y);
			g = sqrt( f*f + 1 );
			if( f < Type(0) )
                g = -g;
			f = ( (x-w)*(x+w) + h*(y/(f+g)-h) ) / x;

			// QR step
			cs = 1;
			sn = 1;
			L1 = L + 1;
			for( int i=L1; i<=k; ++i )
			{
				g = t[i];
				y = S[i];
				h = sn * g;
				g = cs * g;
				w = sqrt( h*h + f*f );
				t[i-1] = w;
				cs = f / w;
				sn = h / w;
				f = x * cs + g * sn;
				g = g * cs - x * sn;
				h = y * sn;
				y = y * cs;

                for( int j=0; j<n; ++j )
                {
                    x = real( V[j][i-1] );
                    w = real( V[j][i] );
                    V[j][i-1] = x*cs + w*sn;
                    V[j][i]   = w*cs - x*sn;
                }

				w = sqrt( h*h + f*f );
				S[i-1] = w;
				cs = f / w;
				sn = h / w;
				f = cs * g + sn * y;
				x = cs * y - sn * g;

                for( int j=0; j<n; ++j )
                {
                    y = real(U[j][i-1]);
                    w = real(U[j][i]);
                    U[j][i-1] = y*cs + w*sn;
                    U[j][i]   = w*cs - y*sn;
                }
			}
			t[L] = 0;
			t[k] = f;
			S[k] = x;
		}

		// convergence
		if( w >= Type(0) )
            continue;
		S[k] = -w;

		for( int j=0; j<n; ++j )
			V[j][k] = -V[j][k];
	}

	// sort dingular values
	for( k=0; k<n; ++k )
	{
		g = -1.0;
		int j = k;
		for( int i=k; i<n; ++i )
		{
			if( S[i] <= g )
                continue;
			g = S[i];
			j = i;
		}

		if( j == k )
            continue;
		S[j] = S[k];
		S[k] = g;

        for( int i=0; i<n; ++i )
        {
            q = V[i][j];
            V[i][j] = V[i][k];
            V[i][k] = q;
        }

        for( int i=0; i<n; ++i )
        {
            q = U[i][j];
            U[i][j] = U[i][k];
            U[i][k] = q;
        }
	}

	// back transformation
    for( k=nM1; k>=0; --k )
    {
        if( b[k] == Type(0) )
            continue;
        q = -B[k][k] / abs(B[k][k]);

        for( int j=0; j<n; ++j )
            U[k][j] *= q;

        for( int j=0; j<n; ++j )
        {
            q = 0;
            for( int i=k; i<m; ++i )
                q += conj(B[i][k]) * U[i][j];
            q /= abs(B[k][k]) * b[k];

            for( int i=k; i<m; ++i )
                U[i][j] -= q * B[i][k];
        }
    }

	if( n > 1 )
	{
		for( k = n-2; k>=0; --k )
		{
			k1 = k + 1;
			if( c[k1] == Type(0) )
                continue;
			q = -conj(B[k][k1]) / abs(B[k][k1]);

			for( int j=0; j<n; ++j )
				V[k1][j] *= q;

			for( int j=0; j<n; ++j )
			{
				q = 0;
				for( int i=k1; i<n; ++i )
					q += B[k][i] * V[i][j];
				q /= abs(B[k][k1]) * c[k1];

				for( int i=k1; i<n; ++i )
					V[i][j] -= q * conj(B[k][i]);
			}
		}
	}
}


/**
 * Get the left singular vectors.
 */
template<typename Type>
inline Matrix< complex<Type> > CSVD<Type>::getU() const
{
    return U;
}


/**
 * Get the singular values matrix.
 */
template<typename Type>
inline Matrix<Type> CSVD<Type>::getSM()
{
    int N = S.size();
    Matrix<Type> tmp( N, N );
    for( int i=0; i<N; ++i )
        tmp[i][i] = S[i];

    return tmp;
}


/**
 * Get the singular values vector.
 */
template<typename Type>
inline Vector<Type> CSVD<Type>::getSV() const
{
    return S;
}


/**
 * Get the right singular vectors.
 */
template<typename Type>
inline Matrix< complex<Type> > CSVD<Type>::getV() const
{
    return V;
}


/**
 * Two norm (max(S)).
 */
template <typename Type>
inline Type CSVD<Type>::norm2() const
{
    return S[0];
}


/**
 * Two norm of condition number (max(S)/min(S)).
 */
template <typename Type>
inline Type CSVD<Type>::cond() const
{
    return ( S[0] / S(S.size()) );
}


/**
 * Effective numerical matrix rank.
 */
template <typename Type>
int CSVD<Type>::rank()
{
    int N = S.size();
    double tol = N * S[0] * EPS;
    int r = 0;

    for( int i=0; i<N; ++i )
        if( S[i] > tol )
            r++;

    return r;
}
