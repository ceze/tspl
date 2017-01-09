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
 *                               matrix-impl.h
 *
 * Implementation for Matrix class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * initialize
 */
template <typename Type>
void Matrix<Type>::init( int rows, int columns )
{
	nRow = rows;
	nColumn = columns;
	nTotal = nRow * nColumn;

	pv0 = new Type[nTotal];
	prow0 = new Type*[nRow];
	prow1 = new Type*[nRow];

	assert( pv0 != NULL );
	assert( prow0 != NULL );
	assert( prow1 != NULL );

	Type *p = pv0;
	pv1 = pv0 - 1;
	for( int i=0; i<nRow; ++i )
	{
		prow0[i] = p;
		prow1[i] = p-1;
		p += nColumn;
	}

	prow1--;
}


/**
 * copy matrix from normal array
 */
template <typename Type>
inline void Matrix<Type>::copyFromArray( const Type *v )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = v[i];
}


/**
 * set matrix by a scalar
 */
template <typename Type>
inline void Matrix<Type>::setByScalar( const Type &x )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = x;
}


/**
 * destroy the matrix
 */
template <typename Type>
void Matrix<Type>::destroy()
{
	if( pv0 == NULL )
		return ;
	else
		delete []pv0;

	if( prow0 != NULL )
		delete []prow0;

	prow1++;
	if( prow1 != NULL )
		delete []prow1;
}


/**
 * constructors and destructor
 */
template <typename Type>
Matrix<Type>::Matrix()
: pv0(0), pv1(0), prow0(0), prow1(0), nRow(0), nColumn(0), nTotal(0)
{
}

template <typename Type>
Matrix<Type>::Matrix( const Matrix<Type> &A )
{
	init( A.nRow, A.nColumn );
	copyFromArray( A.pv0 );
}

template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type &x )
{
	init( rows,columns );
	setByScalar(x);
}

template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type *arrays )
{
	init( rows,columns );
	copyFromArray( arrays );
}

template <typename Type>
Matrix<Type>::~Matrix()
{
	destroy();
}


/**
 * overload evaluate operator = from matrix to matrix
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator=( const Matrix<Type> &A )
{
	if( pv0 == A.pv0 )
		return *this;

	if( nRow == A.nRow && nColumn == A.nColumn )
		copyFromArray( A.pv0 );
	else
	{
		destroy();
		init( A.nRow, A.nColumn );
		copyFromArray( A.pv0 );
	}

	return *this;
}


/**
 * overload evaluate operator = from scalar to matrix
 */
template <typename Type>
inline Matrix<Type>& Matrix<Type>::operator=( const Type &x )
{
	setByScalar( x );

	return *this;
}


/**
 * overload operator [] for 0-offset access
 */
template <typename Type>
inline Type* Matrix<Type>::operator[]( int i )
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return prow0[i];
}

template <typename Type>
inline const Type* Matrix<Type>::operator[]( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return prow0[i];
}


/**
 * overload operator () for 1-offset access
 */
template <typename Type>
inline Type& Matrix<Type>::operator()( int row, int column )
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= nRow ) ;
	assert( 1 <= column);
	assert( column <= nColumn );
#endif

	return  prow1[row][column];
}

template <typename Type>
inline const Type& Matrix<Type>::operator()( int row, int column ) const
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= nRow ) ;
	assert( 1 <= column);
	assert( column <= nColumn );
#endif

	return  prow1[row][column];
}


/**
 * type conversion functions
 */
template <typename Type>
inline Matrix<Type>::operator Type*()
{
	return pv0;
}

template <typename Type>
inline Matrix<Type>::operator const Type*() const
{
	return pv0;
}

//template <typename Type>
//inline Matrix<Type>::operator Type**()
//{
//	return prow0;
//}
//
//template <typename Type>
//inline Matrix<Type>::operator const Type**() const
//{
//	return prow0;
//}


/**
 * get the matrix's size
 */
template <typename Type>
inline long Matrix<Type>::size() const
{
	return nTotal;
}


/**
 * get the matrix's dimension
 */
template <typename Type>
int Matrix<Type>::dim( int dimension ) const
{
#ifdef BOUNDS_CHECK
	assert( dimension >= 1);
	assert( dimension <= 2);
#endif

	if( dimension == 1 )
		return nRow;
	else if( dimension == 2 )
		return nColumn;
	else
		return 0;
}

template <typename Type>
inline int Matrix<Type>::rows() const
{
    return nRow;
}

template <typename Type>
inline int Matrix<Type>::cols() const
{
    return nColumn;
}


/**
 * reallocate matrix's size
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::resize( int rows, int columns )
{
	if(  rows == nRow && columns == nColumn )
		return *this;

	destroy();
	init( rows, columns );

	return *this;
}


/**
 * get the matrix's row vector
 */
template <typename Type>
Vector<Type> Matrix<Type>::getRow( int row ) const
{
#ifdef BOUNDS_CHECK
	assert( row >= 0 );
	assert( row < nRow );
#endif

	Vector<Type> tmp( nColumn );
	for( int j=0; j<nColumn; ++j )
		tmp[j] = prow0[row][j];

	return tmp;
}


/**
 * get the matrix's column vector
 */
template <typename Type>
Vector<Type> Matrix<Type>::getColumn( int column ) const
{
#ifdef BOUNDS_CHECK
	assert( column >= 0 );
	assert( column < nColumn );
#endif

	Vector<Type> tmp( nRow );
	for( int i=0; i<nRow; ++i )
		tmp[i] = prow0[i][column];

	return tmp;
}


/**
 * set the matrix's row vector
 */
template <typename Type>
void Matrix<Type>::setRow( const Vector<Type> &v, int row )
{
#ifdef BOUNDS_CHECK
	assert( row >= 0 );
	assert( row < nRow );
	assert( v.dim() == nColumn );
#endif

	for( int j=0; j<nColumn; ++j )
		prow0[row][j] = v[j];
}


/**
 * set the matrix's column vector
 */
template <typename Type>
void Matrix<Type>::setColumn( const Vector<Type> &v, int column )
{
#ifdef BOUNDS_CHECK
	assert( column >= 0 );
	assert( column < nColumn );
	assert( v.dim() == nRow );
#endif

	for( int i=0; i<nRow; ++i )
		prow0[i][column] = v[i];
}


/**
 * compound assignment operators +=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ += x;
    }

	return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Matrix<Type> &rhs )
{
    assert( nRow == rhs.rows() );
    assert( nColumn == rhs.cols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ += *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators -=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ -= x;
    }

	return *this;
}

template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Matrix<Type> &rhs )
{
    assert( nRow == rhs.rows() );
    assert( nColumn == rhs.cols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ -= *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators *=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ *= x;
    }

	return *this;
}

// WARNING: this is element-by-element multiplication
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Matrix<Type> &rhs )
{
    assert( nRow == rhs.rows() );
    assert( nColumn == rhs.cols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ *= *colPtrR++;
    }

	return *this;
}


/**
 * compound assignment operators /=
 */
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ /= x;
    }

	return *this;
}

// WARNING: this is element-by-element division
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Matrix<Type> &rhs )
{
    assert( nRow == rhs.rows() );
    assert( nColumn == rhs.cols() );

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ /= *colPtrR++;
    }

	return *this;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
ostream& operator<<( ostream &out, const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}


/**
 * Overload the intput stream function.
 */
template <typename Type>
istream& operator>>( istream &in, Matrix<Type> &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.rows() && columns == A.cols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			in >> A[i][j];

	return in;
}


/**
 * get negative matrix
 */
template<typename Type>
Matrix<Type> operator-( const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = -A[i][j];

	return tmp;
}


/**
 * matrix-scalar addition
 */
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp += x;
}

template<typename Type>
inline Matrix<Type> operator+( const Type &x, const Matrix<Type> &A )
{
	return A + x;
}


/**
 * matrix-matrix addition
 */
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp += A2;
}


/**
 * matrix-scalar subtraction
 */
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp -= x;
}

template<typename Type>
inline Matrix<Type> operator-( const Type &x, const Matrix<Type> &A )
{
	Matrix<Type> tmp( A );
	return -tmp += x;
}


/**
 * matrix-matrix subtraction
 */
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp -= A2;
}


/**
 * matrix-scaling multiplication
 */
template <typename Type>
inline Matrix<Type> operator*( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp *= x;
}

template <typename Type>
inline Matrix<Type> operator*( const Type &x, const Matrix<Type> &A )
{
	return A * x;
}


/**
 * matrix-matrix multiplication
 */
template <typename Type>
Matrix<Type> operator*( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.cols() == A2.rows() );

	int rows = A1.rows();
	int columns = A2.cols();
//	int K = A1.cols();

	Matrix<Type> tmp( rows, columns );
//	for( int i=0; i<rows; ++i )
//		for( int j=0; j<columns; ++j )
//		{
//            tmp[i][j] = 0;
//			for( int k=0; k<K; ++k )
//			    tmp[i][j] += A1[i][k] * A2[k][j];
//		}

    mult( A1, A2, tmp );

	return tmp;
}


/**
 * matrix-vector multiplication
 */
template <typename Type>
Vector<Type> operator*( const Matrix<Type> &A, const Vector<Type> &b )
{
	assert( A.cols() == b.dim() );

	int rows = A.rows();
//	int columns = A.cols();

	Vector<Type> tmp(rows);
//	for( int i=0; i<rows; ++i )
//	{
//		Type sum = 0;
//		for( int j=0; j<columns; ++j )
//			sum += A[i][j] * v[j];
//		tmp[i] = sum;
//	}

    mult( A, b, tmp );

	return tmp;
}


/**
 * matrix-scalar division
 */
template <typename Type>
inline Matrix<Type> operator/( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp /= x;
}

template <typename Type>
Matrix<Type> operator/( const Type &x, const Matrix<Type> &A )
{
	int rows = A.rows();
	int clumns = A.cols();

	Matrix<Type> tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}


/**
 * This is an optimized version of matrix multiplication,
 * where the destination matrix has already been allocated.
 */
template <typename Type>
Matrix<Type>& mult( const Matrix<Type> &A, const Matrix<Type> &B,
                    Matrix<Type> &C )
{
    int M = A.rows();
    int N = B.cols();
    int K = A.cols();

    assert( B.rows() == K );

    C.resize( M, N );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
        for( int j=0; j<N; ++j )
        {
            pRow  = &A[i][0];
            pCol  = &B[0][j];
            sum = 0;

            for( int k=0; k<K; ++k )
            {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol += N;
            }
            C[i][j] = sum;
        }
    return C;
}


/**
 * This is an optimized version of matrix and vector multiplication,
 * where the destination vector has already been allocated.
 */
template <typename Type>
Vector<Type>& mult( const Matrix<Type> &A, const Vector<Type> &b,
                    Vector<Type> &c )
{
    int M = A.rows();
    int N = A.cols();

    assert( b.size() == N );

    c.resize( M );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
    {
        pRow  = &A[i][0];
        pCol  = &b[0];
        sum = 0;

        for( int j=0; j<N; ++j )
        {
            sum += (*pRow) * (*pCol);
            pRow++;
            pCol++;
        }
        c[i] = sum;
    }
    return c;
}


/**
 * matrix-matrix elementwise multiplication
 */
template<typename Type>
inline Matrix<Type> elemMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp *= A2;
}

template <typename Type>
inline Matrix<Type>& elemMultEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 *= A2;
}


/**
 * matrix-matrix elementwise division
 */
template <typename Type>
inline Matrix<Type> elemDivd( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp /= A2;
}

template <typename Type>
inline Matrix<Type>& elemDivdEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 /= A2;
}


/**
 * matrix tranpose
 */
template <typename Type>
Matrix<Type> trT( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = A[j][i];

	return tmp;
}


/**
 * matrix conjugate tranpose
 */
template <typename Type>
Matrix<Type> trH( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = conj(A[j][i]);

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A^T * B.
 */
template <typename Type>
Matrix<Type> trMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.rows() == A2.rows() );

	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<Type> tmp( rows, columns );
//	for( int i=0; i<rows; ++i )
//		for( int j=0; j<columns; ++j )
//		{
//			Type sum = 0;
//			for( int k=0; k<K; ++k )
//			   sum += A1[k][i] * A2[k][j];
//			tmp[i][j] = sum;
//		}
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[k][i] * A2[k][j];

	return tmp;
}


/**
 * matrix-vector tranpose multiplication: A^T * b.
 */
template <typename Type>
Vector<Type> trMult( const Matrix<Type> &A, const Vector<Type> &v )
{
	assert( A.rows() == v.dim() );

	int rows = A.rows();
	int columns = A.cols();

	Vector<Type> tmp( columns );
//	for( int i=0; i<columns; ++i )
//	{
//		Type sum = 0;
//		for( int j=0; j<rows; ++j )
//			sum += A[j][i] * v[j];
//		tmp[i] = sum;
//	}
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += A[j][i] * v[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A * B^T.
 */
template <typename Type>
Matrix<Type> multTr( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	assert( A1.cols() == A2.cols() );

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<Type> tmp( rows, columns );
//	for( int i=0; i<rows; ++i )
//		for( int j=0; j<columns; ++j )
//		{
//			Type sum = 0;
//			for( int k=0; k<K; ++k )
//			   sum += A1[i][k] * A2[j][k];
//			tmp[i][j] = sum;
//		}
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * A2[j][k];

	return tmp;
}


/**
 * vector-vector tranpose multiplication: a * b^T.
 */
template <typename Type>
Matrix<Type> multTr( const Vector<Type> &a, const Vector<Type> &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*b[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A^H * B.
 */
template <typename Type>
Matrix<complex<Type> > trMult( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{
	assert( A1.rows() == A2.rows() );

	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<complex<Type> > tmp( rows, columns );
//	for( int i=0; i<rows; ++i )
//		for( int j=0; j<columns; ++j )
//		{
//			Type sum = 0;
//			for( int k=0; k<K; ++k )
//			   sum += A1[k][i] * A2[k][j];
//			tmp[i][j] = sum;
//		}
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += conj(A1[k][i]) * A2[k][j];

	return tmp;
}


/**
 * matrix-vector tranpose multiplication: A^H * b.
 */
template <typename Type>
Vector<complex<Type> > trMult( const Matrix<complex<Type> > &A,
                               const Vector<complex<Type> > &v )
{
	assert( A.rows() == v.dim() );

	int rows = A.rows();
	int columns = A.cols();

	Vector<complex<Type> > tmp( columns );
//	for( int i=0; i<columns; ++i )
//	{
//		Type sum = 0;
//		for( int j=0; j<rows; ++j )
//			sum += A[j][i] * v[j];
//		tmp[i] = sum;
//	}
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += conj(A[j][i]) * v[j];

	return tmp;
}


/**
 * matrix-matrix tranpose multiplication: A * B^H.
 */
template <typename Type>
Matrix<complex<Type> > multTr( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{
	assert( A1.cols() == A2.cols() );

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<complex<Type> > tmp( rows, columns );
//	for( int i=0; i<rows; ++i )
//		for( int j=0; j<columns; ++j )
//		{
//			Type sum = 0;
//			for( int k=0; k<K; ++k )
//			   sum += A1[i][k] * A2[j][k];
//			tmp[i][j] = sum;
//		}
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * conj(A2[j][k]);

	return tmp;
}


/**
 * vector-vector tranpose multiplication: a * b^H.
 */
template <typename Type>
Matrix<complex<Type> > multTr( const Vector<complex<Type> > &a,
                               const Vector<complex<Type> > &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<complex<Type> > tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*conj(b[j]);

	return tmp;
}


/**
 * Generate the identity matrix.
 */
template <typename Type>
Matrix<Type> eye( int N, const Type &x )
{
    Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = x;

	return tmp;
}


/**
 * Get the diagonal entries of matrix.
 */
template <typename Type>
Vector<Type> diag( const Matrix<Type> &A )
{
	int nColumn = A.rows();
	if( nColumn > A.cols() )
		nColumn = A.cols();

	Vector<Type> tmp( nColumn );
	for( int i=0; i<nColumn; ++i )
		tmp[i] = A[i][i];

	return tmp;
}


/**
 * Generate the diagonal of matrix by given its diagonal elements.
 */
template <typename Type>
Matrix<Type> diag( const Vector<Type> &d )
{
	int N = d.size();

	Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = d[i];

	return tmp;
}


/**
 * Compute Frobenius norm of matrix.
 */
template <typename Type>
Type norm( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += A(i,j) * A(i,j);

	return sqrt(sum);
}

template <typename Type>
Type norm( const Matrix<complex<Type> > &A )
{
	int m = A.rows();
	int n = A.cols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += norm(A(i,j));

	return sqrt(sum);
}


/**
 * Swap two matrixes.
 */
template <typename Type> void swap( Matrix<Type> &lhs, Matrix<Type> &rhs )
{
    int m = lhs.rows();
	int n = lhs.cols();

	assert( m == rhs.rows() );
	assert( n == rhs.cols() );

	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            swap( lhs(i,j), rhs(i,j) );
}


/**
 * Matrix's column vecotrs sum.
 */
template <typename Type>
Vector<Type> sum( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
		for( int i=1; i<=m; ++i )
            sum(j) += A(i,j);

	return sum;
}


/**
 * Minimum of matrix's column vecotrs.
 */
template <typename Type>
Vector<Type> min( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<m; ++i )
            if( tmp > A(i,j) )
                tmp = A(i,j);
        sum(j) = tmp;
	}

	return sum;
}


/**
 * Maximum of matrix's column vecotrs.
 */
template <typename Type>
Vector<Type> max( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	Vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<m; ++i )
            if( tmp < A(i,j) )
                tmp = A(i,j);
        sum(j) = tmp;
	}

	return sum;
}


/**
 * Matrix's column vecotrs mean.
 */
template <typename Type>
inline Vector<Type> mean( const Matrix<Type> &A )
{
	return sum(A) / Type(A.rows());
}


/**
 * Convert real matrix to complex matrix.
 */
template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &rA )
{
	int rows = rA.rows();
	int columns = rA.cols();

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = rA[i][j];

    return cA;
}

template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &mR,
                                      const Matrix<Type> &mI )
{
	int rows = mR.rows();
	int columns = mR.cols();

	assert( rows == mI.rows() );
	assert( columns == mI.cols() );

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = complex<Type>( mR[i][j], mI[i][j] );

    return cA;
}


/**
 * Get magnitude of a complex matrix.
 */
template <typename Type>
Matrix<Type> abs( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = abs( A[i][j] );

    return tmp;
}


/**
 * Get angle of a complex matrix.
 */
template <typename Type>
Matrix<Type> arg( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = arg( A[i][j] );

    return tmp;
}


/**
 * Get real part of a complex matrix.
 */
template <typename Type>
Matrix<Type> real( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].real();

    return tmp;
}


/**
 * Get imaginary part of a complex matrix.
 */
template <typename Type>
Matrix<Type> imag( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].imag();

    return tmp;
}
