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
 *                               vector-impl.h
 *
 * Implementation for Vector class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * initialize
 */
template <typename Type>
void Vector<Type>::init( int length )
{
	assert( pv0 == NULL );
	pv0 = new Type[length];

	assert( pv0 != NULL );
	pv1 = pv0 - 1;
	nRow = length;
}


/**
 * copy vector from normal array
 */
template <typename Type>
inline void Vector<Type>::copyFromArray( const Type *v )
{
	for( int i=0; i<nRow; ++i )
		pv0[i] = v[i];
}


/**
 * set vector by a scalar
 */
template <typename Type>
inline void Vector<Type>::setByScalar( const Type &x )
{
	for( int i=0; i<nRow; ++i )
		pv0[i] = x;
}


/**
 * destroy the vector
 */
template <typename Type>
void Vector<Type>::destroy()
{
	if( pv0 == NULL )
		return;

	delete []pv0;

	pv0 = NULL;
	pv1 = NULL;
}


/**
 * constructors and destructor
 */
template <typename Type>
Vector<Type>::Vector()
: pv0(0), pv1(0), nRow(0)
{
}

template <typename Type>
Vector<Type>::Vector( const Vector<Type> &v )
: pv0(0), pv1(0), nRow(0)
{
	init( v.nRow );
	copyFromArray( v.pv0 );
}

template <typename Type>
Vector<Type>::Vector( int length, const Type &x )
:  pv0(0), pv1(0), nRow(0)
{
	init( length );
	setByScalar( x );
}

template <typename Type>
Vector<Type>::Vector( int length, const Type *array )
:  pv0(0), pv1(0), nRow(0)
{
	init( length );
	copyFromArray( array );
}

template <typename Type>
Vector<Type>::~Vector()
{
	destroy();
}


/**
 * overload evaluate operator= from vector to vector
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator=( const Vector<Type> &v )
{
	if( pv0 == v.pv0 )
		return *this;

	if( nRow == v.nRow )
		copyFromArray( v.pv0 );
	else
	{
		destroy();
		init( v.nRow );
		copyFromArray( v.pv0 );
	}

	return *this;
}


/**
 * overload evaluate operator= from scalar to vector
 */
template <typename Type>
inline Vector<Type>& Vector<Type>::operator=( const Type &x )
{
	setByScalar( x );

	return *this;
}


/**
 * overload operator [] for 0-offset access
 */
template <typename Type>
inline Type& Vector<Type>::operator[]( int i )
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return pv0[i];
}

template <typename Type>
inline const Type& Vector<Type>::operator[]( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return pv0[i];
}


/**
 * overload operator () for 1-offset access
 */
template <typename Type>
inline Type& Vector<Type>::operator()( int i )
{
#ifdef BOUNDS_CHECK
	assert( 1 <= i );
	assert( i <= nRow );
#endif

	return pv1[i];
}

template <typename Type>
inline const Type& Vector<Type>::operator()( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 1 <= i );
	assert( i <= nRow );
#endif

	return pv1[i];
}


/**
 * iterators
 */
template <typename Type>
inline typename Vector<Type>::iterator Vector<Type>::begin()
{
    return pv0;
}

template <typename Type>
inline typename Vector<Type>::const_iterator Vector<Type>::begin() const
{
    return pv0;
}

template <typename Type>
inline typename Vector<Type>::iterator Vector<Type>::end()
{
    return pv0 + nRow;
}

template <typename Type>
inline typename Vector<Type>::const_iterator Vector<Type>::end() const
{
    return pv0 + nRow;
}


/**
 * type conversion functions
 */
template <typename Type>
inline Vector<Type>::operator Type*()
{
	return pv0;
}

template <typename Type>
inline Vector<Type>::operator const Type*() const
{
	return pv0;
}


/**
 * get the vector's total size
 */
template <typename Type>
inline int Vector<Type>::size() const
{
	return  nRow;
}


/**
 * get the vector's dimension
 */
template <typename Type>
inline int Vector<Type>::dim() const
{
	return  nRow;
}


/**
 * reallocate vector's size
 */
template <typename Type>
Vector<Type>& Vector<Type>::resize( int length )
{
	if( nRow == length )
		return *this;

	destroy();
	init( length );

	return *this;
}


/**
 * compound assignment operators +=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator+=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ += x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator+=( const Vector<Type> &rhs )
{
    assert( nRow == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ += *itrR++;

	return *this;
}


/**
 * compound assignment operators -=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator-=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ -= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator-=( const Vector<Type> &rhs )
{
    assert( nRow == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ -= *itrR++;

	return *this;
}


/**
 * compound assignment operators *=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator*=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ *= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator*=( const Vector<Type> &rhs )
{
    assert( nRow == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ *= *itrR++;

	return *this;
}


/**
 * compound assignment operators /=
 */
template <typename Type>
Vector<Type>& Vector<Type>::operator/=( const Type &x )
{
    iterator itr = (*this).begin();
    while( itr != (*this).end() )
        *itr++ /= x;

	return *this;
}

template <typename Type>
Vector<Type>& Vector<Type>::operator/=( const Vector<Type> &rhs )
{
    assert( nRow == rhs.dim() );

    iterator itrL = (*this).begin();
    const_iterator itrR = rhs.begin();
    while( itrL != (*this).end() )
        *itrL++ /= *itrR++;

	return *this;
}


/**
 * Overload the output stream function.
 */
template <typename Type>
ostream& operator<<( ostream &out, const Vector<Type> &v )
{
	int N = v.dim();
	out << "size: " << N << " by 1" << "\n";

	for( int i=0; i<N; ++i )
		out << v[i] << " " << "\n";

	return out;
}


/**
 * Overload the input stream function.
 */
template <typename Type>
istream& operator>>( istream &in, Vector<Type> &v )
{
	int N;
	in >> N;

	if( !( N == v.dim() ) )
		v.resize( N );

	for( int i=0; i<N; ++i )
		in >> v[i];

	return in;
}


/**
 * get negative vector
 */
template <typename Type>
Vector<Type> operator-( const Vector<Type> &v )
{
	Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = -(*itrR++);

    return tmp;
}


/**
 * vector-scalar addition.
 */
template <typename Type>
inline Vector<Type> operator+( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp += x;
}

template <typename Type>
inline Vector<Type> operator+( const Type &x, const Vector<Type> &v )
{
	return v+x;
}


/**
 * vector-scalar substraction.
 */
template <typename Type>
inline Vector<Type> operator-( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp -= x;
}

template <typename Type>
inline Vector<Type> operator-( const Type &x, const Vector<Type> &v )
{
	Vector<Type> tmp( v );
	return -tmp += x;
}


/**
 * vector-scalar multiplication.
 */
template <typename Type>
inline Vector<Type> operator*( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp *= x;
}

template <typename Type>
inline Vector<Type> operator*( const Type &x, const Vector<Type> &v )
{
	return v*x;
}


/**
 * vector-scalar division.
 */
template <typename Type>
inline Vector<Type> operator/( const Vector<Type> &v, const Type &x )
{
	Vector<Type> tmp( v );
	return tmp /= x;
}

template <typename Type>
inline Vector<Type> operator/( const Type &x, const Vector<Type> &v )
{
	int N = v.dim();
	Vector<Type> tmp( N );

	for( int i=0; i<N; ++i )
		tmp[i] = x / v[i];

	return tmp;
}


/**
 * vector-vector addition.
 */
template <typename Type>
inline Vector<Type> operator+( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp += v2;
}


/**
 * vector-vector substraction.
 */
template <typename Type>
inline Vector<Type> operator-( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp -= v2;
}


/**
 * vector-vector multiplication.
 */
template <typename Type>
inline Vector<Type> operator*( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp *= v2;
}


/**
 * vector-vector division.
 */
template <typename Type>
inline Vector<Type> operator/( const Vector<Type> &v1, const Vector<Type> &v2 )
{
    Vector<Type> tmp( v1 );
	return tmp /= v2;
}


/**
 * Inner product for vectors.
 */
template <typename Type>
Type dotProd( const Vector<Type> &v1, const Vector<Type> &v2 )
{
	assert( v1.dim() == v2.dim() );

    Type sum = 0;
    typename Vector<Type>::const_iterator itr1 = v1.begin();
    typename Vector<Type>::const_iterator itr2 = v2.begin();

    while( itr1 != v1.end() )
		sum += (*itr1++) * (*itr2++);

	return sum;
}


/**
 * Inner product for vectors.
 */
template <typename Type>
complex<Type> dotProd( const Vector<complex<Type> > &v1,
                       const Vector<complex<Type> > &v2 )
{
	assert( v1.dim() == v2.dim() );

    complex<Type> sum = 0;
    typename Vector<complex<Type> >::const_iterator itr1 = v1.begin();
    typename Vector<complex<Type> >::const_iterator itr2 = v2.begin();

    while( itr1 != v1.end() )
		sum += (*itr1++) * conj(*itr2++);

	return sum;
}


/**
 * Vector's sum.
 */
template <typename Type>
Type sum( const Vector<Type> &v )
{
    Type sum = 0;
    typename Vector<Type>::const_iterator itr = v.begin();

    while( itr != v.end() )
		sum += *itr++;

	return sum;
}


/**
 * Minimum value of vector.
 */
template <typename Type>
Type min( const Vector<Type> &v )
{
    Type m = v[0];
    for( int i=1; i<v.size(); ++i )
        if( m > v[i] )
            m = v[i];

    return m;
}


/**
 * Maximum value of vector.
 */
template <typename Type>
Type max( const Vector<Type> &v )
{
    Type M = v[0];
    for( int i=1; i<v.size(); ++i )
        if( M < v[i] )
            M = v[i];

    return M;
}


/**
 * Vector's norm in Euclidean space.
 */
template <typename Type>
Type norm( const Vector<Type> &v )
{
	Type sum = 0;
	typename Vector<Type>::const_iterator itr = v.begin();

	while( itr != v.end() )
	{
	    sum += (*itr) * (*itr);
	    itr++;
	}

	return Type(sqrt(1.0*sum));
}


/**
 * Vector's norm in Euclidean space.
 */
template <typename Type>
Type norm( const Vector<complex<Type> > &v )
{
	Type sum = 0;
	typename Vector<complex<Type> >::const_iterator itr = v.begin();

	while( itr != v.end() )
	    sum += norm(*itr++);

	return Type(sqrt(1.0*sum));
}


/**
 * return vector's reversion
 */
template <typename Type>
void swap( Vector<Type> &lhs, Vector<Type> &rhs )
{
    typename Vector<Type>::iterator itrL = lhs.begin(),
                                    itrR = rhs.begin();

    while( itrL != lhs.end() )
        std::swap( *itrL++, *itrR++ );
}


/**
 * Generates a vector of n points linearly spaced between and
 * including a and b.
 */
template <typename Type>
Vector<Type> linspace( Type a, Type b, int n )
{
    if( n < 1 )
        return Vector<Type>();
    else if( n == 1 )
        return Vector<Type>( 1, a );
    else
    {
        Type dx = (b-a) / (n-1);

        Vector<Type> tmp(n);
        for( int i=0; i<n; ++i )
            tmp[i] = a + i*dx;

        return tmp;
    }
}


/**
 * Get magnitude of a complex vector.
 */
template <typename Type>
Vector<Type> abs( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = abs(*itrR++);

    return tmp;
}


/**
 * Get angle of a complex vector.
 */
template <typename Type>
Vector<Type> arg( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = arg(*itrR++);

    return tmp;
}


/**
 * Get real part of a complex vector.
 */
template <typename Type>
Vector<Type> real( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = (*itrR++).real();

    return tmp;
}


/**
 * Get imaginary part of a complex vector.
 */
template <typename Type>
Vector<Type> imag( const Vector< complex<Type> > &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector< complex<Type> >::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = (*itrR++).imag();

    return tmp;
}


/**
 * Convert real vector to complex vector.
 */
template <typename Type>
Vector<complex<Type> > complexVector( const Vector<Type> &rv )
{
	int N = rv.dim();

    Vector<complex<Type> > cv( N );
    typename Vector<complex<Type> >::iterator itrL = cv.begin();
    typename Vector<Type>::const_iterator itrR = rv.begin();

    while( itrR != rv.end() )
        *itrL++ = *itrR++;

    return cv;
}

template <typename Type>
Vector<complex<Type> > complexVector( const Vector<Type> &vR,
                                      const Vector<Type> &vI )
{
	int N = vR.dim();

	assert( N == vI.dim() );

    Vector<complex<Type> > cv( N );
    typename Vector<complex<Type> >::iterator itrC = cv.begin();
    typename Vector<Type>::const_iterator itrR = vR.begin(),
                                          itrI = vI.begin();

    while( itrC != cv.end() )
        *itrC++ = complex<Type>( *itrR++, *itrI++ );

    return cv;
}
