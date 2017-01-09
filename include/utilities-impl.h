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
 *                               utilities-impl.h
 *
 * Implementation for utilities.
 *
 * Zhang Ming, 2010-01 (revised 2010-08), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Modulus after division. return a integer in the range of 0 to n-1.
 * e.g. -1%5=-1, mod(-1,5)=4
 */
int mod( int m, int n )
{
    if( n != 0 )
    {
        int r = m % n;
        if( r < 0 )
            r += n;

        return r;
    }
    else
    {
        cerr << "The dividend shouldn't be zero." << endl;
        return 0;
    }
}


/**
 * Rounds the elements of a/b to the nearest integers
 * greater than or equal to a/b.
 * e.g. ceil(10,2) = 5, ceil(10,3)=4.
 */
int ceil( int m, int n )
{
    if( n != 0 )
    {
        int q = m / n;
        if( m%n != 0 )
            q += 1;

        return q;
    }
    else
    {
        cerr << "The dividend shouldn't be zero." << endl;
        return 0;
    }
}


/**
 * Flip vector left to right.
 */
template <typename Type>
Vector<Type> reverse( const Vector<Type> &v )
{
    Vector<Type> tmp( v.size() );
    typename Vector<Type>::iterator ib = tmp.begin();
    typename Vector<Type>::const_iterator ie = v.end();

    while( ib != tmp.end() )
        *ib++ = *--ie;

    return tmp;
}

template <typename Type>
inline Vector<Type> flip( const Vector<Type> &v )
{
    return reverse( v );
}


/**
 * Vector's shift.
 */
template <typename Type>
Vector<Type> shift( const Vector<Type> &v, int shiftsize )
{
    Vector<Type> tmp( v.dim() );

    if( shiftsize >= 0 )
    {
        typename Vector<Type>::iterator itrL = tmp.begin()+shiftsize;
        typename Vector<Type>::const_iterator itrR = v.begin();
        while( itrL != tmp.end() )
            *itrL++ = *itrR++;
    }
    else
    {
        typename Vector<Type>::iterator itrL = tmp.begin();
        typename Vector<Type>::const_iterator itrR = v.begin()-shiftsize;
        while( itrR != v.end() )
            *itrL++ = *itrR++;
    }

    return tmp;
}


/**
 * Vector's circulary shift.
 */
template <typename Type>
Vector<Type> circshift( const Vector<Type> &v, int shiftsize )
{
    Vector<Type> tmp( v.dim() );

    if( shiftsize >= 0 )
    {
        typename Vector<Type>::iterator itrL = tmp.begin()+shiftsize;
        typename Vector<Type>::const_iterator itrR = v.begin();
        while( itrL != tmp.end() )
            *itrL++ = *itrR++;

        itrL = tmp.begin();
        while( itrR != v.end() )
            *itrL++ = *itrR++;
    }
    else
    {
        typename Vector<Type>::iterator itrL = tmp.begin();
        typename Vector<Type>::const_iterator itrR = v.begin()-shiftsize;
        while( itrR != v.end() )
            *itrL++ = *itrR++;

        itrR = v.begin();
        while( itrL != tmp.end() )
            *itrL++ = *itrR++;
    }

    return tmp;
}


/**
 * Vector's fft shift.
 */
template <typename Type>
inline Vector<Type> fftshift( const Vector<Type> &v )
{
    int shiftsize = v.dim() - v.dim()/2 - 1;
    return circshift( v, shiftsize );
}


/**
 * dyadic upsampling
 * w = dyadup(v, evenodd), where v is a vector, returns an extended
 * copy of vector v obtained by inserting zeros. Whether the zeros
 * are inserted as even- or odd-indexed elements of w depends on the
 * value of positive integer evenodd:
 * If evenodd is even, then w[2k]=0, w[2k+1]=v[k], w.size()=2*v.size()+1.
 * If evenodd is odd, then w[2k]=v[k], w[2k+1]=0. w.size()=2*v.size()-1.
 */
template <typename Type>
Vector<Type> dyadUp( const Vector<Type> &v, int evenodd )
{
    int length = v.dim();
    Vector<Type> tmp;

    if( evenodd%2 == 0 )
    {
        tmp.resize( 2*length+1 );
        for( int i=0; i<length; ++i )
        {
            tmp[2*i] = 0;
            tmp[2*i+1] = v[i];
        }
        tmp[2*length] = 0;
    }
    else
    {
        tmp.resize( 2*length-1 );
        for( int i=0; i<length-1; ++i )
        {
            tmp[2*i] = v[i];
            tmp[2*i+1] = 0;
        }
        tmp[2*length-2] = v[length-1];
    }

    return tmp;
}


/**
 * dyadic downsampling
 * w = dyadup(v, evenodd), where v is a vector, returns a version of v
 * that has been downsampled by 2. Whether w contains the even or odd
 * indexed samples of v depends on the value of positive integer evenodd:
 * If evenodd is even, then w[k]=v[2*k], w.size()=(v.size()+1)/2.
 * If evenodd is odd, then w[k]=v[2*k+1], w.size()=v.size()/2.
 */
template <typename Type>
Vector<Type> dyadDown( const Vector<Type> &v, int evenodd )
{
    int length = v.dim();
    Vector<Type> tmp;

    if( evenodd%2 == 0 )
    {
        tmp.resize( (length+1)/2 );
        for( int i=0; i<tmp.dim(); ++i )
            tmp[i] = v[2*i];
    }
    else
    {
        tmp.resize( length/2 );
        for( int i=0; i<tmp.dim(); ++i )
            tmp[i] = v[2*i+1];
    }

    return tmp;
}


/**
 * Real signal interpolation by the method of padding zeros in frequency domain.
 * The interpolation factor should be >= 1.
 */
template <typename Type>
Vector<Type> fftInterp( const Vector<Type> &sn, int factor )
{
    int N = sn.size(),
        halfN = N/2,
        offset = (factor-1)*N;

    Vector< complex<Type> > Sk = fft(sn);
    Vector< complex<Type> > Xk(factor*N);

    for( int i=0; i<=halfN; ++i )
        Xk[i] = Type(factor)*Sk[i];
    for( int i=halfN+1; i<N; ++i )
        Xk[offset+i] = Type(factor)*Sk[i];

    return ifftc2r(Xk);
}


/**
 * Complex signal interpolation by the method of padding zeros in frequency domain.
 * The interpolation factor should be >= 1.
 */
template <typename Type>
Vector< complex<Type> > fftInterp( const Vector< complex<Type> > &sn,
                                   int factor )
{
    int N = sn.size(),
        halfN = N/2,
        offset = (factor-1)*N;

    Vector< complex<Type> > Sk = fft(sn);
    Vector< complex<Type> > Xk(factor*N);

    for( int i=0; i<=halfN; ++i )
        Xk[i] = Type(factor)*Sk[i];
    for( int i=halfN+1; i<N; ++i )
        Xk[offset+i] = Type(factor)*Sk[i];

    return ifft(Xk);
}


/**
 * Keep part of vector.
 * For a vector, w = wkeep(v,L,opt) extracts the vector w from the vector v.
 * The length of w is L. If direction = "center" ("left", "rigth",
 * respectively), w is the central (left, right, respectively) part of v.
 * w = wkeep(x,L) is equivalent to w = wkeep(v,L,"center").
 * w = wkeep(v,L,first) returns the vector v[first] to v[first+L-1].
 */
template <typename Type>
Vector<Type> wkeep( const Vector<Type> &v, int length, int first )
{
    Vector<Type> tmp(length);

    if( ( 0 < length ) && ( length <= v.dim()-first ) )
    {
        for( int i=0; i<length; ++i )
            tmp[i] = v[first+i];

        return tmp;
    }
    else
    {
        cerr << "Invalid length input." << endl;
        return tmp;
    }
}

template <typename Type>
Vector<Type> wkeep( const Vector<Type> &v, int length,
                    const string &direction )
{
    int lv = v.dim();
    Vector<Type> tmp(length);

    if( ( 0 <= length ) && ( length <= lv ) )
    {
        if( direction == "right" )
            for( int i=0; i<length; ++i )
                tmp[i] = v[lv-length+i];
        else if( direction == "left" )
            for( int i=0; i<length; ++i )
                tmp[i] = v[i];
        else
        {
            int first = (lv-length)/2;
            for( int i=0; i<length; ++i )
                tmp[i] = v[first+i];
        }

        return tmp;
    }
    else
    {
        cerr << "Invalid length input." << endl;
        return tmp;
    }
}


/**
 * extend vector
 * The extension types are specified by the string "direction", include
 * "left", "right" and "both". The default type is "both". The valid
 * extension modes, which specified by strint "mode" are: zero padding
 * ("zpd"), periodized extension("ppd") and symetirc extension("sym").
 * The default mode is "zpd".
 */
template <typename Type>
Vector<Type> wextend( const Vector<Type> &v, int extLength,
                      const string &direction, const string &mode )
{
    if( extLength >= 0 )
    {
        Vector<Type> tmp;
        int lv = v.dim();

        if( direction == "right" )
        {
            tmp.resize( lv+extLength );
            for( int i=0; i<lv; ++i )
                tmp[i] = v[i];

            if( mode == "sym" )
                for( int i=0; i<extLength; ++i )
                    tmp[lv+i] = v[lv-1-i];
            else if( mode == "ppd" )
                for( int i=0; i<extLength; ++i )
                    tmp[lv+i] = v[i];
            else
                for( int i=0; i<extLength; ++i )
                    tmp[lv+i] = 0;
        }
        else if( direction == "left" )
        {
            tmp.resize( lv+extLength );

            if( mode == "sym" )
                for( int i=0; i<extLength; ++i )
                    tmp[i] = v[extLength-1-i];
            else if( mode == "ppd" )
                for( int i=0; i<extLength; ++i )
                    tmp[i] = v[lv-extLength+i];
            else
                for( int i=0; i<extLength; ++i )
                    tmp[i] = 0;

            for( int i=0; i<lv; ++i )
                tmp[i+extLength] = v[i];
        }
        else
        {
            tmp.resize( lv+2*extLength );
            for( int i=0; i<lv; ++i )
                tmp[i+extLength] = v[i];

            if( mode == "sym" )
                for( int i=0; i<extLength; ++i )
                {
                    tmp[i] = v[extLength-1-i];
                    tmp[lv+extLength+i] = v[lv-1-i];
                }
            else if( mode == "ppd" )
                for( int i=0; i<extLength; ++i )
                {
                    tmp[i] = v[lv-extLength+i];
                    tmp[lv+extLength+i] = v[i];
                }
            else
                for( int i=0; i<extLength; ++i )
                {
                    tmp[i] = 0;
                    tmp[lv+extLength+i] = 0;
                }
        }
        return tmp;
    }
    else
    {
        cerr << "The extesion length should be greater zero." << endl;
        return Vector<Type>(0);
    }
}
