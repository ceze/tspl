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
 * Implementation for Vector math functions.
 *
 * Zhang Ming, 2010-08, Xi'an Jiaotong University.
 *****************************************************************************/


template <typename Type>
Vector<Type> abs( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = abs(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> cos( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = cos(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> sin( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = sin(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> tan( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = tan(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> acos( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = acos(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> asin( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = asin(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> atan( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = atan(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> exp( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = exp(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> log( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = log(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> log10( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = log10(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> sqrt( const Vector<Type> &v )
{
    Vector<Type> tmp( v.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = v.begin();

    while( itrL != tmp.end() )
        *itrL++ = sqrt(*itrR++);

    return tmp;
}


template <typename Type>
Vector<Type> pow( const Vector<Type> &b, const Vector<Type> &e )
{
    assert( b.dim() == e.dim() );

    Vector<Type> tmp( b.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR1 = b.begin();
    typename Vector<Type>::const_iterator itrR2 = e.begin();

    while( itrL != tmp.end() )
        *itrL++ = pow( *itrR1++, *itrR2++ );

    return tmp;
}


template <typename Type>
Vector<Type> pow( const Vector<Type> &b, const Type &e )
{
    Vector<Type> tmp( b.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = b.begin();

    while( itrL != tmp.end() )
        *itrL++ = pow( *itrR++, e );

    return tmp;
}


template <typename Type>
Vector<Type> pow( const Type &b, const Vector<Type> &e )
{
    Vector<Type> tmp( e.dim() );
    typename Vector<Type>::iterator itrL = tmp.begin();
    typename Vector<Type>::const_iterator itrR = e.begin();

    while( itrL != tmp.end() )
        *itrL++ = pow( b, *itrR++ );

    return tmp;
}


/**
 * Normal distribution with expectation "u" and variance "r".
 */
template <typename Type>
Vector<Type> gauss( const Vector<Type> &x, const Type &u, const Type &r )
{
    Vector<Type> tmp(x);

    tmp = (tmp-u)*(tmp-u) / ( -2*r*r );
    tmp = exp(tmp) / Type( (sqrt(2*PI)*r) );

    return tmp;
}
