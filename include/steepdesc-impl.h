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
 *                               steepdesc-impl.h
 *
 * Implementation for SteepDesc class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Dtype, typename Ftype>
SteepDesc<Dtype, Ftype>::SteepDesc() : LineSearch<Dtype, Ftype>()
{
}

template <typename Dtype, typename Ftype>
SteepDesc<Dtype, Ftype>::~SteepDesc()
{
}


/**
 * Finding the optimal solution. The default tolerance error and maximum
 * iteratin number are "tol=1.0e-6" and "maxItr=100", respectively.
 */
template <typename Dtype, typename Ftype>
void SteepDesc<Dtype, Ftype>::optimize( Ftype &func, Vector<Dtype> &x0,
                                        Dtype tol, int maxItr )
{
    // initialize parameters.
    int k = 0;
    Vector<Dtype> x(x0);
    Dtype fx = func(x);
    this->funcNum++;
    Vector<Dtype> gnorm(maxItr);
    Vector<Dtype> g = func.grad(x);
    gnorm[k++]= norm(g);

    Dtype alpha;
    Vector<Dtype> d(g.dim());

    while( ( gnorm(k) > tol ) && ( k < maxItr ) )
    {
        // descent direction
        d =  - g;

        // one dimension searching
        alpha = this->getStep( func, x, d );
        if( !this->success )
            break;

		// update
        x += alpha*d;
        fx = func(x);
        this->funcNum++;
        g = func.grad(x);
        gnorm[k++] = norm(g);
    }

    xOpt = x;
    fMin = fx;
    gradNorm.resize(k);
    for( int i=0; i<k; ++i )
        gradNorm[i] = gnorm[i];

    if( gradNorm[k-1] > tol )
        this->success = false;
}


/**
 * Get the optimum point.
 */
template <typename Dtype, typename Ftype>
inline Vector<Dtype> SteepDesc<Dtype, Ftype>::getOptValue() const
{
    return xOpt;
}


/**
 * Get the norm of gradient in each iteration.
 */
template <typename Dtype, typename Ftype>
inline Vector<Dtype> SteepDesc<Dtype, Ftype>::getGradNorm() const
{
    return gradNorm;
}


/**
 * Get the minimum value of objective function.
 */
template <typename Dtype, typename Ftype>
inline Dtype SteepDesc<Dtype, Ftype>::getFuncMin() const
{
    return fMin;
}


/**
 * Get the iteration number.
 */
template <typename Dtype, typename Ftype>
inline int SteepDesc<Dtype, Ftype>::getItrNum() const
{
    return gradNorm.dim()-1;
}
