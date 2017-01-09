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
 *                               bfgs-impl.h
 *
 * Implementation for BFGS class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Dtype, typename Ftype>
BFGS<Dtype, Ftype>::BFGS() : LineSearch<Dtype, Ftype>()
{
}

template <typename Dtype, typename Ftype>
BFGS<Dtype, Ftype>::~BFGS()
{
}


/**
 * Finding the optimal solution. The default tolerance error and maximum
 * iteratin number are "tol=1.0e-6" and "maxItr=100", respectively.
 */
template <typename Dtype, typename Ftype>
void BFGS<Dtype, Ftype>::optimize( Ftype &func, Vector<Dtype> &x0,
                                   Dtype tol, int maxItr )
{
    // initialize parameters.
    int k = 0,
        cnt = 0,
        N = x0.dim();

    Dtype ys,
          yHy,
          alpha;
    Vector<Dtype> d(N),
                  s(N),
                  y(N),
                  v(N),
                  Hy(N),
                  gPrev(N);
    Matrix<Dtype> H = eye( N, Dtype(1.0) );

    Vector<Dtype> x(x0);
    Dtype fx = func(x);
    this->funcNum++;
    Vector<Dtype> gnorm(maxItr);
    Vector<Dtype> g = func.grad(x);
    gnorm[k++]= norm(g);

    while( ( gnorm(k) > tol ) && ( k < maxItr ) )
    {
        // descent direction
        d = - H * g;

        // one dimension searching
        alpha = this->getStep( func, x, d );

        // check flag for restart
        if( !this->success )
            if( norm(H-eye(N,Dtype(1.0))) < EPS )
                break;
            else
            {
                H = eye( N, Dtype(1.0) );
                cnt++;
                if( cnt == maxItr )
                    break;
            }
        else
        {
            // update
            s = alpha * d;
            x += s;
            fx = func(x);
            this->funcNum++;
            gPrev = g;
            g = func.grad(x);
            y = g - gPrev;

            Hy = H * y;
            ys = dotProd( y, s );
            yHy = dotProd( y, Hy );
            if( (ys < EPS) || (yHy < EPS) )
                H = eye( N, Dtype(1.0) );
            else
            {
                v = sqrt(yHy) * ( s/ys - Hy/yHy );
                H = H + multTr(s,s)/ys - multTr(Hy,Hy)/yHy + multTr(v,v);
            }
            gnorm[k++] = norm(g);
        }
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
inline Vector<Dtype> BFGS<Dtype, Ftype>::getOptValue() const
{
    return xOpt;
}


/**
 * Get the norm of gradient in each iteration.
 */
template <typename Dtype, typename Ftype>
inline Vector<Dtype> BFGS<Dtype, Ftype>::getGradNorm() const
{
    return gradNorm;
}


/**
 * Get the minimum value of objective function.
 */
template <typename Dtype, typename Ftype>
inline Dtype BFGS<Dtype, Ftype>::getFuncMin() const
{
    return fMin;
}


/**
 * Get the iteration number.
 */
template <typename Dtype, typename Ftype>
inline int BFGS<Dtype, Ftype>::getItrNum() const
{
    return gradNorm.dim()-1;
}
