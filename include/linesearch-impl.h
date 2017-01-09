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
 *                               linesearch-impl.h
 *
 * Implementation for LineSearch class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template <typename Dtype, typename Ftype>
LineSearch<Dtype, Ftype>::LineSearch()
{
    funcNum = 0;
    success = true;
}

template <typename Dtype, typename Ftype>
LineSearch<Dtype, Ftype>::~LineSearch()
{
}


/**
 * Finding the step size at point "xk" in direction of "dk" of function
 * "func". The default heuristics number is "maxItr=10".
 */
template <typename Dtype, typename Ftype>
Dtype LineSearch<Dtype, Ftype>::getStep( Ftype &func, Vector<Dtype> &xk,
                                         Vector<Dtype> &dk, int maxItr )
{
    // Set line search parameters that everyone uses.
    Dtype mu = Dtype(0.001),
          kUp = Dtype(0.5),
          kLow = Dtype(0.1),
          alpha = Dtype(1.0),
          alphaMin,
          alphaMax;

    Dtype fNew,
          fk = func(xk);

    Vector<Dtype> xNew,
                  gk = func.grad(xk);

    Dtype gd = dotProd( gk, dk );

    for( int i=0; i<maxItr; ++i )
    {
        xNew = xk + alpha*dk;
        fNew = func(xNew);
        funcNum++;

        if( fNew < fk+mu*alpha*gd )
        {
            success = true;
            return alpha;
        }
        else
        {
            alphaMin = kLow*alpha;
            alphaMax = kUp*alpha;

            // Compute the step by using quadratic polynomial interpolation.
            alpha = Dtype(-0.5)*alpha*alpha*gd / ( fNew-fk-alpha*gd );

            // bound checking
            if( alpha < alphaMin )
                alpha = alphaMin;
            else if( alpha > alphaMax )
                alpha = alphaMax;
        }
    }

    if( fNew>=fk )
    {
        success = false;
        return Dtype(0.0);
    }
    else
    {
        success = true;
        return alpha;
    }
}


/**
 * Get the number of objective function's calculation.
 */
template <typename Dtype, typename Ftype>
inline int LineSearch<Dtype, Ftype>::getFuncNum() const
{
    return funcNum;
}


/**
 * Judgement whether the optimal solution is found or not.
 */
template <typename Dtype, typename Ftype>
inline bool LineSearch<Dtype, Ftype>::isSuccess() const
{
    return success;
}
