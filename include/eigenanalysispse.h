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
 *                              eigenanalysispse.h
 *
 * Eigenanalysis algorithms for spectrum estimation.
 *
 * The eigenanalysis algorithms perform eigen decomposition of the signal's
 * auto-correlation matrix (ACM) to estimate signal's frequency content. The
 * ACM can be decomposed into signal-subspace and noise-subspace, then use
 * the orthogonality of the tow subspaces to estimate the signal's spectrum.
 *
 * These algorithms, such as Pisarenko, MUSIC (MUltiple SIgnal Classification)
 * and ESPRIT (Estimation of Signal Parameters by Rotational Invariance
 * Techniques) are particularly suitable for signals that are the sum of
 * sinusoids with additive white Gaussian noise. The model order (the number
 * of complex exponential signal) is estimated by minimizing the MDL criterion.
 *
 * This file also provide the Capon's maximum likehood method or minimum
 * variance method for specturm estimation.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef EIGENANALYSISPSE_H
#define EIGENANALYSISPSE_H


#include <toeplitz.h>
#include <levinson.h>
#include <linequs1.h>
#include <svd.h>
#include <evd.h>


namespace splab
{

    template<typename Type> Vector<Type> caponPSE( const Vector<Type>&,
                                                   int, int );
    template<typename Type> Vector<Type> pisarenkoPSE( const Vector<Type>&,
                                                       int, int, int );
    template<typename Type> Vector<Type> musicPSE( const Vector<Type>&,
                                                   int, int, int );
    template<typename Type> Vector<Type> espritPSE( const Vector<Type>&,
                                                    int, int );

    template<typename Type> int orderEst( const Vector<Type>&, int );


    #include <eigenanalysispse-impl.h>

}
// namespace splab


#endif
// EIGENANALYSISPSE_H
