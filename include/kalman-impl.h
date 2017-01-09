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
 *                              kalman-impl.h
 *
 * Implementation for Kalman Filter.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * The simple Kalman filter for one step.
 * A ---> system matrix defining linear dynamic system
 * C ---> measurement matrix defining relationship between system's state
 *        and measurements
 * Q ---> covariance matrix of process noise in system state dynamics
 * R ---> covariance matrix of measurements uncertainty
 * y ---> measurements at time t
 * xPrev ---> previous estimated state vector of the linear dynamic system
 * initDiagV ---> diagonal vector for initializing the covariance matrix of
 *                state estimation uncertainty
 */
template <typename Type>
Vector<Type> kalman( const Matrix<Type> &A, const Matrix<Type> &C,
                     const Matrix<Type> &Q, const Matrix<Type> &R,
                     const Vector<Type> &xPrev, const Vector<Type> &y,
                     const Vector<Type> &initDiagV )
{
    int N = xPrev.size();

    // covariance matrix of state estimation uncertainty
    static Matrix<Type> V = diag(initDiagV);

    // previoused state vector
    Vector<Type> xPred = A * xPrev;

    // inovation
    Vector<Type> alpha = y - C * xPred;

    Matrix<Type> CTran = trT( C );
    Matrix<Type> VPred = A*V*trT(A) + Q;

    // Kalman gain matrix
    Matrix<Type> KGain = VPred*CTran * inv(C*VPred*CTran+R);

    V = ( eye(N,Type(1.0)) - KGain*C ) * VPred;

    // return the estimation of the state vector
    return xPred + KGain * alpha;
}
